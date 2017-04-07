suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("-f", "--file"), 
              type="character",
              help="Path of FASTQ file to process"),
  make_option(c("-t", "--trim"),
              type="integer",
              default=0,
              help="Pre-trim all reads to this length before processing"),
  make_option(c("-k", "--keep"),
              type="integer",
              default=18,
              help="Minimum number of bases required after barcode to keep read"),
  make_option(c("-b", "--barcode"),
              type="character",
              default="CTGA",
              help="Barcode sequences (comma-separated) that follow random barcode"),
  make_option(c("-r", "--randombarcode"),
             type="integer",
             default=5,
             help="Number of bases at the start of each read used for random barcode"),
  make_option(c("-c", "--chunksize"),
              type="integer",
              default=1000,
              help="Number of reads to process at once (in thousands)"),
  make_option(c("-o", "--output"),
              type="character",
              help="Output FASTQ file (gzip compressed)"),
  make_option(c("-p", "--processors"),
              type="character",
              default=2,
              help="Number of simultaneous processing cores to utilize"))
        

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(stringr))

pn <- function(value) {
  prettyNum(value, big.mark=",")
}

# From HTSeqGenie

##' Scheduled parallel processing
##'
##' @param inext A function (without argument) returning an object to process; NULL if none left; this function is run in the main thread
##' @param fun Function to process the object returned by inext; this function is run in children thread
##' @param max.parallel.jobs Number of jobs to start in parallel 
##' @param ... Further arguments passed to fun
##' @param stop.onfail Throw error if one 
##' @param tracefun Callback function that will be executed in a separate thread 
##' @param tracefun.period Time intervall between calls to tracefun
##' @return Return value of applied function
##' @export
##' @keywords internal
sclapply <- function(inext, fun, max.parallel.jobs, ..., stop.onfail=TRUE, tracefun=NULL, tracefun.period=60) {
  ## initialise scheduler
  sjobs <- character(0) ## jobs (state)
  rjobs <- list() ## jobs (result)
  pnodes <- vector(mode="list", max.parallel.jobs) ## nodes (process)
  jnodes <- rep(NA, max.parallel.jobs) ## nodes (job id)

  ## cleanup procedure, based on mclapply
  cleanup <- function() {
    z <- pnodes[!sapply(pnodes, is.null)]
    if (length(z)>0) mccollect(z, wait=FALSE, timeout=4)
    z <- pnodes[!sapply(pnodes, is.null)]
    if (length(z)>0) mccollect(z, wait=FALSE, timeout=4)
    z <- pnodes[!sapply(pnodes, is.null)]
    if (length(z)>0) try(parallel:::mckill(z, tools::SIGTERM), silent=TRUE)
    z <- pnodes[!sapply(pnodes, is.null)]
    if (length(z)>0) mccollect(z, wait=FALSE, timeout=4)
    z <- pnodes[!sapply(pnodes, is.null)]
    if (length(z)>0) mccollect(z, wait=FALSE, timeout=4)
  }
  on.exit(cleanup())

  ## start scheduler
  collect.timeout <- 1 ## 1 seconds wait between each iteration
  inextdata <- NULL
  i <- 0
  repeat {
    ## is there a new job to process?
    if (is.null(inextdata)) inextdata <- inext()
    
    ## are all the jobs done?
    if (is.null(inextdata)) {
      if (length(sjobs)==0) break ## no jobs have been run
      if (all(sjobs=="done")) break
    }
    
    ## fire fun(inextdata) on node i
    repeat {
      ## periodic trace
      if (!is.null(tracefun)) {
        current <- proc.time()["elapsed"]
        if (!exists("last.tracefun") || (current-last.tracefun)>tracefun.period) {
          tracefun(type="timer", sjobs=sjobs, jnodes=jnodes)
          last.tracefun <- current
        }
      }
      
      i <- (i %% length(pnodes))+1
      process <- pnodes[[i]]
      if (!is.null(process)) {
        status <- mccollect(process, wait=FALSE, timeout=collect.timeout) ## wait collect.timeout seconds        
        if (is.null(status)) {
          ## node busy
          fire <- FALSE
        }
        else {
          ## node done: save results and fire a new job
          mccollect(process) ## kill job
          ## did an error occur?
          if (class(status[[1]])=="try-error") { 
            if (!is.null(tracefun)) {
              tracefun("error", sjobs=sjobs, jnodes=jnodes, chunkid=jnodes[i], error=status[[1]])
            }
            if (stop.onfail) {
              stop(paste("tools.R/sclapply: error in chunkid=", jnodes[i], ": ", status[[1]], sep=""))
            }
          }
          if (!is.null(tracefun)) tracefun("done", chunkid=jnodes[i])
          rjobs[jnodes[i]] <- status
          sjobs[jnodes[i]] <- "done"
          fire <- TRUE
        }
      } else {
        ## virgin node
        fire <- TRUE
      }
      
      ## fire a new job
      if (fire && !is.null(inextdata)) {
        jnodes[i] <- length(sjobs)+1
        sjobs[jnodes[i]] <- "running"
        if (!is.null(tracefun)) tracefun("start", chunkid=jnodes[i])
        pnodes[[i]] <- mcparallel(fun(inextdata, ...), mc.set.seed=TRUE)
        inextdata <- NULL
        break
      }

      ## no more job
      if (is.null(inextdata)) break
    }
  }
  
  rjobs
}

process_chunk <- function(fq_chunk, opt) {
    
    if(opt$trim > 0) {
      fq_chunk <- narrow(fq_chunk, start=1, end=pmin(width(fq_chunk), opt$trim))
    }
  
  output_file <- opt$output
  
  barcodes <- strsplit(opt$barcode, split=",")[[1]]
  if(length(unique(str_length(barcodes))) != 1) {
    stop("Fixed barcodes must all be the same length: ", barcodes)
  }
  
  barcode_start <- opt$randombarcode + 1
  barcode_end   <- barcode_start + str_length(barcodes[1]) - 1
  
  fq.lockfile <- paste0(output_file, ".lock")
  
  barcode_reads <- narrow(sread(fq_chunk), start=barcode_start, end=barcode_end)

  bc_matches <- as.list(rep(NA, times=length(barcodes)))
  names(bc_matches) <- barcodes
  
  for(barcode in barcodes) {
    matches <- elementLengths(vmatchPattern(barcode, barcode_reads, fixed=FALSE)) == 1
    n_count <- elementLengths(vmatchPattern("N", barcode_reads, fixed=TRUE))
    bc_matches[[barcode]] <- which(matches == TRUE & n_count <= 1)
  }
  
  # don't allow a read to match multiple barcodes
  previous_matches <- c()
  for(i in seq_along(bc_matches)) {
    bc_matches[[i]] <- bc_matches[[i]][!bc_matches[[i]] %in% previous_matches]
    previous_matches <- c(previous_matches, bc_matches[[i]])
  }
  
  total_matches <- sum(elementLengths(bc_matches))
  

  message("[", opt$file, "] ", pn(length(fq_chunk)), " reads with ", pn(total_matches), " barcode matches")
  
  if(total_matches > 0) {
    for(barcode in barcodes) {
      fq.matched <- fq_chunk[bc_matches[[barcode]]]
      if(length(fq.matched) == 0) next
        
      # Reject reads that are too short
      fq.matched <- fq.matched[width(sread(fq.matched)) >= barcode_end + opt$keep]

      # Keep random barcode
      random_bc  <- substr(sread(fq.matched), 1, opt$randombarcode)

      fq.matched <- narrow(fq.matched, start=barcode_end + 1, width=width(fq.matched) - (barcode_end + 1))

      fq.new <- ShortReadQ(sread   = sread(fq.matched),
                           quality = quality(fq.matched),
                           id      = BStringSet(paste0(random_bc, "_", barcode)))

      lock_status <- system(paste("lockfile", "-1", fq.lockfile), intern=FALSE)
      if(lock_status > 0) stop("lockfile command failed.")
      output_mode <- ifelse(file.exists(output_file), "a", "w")
      writeFastq(fq.new, file=output_file, compress=TRUE, mode=output_mode)
      file.remove(fq.lockfile)
    }
  }
  
  TRUE
}

yieldHelper <- function() {
  fq <- yield(fqstream, withIds=FALSE)
  if(length(fq) > 0) {
    fq
  } else {
    NULL
  }
}

# Restrict FastqStreamer threading
nothing <- .Call(ShortRead:::.set_omp_threads, 1L)

if(file.exists(opt$output)) stop("Output file ", opt$output, " already exists.")

fqstream <- FastqStreamer(opt$file, n=opt$chunksize * 1000)

results <- sclapply(yieldHelper, process_chunk, max.parallel.jobs=opt$processors, opt)
close(fqstream)

