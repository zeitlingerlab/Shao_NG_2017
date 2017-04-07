suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

# bamtoolsr.r - some commands for processing BAM files

option_list <- list(
  make_option(c("-f", "--file"), 
              type="character",
              default=NA,
              help="Path of BAM file to process"),			  
  make_option(c("-e", "--extension"),
              type="character",
              default="native",
              help="Extension length ('native' for no extension, 'auto' to estimate fragment size)"),			  
  make_option(c("-n", "--name"),
              type="character",
              help="R variable name for resulting coverage and ranges objects"),			               
  make_option(c("-p", "--paired"),
              action="store_true",
              default=FALSE,
              help="Paired-end data"),
  make_option(c("-s", "--skipfilter"),
              action="store_true",
              default=FALSE,
              help="Skip artifact filter"),
  make_option(c("-c", "--cores"),
              type="integer",
              default=2,
              help="Number of processor cores to use when performing artifact filtering")
  )


readBAM <- function(filepath) {
  bam.gr <- granges(readGAlignments(filepath))
  bam.gr
}

readBAM_pe <- function(filepath) {
  bam.gr <- granges(readGAlignmentPairs(filepath))
  strand(bam.gr) <- "*"
  bam.gr
}

artifact_filter <- function(gr, dup_count_limit=10, ext_length, cores=2) {
  message(id, "artifact_filter: Starting with ", pn(length(gr)), " reads")
  
  grl <- split(gr, seqnames(gr))
  grl <- mclapply(grl[elementLengths(grl) > 0], artifact_filter_base, dup_count_limit, ext_length, mc.cores=cores, mc.preschedule=FALSE)
  gr <- do.call(c, unname(grl))
  message(id, "artifact_filter: Returning ", pn(length(gr)), " reads after filtering")
  gr
}

artifact_filter_base <- function(gr, dup_count_limit=10, ext_length) {
  current_chr <- seqnames(gr)[1]
  message(id, "  ", current_chr, " with ", pn(length(gr)), " reads")
  uniq.gr <- unique(gr)
  uniq.gr$count <- countOverlaps(uniq.gr, gr, type="equal")
  
  uniq.above <- subset(uniq.gr, count >= .(dup_count_limit))
  uniq.below <- subset(uniq.gr, count <  .(dup_count_limit))

  message(id, "  ", current_chr, " above limit: ", pn(length(uniq.above)))
  message(id, "  ", current_chr, " below limit: ", pn(length(uniq.below)))

  if(length(uniq.above) > 0) {

    i_pos <- which(strand(uniq.above) == "+")
    i_neg <- which(strand(uniq.above) == "-")

    gr <- trim(resize(gr, 1))

    shift_value <- rep(ext_length, times=length(uniq.above))
    if(length(i_neg) > 0) shift_value[i_neg] <- shift_value[i_neg] * -1

    count_regions.gr <- trim(shift(resize(uniq.above, 1), shift_value))
    count_regions.gr <- trim(resize(count_regions.gr, 51, fix="center"))
  
  
    uniq.above$paired_count <- 0
  
    if(length(i_pos)) uniq.above$paired_count[i_pos] <- countOverlaps(count_regions.gr[i_pos], gr[strand(gr) == "-"], type="any", ignore.strand=TRUE)
    if(length(i_neg)) uniq.above$paired_count[i_neg] <- countOverlaps(count_regions.gr[i_neg], gr[strand(gr) == "+"], type="any", ignore.strand=TRUE)

    uniq.above$new_count <- with(mcols(uniq.above), pmax(pmin(count, paired_count), 1))
  
    expand.above <- GRanges(ranges     = IRanges(start = rep(start(uniq.above), times=uniq.above$new_count),
                                                 end   = rep(end(uniq.above),   times=uniq.above$new_count)),
                            strand     = rep(strand(uniq.above),   times=uniq.above$new_count),
                            seqnames   = rep(seqnames(uniq.above), times=uniq.above$new_count),
                            seqlengths = seqlengths(gr))

    expand.below <- GRanges(ranges     = IRanges(start=rep(start(uniq.below), times=uniq.below$count),
                                                 end=rep(end(uniq.below),     times=uniq.below$count)),
                            strand     = rep(strand(uniq.below),   times=uniq.below$count),
                            seqnames   = rep(seqnames(uniq.below), times=uniq.below$count),
                            seqlengths = seqlengths(gr))

    filtered.gr <- c(expand.above, expand.below)
    message(id, "  ", current_chr, " returning ", pn(length(filtered.gr)), " reads")
    filtered.gr
  } else {
    message(id, "  ", current_chr, " returning ", pn(length(gr)), " reads")
    return(gr)
  }
}

pn <- function(value) {
  prettyNum(value, big.mark=",")
}

opt <- parse_args(OptionParser(option_list=option_list))

if(is.na(opt$file)) {
  message("No BAM file specified. Use --help to list available options.")
  q(status=1)
}

suppressPackageStartupMessages(library(Rsamtools, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(chipseq, warn.conflicts=F))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts=F))
suppressPackageStartupMessages(library(GenomicAlignments))

id <- "[unknown]"

bam_file      <- opt$file
ext_length    <- opt$extension
var_name      <- opt$name
run_filter    <- !opt$skipfilter

id <- paste("[", var_name, "] ", sep="")
paired        <- opt$paired

est.frag.size <- NULL

ranges_name <- paste0(var_name, ".ranges")

message(id, "Converting BAM to ranges object:")
message(id, "Input BAM: ", bam_file)
message(id, "Object name: ", ranges_name)

if(!file.exists(bam_file)) {
	stop("Could not open BAM file: ", bam_file)
}

if(paired) {
  bam.gr <- readBAM_pe(bam_file)
} else {
  bam.gr <- readBAM(bam_file)
}

if(paired == FALSE & (ext_length == "auto" | run_filter == TRUE)) {
  message(id, "Estimating fragment length...")
  est.frag.size <- median(estimate.mean.fraglen(bam.gr, method="coverage"))
  message(id, "Fragment length estimate: ", est.frag.size)
}

if(paired == TRUE | ext_length == "native") {
	ext_length <- NULL
} else {
	if(ext_length == "auto")
	  ext_length <- est.frag.size
	else
	  ext_length <- as.integer(ext_length)
}

if(run_filter == TRUE & !paired) {
  bam.gr <- artifact_filter(bam.gr, ext_length=est.frag.size, cores=opt$cores)
}

message(id, "Extension length: ", ifelse(is.null(ext_length), "native", ext_length))

if(!is.null(ext_length)) bam.gr <- trim(resize(bam.gr, ext_length))

message(id, "Saving ranges object...")
bam.gr <- bam.gr[order(bam.gr)]
assign(ranges_name, bam.gr)
save(list=ranges_name, file=paste(ranges_name, ".RData", sep=""))
	
message(id, "Generating coverage object...")
sample.cov <- coverage(bam.gr)
rm(bam.gr)
nothing <- gc()

cov_object_name <- paste(var_name, ".cov", sep="")
assign(cov_object_name, sample.cov)
message(id, "Saving coverage object...")
save(list=cov_object_name, file=paste(cov_object_name, ".RData", sep=""))
