suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts=F, quietly=T))

option_list <- list(
		make_option(c("-f", "--file"),
              type="character",
              help="Path to granges file"),
	  	make_option(c("-s", "--file2"),
	          type="character",
			  default= NA,
	          help="Path to the second granges file"),
		make_option(c("-n", "--name"),
			  type="character",
			  help="name of the output data"),
	  	make_option(c("-m", "--normalization"),
			  default = "T",
	  		  type="character",
	  		  help="parameter to indicate if normalize to the total reads, default is true"),
		make_option(c("-t", "--type"),
			  default = "chipnexus",
			  type="character",
			  help="experiment type, chipseq or chipnexus"))

opt <- parse_args(OptionParser(option_list=option_list))


total_reads_normalization <- function(gr, type){
	if(type == "chipnexus"){
		total_reads <- length(gr)
		gr.pos <- gr[strand(gr) == "+"]
		gr.neg <- gr[strand(gr) == "-"]
		gr.pos.cov <- coverage(gr.pos) / total_reads * 1000000
		gr.neg.cov <- coverage(gr.neg) / total_reads * 1000000 *(-1)
		cov.list <- list(pos =gr.pos.cov, neg = gr.neg.cov )
	}
	if(type == "chipseq"){
		total_reads <- length(gr)
		cov <- coverage(gr) / total_reads *1000000
		cov.list <- list(pos = cov)
	}
	cov.list
}

if(opt$type == "chipnexus"){
	message("sample being processed is chipnexus")
	message("reading the granges file")
	gr  <- readRDS(opt$file)
	gr <- resize(gr, 1, "start")

	if(opt$normalization == "T"){
		message("generating the normalized coverage files")
		cov.list <- total_reads_normalization(gr, "chipnexus")
	}else{
		message("generating the non-normalized coverage files")
		gr.pos <- gr[strand(gr) == "+"]
		gr.neg <- gr[strand(gr) == "-"]
		gr.pos.cov <- coverage(gr.pos)
		gr.neg.cov <- coverage(gr.neg) * (-1)
		cov.list <- list(pos =gr.pos.cov, neg = gr.neg.cov )
	}

	if(!is.na(opt$file2)){
		message("reading the second granges file")
		gr2  <- readRDS(opt$file2)
		gr2 <- resize(gr2, 1, "start")

		if(opt$normalization == "T"){
			message("generating the second normalized coverage files")
			cov.list2 <- total_reads_normalization(gr2, "chipnexus")
			cov.list <- list(pos =cov.list$pos + cov.list2$pos, neg = cov.list$neg + cov.list2$neg )
		}else{
			message("generating the second non-normalized coverage files")
			gr.pos2 <- gr[strand(gr2) == "+"]
			gr.neg2 <- gr[strand(gr2) == "-"]
			gr.pos.cov2 <- coverage(gr.pos2)
			gr.neg.cov2 <- coverage(gr.neg2) * (-1)
			cov.list2 <- list(pos =gr.pos.cov2, neg = gr.neg.cov2 )
			cov.list <- list(pos =cov.list$pos + cov.list2$pos, neg = cov.list$neg + cov.list2$neg )
		}

	}
	message("exporting bigwig files")
	export(cov.list$pos, paste0(opt$name, "_positive.bw"))
	export(cov.list$neg, paste0(opt$name, "_negative.bw"))

}

if(opt$type == "chipseq"){
	message("sample being processed is chipseq")
	message("reading the granges file")
	gr <- get(load(opt$file))
	if(opt$normalization == "T"){
		message("generating the normalized coverage files")
		cov.list <- total_reads_normalization(gr, "chipseq")
	}else{
		message("generating the non-normalized coverage files")
		cov.list <- list(pos=coverage(gr))
	}

	if(!is.na(opt$file2)){
		message("reading the second granges")
		gr2 <- get(load(opt$file2))
		if(opt$normalization == "T"){
			message("generating the normalized coverage files")
			cov.list2 <- total_reads_normalization(gr2, "chipseq")
			cov.list <- list(pos =cov.list$pos + cov.list2$pos)
		}else{
			message("generating the non-normalized coverage files")
			cov.list2 <- list(pos=coverage(gr2))
			cov.list <- list(pos =cov.list$pos + cov.list2$pos)
		}
	}
	message("exporting bigwig files")
	export(cov.list$pos, paste0(opt$name, ".bw"))
}
