suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts=F, quietly=T))

option_list <- list(
		make_option(c("-f", "--file"), 
              type="character",
              help="Path to granges file"), 
      	make_option(c("-d", "--file2"), 
              type="character",
              default=NA,
              help="Path to the second granges file"),       
	  	make_option(c("-e", "--experiment"), 
	          type="character",
			  default= "dm3",
	          help="experiment genome"), 
	  	make_option(c("-s", "--spikein"), 
	  	      type="character",
	  		  default= "hg19",
	  	      help="spikein genome"), 
		make_option(c("-n", "--name"),
			  type="character",
			  help="name of the output data"))
			  
opt <- parse_args(OptionParser(option_list=option_list))

filter_gr <- function(genome, gr){
	message("filtering granges")
	chr <- grep(genome, seqlevels(gr), value=T)
	filtered_gr <- gr
	seqlevels(filtered_gr, force=T) <- chr
	seqlevels(filtered_gr, force=T) <- gsub(".*chr", "chr", seqlevels(filtered_gr)) 
	filtered_gr
}

nromalize_gr <- function(exp_gr, spikein_gr){
	message("normalization")
	exp_gr <- resize(exp_gr, 1, "start")
	exp_gr_pos <- exp_gr[strand(exp_gr) == "+"]
	exp_gr_neg <- exp_gr[strand(exp_gr) == "-"]
	
	exp_cov_pos <- coverage(exp_gr_pos) / length(spikein_gr) * 1000000
	exp_cov_neg <- coverage(exp_gr_neg) / length(spikein_gr) * 1000000 *(-1)
	cov.list <- list(pos=exp_cov_pos, neg =exp_cov_neg )
	cov.list
}


message("read granges file")
gr <- readRDS(opt$file)
exp_gr <- filter_gr(opt$experiment, gr)
spikein_gr <- filter_gr(opt$spikein, gr)
cov.list <- nromalize_gr(exp_gr, spikein_gr)

if(!is.na(opt$file2)){
    gr2 <- readRDS(opt$file2)
    exp_gr2 <- filter_gr(opt$experiment, gr2)
    spikein_gr2 <- filter_gr(opt$spikein, gr2)
    cov.list2 <- nromalize_gr(exp_gr2, spikein_gr2)
    cov.list <-list(pos = cov.list$pos + cov.list2$pos, neg=cov.list$neg + cov.list2$neg)
}

message("exporting bigwig file")
export(cov.list$pos, paste0(opt$name, "_positive.bw"))
export(cov.list$neg, paste0(opt$name, "_negative.bw"))
