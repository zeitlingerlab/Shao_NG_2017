library(GenomicRanges)
library(magrittr)

align_and_score_chipseq_peaks <- function(gr, cov) {
  gr$summit <- regionWhichMaxs(gr, chipseq.cov)
  start(gr) <-gr$summit
  end(gr) <- start(gr)
  gr$chipseq_signal <- resize(gr, width=200, fix="center") %>% regionSums(., cov)
  gr
}

align_and_score_chipnexus_peaks <- function(gr, sample, type, width=80) {
  if(type=="center"){
    gr$summit <- floor((regionWhichMaxs(gr, sample$pos)+ regionWhichMins(gr, sample$neg))/2)
  }
  if(type=="positive"){
    gr_p <- gr[strand(gr)== "+"]
    gr_n <- gr[strand(gr)== "-"]
    gr_p$summit <- regionWhichMaxs(gr_p, sample$pos)
    gr_n$summit <- regionWhichMins(gr_n, sample$neg)
    gr <- c(gr_p, gr_n)
  }
  if(type=="negative"){
    gr_p <- gr[strand(gr)== "+"]
    gr_n <- gr[strand(gr)== "-"]
    gr_p$summit <- regionWhichMins(gr_p, sample$neg)
    gr_n$summit <- regionWhichMaxs(gr_n, sample$pos)
    gr <- c(gr_p, gr_n)
  }
  start(gr) <- gr$summit
  end(gr) <- start(gr)
  gr$nexus_signal <- regionSums(resize(gr, width, "center"), sample$pos) + abs(regionSums(resize(gr, width, "center"), sample$neg))
  gr.o <-gr[order(gr$nexus_signal, decreasing=T)]
  gr.o
}