library(plyr)
library(BSgenome.Dmelanogaster.UCSC.dm3)
library(magrittr)

all_tss <- get(load("/data/analysis_code/rdata/dme_mrna_unique_tss.RData"))
promoter_table <- read.table("/data/analysis_code/promoter_elements.txt", header=T)

find_motif <- function(motif_name, fb_t_id, group_name, mismatch=0) {
    
    motif_info <- subset(promoter_table, name == motif_name)
    motif <- DNAString(motif_info$motif)
    up_dis <- motif_info$window_start
    down_dis <- motif_info$window_end
    
    gene_tss <- all_tss[all_tss$fb_t_id %in% fb_t_id]
    
    if(up_dis >= 0 & down_dis >=0){
      tss_r <- resize(gene_tss, down_dis, "start") %>%
               resize(., down_dis - up_dis, "end")
    }
    if(up_dis < 0 & down_dis >=0){
      tss_r <- resize(gene_tss, down_dis, "start") %>%
               resize(., abs(up_dis)+down_dis, "end")
    }
    if(up_dis < 0 & down_dis <0){
      tss_r <- resize(gene_tss, abs(up_dis), "end") %>%
               resize(., abs(up_dis)-abs(down_dis), "start")
    }
    
    promoter_seq <- getSeq(Dmelanogaster, tss_r)
    names(promoter_seq) <- tss_r$fb_t_id
    
    count_df <- vcountPattern(motif, promoter_seq, fixed = FALSE, min.mismatch = 0, max.mismatch = mismatch) %>%
                data.frame(fb_t_id = fb_t_id, count =.,motif_name=motif_name, group=group_name) 

    motif_enrich <- ddply(count_df, .(motif_name, group), summarize, with_motif = sum(count != 0), without_motif = sum(count == 0))
    motif_enrich
}


proportion_test <- function(values, enrichment = TRUE) {
    stopifnot(is.logical(enrichment))
    m <- matrix(values, nrow = 2, byrow = TRUE)
    alt.test <- ifelse(enrichment, "greater", "less")
    prop.test(m, alternative = alt.test)$p.value
}

single_motif_test <- function(row.df, enrichment) {
    proportion_test(as.numeric(c(row.df$s1_W, row.df$s1_WO, row.df$s2_W, row.df$s2_WO)), 
        enrichment)
}

motif_count_comparison <- function(set1, set2) {
    set1_df <- set1[, -2]
    set2_df <- set2[, -2]

    names(set1_df)[2:3] <- c("s1_W", "s1_WO")
    names(set2_df)[2:3] <- c("s2_W", "s2_WO")

    set_df <- merge(set1_df, set2_df)

    test_df <- subset(set_df, s1_W > 0 | s2_W > 0)

    e_set_test.df <- test_df
    e_set_test.df$test_type <- "enrichment"

    e_set_test.df$pvalue <- as.vector(by(e_set_test.df, 1:nrow(e_set_test.df), 
        single_motif_test, enrichment = TRUE))

    d_set_test.df <- test_df
    d_set_test.df$test_type <- "depletion"

    d_set_test.df$pvalue <- as.vector(by(d_set_test.df, 1:nrow(d_set_test.df), 
        single_motif_test, enrichment = FALSE))

    set_test.df <- rbind(e_set_test.df, d_set_test.df)

    set_notest.df <- subset(set_df, s1_W == 0 & s2_W == 0)
    if(nrow(set_notest.df) > 0) {
        set_notest.df$pvalue <- 1
        e_set_notest.df <- set_notest.df
        e_set_notest.df$test_type <- "enrichment"
        
        d_set_notest.df <- set_notest.df
        d_set_notest.df$pvalue <- 1
        d_set_notest.df$test_type <- "depletion"

        set.df <- rbind(set_test.df, e_set_notest.df, d_set_notest.df)
    } else {
        set.df <- set_test.df
    }

    set.df <- transform(set.df, enrichment = (s1_W/(s1_W + s1_WO))/(s2_W/(s2_W + s2_WO)))
    set.df <- transform(set.df, enrichment = ifelse(test_type == "enrichment", enrichment, 1/enrichment))
    set.df
    
    set.df$group <- set1$group
    set.df$enrichment <- ifelse(set.df$test_type == "enrichment", log2(set.df$enrichment), log2(set.df$enrichment))
    set.df
}


get_motif_enrichment <- function(test_tx,control_tx, group, motif, mismatch =0){
  test_motif <- find_motif(motif_name = motif,fb_t_id = test_tx,  group_name = group, mismatch = mismatch)
  control_motif <- find_motif(motif_name = motif,fb_t_id = control_tx,  group_name = group, mismatch = mismatch)
  motif_enrich_df <- motif_count_comparison(test_motif, control_motif) 
  motif_enrich_df
}

clean_motif_table <- function(results.df,pval_cutoff = 0.05){
  
  results.df$padj <- p.adjust(results.df$pvalue, method = "BH", n = nrow(results.df))
  
  results.df <- results.df[order(results.df$padj), ] %>%
                .[!duplicated(paste(.$group, .$motif_name)),]
  
  results.df$star <- with(results.df, ifelse(padj < pval_cutoff, "*", ""))
  results.df$enrichment <- ifelse(results.df$test_type == "enrichment",results.df$enrichment, -1 * results.df$enrichment)
  results.df$enrichment <- ifelse(is.infinite(results.df$enrichment), NA, results.df$enrichment)
  
  results.df$enrichment[is.na(results.df$enrichment) & results.df$padj < pval_cutoff & 
                          results.df$test_type == "depletion"] <- min(results.df$enrichment, na.rm = TRUE)
  
  results.df$enrichment[is.na(results.df$enrichment) & results.df$padj < pval_cutoff & 
                          results.df$test_type == "enrichment"] <- max(results.df$enrichment, na.rm = TRUE)

  results.df
  }




