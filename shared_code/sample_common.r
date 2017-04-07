sample_list <- read.csv("./sample_list.csv")

load_bigwig <- function(name, sample_format = "default", gr=NULL){
  if(!(name %in% sample_list$sample_name) &!(name %in% sample_list$short_name)){
    message("No sample found, please use the sample_list.csv spreadsheet to get the correct sample_name or short_name")
  }else{
    if(sample_format == "default"){
      if(name %in% sample_list$sample_name){
        sample.path <- list(pos =  unique(as.character(subset(sample_list, sample_name == name)$bigwig_positive)),
                            neg =  unique(as.character(subset(sample_list, sample_name ==name)$bigwig_negative)))
      }
      if(name %in% sample_list$short_name){
        sample.path <- list(pos =  unique(as.character(subset(sample_list, short_name == name)$bigwig_positive)),
                            neg =  unique(as.character(subset(sample_list, short_name ==name)$bigwig_negative)))
      }
      if(nchar(sample.path$neg) < 2){
        sample.path <-sample.path$pos
      }
    }
    
    if(sample_format == "separate"){
      if(name %in% sample_list$short_name){
        sample.path1 <- list(pos =  unique(as.character(subset(sample_list, short_name == name)$bigwig_positive))[1],
                             neg =  unique(as.character(subset(sample_list, short_name ==name)$bigwig_negative))[1])
        sample.path2 <- list(pos =  unique(as.character(subset(sample_list, short_name == name)$bigwig_positive))[2],
                             neg =  unique(as.character(subset(sample_list, short_name ==name)$bigwig_negative))[2])
        sample.path <- list(sample.path1, sample.path2)
      }else{
        message("Separate loading can only be used when sample short name is provided")
      }
    }
    
    if(sample_format == "data"){
      if(name %in% sample_list$sample_name){
        sample.path <- list(pos =  unique(as.character(subset(sample_list, sample_name == name)$bigwig_positive)),
                            neg =  unique(as.character(subset(sample_list, sample_name ==name)$bigwig_negative)))
        sample.path <- list(pos=check_coverage_argument(sample.path$pos, regions=gr), 
                            neg=check_coverage_argument(sample.path$neg, regions=gr))
      }
      if(name %in% sample_list$short_name){
        if(length(grep("spikein", name)) == 0 & length(grep("chipseq", name)) == 0){
          sample.path <- list(pos =  unique(as.character(subset(sample_list, short_name == name)$bigwig_positive)),
                              neg =  unique(as.character(subset(sample_list, short_name ==name)$bigwig_negative)))
          sample.path <- list(pos=check_coverage_argument(sample.path$pos, regions=gr), 
                              neg=check_coverage_argument(sample.path$neg, regions=gr))
        }
        if(length(grep("spikein", name)) == 1){
          sample.path1 <- list(pos =  unique(as.character(subset(sample_list, short_name == name)$bigwig_positive))[1],
                              neg =  unique(as.character(subset(sample_list, short_name ==name)$bigwig_negative))[1])
          sample.path2 <- list(pos =  unique(as.character(subset(sample_list, short_name == name)$bigwig_positive))[2],
                               neg =  unique(as.character(subset(sample_list, short_name ==name)$bigwig_negative))[2])
          sample.path <- list(pos=check_coverage_argument(sample.path1$pos, regions=gr) +check_coverage_argument(sample.path2$pos, regions=gr), 
                              neg=check_coverage_argument(sample.path1$neg, regions=gr) + check_coverage_argument(sample.path2$neg, regions=gr))
        }
      }
    }
    
    sample.path
  }
}

