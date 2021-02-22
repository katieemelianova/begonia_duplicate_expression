library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(magrittr)
library(wCorr)



# give full dataframe of expression and an orthogroup to plot
plot_paralogs<-function(df, og){
  test_set<-df %>% filter(orthogroup == og)
  c1<-test_set$counts[1:6]
  c2<-test_set$counts[7:12]
  p1<-test_set$counts[13:18]
  p2<-test_set$counts[19:24]
  up_lim<-max(test_set$counts) +1
  low_lim<-min(test_set$counts) -1
  plot(c1, col="red", type='l', xaxt="n", ylim = c(low_lim, up_lim), pch=16)
  axis(1, at=1:6, labels=c("fflower", "leaf", "mflower", "petiole", "root", "vegbud"))
  lines(c2, col="red", pch=16)
  lines(p1, col="blue", pch=16)
  lines(p2, col="blue", pch=16)s
}

ortho<-read_tsv("/Users/katie/Desktop/Bg/begonia_duplicate_expression/Orthogroups/Orthogroups.tsv")
colnames(ortho)<-c("orthogroup", "con", "ple")

# filter out only rows where two copies present in both CON and PLE
# separate the string of two contigs into two separate columns
x<-ortho %>%
  filter(str_count(con, "CON") == 2) %>% 
  filter(str_count(ple, "PLE") == 2) %>%
  separate(con, c("con1", "con2"), sep=", ") %>% 
  separate(ple, c("ple1", "ple2"), sep=", ")

# edit names of trinity contigs to match the ones in the counts files
x$con1<-x$con1 %>% str_remove(".p1")
x$con2<-x$con2 %>% str_remove(".p1")
x$ple1<-x$ple1 %>% str_remove(".p1")
x$ple2<-x$ple2 %>% str_remove(".p1")

con_counts<-read.table("con_avg_cpm")
ple_counts<-read.table("ple_avg_cpm")

# make the rownames an explicit variable
con_counts<-con_counts %>% tibble::rownames_to_column()
ple_counts<-ple_counts %>% tibble::rownames_to_column()

con_counts$rowname <- con_counts$rowname %>% str_remove(".p1")
ple_counts$rowname <- ple_counts$rowname %>% str_remove(".p1")


##########################################################################################
##########################################################################################



getit_realthing<-function(l){
  c1<-con_counts %>% filter(rowname %in% l$con1)
  c2<-con_counts %>% filter(rowname %in% l$con2)
  
  p1<-ple_counts %>% filter(rowname %in% l$ple1)
  p2<-ple_counts %>% filter(rowname %in% l$ple2)
  to_return_list<-list(con=c(), ple=c(), con1_name=c(), con2_name=c(), ple1_name=c(), ple2_name=c())
  
  # assert we have a complete set of 4 homologs
  if (length(c(c1$rowname, c2$rowname, p1$rowname, p2$rowname)) == 4){
    to_return_list[["con1_name"]]<-c1$rowname
    to_return_list[["con2_name"]]<-c2$rowname
    to_return_list[["ple1_name"]]<-p1$rowname
    to_return_list[["ple2_name"]]<-p2$rowname
    to_return_list[['orthogroup']]<-l$orthogroup
    
    to_return_list[["con1_expr"]]<-c1[,2:7]
    to_return_list[["con2_expr"]]<-c2[,2:7]
    to_return_list[["ple1_expr"]]<-p1[,2:7]
    to_return_list[["ple2_expr"]]<-p2[,2:7]
    
    c1<-c1[,2:7]
    c2<-c2[,2:7]
    c1<-c1 %>% as.character() %>% as.numeric()
    c2<-c2 %>% as.character() %>% as.numeric()
    p1<-p1[,2:7]
    p2<-p2[,2:7]
    p1<-p1 %>% as.character() %>% as.numeric()
    p2<-p2 %>% as.character() %>% as.numeric()
    
    # get difference between paralogs at each tissue
    con_diff<- (abs(c1 - c2))^2
    ple_diff<- (abs(p1 - p2))^2
    
    con_diff_sum <- sum(con_diff)
    ple_diff_sum <- sum(ple_diff)
    con_maxdiff <- max(con_diff)
    ple_maxdiff <- max(ple_diff)
    
    to_return_list[["con_maxdiff"]] <- con_maxdiff
    to_return_list[["ple_maxdiff"]] <- ple_maxdiff
    
    # get ratio of largest sum of differences to smallest sum of differences
    sums_ratio<-max(c(con_diff_sum, ple_diff_sum))/min(c(con_diff_sum, ple_diff_sum))
    
    # get ratio of largest maximum difference to smallest maximum difference
    max_diff_ratio<-max(c(con_maxdiff, ple_maxdiff))/min(c(con_maxdiff, ple_maxdiff))
    
    to_return_list[["max_diff_ratio"]] <- max_diff_ratio
    to_return_list[["sums_ratio"]] <- sums_ratio
    
    to_return_list[["con_sumdiffs"]] <- con_diff_sum
    to_return_list[["ple_sumdiffs"]] <- ple_diff_sum
    
  }
  return(to_return_list)
}


for (i in 1:length(rownames(x))){
  con_ple_list<-getit_realthing(x[i,])

  if (!is.null(con_ple_list[["con_sumdiffs"]]) & !is.null(con_ple_list[["ple_sumdiffs"]])) { 
    a<-cbind(con_ple_list[["con1_expr"]], con_ple_list[["con2_expr"]], con_ple_list[["ple1_expr"]], con_ple_list[["ple2_expr"]])
    a<-a %>% as.character() %>% as.numeric()
    
    cint<-con_ple_list[["con_sumdiffs"]] %>% round(0)
    pint<-con_ple_list[["ple_sumdiffs"]] %>% round(0)
    
    crange<-con_ple_list[["con_maxdiff"]] %>% round(0)
    prange<-con_ple_list[["ple_maxdiff"]] %>% round(0)
    
    max_diff_ratio<-con_ple_list[["max_diff_ratio"]] %>% round(1)
    sums_ratio<-con_ple_list[["sums_ratio"]] %>% round(1)

    ortho<-con_ple_list[["orthogroup"]]
    
    
    diff_exspr<-paste(sums_ratio, max_diff_ratio, collapse="      ")
    colours<-c(rep("red", 12), rep("blue", 12))
    up_lim<-max(a) +1
    low_lim<-min(a) -1
    #pdf(plot_name)
    plot(a[1:6], pch=8, type="l", xaxt='n', main=diff_exspr, col="red", ylim = c(low_lim, up_lim), xlab="Tissue", ylab="FPKM")
    axis(1, at=1:6, labels=c("fflower", "leaf", "mflower", "petiole", "root", "vegbud"))
    #legend("bottomright", legend=c("conchifolia", "plebeja"), col=c("red", "blue"), lty=1:2, cex=0.8)
    lines(a[7:12], col="red")
    lines(a[13:18], lty="dashed", col="blue")
    lines(a[19:24], lty="dashed", col="blue") 
    #dev.off()
    Sys.sleep(6)
    #plot_number<-plot_number+1
  }
}


plot_paralogs(total_df, "OG0000587")
plot_paralogs(total_df, "OG0000598")
plot_paralogs(total_df, "OG0000624")
plot_paralogs(total_df, "OG0000643")
plot_paralogs(total_df, "OG0000654")





