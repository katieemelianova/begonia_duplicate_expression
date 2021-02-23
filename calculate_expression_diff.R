library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(magrittr)
library(wCorr)





con_counts %>% 
  filter(rowname == as.character("CON_TRINITY_DN6298_c1_g2_i2")) %>% 
  select(CONfemaleFlower, CONleaf, CONmaleFlower, CONpetiole, CONroot, CONvegBud) %>% 
  as.character() %>%
  as.numeric()


con_counts %>% filter(rowname == "CON_TRINITY_DN6687_c3_g3_i1") %>% select(CONfemaleFlower, CONleaf, CONmaleFlower, CONpetiole, CONroot, CONvegBud)


ortho<-read_tsv("/Users/katie/Desktop/Bg/begonia_duplicate_expression/Orthogroups/Orthogroups.tsv")
colnames(ortho)<-c("orthogroup", "con", "ple")

# filter out only rows where two copies present in both CON and PLE
# separate the string of two contigs into two separate columns
orthogroups<-ortho %>%
  filter(str_count(con, "CON") == 2) %>% 
  filter(str_count(ple, "PLE") == 2) %>%
  separate(con, c("con1", "con2"), sep=", ") %>% 
  separate(ple, c("ple1", "ple2"), sep=", ")

# edit names of trinity contigs to match the ones in the counts files
orthogroups$con1<-orthogroups$con1 %>% str_remove(".p1")
orthogroups$con2<-orthogroups$con2 %>% str_remove(".p1")
orthogroups$ple1<-orthogroups$ple1 %>% str_remove(".p1")
orthogroups$ple2<-orthogroups$ple2 %>% str_remove(".p1")

con_counts<-read.table("con_avg_cpm") %>% tibble::rownames_to_column()
ple_counts<-read.table("ple_avg_cpm") %>% tibble::rownames_to_column()

con_counts$rowname <- con_counts$rowname %>% str_remove(".p1")
ple_counts$rowname <- ple_counts$rowname %>% str_remove(".p1")


plot_paralogs<-function(og){
  con1<-orthogroups %>% filter(orthogroup == og) %>% select(con1)
  con2<-orthogroups %>% filter(orthogroup == og) %>% select(con2)
  ple1<-orthogroups %>% filter(orthogroup == og) %>% select(ple1)
  ple2<-orthogroups %>% filter(orthogroup == og) %>% select(ple2)
  
  con1_ex <- con_counts %>% filter(rowname == as.character(con1)) %>% select(CONfemaleFlower, CONleaf, CONmaleFlower, CONpetiole, CONroot, CONvegBud) %>% as.character() %>% as.numeric() 
  con2_ex <- con_counts %>% filter(rowname == as.character(con2)) %>% select(CONfemaleFlower, CONleaf, CONmaleFlower, CONpetiole, CONroot, CONvegBud) %>% as.character() %>% as.numeric() 
  ple1_ex <- ple_counts %>% filter(rowname == as.character(ple1)) %>% select(PLEfemaleFlower, PLEleaf, PLEmaleFlower, PLEpetiole, PLEroot, PLEvegBud) %>% as.character() %>% as.numeric() 
  ple2_ex <- ple_counts %>% filter(rowname == as.character(ple2)) %>% select(PLEfemaleFlower, PLEleaf, PLEmaleFlower, PLEpetiole, PLEroot, PLEvegBud) %>% as.character() %>% as.numeric() 
  
  up_lim<-max(c(con1_ex, con2_ex, ple1_ex, ple2_ex)) +1
  low_lim<-min(c(con1_ex, con2_ex, ple1_ex, ple2_ex)) -1
  
  plot(con1_ex, col="red", type='l', xaxt="n", ylim = c(low_lim, up_lim), pch=16, xlab="Tissue", ylab="FPKM")
  axis(1, at=1:6, labels=c("fflower", "leaf", "mflower", "petiole", "root", "vegbud"))
  lines(ple1_ex, col="red", pch=16)
  lines(con2_ex, col="blue", pch=16)
  lines(ple2_ex, col="blue", pch=16)
}


##########################################################################################
##########################################################################################



calculate_difference<-function(l){
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
    
    c1<-c1[,2:7] %>% as.character() %>% as.numeric()
    c2<-c2[,2:7] %>% as.character() %>% as.numeric()
    p1<-p1[,2:7] %>% as.character() %>% as.numeric()
    p2<-p2[,2:7] %>% as.character() %>% as.numeric()

    # get difference between paralogs at each tissue
    con_diff<- (abs(c1 - c2))^2
    ple_diff<- (abs(p1 - p2))^2
    
    con_diff_sum <- sum(con_diff)
    ple_diff_sum <- sum(ple_diff)
    con_maxdiff <- max(con_diff)
    ple_maxdiff <- max(ple_diff)

    # get ratio of largest sum of differences to smallest sum of differences
    sums_ratio<-max(c(con_diff_sum, ple_diff_sum))/min(c(con_diff_sum, ple_diff_sum))
    
    # get ratio of largest maximum difference to smallest maximum difference
    max_diff_ratio<-max(c(con_maxdiff, ple_maxdiff))/min(c(con_maxdiff, ple_maxdiff))
    
    to_return_list[["max_diff_ratio"]] <- max_diff_ratio
    to_return_list[["sums_ratio"]] <- sums_ratio
  }
  return(to_return_list)
}



diffs_df<-data.frame(orthogroup=c(), sumdiffs=c(), maxdiffs<-c())
plot(log(diffs_df$sumdiffs), log(diffs_df$maxdiffs))




# orthogroups whete maxdiff is > 3
plot_paralogs("OG0000529")
plot_paralogs("OG0000553")
plot_paralogs("OG0000587")
plot_paralogs("OG0000598")
plot_paralogs("OG0000605")
plot_paralogs("OG0000606")
plot_paralogs("OG0000619")
plot_paralogs("OG0000622")
plot_paralogs("OG0000624")
plot_paralogs("OG0000625")
plot_paralogs("OG0000643")
plot_paralogs("OG0000651")
plot_paralogs("OG0000654")
plot_paralogs("OG0000666")
plot_paralogs("OG0000672")
plot_paralogs("OG0000682")
plot_paralogs("OG0000685")
plot_paralogs("OG0000689")
plot_paralogs("OG0000697")
plot_paralogs("OG0000703")
plot_paralogs("OG0000712")
plot_paralogs("OG0000733")
plot_paralogs("OG0000737")
plot_paralogs("OG0000747")
plot_paralogs("OG0000749")
plot_paralogs("OG0000753")
plot_paralogs("OG0000765")
plot_paralogs("OG0000769")
plot_paralogs("OG0000784")
plot_paralogs("OG0000799")
plot_paralogs("OG0000805")
plot_paralogs("OG0000807")
plot_paralogs("OG0000837")
plot_paralogs("OG0000861")
plot_paralogs("OG0000862")
plot_paralogs("OG0000887")
plot_paralogs("OG0000908")
plot_paralogs("OG0000915")
plot_paralogs("OG0000933")
plot_paralogs("OG0000958")
plot_paralogs("OG0000961")
plot_paralogs("OG0000963")
plot_paralogs("OG0000966")
plot_paralogs("OG0000998")
plot_paralogs("OG0001005")
plot_paralogs("OG0001014")
plot_paralogs("OG0001015")
plot_paralogs("OG0001018")
plot_paralogs("OG0001030")
plot_paralogs("OG0001035")
plot_paralogs("OG0001045")
plot_paralogs("OG0001048")
plot_paralogs("OG0001052")
plot_paralogs("OG0001053")
plot_paralogs("OG0001060")
plot_paralogs("OG0001071")
plot_paralogs("OG0001074")
plot_paralogs("OG0001087")
plot_paralogs("OG0001091")
plot_paralogs("OG0001095")
plot_paralogs("OG0001098")
plot_paralogs("OG0001105")
plot_paralogs("OG0001106")
plot_paralogs("OG0001115")
plot_paralogs("OG0001119")
plot_paralogs("OG0001121")
plot_paralogs("OG0001137")
plot_paralogs("OG0001141")
plot_paralogs("OG0001143")
plot_paralogs("OG0001148")
plot_paralogs("OG0001156")
plot_paralogs("OG0001171")
plot_paralogs("OG0001179")
plot_paralogs("OG0001184")
plot_paralogs("OG0001186")
plot_paralogs("OG0001196")
plot_paralogs("OG0001200")
plot_paralogs("OG0001202")
plot_paralogs("OG0001205")
plot_paralogs("OG0001209")
plot_paralogs("OG0001220")
plot_paralogs("OG0001229")
plot_paralogs("OG0001247")
plot_paralogs("OG0001272")
plot_paralogs("OG0001280")
plot_paralogs("OG0001294")
plot_paralogs("OG0001297")
plot_paralogs("OG0001301")
plot_paralogs("OG0001304")
plot_paralogs("OG0001308")
plot_paralogs("OG0001309")
plot_paralogs("OG0001337")
plot_paralogs("OG0001338")
plot_paralogs("OG0001340")
plot_paralogs("OG0001355")


diffs_df %>% filter(maxdiffs > 3) %>% select(orthogroup)

for (i in 1:length(rownames(orthogroups))){
  con_ple_list<-calculate_difference(orthogroups[i,])

  if (!is.null(con_ple_list[["max_diff_ratio"]]) & !is.null(con_ple_list[["sums_ratio"]])) { 
    
    con1_exp<-con_ple_list[["con1_expr"]] %>% as.character() %>% as.numeric()
    con2_exp<-con_ple_list[["con2_expr"]] %>% as.character() %>% as.numeric()
    ple1_exp<-con_ple_list[["ple1_expr"]] %>% as.character() %>% as.numeric()
    ple2_exp<-con_ple_list[["ple2_expr"]] %>% as.character() %>% as.numeric()
    orthogroup<-con_ple_list[["orthogroup"]]                      
    max_diff_ratio<-con_ple_list[["max_diff_ratio"]] %>% round(1)
    sums_ratio<-con_ple_list[["sums_ratio"]] %>% round(1)

    mydiffs<-data.frame(orthogroup=orthogroup, sumdiffs=sums_ratio, maxdiffs=max_diff_ratio)
    diffs_df<-rbind(diffs_df, mydiffs)

    #diff_exspr<-paste(sums_ratio, max_diff_ratio, collapse="      ")
    #up_lim<-max(c(con1_exp, con2_exp, ple1_exp, ple2_exp)) +1
    #low_lim<-min(c(con1_exp, con2_exp, ple1_exp, ple2_exp)) -1
    #plot(con1_exp, pch=8, type="l", xaxt='n', main=diff_exspr, col="red", ylim = c(low_lim, up_lim), xlab="Tissue", ylab="FPKM")
    #axis(1, at=1:6, labels=c("fflower", "leaf", "mflower", "petiole", "root", "vegbud"))
    #lines(con2_exp, col="red")
    #lines(ple1_exp, lty="dashed", col="blue")
    #lines(ple2_exp, lty="dashed", col="blue") 
    #Sys.sleep(3)
  }
}
 




