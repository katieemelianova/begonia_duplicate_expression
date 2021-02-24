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
  
  
  print(con1)
  print(con2)
  print(ple1)
  print(ple2)
  print(con1_ex)
  print(con2_ex)
  print(ple1_ex)
  print(ple2_ex)
  
  up_lim<-max(c(con1_ex, con2_ex, ple1_ex, ple2_ex)) +1
  low_lim<-min(c(con1_ex, con2_ex, ple1_ex, ple2_ex)) -1
  
  print(up_lim)
  print(low_lim)
  
  plot(con1_ex, col="red", type='l', xaxt="n", ylim = c(low_lim, up_lim), pch=16, xlab="Tissue", ylab="FPKM")
  axis(1, at=1:6, labels=c("fflower", "leaf", "mflower", "petiole", "root", "vegbud"))
  lines(con2_ex, col="red", pch=16)
  lines(ple1_ex, col="blue", pch=16)
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

    
    c1c2_diff<- (abs(c1 - c2))
    p1p2_diff<- (abs(p1 - p2))
    c1p1_diff<- (abs(c1 - p1))
    c1p2_diff<- (abs(c1 - p2))
    c2p1_diff<- (abs(c2 - p1))
    c2p2_diff<- (abs(c2 - p2))
    
    to_return_list[["c1c2_diffs"]] <- c1c2_diff
    to_return_list[["p1p2_diffs"]] <- p1p2_diff
    to_return_list[["c1p1_diffs"]] <- c1p1_diff
    to_return_list[["c1p2_diffs"]] <- c1p2_diff
    to_return_list[["c2p1_diffs"]] <- c2p1_diff
    to_return_list[["c2p2_diffs"]] <- c2p2_diff
    
    c1c2_diff_sum<- sum(c1c2_diff)
    p1p2_diff_sum <- sum(p1p2_diff)
    c1p1_diff_sum <- sum(c1p1_diff)
    c1p2_diff_sum <- sum(c1p2_diff)
    c2p1_diff_sum <- sum(c2p1_diff)
    c2p2_diff_sum <- sum(c2p2_diff)
    
    c1c2_maxdiff <- max(c1c2_diff)
    p1p2_maxdiff <- max(p1p2_diff)
    c1p1_maxdiff <- max(c1p1_diff)
    c1p2_maxdiff <- max(c1p2_diff)
    c2p1_maxdiff <- max(c2p1_diff)
    c2p2_maxdiff <- max(c2p2_diff)

    MDR_orth1 <- max(c(c1c2_maxdiff, p1p2_maxdiff))/min(c(c1c2_maxdiff, p1p2_maxdiff))
    MDR_para1 <- max(c(c1p1_maxdiff, c2p2_maxdiff))/min(c(c1p1_maxdiff, c2p2_maxdiff))
    MDR_para2 <- max(c(c1p2_maxdiff, c2p1_maxdiff))/min(c(c1p2_maxdiff, c2p1_maxdiff))
    
    SR_orth1 <- max(c(c1c2_diff_sum, p1p2_diff_sum))/min(c(c1c2_diff_sum, p1p2_diff_sum))
    SR_para1 <- max(c(c1p1_diff_sum, c2p2_diff_sum))/min(c(c1p1_diff_sum, c2p2_diff_sum))
    SR_para2 <- max(c(c1p2_diff_sum, c2p1_diff_sum))/min(c(c1p2_diff_sum, c2p1_diff_sum))
    
    to_return_list[["MDR_orth1"]] <- MDR_orth1
    to_return_list[["MDR_para1"]] <- MDR_para1
    to_return_list[["MDR_para2"]] <- MDR_para2
    
    to_return_list[["SR_orth1"]] <- SR_orth1
    to_return_list[["SR_para1"]] <- SR_para1
    to_return_list[["SR_para2"]] <- SR_para2

  }
  return(to_return_list)
}







plot(log(diffs_df$sumdiffs), log(diffs_df$maxdiffs))

# of interest
"OG0000844"
"OG0000819"
"OG0000888"
test2<-orthogroups %>% filter(orthogroup == "OG0000888")
calculate_difference(test2)

diffs_df %>% filter(maxdiffs < 1.5) %>% select(orthogroup)
plot_paralogs("OG0000888")
# maxdiffs <1.5


diffs_df<-data.frame(orthogroup=c(), 
                     MDR_orth1=c(), 
                     MDR_para1=c(), 
                     MDR_para2=c(), 
                     SR_orth1=c(), 
                     SR_para1=c(), 
                     SR_para2=c())

# should all be pretty conserved
diffs_df %>% filter(MDR_orth1 < 1.5 & MDR_para1 < 1.5 & MDR_para2 < 1.5)
plot_paralogs("OG0001149")
plot_paralogs("OG0000771")
plot_paralogs("OG0001169")



# at least one ortho pair is v different in ratio and paralogs are low ratio
# both species are similarly different but at least one duplicate is different
# change tis to stipulate one then the other to make sure one ortholog pair is conserved and the other not
diffs_df %>% filter(MDR_orth1 < 1.5) %>% filter(MDR_para1 > 2 | MDR_para1 > 2)


# OG0001071  6.121649  2.902226  3.528464 2.926645  2.572447  2.604369
plot_paralogs("OG0001071")

# OG0000606  1.854607  1.854728  4.678835 1.964908  2.021434  6.197488
# where two orthologs are similar and two orthologs are not
plot_paralogs("OG0000606")

# OG0000682  3.450186  3.856872 30.261497 2.432495  2.460286 18.279981
# DITTO where two orthologs are similar and two orthologs are not
plot_paralogs("OG0000682")

# OG0000963  1.885552  1.973281 31.012754 1.822713  1.913121 19.819666
# DITTO
plot_paralogs("OG0000963")


# OG0001083  1.196543  1.112929  1.135611 1.312299  1.087589  1.036780
# similar - all less than 1.5ish
plot_paralogs("OG0001083")



for (i in 1:length(rownames(orthogroups))){
  con_ple_list<-calculate_difference(orthogroups[i,])

  if (!is.null(con_ple_list[["MDR_orth1"]])) { 
    
    con1_exp<-con_ple_list[["con1_expr"]] %>% as.character() %>% as.numeric()
    con2_exp<-con_ple_list[["con2_expr"]] %>% as.character() %>% as.numeric()
    ple1_exp<-con_ple_list[["ple1_expr"]] %>% as.character() %>% as.numeric()
    ple2_exp<-con_ple_list[["ple2_expr"]] %>% as.character() %>% as.numeric()
    orthogroup<-con_ple_list[["orthogroup"]]                      
    

    MDR_orth1 <- con_ple_list[["MDR_orth1"]]
    MDR_para1 <- con_ple_list[["MDR_para1"]]
    MDR_para2 <- con_ple_list[["MDR_para2"]]
    
    SR_orth1 <- con_ple_list[["SR_orth1"]]
    SR_para1 <- con_ple_list[["SR_para1"]]
    SR_para2 <- con_ple_list[["SR_para2"]]
    

    mydiffs<-data.frame(orthogroup=orthogroup, 
                        MDR_orth1=MDR_orth1, 
                        MDR_para1=MDR_para1, 
                        MDR_para2=MDR_para2, 
                        SR_orth1=SR_orth1, 
                        SR_para1=SR_para1, 
                        SR_para2=SR_para2)
    

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
 




