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

    
    c1c2_diff_sum <- sum(c1c2_diff)
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
    
    to_return_list[["c1c2_diff_sum"]] <- c1c2_diff_sum
    to_return_list[["p1p2_diff_sum"]] <- p1p2_diff_sum
    to_return_list[["c1p1_diff_sum"]] <- c1p1_diff_sum
    to_return_list[["c1p2_diff_sum"]] <- c1p2_diff_sum
    to_return_list[["c2p1_diff_sum"]] <- c2p1_diff_sum
    to_return_list[["c2p2_diff_sum"]] <- c2p2_diff_sum

    MDR_orth1 <- max(c(c1c2_maxdiff, p1p2_maxdiff))/min(c(c1c2_maxdiff, p1p2_maxdiff))
    MDR_para1 <- max(c(c1p1_maxdiff, c2p2_maxdiff))/min(c(c1p1_maxdiff, c2p2_maxdiff))
    MDR_para2 <- max(c(c1p2_maxdiff, c2p1_maxdiff))/min(c(c1p2_maxdiff, c2p1_maxdiff))
    MDR_para3 <- max(c(c1p1_maxdiff, c1p2_maxdiff))/min(c(c1p1_maxdiff, c1p2_maxdiff))
    MDR_para4 <- max(c(c2p2_maxdiff, c2p1_maxdiff))/min(c(c2p2_maxdiff, c2p1_maxdiff))
    
    SR_orth1 <- max(c(c1c2_diff_sum, p1p2_diff_sum))/min(c(c1c2_diff_sum, p1p2_diff_sum))
    SR_para1 <- max(c(c1p1_diff_sum, c2p2_diff_sum))/min(c(c1p1_diff_sum, c2p2_diff_sum))
    SR_para2 <- max(c(c1p2_diff_sum, c2p1_diff_sum))/min(c(c1p2_diff_sum, c2p1_diff_sum))
    SR_para3 <- max(c(c1p1_diff_sum, c1p2_diff_sum))/min(c(c1p1_diff_sum, c1p2_diff_sum))
    SR_para4 <- max(c(c2p2_diff_sum, c2p1_diff_sum))/min(c(c2p2_diff_sum, c2p1_diff_sum))

    
    
    to_return_list[["MDR_orth1"]] <- MDR_orth1
    to_return_list[["MDR_para1"]] <- MDR_para1
    to_return_list[["MDR_para2"]] <- MDR_para2
    to_return_list[["MDR_para3"]] <- MDR_para3
    to_return_list[["MDR_para4"]] <- MDR_para4
    
    to_return_list[["SR_orth1"]] <- SR_orth1
    to_return_list[["SR_para1"]] <- SR_para1
    to_return_list[["SR_para2"]] <- SR_para2
    to_return_list[["SR_para3"]] <- SR_para3
    to_return_list[["SR_para4"]] <- SR_para4

    
  }
  return(to_return_list)
}







plot(log(diffs_df$sumdiffs), log(diffs_df$maxdiffs))

# of interest
"OG0000844"
"OG0000819"
"OG0000888"
test2<-orthogroups %>% filter(orthogroup == "OG0000973")
calculate_difference(test2)

diffs_df %>% filter(maxdiffs < 1.5) %>% select(orthogroup)
plot_paralogs("OG0000888")
# maxdiffs <1.5




# all ratios are low is low
diffs_df %>% filter(SR_orth1 < 1.3 & SR_para1 < 1.3 & SR_para2 < 1.3 & MDR_orth1 < 1.5 & MDR_para1 < 1.5 & MDR_para2 < 1.5) %>% select(orthogroup)



test<-orthogroups %>% filter(orthogroup == "OG0000888")
calculate_difference(test)



# at least one ortho pair is v different in ratio and paralogs are low ratio
# both species are similarly different but at least one duplicate is different
# change tis to stipulate one then the other to make sure one ortholog pair is conserved and the other not
diffs_df %>% filter(MDR_orth1 < 1.5) %>% filter(MDR_para1 > 2 | MDR_para1 > 2)



diffs_df<-data.frame(orthogroup=c(), 
                     MDR_orth1=c(), 
                     MDR_para1=c(), 
                     MDR_para2=c(),
                     MDR_para3=c(), 
                     MDR_para4=c(),
                     SR_orth1=c(), 
                     SR_para1=c(), 
                     SR_para2=c(),
                     SR_para3=c(), 
                     SR_para4=c())


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
    MDR_para3 <- con_ple_list[["MDR_para3"]]
    MDR_para4 <- con_ple_list[["MDR_para4"]]
    
    SR_orth1 <- con_ple_list[["SR_orth1"]]
    SR_para1 <- con_ple_list[["SR_para1"]]
    SR_para2 <- con_ple_list[["SR_para2"]]
    SR_para3 <- con_ple_list[["SR_para3"]]
    SR_para4 <- con_ple_list[["SR_para4"]]
    
    print(MDR_para3)
    print(MDR_para4)
    print(SR_para3)
    print(SR_para4)
    

    mydiffs<-data.frame(orthogroup=orthogroup, 
                        MDR_orth1=MDR_orth1, 
                        MDR_para1=MDR_para1, 
                        MDR_para2=MDR_para2,
                        MDR_para3=MDR_para3, 
                        MDR_para4=MDR_para4,
                        SR_orth1=SR_orth1, 
                        SR_para1=SR_para1, 
                        SR_para2=SR_para2,
                        SR_para3=SR_para3, 
                        SR_para4=SR_para4)
    

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
 

x<-diffs_df %>% select(orthogroup, MDR_orth1, MDR_para1, MDR_para2, MDR_para3, MDR_para4) %>% 
  rowwise() %>% 
  mutate(max=max(c(MDR_orth1, MDR_para1, MDR_para2, MDR_para3, MDR_para4)), 
         min=min(c(MDR_orth1, MDR_para1, MDR_para2, MDR_para3, MDR_para4))) %>% 
  mutate(MDRo1_min=min/MDR_orth1, MDR_p1_min=min/MDR_para1, MDR_p2_min=min/MDR_para2, MDR_p3_min=min/MDR_para3, MDR_p4_min=min/MDR_para4) %>%
  mutate(MDRo1_max=MDR_orth1/max, MDR_p1_max=MDR_para1/max, MDR_p2_max=MDR_para2/max, MDR_p3_max=MDR_para3/max, MDR_p4_max=MDR_para4/max) %>%
  select(-min, -max)
  
  


x %>% data.frame()
#OG0000654  3.589369  2.959965 12.919118  2.103268 18.181292 2.923557  2.416302  6.770899  1.558232 10.499425
plot_paralogs("OG0000654")

#OG0000682  3.450186  3.856872 30.261497 10.942561  1.394645 2.432495  2.460286 18.279981 10.409343  1.400984
plot_paralogs("OG0000682")

#OG0000689  3.247325 41.444500  3.156518  2.227006 58.742678 2.081218 47.051260  2.023502  1.073839 88.661654
plot_paralogs("OG0000689")

#OG0000749  5.303714  3.082321  6.074614  1.265988  2.495000 5.944943  3.434065  5.732187  1.178734  1.967558
plot_paralogs("OG0000749")

#OG0000784  3.237085  3.601510 10.141053  2.515502 14.519214 4.217085  3.383158 10.581139  2.543004 14.076923
plot_paralogs("OG0000784")

#OG0000760  1.544114  1.647049  1.535229  3.215405  1.271617 1.511038  2.039850  2.060357  3.319786  1.265991
plot_paralogs("OG0000760")

#OG0000765  1.893011  3.356271  2.177789  1.288507  1.196065 1.179970  1.380102  2.687931  1.225422  2.386671
plot_paralogs("OG0000765")

#OG0000805  1.953466  5.103871  1.640552  2.311589  1.345857 1.402595  7.704918  1.370164  3.365170  1.671047
plot_paralogs("OG0000805")

#OG0000819  1.021309  4.224744  1.083244  6.643145  1.703333 1.311555  5.600788  1.569478  8.962190  2.511426
plot_paralogs("OG0000819")



#OG0000643  4.648703  1.000861  1.164995  1.388293  1.190648
plot_paralogs("OG0000643")

1 ortholog red/blue similar = 




