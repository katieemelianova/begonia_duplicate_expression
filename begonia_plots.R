library(ggplot2)
library(reshape)
library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)
library(UpSetR)
library(ggpubr)
library(pheatmap)
library(edgeR)
library(GO.db)
library(seqinr)
library(readr)

#########################################################
###           prepare normalised counts                ##
###     Do not filter by cpm to include low counts     ##
#########################################################

#read in featurecounts output
con_counts<-read.table("/Users/katie/Desktop/Bg/begonia_duplicate_expression/featurecount/con/con_counts", header=TRUE)
ple_counts<-read.table("/Users/katie/Desktop/Bg/begonia_duplicate_expression/featurecount/ple/ple_counts", header=TRUE)

# prepend them with CON_ and PLE_ to match that of Orthogroups
con_names<-paste("CON", con_counts$Geneid, sep="_")
ple_names<-paste("PLE", ple_counts$Geneid, sep="_")

# remove all the path names from the column names to make them shorter
con_cols<-colnames(con_counts) %>% 
  str_remove("X.home.kemelianova.to_download.trimmed_decontaminated.star.") %>% 
  str_remove("Aligned.sortedByCoord.out.bam")

ple_cols<-colnames(ple_counts) %>% 
  str_remove("X.home.kemelianova.to_download.trimmed_decontaminated.star.") %>% 
  str_remove("Aligned.sortedByCoord.out.bam")

# rename columns using the shorter version
colnames(con_counts)<-con_cols
colnames(ple_counts)<-ple_cols

con_lengths<-con_counts %>% dplyr::select(Length) %>% set_rownames(con_names)
ple_lengths<-ple_counts %>% dplyr::select(Length) %>% set_rownames(ple_names)

# select columns of the counts, excluding chr start end etc info at beginning
con_counts<-con_counts[,7:23]
row.names(con_counts)<-con_names
ple_counts<-ple_counts[,7:23]
row.names(ple_counts)<- ple_names

# assign groups to the columns in con and ple
# missing letters are for Con root 2 and Ple petiole 3
con_groups<-c('A','A','A','B','B','B','C','C','C','D','D','D','E','E','F','F','F')
ple_groups<-c('G','G','G','H','H','H','I','I','I','J','J','K','K','K','L','L','L')

# create DGE list
# keep only rows where counts per million is above 100 in at least 2 samples
# calculate normalisation factors
# estimate tagwise and common dispersal
con_d<-DGEList(counts=con_counts,group=factor(con_groups))
con_d_lengths<-con_lengths %>% data.frame()
con_d$samples$lib.size<-colSums(con_d$counts)
con_d<-calcNormFactors(con_d, method="TMM")

ple_d<-DGEList(counts= ple_counts,group=factor(ple_groups))
ple_d_lengths<-ple_lengths %>% data.frame()
rownames(ple_d_lengths)<-rownames(ple_d$counts)
ple_d$samples$lib.size<-colSums(ple_d$counts)
ple_d<-calcNormFactors(ple_d, method="TMM")



#########################################################
###           get average FPKM per tissue              ##
#########################################################


CONfemaleFlower<-con_d$counts %>% data.frame() %>% dplyr::select(CONfemaleFlower1, CONfemaleFlower2, CONfemaleFlower3) %>% rpkmByGroup(gene.length=con_d_lengths)  %>% data.frame()
CONleaf<-con_d$counts %>% data.frame() %>% dplyr::select(CONleaf1, CONleaf2, CONleaf3) %>% rpkmByGroup(gene.length=con_d_lengths) %>% data.frame()
CONmaleFlower<-con_d$counts %>% data.frame() %>% dplyr::select(CONmaleFlower1, CONmaleFlower2, CONmaleFlower3) %>% rpkmByGroup(gene.length=con_d_lengths) %>% data.frame()
CONpetiole<-con_d$counts %>% data.frame() %>% dplyr::select(CONpetiole1, CONpetiole2, CONpetiole3) %>% rpkmByGroup(gene.length=con_d_lengths) %>% data.frame()
CONroot<-con_d$counts %>% data.frame() %>% dplyr::select(CONroot1, CONroot3, CONvegBud1) %>% rpkmByGroup(gene.length=con_d_lengths) %>% data.frame()
CONvegBud<-con_d$counts %>% data.frame() %>% dplyr::select(CONvegBud1, CONvegBud2, CONvegBud3) %>% rpkmByGroup(gene.length=con_d_lengths) %>% data.frame()

PLEfemaleFlower<-ple_d$counts %>% data.frame() %>% dplyr::select(PLEfemaleFlower1, PLEfemaleFlower2, PLEfemaleFlower3) %>% rpkmByGroup(gene.length=ple_d_lengths)  %>% data.frame()
PLEleaf<-ple_d$counts %>% data.frame() %>% dplyr::select(PLEleaf1, PLEleaf2, PLEleaf3) %>% rpkmByGroup(gene.length=ple_d_lengths)  %>% data.frame()
PLEmaleFlower<-ple_d$counts %>% data.frame() %>% dplyr::select(PLEmaleFlower1, PLEmaleFlower2, PLEmaleFlower3) %>% rpkmByGroup(gene.length=ple_d_lengths)  %>% data.frame()
PLEpetiole<-ple_d$counts %>% data.frame() %>% dplyr::select(PLEpetiole1, PLEpetiole2) %>% rpkmByGroup(gene.length=ple_d_lengths)  %>% data.frame()
PLEroot<-ple_d$counts %>% data.frame() %>% dplyr::select(PLEroot1, PLEroot2, PLEroot3) %>% rpkmByGroup(gene.length=ple_d_lengths) %>% data.frame()
PLEvegBud<-ple_d$counts %>% data.frame() %>% dplyr::select(PLEvegBud1, PLEvegBud2, PLEvegBud3) %>% rpkmByGroup(gene.length=ple_d_lengths)  %>% data.frame()

# bind con and ple together per species
con_avg_cpm<-cbind(CONfemaleFlower, CONleaf, CONmaleFlower, CONpetiole, CONroot, CONvegBud)
ple_avg_cpm<-cbind(PLEfemaleFlower, PLEleaf, PLEmaleFlower, PLEpetiole, PLEroot, PLEvegBud)

# transfer row names
rownames(con_avg_cpm)<-rownames(con_d$counts)
rownames(ple_avg_cpm)<-rownames(ple_d$counts)

# assign colnames
colnames(con_avg_cpm)<-c("CONfemaleFlower", "CONleaf", "CONmaleFlower", "CONpetiole", "CONroot", "CONvegBud")
colnames(ple_avg_cpm)<-c("PLEfemaleFlower", "PLEleaf", "PLEmaleFlower", "PLEpetiole", "PLEroot", "PLEvegBud")


############################################################################
###     filter features which are expressed > 1 FPKM in *one* tissue       #
###    all other tissues are expressed at < 1 (i.e. unique to one tissue)  #
############################################################################

CONfemaleFlower_specific<-con_avg_cpm %>% filter(CONfemaleFlower > 1 & CONroot < 1 & CONpetiole < 1 & CONvegBud < 1 & CONmaleFlower < 1 & CONleaf < 1)
CONmaleFlower_specific<-con_avg_cpm %>% filter(CONmaleFlower > 1 & CONroot < 1 & CONpetiole < 1 & CONvegBud < 1 & CONleaf < 1 & CONfemaleFlower < 1)
CONleaf_specific<-con_avg_cpm %>% filter(CONleaf > 1 & CONmaleFlower < 1 & CONroot < 1 & CONpetiole < 1 & CONvegBud < 1 & CONfemaleFlower < 1)
CONpetiole_specific<-con_avg_cpm %>% filter(CONpetiole > 1 & CONleaf < 1 & CONmaleFlower < 1 & CONroot < 1 & CONvegBud < 1 & CONfemaleFlower < 1)
CONroot_specific<-con_avg_cpm %>% filter(CONroot > 1 & CONpetiole < 1 & CONleaf < 1 & CONmaleFlower < 1 & CONvegBud < 1 & CONfemaleFlower < 1)
CONvegBud_specific<-con_avg_cpm %>% filter(CONvegBud > 1 & CONroot < 1 & CONpetiole < 1 & CONleaf < 1 & CONmaleFlower < 1 & CONfemaleFlower < 1)

PLEfemaleFlower_specific<-ple_avg_cpm %>% filter(PLEfemaleFlower > 1 & PLEroot < 1 & PLEpetiole < 1 & PLEvegBud < 1 & PLEmaleFlower < 1 & PLEleaf < 1)
PLEmaleFlower_specific<-ple_avg_cpm %>% filter(PLEmaleFlower > 1 & PLEroot < 1 & PLEpetiole < 1 & PLEvegBud < 1 & PLEleaf < 1 & PLEfemaleFlower < 1)
PLEleaf_specific<-ple_avg_cpm %>% filter(PLEleaf > 1 & PLEmaleFlower < 1 & PLEroot < 1 & PLEpetiole < 1 & PLEvegBud < 1 & PLEfemaleFlower < 1)
PLEpetiole_specific<-ple_avg_cpm %>% filter(PLEpetiole > 1 & PLEleaf < 1 & PLEmaleFlower < 1 & PLEroot < 1 & PLEvegBud < 1 & PLEfemaleFlower < 1)
PLEroot_specific<-ple_avg_cpm %>% filter(PLEroot > 1 & PLEpetiole < 1 & PLEleaf < 1 & PLEmaleFlower < 1 & PLEvegBud < 1 & PLEfemaleFlower < 1)
PLEvegBud_specific<-ple_avg_cpm %>% filter(PLEvegBud > 1 & PLEroot < 1 & PLEpetiole < 1 & PLEleaf < 1 & PLEmaleFlower < 1 & PLEfemaleFlower < 1)
#Removeing this one because it has a length one 1 and the GO term name is v long
PLEvegBud_specific_goterms<-PLEvegBud_specific_goterms[!PLEvegBud_specific_goterms %in% c('GO:0000447')]


##############################################################################
###                  Loop through uniquely expressed genes                   #
###  use to find GO terms associated with that gene (from trinotate output)  #
##     append these GO terms to GO terms foer that tissue specific GO terms  #
##############################################################################

con_annotation<-read.table("con_go_annotations_trans.txt")
ple_annotation<-read.table("ple_go_annotations_trans.txt")

CONfemaleFlower_specific_goterms<-c()
CONmaleFlower_specific_goterms<-c()
CONleaf_specific_goterms<-c()
CONpetiole_specific_goterms<-c()
CONroot_specific_goterms<-c()
CONvegBud_specific_goterms<-c()

PLEfemaleFlower_specific_goterms<-c()
PLEmaleFlower_specific_goterms<-c()
PLEleaf_specific_goterms<-c()
PLEpetiole_specific_goterms<-c()
PLEroot_specific_goterms<-c()
PLEvegBud_specific_goterms<-c()


for (i in rownames(CONfemaleFlower_specific)){
  seqname<-str_split(i, "[.]")[[1]][1]
  seqname<-str_remove(seqname, "CON_")
  line<-con_annotation[con_annotation$V1 == seqname,]
  go_terms<-str_split(line$V2, ",")
  if(length(go_terms) != 0){
    go_terms = go_terms[[1]]
    CONfemaleFlower_specific_goterms <-c(CONfemaleFlower_specific_goterms, go_terms)
  }
}

for (i in rownames(CONmaleFlower_specific)){
  seqname<-str_split(i, "[.]")[[1]][1]
  seqname<-str_remove(seqname, "CON_")
  line<-con_annotation[con_annotation$V1 == seqname,]
  go_terms<-str_split(line$V2, ",")
  if(length(go_terms) != 0){
    go_terms = go_terms[[1]]
    CONmaleFlower_specific_goterms <-c(CONmaleFlower_specific_goterms, go_terms)
  }
}

for (i in rownames(CONleaf_specific)){
  seqname<-str_split(i, "[.]")[[1]][1]
  seqname<-str_remove(seqname, "CON_")
  line<-con_annotation[con_annotation$V1 == seqname,]
  go_terms<-str_split(line$V2, ",")
  if(length(go_terms) != 0){
    go_terms = go_terms[[1]]
    CONleaf_specific_goterms <-c(CONleaf_specific_goterms, go_terms)
  }
}

for (i in rownames(CONpetiole_specific)){
  seqname<-str_split(i, "[.]")[[1]][1]
  seqname<-str_remove(seqname, "CON_")
  line<-con_annotation[con_annotation$V1 == seqname,]
  go_terms<-str_split(line$V2, ",")
  if(length(go_terms) != 0){
    go_terms = go_terms[[1]]
    CONpetiole_specific_goterms <-c(CONpetiole_specific_goterms, go_terms)
  }
}

for (i in rownames(CONroot_specific)){
  seqname<-str_split(i, "[.]")[[1]][1]
  seqname<-str_remove(seqname, "CON_")
  line<-con_annotation[con_annotation$V1 == seqname,]
  go_terms<-str_split(line$V2, ",")
  if(length(go_terms) != 0){
    go_terms = go_terms[[1]]
    CONroot_specific_goterms <-c(CONroot_specific_goterms, go_terms)
  }
}

for (i in rownames(CONvegBud_specific)){
  seqname<-str_split(i, "[.]")[[1]][1]
  seqname<-str_remove(seqname, "CON_")
  line<-con_annotation[con_annotation$V1 == seqname,]
  go_terms<-str_split(line$V2, ",")
  if(length(go_terms) != 0){
    go_terms = go_terms[[1]]
    CONvegBud_specific_goterms <-c(CONvegBud_specific_goterms, go_terms)
  }
}

########

for (i in rownames(PLEfemaleFlower_specific)){
  seqname<-str_split(i, "[.]")[[1]][1]
  seqname<-str_remove(seqname, "PLE_")
  line<-ple_annotation[ple_annotation$V1 == seqname,]
  go_terms<-str_split(line$V2, ",")
  if(length(go_terms) != 0){
    go_terms = go_terms[[1]]
    PLEfemaleFlower_specific_goterms <-c(PLEfemaleFlower_specific_goterms, go_terms)
  }
}

for (i in rownames(PLEmaleFlower_specific)){
  seqname<-str_split(i, "[.]")[[1]][1]
  seqname<-str_remove(seqname, "PLE_")
  line<-ple_annotation[ple_annotation$V1 == seqname,]
  go_terms<-str_split(line$V2, ",")
  if(length(go_terms) != 0){
    go_terms = go_terms[[1]]
    PLEmaleFlower_specific_goterms <-c(PLEmaleFlower_specific_goterms, go_terms)
  }
}

for (i in rownames(PLEleaf_specific)){
  seqname<-str_split(i, "[.]")[[1]][1]
  seqname<-str_remove(seqname, "PLE_")
  line<-ple_annotation[ple_annotation$V1 == seqname,]
  go_terms<-str_split(line$V2, ",")
  if(length(go_terms) != 0){
    go_terms = go_terms[[1]]
    PLEeaf_specific_goterms <-c(PLEleaf_specific_goterms, go_terms)
  }
}

for (i in rownames(PLEpetiole_specific)){
  seqname<-str_split(i, "[.]")[[1]][1]
  seqname<-str_remove(seqname, "PLE_")
  line<-ple_annotation[ple_annotation$V1 == seqname,]
  go_terms<-str_split(line$V2, ",")
  if(length(go_terms) != 0){
    go_terms = go_terms[[1]]
    PLEpetiole_specific_goterms <-c(PLEpetiole_specific_goterms, go_terms)
  }
}

for (i in rownames(PLEroot_specific)){
  seqname<-str_split(i, "[.]")[[1]][1]
  seqname<-str_remove(seqname, "PLE_")
  line<-ple_annotation[ple_annotation$V1 == seqname,]
  go_terms<-str_split(line$V2, ",")
  if(length(go_terms) != 0){
    go_terms = go_terms[[1]]
    PLEroot_specific_goterms <-c(PLEroot_specific_goterms, go_terms)
  }
}

for (i in rownames(PLEvegBud_specific)){
  seqname<-str_split(i, "[.]")[[1]][1]
  seqname<-str_remove(seqname, "PLE_")
  line<-ple_annotation[ple_annotation$V1 == seqname,]
  go_terms<-str_split(line$V2, ",")
  if(length(go_terms) != 0){
    go_terms = go_terms[[1]]
    PLEvegBud_specific_goterms <-c(PLEvegBud_specific_goterms, go_terms)
  }
}


# a function to plot the top 25 GO terms from tissue specific GO terms
# splits the GO terms into different ontologies (MF, CC, BP)
plot_go_graph<-function(input_terms){
  df<-data.frame(table(input_terms), Term(names(table(input_terms))), Ontology(names(table(input_terms))))
  colnames(df)<-c("go_term", "count", "term_name", "ontology")
  df_mf <- df %>% filter(ontology == "MF")
  df_bp <- df %>% filter(ontology == "BP")
  df_cc <- df %>% filter(ontology == "CC")
  df_mf$term_name<-factor(df_mf$term_name, levels = df_mf$term_name[order(df_mf$count, decreasing = FALSE)])
  df_bp$term_name<-factor(df_bp$term_name, levels = df_bp$term_name[order(df_bp$count, decreasing = FALSE)])
  df_cc$term_name<-factor(df_cc$term_name, levels = df_cc$term_name[order(df_cc$count, decreasing = FALSE)])
  df_toplot<-rbind(df_mf[1:25,], df_bp[1:25,], df_cc[1:25,])
  df_toplot<-df_toplot %>% drop_na()
  ggplot(df_toplot, aes(x=term_name, y=count, colour=ontology, fill=ontology)) + 
    geom_bar(stat="identity") + 
    coord_flip()  + 
    theme(axis.text.y=element_text(size=12), 
          axis.title.x=element_blank(), 
          axis.title.y=element_blank())
}


# plot each tissue specific GO term graph for B. conchifolia
con_a<-plot_go_graph(CONfemaleFlower_specific_goterms)
con_b<-plot_go_graph(CONmaleFlower_specific_goterms)
con_c<-plot_go_graph(CONleaf_specific_goterms)
con_d<-plot_go_graph(CONpetiole_specific_goterms)
con_e<-plot_go_graph(CONroot_specific_goterms)
con_f<-plot_go_graph(CONvegBud_specific_goterms)

# plot each tissue specific GO term graph for B. plebeja
ple_a<-plot_go_graph(PLEfemaleFlower_specific_goterms)
ple_b<-plot_go_graph(PLEmaleFlower_specific_goterms)
ple_c<-plot_go_graph(PLEleaf_specific_goterms)
#ple_d<-plot_go_graph(PLEpetiole_specific_goterms)
ple_e<-plot_go_graph(PLEroot_specific_goterms)
ple_f<-plot_go_graph(PLEvegBud_specific_goterms)




# check how long each Go term list is for tables in manuscript
length(CONfemaleFlower_specific_goterms)
length(CONmaleFlower_specific_goterms)
length(CONleaf_specific_goterms)
length(CONpetiole_specific_goterms)
length(CONroot_specific_goterms)
length(CONvegBud_specific_goterms)

length(PLEfemaleFlower_specific_goterms)
length(PLEmaleFlower_specific_goterms)
length(PLEleaf_specific_goterms)
length(PLEpetiole_specific_goterms)
length(PLEroot_specific_goterms)
length(PLEvegBud_specific_goterms)

# arrange all figures together
# to go in supplemental due to size
figure_con <- ggarrange(con_a, con_b, con_c, con_d, con_e, con_f,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    ncol = 3, nrow = 2)


figure_ple <- ggarrange(ple_a, ple_b, ple_d, ple_e, ple_f,
                        labels = c("A", "B", "C", "D", "E"),
                        ncol = 3, nrow = 2)


pdf("supplementary_fig3_tissue_specific_goterms_con.pdf", width = 30, height = 25)
figure_con
dev.off()

pdf("supplementary_fig4_tissue_specific_goterms_ple.pdf", width = 35, height = 25)
figure_ple
dev.off()




#########################################################
###     uniquely expressed transcripts upset plot      ##
#########################################################

# set FPKM threshold
threshold <- 1

# use threshold to get list of transcript names which are expressed more than
con_female_flower_expressed<-rownames(con_avg_cpm)[con_avg_cpm[,1] >=threshold]
con_leaf_expressed<-rownames(con_avg_cpm)[con_avg_cpm[,2] >=threshold]
con_male_flower_expressed<-rownames(con_avg_cpm)[con_avg_cpm[,3] >=threshold]
con_petiole_expressed<-rownames(con_avg_cpm)[con_avg_cpm[,4] >=threshold]
con_root_expressed<-rownames(con_avg_cpm)[con_avg_cpm[,5] >=threshold]
con_vegbud_expressed<-rownames(con_avg_cpm)[con_avg_cpm[,6] >=threshold]

ple_female_flower_expressed<-rownames(ple_avg_cpm)[ple_avg_cpm[,1] >=threshold]
ple_leaf_expressed<-rownames(ple_avg_cpm)[ple_avg_cpm[,2] >=threshold]
ple_male_flower_expressed<-rownames(ple_avg_cpm)[ple_avg_cpm[,3] >=threshold]
ple_petiole_expressed<-rownames(ple_avg_cpm)[ple_avg_cpm[,4] >=threshold]
ple_root_expressed<-rownames(ple_avg_cpm)[ple_avg_cpm[,5] >=threshold]
ple_vegbud_expressed<-rownames(ple_avg_cpm)[ple_avg_cpm[,6] >=threshold]

# use those lists as input to UpSet plots
con_listInput <- list(`Female Flower` = con_female_flower_expressed,
                      `Leaf` = con_leaf_expressed,
                      `Male Flower` = con_male_flower_expressed,
                      `Petiole` = con_petiole_expressed,
                      `Root` = con_root_expressed,
                      `Vegetative Bud` = con_vegbud_expressed)

ple_listInput <- list(`Female Flower` = ple_female_flower_expressed,
                      `Leaf` = ple_leaf_expressed,
                      `Male Flower` = ple_male_flower_expressed,
                      `Petiole` = ple_petiole_expressed,
                      `Root` = ple_root_expressed,
                      `Vegetative Bud` = ple_vegbud_expressed)


pdf("figure4_shared_tissue_transcripts_con.pdf", width = 25, height = 15)
upset(fromList(con_listInput), order.by = "freq", nsets=6, text.scale = c(3, 3, 3, 3, 3, 2.5), point.size = 4)
dev.off()

pdf("figure5_shared_tissue_transcripts_ple.pdf", width = 25, height = 15)
upset(fromList(ple_listInput), order.by = "freq", nsets=6, text.scale = c(3, 3, 3, 3, 3, 2.5), point.size = 4)
dev.off()



#########################################################
###             sequence length histograms             ##
#########################################################


decontaminated_con<-"/Volumes/BACKUP3/trimmed_decontaminated/trinity/con/trinity_out_dir/Trinity.fasta"
contaminated_con<-"/Volumes/BACKUP3/balrog_download/blobplot/con_trimmed_trinity.fasta"

decontaminated_ple<-"/Volumes/BACKUP3/trimmed_decontaminated/trinity/ple/trinity_out_dir/Trinity.fasta"
contaminated_ple<-"/Volumes/BACKUP3/balrog_download/blobplot/ple_trimmed_trinity.fasta"


fasta_decon_con <- read.fasta(file = decontaminated_con)
lengths_decon_con<-getLength(fasta_decon_con) 
fasta_decon_ple <- read.fasta(file = decontaminated_ple)
lengths_decon_ple<-getLength(fasta_decon_ple) 

fasta_contaminated_con <- read.fasta(file = contaminated_con)
lengths_contaminated_con<-getLength(fasta_contaminated_con) 
fasta_contaminated_ple <- read.fasta(file = contaminated_ple)
lengths_contaminated_ple<-getLength(fasta_contaminated_ple) 

lengths_cont_df_con<-data.frame(lengths=lengths_contaminated_con, assembly="contaminated", species="B. conchifolia")
lengths_cont_df_ple<-data.frame(lengths=lengths_contaminated_ple, assembly="contaminated", species="B. plebeja")
lengths_decon_df_con<-data.frame(lengths=lengths_decon_con, assembly="decontaminated", species="B. conchifolia")
lengths_decon_df_ple<-data.frame(lengths=lengths_decon_ple, assembly="decontaminated", species="B. plebeja")

# print summary of lengths
summary(lengths_cont_df_con$lengths)
summary(lengths_cont_df_ple$lengths)
summary(lengths_decon_df_con$lengths)
summary(lengths_decon_df_ple$lengths)

# bind everything together to plot
lengths_all<-rbind(lengths_cont_df_con, lengths_cont_df_ple, lengths_decon_df_con, lengths_decon_df_ple)


pdf("figure3_sequence_length_histogram.pdf", width = 10, height = 7)
ggplot(lengths_all, aes(x=lengths, colour=assembly, fill=assembly)) + 
  geom_histogram(alpha=0.8, position = "dodge") + facet_grid(rows = vars(species)) +
  ylab("Number of transcripts") + 
  xlab("Transcript length") + 
  theme(legend.text=element_text(size=15),
        legend.title = element_blank(),
        axis.title.x=element_text(size=15), 
        axis.title.y=element_text(size=15),
        axis.text.x= element_text(size=15),
        axis.text.y= element_text(size=15),
        strip.text.y = element_text(size = 15, face="italic"))
dev.off()


#########################################################
###             CHS expression heatmap                 ##
#########################################################


# these are the expression values grepped out of avg_fpkm files for each sequence used in the phylogenetic tree

chs_expression<-data.frame(CON_TRINITY_DN7494_c0_g2_i3 = c(1079.21105353443,	1355.12934525355,	1578.51236944914,	530.944363270742,	2050.25822525919,	2350.4780647317),
                 PLE_TRINITY_DN4929_c0_g1_i2 = c(924.596399332199,	791.108169198113,	2903.54928095096,	408.091974087824,	1889.47347541453,	1911.29304864381),
                 CON_TRINITY_DN7494_c0_g1_i13 = c(125.100449245035,	198.371872561107,	27.5750096387708,	87.9648832920222,	576.921038683795,	497.63994605491),
                 PLE_TRINITY_DN4929_c0_g5_i2 = c(22.9793092298,	19.0577362046377,	50.1597318327533,	11.9099250336282,	635.852980027475,	139.94768603552),
                 CON_TRINITY_DN7494_c0_g4_i1 = c(0.157425406751893,	0.198647581483009,	0.112911650006348,	0.187629954492802,	2.42180462748974,	0.270607978612891),
                 PLE_TRINITY_DN4929_c0_g4_i3 = c(6.52592654631693,	2.06438283136467,	8.56467592589611,	9.50761752902403,	14.9072008153857,	6.77871413195937),
                 CON_TRINITY_DN7494_c0_g5_i1 = c(132.759671977685,	185.082362877332,	17.1353878223157,	13.6595518410768,	278.794311267313,	365.584617551438),
                 PLE_TRINITY_DN9097_c0_g2_i1 = c(193.225158625332,	38.9603206234101,	12.264579202113,	13.6518806166165,	493.23784385718,	606.763808065943)
                 )

chs_expression<-t(chs_expression)
colnames(chs_expression) <- c("Female flower"	,"Leaf",	"Male flower"	,"Petiole",	"Root"	,"Vegetative bud")

transcript_ids<-c("CON_TRINITY_DN7494_c0_g2_i3", 
                  "PLE_TRINITY_DN4929_c0_g1_i2", 
                  "CON_TRINITY_DN7494_c0_g1_i13",
                  "PLE_TRINITY_DN4929_c0_g5_i2", 
                  "CON_TRINITY_DN7494_c0_g4_i1", 
                  "PLE_TRINITY_DN4929_c0_g4_i3", 
                  "CON_TRINITY_DN7494_c0_g5_i1", 
                  "PLE_TRINITY_DN9097_c0_g2_i1")

group_info<-data.frame(paralog_group = c(rep(c("group1"), 2), rep(c("group2"), 2), rep(c("group3"), 2), rep(c("group4"), 2)) )
rownames(group_info) <-transcript_ids
group_info_cols<-list(paralog_group = c(group1 = "lightblue3", group2 = "pink2", group3 = "lightgoldenrod2", group4 = "darkseagreen2"))

pdf("figure8_chs_heatmap.pdf", width = 21, height = 15)
out<-pheatmap(log2(chs_expression), 
         cluster_rows = FALSE, 
         fontsize = 25, 
         angle_col = 45, 
         fontsize_col = 25, 
         fontsize_row = 20, 
         annotation_row=group_info, 
         annotation_names_row = FALSE,
         annotation_colors = group_info_cols
         )
dev.off()


#########################################################
###         annotation statistics barplot              ##
#########################################################


con_annot<-read_tsv("con_trinotate_annotation_report.xls")
ple_annot<-read_tsv("ple_trinotate_annotation_report.xls")


x<-ple_annot %>% 
  rowwise() %>% 
  filter(length(which(c(sprot_Top_BLASTX_hit, 
                        sprot_Top_BLASTP_hit, 
                        Pfam, 
                        eggnog, 
                        Kegg, 
                        gene_ontology_BLASTX, 
                        gene_ontology_BLASTP, 
                        gene_ontology_Pfam) == ".")) == 8)

completely_unannotated<-x %>% dplyr::select(transcript_id) %>% pull() %>% unique() %>% length()



con_blastx_annot<-con_annot %>% filter(sprot_Top_BLASTX_hit != ".") %>% dplyr::select(transcript_id) %>% pull() %>% unique()
con_blastp_annot<-con_annot %>% filter(sprot_Top_BLASTP_hit != ".") %>% dplyr::select(transcript_id) %>% pull() %>% unique()
con_pfam_annot<-con_annot %>% filter(Pfam != ".") %>% dplyr::select(transcript_id) %>% pull() %>% unique()
con_eggnog_annot<-con_annot %>% filter(eggnog != ".") %>% dplyr::select(transcript_id) %>% pull() %>% unique()
con_kegg_annot<-con_annot %>% filter(Kegg != ".") %>% dplyr::select(transcript_id) %>% pull() %>% unique()
con_GOblastx_annot<-con_annot %>% filter(gene_ontology_BLASTX != ".") %>% dplyr::select(transcript_id) %>% pull() %>% unique()
con_GOblastp_annot<-con_annot %>% filter(gene_ontology_BLASTP != ".") %>% dplyr::select(transcript_id) %>% pull() %>% unique()
con_GOpfam_annot<-con_annot %>% filter(gene_ontology_Pfam != ".") %>% dplyr::select(transcript_id) %>% pull() %>% unique()


ple_blastx_annot<-ple_annot %>% filter(sprot_Top_BLASTX_hit != ".") %>% dplyr::select(transcript_id) %>% pull() %>% unique()
ple_blastp_annot<-ple_annot %>% filter(sprot_Top_BLASTP_hit != ".") %>% dplyr::select(transcript_id) %>% pull() %>% unique()
ple_pfam_annot<-ple_annot %>% filter(Pfam != ".") %>% dplyr::select(transcript_id) %>% pull() %>% unique()
ple_eggnog_annot<-ple_annot %>% filter(eggnog != ".") %>% dplyr::select(transcript_id) %>% pull() %>% unique()
ple_kegg_annot<-ple_annot %>% filter(Kegg != ".") %>% dplyr::select(transcript_id) %>% pull() %>% unique()
ple_GOblastx_annot<-ple_annot %>% filter(gene_ontology_BLASTX != ".") %>% dplyr::select(transcript_id) %>% pull() %>% unique()
ple_GOblastp_annot<-ple_annot %>% filter(gene_ontology_BLASTP != ".") %>% dplyr::select(transcript_id) %>% pull() %>% unique()
ple_GOpfam_annot<-ple_annot %>% filter(gene_ontology_Pfam != ".") %>% dplyr::select(transcript_id) %>% pull() %>% unique()



con_listInput<-list(blastx = con_blastx_annot,
                    blasp = con_blastp_annot,
                    pfam = con_pfam_annot,
                    eggnog = con_eggnog_annot,
                    kegg = con_kegg_annot,
                    GOblastx = con_GOblastx_annot,
                    GOblastp = con_GOblastp_annot,
                    GOpfam = con_GOpfam_annot)


ple_listInput<-list(blastx = ple_blastx_annot,
                    blasp = ple_blastp_annot,
                    pfam = ple_pfam_annot,
                    eggnog = ple_eggnog_annot,
                    kegg = ple_kegg_annot,
                    GOblastx = ple_GOblastx_annot,
                    GOblastp = ple_GOblastp_annot,
                    GOpfam = ple_GOpfam_annot)


pdf("supplementary_fig1_shared_annotation_sources_con.pdf", width = 25, height = 15)
out<-upset(fromList(con_listInput), order.by = "freq", nsets=8, text.scale = c(3, 3, 3, 3, 3, 2.5), point.size = 4)
out
dev.off()

pdf("supplementary_fig2_shared_annotation_sources_ple.pdf", width = 25, height = 15)
out<-upset(fromList(ple_listInput), order.by = "freq", nsets=8, text.scale = c(3, 3, 3, 3, 3, 2.5), point.size = 4)
out
dev.off()


########
#GO term annotations plots

con_total<-con_annot %>% dplyr::select(transcript_id) %>% pull() %>% unique()
ple_total<-ple_annot %>% dplyr::select(transcript_id) %>% pull() %>% unique()


length(con_total)
length(con_blastx_annot)
length(con_blastp_annot)
length(con_pfam_annot)
length(con_eggnog_annot)
length(con_kegg_annot)
length(con_GOblastx_annot)
length(con_GOblastp_annot)
length(con_GOpfam_annot)

length(ple_total)
length(ple_blastx_annot)
length(ple_blastp_annot)
length(ple_pfam_annot)
length(ple_eggnog_annot)
length(ple_kegg_annot)
length(ple_GOblastx_annot)
length(ple_GOblastp_annot)
length(ple_GOpfam_annot)
  


con_annotation_stats<-data.frame(category=c("total",
                                            "sprot_Top_BLASTX_hit", 
                                            "sprot_Top_BLASTP_hit",
                                            "Pfam",
                                            "eggnog",
                                            "Kegg",
                                            "gene_ontology_BLASTX", 
                                            "gene_ontology_BLASTP", 
                                            "gene_ontology_Pfam"
                                            ),
                                 num_transcripts=c(length(con_total),
                                                   length(con_blastx_annot),
                                                   length(con_blastp_annot),
                                                   length(con_pfam_annot),
                                                   length(con_eggnog_annot),
                                                   length(con_kegg_annot),
                                                   length(con_GOblastx_annot),
                                                   length(con_GOblastp_annot),
                                                   length(con_GOpfam_annot)),
                                 species="B. conchifolia")




ple_annotation_stats<-data.frame(category=c("total",
                                            "sprot_Top_BLASTX_hit", 
                                            "sprot_Top_BLASTP_hit",
                                            "Pfam",
                                            "eggnog",
                                            "Kegg",
                                            "gene_ontology_BLASTX", 
                                            "gene_ontology_BLASTP", 
                                            "gene_ontology_Pfam"),
                                 num_transcripts=c(length(ple_total),
                                                   length(ple_blastx_annot),
                                                   length(ple_blastp_annot),
                                                   length(ple_pfam_annot),
                                                   length(ple_eggnog_annot),
                                                   length(ple_kegg_annot),
                                                   length(ple_GOblastx_annot),
                                                   length(ple_GOblastp_annot),
                                                   length(ple_GOpfam_annot)),
                                 species="B. plebeja")




annotation_stats<-rbind(con_annotation_stats, ple_annotation_stats)
pdf("figure6_annotation_category_barplot.pdf", height = 10, width = 15)
ggplot(annotation_stats, aes(x=category, y=num_transcripts, colour=species, fill=species)) + 
  geom_bar(stat="identity", position = "dodge") + 
  coord_flip() +
  labs(y = "Number of unigenes", x = "Annotation Category") +
  theme(legend.text=element_text(size=20, face="italic"),
        legend.title = element_blank(),
        axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20))
dev.off()





