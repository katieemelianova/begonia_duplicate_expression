library(ggplot2)
library(reshape)



###################### orthogroup histogram and bubble plot
######################

ortho<-read.table("/Users/katie/Desktop/Bg/begonia_duplicate_expression/orthogroup_size_table.txt", header = FALSE)
colnames(ortho) <- c("og", "c", "p")

a<-ortho %>% dplyr::select(c, p) %>% table() %>% data.frame()

# bubble plot of gene family sizes
ggplot(a, aes(x=c, y=p, size=Freq)) +
  geom_point(alpha=0.5) +
  scale_size(range=c(0, 30))


hist(ortho$c)
hist(ortho$p)

con_orthogroups<-data.frame(count=ortho$c, species="con")
ple_orthogroups<-data.frame(count=ortho$p, species="ple")
orthogroups<-rbind(test1, test2)


ggplot(orthogroups, aes(x=count, colour=species, fill=species)) + geom_histogram(alpha=0.8, position = "dodge", bins=11)

###################### GO plots
######################
library(edgeR)
#library(baySeq)

#read in featurecounts output
con_counts<-read.table("/Users/katie/Desktop/Bg/begonia_duplicate_expression/featurecount/con/con_counts", header=TRUE)
ple_counts<-read.table("/Users/katie/Desktop/Bg/begonia_duplicate_expression/featurecount/ple/ple_counts", header=TRUE)


# gen the gene IDs (these are unique, unlike the Chr column)
# make them shorter, chopping off front diuplicated naming pattern
con_geneid<-sapply(con_counts$Geneid, function(x) str_split(x, "~~")[[1]][2]) %>% as.vector()
ple_geneid<-sapply(ple_counts$Geneid, function(x) str_split(x, "~~")[[1]][2]) %>% as.vector()

# prepend them with CON_ and PLE_ to match that of Orthogroups
con_names<-paste("CON", con_geneid, sep="_")
ple_names<-paste("PLE", ple_geneid, sep="_")


# remove all the path names fro  the column names to make them shorter
con_cols<-colnames(con_counts) %>% 
  str_remove("X.home.kemelianova.to_download.trimmed_decontaminated.star.") %>% 
  str_remove("Aligned.sortedByCoord.out.bam")

ple_cols<-colnames(ple_counts) %>% 
  str_remove("X.home.kemelianova.to_download.trimmed_decontaminated.star.") %>% 
  str_remove("Aligned.sortedByCoord.out.bam")

# rename columns using the shorter version
colnames(con_counts)<-con_cols
colnames(ple_counts)<-ple_cols

con_lengths<-con_counts %>% select(Length) %>% set_rownames(con_names)
ple_lengths<-ple_counts %>% select(Length) %>% set_rownames(ple_names)

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
plotMDS(con_d, method="bcv", col=as.numeric(con_d$samples$group))



ple_d<-DGEList(counts= ple_counts,group=factor(ple_groups))
ple_d_lengths<-ple_lengths %>% data.frame()
rownames(ple_d_lengths)<-rownames(ple_d$counts)
ple_d$samples$lib.size<-colSums(ple_d$counts)
ple_d<-calcNormFactors(ple_d, method="TMM")
plotMDS(ple_d, method="bcv", col=as.numeric(ple_d$samples$group))


CONfemaleFlower<-con_d$counts %>% data.frame() %>% select(CONfemaleFlower1, CONfemaleFlower2, CONfemaleFlower3) %>% rpkmByGroup(gene.length=con_d_lengths)  %>% data.frame()
CONleaf<-con_d$counts %>% data.frame() %>% select(CONleaf1, CONleaf2, CONleaf3) %>% rpkmByGroup(gene.length=con_d_lengths) %>% data.frame()
CONmaleFlower<-con_d$counts %>% data.frame() %>% select(CONmaleFlower1, CONmaleFlower2, CONmaleFlower3) %>% rpkmByGroup(gene.length=con_d_lengths) %>% data.frame()
CONpetiole<-con_d$counts %>% data.frame() %>% select(CONpetiole1, CONpetiole2, CONpetiole3) %>% rpkmByGroup(gene.length=con_d_lengths) %>% data.frame()
CONroot<-con_d$counts %>% data.frame() %>% select(CONroot1, CONroot3, CONvegBud1) %>% rpkmByGroup(gene.length=con_d_lengths) %>% data.frame()
CONvegBud<-con_d$counts %>% data.frame() %>% select(CONvegBud1, CONvegBud2, CONvegBud3) %>% rpkmByGroup(gene.length=con_d_lengths) %>% data.frame()

PLEfemaleFlower<-ple_d$counts %>% data.frame() %>% select(PLEfemaleFlower1, PLEfemaleFlower2, PLEfemaleFlower3) %>% rpkmByGroup(gene.length=ple_d_lengths)  %>% data.frame()
PLEleaf<-ple_d$counts %>% data.frame() %>% select(PLEleaf1, PLEleaf2, PLEleaf3) %>% rpkmByGroup(gene.length=ple_d_lengths)  %>% data.frame()
PLEmaleFlower<-ple_d$counts %>% data.frame() %>% select(PLEmaleFlower1, PLEmaleFlower2, PLEmaleFlower3) %>% rpkmByGroup(gene.length=ple_d_lengths)  %>% data.frame()
PLEpetiole<-ple_d$counts %>% data.frame() %>% select(PLEpetiole1, PLEpetiole2) %>% rpkmByGroup(gene.length=ple_d_lengths)  %>% data.frame()
PLEroot<-ple_d$counts %>% data.frame() %>% select(PLEroot1, PLEroot2, PLEroot3) %>% rpkmByGroup(gene.length=ple_d_lengths) %>% data.frame()
PLEvegBud<-ple_d$counts %>% data.frame() %>% select(PLEvegBud1, PLEvegBud2, PLEvegBud3) %>% rpkmByGroup(gene.length=ple_d_lengths)  %>% data.frame()


# bind con and ple together per species
con_avg_cpm<-cbind(CONfemaleFlower, CONleaf, CONmaleFlower, CONpetiole, CONroot, CONvegBud)
ple_avg_cpm<-cbind(PLEfemaleFlower, PLEleaf, PLEmaleFlower, PLEpetiole, PLEroot, PLEvegBud)

# transfer row names
rownames(con_avg_cpm)<-rownames(con_d$counts)
rownames(ple_avg_cpm)<-rownames(ple_d$counts)

# assign colnames
colnames(con_avg_cpm)<-c("CONfemaleFlower", "CONleaf", "CONmaleFlower", "CONpetiole", "CONroot", "CONvegBud")
colnames(ple_avg_cpm)<-c("PLEfemaleFlower", "PLEleaf", "PLEmaleFlower", "PLEpetiole", "PLEroot", "PLEvegBud")


#rownames(con_avg_cpm) <- str_remove(rownames(con_avg_cpm), ".p1")

CONfemaleFlower_specific<-con_avg_cpm %>% filter(CONfemaleFlower > 1 & CONroot < 1 & CONpetiole < 1 & CONvegBud < 1 & CONmaleFlower < 1 & CONleaf < 1)
CONmaleFlower_specific<-con_avg_cpm %>% filter(CONmaleFlower > 1 & CONroot < 1 & CONpetiole < 1 & CONvegBud < 1 & CONleaf < 1 & CONfemaleFlower < 1)
CONleaf_specific<-con_avg_cpm %>% filter(CONleaf > 1 & CONmaleFlower < 1 & CONroot < 1 & CONpetiole < 1 & CONvegBud < 1 & CONfemaleFlower < 1)
CONpetiole_specific<-con_avg_cpm %>% filter(CONpetiole > 1 & CONleaf < 1 & CONmaleFlower < 1 & CONroot < 1 & CONvegBud < 1 & CONfemaleFlower < 1)
CONroot_specific<-con_avg_cpm %>% filter(CONroot > 1 & CONpetiole < 1 & CONleaf < 1 & CONmaleFlower < 1 & CONvegBud < 1 & CONfemaleFlower < 1)
CONvegBud_specific<-con_avg_cpm %>% filter(CONvegBud > 1 & CONroot < 1 & CONpetiole < 1 & CONleaf < 1 & CONmaleFlower < 1 & CONfemaleFlower < 1)


con_annotation<-read.table("/Users/katie/Desktop/Bg/begonia_duplicate_expression/download_updated/trinotate/con/con_go_annotations_trans.txt")

CONfemaleFlower_specific_goterms<-c()
CONmaleFlower_specific_goterms<-c()
CONleaf_specific_goterms<-c()
CONpetiole_specific_goterms<-c()
CONroot_specific_goterms<-c()
CONvegBud_specific_goterms<-c()

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



library(GO.db)
Term("GO:0016021")
Ontology("GO:0016021")

plot_go_graph<-function(input_terms){
  df<-data.frame(table(input_terms), Term(names(table(input_terms))), Ontology(names(table(input_terms))))
  colnames(df)<-c("go_term", "count", "term_name", "ontology")
  df_mf <- df %>% filter(ontology == "MF")
  df_bp <- df %>% filter(ontology == "BP")
  df_cc <- df %>% filter(ontology == "CC")
  df_mf$term_name<-factor(df_mf$term_name, levels = df_mf$term_name[order(df_mf$count, decreasing = FALSE)])
  df_bp$term_name<-factor(df_bp$term_name, levels = df_bp$term_name[order(df_bp$count, decreasing = FALSE)])
  df_cc$term_name<-factor(df_cc$term_name, levels = df_cc$term_name[order(df_cc$count, decreasing = FALSE)])
  df_toplot<-rbind(df_mf[1:30,], df_bp[1:30,], df_cc[1:30,])
  ggplot(df_toplot, aes(x=term_name, y=count, colour=ontology, fill=ontology)) + geom_bar(stat="identity") + coord_flip() 
}

plot_go_graph(CONfemaleFlower_specific_goterms)
plot_go_graph(CONmaleFlower_specific_goterms)
plot_go_graph(CONleaf_specific_goterms)
plot_go_graph(CONpetiole_specific_goterms)
plot_go_graph(CONroot_specific_goterms)
plot_go_graph(CONvegBud_specific_goterms)



###################### 
######################

#GO term annotations plots
# only problem is that the report generated by trinonotate doesnt provide info about sharing of annotations sources
# no idea if these overlap


# CONCHIFOLIA
#gene_id	20048
#transcript_id	20048
#prot_coords	16756
#prot_id	16756
#sprot_Top_BLASTX_hit	16522
#gene_ontology_BLASTX	16211
#Kegg	15060
#sprot_Top_BLASTP_hit	14233
#gene_ontology_BLASTP	13919
#Pfam	13624
#gene_ontology_Pfam	8953
#eggnog	278


# PLEBEJA
#transcript_id	22241
##gene_id	22241
#prot_id	17771
#prot_coords	17771
#sprot_Top_BLASTP_hit	15269
#gene_ontology_BLASTP	14970
#Kegg	14828
#Pfam	14305
#gene_ontology_Pfam	9357
#sprot_Top_BLASTX_hit	3138
#gene_ontology_BLASTX	3073
#eggnog	281




con_annotation_stats<-data.frame(category=c("gene_id", 
                                            "transcript_id", 
                                            "prot_coords", 
                                            "prot_id", 
                                            "sprot_Top_BLASTX_hit", 
                                            "gene_ontology_BLASTX", 
                                            "Kegg", 
                                            "sprot_Top_BLASTP_hit", 
                                            "gene_ontology_BLASTP", 
                                            "Pfam", 
                                            "gene_ontology_Pfam", 
                                            "eggnog"),
                                 num_transcripts=c(20048, 
                                                   20048, 
                                                   16756, 
                                                   16756, 
                                                   16522, 
                                                   16211, 
                                                   15060, 
                                                   14233, 
                                                   13919, 
                                                   13624, 
                                                   8953, 
                                                   278))




ple_annotation_stats<-data.frame(category=c("transcript_id", 
                                            "gene_id", 
                                            "prot_id", 
                                            "prot_coords", 
                                            "sprot_Top_BLASTP_hit", 
                                            "gene_ontology_BLASTP", 
                                            "Kegg", 
                                            "Pfam", 
                                            "gene_ontology_Pfam", 
                                            "sprot_Top_BLASTX_hit", 
                                            "gene_ontology_BLASTX", 
                                            "eggnog"),
                                 num_transcripts=c(22241, 
                                                   22241, 
                                                   17771, 
                                                   17771, 
                                                   15269, 
                                                   14970, 
                                                   14828, 
                                                   14305, 
                                                   9357, 
                                                   3138, 
                                                   3073, 
                                                   281))




annotation_stats<-inner_join(con_annotation_stats, ple_annotation_stats, by="category")
annotation_stats<-melt(annotation_stats)
colnames(annotation_stats)<-c("category", "species", "count")
ggplot(annotation_stats, aes(x=category, y=count, colour=species, fill=species)) + geom_bar(stat="identity", position = "dodge") + coord_flip() 







