
# this table was generated using count_gene_families.py (same dir)
ortho_sizes<-read.table("/Users/katie/Desktop/Bg/begonia_duplicate_expression/orthogroup_size_table.txt", header = FALSE)
colnames(ortho_sizes) <- c("og", "c", "p")

ortho<-read_tsv("/Users/katie/Desktop/Bg/begonia_duplicate_expression/Orthogroups/Orthogroups.tsv")
colnames(ortho)<-c("orthogroup", "con", "ple")

# filter orthogroups where conchifolia is not expanded and plebeja is
expanded_in_ple<-ortho_sizes %>% filter(c == 1 & p > 5) %>% dplyr::select(og) %>% pull()

# filter out orthogroup memeber info for these expandded orthogroups
orthogroups_ple_expanded<-ortho %>% filter(orthogroup %in% expanded_in_ple)


split_contigs<-str_split(orthogroups_ple_expanded$ple, ", ", n = Inf, simplify = FALSE)
contig_per_orthogroup<-sapply(split_contigs, function(x) x[[1]][1] %>% str_remove("PLE_"))

ple_annotation<-read.table("/Users/katie/Desktop/Bg/begonia_duplicate_expression/download_updated/trinotate/ple/ple_go_annotations_trans.txt")
go_terms_ple_expanded<-ple_annotation %>% filter(V1 %in% contig_per_orthogroup)
go_terms_ple_expanded<-str_split(go_terms_ple_expanded$V2, ",")
sapply(go_terms_ple_expanded, function(x) print(Term(x[[1]])))


#library(GO.db)
#Term("GO:0016021")
#Ontology("GO:0016021")