library("phyloseq") 
packageVersion("phyloseq")
library(phyloseq)
library(ggplot2)
library("DESeq2")
packageVersion("DESeq2")
library("ggplot2")

#otu_A
setwd("/Users/shadieshghi/Desktop/new_brendan/kungsängen_output/phyloseq_otu")
phyloseq_f_otua=readRDS("/Users/shadieshghi/Desktop/new_brendan/kungsängen_output/phyloseq_otu/phyloseq_otuA_fungi.rds")
phyloseq_p_otua=readRDS("/Users/shadieshghi/Desktop/new_brendan/kungsängen_output/phyloseq_otu/phyloseq_otuA_protists.rds")
phyloseq_otua = merge_phyloseq(phyloseq_f_otua, phyloseq_p_otua)
phyloseq_otua
phyloseq_otua=subset_samples(phyloseq_otua, Sites != "None")
cds_otua = phyloseq_to_deseq2(phyloseq_otua,  ~Sites)
cds_otua= DESeq(cds_otua, test = "Wald", fitType = "parametric")
head(sample_data(phyloseq_otua))
res_otua= results(cds_otua, cooksCutoff = TRUE)

#res1= results(cds_pro, cooksCutoff = FALSE)
alpha= 0.05
sigtab_otua= res_otua[which(res_otua$padj < alpha),]
sigtab_otua = cbind(as(sigtab_otua, "data.frame"), as(tax_table(phyloseq_otua)[rownames(sigtab_otua), ], "matrix"))
head(sigtab_otua)
dim(sigtab_otua)
write.table(sigtab_otua, "/Users/shadieshghi/Desktop/new_brendan/kungsängen_output/phyloseq_otu/sigtab_otua.txt", row.names =T , col.names = T, sep = "\t")


library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab_otua$log2FoldChange, sigtab_otua$phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_otua$phylum = factor(as.character(sigtab_otua$phylum), levels=names(x))
# Genus order
x = tapply(sigtab_otua$log2FoldChange, sigtab_otua$order, function(x) max(x))
x = sort(x, TRUE)
sigtab_otua$order = factor(as.character(sigtab_otua$order), levels=names(x))
ggplot(sigtab_otua, aes(x=order, y=log2FoldChange, color=phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(text = element_text(size = 30))+
  labs(title = "OTU_A")

########################################################################
########################################################################
#otu_s
setwd("/Users/shadieshghi/Desktop/new_brendan/kungsängen_output/phyloseq_otu")
phyloseq_f_otus=readRDS("/Users/shadieshghi/Desktop/new_brendan/kungsängen_output/phyloseq_otu/phyloseq_otuS_fungi.rds")
phyloseq_p_otus=readRDS("/Users/shadieshghi/Desktop/new_brendan/kungsängen_output/phyloseq_otu/phyloseq_otuS_protists.rds")
phyloseq_otus = merge_phyloseq(phyloseq_f_otus, phyloseq_p_otus)
phyloseq_otus
phyloseq_otus=subset_samples(phyloseq_otus, Sites != "None")
cds_otus = phyloseq_to_deseq2(phyloseq_otus,  ~Sites)
cds_otus= DESeq(cds_otus, test = "Wald", fitType = "parametric")

res_otus= results(cds_otus, cooksCutoff = TRUE)

#res1= results(cds_pro, cooksCutoff = FALSE)
alpha= 0.05
sigtab_otus= res_otus[which(res_otus$padj < alpha),]
sigtab_otus = cbind(as(sigtab_otus, "data.frame"), as(tax_table(phyloseq_otus)[rownames(sigtab_otus), ], "matrix"))
head(sigtab_otus)
sigtab_otus
dim(sigtab_otus)
write.table(sigtab_otus, "/Users/shadieshghi/Desktop/new_brendan/kungsängen_output/phyloseq_otu/sigtab_otus.txt", row.names =T , col.names = T, sep = "\t")


library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab_otus$log2FoldChange, sigtab_otus$phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_otus$phylum = factor(as.character(sigtab_otus$phylum), levels=names(x))
# Genus order
x = tapply(sigtab_otus$log2FoldChange, sigtab_otus$order, function(x) max(x))
x = sort(x, TRUE)
sigtab_otus$order = factor(as.character(sigtab_otus$order), levels=names(x))
ggplot(sigtab_otus, aes(x=order, y=log2FoldChange, color=phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(text = element_text(size = 30))+
  labs(title = "OTU_S")

##################################################
# Phylum order
x = tapply(sigtab_otus$log2FoldChange, sigtab_otus$kingdom, function(x) max(x))
x = sort(x, TRUE)
sigtab_otus$kingdom = factor(as.character(sigtab_otus$kingdom), levels=names(x))
# Genus order
x = tapply(sigtab_otus$log2FoldChange, sigtab_otus$order, function(x) max(x))
x = sort(x, TRUE)
sigtab_otus$order = factor(as.character(sigtab_otus$order), levels=names(x))
ggplot(sigtab_otus, aes(x=phylum, y=log2FoldChange, color=kingdom)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(text = element_text(size = 30))+
  labs(title = "OTU_S")

#################################################


x = tapply(sigtab_otus$log2FoldChange, sigtab_otus$kingdom, function(x) max(x))
x = sort(x, TRUE)
sigtab_otus$kingdom = factor(as.character(sigtab_otus$kingdom), levels=names(x))
# Genus order
x = tapply(sigtab_otus$log2FoldChange, sigtab_otus$order, function(x) max(x))
x = sort(x, TRUE)
sigtab_otus$order = factor(as.character(sigtab_otus$order), levels=names(x))
ggplot(sigtab_otus, aes(x=order, y=log2FoldChange, color=kingdom)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(text = element_text(size = 30))+
  labs(title = "OTU_S")

#otu_c
setwd("/Users/shadieshghi/Desktop/new_brendan/kungsängen_output/phyloseq_otu")
phyloseq_f_otuc=readRDS("/Users/shadieshghi/Desktop/new_brendan/kungsängen_output/phyloseq_otu/phyloseq_otuC_fungi.rds")
phyloseq_p_otuc=readRDS("/Users/shadieshghi/Desktop/new_brendan/kungsängen_output/phyloseq_otu/phyloseq_otuC_protists.rds")
phyloseq_otuc = merge_phyloseq(phyloseq_f_otuc, phyloseq_p_otuc)
phyloseq_otuc
phyloseq_otuc=subset_samples(phyloseq_otuc, Sites != "None")
cds_otuc = phyloseq_to_deseq2(phyloseq_otuc,  ~Sites)
cds_otuc= DESeq(cds_otuc, test = "Wald", fitType = "parametric")

res_otuc= results(cds_otuc, cooksCutoff = TRUE)

#res1= results(cds_pro, cooksCutoff = FALSE)
alpha= 0.05
sigtab_otuc= res_otuc[which(res_otuc$padj < alpha),]
sigtab_otuc = cbind(as(sigtab_otuc, "data.frame"), as(tax_table(phyloseq_otuc)[rownames(sigtab_otuc), ], "matrix"))
head(sigtab_otuc)
dim(sigtab_otuc)
write.table(sigtab_otuc, "/Users/shadieshghi/Desktop/new_brendan/kungsängen_output/phyloseq_otu/sigtab_otuc.txt", row.names =T , col.names = T, sep = "\t")


library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab_otuc$log2FoldChange, sigtab_otuc$phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_otuc$phylum = factor(as.character(sigtab_otuc$phylum), levels=names(x))
# Genus order
x = tapply(sigtab_otuc$log2FoldChange, sigtab_otuc$order, function(x) max(x))
x = sort(x, TRUE)
sigtab_otuc$order = factor(as.character(sigtab_otuc$order), levels=names(x))
ggplot(sigtab_otuc, aes(x=order, y=log2FoldChange, color=phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(text = element_text(size = 30))+
  labs(title = "OTU_C")