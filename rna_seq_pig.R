#Adapted largely from:
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4934518/pdf/f1000research-5-9996.pdf
  
  #load("RNASeq.RData_pig")
  #install 'Rsubread' package. DO ONCE!
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
  
  #BiocManager::install("Rsubread")
library(Rsubread) #For alignment

#install 'edgeR'
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("edgeR")
library(edgeR) #For DE analysis

#Download reference fasta and annotation (gtf) from ensemble, ...
#here: //ftp.ensembl.org/pub/release-97/fasta/sus_scrofa/dna/
#and here: ftp://ftp.ensembl.org/pub/release-97/gtf/sus_scrofa/

setwd("/users/PAS1394/osu10255/")
  
#Step 1: Building an index
#Must provide a single FASTA file (eg. “genome.fa”)
#build index
#buildindex(basename="ssRsubread_index",reference="/users/PAS1394/osu10255/sus_scrofa_genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa", memory=128000)

#Load files: 
#fastq.files.R1<-list.files(path = "TrimmedReads_pig", pattern = "Ss.*R1_001.fq.gz$", full.names = TRUE)
#fastq.files.R2<-list.files(path = "TrimmedReads_pig", pattern = "Ss.*R2_001.fq.gz$", full.names = TRUE)


#Step 2: Aligning the reads 
#NOTE:Consider mapping with "Subjunc" function, enables exon-spanning and alternative splicing alignment 
#Map paired-end reads:

#Aligning all reads
subjunc(index="ssRsubread_index",
        readfile1 = fastq.files.R1, readfile2 = fastq.files.R2,
        nthreads = 28)

#Aligning only one pair of reads
subjunc(index="ssRsubread_index",
        readfile1 = "TrimmedReads_pig/Ss_Dc_3_CKDL190143536-1a-12_H723FBBXX_L008_R1_001.fq.gz", 
        readfile2 = "TrimmedReads_pig/Ss_Dc_3_CKDL190143536-1a-12_H723FBBXX_L008_R2_001.fq.gz",
        nthreads = 28)

#Check parameters used in alignment: 
args(subjunc)

##Summary of proportion of read alignment: 
bam.files <- list.files(path = "TrimmedReads_pig", pattern = "Ss.*BAM$", full.names = TRUE)
#props<-propmapped(bam.files[1], properlyPaired=TRUE)
#props<-propmapped(bam.files, properlyPaired=TRUE)
#write.table(props,"PhalignmentProportionsRsubread.txt", sep = "\t")


fc <- featureCounts(bam.files, annot.ext = "/users/PAS1394/osu10255/sus_scrofa_genome/Sus_scrofa.Sscrofa11.1.97.gtf.gz", 
                    isGTFAnnotationFile = TRUE, nthreads=28, isPairedEnd=TRUE, 
                    GTF.featureType = "gene")
# See what slots are stored in fc
names(fc)

## Take a look at the featurecounts stats
fc$stat
fc$annotation
fc$counts
fc$targets
annotation<-(fc$annotation)

##Counts 
head(fc$counts)


##Remane to reduce the filenames, for aesthetics:
library(stringr)
samplenames<-bam.files
samplenames<-word(samplenames, 2, sep = fixed('/')) #change the number, and the character  
samplenames<-word(samplenames, 1, sep = fixed('_CKDL190143536'))
samplenames
#samplenames<-gsub(paste("_0_illumina12index"), "", samplenames)
#samplenames<-gsub(paste("HMH2YBGX9_"), "", samplenames)
colnames(fc$counts)
colnames(fc$counts)<-samplenames
colnames(fc$counts)
colnames(fc$stat)[-1]<-samplenames
fc$stat

#transfrom fc to edgeR DGEList object (fc2)
fc2<-DGEList(counts = fc$counts, lib.size = colSums(fc$counts),
             norm.factors = rep(1,ncol(fc$counts)), samples = NULL,
             group = NULL, genes = fc$annotation, remove.zeros = FALSE)
fc2$counts
fc2$samples
fc2$genes

#Organizing Sus sascrofa gene ID/Entrez IDs

library(biomaRt)
######Annotation of DE genes:
library("biomaRt")
#Load protist database
ensembl<-useMart(biomart = "protists_mart", host = "protists.ensembl.org")
ensembl <- useMart("ensembl")
ensembl
#List datasets:
listDatasets(ensembl)
#Select dataset to use: sscrofa_gene_ensembl/
#ensembl =  useDataset("sscrofa_gene_ensembl", mart = ensembl)
ensembl <- useMart("ensembl",dataset="sscrofa_gene_ensembl")
#list attributes:
listAttributes(ensembl)
head(listAttributes(ensembl))
#getBM(attributes = c("broad_p_infestans","description"), filters = "broad_p_infestans",values = rownames(DEgenesPAMA3$table), mart = ensembl)
#getBM(attributes = c("broad_p_infestans","description"), filters = "broad_p_infestans",values = , mart = ensembl)


#or
features.RA <- getBM(attributes = c("ensembl_gene_id","go_id",
                                    "description"),
                     filters = c("ensembl_gene_id"),
                     values = rownames(head(fc2$genes)),
                     mart = ensembl)

features.RA <- getBM(attributes = c("ensembl_gene_id","go_id",
                                    "description"),
                     filters = c("ensembl_gene_id"),
                     values = rownames(toptags.RA),
                     mart = ensembl)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("org.Ss.eg.db")
library(org.Ss.eg.db)
columns(org.Ss.eg.db)
keytypes(org.Ss.eg.db)


#if (!"BiocManager" %in% rownames(installed.packages()))
#  install.packages("BiocManager")
#BiocManager::install(c("AnnotationHub", "Organism.dplyr"))


geneid <- rownames(fc2)
keytypes(Homo.sapiens) ##see available keytypes
columns(Homo.sapiens) ## see available columns/descriptions
genes <- AnnotationDbi::select(Homo.sapiens, keys=geneid, columns=c("SYMBOL","ENTREZID", "TXCHROM", "DEFINITION"),
                               keytype="ENSEMBL")
dim(genes)

head(genes)
#remove duplicates
genes <- genes[!duplicated(genes$ENSEMBL),]
fc2$genes<-genes


###Organising gene annotations
#install the 'Homo.sapiens' annotation package
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("Homo.sapiens")
library(Homo.sapiens) 
geneid <- rownames(fc2)
keytypes(Homo.sapiens) ##see available keytypes
columns(Homo.sapiens) ## see available columns/descriptions
genes <- AnnotationDbi::select(Homo.sapiens, keys=geneid, columns=c("SYMBOL","ENTREZID", "TXCHROM", "DEFINITION"),
                               keytype="ENSEMBL")
dim(genes)

head(genes)
#remove duplicates
genes <- genes[!duplicated(genes$ENSEMBL),]
fc2$genes<-genes

##SAMPLE INFO
#Add group info:
fc2$samples
treatment<-as.factor(rep(c("UTR","TR"), c(3,3)))
fc2$samples$group<-treatment
fc2$samples
#Add lane info:
#lane <- as.factor(rep(c(1,2,3,4), c(8,8,8,8)))
#lane <- as.factor(rep(c(1,2,3,4), c(6,6,6,6))) ##FOR SECOND ITERATION!!!!!!
#fc2$samples$lane<-lane

#double check the whole object:
fc2
#make a copy:
fc2_orig<-fc2

#####DATA PRE-PROCESSING:
##Transformations from the raw-scale
##CPM normalizaton:
#raw counts are converted to CPM and log-CPM values using the cpm function
cpm <- cpm(fc2)
lcpm <- cpm(fc2, log=TRUE)

L <- mean(fc2$samples$lib.size) * 1e-6 #average library size in Millions
M <- median(fc2$samples$lib.size) * 1e-6 #median lib size 

summary(lcpm)
summary(cpm)[,1:3]

###Removing genes that are lowly expressed
#get # of genes with zero counts across all samples:
table(rowSums(fc2$counts==0)==ncol(fc2$counts)) #~52% of genes had count zero across all samples!

###OPTION 1:
#Use the "filterByExpr" edgeR function. By default, the fxn keeps genes with about 10 read counts or more in a 
#minimum number of samples, where the number of samples is chosen according to the minimum group sample size. 
#The actual filtering uses CPM values rather than counts in order to avoid giving preference to samples with large library sizes. 
#For example, if the median library size is about 51 million and 10/51 ≈ 0.2, so the filterByExpr function keeps genes that have a CPM of 0.2 or
#more in at least three samples. 

#keep.exprs <- filterByExpr(fc2, group=group)
keep.exprs <- filterByExpr(fc2)
keep.exprs <- filterByExpr(fc2, min.count = 5, min.total.count = 5)
fc2 <- fc2[keep.exprs,, keep.lib.sizes=FALSE]
dim(fc2) ##18270 genes kept. 

#NOTE: Consider adjusting the read count cutoff. By default, the function keeps genes with about 10 read counts or more in a minimum number of samples 
keep.exprs <- filterByExpr(fc2, min.count = 2, min.total.count = 5)
fc2 <- fc2[keep.exprs,, keep.lib.sizes=FALSE]
dim(fc2) ##18270 genes kept. 

###OPTION 2:
# get # of genes in this dataset that have zero counts across all 6?? samples.
#fc2<-fc2_orig
#cpm <- cpm(fc2)
#lcpm <- cpm(fc2, log=TRUE)
#table(rowSums(fc2$counts==0)==6) ##change 9 to the number of samples you have in your dataset
#keep.exprs <- rowSums(cpm>1)>=3 #Filter by cpm > 1 in 3? or more samples 
#fc2 <- fc2[keep.exprs,, keep.lib.sizes=FALSE]
#dim(fc2) #14007 genes kept. Option 2 is more stringent. 

#Go with option 1
fc2<-fc2_orig
cpm <- cpm(fc2)
lcpm <- cpm(fc2, log=TRUE)

L <- mean(fc2$samples$lib.size) * 1e-6 #average library size in Millions
M <- median(fc2$samples$lib.size) * 1e-6 #median lib size


keep.exprs <- filterByExpr(fc2)
keep.exprs <- filterByExpr(fc2, min.count = 5, min.total.count = 5)
fc2 <- fc2[keep.exprs,, keep.lib.sizes=FALSE]
dim(fc2)
#make plots
#remove sample names, too many
library(RColorBrewer)
lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(fc2)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="") #density with lcpm from unfiltered data
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("topright", samplenames, text.col=col, bty="n")
lcpm2 <- cpm(fc2, log=TRUE) #NOTE: reculculating lcpm from filtered data!!! 
plot(density(lcpm2[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm2[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

##Normalization:
fc2_unNorm<-fc2 #create a copy of unnormalized data
fc2 <- calcNormFactors(fc2, method = "TMM")
fc2$samples
fc2_unNorm$samples

lcpm <- cpm(fc2_unNorm, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="Unnormalized data",ylab="Log-cpm")

lcpm <- cpm(fc2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="Normalized data",ylab="Log-cpm")

###MDS plots 
#To create the MDS plots, we assign different colours to the factors of interest.
#Dimensions 1 and 2 are examined using the colour grouping (i.e treatment).
lcpm <- cpm(fc2, log=TRUE)
par(mfrow=c(1,2))
col.group <- fc2$samples$group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
#col.lane <- lane
#levels(col.lane) <- brewer.pal(nlevels(col.lane), "Set2")
#col.lane <- as.character(col.lane)
plotMDS(lcpm, col=col.group)
plotMDS(lcpm, labels=fc2$samples$group, col=col.group)

#title(main="A. Sample groups")
#Dimensions 3 and 4 are examined using the colour grouping defined by sequencing lanes (batch).
#plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
#title(main="B. Sequencing lanes")
dev.off()


####PCA:
# transpose the data to have variables (genes) as columns
data_for_PCA <- t(fc2$counts)
dim(data_for_PCA)

## calculate MDS (matrix of dissimilarities)
mds <- cmdscale(dist(data_for_PCA), k=3, eig=TRUE)  
# k = the maximum dimension of the space which the data are to be represented in
# eig = indicates whether eigenvalues should be returned
mds$eig

##How many components can explain the variabilty?
# transform the Eigen values into percentage
eig_pc <- mds$eig * 100 / sum(mds$eig)
# plot the PCA
#png(file="~/PCA_PropExplainedVariance.png")
barplot(eig_pc,
        las=1,
        xlab="Dimensions", 
        ylab="Proportion of explained variance (%)", y.axis=NULL,
        col="darkgrey", names.arg = c(1,2,3,4,5,6))

#dev.off()
## calculate MDS
#mds <- cmdscale(dist(data_for_PCA)) # Performs MDS analysis 
#Samples representation
#png(file="~/PCA_Dim1vsDim2.png")
#plot(mds[,1], -mds[,2], type="n", xlab="Dimension 1", ylab="Dimension 2", main="")
#text(mds[,1], -mds[,2], rownames(mds), cex=0.8)

#CONCLUSSION: Both MDS plot & PCA support DE analysis between TR and UTR!!!!

###########################################
###STATISTICS!
###########################################
##FITTING MODELS TO DATA:
##Fitting genewise negative binomial generalized linear models (unlike limma lmFit that uses logcpm--assuming Normality--after adding 'voom' weights)
# create design matrix for differential expression analysis;
# if you wanted to account for batch here, you could simply include a batch
# term in the linear model at this step, e.g.:
#mod <- model.matrix(~0 + fc2$samples$lane + fc2$samples$group ) ##IF WE NEEDED TO ACCOUNT FOR LANE EFFECTS, IN THIS CASE NO

#Create a design matrix fir the linear models
#No reference, and no intercept, easier for pair-wise comparisons with contrast
mod <- model.matrix(~0+fc2$samples$group)
mod
rownames(mod) <- colnames(fc2)
colnames(mod)<-c("TR","UTR")
mod

#Estimate dispersion: 
fc3 <- estimateDisp(fc2, design = mod )
#plot biological coefficient of variation
plotBCV(fc3)

#make contrasts:
#Below, a +ve log FC denotes gene upregulated TR relative to UTR baseline  
contr.matrix <- makeContrasts(
  TRvsUTR = TR - UTR, 
  levels = colnames(mod))

contr.matrix

#Fit GLM or GLM QL:
#fit <- glmQLFit(fc3, mod)
fit <- glmFit(fc3, mod)

#Perform Hypothesis testing (is DE significant?) using Likelihood ratio test
#glf_TRvsUTR <-glmQLFTest(fit, contrast = contr.matrix[,"TRvsUTR"])
glf_TRvsUTR <-glmLRT(fit, contrast = contr.matrix[,"TRvsUTR"])

#perform multiple testing correction/p value adjustment.
#Why? Here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6099145/
#First, see the top DE genes:
topTags(glf_TRvsUTR)
#Get all DE genes with adjusted p.value cutoff
pvals_TRvsUTR<-topTags(glf_TRvsUTR, n = "Inf", adjust.method = "BH", sort.by = "PValue", p.value = 1 )
pvals_TRvsUTR<-pvals_TRvsUTR[[1]]

##Modified code for this Pig Dataset
##count genes with FDF<=0.5, but consider different FDRs 
#use the data (science) manipulation package 'dplyr'
#install.packages("dplyr")
library(dplyr)
#see top 25 genes, with FDR cutoff of 0.5
pvals_TRvsUTR %>%
  filter(FDR <= 0.5207219) %>%
  head(25)

##What is the P-value less than 0.05 that is associated with FDR
pvals_TRvsUTR %>% 
  filter(PValue <= 0.05) %>%
  tail()
#bottom genes:
pvals_TRvsUTR %>%
  filter(FDR <= 0.5207219) %>%
  count() #1134 genes
#get the number of genes with FDR cutoff of 0.5
pvals_TRvsUTR %>%
  filter(FDR <= 0.5207219) %>%
  count() #1134 genes

#GET TABLE OF DE WITH GENES OF INTEREST
de.genes.table<-pvals_TRvsUTR %>%
  filter(FDR <= 0.5207219)
genes.interest<-c("ENSSSCG00000009765", "ENSSSCG00000011570","ENSSSCG00000021203")
de.genes.interest.table<-de.genes.table[de.genes.table$GeneID %in% genes.interest, ]

de.genes.interest.table

de.genes.interest.table$PValue
###Useful graphical representations of differential expression results
#summary(decideTests(glf_TRvsUTR))
#Get number of up and down regulated genes:
summary(decideTests(glf_TRvsUTR, adjust.method = "BH", p.value = 0.5207219 ))
dt<-decideTests(glf_TRvsUTR, adjust.method = "BH", p.value = 0.5207219)
summary(dt)

#Get gene IDs (ENSEMBLE) of DE genes
de.genes <- which(dt[,1]!=0) #DE genes, 0 represents not DE genes
de.genes
length(de.genes)
##Get the DE genes from the hypothesis testing object
#head(glf_TRvsUTR$genes$SYMBOL[de.genes], n=20)
#head(glf_TRvsUTR$genes$ENSEMBL[de.genes], n=20)
#head(glf_TRvsUTR$genes$DEFINITION[de.genes], n=20)

#The magnitude of the differential expression changes can be visualized with a fitted model MD plot:
#Mean-difference plot
plotMD(glf_TRvsUTR, main = "TR vs UTR") #Significantly DE genes at a FDR of 5% are highlighted

library(ggplot2)
#create a dataframe with 'sign' column showing Up, Down, or Non.sig
sign.dat<-glf_TRvsUTR$table %>%
  dplyr::mutate(sign=case_when(logFC>0 & PValue < 0.05 ~ "Up",
                               logFC<0 & PValue < 0.05 ~ "Down",
                               PValue > 0.05 ~ "Non.sigf")
  )

head(sign.dat)
ggplot(sign.dat, aes(x = logCPM, y=logFC,col=sign)) +
  #  geom_point(alpha=0.4) + 
  geom_point() + 
  scale_colour_manual(values=c("blue","black","red")) + 
  labs(x = "Average log CPM", y = "log-fold-change") +
  theme(legend.position = c(0.9, 0.9), legend.title = element_blank()) +
  theme(
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14), 
    axis.text.x = element_text(size = 14), 
    axis.text.y = element_text(size = 14)
  ) + ylim(-3,12) ####ggplot changing y-axis

######interactive MD plot
#install Glimma package
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Glimma")

library(Glimma)
glMDPlot(glf_TRvsUTR, coef=1, status=dt, main=colnames(glf_TRvsUTR)[1],
         side.main="ENTREZID", counts=lcpm, groups=glf_TRvsUTR$samples$group, launch=TRUE)


#We use glmTreat to narrow down the list of DE genes and focus on genes that are more
#biologically meaningful. We test whether the differential expression is significantly above a
#log2-fold-change of log2 1.2, i.e., a fold-change of 1.2

tr <- glmTreat(fit, contrast=contr.matrix, lfc=log2(1.2))
topTags(tr)
summary(decideTests(tr))
plotMD(tr) 


###Heatmap of the top 100 DE genes
#install.packages("gplots")
library(gplots)
tr.vs.utr.topgenes <- pvals_TRvsUTR$GeneID[1:100]
i <- which(fc2$genes$GeneID %in% tr.vs.utr.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
#with gene names
heatmap.2(lcpm[i,], scale="row", margins = c(5,10),
          labRow=fc2$genes$GeneID[i], labCol=fc2$samples$group,
          col=mycol, trace="none", density.info="none", dendrogram="column")
heatmap.2(lcpm[i,], scale="row", margins = c(5,10),
          labRow=fc2$genes$GeneID[i], labCol=fc2$samples$group,
          col=mycol, trace="none", density.info="none", dendrogram="row")

###heatmap filtering one sample
heatmap.2(lcpm_filter, scale = "row", col = mycol, 
          labRow=fc2$genes$GeneID[i], labCol=fc2$samples$group[-1], 
          trace="none", density.info="none", dendrogram="column" )


#Note: All are up regulated in TR (DC infection)
#dev.off()


#####Pathway Analysis
##Gene set testing with camera
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata"))

#idx <- ids2indices(Hs.c2,id=rownames(fc3))
idx <- ids2indices(Hs.c2,id=fc3$genes$ENTREZID)
cam.TRvsUTR <- camera(fc3,idx,mod,contrast=contr.matrix[,1])
head(cam.TRvsUTR,200) 
#We can see lots of cancer genes 
#Also:
#DEBIASI_APOPTOSIS_BY_REOVIRUS_INFECTION_UP
#RAMJAUN_APOPTOSIS_BY_TGFB1_VIA_SMAD4_UP
barcodeplot(glf_TRvsUTR$table$logFC,
            index=idx[["DEBIASI_APOPTOSIS_BY_REOVIRUS_INFECTION_UP"]],
            index2=idx[["DEBIASI_APOPTOSIS_BY_REOVIRUS_INFECTION_DN"]],
            labels=c("UTR","TR"),
            main="DEBIASI_APOPTOSIS_BY_REOVIRUS_INFECTION",
            alpha=1)


######Pathway Analysis with Pathview
#install pathview:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("pathview")
library(pathview)
#Install gage
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("gage")
library(gage)

####Getting KEGG data 

####1. Finding your species in KEGG
#First, we need to make sure the species we work with is found in KEGG database. 
#You can go the following webpage and see the organism is list:
#http://rest.kegg.jp/list/organism

#In addition, ```pathview``` also carry a data matrix with KEGG supported species.  You can explore this matrix as well.
data(korg)
head(korg)




######PATHWAY ANALYSIS PIG:
########Combine GAGE and EdgeR analysis 
DE_entrez_pig <- read.delim("~/DE_entrez_pig", header=FALSE)
colnames(DE_entrez_pig)<-c("ensemble","entrez")
DEnames_entrez <- DE_entrez_pig$entrez # get the DE ENTREZ names 
DEnames_ensembl <- DE_entrez_pig$ensemble # get the DE ENS names 

fc_data<-glf_TRvsUTR$table #get all FC data
fc_data$ensemble<-rownames(fc_data)
head(fc_data)
fc_data_DE<-fc_data[fc_data$ensemble %in% DEnames_ensembl,] #get FC for DE genes only
head(fc_data_DE)
nrow(fc_data_DE)
#now merge with the convertion data.frame:
fc_data_DE_ensembl_entrez<-left_join(fc_data_DE,DE_entrez_pig,  by = 'ensemble')
rownames(fc_data_DE_ensembl_entrez)<-fc_data_DE_ensembl_entrez$entrez  
head(fc_data_DE_ensembl_entrez)
head(fc_data_DE_ensembl_entrez$logFC)
DE_foldchange_entrez<-fc_data_DE_ensembl_entrez$logFC
DE_foldchange_entrez<-data.frame(DE_foldchange_entrez)
rownames(DE_foldchange_entrez)<-rownames(fc_data_DE_ensembl_entrez)
colnames(DE_foldchange_entrez)<-"FC"
head(DE_foldchange_entrez)



#Get short species name
org <- "Sus suscrofa"
species <- unlist(sapply(1:ncol(korg), function(i) {
  agrep(org, korg[, i])
}))
korg[species, 1, drop = F]
korg[species, 3, drop = F] #ssc
korg[species, 5, drop = F]

#####2. Creating the KEGG dataset for GAGE analysis
kegg_hs <- kegg.gsets("ssc")
#kegg_hs <- kegg.gsets("hsa", id.type = "entrez")
kegg.gs <- kegg_hs$kg.sets[kegg_hs$sigmet.idx]


### GAGE analysis 
#fc.kegg.p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL) REALLYYY???
##OR
fc.kegg.p <- gage(DE_foldchange_entrez, gsets = kegg.gs, ref = NULL, samp = NULL)
head(fc.kegg.p$greater[,1:5], 20)
head(fc.kegg.p$less[,1:5], 20)

##Look for some immunity pathways, such as:
#ssc04650 Natural killer cell mediated cytotoxicity                                    
#ssc04657 IL-17 signaling pathway                                                       
#ssc04659 Th17 cell differentiation                                                      
#ssc04660 T cell receptor signaling pathway                                             
#ssc04662 B cell receptor signaling pathway
#ssc04210 Apoptosis

#visualize ssc04660 T cell receptor signaling pathway:
pv.out.list <- sapply("ssc04660", function(pid) pathview(
  gene.data =  DE_foldchange_entrez, pathway.id = pid,
  species = "ssc", out.suffix="utr_tr",
  gene.idtype="KEGG"))

#ssc04210 Apoptosis
pv.out.list <- sapply("ssc04210", function(pid) pathview(
  gene.data =  DE_foldchange_entrez, pathway.id = pid,
  species = "ssc", out.suffix="utr_tr",
  gene.idtype="KEGG"))

#ssc04650 Natural killer cell mediated cytotoxicity 
pv.out.list <- sapply("ssc04650", function(pid) pathview(
  gene.data =  DE_foldchange_entrez, pathway.id = pid,
  species = "ssc", out.suffix="utr_tr",
  gene.idtype="KEGG"))

#ssc04060 Cytokine-cytokine receptor interaction 
pv.out.list <- sapply("ssc04060", function(pid) pathview(
  gene.data =  DE_foldchange_entrez, pathway.id = pid,
  species = "ssc", out.suffix="utr_tr",
  gene.idtype="KEGG"))

#visualize the Jak-STAT pathway: ssc04630
pv.out.list <- sapply("ssc04630", function(pid) pathview(
  gene.data =  DE_foldchange_entrez, pathway.id = pid,
  species = "ssc", out.suffix="utr_tr",
  gene.idtype="KEGG"))

#Significance test: NOT WORKING!!! 
utr_tr_sig_kegg <-sigGeneSet(fc.kegg.p, outname="utr_tr_pig.kegg")

###Visualize the data using Pathview 

##1. Selecting the path ids for the upregulated set. 
greater_set <- native_kegg_fc$greater[, "q.val"] < 0.1 &
  !is.na(native_kegg_fc$greater[, "q.val"])
greater_ids <- rownames(native_kegg_fc$greater)[greater_set]
head(greater_ids)

##2. Selecting path ids for down-regulated set
less_set <- native_kegg_fc$less[, "q.val"] < 0.1 &
  !is.na(native_kegg_fc$less[,"q.val"])
less_ids <- rownames(native_kegg_fc$less)[less_set]
less_ids

##3. Combine up and down-regulated path ids. 
c(greater_ids,less_ids)
combine_ids <- substr(c(greater_ids, less_ids), 1, 8)
head(combine_ids)

##4. Visualization 
#Here we are going to get the first three pathways and visualize using pathview
pv.out.list <- sapply(combine_ids[1:3], function(pid) pathview(
  gene.data =  DE_foldchange_entrez, pathway.id = pid,
  species = "hsa", out.suffix="utr_tr",
  gene.idtype="KEGG"))

#visualize the apoptosis pathways: hsa04210 and hsa04215??
#Note: most of the genes in the pathway are upregulated!
pv.out.list <- sapply("hsa04210", function(pid) pathview(
  gene.data =  DE_foldchange_entrez, pathway.id = pid,
  species = "hsa", out.suffix="utr_tr",
  gene.idtype="KEGG"))

#visualize the NF-kappa B signaling pathway: hsa04064
pv.out.list <- sapply("hsa04064", function(pid) pathview(
  gene.data =  DE_foldchange_entrez, pathway.id = pid,
  species = "hsa", out.suffix="utr_tr",
  gene.idtype="KEGG"))

#####5. Mapping DE genes to individual pathways
#We can map genes that show statistically significant changes to individual pathways. 
#As an example, we can take apoptosis pathway and see which of the DE genes from EdgeR is present. 

pv_apoptosis <- pathview(gene.data = DE_foldchange_entrez, 
                         gene.idtype = "KEGG", 
                         pathway.id = "hsa04210", 
                         species = "hsa", 
                         out.suffix = "Apoptosis", 
                         keys.align = "y", 
                         kegg.native = T, 
                         match.data = T, 
                         key.pos = "topright")

#See how the output look like. 
head(pv_apoptosis)

#save.image(file = "RNASeq.RData")

#Converting an object to a table

write.table(de.genes, file = "DE_genes.txt", sep= "\t", row.names = TRUE, col.names = NA)
read.table("DE_genes.txt")

#Subsetting

DE_Ensembl_Genes<-DE_genes$X

#Converting this object to a txt file

write.table(DE_Ensembl_Genes, file = "DE_Ensembl_DAVID.txt", sep= "\t", row.names = TRUE, col.names = NA)
read.table("DE_Ensembl_DAVID.txt")

#Looking for the P-value

#GET TABLE OF DE WITH GENES OF INTEREST
de.genes.table<-pvals_TRvsUTR %>% 
  filter(FDR <= 0.5207219)
genes.interest<-c("ENSSSCG00000012077", "ENSSSCG00000008648","ENSSSCG00000017705","ENSSSCG00000007682","ENSSSCG00000016850","ENSSSCG00000005587","ENSSSCG00000012520","ENSSSCG00000017416","ENSSSCG00000005311","ENSSSCG00000012959","ENSSSCG00000020740","ENSSSCG00000030368","ENSSSCG00000027997","ENSSSCG00000030850","ENSSSCG00000006648")
de.genes.interest.table<-de.genes.table[de.genes.table$GeneID %in% genes.interest, ]

de.genes.interest.table

write.table(de.genes.interest.table, file = "DE_TopGenes_Pig.txt", sep= "\t", row.names = TRUE, col.names = NA)
read.table("DE_TopGenes_Pig.txt")

#Subsetting
de.genes.interest.table$PValue

#Converting the whole DEFs object to a table

write.table(de.genes.table, file = "AllDE_genes_Pig.txt", sep= "\t", row.names = TRUE, col.names = NA)
read.table("AllDE_genes_Pig.txt")

#Visualization of Functional Enrichment Result
#library(enrichplot)
#library(DOSE)
library(clusterProfiler)
# KEGG over-representation test: ONLY DE gnes list!!!
data(geneList, package="DOSE")
#gene <- names(geneList)[abs(geneList) > 2]
#gene<-names(DE_foldchange_entrez)
gene<-row.names(DE_foldchange_entrez)
kk <- enrichKEGG(gene         = gene,
                 organism     = 'ssc',
                 pvalueCutoff = 0.5)
kk

head(kk)

#dot plot:
dotplot(kk, showCategory=60) + ggtitle("dotplot for KEGG overrepresentaion")

#KEGG Gene Set Enrichment Analysis: ##Ranked (ALL?) gene list
#convert the edger.fc genelist ENSEMBLE IDs to ENTREZID

edger.fc.df<-as.data.frame(edger.fc)
head(edger.fc.df)
edger.fc.df$ENSEMBL<-rownames(edger.fc.df)
head(edger.fc.df)
#Now join:
head(ann)
head(edger.fc.df)
edger.fc.df_ann<-left_join(edger.fc.df,ann,  by = 'ENSEMBL') %>%
  na.omit() %>%
  distinct(ENTREZID, .keep_all= TRUE)
head(edger.fc.df_ann)
rownames(edger.fc.df_ann)<-edger.fc.df_ann$ENTREZID
head(edger.fc.df_ann)
#Now convert to an entrez gene list!
edger.fc.geneList.entrez<-edger.fc.df_ann$edger.fc 
names(edger.fc.geneList.entrez)<-rownames(edger.fc.df_ann)
head(edger.fc.geneList.entrez)
##NOTE: geneList is like vector of entrez gene ids with their Fold change values?
length(edger.fc.geneList.entrez)
edger.fc.geneList.entrez.sorted<-sort(edger.fc.geneList.entrez, decreasing = T)
kk2 <- gseKEGG(geneList     = edger.fc.geneList.entrez.sorted,
               organism     = 'hsa',
               nPerm        = 10000,
               minGSSize    = 5,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
kk2[,1:3]
dotplot(kk2, showCategory=60) + ggtitle("dotplot for KEGG enrichment")
