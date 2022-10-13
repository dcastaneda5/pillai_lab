library(BiocManager)
library(dada2)

#This script has been adapted from the DADA2 tutorial. In here, I've added the parameters separately in variables I can easily change myself. It is meant for
#ama1d2 and ama1d3 data for amplicon deep sequencing on Plasmodium falciparum data. Additionally, it creates a new file in the output directory where it
#adds up the number of reads, infers the proportion of each haplotype.

path<-"/Users/danielcm/Desktop/Gondar/illumina/raw_fastq/MS3194860-600V3-19-Sep-2022-PillaiLab-3447444/ama1d3_fastq/trimmomatic_output_4_20/"
out_path<-"/Users/danielcm/Desktop/Gondar/illumina/raw_fastq/MS3194860-600V3-19-Sep-2022-PillaiLab-3447444/ama1d3_fastq/trimmomatic_output_4_20/"

#My parameters below ------------------
ama1d3_length = 516
primers_ama1d3_length = c(13,12)
primers_ama1d3_sum = 25
primers_ama1d2_sum = 42
primers_ama1d2_length = c(22,20)
ama1d2_length = 479
myMinOverlap = 12
myTruncLen_ama1d2 = c(280, ama1d2_length+myMinOverlap-280)
myTruncLen_ama1d3 = c(280, ama1d3_length+myMinOverlap-280)
myMaxEE = c(2,2)
myTruncQ = 5
myMaxConsistF = 10
myMaxConsistR = 10
myPooling = TRUE
myOmegaA = 1e-40
myChimeraMethod1 = "consensus"
myChimeraMethod2 = "pooled"
myChimeraMethod3 = "per-sample"

#-------------------

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
fnFs
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 3) #the 4 tells where to split the name of the file, which in this case is the 4th underscore.
sample.names
#plotQualityProfile(fnFs[3:4])
#plotQualityProfile(fnRs[3:4])
#dev.off()

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtFs
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=myTruncLen_ama1d3,
                     maxN=0, maxEE=myMaxEE, truncQ=myTruncQ, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE,trimLeft = primers_ama1d3_length)
out

errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST = myMaxConsistF, verbose = TRUE)
errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST = myMaxConsistR, verbose = TRUE)
plotErrors(errF, nominalQ=TRUE)

exists <- file.exists(filtFs) & file.exists(filtRs) #only working with the files that passed the filters
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]

dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool = myPooling, OMEGA_A=myOmegaA)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool = myPooling, OMEGA_A=myOmegaA)

dadaFs[[1]]$denoised #sample 1 forward
dadaRs[[1]] #sample 1 reverse

#dadaFs$NP

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, minOverlap = myMinOverlap)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab<-makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab2<-seqtab[,nchar(colnames(seqtab)) %in% ama1d3_length-primers_ama1d3_sum] #479 bp only for ama1d2! change for ama1d3 (516 bp)!

seqtab.nochim <- removeBimeraDenovo(seqtab2, method=myChimeraMethod1, multithread=TRUE, verbose=TRUE)
#seqtab.nochim <- removeBimeraDenovo(seqtab2, method=myChimeraMethod2, multithread=TRUE, verbose=TRUE)
#seqtab.nochim <- removeBimeraDenovo(seqtab2, method=myChimeraMethod3, multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#View(seqtab.nochim)
seqtab_df = as.data.frame(seqtab.nochim)
sum_reads = rowSums(seqtab_df) #adds up the reads of each row into the vector
total_haplotypes = ncol(seqtab_df) #prints the number of columns
columns_df = ncol(seqtab_df)
columns_df #number of haplotypes detected
length(sum_reads) #number of samples
length(seqtab_df[,1]) #number of samples
View(seqtab_df)
#----------------------------------------
#This function looks for the highest haplotype ID in the original order from dada2's output
maxn<-function(n) function(x) order(x,decreasing=TRUE)[n] #function that finds the highest value index of given columns
for(i in 1:8){
  index_top_hap = (apply(seqtab_df,1,maxn(i)))
  seqtab_df[,ncol(seqtab_df)+1]<-index_top_hap
  colnames(seqtab_df)[ncol(seqtab_df)]<-paste0("hap_id_top_",i)
}
#----------------------------------------
#View(seqtab_df) #sanity check

#----------------------------------------
#This for loop looks at the proportion of the haplotypes given the max number of reads
for (i in 1:columns_df){
  column = seqtab_df[,i] #gets the value of reads of the ith haplotype
  proportion_hap = column/sum_reads #gets the proportion of reads for the ith haplotype
  seqtab_df[,ncol(seqtab_df)+1]<-proportion_hap #creates a new column with the proportion
  colnames(seqtab_df)[ncol(seqtab_df)]<-paste0("haplotype_proportion_",i) #labels the new column
}
seqtab_df$total_reads<-(sum_reads)
View(seqtab_df)
seqtab_df

#Writing the output
out_file = paste0(out_path,"seqtab_nochim31_preprocessed_ama1d3.csv")
write.csv(seqtab_df,out_file)


table(nchar(getSequences(seqtab)))
