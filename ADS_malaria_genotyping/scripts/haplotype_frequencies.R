#This code was made to create the figures of msp MOI in regards of its frequency across the 
#20 samples that have been analyzed.

remove.packages(c("ggplot2", "data.table"))
install.packages('Rcpp', dependencies = TRUE)
install.packages('ggplot2', dependencies = TRUE)
install.packages('data.table', dependencies = TRUE)

install.packages("plotrix")
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}
install.packages("ggplot2")
library(ggplot2)
library(dplyr)
library(grDevices)
library(plotrix)
msp1_data = as.matrix(read.table("/Users/Danniel/Desktop/msp_corrections_study/msp_haplotype_frequency_patient_data.csv",header=TRUE,sep=",",
                                 row.names=1,
                                 as.is=TRUE))
msp1_data
colnames=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10",
                      "S11","S12","S13","S14","S15")
rownames=c("d046","824d","b98b","7523","f20b","5814","7da4","4991","7fab","0d0b","6bbd","fe03",
                      "2110","286d","9a3e","f684","f238","7159","a127","5ca0","b366","0e28","8526","37b9",
                      "ee0a","879d","a1b1","416c","9893","eb98","c1f0","5655","e8e3")

msp2_data = as.matrix(read.table("/Users/Danniel/Desktop/msp_corrections_study/msp_haplotype_frequency_msp2_patients.csv",header=TRUE,sep=",",
                                 row.names=1,
                                 as.is=TRUE))
colnames=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10",
          "S11","S12","S13","S14","S15")
rownames=c("d046","824d","b98b","7523","5814","7da4","4991","7fab","0d0b","6bbd","fe03",
           "2110","286d","9a3e","f684","f238","7159","a127","5ca0","b366","0e28","8526","37b9",
           "ee0a","879d","a1b1","416c","9893","eb98","c1f0","5655","e8e3")

msp2_data

coul=c("coral3","yellow","cadetblue4","burlywood4","blueviolet","blue4","blanchedalmond","lightslateblue","chocolate4","bisque4","azure4","aquamarine4",
       "brown4","deeppink3","darkseagreen4","darksalmon","red","darkolivegreen1","gold",
       "blue","white","darkmagenta","steelblue","orange","black","green","lightblue","cyan","magenta")
barplot(msp2_data,col=coul,border="black",xlab="Sample",
        density=c(300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,
                  300,300,300,300,300,300,300,300,35,35,35),angle=c(45),ylab="Frequency",xlim=c(0,35),main="msp2",legend=TRUE,args.legend = list(x=35,y=1.1,ncol=2,text.width=2.5))

msp2_data_control=as.matrix(read.table("/Users/Danniel/Desktop/msp_corrections_study/msp2_haplotype_frequency_control_data.csv",header=TRUE,sep=",",
                                       row.names=1,
                                       as.is=TRUE))
colnames=c("C16","C17","C18","C19","C20")
rownames=c("d046","824d","b98b","f20b","7523")

msp2_data_control
barplot(msp2_data_control,col=coul,border="black",xlab="Sample",
        density=c(300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,
                  300,300,300,300,300,300,300,300,35,35,35,35),angle=c(45,90,180),ylab="Frequency",xlim=c(0,15),main="msp2",legend=TRUE,args.legend = list(x=10,y=1,ncol=1,text.width=2))


moi_msp2 = matrix(c(1,2,3,7,12,1),ncol=2)
slices=c(7,12,1)
perncentage = round(slices/sum(slices)*100,digits = 4)
labels_pie=c("MOI = 1","MOI = 2","MOI = 3")
lbls=paste(labels_pie,"\n",perncentage)
lbls=paste(lbls,"  %",sep="")
pos=pie3D(slices,labels=lbls,explode=0.1,main="msp1",labelcex=0.7,theta=1.2,col=c("purple","violet","deepskyblue","white"),
      radius=0.9,labelpos=c(0.69,3,6.1))

msp2_shared = read.csv("/Users/Danniel/Desktop/shared_haplotypes_msp.csv")
coul=brewer.pal(4,"Pastel2")

p = ggplot(msp2_shared,aes(x=msp2_shared$Number.of.samples,y=msp2_shared$Haplotypes))+geom_bar(stat="identity",color="black",fill=rgb(0.1,0.4,0.5,0.7))+theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
p+xlab("Number of samples")+ylab("Haplotypes")+ggtitle("                                                                    msp2")+scale_y_discrete(limits=c(0,3,6,9,12,15,18,21,24,27,30))

#
#
#
#msp1 haplotyping now:

dev.off()
msp1_data = as.matrix(read.table("/Users/Danniel/Desktop/msp_corrections_study/msp_haplotype_frequency_patient_data.csv",header=TRUE,sep=",",
                                 row.names=1,
                                 as.is=TRUE))

colnames_msp1=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10",
           "S11","S12","S13","S14","S15")
rownames_msp1=c("7bb6","4f17","3a96","3624","05ed","c682","3e82","0f56","12b6","047c","6b73","3bbd","52de")

msp1_data

coul=c("coral3","yellow","cadetblue4","burlywood4","blueviolet","blue4","blue","orange","chocolate4","cyan","white","green",
       "red","deeppink3","darkseagreen4","darksalmon")
barplot(msp1_data,col=coul,border="black",xlab="Sample",
        density=c(300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,
                  300,300,300,300,300,300,300),angle=c(45,90,180),ylab="Frequency",xlim=c(0,29),main="msp1",legend=TRUE,args.legend = list(x=29,y=1.1,ncol=2,text.width=2.5))

msp1_controls = as.matrix(read.table("/Users/Danniel/Desktop/msp_corrections_study/msp1_haplotype_frequency_control_data.csv",header=TRUE,sep=",",
                                     row.names=1,
                                     as.is=TRUE))

#msp1_controls
colnames_msp1_control=c("C16","C17","C18","C19","C20")
rownames_msp1_control=c("2aa6","ff79","c00c","3a96","3624")
msp1_controls

coul2=c("deeppink3","darkseagreen4","darksalmon","cadetblue4","burlywood4")
barplot(msp1_controls,col=coul2,border="black",xlab="Sample",
        density=c(300,300,300,300,300),angle=c(45,90,180),ylab="Frequency",xlim=c(0,15),main="msp1",legend=TRUE,args.legend = list(x=10,y=1.0,ncol=1,text.width=2.0))


#pie chart of msp1 patient data
moi_msp2
slices=c(1,2,2)
perncentage = round(slices/sum(slices)*100,digits = 4)
labels_pie=c("COI = 1","COI = 2","COI = 3")
lbls=paste(labels_pie,"\n",perncentage)
lbls=paste(lbls,"  %",sep="")
pos=pie3D(slices,labels=lbls,explode=0.1,main="msp2",labelcex=0.7,theta=1.2,col=c("purple","violet","deepskyblue","white"),
          radius=0.9,labelpos=c(1,3,5.5,6))

#pie chart of msp1 control data
moi_msp1_controls = matrix(c(1,2,1,4),ncol=2)
moi_msp1_controls
slices=c(1,4)
perncentage = round(slices/sum(slices)*100,digits = 4)
labels_pie=c("MOI = 1","MOI = 2")
lbls=paste(labels_pie,"\n",perncentage)
lbls=paste(lbls,"  %",sep="")
pos=pie3D(slices,labels=lbls,explode=0.1,main="msp1",labelcex=0.7,theta=1.2,col=c("purple","violet","deepskyblue","white"),
          radius=0.9,labelpos=c(0.69,3,6.1))

#bar chart of shared moi among patients msp1 or msp2
msp1_shared = read.csv("/Users/Danniel/Desktop/msp_corrections_study/shared_haplotypes_msp2_patients.csv",sep = ',')
coul=brewer.pal(4,"Pastel2")
msp1_shared
p = ggplot(msp1_shared,aes(x=msp1_shared$Number.of.samples,y=msp1_shared$Haplotypes))+geom_bar(stat="identity",color="black",fill=rgb(0.1,0.4,0.5))+theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
p+xlab("Number of samples")+ylab("Haplotypes")+ggtitle("                                              msp2")+scale_y_discrete(limits=c(0,3,6,9,12,15,18,21,24,27,30))          


#bar chart of shared moi among controls of msp1 or msp2
msp1_shared = read.csv("/Users/Danniel/Desktop/msp_corrections_study/shared_haplotypes_msp2_controls.csv",sep = ',')
coul=brewer.pal(4,"Pastel2")
msp1_shared
p = ggplot(msp1_shared,aes(x=msp1_shared$Number.of.samples,y=msp1_shared$Haplotypes))+geom_bar(stat="identity",color="black",fill=rgb(0.1,0.4,0.5))+theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
p+xlab("Number of samples")+ylab("Haplotypes")+ggtitle("                                              msp2")+scale_y_discrete(limits=c(0,3,6,9,12,15,18,21,24,27,30))          




#agreement between msp1 and msp2 MOI
msp1_haplotypes=c(2,2,2,2,2,2,2,3,1,2,1,1,2,1,1,1,1,2,2,2)
msp2_haplotypes=c(2,2,2,2,2,2,3,7,2,1,2,1,1,3,2,2,1,2,3,3)
wilcox.test(msp1_haplotypes, conf.int = TRUE, conf.level = 0.95)
wilcox.test(msp2_haplotypes, conf.int = TRUE, conf.level = 0.95)


shapiro.test(msp1_haplotypes)
shapiro.test(msp2_haplotypes)
t.test(x=msp1_haplotypes,y=msp2_haplotypes,paired = FALSE)
wilcox.test(msp1_haplotypes,msp2_haplotypes,paired=TRUE)
msp1_control=c(1,1,2,2,2)
msp2_control=c(2,1,2,3,3)
expected_control=c(1,1,1,3,3)
wilcox.test(msp1_control,expected_control,paired=TRUE)
wilcox.test(msp2_control,expected_control,paired=TRUE)
t.test(msp1_control,expected_control,paired=TRUE)
t.test(msp2_control,expected_control,paired=TRUE)
        
