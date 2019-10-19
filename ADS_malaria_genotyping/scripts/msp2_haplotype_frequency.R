install.packages("AER")
data(Affairs, package="AER")
library(ggplot2)
library(dplyr)
library(grDevices)
ggplot2::ggplot(Affairs)
data.frame(Affairs)
msp2_data = as.matrix(read.table("/Users/Danniel/Desktop/msp_haplotype_frequency.csv",header=TRUE,sep=",",
                                 row.names=1,
                                 as.is=TRUE))
msp2_data
colnames=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10",
                      "S11","S12","S13","S14","S15","C1","C2","C3","C4","C5")
rownames=c("d046","824d","b98b","7523","f20b","5814","7da4","4991","7fab","0d0b","6bbd","fe03",
                      "2110","286d","9a3e","f684","f238","7159","a127","5ca0","b366","0e28","8526","37b9",
                      "ee0a","879d","a1b1","416c","9893","eb98","c1f0","5655","e8e3")

msp2_data

coul=c("coral3","yellow","cadetblue4","burlywood4","blueviolet","blue4","blanchedalmond","lightslateblue","chocolate4","bisque4","azure4","aquamarine4",
       "brown4","deeppink3","darkseagreen4","darksalmon","red","darkolivegreen1","gold",
       "blue","white","darkmagenta","steelblue","orange","black","green","lightblue","magenta","cyan")
barplot(msp2_data,col=coul,border="black",xlab="Samples",
        density=c(300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,
                  300,300,300,300,300,300,300,300,35,35,35,35),angle=c(45,90,180),ylab="Frequency",xlim=c(0,35),legend=TRUE,args.legend = list(x=35,y=1.1,ncol=2,text.width=2.5))

        