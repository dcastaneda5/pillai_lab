#Boxplot for msp1 and msp2 markers (swarm)
#against srst2 and pcr
library("ggpubr")
library(ggsignif)
moi_msp1_swarm = c(2,2,2,2,2,2,2,3,1,2,1,1,2,1,1,1,1,2,2,2)
moi_msp2_swarm = c(2,2,2,2,2,2,3,7,2,1,2,1,1,3,2,2,1,2,3,3)
moi_srst2 = c(2,2,2,2,2,3,3,3,2,3,2,2,2,2,0,1,2,2,3,3)
moi_pcr = c(1,1,1,1,1,3,3,0,1,2,1,2,1,1,1,1,1,1,3,2)

shapiro.test(moi_msp1_swarm)
shapiro.test(moi_msp2_swarm)
shapiro.test(moi_srst2)
shapiro.test(moi_pcr)

#NONE OF MY DATA HAS A NORMAL DISTRIBUTION
length(moi_msp1_swarm)
length(moi_msp2_swarm)
length(moi_srst2)
length(moi_pcr)

mean1=mean(moi_msp1_swarm)
mean2=mean(moi_msp2_swarm)
mean3=mean(moi_srst2)
mean4=mean(moi_pcr)
sd1=sd(moi_msp1_swarm)
sd2=sd(moi_msp2_swarm)
sd3=sd(moi_srst2)
sd4=sd(moi_pcr)
se1=sd1/(sqrt(20))
se2=sd2/(sqrt(20))
se3=sd3/(sqrt(20))
se4=sd4/(sqrt(20))

MOI=c(mean1,mean2,mean3,mean4)
sds=c(sd1,sd2,sd3,sd4)
ses=c(se1,se2,se3,se4)

my_comparisons=list(c("msp1 swarm","msp2 swarm"),c("msp2 swarm","SRST2"))
matrix_analysis = matrix(c("msp1 swarm","msp2 swarm","SRST2","Nested PCR",mean1,mean2,mean3,mean4),ncol=2)
matrix_analysis
df_msp = as.data.frame(matrix_analysis)
Method=c("msp1 swarm","msp2 swarm","SRST2","Nested PCR")
ggplot(df_msp) + 
  geom_bar( aes(x=Method, y=MOI), stat='identity', fill="skyblue") +
  geom_errorbar( aes(x=Method,ymin=MOI-ses,ymax=MOI+ses),width=0.1) +
  theme_bw() + ylim(0,2.9)
  
 
wilcox.test(x=moi_msp1_swarm,y = moi_msp2_swarm,paired = TRUE)
wilcox.test(x=moi_msp1_swarm,y = moi_pcr, paired=TRUE)
wilcox.test(x=moi_msp1_swarm,y=moi_srst2,paired=TRUE)
wilcox.test(x=moi_msp2_swarm,y=moi_pcr,paired=TRUE)
wilcox.test(x=moi_msp2_swarm,y=moi_srst2,paired=TRUE)
wilcox.test(x=moi_pcr,y=moi_srst2)

t.test(x=moi_msp1_swarm,y = moi_msp2_swarm)
t.test(x=moi_msp1_swarm,y = moi_pcr)
t.test(x=moi_msp1_swarm,y=moi_srst2)
t.test(x=moi_msp2_swarm,y=moi_pcr)
t.test(x=moi_msp2_swarm,y=moi_srst2)
t.test(x=moi_pcr,y=moi_srst2)

moi_pcr=c(1,1,1,1,1,3,3,0,1,2,1,2,1,1,1,1,1,1,3,2)
moi_srs=c(1,1,2,1,2,3,3,3,1,3,2,2,2,1,0,1,2,2,3,2)
mean(moi_pcr)
mean(moi_srs)
sd_pcr=sd(moi_pcr)
sd_srs=sd(moi_srs)
sd_pcr
sd_srs
se_pcr = sd_pcr/(sqrt(20))
se_srs = sd_srs/(sqrt(20))
se_pcr
se_srs

moi_msp1 = c(1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,3)
moi_msp2 = c(1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,7)
mean(moi_msp1)
sd_msp1 =sd(moi_msp1)
se_msp1 = sd_msp1/(sqrt(20))
se_msp1
mean(moi_msp2)
sd_msp2 = sd(moi_msp2)
sd_msp2
se_msp2 = sd_msp2/(sqrt(20))
se_msp2
