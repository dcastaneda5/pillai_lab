#Boxplot for msp1 and msp2 markers (swarm)
#against DSCH and pcr
library("ggpubr")
library(ggsignif)

dev.off()
#defining vector variables:
nested_pcr_patients=c(1,1,1,1,1,3,3,0,1,2,1,2,1,1,1)
dsch_patients=c(2,3,2,2,2,3,3,3,3,3,2,2,2,2,2)
dscvh_msp1_patients=c(2,2,2,2,2,2,2,3,1,2,1,1,2,1,1)
dscvh_msp2_patients=c(2,2,2,2,2,2,3,7,2,1,2,1,1,3,2)
dshr_msp1_patients=c(3,3,4,6,2,5,8,5,3,6,4,3,2,3,2)
dshr_msp2_patients=c(5,2,6,1,8,4,5,9,2,9,3,2,1,1,6)

#confirming length of 15 across variables
length(nested_pcr_patients)
length(dsch_patients)
length(dscvh_msp1_patients)
length(dscvh_msp2_patients)
length(dshr_msp1_patients)
length(dshr_msp2_patients)

#looking at its distribution (normal or not-normal)
shapiro.test(nested_pcr_patients)
shapiro.test(dsch_patients)
shapiro.test(dscvh_msp1_patients)
shapiro.test(dscvh_msp2_patients)
shapiro.test(dshr_msp1_patients) #normally distributed
shapiro.test(dshr_msp2_patients) #normally distributed

mean1=mean(dshr_msp1_patients)
mean2=mean(dshr_msp2_patients)
mean3=mean(dscvh_msp1_patients)
mean4=mean(dscvh_msp2_patients)
mean5=mean(dsch_patients)
mean6=mean(nested_pcr_patients)

sd1=sd(dshr_msp1_patients)
sd2=sd(dshr_msp2_patients)
sd3=sd(dscvh_msp1_patients)
sd4=sd(dscvh_msp2_patients)
sd5=sd(dsch_patients)
sd6=sd(nested_pcr_patients)

se1=sd1/(sqrt(15))
se2=sd2/(sqrt(15))
se3=sd3/(sqrt(15))
se4=sd4/(sqrt(15))
se5=sd5/(sqrt(15))
se6=sd6/(sqrt(15))

#grouping and printing variables to see the order
COI=c(mean1,mean2,mean3,mean4,mean5,mean6)
COI
sds=c(sd1,sd2,sd3,sd4,sd5,sd6)
sds
ses=c(se1,se2,se3,se4,se5,se6)
ses

matrix_analysis = matrix(c(" DSHR msp1"," DSHR msp2","DSCVH msp1","DSCVH msp2","DSCH","Nested PCR",mean1,mean2,mean3,mean4,mean5,mean6),ncol=2)
matrix_analysis
df_msp = as.data.frame(matrix_analysis)
Method=c(" DSHR msp1"," DSHR msp2","DSCVH msp1","DSCVH msp2","DSCH","Nested PCR")
ggplot(df_msp) + 
  geom_bar( aes(x=Method, y=COI), stat='identity', fill="skyblue") +
  geom_errorbar( aes(x=Method,ymin=COI-ses,ymax=COI+ses),width=0.1) +
  theme_bw() + ylim(0,6)
  

#------------------------------------------------------------------------------#
#Determining if msp1 and msp2 are significantly different in DSCH and nested PCR
DSCH_msp1 = c(2,3,2,1,2,3,3,3,3,3,2,2,2,1,2)
DSCH_msp2 = c(2,2,2,2,2,2,2,2,2,2,1,2,2,2,2)
wilcox.test(x = DSCH_msp1,y = DSCH_msp2,paired = FALSE)
#NON SIGNIFICANT
pcr_msp1 = c(0,1,1,1,1,3,3,0,1,2,1,2,1,1,1)
pcr_msp2 = c(1,1,1,1,1,1,2,0,0,1,1,1,1,0,1)
wilcox.test(x=pcr_msp1,y=pcr_msp2,paired=FALSE)
#NON SIGNIFICANT
t.test(x = dshr_msp1_patients,y = dshr_msp2_patients,paired = FALSE,var.equal = FALSE)

my_comparisons_DSCH = c("msp1","msp2")
mean_DSCH_msp1 = mean(DSCH_msp1)
mean_DSCH_msp2 = mean(DSCH_msp2)
sd_DSCH_msp1 = sd(DSCH_msp1)
sd_DSCH_msp2 = sd(DSCH_msp2)
se_DSCH_msp1 = sd_DSCH_msp1/(sqrt(20))
se_DSCH_msp2 = sd_DSCH_msp2/(sqrt(20))

MOI_DSCH_msp1_msp2 = c(mean_DSCH_msp1,mean_DSCH_msp2)
sds_DSCH = c(sd_DSCH_msp1,sd_DSCH_msp2)
ses_DSCH = c(se_DSCH_msp1,se_DSCH_msp2)
matrix_analysis_DSCH = matrix(c("msp1","msp2",mean_DSCH_msp1,mean_DSCH_msp2),ncol=2)
matrix_analysis_DSCH
df_DSCH = as.data.frame(matrix_analysis_DSCH)
Method=c("msp1","msp2")
ggplot(df_DSCH) + 
  geom_bar( aes(x=Method, y=MOI_DSCH_msp1_msp2), stat='identity', fill="skyblue") +
  geom_errorbar( aes(x=Method,ymin=MOI_DSCH_msp1_msp2-ses_DSCH,ymax=MOI_DSCH_msp1_msp2+ses_DSCH),width=0.1) +
  theme_bw() + ylim(0,2.9) + xlab("DSCH gene") + ylab("COI")

my_comparisons_pcr = c("msp1","msp2")
mean_pcr_msp1 = mean(pcr_msp1)
mean_pcr_msp2 = mean(pcr_msp2)
sd_pcr_msp1 = sd(pcr_msp1)
sd_pcr_msp2 = sd(pcr_msp2)
se_pcr_msp1 = sd_pcr_msp1/(sqrt(20))
se_pcr_msp2 = sd_pcr_msp2/(sqrt(20))

MOI_pcr_msp1_msp2 = c(mean_pcr_msp1,mean_pcr_msp2)
sds_pcr = c(sd_pcr_msp1,sd_pcr_msp2)
ses_pcr = c(se_pcr_msp1,se_pcr_msp2)
matrix_analysis_pcr = matrix(c("msp1","msp2",mean_pcr_msp1,mean_pcr_msp2),ncol=2)
matrix_analysis_pcr
df_pcr = as.data.frame(matrix_analysis_pcr)
ggplot(df_pcr) + 
  geom_bar( aes(x=Method, y=MOI_pcr_msp1_msp2), stat='identity', fill="skyblue") +
  geom_errorbar( aes(x=Method,ymin=MOI_pcr_msp1_msp2-ses_pcr,ymax=MOI_pcr_msp1_msp2+ses_pcr),width=0.1) +
  theme_bw() + ylim(0,2.9) + xlab("Nested PCR gene") + ylab("COI")

my_comparisons_dshr=c("msp1","msp2")
MOI_dshr_msp1_msp2 = c(mean1,mean2)
sds_dshr = c(sd1,sd2)
ses_dshr = c(se1,se2)
matrix_analysis_dshr = matrix(c("msp1","msp2",mean1,mean2),ncol=2)
matrix_analysis_dshr
df_dshr = as.data.frame(matrix_analysis_dshr)
ggplot(df_dshr) + 
  geom_bar( aes(x=Method, y=MOI_dshr_msp1_msp2), stat='identity', fill="skyblue") +
  geom_errorbar( aes(x=Method,ymin=MOI_dshr_msp1_msp2-ses_dshr,ymax=MOI_dshr_msp1_msp2+ses_dshr),width=0.1) +
  theme_bw() + ylim(0,6) + xlab("DSHR gene") + ylab("COI")

#test combinations
wilcox.test(x=nested_pcr_patients,y=dsch_patients,paired=FALSE)
wilcox.test(x=nested_pcr_patients,y=dscvh_msp1_patients,paired=FALSE)
wilcox.test(x=nested_pcr_patients,y=dscvh_msp2_patients,paired=FALSE)
wilcox.test(x=nested_pcr_patients,y=dshr_msp1_patients,paired=FALSE)
wilcox.test(x=nested_pcr_patients,y=dshr_msp2_patients,paired=FALSE)
wilcox.test(x=dsch_patients,y=dscvh_msp1_patients,paired=FALSE)
wilcox.test(x=dsch_patients,y=dscvh_msp2_patients,paired=FALSE)
wilcox.test(x=dsch_patients,y=dshr_msp1_patients,paired=FALSE)
wilcox.test(x=dsch_patients,y=dshr_msp2_patients,paired=FALSE)
wilcox.test(x=dscvh_msp1_patients,y=dscvh_msp2_patients,paired=FALSE)
wilcox.test(x=dscvh_msp1_patients,y=dshr_msp1_patients,paired=FALSE)
wilcox.test(x=dscvh_msp1_patients,y=dshr_msp2_patients,paired=FALSE)
wilcox.test(x=dscvh_msp2_patients,y=dshr_msp1_patients,paired=FALSE)
wilcox.test(x=dscvh_msp2_patients,y=dshr_msp2_patients,paired=FALSE)

dscvh_msp2_patients
dshr_msp2_patients
wilcox.test(x=dshr_msp2_patients,y=dscvh_msp2_patients,paired=FALSE)
