library(boot) 

# function to obtain the mean
Bmean <- function(data, indices) {
  d <- data[indices] # allows boot to select sample 
  return(mean(d))
} 

# bootstrapping with 1000 replications 
results <- boot(data=data0, statistic=Bmean, R=1000)

# view results
results 
plot(results)

# get 95% confidence interval 
boot.ci(results, type=c("norm", "basic", "perc", "bca"))

#function for 95% confidence interval
conf_interval <- function(msp_vector){
  msp_results = boot(data=msp_vector, statistic=Bmean, R=1000)
  final=boot.ci(msp_results,type="norm")
  print(final)
}

#dataset for DSCVH results for msp2
msp2_d1_5=c(2,2,2,1,2,2,2,2,1,1,2,1,1,2,2)
msp2_d1_1=c(2,3,4,1,3,3,3,4,2,2,2,3,2,3,5)
msp2_d1_f_5=c(2,2,2,1,2,2,2,2,1,1,2,1,1,2,2)
msp2_d1_f_1=c(2,4,4,1,3,3,3,4,2,2,2,3,2,3,5)
msp2_d2_5=c(2,2,2,2,2,2,3,7,2,1,2,1,1,3,2)
msp2_d2_1=c(3,5,5,2,3,6,7,16,3,3,2,4,2,4,5)
msp2_d3_5=c(2,2,2,2,2,3,3,12,2,1,2,2,1,3,2)
msp2_d3_1=c(3,5,5,2,3,7,8,22,4,3,2,6,2,4,5)
msp2_control_d1_5=c(2,1,2,3,3)
msp2_control_d1_1=c(3,1,4,3,5)
msp2_control_d1_f_5=c(2,1,2,3,3)
msp2_control_d1_f_1=c(3,1,4,3,5)
msp2_control_d2_5=c(2,1,2,3,3)
msp2_control_d2_1=c(7,1,5,4,6)
msp2_control_d3_5=c(2,1,2,3,3)
msp2_control_d3_1=c(4,1,5,4,6)

#running ci intervals for DSCVH results for msp2
conf_interval(msp2_d1_5)
conf_interval(msp2_d1_1)
conf_interval(msp2_d1_f_5)
conf_interval(msp2_d1_f_1)
conf_interval(msp2_d2_5)
conf_interval(msp2_d2_1)
conf_interval(msp2_d3_5)
conf_interval(msp2_d3_1)
#
conf_interval(msp2_control_d1_5)
conf_interval(msp2_control_d1_1)
conf_interval(msp2_control_d1_f_5)
conf_interval(msp2_control_d1_f_1)
conf_interval(msp2_control_d2_5)
conf_interval(msp2_control_d2_1)
conf_interval(msp2_control_d3_5)
conf_interval(msp2_control_d3_1)

#dataset for DSCVH results for msp1
msp1_d1_5=c(1,1,1,1,1,1,1,2,1,1,1,1,1,1,1)
msp1_d1_1=c(2,2,2,4,2,4,5,2,7,2,4,3,1,2,1)
msp1_d1_f_5=c(1,1,1,1,1,1,1,2,1,1,1,1,1,1,1)
msp1_d1_f_1=c(2,2,2,4,3,4,5,2,8,2,4,3,1,2,38)
msp1_d2_5=c(2,2,2,2,2,2,2,3,1,2,1,1,2,1,1)
msp1_d2_1=c(3,3,3,5,5,5,6,3,12,3,5,4,2,4,39)
msp1_d3_5=c(2,2,2,2,2,2,2,3,1,2,1,1,1,1,1)
msp1_d3_1=c(3,5,4,6,6,5,7,4,9,4,5,4,2,4,47)
msp1_control_d1_5=c(1,1,1,1,1)
msp1_control_d1_1=c(2,5,2,2,2)
msp1_control_d1_5_f=c(1,1,1,1,1)
msp1_control_d1_1_f=c(2,5,2,2,2)
msp1_control_d2_5=c(1,1,2,2,2)
msp1_control_d2_1=c(2,6,3,3,4)
msp1_control_d3_5=c(1,2,2,2,2)
msp1_control_d3_1=c(2,6,4,4,4)

#running ci intervals for DSCVH results for msp1
conf_interval(msp1_d1_5)
conf_interval(msp1_d1_1)
conf_interval(msp1_d1_f_5)
conf_interval(msp1_d1_f_1)
conf_interval(msp1_d2_5)
conf_interval(msp1_d2_1)
conf_interval(msp1_d3_5)
conf_interval(msp1_d3_1)
#
conf_interval(msp1_control_d1_5)
conf_interval(msp1_control_d1_1)
conf_interval(msp1_control_d1_5_f)
conf_interval(msp1_control_d1_1_f)
conf_interval(msp1_control_d2_5)
conf_interval(msp1_control_d2_1)
conf_interval(msp1_control_d3_5)
conf_interval(msp1_control_d3_1)

#nested_pcr and dsch
msp1_nested=c(0,1,1,1,1,3,3,0,1,2,1,2,1,1,1)
msp2_nested=c(1,1,1,1,1,1,2,0,0,1,1,1,1,0,1)
dsch_msp1=c(2,3,2,1,2,3,3,3,3,3,2,2,1,2,2)
dsch_msp2=c(2,2,2,2,2,2,2,2,2,2,1,2,2,2,2)
coi_nested=c(1,1,1,1,1,3,3,0,1,2,1,2,1,1,1)
coi_dsch=c(2,3,2,2,2,3,3,3,3,3,2,2,2,2,2)
conf_interval(msp1_nested)
conf_interval(msp2_nested)
conf_interval(dsch_msp1)
conf_interval(dsch_msp2)
conf_interval(coi_nested)
conf_interval(coi_dsch)

#nested_pcr and dsch controls
msp1_nested_c=c(1,1,1,3,2)
msp2_nested_c=c(1,1,1,2,1)
msp1_dsch_c=c(2,1,1,3,3)
msp2_dsch_c=c(2,2,2,2,2)
coi_nested_c=c(1,1,1,3,2)
coi_dsch_c=c(2,2,2,3,3)
conf_interval(msp1_nested_c)
conf_interval(msp2_nested_c)
conf_interval(msp1_dsch_c)
conf_interval(msp2_dsch_c)
conf_interval(coi_nested_c)
conf_interval(coi_dsch_c)
