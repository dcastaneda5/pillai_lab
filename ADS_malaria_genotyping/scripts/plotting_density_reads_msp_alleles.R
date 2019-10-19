#plotting allele curves vs number of reads by using SRST2 to call specific genes
#Daniel Castaneda Mogollon

library("coda")
require("logspline")

file_srst2 = '/Volumes/Seagate/msp_project/files/filtered/srst2_trimmed_output__fullgenes__srst2_whole_db__results.xls'
df_srst2 = readxl::read_xls(path = file_srst2, sheet=2)
print(df_srst2)
class(df_srst2)
k1_data = subset(x = df_srst2, allele=="K1",select=c("allele","depth"))
k1_data
mad20_data = subset(x = df_srst2, allele=="MAD20",select=c("allele","depth"))
mad20_data
ro33_data = subset(x = df_srst2, allele=="RO33",select=c("allele","depth"))
fc27_data = subset(x = df_srst2, allele=="FC27",select=c("allele","depth"))
d7_data = subset(x = df_srst2, allele=="3D7",select=c("allele","depth"))
#k1_density2 = oldlogspline(uncensored = k1_data$allele, lbound = 0)
k1_density = density(x = k1_data$depth, kernel = "g", from = 0, to = 10000)
#mad20_density2 = logspline(mad20_data)
mad20_density = density(x = mad20_data$depth, kernel = "g", from = 0, to = 10000)
ro33_density = density(x=ro33_data$depth, kernel = "g", from = 0, to = 10000)
fc27_density = density(x=fc27_data$depth, kernel = 'g', from = 0, to = 10000)
#d7_density2 = density(x=d7_data$depth, kernel = 'g')
d7_density = density(x=d7_data$depth, kernel = 'g', from=0,to=10000)
plot(x = d7_density,ylab = "Density", xlab = "Mapped DNA reads",main="") 
#hist(d7_data$depth,breaks = 50,probability = TRUE)
#lines(density(d7_data$depth, kernel='g', from=0, to=10000),col='red')
lines(x= mad20_density,col="blue")
lines(x = ro33_density, col="red")
lines(x= fc27_density,col="purple")
lines(x= k1_density, col='green')
axis(2,at=y,labels=scientific_10x(y))
legend(x="topright",y=0.92, legend=c("3D7","MAD20","RO33","FC27","K1"),
       col=c("black","blue","red","purple","green"),cex = 0.8,
       pch=c("-","-", "-", "-","-"))

#all reads density
all_alleles = density(x = df_srst2$depth, kernel='g', from = 0, to = 10000)


#getting the 0.05 lower bound
dn_k1 = cumsum(k1_density$y/sum(k1_density$y))
li_k1 = which(dn_k1>=0.05)[1]
k1_density$x[c(li_k1)]

dn_mad20 = cumsum(mad20_density$y/sum(mad20_density$y))
li_mad20 = which(dn_mad20>=0.05)[1]
mad20_density$x[c(li_mad20)]

dn_ro33 = cumsum(ro33_density$y/sum(ro33_density$y))
li_ro33 = which(dn_ro33>=0.05)[1]
ro33_density$x[c(li_ro33)]

dn_d7 = cumsum(d7_density$y/sum(d7_density$y))
li_d7 = which(dn_d7>=0.05)[1]
d7_density$x[c(li_d7)]

dn_fc27 = cumsum(fc27_density$y/sum(fc27_density$y))
li_fc27 = which(dn_fc27>=0.05)[1]
fc27_density$x[c(li_fc27)]

dn_all = cumsum(all_alleles$y/sum(all_alleles$y))
li_all = which(dn_all>=0.1)[1]
all_alleles$x[c(li_all)]

#dn_d72 = cumsum(d7_density2$y/sum(d7_density2$y))
#li_d72 = which(dn_d72>=0.05)[1]
#d7_density2$x[c(li_d72)]

plot(all_alleles, ylab='Density', xlab='DNA read size', main='')




