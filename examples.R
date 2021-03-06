#library('rgdal')
library('rasterVis')

#run the validation tools, by default using the 150 by 150 validation data and filter
MODISCRS <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
UTM21S <- "+proj=utm +zone=21 +south"

data(prodespoints00) # mask
data(pdd)           # reference data
groundtruth<-pdd
data(tssarar2) # corrected st model with original array cusum
data(tssarar3)#sar cusum
#load('tssarar4.Rdata')#sar mosum
#load('tssarar5.Rdata')#ar cusum
#load('tssarar6.Rdata')#ar mosum

#####################################
#example: validate bfast results
#load("C:\\Users\\m_lu0002\\Dropbox\\mengluchu\\bfast\\a150p05t.Rdata")
#data(prodespoints00) # mask
#data(pdd)           # reference data
#groundtruth<-pdd #
#c1<-generatepchange(a150p05t)
#generatemapRAS1(c1,groundtruth,prodespoints00)
#color2=rep(c('skyblue','red','gold','yellow','blue','antiquewhite2','darkolivegreen1','darkgreen','purple','pink','chartreuse1','deepskyblue4'),1000)
##########################



#load("C:/Users/m_lu0002/Dropbox/mengluchu/bfast/SpatialPolygons.R")
#load("C:/Users/m_lu0002/Dropbox/mengluchu/bfast/tdall2.Rdata")
#load("C:/Users/m_lu0002/Dropbox/mengluchu/bfast/pp0.Rdata")
#load("C:/Users/m_lu0002/Dropbox/mengluchu/bfast/rpp.Rdata")
#load("C:/Users/m_lu0002/Dropbox/mengluchu/bfast/xyc.Rdata")  # whole 22500 points

#load("C:/Users/m_lu0002/Dropbox/mengluchu/bfast2/bt.Rdata")
#load("C:/Users/m_lu0002/Dropbox/mengluchu/bfast2/t3darrbfamul2.Rdata")
#load("C:/Users/m_lu0002/Dropbox/mengluchu/bfast/edivisive1.Rdata")
#load("C:/Users/m_lu0002/Dropbox/mengluchu/bfast/edivisive2.Rdata")
#load("C:/Users/m_lu0002/Dropbox/mengluchu/bfast/bfast1.Rdata")


#SpatialPolygon <- spatialPolygons
#tdall <- tdall2
############################
##150by 150 # coordinate of MODIS array
x<-c(58930:59079)
y<-c(48210:48359)
#load("edivisive1.Rdata")


######### generate confusion matrix: we use pontias' producer's accuracy ###########################

### generate confusion matrix from result array, can use a mask to mask only interest area (sppoints) and set p-value#######

### generate confusion matrix from p-value array, can use a mask to mask only interest area (sppoints) and set p-value#######

#load('prodespoints00.Rdata') # mask
#load('pdd.Rdata')
#groundtruth<-pdd
#cpe2<- generatepchange(edivisive2,mask=prodespoints00,reference.sppoints=pdd,pv=0.05,x=c(58930:59079),y=c(48210:48359),crs= CRS("+proj=utm +zone=21 +south"))
#cpe<- generatepchange(edivisive1,mask=prodespoints00,reference.sppoints=pdd,pv=0.05,x=c(58930:59079),y=c(48210:48359),crs= CRS("+proj=utm +zone=21 +south"))
#cpb<-generatepchange(bfast1,mask=prodespoints00,reference.sppoints=pdd,pv=0.05,x=c(58930:59079),y=c(48210:48359),crs= CRS("+proj=utm +zone=21 +south"))
#
#cpe2cm<-generatecmchange(result.array=edivisive2 ,reference.sppoints=pdd,mask=prodespoints00,pv=0.05,x=c(58930:59079),y=c(48210:48359))
#cpe1cm<-generatecmchange(result.array=edivisive1 ,reference.sppoints=pdd,mask=prodespoints00,pv=0.05,x=c(58930:59079),y=c(48210:48359))
#cpbcm<-generatecmchange(result.array=bfast1 ,reference.sppoints=pdd,mask=prodespoints00,pv=0.05,x=c(58930:59079),y=c(48210:48359))
#cts<-rbind(cpe2cm, cpe1cm,cpbcm)
#
### do things
e1[is.na(e1)]<-1
summary(l1)
ptssarar1<-generateppvalue(result.array=e1,pv=0.05)
ptssarar11<-generateppvalue(result.array=l1,pv=0.005)
#ptssarar1<-generateppvalue(tssarar1,pv=0.05)
#ptssarar11<-generateppvalue(tssarar1,pv=0.005)
ptssarar2<-generateppvalue(tssarar2,pv=0.05)
ptssarar22<-generateppvalue(tssarar2,pv=0.025)
#ptssarar3<-generateppvalue(tssarar3,pv=0.05)
#ptssarar4<-generateppvalue(tssarar4,pv=0.05)
#ptssarar5<-generateppvalue(tssarar5,pv=0.05)
#ptssarar6<-generateppvalue(tssarar6,pv=0.05)

#ttssarar1<-generatecmpvalue(tssarar1,pdd,pv=0.05)
#ttssarar11<-generatecmpvalue(result.array=tssarar1 ,pv=0.2)
ttssarar2<-generatecmpvalue(tssarar2,pdd,pv=0.05)
#ttssarar22<-generatecmpvalue(tssarar2,pdd,pv=0.1)
#ttssarar3<-generatecmpvalue(tssarar3,pdd,pv=0.05)
#ttssarar4<-generatecmpvalue(tssarar4,pdd,pv=0.05)
#ttssarar5<-generatecmpvalue(tssarar5,pdd,pv=0.05)
#ttssarar6<-generatecmpvalue(tssarar6,pdd,pv=0.05)

#cts<-rbind(ttssarar4,ttssarar6,ttssarar2,ttssarar22，ttssarar3，ttssarar5，ttssarar1，ttssarar11)

#jpeg('Pontius Producer Accuracy ts.jpg', height=6, width=10, res=400,unit="in")
#names1<-c("OLS-MOSUM p-value: 0.05","AR(1) OLS-MOSUM pvalue: 0.05","SAR OLS-MOSUM p-value: 0.05", "SAR OLS-MOSUM p-value: 0.1","OLS-CUSUM p-value: 0.05","AR(1) OLS-CUSUM p-value: 0.05","SAR OLS-CUSUM p-value: 0.05","SAR OLS-CUSUM p-value: 0.2")
#tex1='Confusion Matrix of Different Methods'

#barplotcm(cts, n=8)
generatemapRAS1(ptssarar22,groundtruth,prodespoints00)
generatemapRAS1(ptssarar1,groundtruth,prodespoints00)

#generatemapRAS1(cpe2,groundtruth,prodespoints00)
#generatemapRAS1(cpb,groundtruth,prodespoints00)

#take 1 object for comparison, plot the raster
#return the object
#r1<-generateRAS(pva=cpe,groundtruth,prodespoints00)
#rb<-generateRAS(pva=cpb,groundtruth,prodespoints00)
#r2<-generateRAS(pva=cpe2,groundtruth,prodespoints00)
#rasterStack1 = stack(r1,r2, rb)
#names(rasterStack1) = c("e.divisive","e.divisive dst", "BFAST" )

#l<- levelplot(rasterStack1
#              , pretty=TRUE, margin=F
#              , xlab="Latitude (Sinusoidal Projection)", ylab="Longitude (Sinusoidal Projection)",
#              col.regions=cols,names.attr=c("e.divisive(monthly)", "e.divisive dst(monthly)", "BFAST (monthly)" ))
#plot(l)
#take all four
