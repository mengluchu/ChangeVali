#library('rgdal')
#library('raster')

#load('prodespoints00.Rdata') # mask
#load('pdd.Rdata')           # reference data
#groundtruth<-pdd
#load('tssarar1.Rdata') # corrected st model with original array cusum
#load('tssarar2.Rdata')#mosum
#load('tssarar3.Rdata')#sar cusum
#load('tssarar4.Rdata')#sar mosum
#load('tssarar5.Rdata')#ar cusum
#load('tssarar6.Rdata')#ar mosum

#color2=rep(c('skyblue','red','gold','yellow','blue','antiquewhite2','darkolivegreen1','darkgreen','purple','pink','chartreuse1','deepskyblue4'),1000)

##150by 150 # coordinate of MODIS array
#x<-c(58930:59079)
#y<-c(48210:48359)
#load("edivisive1.Rdata")


######### generate confusion matrix: we use pontias' producer's accuracy ###########################
gcm<-function(changepoint, reference.sppoint=pdd,totalp=19167){
  #changepoint<-pvpoint

  modis.mt52<-changepoint
  deterpoints<-reference.sppoint
  result<-data.frame()
  deterpoints2<-c()


  b9<-length(modis.mt52)
  b10<-length(unique(modis.mt52@coords))/2

  b1<-length( deterpoints) #1898
  b2<-length(which(!is.na(over(deterpoints,modis.mt52)))) #997
  b3<-length(unique(deterpoints[is.na(over(deterpoints,modis.mt52))]@coords))/2

  b4<-length(unique(deterpoints[!is.na(over(deterpoints,modis.mt52))]@coords))/2
  b5<-length(unique(modis.mt52[is.na(over(modis.mt52,deterpoints))]@coords))/2 #5381 /4061
  b6<-length(unique(modis.mt52[!is.na(over(modis.mt52,deterpoints))]@coords))/2
  b7<-length(which(!is.na(over(modis.mt52,deterpoints)))) #1400/997
  b8<-length(which(is.na(over(modis.mt52,deterpoints))))

  TP = b7
  FN=  b3
  FP=  b8
  P=TP+FN
  N=totalp-P
  TN=N-FP
  ALL=totalp

  PRODUCER =TP/(TP+FP+FN)

  result<-cbind(TP,FN,FP,TN,PRODUCER)
  names(result)<-c("TP","FN","FP","TN","producer's accuracy")

  p1<-paste('total reference points:', b1,    'unique reference points not in result changepoints:',b3,
            'unique referece points in result changepoints',b4, 'unique reference points not in result changepoints', b5, 'unique result changepoints in reference points',b6,
            'unique total result changepoints',b10, sep='        ')

  p2<-paste( 'TRUE POSITIVE',TP ,
             'FALSE NEGATIVE',FN ,
             'FALSE POSITIVE', FP ,
             'PRODES POSITIVE',P ,
             'PRODES NEGATIVE',N ,
             'TRUE NEGATIVE',TN   ,
             'Precision', PRODUCER ,
             sep=', ')
  print(p1)
  print(p2)
  return(result)

}

### generate confusion matrix from result array, can use a mask to mask only interest area (sppoints) and set p-value#######

generatecmchange<-function( result.array,mask=prodespoints00,reference.sppoints=pdd,pv=0.05,x=c(58930:59079),y=c(48210:48359),crs=NULL,totalp=19167) #  put functions together and generate confusion matrix (filtered with a mask)
{
  pvpoint<- array2sp(changearray=result.array,x=x,y=y,crs=crs )
  if(!is.null(mask))
    pvpoint<-pvpoint[mask,]
  showre<-gcm(pvpoint,reference.sppoints,totalp)
  return(showre)
}

### generate confusion matrix from p-value array, can use a mask to mask only interest area (sppoints) and set p-value#######
generatecmpvalue<-function(result.array=pvaluemx ,reference.sppoints=pdd,mask=prodespoints00,pv=0.05,x=c(58930:59079),y=c(48210:48359),crs=NULL) #  put functions together and generate confusion matrix (filtered with a mask)
{

  pvpoint<-pvaluearray2sp(parray=result.array,x=x,y=y,1,1,pv,crs=crs)
  pvpoint<-pvpoint[mask,]
  showre<-gcm(pvpoint,reference.sppoints)
  return(showre)
}


generateppvalue<-function(result.array=pvaluemx ,reference.sppoints=pdd,mask=prodespoints00,pv=0.05,x=c(58930:59079),y=c(48210:48359),crs=NULL)
{
  pvpoint<-pvaluearray2sp(parray=result.array,x=x,y=y,xoff=1,yoff=1,pv=pv,crs=crs)
  if(!is.null(mask))
  pvpoint<-pvpoint[mask,]
  return(pvpoint)
}

generatepchange<-function( result.array,mask=prodespoints00,reference.sppoints=pdd,pv=0.05,x=c(58930:59079),y=c(48210:48359),crs=NULL) #  put functions together and generate confusion matrix (filtered with a mask)
{
  pvpoint<- array2sp(changearray=result.array,x=x,y=y,crs=crs )
  if(!is.null(mask))
    pvpoint<-pvpoint[mask,]

  return(pvpoint)
}

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
#ptssarar1<-generateppvalue(tssarar1,pv=0.05)
#ptssarar11<-generateppvalue(tssarar1,pv=0.005)
#ptssarar2<-generateppvalue(tssarar2,pv=0.05)
#ptssarar22<-generateppvalue(tssarar2,pv=0.025)
#ptssarar3<-generateppvalue(tssarar3,pv=0.05)
#ptssarar4<-generateppvalue(tssarar4,pv=0.05)
#ptssarar5<-generateppvalue(tssarar5,pv=0.05)
#ptssarar6<-generateppvalue(tssarar6,pv=0.05)

#ttssarar1<-generatecmpvalue(tssarar1,pdd,pv=0.05)
#ttssarar11<-generatecmpvalue(result.array=tssarar1 ,pv=0.2)
#ttssarar2<-generatecmpvalue(tssarar2,pdd,pv=0.05)
#ttssarar22<-generatecmpvalue(tssarar2,pdd,pv=0.1)
#ttssarar3<-generatecmpvalue(tssarar3,pdd,pv=0.05)
#ttssarar4<-generatecmpvalue(tssarar4,pdd,pv=0.05)
#ttssarar5<-generatecmpvalue(tssarar5,pdd,pv=0.05)
#ttssarar6<-generatecmpvalue(tssarar6,pdd,pv=0.05)

#cts<-rbind(ttssarar4,ttssarar6,ttssarar2,ttssarar22，ttssarar3，ttssarar5，ttssarar1，ttssarar11)

#jpeg('Pontius Producer Accuracy ts.jpg', height=6, width=10, res=400,unit="in")
#names1<-c("OLS-MOSUM p-value: 0.05","AR(1) OLS-MOSUM pvalue: 0.05","SAR OLS-MOSUM p-value: 0.05", "SAR OLS-MOSUM p-value: 0.1","OLS-CUSUM p-value: 0.05","AR(1) OLS-CUSUM p-value: 0.05","SAR OLS-CUSUM p-value: 0.05","SAR OLS-CUSUM p-value: 0.2")
#tex1='Confusion Matrix of Different Methods'

barplotcm<- function(cts,  names1=c("Ediv dst","ediv","BFAST"),n=8,title="Confusion Matrix")
{
par(mfrow=c(1,2))
barplot(cts[,1:4]/22500,beside=TRUE, main=title,
        legend.text=names1,

        cex.main=1,
        font.main=1,
        col=bpy.colors(n),
        #=rainbow(8,s = 1, v = 1, start = 0.05, end = 0.35),
        #legend.text=c("MOSUM","CUSUM","AR(1) MOSUM 0.05","AR(1) CUSUM 0.05","SAR MOSUM 0.05", "SAR MOSUM 0.025","SAR CUSUM 0.05","SAR CUSUM 0.005"),
        args.legend = list(x = "topleft", bty = "n",cex=0.8))

 barplot(cts[,5],beside=TRUE, main='Pontius Producer\'s Accuracy',
         col=bpy.colors(n),
         font.main=1,
         cex.main=1,
        #legend.text=c("MOSUM","CUSUM","SAR MOSUM 0.05", "SAR MOSUM 0.025","SAR CUSUM 0.05","SAR CUSUM 0.005"),
        #args.legend = list(x = "topleft", bty = "n",cex=0.6))
 )
}
#barplotcm(cts, n=8)
#generatemapRAS1(cpe,groundtruth,prodespoints00)
#generatemapRAS1(cpe2,groundtruth,prodespoints00)
#generatemapRAS1(cpb,groundtruth,prodespoints00)

generatemapRAS1<-function (pva, groundtruth=groundtruth,prodespoints00=prodespoints00,cols=c("#BEBADA","#FFFFB3",  "#8DD3C7", "#FB8072")) {
groundtruth<-spTransform(groundtruth,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))
prodespoints00<-spTransform(prodespoints00,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))
pva<-spTransform(pva,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))

FP<-pva[which(is.na(over(pva,groundtruth)))]
FN<-groundtruth[which(is.na(over(groundtruth,pva)))]
TP<-pva[groundtruth,]
PRODESNE<-prodespoints00[ which(is.na(over(prodespoints00,groundtruth))),]
BFASTNE<-prodespoints00[ which(is.na(over(prodespoints00,pva))),]
TN<-PRODESNE[BFASTNE,]

TPV<-rep(1,length(data.frame(TP)[,1]))
TNV<-rep(2,length(data.frame(TN)[,1]))
FPV<-rep(3,length(data.frame(FP)[,1]))
FNV<-rep(4,length(data.frame(FN)[,1]))
#plot(raster(TPdf))
TPdf<-SpatialPointsDataFrame( coordinates(data.frame(TP)),data=data.frame(TPV))
TNdf<-SpatialPointsDataFrame( coordinates(data.frame(TN)),data=data.frame(TNV))
FPdf<-SpatialPointsDataFrame( coordinates(data.frame(FP)),data=data.frame(FPV))
FNdf<-SpatialPointsDataFrame( coordinates(data.frame(FN)),data=data.frame(FNV))
a<-data.frame(TPdf)
b<-data.frame(FPdf)
c<-data.frame(FNdf)
d<-data.frame(TNdf)

names(c)<-names(a)
names(b)<-names(a)
names(d)<-names(a)
abcd<-rbind(a,b,c,d)

coordinates(abcd)<-~x+y
gridded(abcd)<-TRUE

r1<-raster(abcd)

raty = ratify(r1)
rat <- levels(raty)[[1]]
rat$class <- c('True Positive', 'True Negative', 'False Positive', 'False Negative')
levels(raty) <- rat
#cols = c("#ffff99", "#7fc97f", "#1f78b4", "#fdc086")
 
#c("#BEBADA","#FFFFB3",  "#8DD3C7", "#FB8072") #the sequence is reversed

ty<-levelplot(raty,col.regions =cols,names.attr=c("edivisive"))

plot(ty)
}
 

#take 1 object for comparison, plot the raster
generateRAS<-function (pva , groundtruth=groundtruth,prodespoints00=prodespoints00 ) {
  groundtruth<-spTransform(groundtruth,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))
  prodespoints00<-spTransform(prodespoints00,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))
  pva<-spTransform(pva,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))

  FP<-pva[which(is.na(over(pva,groundtruth)))]
  FN<-groundtruth[which(is.na(over(groundtruth,pva)))]
  TP<-pva[groundtruth,]
  PRODESNE<-prodespoints00[ which(is.na(over(prodespoints00,groundtruth))),]
  BFASTNE<-prodespoints00[ which(is.na(over(prodespoints00,pva))),]
  TN<-PRODESNE[BFASTNE,]

  TPV<-rep(1,length(data.frame(TP)[,1]))
  TNV<-rep(2,length(data.frame(TN)[,1]))
  FPV<-rep(3,length(data.frame(FP)[,1]))
  FNV<-rep(4,length(data.frame(FN)[,1]))
  #plot(raster(TPdf))
  TPdf<-SpatialPointsDataFrame( coordinates(data.frame(TP)),data=data.frame(TPV))
  TNdf<-SpatialPointsDataFrame( coordinates(data.frame(TN)),data=data.frame(TNV))
  FPdf<-SpatialPointsDataFrame( coordinates(data.frame(FP)),data=data.frame(FPV))
  FNdf<-SpatialPointsDataFrame( coordinates(data.frame(FN)),data=data.frame(FNV))
  a<-data.frame(TPdf)
  b<-data.frame(FPdf)
  c<-data.frame(FNdf)
  d<-data.frame(TNdf)

  names(c)<-names(a)
  names(b)<-names(a)
  names(d)<-names(a)
  abcd<-rbind(a,b,c,d)

  coordinates(abcd)<-~x+y
  gridded(abcd)<-TRUE

  r1<-raster(abcd)

  raty = ratify(r1)
  rat <- levels(raty)[[1]]
  rat$class <- c('True Positive', 'True Negative', 'False Positive', 'False Negative')
  levels(raty) <- rat
  return(raty)

}
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
