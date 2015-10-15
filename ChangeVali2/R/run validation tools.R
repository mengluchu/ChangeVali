#example: TP of MOSUM and SAR MOSUM 
#cmsarm<-tfpn(psarmosum,reference.sppoint = pmosum) #TP, FN, FP, TN
#TPts1<-stfdfevi8c[cmsarm[[1]],]# sar mosum and mosum agreed TP
# example to run:
#load("/Users/lumeng/Dropbox/mengluchu/bfast/tssarar2.Rdata")
#load("/Users/lumeng/Dropbox/mengluchu/bfast/pdd.Rdata")
#sptssarar2<-pvaluepoint2sp(tssarar2,x=c(58930:59079),y=c(48210:48359),xoff=1,yoff=1,pvalue=0.05,crs=CRS("+proj=utm +zone=21 +south")) #map changes from an array of p-values
#tfpn1<-tfpn(sptssarar2,pdd)
#TPsp<-tfpn1[[1]]
#FNsp<-tfpn1[[2]]
#FPsp<-tfpn1[[3]]
#TNsp<-tfpn1[[4]]

#example:
#get the time 
#timeroccu<-time(mochangepointcu@time[mochangepointcu@index[,2]] )
#timestablecu<-time(mochangepoint2cu@time[mochangepoint2cu@index[,2]] )
#plotting the time histogram 
#jpeg("time cusum roc stable.jpeg ",width = 480,height = 480)
#par(mfrow=c(1,2))
#par(mar=c(5,2,2,2))
#hist.time(timeroccu,title=" time of cusum roc detected change")
#hist.time(timestablecu,title=" time of cusum stable detected change")
#dev.off()
#
#time differences
#bfa.de<-timedif(sts1=changests,sts2=deterpoinf)
#edi.de<-timedif(sts1=stfdfedivi,sts2=deterpoinf)
#bfa.edi<-timedif(sts1=changests,sts2=stfdfedivi)
#bfa.de1<-bfa.de/365
#edi.de1<-edi.de/365
#bfa.edi1<-bfa.edi/3600/24/365
#hist(as.numeric(bfa.de1 ),,breaks=10,main="Differences in time: Bfast Vs.Deter",xlab="year")
#v.hist$counts <- v.hist$counts/sum(v.hist$counts)
#plot(v.hist,ylab="probability",main="Differences in time: Bfast Vs.DETER",xlab="year")
