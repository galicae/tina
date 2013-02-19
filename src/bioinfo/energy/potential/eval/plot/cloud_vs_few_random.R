#file containing data:
#ids,nativescore,shuffle0score,shuffle1score,etc...


file <- "Desktop/tmp.tmp";
a <- read.delim(file,header=TRUE);
aa <- na.omit(a[,c(2:length(a[1,]))]);
bucketnumber <- 50
main=paste0("histogram comparing native vs ",length(a[1,])-2," random sequences on ",length(a[,1])," structures")
xlab="score"
ylab="frequency"
names <- names(a)[c(2:length(a[1,]))]

borders <- c(floor(range(aa)[1]),ceiling(range(aa)[2]));
range <- r[2]-r[1];
bucketsize <- ceiling(range/bucketnumber)
bucketborders <- seq(borders[1],borders[2]+bucketsize,bucketsize);
             
for(i in 1:length(aa[1,])){
  h <- hist(aa[,i],bucketborders,plot=FALSE);
  if(i == 1){
    d <- h$counts;
  }else{
    d <- cbind(d, h$counts);
  }
}
col <- c(2:(length(d[1,])+1));
barplot(t(d),beside=TRUE,col=col,axisnames=TRUE,names.arg=bucketborders[c(2:length(bucketborders))],cex.names=0.5,las=2,main=main,xlab=xlab,ylab=ylab);
legend((bucketnumber-1)*4,max(d),lwd=6,xjust=TRUE,bty="n",legend=names,lty=1, col=col)