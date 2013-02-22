
a <- read.delim("Desktop/voroeval/randomsamplingPerPeptide/sampling.zscores",header=FALSE)
b <- read.delim("Desktop/voroeval/randomsamplingPerPeptide/sampling.zscores",header=FALSE)
bucketnum <- 50

data <- cbind(a[,2],b[,2])
names <- c("random1","random2")

r <- range(data)
bucketsize <- ceiling((r[2]-r[1])/bucketnum)
buckets <- seq(floor(r[1]),ceiling(r[2])+bucketsize,bucketsize)


for(d in 1:length(data[1,])){
  scores <- data[,d]
  h <- hist(scores,breaks=buckets,plot=FALSE)
  if(d == 1){
    hdata <- h$counts
  }else{
    hdata <- rbind(hdata,h$counts) 
  }
}

col <- c(2:(length(data[1,])+1))
barplot(hdata,beside=TRUE,col=col,names.arg=buckets[c(2:length(buckets))],las=2)
legend(x="topright",bty="n",legend=names,col=col,lty=1,lwd=4)

