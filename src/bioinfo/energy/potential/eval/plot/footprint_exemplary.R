file <- "Desktop/tmp.tmp";
a <- read.delim(file,header=TRUE);
main <- paste0("energy footprint of native vs ",length(a[1,])-1," shuffles")

cummulativeMean <- function(a){
  for(i in 1:length(a)){
    if(i == 1){
      c <- a[i];
    }else{
      c <- rbind(c,sum(a[c(1:i)])/i);
    }
  }
  return(c);
}

names <- c(names(a),"cummulative");
colors = c(c(1:length(a)),"lightgrey")
linetype <- c(1,1,1,2);
plot(c(1,length(a[,1])),c(range(a)[1],range(a)[2]),col="white",ylab="energy values (cummulative)",xlab="peptide length",main=main);
abline(h=0,col="lightgrey")
for(i in 1:length(a)){
    tmp <- cummulativeMean(a[,i]);
    lines(tmp,col=i,lty=2);
    lines(a[,i],col=i,lty=1);
}
legend(xjust=1,bty="n",x=length(a[,1]),y=range(a)[2],legend=names,col=colors,lty=linetype,cex=0.5);