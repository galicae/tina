#file containing data as
#id   score
#with firstline native


file <- "Desktop/tmp.tmp";
a <- read.delim(file,header=TRUE);
native <- a[1,2];
aa <- na.omit(a[c(2:length(a[,1])),2]);
bucketnumber <- 50;
main=paste0("histogram comparing native vs ",length(aa)," random sequences")
xlab="score";
ylab="frequency";
names <- c("native","shuffles");

hist(aa,breaks=bucketnumber,col="lightgrey",main=main,xlab=xlab,ylab=ylab);
abline(v=native,col="red",lwd=5);

h <- hist(aa,plot=FALSE);
lty <- c(1,1);
legend <- c("shuffles","native");
col <- c("lightgrey","red");
lwd <- 6;
legend(max(aa),max(h$breaks),xjust=TRUE,bty="n",legend=legend,lty=lty, col=col,lwd=lwd);
dis <- round(100*(length(aa[aa <= native])/length(aa)),digits=1)
legend(min(aa),max(h$breaks),bty="n",legend=c(paste0(dis,"% < native")),col=0,lty=0)
