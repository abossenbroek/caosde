b02<-readMat("benefits_a02_p02.mat")
b<-readMat("benefits_a02_p0.mat")
b1<-readMat("benefits_a02_p1.mat")

benlineplot <- function(benefit, color) {
	bad <- sapply(benefit, function(x) is.nan(x) || x < 0 || x > 2e4)
	lines(density(benefit[!bad,]), col=color);
}



jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", 
		"cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
		
cm <- jet.colors(5)

ben <- b0$benefit2.1[,]
bad <- sapply(ben, function(x) is.nan(x) || x < 0 || x > 2e4)
plot(density(ben[!bad]), ylim=c(0, 0.00035), xlim=c(1000, 20000),main="Density plots of a 10% Level Click Fund profit",xlab="Terminal Level Click Fund Profit\n alpha = 0.2  p = 1", col=cm[1])

benlineplot(b02$benefit2.1, cm[2])
benlineplot(b1$benefit2.1, cm[3])

legend(10000, 0.0003, c("Black Scholes market", "p = 0.2", "p = 1"), lty=c(1,1,1), col=c(cm[1], cm[2], cm[3]))