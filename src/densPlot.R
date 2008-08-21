require(RColorBrewer)

b<-readMat("benefits2.mat")

benlineplot <- function(benefit, color) {
	bad <- sapply(benefit, function(x) is.nan(x) || x < 0 || x > 2e4)
	lines(density(benefit[!bad,]), col=color);
}



jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", 
		"cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
		
cm <- jet.colors(8)
b<-readMat("benefits2_3.mat")
ben <- b$benefit2.0[,]
bad <- sapply(ben, function(x) is.nan(x) || x < 0 || x > 2e4)
plot(density(ben[!bad]), ylim=c(0, 0.00035), xlim=c(1000, 20000),main="Density plots of the Level Click Fund profit",xlab="Terminal Level Click Fund Profit\n alpha = 0.2  p = 1")

#plot(density(ben[!bad]), main="Density plots of the Level Click Fund profit",xlab="Terminal Level Click Fund Profit\n alpha = 0.2  p = 0")

benlineplot(b$benefit2.01,cm[1])
benlineplot(b$benefit2.05,cm[2])
benlineplot(b$benefit2.1,cm[3])
benlineplot(b$benefit2.2,cm[4])
benlineplot(b$benefit2.5,cm[5])
benlineplot(b$benefit2.10,cm[6])
benlineplot(b$benefit2.15,cm[7])

legend(10000, 0.0003, c("no click fund", "l = 1%", "l = 5%", "l = 10%", "l = 20%", "l = 50%", "l = 100%", "l = 150%"), lty=c(1, 1, 1, 1, 1, 1, 1, 1), col=c("black", cm[1], cm[2], cm[3], cm[4], cm[5], cm[6], cm[7]))

