require(R.matlab)
b<-readMat("benefits2.mat")


ben <- b$benefit2.0[,]
bad <- sapply(ben, function(x) is.nan(x))
hdr.den(ben[!bad], den=density(ben[!bad]))
rug(ben[!bad])

ben <- b$benefit2.01
bad <- sapply(ben, function(x) is.nan(x))
#lines(density(ben[!bad,], bw=bw), col="red")
hdr.den(ben[!bad,], den=density(ben[!bad,]))
rug(ben[!bad,])

ben <- b$benefit.05
bad <- sapply(ben, function(x) is.nan(x))
#lines(density(ben[!bad,], bw=bw), col="red")
hdr.den(ben[!bad,], den=density(ben[!bad,]))
rug(ben[!bad,])

ben <- b$benefit2.1
bad <- sapply(ben, function(x) is.nan(x))
#lines(density(ben[!bad,], bw=bw), col="red")
hdr.den(ben[!bad,], den=density(ben[!bad,]))
rug(ben[!bad,])

ben <- b$benefit2.15
bad <- sapply(ben, function(x) is.nan(x) || x < 0 || x > 2e4)
#lines(density(ben[!bad,], bw=bw), col="red")
hdr.den(ben[!bad,], den=density(ben[!bad,]))
rug(ben[!bad,])


ben <- b$benefit.5
bad <- sapply(ben, function(x) is.nan(x))
#lines(density(ben[!bad,], bw=bw), col="red")
hdr.den(ben[!bad,], den=density(ben[!bad,]))
rug(ben[!bad])


ben <- b$benefit2.10
bad <- sapply(ben, function(x) is.nan(x))
plot(density(ben[!bad,]), col="blue")
hdr.den(ben[!bad,], den=density(ben[!bad,]))
rug(ben[!bad,])

hist(ben[!bad,])
ben <- b$benefit2.0
bad <- sapply(ben, function(x) is.nan(x))
#lines(density(ben[!bad,], bw=bw), col="green")
