require(hdrcde)

benefitplot <- function(benefit, name) {
	bad <- sapply(benefit, function(x) is.nan(x) || x < 0 || x > 2e4)
	hdr.den(benefit[!bad,], den=density(benefit[!bad,]), main=name, xlab="Fund profit")
	rug(benefit[!bad,])
	}