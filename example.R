library(languageR)
library(robustbase)


source('tautests-public.txt')

data(selfPacedReadingHeid)
x <- selfPacedReadingHeid
x2 <- subset(x, subset = (Condition == "baseheid"))

n <- nrow(x2)

# set the design matrix

# Are here any outliers in the data?

# Build the design matrix
xa <- model.matrix(RT ~ RT4WordsBack + RT3WordsBack + 
RT1WordBack + RTtoPrime +
RT2WordsLater + LengthInLetters + RT2WordsBack + 
FamilySize + RT1WordLater +  Rating + NumberOfSynsets +
BaseFrequency + RootFrequency, data=x2)

# Compute the tau-regression estimators
a <- fasttau_Opt(x=xa, y = x2$RT)

y <- x2$RT
# Compute robust leverage measures
# (Mahalanobis distances computed in the design space)
library(rrcov)
b2 <- .1278
res <- as.vector(y - xa %*% a$beta) / (a$tauscale/sqrt(b2))
x.cov <- CovMMest(x=xa[,-1])
mu <- getCenter(x.cov)
x.cov <- getCov(x.cov)
x.ma <- mahalanobis(x=xa[,-1], center=mu, cov=x.cov)

# Plot standardized residuals versus Mahalanobis distances
plot(res ~ x.ma, pch=19, col='gray', cex=1.3, ylab='Standardized residuals',
xlab='Robust Distances')
# +/- 3.5 standard deviations
abline(h=3.5, lwd=2, col='gray40', lty=2)
abline(h=-3.5, lwd=2, col='gray40', lty=2)
abline(h=0, lwd=2)
# Cut-off points for high-leverage points
abline(v=qchisq(.995, df=13), lwd=2, lty=2, col='gray40')
# Highlight outliers
oo <- (abs(res)>3.5)
points(res ~ x.ma, pch=19, col='gray30', cex=1.3, subset=oo)

# Now show the estimated coefficients and
# p-values for each individual t-test
# using an asymptotic normal approximation for
# the distribution of the regression estimators
uu <- abs(a$beta/sqrt(diag(a$cov)))
round(2*(1 - pt(uu, df=nrow(xa)-ncol(xa))), 3)

#(Intercept)     0.000
#RT4WordsBack    0.000
#RT3WordsBack    0.000
#RT1WordBack     0.000
#RTtoPrime       0.000
#RT2WordsLater   0.000
#LengthInLetters 0.014
#RT2WordsBack    0.045
#FamilySize      0.223
#RT1WordLater    0.394
#Rating          0.842
#NumberOfSynsets 0.188
#BaseFrequency   0.619
#RootFrequency   0.860

# now using the FRB 
b2 <- 0.1278
fi <- as.vector(xa %*% a$beta)
rr <- y - fi
tmp2 <- tautestRBPairs(x=xa, fi=fi, rr=rr, beta=a$beta, sigma=a$scale, R=1000)$betastR
a$cov <- var(tmp2)
uu <- abs(a$beta/sqrt(diag(a$cov)))
round(2*(1 - pt(uu, df=nrow(xa)-ncol(xa))), 3)

#(Intercept)     0.000
#RT4WordsBack    0.002
#RT3WordsBack    0.000
#RT1WordBack     0.000
#RTtoPrime       0.012
#RT2WordsLater   0.004
#LengthInLetters 0.040
#RT2WordsBack    0.240
#FamilySize      0.171
#RT1WordLater    0.470
#Rating          0.896
#NumberOfSynsets 0.354
#BaseFrequency   0.714
#RootFrequency   0.849


# remove all with (rounded) p-values > 0.05

# Build design matrices to compute the FRB p-values
# for the robust tests of the hypothesis
# that all the above betas are indeed zero

xa <- model.matrix(RT ~ RT4WordsBack + RT3WordsBack + 
RT1WordBack + RTtoPrime +
RT2WordsLater + LengthInLetters + RT2WordsBack + 
FamilySize + RT1WordLater +  Rating + NumberOfSynsets +
BaseFrequency + RootFrequency, data=x2)

x0 <- model.matrix(RT ~ RT4WordsBack + RT3WordsBack + 
RT1WordBack + RTtoPrime +
RT2WordsLater + LengthInLetters, data=x2)

y <- x2$RT

# Compute the FRB p-values

set.seed(123)
tau.t <- tautest(x0=x0, xa=xa, y = y, R=1000, Nsamp=2500)
round(unlist(tau.t), 4)

# XXXtest -> test statistic value
# XXXpvalue -> p-value estimated with the asymptotic approximation
# XXXBpvalue -> p-value estimated with the FRB

#        tauHa         tauH0     ratiotest    Sratiotest      ASconstS 
#       0.3168        0.3193        0.0156        0.0315        0.4530 
# Sratiopvalue      SscaleH0      SscaleHa         Wtest         Stest 
#       0.2288        0.2535        0.2496        8.0959       42.3255 
#      Wpvalue       Spvalue   ratiopvalue       ASconst       LRTtest 
#       0.3242        0.0000        0.2211        0.9264        9.5116 
#    LRTpvalue  ratioBpvalue      WBpvalue      SBpvalue    LRTBpvalue 
#       0.2180        0.3154        0.4072        0.2555        0.2645 
#SratioBpvalue 
#       0.1018 


# the bootstrap p-values are systematically larger
# indicating not enough evidence to reject the smaller
# model. In fact, the smaller model gives better
# robust predictions!


# prediction powers via 5-fold CV for different trimming values



# trimmed MSE function
# returns the average of the (1-alpha)100%
# smallest elements in "x", each of them squared
tm <- function(x, alpha) {
n <- length(x)
n0 <- floor(alpha * n)
n <- n - n0
return( mean( (sort(x^2))[1:n] ) ) 
}


n <- dim(x2)[1]
# Number of 5-fold CV runs 
N <- 500 
# trimming proportions
alp.tr <- c(0, 0.01, 0.02, 0.05, 0.10, .15, .2)
la <- length(alp.tr)
# store the TMSPE
mse.r1 <- mse.r2 <- matrix(0, N, la)
set.seed(123)
# 5-fold CV
ii <- c(1, rep(1:5, each=131)) #ii <- c(1:4, rep(1:5, each=124))
for(i in 1:N) {
	ii <- sample(ii)
	pr.r1 <- pr.r2 <- rep(0, n)
	for(j in 1:5) {
		# fit full model 
		a <- fasttau_Opt(x=xa[ii!=j,], y = y[ii!=j])
		# fit reduced model
		b <- fasttau_Opt(x=x0[ii!=j,], y = y[ii!=j])
		# predictions (full and reduced models)
		pr.r1[ii==j] <- as.vector(xa[ii==j,] %*% a$beta)
		pr.r2[ii==j] <- as.vector(x0[ii==j,] %*% b$beta)
	}
	for(j in 1:la) {
		# TMSPE  for full and reduced models
		mse.r1[i, j] <- tm( (y - pr.r1), alpha=alp.tr[j] ) # 0.05
		mse.r2[i, j] <- tm( (y - pr.r2), alpha=alp.tr[j] )
	}
	print(c(i, mean(mse.r1[1:i,4]), mean(mse.r2[1:i,4])))
}

# show the results
boxplot(mse.r1[,1], mse.r2[,1], 
mse.r1[,2], mse.r2[,2],
mse.r1[,3], mse.r2[,3],
mse.r1[,4], mse.r2[,4],
mse.r1[,5], mse.r2[,5],
mse.r1[,6], mse.r2[,6],
mse.r1[,7], mse.r2[,7],
labels=rep(la, each=2), col=rep(c('gray90', 'lightblue'), 5))

# show the results
par(mfrow=c(2, 3))
for(j in 2:7) {
boxplot(mse.r1[,j], mse.r2[,j], col=c('gray90', 'lightblue'), 
main=paste('Trimming: ', alp.tr[j], sep=''),
names = c('Full', 'Reduced'), main='', ylab='TMSPE', cex.axis=1.3, cex.lab=1.3)
}

# 10% trimmed MSPE
boxplot(mse.r1[,5], mse.r2[,5], col=c('gray90', 'lightblue'), 
#main=paste('Trimming: ', alp.tr[j], sep=''),
names = c('Full', 'Reduced'), main='', ylab='TMSPE', cex.axis=1.3, cex.lab=1.3)

