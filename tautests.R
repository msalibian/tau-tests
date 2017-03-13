
tautest <- function(x0, xa, y, R, Nsamp){

oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv());
on.exit(assign(".Random.seed", oldSeed, envir=globalenv()));

#R = 2000; % bootstrap samples
#
#b2 = 0.1278;
RWLStol = 1e-5;
Mscaletol = 1e-5;
#Nsamp = 100;

c1 = 0.4046;
c2 = 1.090;
b1 = 0.500;
b2 = 0.1278;

# control for fastS()
# control.N : number of random subsamples (e.g. 500)
#        .k : number of initial IRWLS steps on each candidate (e.g. 2)
#        .t : number of best candidates to fully improve (e.g. 5)
#        .r : number of iterations in scale approximation in case approx=1 (e.g. 2)
#        .approx : if 0, fully compute S-scale, otherwise approximate

n <- nrow(xa)
p <- ncol(xa)
q <- p - ncol(x0)
p0 <-  ncol(x0)

fasttauresHa <- fasttau_Opt(xa, y, N=Nsamp, RWLStol=RWLStol, Mscaletol=Mscaletol, b=b1, c1=c1, c2=c2)
betaHa <- fasttauresHa$beta
tauHa <- fasttauresHa$tauscale/sqrt(b2)
sigmaHa <- fasttauresHa$scale
fasttauresH0 <- fasttau_Opt(x0, y, Nsamp, RWLStol=RWLStol, Mscaletol=Mscaletol, b=b1, c1=c1, c2=c2)
betaH0 <- fasttauresH0$beta
tauH0 <- fasttauresH0$tauscale/sqrt(b2)
sigmaH0 <- fasttauresH0$scale

result.tauH0 <- tauH0
result.tauHa <- tauHa

#Ratio of scales test statistic
Ratiotest <- (tauH0^2 - tauHa^2) / tauHa^2

#Ratio of S-scales test statistic
S0 <- fasttau_Opt(x0, y, Nsamp, RWLStol=RWLStol, Mscaletol=Mscaletol, b=b1, c1=c1, c2=c1)
Sa <- fasttau_Opt(xa, y, Nsamp, RWLStol=RWLStol, Mscaletol=Mscaletol, b=b1, c1=c1, c2=c1)
SRatiotest <- (S0$scale^2 - Sa$scale^2 )/ (Sa$scale^2)
#print(Sa)
SresA <- y - xa %*% Sa$beta
#Ratio of S-scales test
B1 <- mean(psiOptfun(SresA/Sa$scale, c1) * (SresA/Sa$scale))
D0S <- mean(psiderOptfun(SresA/Sa$scale, c1))
H0S <- mean(psiOptfun(SresA/Sa$scale, c1)^2)

result.ASconstS <- B1 %*% D0S / H0S

result.Sratiopvalue <- 1 - pchisq(SRatiotest*n*result.ASconstS, q)

result.SscaleH0 <- S0$scale
result.SscaleHa <- Sa$scale

#Score and Wald and alternative LRT tests
WSteststats <- tautestWaldScores(y, x0, xa, betaH0, betaHa, sigmaHa)
result.Wtest <- WSteststats$W
result.Wpvalue <- 1-pchisq(WSteststats$W, q)
result.Stest <- WSteststats$R
result.Spvalue <- 1-pchisq(WSteststats$R, q)

#Ratio of scales test
result.ratiopvalue <- 1-pchisq(Ratiotest*n*WSteststats$ratio, q)
result.ASconst <- WSteststats$ratio

#LRT tests
result.LRTtest <- WSteststats$LRT
result.LRTpvalue <- 1-pchisq(WSteststats$LRT, q)

#Bootstrap tests

#Transform data to H0 scheme
rr <- y - xa %*% betaHa     #residuals
fi <- x0 %*% betaH0
yy <- fi + rr       #y's under H0

#Tau estimates for data under H0
betaHa <- c(betaH0, rep(0,q))
tauH0 <- tauHa
sigmaH0 <- sigmaHa

RBresultH0 <- tautestRBPairs(x0, fi, rr, betaH0, sigmaH0, R)
taus.H0 <- RBresultH0$taustE
RBresultHa <- tautestRBPairs(xa, fi, rr, betaHa, sigmaHa, R)
taus.Ha <- RBresultHa$taustE

#Bootstrap version of ratio of scales test
RatioBtests <- (taus.H0^2 - taus.Ha^2) / taus.Ha^2
result.ratioBpvalue <- (sum(RatioBtests > Ratiotest)+1) / (R+2)

# Transform data to H0 scheme
rr <- y - xa %*% Sa$beta   # residuals
fi <- x0 %*% S0$beta
yy <- fi + rr   #y's under H0

# S-estimates for data under H0
SbetaH0 <- S0$beta
SbetaHa <- append(SbetaH0, rep(0,q))
SscaleHa <- Sa$scale
SscaleH0 <- SscaleHa

SRBresultH0 <- fastSRBPairs(x0, fi, rr, SbetaH0, SscaleH0, RBresultH0$bootmatrix)
SBscales.H0 <- SRBresultH0$sigmastE
SRBresultHa <- fastSRBPairs(xa, fi, rr, SbetaHa, SscaleHa, RBresultH0$bootmatrix)
SBscales.Ha <- SRBresultHa$sigmastE

# Bootstrap version of ratio of S-scales test
SRatioBtests <- (SBscales.H0^2 - SBscales.Ha^2) / (SBscales.Ha^2)
result.SratioBpvalue <- (sum(SRatioBtests > SRatiotest)+1) / (R+2)

#Transform data to H0 scheme
rr <- y - xa %*% fasttauresHa$beta     #residuals
fi <- x0 %*% fasttauresH0$beta
yy <- fi + rr       #y's under H0

#Bootstrap version of the other tests
WaldBtests <- rep(0, R)
# Wtmp <- solve(WSteststats$V22, t(RBresultHa$betastR[,(p-q+1):p]))  #CHECK UNSURE

ScoreBtests <- rep(0, R)
LRTBtests <- rep(0, R)

for(r in 1:R){

	#Calculate residuals of both models for rth bootstrap sample
	binds <- t(RBresultH0$bootmatrix[,r])
	yyb <- yy[binds]
	Bres <- yyb - x0[binds,,drop=F] %*% RBresultH0$betastR[r,]
	BresA <- yyb - xa[binds,,drop=F] %*% RBresultHa$betastR[r,]


	#Calculate M
	#can't use square root weights for this following one (weights can be negative!)
	weightsM <- ((WSteststats$Wn * psiderOptfun(BresA/RBresultHa$sigmastR[r], c1) + psiderOptfun(BresA/RBresultHa$sigmastR[r], c2)) / RBresultHa$sigmastR[r])
	WM <- weightsM %*% rep(1, p)
	xwM <- xa[binds,] * WM
	M <- t(xa[binds,]) %*% xwM / n
	
	#Calculate Q
	weightsQ <- ((WSteststats$Wn * psiOptfun(BresA/RBresultHa$sigmastR[r], c1) + psiOptfun(BresA/RBresultHa$sigmastR[r], c2))^2)
	WQ <- weightsQ %*% rep(1, p)
	xwQ <- xa[binds,] * WQ
	Q <- t(xa[binds,]) %*% xwQ / n

  # Calculate V
	Minv <- solve(M)
  V <- Minv %*% Q %*% Minv
	result.V22 <- V[(p0+1):p,(p0+1):p]
  
  
  #Calculate Wald test statistic for rth bootstrap sample
	WaldBtests[r] <- RBresultHa$betastR[r,(p-q+1):p] %*% solve(result.V22, as.vector(RBresultHa$betastR[r,(p-q+1):p]))
	#Calculate Score test statistic for rth bootstrap sample
	weightsZ <- (WSteststats$Wn * psiOptfun(Bres/RBresultHa$sigmastR[r], c1) + psiOptfun(Bres/RBresultHa$sigmastR[r], c2))
	Z <- t(xa[binds,(p-q+1):p]) %*% weightsZ / n
#	ScoreBtests[r] <- n * t(Z) %*% solve(WSteststats$C, Z)

	result.C <- Q[(p0+1):p,(p0+1):p]-M[(p0+1):p,1:p0]%*% solve(M[1:p0,1:p0,drop=F], Q[1:p0,(p0+1):p,drop=F])-
	  Q[(p0+1):p,1:p0]%*% solve(M[1:p0,1:p0,drop=F], M[1:p0,(p0+1):p,drop=F])+
	  M[(p0+1):p,1:p0]%*% solve(M[1:p0,1:p0,drop=F], Q[1:p0,1:p0,drop=F])%*%solve(M[1:p0,1:p0,drop=F], M[1:p0,(p0+1):p,drop=F])  

  ScoreBtests[r] <- n * t(Z) %*% solve(result.C, Z)
  
  
  #Calculate LRT test statistic for rth bootstrap sample
	obj0 <- WSteststats$Wn * rhoOptfun(Bres/RBresultHa$sigmastR[r], c1) + rhoOptfun(Bres/RBresultHa$sigmastR[r], c2)
	objA <- WSteststats$Wn * rhoOptfun(BresA/RBresultHa$sigmastR[r], c1) + rhoOptfun(BresA/RBresultHa$sigmastR[r], c2)
	LRTBtests[r] <- 2 * sum(obj0 - objA) * WSteststats$D0 / WSteststats$H0

}

#Bootstrap version of the Wald test
result.WBpvalue <- (sum(WaldBtests > as.numeric(WSteststats$W/n))+1) / (R+2)
#Bootstrap version of the Score test
result.SBpvalue <- (sum(ScoreBtests > as.numeric(WSteststats$R))+1) / (R+2)
#Bootstrap version of the LRT test
result.LRTBpvalue <- (sum(LRTBtests > as.numeric(WSteststats$LRT))+1) / (R+2)

return(list(tauHa=result.tauHa, tauH0=result.tauH0, ratiotest=Ratiotest, Sratiotest=SRatiotest, ASconstS=result.ASconstS, 
Sratiopvalue=result.Sratiopvalue, SscaleH0=result.SscaleH0, SscaleHa=result.SscaleHa, Wtest=result.Wtest, 
Stest=result.Stest, Wpvalue=result.Wpvalue, Spvalue=result.Spvalue, ratiopvalue=result.ratiopvalue, ASconst=result.ASconst, 
LRTtest=result.LRTtest, LRTpvalue=result.LRTpvalue, ratioBpvalue=result.ratioBpvalue, WBpvalue=result.WBpvalue, 
SBpvalue=result.SBpvalue, LRTBpvalue=result.LRTBpvalue, 
SratioBpvalue = result.SratioBpvalue))

}


#-------------- Functions to do Fast Robust Bootstrap ---------------------

tautestRBPairs <- function(x, fi, rr, beta, sigma, R){
  
  x <- as.matrix(x)
  
  n <- length(fi)
  p <- ncol(x)
  
  c1 = 0.4046
  c2 = 1.090
  b1 = 0.500
  b2 = 0.1278
  
  y <- fi + rr
  
  res <- as.vector(y - x %*% beta)
  #tau <- sqrt(sigma^2 * mean( rhoOptfun(res/sigma,c2) ) / b2)
  
  Wn_1 <- sum(2*rhoOptfun(res/sigma, c2) - psiOptfun(res/sigma, c2) * (res/sigma))
  Wn_2 <- sum(psiOptfun(res/sigma, c1) * (res/sigma))
  Wn <- Wn_1/Wn_2
  
  weights <- (Wn * fw(res/sigma, c1) + fw(res/sigma, c2))
  scaleweights <- rhoOptfun(res/sigma, c1)
  
  sqweights <- sqrt(weights/sigma)
  #sqW <- sqweights %*% rep(1,p)
  xw <- x * sqweights #sqW
  yw <- y * sqweights
  A <- t(xw) %*% xw
  Ainv <- solve(A)
  
  #can't use square root weights for this following one (weights can be negative!)
  weightsB <- ((Wn * psiderOptfun(res/sigma, c1) + psiderOptfun(res/sigma, c2)) / sigma)
  WB <- weightsB # %*% rep(1,p)
  xwB <- x * WB
  B <- t(x) %*% xwB
  
  weightsv <- ((Wn * psiderOptfun(res/sigma, c1) + psiderOptfun(res/sigma, c2)) * res / sigma^2)
  Wv <- weightsv # %*% rep(1,p)
  #v <- colSums(x * Wv)
  v <- as.vector( t(x) %*% Wv )
  
  weightsb <- psiOptfun(res/sigma, c1)
  Wb <- weightsb # %*% rep(1,p)
  #b <- colSums(x * Wb)
  b <- as.vector( t(x) %*% Wb )
  
  corr11 <- Ainv %*% B
  corr12 <- Ainv %*% v
  #corr21 <- matrix(1 / b1 / n * b, 1, p)
  corr21 <- b / b1 / n 
  corr22 <- mean(psiOptfun(res/sigma, c1) * (res/sigma)) / b1
  
  #corrmat <- solve(rbind(cbind(corr11, corr12), cbind(corr21, corr22)))
  corrmat1 <- cbind(corr11, corr12)
  corrmat2 <- c(corr21, corr22)
  corrmat <- solve( rbind(corrmat1, corrmat2) )
  
  #set.seed(123)
  #To fix the seed so that we can reproduce the results
  set.seed(123)
  bootmatrix <- matrix(sample(1:n, R*n, repl=TRUE), n, R)
  
  betast <- matrix(0,R,p)
  sigmast <- rep(0,R)
  betastR <- matrix(0,R,p)
  sigmastR <- rep(0,R)
  sigmastE <- rep(0,R)
  taustE <- rep(0,R)
  
  #xaux <- solve(t(xw) %*% xw)
  #xaux <- Ainv
  
  for(r in 1:R){
    
    binds <- as.vector(bootmatrix[,r])
    
    yy <- y[binds]
    xxw <- xw[binds,]
    yyw <- yw[binds]
    
    betast[r,] <- solve(t(xxw) %*% xxw, t(xxw) %*% yyw)
    sigmast[r] <- sigma * mean(scaleweights[binds]) / b1
    
    RBres <- c(beta, sigma) + as.vector( corrmat %*% c(betast[r,] - beta, sigmast[r] - sigma) )
    
    betastR[r,] <- RBres[1:p]
    sigmastR[r] <- RBres[p+1]

    #Exact calculation of scale corresponding to approximated beta

    sigmastE[r] <- Mscale(yy - x[binds,] %*% betastR[r,], b1, c1, RBres[p+1], 1e-5)
    
    taustE[r] <- sqrt(sigmastE[r]^2 * mean(rhoOptfun((yy - x[binds,] %*% betastR[r,]) / 
      sigmastE[r], c2)) / b2)
    
  }
  
  return(list(betastR = betastR, sigmastR = sigmastR, sigmastE = sigmastE, taustE = taustE, 
              bootmatrix = bootmatrix))
  
}

# ----------------------------------------------------------------------------------------------
fastSRBPairs <- function(x, fi, rr, beta, sigma, bootmatrix){
  
  n <- length(fi)
  p <- ncol(x)
  
  R <- ncol(bootmatrix)
  
  #rand('seed',31);
  # To fix the seed so that we can reproduce the results
  #s <- RandStream('swb2712','Seed',31);
  #RandStream.setDefaultStream(s);
  
  c1 <- 0.4046
  c2 <- 1.090
  b1 <- 0.500
  b2 <- 0.1278
  
  betast <- matrix(rep(0,R*p),R,p)
  sigmast <- matrix(rep(0,R),R,1)
  betastR <- matrix(rep(0,R),R,p)
  sigmastE <- matrix(rep(0,R),R,1)
  
  y <- fi + rr
  # res <- y - x %*% beta
  res <- rr
  
  weights <- fw(res/sigma, c1)
  scaleweights <- rhoOptfun(res/sigma, c1)
  
  sqweights <- (weights/sigma)^(1/2)
  sqW <- sqweights %*% rep(1,p)
  xw <- x * sqW
  yw <- y * sqweights
  A <- t(xw) %*% xw
  
  weightsM <- psiderOptfun(res/sigma, c1)
  xwM <- x * (weightsM %*% rep(1,p))
  M <- sigma * solve((t(x) %*% xwM), A)
  
  an <- mean( psiOptfun(res/sigma, c1)* res / sigma ) / b1;
  
  d <- solve((t(x) %*% xwM), (t(x) %*% (res * weightsM / sigma))) / an
  
  for(r in 1:R){
    
    binds <- as.vector(bootmatrix[,r])
    
    yy <- y[binds]
    xxw <- xw[binds,]
    yyw <- yw[binds]
    
    betast[r,] <- solve(t(xxw) %*% xxw, t(xxw) %*% yyw)
    
    sigmast[r] <- sigma / b1 * mean(scaleweights[binds])
    
    betastR[r,] <- t(beta + M %*% (betast[r,] - beta) + d * (sigmast[r] - sigma))
    
    # Exact calculation of scale corresponding to approximated beta
    rb <- yy - x[binds,,drop=F] %*% betastR[r,]
    sigmastE[r] <- Mscale(rb,b1,c1,sigmast[r],1e-5)
  }
  
  return(list(bootmatrix=bootmatrix, betastR=betastR, sigmastE=sigmastE))
}


#----------------- Functions to calculate test statistics -----------------------------------

tautestWaldScores <- function(y, x0, xa, beta0, betaA, sigmaA){

xa <- as.matrix(xa)
x0 <- as.matrix(x0)

n <- nrow(xa)
pa <- ncol(xa)
p0 <- ncol(x0)

c1 = 0.4046
c2 = 1.090
b1 = 0.500
b2 = 0.1278

#Calculate M, Q, V based on the largest model
resA <- y - xa %*% betaA
res <- y - x0 %*% beta0

Wn_1 <- colSums(2*rhoOptfun(resA/sigmaA, c2) - psiOptfun(resA/sigmaA, c2) * (resA/sigmaA))
Wn_2 <- colSums(psiOptfun(resA/sigmaA, c1) * (resA/sigmaA))
Wn <- as.numeric(Wn_1/Wn_2)

#Calculate M
#can't use square root weights for this following one (weights can be negative!)
weightsM <- ((Wn * psiderOptfun(resA/sigmaA, c1) + psiderOptfun(resA/sigmaA, c2)) / sigmaA)
WM <- weightsM %*% rep(1, pa)
xwM <- xa * WM
M <- t(xa) %*% xwM / n

#Calculate Q
weightsQ <- ((Wn * psiOptfun(resA/sigmaA, c1) + psiOptfun(resA/sigmaA, c2))^2)
WQ <- weightsQ %*% rep(1, pa)
xwQ <- xa * WQ
Q <- t(xa) %*% xwQ / n

#Calculate V
Minv <- solve(M)
V <- Minv %*% Q %*% Minv

#Calculate Wald test statistic

result.V22 <- V[(p0+1):pa,(p0+1):pa]
result.W <- n * betaA[(p0+1):pa] %*% solve(result.V22, betaA[(p0+1):pa])

#Calculate Score test statistic
weightsZ <- (Wn * psiOptfun(res/sigmaA, c1) + psiOptfun(res/sigmaA, c2))
Z <- t(xa[,(p0+1):pa]) %*% weightsZ / n
#print(M[1:p0,1:p0])
#print(M[1:p0,(p0+1):pa])
#M221 <- M[(p0+1):pa,(p0+1):pa] - M[(p0+1):pa,1:p0] %*% solve(M[1:p0,1:p0,drop=F], M[1:p0,(p0+1):pa,drop=F])
#result.C <- M221 %*% V[(p0+1):pa,(p0+1):pa] %*% t(M221)


result.C <- Q[(p0+1):pa,(p0+1):pa]-M[(p0+1):pa,1:p0]%*% solve(M[1:p0,1:p0,drop=F], Q[1:p0,(p0+1):pa,drop=F]) - 
Q[(p0+1):pa,1:p0]%*% solve(M[1:p0,1:p0,drop=F], M[1:p0,(p0+1):pa,drop=F])+
  M[(p0+1):pa,1:p0]%*% solve(M[1:p0,1:p0,drop=F], Q[1:p0,1:p0,drop=F])%*%
  solve(M[1:p0,1:p0,drop=F], M[1:p0,(p0+1):pa,drop=F])




result.R <- n * t(Z) %*% solve(result.C, Z)

#Calculate ratio of scales test statistic

M2 <- mean(rhoOptfun(resA/sigmaA, c2))

psi0 <- Wn * psiOptfun(resA/sigmaA, c1) + psiOptfun(resA/sigmaA, c2)
H0 <- mean(psi0^2)
#result.H0 <- H0

psi0der <- Wn * psiderOptfun(resA/sigmaA, c1) + psiderOptfun(resA/sigmaA, c2)
D0 <- mean(psi0der)
#result.D0 <- D0

result.ratio <- 2 * D0 %*% M2 / H0

#Calculate LR Type test statistic
obj0 <- Wn * rhoOptfun(res/sigmaA, c1) + rhoOptfun(res/sigmaA, c2)
objA <- Wn * rhoOptfun(resA/sigmaA, c1) + rhoOptfun(resA/sigmaA, c2)

result.LRT <- 2 * colSums(obj0 - objA) %*% D0 / H0

return(list(Wn = Wn, V22 = result.V22, W = result.W, C = result.C, R = result.R, H0 = H0, 
D0 = D0, ratio = result.ratio, LRT = result.LRT))

}

#--------------------------------------------------------------------------------------------

rhoOptfun <- function(x, c){

#Computes optimal rho function

tmp <- x^2 / 2 / (3.25*c^2)
tmp2 <- (1.792 - 0.972 * x^2 / c^2 + 0.432 * x^4 / c^4 - 0.052 * x^6 / c^6 + 0.002 * x^8 / c^8) / 3.25
tmp[abs(x) > 2*c] <- tmp2[abs(x) > 2*c]
tmp[abs(x) > 3*c] <- 1
return(tmp)

}

#--------------------------------------------------------------------------------------------

psiOptfun <- function(x, c){

#Computes optimal rho function's first derivative

tmp <- x / (3.25*c^2)
tmp2 <- (-1.944 * x / c^2 + 1.728 * x^3 / c^4 - 0.312 * x^5 / c^6 + 0.016 * x^7 / c^8) / 3.25
tmp[abs(x) > 2*c] <- tmp2[abs(x) > 2*c]
tmp[abs(x) > 3*c] <- 0
return(tmp)

}

#--------------------------------------------------------------------------------------------

fw <- function(x, c){

# weight function = psi(x)/x

tmp <- x^0 / (3.25*c^2)
tmp2 <- (-1.944  / c^2 + 1.728 * x^2 / c^4 - 0.312 * x^4 / c^6 + 0.016 * x^6 / c^8) / 3.25
tmp[abs(x) > 2*c] <- tmp2[abs(x) > 2*c]
tmp[abs(x) > 3*c] <- 0
tmp[abs(tmp) < 10e-10] <- 0
return(tmp)

}

#--------------------------------------------------------------------------------------------

psiderOptfun <- function(x,c){

tmp <- x^0 / (3.25*c^2)
tmp2 <- (-1.944 / c^2 + 5.184 * x^2 / c^4 - 1.560 * x^4 / c^6 + 0.112 * x^6 / c^8) / 3.25
tmp[abs(x) > 2*c] <- tmp2[abs(x) > 2*c]
tmp[abs(x) > 3*c] <- 0
return(tmp)

}

fasttau_Opt <- function(x, y, N=500, kk=2, tt=5, rr=2, approximate=0, RWLStol=1e-5, Mscaletol=1e-5, b1=0.5, c1=0.4046, c2=1.09)
{
# fast-tau algorithm for linear regression
#
# tau-estimate is tuned to have 95% efficiency, and 50% bdp,
# using Optimal rho-function
#
# INPUT:
# 	y : response vector (n x 1)
# 	x : covariates matrix (n x p), possibly including intercept column
# 	N : number of elemental starts, e.g. 500 
#	  kk : number of initial IRWLS steps on each candidate
#	  tt : number of best solutions to RWLS-iterate until convergence
#   rr : number of iterations in scale approximation in case approximate=1
#   approximate : if 0, fully compute S-scale in objective function evaluation, 
#               otherwise approximate 
# OUTPUT:
# 	res$beta : tau-estimate of regression coefficients
#	res$scale : tau-estimate of residual scale

if (tt<1) stop("parameter tt should be at least 1")

x <- as.matrix(x)

n <- nrow(x)
p <- ncol(x)

#c1 <- .4046
#b1 <- .5
#c2 <- 1.09
#b2 <- .1278

bestbetas <- matrix(0, p, tt)
bestscales <- 1e20 * rep(1, tt)
besttauscales <- 1e20 * rep(1, tt)
worsti <- 1
rworst <- y

for (i in 1:N) {
    # find a p-subset in general position.
    singular <- 1; itertest <- 1
    while (singular==1 && itertest<100) {
        ranset <- sample(n, p) #randomset(n,p)
        xj <- x[ranset,]
        yj <- y[ranset] 
        bj <- as.matrix(qr.coef(qr(xj),yj))
        singular <- any(!is.finite(bj))
	      itertest <- itertest + 1
    }
    if (itertest==100) stop("Too many degenerate subsamples")

    # perform kk steps of IRLS on elemental start
    if (kk > 0) {
        tmp <- IWLSiteration(x, y, bj, 0, kk, RWLStol, b1, c1, c2)
        betarw <- tmp$betarw
        resrw <- y - x %*% betarw
        scalerw <- tmp$scalerw
    }
    else {
        betarw <- bj
        resrw <- y - x %*% betarw
        scalerw <- median(abs(resrw))/.6745
    }
    
  # long-term memory vector, for finding a special extra candidate at the end :
    if (i > 1) LTMvec = LTMvec + abs(resrw)
    else LTMvec = abs(resrw)
    
    
    # check whether new subsample yields one of best t tau-objective values
    
    if (!approximate) { # compute actual scale, but use tau-conditions!
      scaletest1 <- mean(rhoOpt(resrw / bestscales[worsti],c1)) < b1
      scaletest2 <- sum(rhoOpt(resrw / bestscales[worsti],c2)) < sum(rhoOpt(rworst/bestscales[worsti],c2))
      if (scaletest1 || scaletest2) {
			   # if conditions fulfulled, compute objective value
            snew <- Mscale(resrw, b1, c1, scalerw, Mscaletol)
            taunew <- snew * sqrt(mean(rhoOpt(resrw/snew,c2)))
            if (taunew < besttauscales[worsti]) {
				    # if candidate has indeed better tau than the worst of the tt best until now, keep it. 
              besttauscales[worsti] <- taunew
              bestscales[worsti] <- snew
              bestbetas[,worsti] <- betarw
              worsti <- which.max(besttauscales) 
              rworst <- y - x %*% bestbetas[,worsti]
            }
      }
    }
    else { # or just compute approximations (and don't bother with the conditions)
       snew = scalerw;
       if (rr>0) {
          for (kstep in 1:rr) { 
              snew <- sqrt( snew^2 * mean( rhoOpt(resrw/snew,c1) ) / b1 )
          }
       }
       taunew <- snew * sqrt(mean(rhoOpt(resrw/snew,c2)))
       if (taunew < besttauscales[worsti]) {
	     # if candidate has indeed better tau than the worst of the tt best until now, keep it. 
          besttauscales[worsti] <- taunew
          bestscales[worsti] <- snew
          bestbetas[,worsti] <- betarw
          worsti <- which.max(besttauscales) 
          rworst <- y - x %*% bestbetas[,worsti]
       }
    }
}

# consider an extra subsample, made up of badly fit observations

IXLTM <- order(LTMvec, decreasing=T)
singular <- 1 
extrasize <- p
while (singular==1) {
    xs <- x[IXLTM[1:extrasize],]
    ys <- y[IXLTM[1:extrasize]]
	  bbeta <- as.matrix(qr.coef(qr(xs),ys))
	  singular <- any(!is.finite(bbeta))
    extrasize <- extrasize + 1
}

# perform kk steps of IRLS on elemental start
if (kk > 0) {
    tmp <- IWLSiteration(x, y, bbeta, 0, kk, RWLStol, b1, c1, c2)
    betarw <- tmp$betarw
    resrw <- y - x %*% betarw
    scalerw <- tmp$scalerw
}
else {
    betarw <- bbeta
    resrw <- y - x %*% betarw
    scalerw <- median(abs(resrw))/.6745
}

# check whether this candidate yields one of best t tau-objective values

if (!approximate) { # compute actual scale, but use tau-conditions!
  scaletest1 <- mean(rhoOpt(resrw / bestscales[worsti],c1)) < b1
  scaletest2 <- sum(rhoOpt(resrw / bestscales[worsti],c2)) < sum(rhoOpt(rworst/bestscales[worsti],c2))
  if (scaletest1 || scaletest2) {
	   # if conditions fulfulled, compute objective value
        snew <- Mscale(resrw, b1, c1, scalerw, Mscaletol)
        taunew <- snew * sqrt(mean(rhoOpt(resrw/snew,c2)))
        if (taunew < besttauscales[worsti]) {
		    # if candidate has indeed better tau than the worst of the tt best until now, keep it. 
          besttauscales[worsti] <- taunew
          bestscales[worsti] <- snew
          bestbetas[,worsti] <- betarw
          worsti <- which.max(besttauscales) 
          rworst <- y - x %*% bestbetas[,worsti]
        }
  }
}
else { # or just compute approximations (and don't bother with the conditions)
   snew = scalerw;
   if (rr>0) {
     for (kstep in 1:rr) { 
        snew <- sqrt( snew^2 * mean( rhoOpt(resrw/snew,c1) ) / b1 )
     }
   }  
   taunew <- snew * sqrt(mean(rhoOpt(resrw/snew,c2)))
   if (taunew < besttauscales[worsti]) {
   # if candidate has indeed better tau than the worst of the tt best until now, keep it. 
      besttauscales[worsti] <- taunew
      bestscales[worsti] <- snew
      bestbetas[,worsti] <- betarw
      worsti <- which.max(besttauscales) 
      rworst <- y - x %*% bestbetas[,worsti]
   }
}

superbesttauscale <- 1e20

# RWLS-iterate each of the best tt candidates until convergence, and retain the best result
for (i in 1:tt) {
    tmp <- IWLSiteration(x, y, bestbetas[,i], bestscales[i], 500, RWLStol, b1, c1, c2)
    resrw <- y - x %*% tmp$betarw
    tauscalerw <- tmp$scalerw * sqrt(mean(rhoOpt(resrw/tmp$scalerw,c2)))
    if (tauscalerw < superbesttauscale) {
        superbesttauscale <- tauscalerw
        superbestbeta <- tmp$betarw
        superbestscale <- tmp$scalerw
    }
}

superbestscale <- Mscale(y - x%*%superbestbeta, b1, c1, superbestscale, Mscaletol)
superbesttauscale <- superbestscale * sqrt(mean(rhoOpt((y - x%*%superbestbeta)/superbestscale,c2)))

# Olive's two extra candidates:

# add LS candidate
betaLS <- as.matrix(qr.coef(qr(x),y))
resLS <- y - x %*% betaLS
scaleLS <- median(abs(resLS))/.6745  
scaletest1 <- mean(rhoOpt(resLS / superbestscale,c1)) < b1
scaletest2 <- sum(rhoOpt(resLS / superbestscale,c2)) < sum(rhoOpt((y - x%*%superbestbeta)/superbestscale,c2))
if (scaletest1 || scaletest2) {
    snew <- Mscale(resLS, b1, c1, scaleLS, Mscaletol)
    taunew <- snew * sqrt(mean(rhoOpt(resLS/snew,c2)))
    if (taunew < superbesttauscale) {
        superbestscale <- snew
        superbestbeta <- betaLS
        superbesttauscale <- taunew
    }   
}

# add HB candidate
IXmed <- order(abs(y - median(y)))
xhalf <- x[IXmed[1:floor(n/2)],]
yhalf <- y[IXmed[1:floor(n/2)]]
bbeta <- as.matrix(qr.coef(qr(xhalf),yhalf))
# + 10 C-steps
tmp <- IWLSiteration(x, y, bbeta, 0, 10, RWLStol, b1, c1, c2)
betaHB <- tmp$betarw
resHB <- y - x %*% betaHB
scaleHB <- tmp$scalerw
scaletest1 <- mean(rhoOpt(resHB / superbestscale,c1)) < b1
scaletest2 <- sum(rhoOpt(resHB / superbestscale,c2)) < sum(rhoOpt((y - x%*%superbestbeta)/superbestscale,c2))
if (scaletest1 || scaletest2) {
    snew <- Mscale(resHB, b1, c1, scaleHB, Mscaletol)
    taunew <- snew * sqrt(mean(rhoOpt(resHB/snew,c2)))
    if (taunew < superbesttauscale) {
        superbestbeta <- betaHB
        superbesttauscale <- taunew
    }   
}

b2 <- .1278

res <- as.vector(y - x %*% superbestbeta)
#sigma <- superbestscale
sigma <- sqrt( superbestscale^2 * mean( rhoOpt(res/superbestscale,c2) ) / b2 )
Wn_1 <- sum(2*rhoOptfun(res/sigma, c2) - psiOptfun(res/sigma, c2) * (res/sigma))
Wn_2 <- sum(psiOptfun(res/sigma, c1) * (res/sigma))
Wn <- Wn_1/Wn_2
psider0 <- ((Wn * psiderOptfun(res/sigma, c1) + psiderOptfun(res/sigma, c2))  / sigma)
psi0 <- Wn * psiOptfun(res/sigma, c1) + psiOptfun(res/sigma, c2)
A <- solve(t( x ) %*% x  ) * mean(psi0^2) / (mean(psider0))^2 

B0 <- solve( t(x) %*% (x*psider0) )
A0 <- ( t(x * psi0) %*% (x*psi0) )
A1 <- B0 %*% A0 %*% B0

return(list( beta = superbestbeta, tauscale = superbesttauscale, scale = superbestscale, cov=A, cov2=A1 ))

}


# -------------------------------------------------------------------

IWLSiteration <- function(x, y, inib, iniscale, maxiter, tol, b, c1, c2)
{  
# approximate IRWLS iteration; pass maxiter=500, say, if convergence is desired
# e.g. tol = 1e-11

n <- nrow(x)
p <- ncol(x)

res <- y - x %*% inib
if (iniscale == 0)
    scale <- median(abs(res))/.6745
else
    scale <- iniscale

oldbeta <- inib

betadiff <- 2*tol
iter <- 0
while ((betadiff > tol) && (iter < maxiter)) {
    scale <- sqrt( scale^2 * mean( rhoOpt(res/scale,c1) ) / b )
    scaledres <- res/scale
    Wn.teller <- sum(WtellerOpt(scaledres,c2))
    Wn.noemer <- sum(psixOpt(scaledres,c1))
    Wn <- Wn.teller / Wn.noemer
    weights <- (Wn * fwOpt(scaledres,c1) + fwOpt(scaledres,c2))
    sqweights <- sqrt(weights)
    xw <- x * as.vector(sqweights)
    yw <- y * sqweights
    newbeta <- qr.coef(qr(xw),yw)
    if (any(!is.finite(newbeta))) {
        newbeta <- inib
        scale <- iniscale
        break
    }
    betadiff <- sqrt(sum((oldbeta - newbeta)^2))
    res <- y - x %*% newbeta
    oldbeta <- newbeta
    iter <- iter + 1
}

return( list( betarw = newbeta, scalerw = scale ) )

}

#--------------------------------------------------------------------------  

Mscale <- function(u, b, c, initialsc, Mscaletol) 
{
# from Kristel's fastSreg
if (initialsc==0)
    initialsc = median(abs(u))/.6745
maxit <- 100
sc <- initialsc
i <- 0 
err <- 1
while  (( i < maxit ) & (err > Mscaletol)) {
    sc2 <- sqrt( sc^2 * mean(rhoOpt(u/sc,c)) / b)
    err <- abs(sc2/sc - 1)
    sc <- sc2
    i <- i+1
}

return(sc)

}

#---------------------------------------------------------------------------------------

randomset <- function(tot,nel) {
ranset <- rep(0,nel)
for (j in 1:nel) {
 	num <- ceiling(runif(1)*tot)
   	if (j > 1) {
     		while (any(ranset==num)) 
     			num <- ceiling(runif(1)*tot)
	}
	ranset[j] <- num
}
return(ranset)
}

# --------------------------------------------------------------------

rhoOpt <- function(x, cc)
{
	tmp <- x^2 / 2 / (3.25*cc^2)
	tmp2 <- (1.792 - 0.972 * x^2 / cc^2 + 0.432 * x^4 / cc^4 - 0.052 * x^6 / cc^6 + 0.002 * x^8 / cc^8) / 3.25
	tmp[abs(x) > 2*cc] <- tmp2[abs(x) > 2*cc]
	tmp[abs(x) > 3*cc] <- 1
	tmp
	
}

# --------------------------------------------------------------------

fwOpt <- function(x, cc)
{
	tmp <- (-1.944 / cc^2 + 1.728 * x^2 / cc^4 - 0.312 * x^4 / cc^6 + 0.016 * x^6 / cc^8) / 3.25
	tmp[abs(x) < 2*cc] <- 1 / (3.25*cc^2)
	tmp[abs(x) > 3*cc] <- 0
	tmp[abs(tmp) < 10e-10] <- 0
	tmp
	
}

# --------------------------------------------------------------------

psiOpt <- function(x, cc)
{
	tmp <- x / (3.25*cc^2)
	tmp2 <- (-1.944 * x / cc^2 + 1.728 * x^3 / cc^4 - 0.312 * x^5 / cc^6 + 0.016 * x^7 / cc^8) / 3.25
	tmp[abs(x) > 2*cc] <- tmp2[abs(x) > 2*cc]
	tmp[abs(x) > 3*cc] <- 0
	tmp
	
}

# --------------------------------------------------------------------

psixOpt <- function(x, cc)
{
	tmp <- x^2 / (3.25*cc^2)
	tmp2 <- (-1.944 * x^2 / cc^2 + 1.728 * x^4 / cc^4 - 0.312 * x^6 / cc^6 + 0.016 * x^8 / cc^8) / 3.25
	tmp[abs(x) > 2*cc] <- tmp2[abs(x) > 2*cc]
	tmp[abs(x) > 3*cc] <- 0
	tmp
	
}
# --------------------------------------------------------------------

WtellerOpt <- function(x, cc)
{
	tmp <- (3.584 - 0.864 * x^4 / cc^4 + 0.208 * x^6 / cc^6 - 0.012 * x^8 / cc^8) / 3.25
	tmp[abs(x) < 2*cc] <- 0
	tmp[abs(x) > 3*cc] <- 2
	tmp
	
}

# --------------------------------------------------------------------

psiOpt <- function(x, cc)
{
	tmp <- x / (3.25*cc^2)
	tmp2 <- (-1.944 * x / cc^2 + 1.728 * x^3 / cc^4 - 0.312 * x^5 / cc^6 + 0.016 * x^7 / cc^8) / 3.25
	tmp[abs(x) > 2*cc] <- tmp2[abs(x) > 2*cc]
	tmp[abs(x) > 3*cc] <- 0
	tmp
	
}


psiderOptfun <- function(x,cc){
tmp <- x^0 / (3.25*cc^2)
tmp2 <- (-1.944 / cc^2 + 5.184 * x^2 / cc^4 - 1.560 * x^4 / cc^6 + 0.112 * x^6 / cc^8) / 3.25
tmp[abs(x) > 2*cc] <- tmp2[abs(x) > 2*cc]
tmp[abs(x) > 3*cc] <- 0
return(tmp)
}



#--------------------- Functions to do Fast Robust Bootstrap --------------------------------

tautestRBPairs <- function(x, fi, rr, beta, sigma, R){

x <- as.matrix(x)

n <- length(fi)
p <- ncol(x)

c1 <- 0.4046
c2 <- 1.090
b1 <- 0.500
b2 <- 0.1278

y <- fi + rr

res <- y - x %*% beta
#tau <- sqrt(sigma^2 * mean( rhoOptfun(res/sigma,c2) ) / b2)

Wn_1 <- colSums(2*rhoOptfun(res/sigma, c2) - psiOptfun(res/sigma, c2) * (res/sigma))
Wn_2 <- colSums(psiOptfun(res/sigma, c1) * (res/sigma))
Wn <- Wn_1/Wn_2

weights <- (Wn * fw(res/sigma, c1) + fw(res/sigma, c2))
scaleweights <- rhoOptfun(res/sigma, c1)

sqweights <- (weights/sigma)^(1/2)
sqW <- sqweights %*% rep(1,p)
xw <- x * sqW
yw <- y * sqweights
A <- t(xw) %*% xw
Ainv <- solve(A)

#can't use square root weights for this following one (weights can be negative!)
weightsB <- ((Wn * psiderOptfun(res/sigma, c1) + psiderOptfun(res/sigma, c2)) / sigma)
WB <- weightsB %*% rep(1,p)
xwB <- x * WB
B <- t(x) %*% xwB

weightsv <- ((Wn * psiderOptfun(res/sigma, c1) + psiderOptfun(res/sigma, c2)) * res / sigma^2)
Wv <- weightsv %*% rep(1,p)
v <- colSums(x * Wv)

weightsb <- psiOptfun(res/sigma, c1)
Wb <- weightsb %*% rep(1,p)
b <- colSums(x * Wb)

corr11 <- Ainv %*% B
corr12 <- Ainv %*% v
corr21 <- matrix(1 / b1 / n * b, 1, p)
corr22 <- 1 / b1 * mean(psiOptfun(res/sigma, c1) * (res/sigma))

corrmat <- solve(rbind(cbind(corr11, corr12), cbind(corr21, corr22)))

#set.seed(123)
#To fix the seed so that we can reproduce the results
set.seed(123)
bootmatrix <- matrix(sample(1:n, R*n, repl=TRUE), n, R)

betast <- matrix(rep(0,R*p),R,p)
sigmast <- matrix(rep(0,R),R,1)
#taust <- matrix(rep(0,R),R,1)
betastR <- matrix(rep(0,R),R,p)
sigmastR <- matrix(rep(0,R),R,1)
sigmastE <- matrix(rep(0,R),R,1)
taustE <- matrix(rep(0,R),R,1)

#xaux <- solve(t(xw) %*% xw)
#xaux <- Ainv

for(r in 1:R){

	binds <- as.vector(bootmatrix[,r])

	yy <- y[binds]
	xxw <- xw[binds,]
	#xaux <- solve(t(xw) %*% xw)
	yyw <- yw[binds]

	#betast[r,] <- t(solve(t(xw[binds,]) %*% xw[binds,], t(xw[binds,]) %*% yw[binds]))
	betast[r,] <- solve(t(xxw) %*% xxw, t(xxw) %*% yyw)

	sigmast[r] <- sigma / b1 * mean(scaleweights[binds])

	#RBres <- t(cbind(t(beta), sigma)) + corrmat %*% t(cbind(betast[r,] - t(beta), sigmast[r] - sigma))
	RBres <- c(beta, sigma) + corrmat %*% c(betast[r,] - beta, sigmast[r] - sigma)

	
	betastR[r,] <- RBres[1:p]
	sigmastR[r] <- RBres[p+1]
	#Exact calculation of scale corresponding to approximated beta
	sigmastE[r] <- solveMscale(yy - x[binds,] %*% betastR[r,], RBres[p+1], 1e-5, b1, c1)

	#taust[r] <- sqrt(sigmast[r]^2 %*% mean(rhoOptfun((yy - x%*%t(betast[r,]))/sigmast[r], c2))/b2
	taustE[r] <- sqrt(sigmastE[r]^2 * mean(rhoOptfun((yy - x[binds,] %*% betastR[r,]) / 
sigmastE[r], c2)) / b2)

}

return(list(betastR = betastR, sigmastR = sigmastR, sigmastE = sigmastE, taustE = taustE, 
bootmatrix = bootmatrix))

}

#------------------------------------------------------------------------
#
#rhoOptfun <- function(x, c){
#
##Computes optimal rho function
#
#tmp <- x^2 / 2 / (3.25*c^2)
#tmp2 <- (1.792 - 0.972 * x^2 / c^2 + 0.432 * x^4 / c^4 - 0.052 * x^6 / c^6 +
# 0.002 * x^8 / c^8) / 3.25 #tmp[abs(x) > 2*c] <- tmp2[abs(x) > 2*c]
#tmp[abs(x) > 3*c] <- 1
#return(tmp)
#
#}
#
##---------------------------------------------------------------------
#
#psiOptfun <- function(x, c){
#
##Computes optimal rho function's first derivative
#
#tmp <- x / (3.25*c^2)
#tmp2 <- (-1.944 * x / c^2 + 1.728 * x^3 / c^4 - 0.312 * x^5 / c^6 + 
#0.016 * x^7 / c^8) / 3.25
#tmp[abs(x) > 2*c] <- tmp2[abs(x) > 2*c]
#tmp[abs(x) > 3*c] <- 0
#return(tmp)
#
#}
#
##-------------------------------------------------------------------------
#
#fw <- function(x, c){
#
## weight function = psi(x)/x
#
#tmp <- x^0 / (3.25*c^2)
#tmp2 <- (-1.944  / c^2 + 1.728 * x^2 / c^4 - 0.312 * x^4 / c^6 + 
#0.016 * x^6 / c^8) / 3.25
#tmp[abs(x) > 2*c] <- tmp2[abs(x) > 2*c]
#tmp[abs(x) > 3*c] <- 0
#tmp[abs(tmp) < 10e-10] <- 0
#return(tmp)
#
#}

#----------------------------------------------------------------------------

solveMscale <- function(x, initialscale=0, tol, b, c){

#M-estimator of scale using the biweight function

x <- as.matrix(x)

maxiter <- 100
n <- nrow(x)
if (initialscale==0){
	s <- median(abs(x))/0.6745
}else{
	s <- initialscale
}

rhoold <- mean(rhoOptfun(x/s,c)) - b
iter <- 0
while (abs(rhoold) > tol & iter < maxiter){
	delta <- rhoold / mean( psiOptfun(x/s,c) * (x/(s^2)))
	isqu <- 1; ok <- 0
	while (isqu < 30 & ok != 1){
		rhonew <- mean( rhoOptfun(x/(s+delta),c)) - b
		if( abs(rhonew) < abs(rhoold)){
			s <- s + delta; ok <- 1
		}else{
			delta <- delta/2; isqu <- isqu + 1
		}
	}
	if(isqu == 30){
		maxiter <- iter  #we tell it to stop, but we keep the iter for info
	}
	rhoold <- rhonew
	iter <- iter + 1
}
return(abs(s))

}

#---------------------------------------------------------------------------
#
#psiderOptfun <- function(x,c){
#
#tmp <- x^0 / (3.25*c^2)
#tmp2 <- (-1.944 / c^2 + 5.184 * x^2 / c^4 - 1.560 * x^4 / c^6 + 0.112 * x^6 / c^8) / 3.25
#tmp[abs(x) > 2*c] <- tmp2[abs(x) > 2*c]
#tmp[abs(x) > 3*c] <- 0
#return(tmp)
#
#}
