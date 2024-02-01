# SCENT
Gaussian graphical models are widely used to study the dependence structure among variables. When samples are obtained from multiple conditions or populations, joint analysis of multiple graphical models is desired due to its capacity to borrow strength across populations. In this package, the function simultaneously clusters and estimates graphical models from multiple populations. Precision matrices from different populations are uniquely organized as a three-way tensor array, and a low-rank sparse model is used for joint population clustering and network estimation. 

This package implements an ADMM algorithm to conduct Simultaneous Clustering and Estimation of Networks.

To install the package, run the following code in R:

library(devtools)
install_github("reagan0323/SCENT")

A toy example of the SCENT function:

library(SCENT)
SS=rWishart(6,20,0.01*diag(20))
out=SCENT(SS,r=2,nn=rep(100,6),
              rho=50)
par(mfrow=c(1,2))
plot(out$primal)
plot(out$dual)
out$U
