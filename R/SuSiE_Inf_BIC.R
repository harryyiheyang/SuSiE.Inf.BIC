#' Estimation of Causal Variants Using SuSiE-Inf and Selection of Optimal Number of Single Effects using BIC
#'
#' This function utilizes the SuSiE (Sum of Single Effects) method for fine-mapping to identify causal variants with significant effect sizes. It used BIC (Bayesian Informative Criterion) to select the optimal L. It employs REML (Restricted Maximum Likelihood) for the estimation of infinitesimal effects and for determining the optimal variance of these effects. Additionally, the function utilizes a score test to assess whether the variance of the infinitesimal effect is zero. The function has been updated to include parameters for controlling the precision of the estimation process and the option to estimate residual variance.
#'
#' @param z Z scores of GWAS effect sizes.
#' @param R LD (Linkage Disequilibrium) matrix of variants.
#' @param n Sample size of GWAS data.
#' @param Lvec A vector of numbers of single effects to consider, default is c(0:9).
#' @param pip.cred.thres A threshold to define the credible set, default is 0.95.
#' @param max.iter Maximum number of iterations for estimating the infinitesimal effect, default is 15.
#' @param max.eps The maximum epsilon for convergence in the iterative process, default is 0.001.
#' @param reml.iter Number of iterations for the inner loop of the REML, default is 10.
#' @param susie.iter Number of iterations for the inner loop of the SuSiE, default is 500.
#' @param score.test Perform score test of variance component or not, default to F.
#' @param pv.thres In score test, a threshold in score test of variance, default to 0.05. If the P-value of score test is larger than this threshold, the infinitesimal effect is removed.
#' @param pip.thres In score test, a threshold to determine which variants are considered causal, default is 0.5. This will only be used in score test as a fixed effect adjustment.
#' @param eigen.thres The threshold of eigenvalues for modelling the infinitesimal effect. Default is 1.
#' @param varinf.upper.boundary The upper boundary for the prior variance of infinitesimal effects, multiplied by var(y) to adapt to different locus variances. Default is 0.25.
#' @return A list containing the following elements:
#'   - \code{eta}: Linear predictor, which is the sum of the causal effect (\code{beta}) and the infinitesimal effect (\code{alpha}). It represents the total genetic effect of variants on the trait, combining the direct effects of causal variants and the background genetic effects.
#'   - \code{beta}: Causal effect. These are the effect estimates for variants identified with significant effect sizes through the SuSiE method, representing the direct impact of these variants on the trait.
#'   - \code{alpha}: Infinitesimal effect. This represents the aggregate effect of all other genetic variants in the background, excluding the identified causal variants. It can be considered as the unexplained genetic variability in the model.
#'   - \code{var.inf}: Variance of the infinitesimal effect. Estimated via REML, this reflects the degree of variation in the infinitesimal effect across the genetic variation.
#'   - \code{pv}: P-value from the score test. This value is used to assess whether the variance of the infinitesimal effect is significantly non-zero, serving as an important metric to test if the model adequately explains the genetic variability.
#'   - \code{pip}: Posterior Inclusion Probabilities for the estimates of causal effects. These are the probabilities calculated by the SuSiE method for each variant, assessing the likelihood that a variant is causal.
#'   - \code{fit.susie}: Output from the SuSiE method. This is the complete result object from the SuSiE method after performing fine-mapping, containing detailed information about the selected model and estimated parameters, allowing for further analysis and interpretation of the results.
#' @import CppMatrix
#' @importFrom susieR susie_rss
#' @importFrom Matrix bdiag
#' @export
#'
SuSiE_Inf_BIC <- function(z, R, n, Lvec = c(0:15), pip.cred.thres = 0.95, max.iter = 50, max.eps = 0.001, susie.iter = 500, reml.iter = 10, score.test = T, pip.thres = 0.5, pv.thres = 0.05, eigen.thres = 1, varinf.upper.boundary=0.25) {

var.inf=0.5
alpha=beta=z*0
Theta=matrixGeneralizedInverse(R)
m=length(z)
LD=R
varinf.upper.boundary=varinf.upper.boundary*mean(z*(Theta%*%z))

fiteigen=matrixEigen(R)
Umat=fiteigen$vectors
Dvec=fiteigen$values
Kthres=eigen_cumsum(Dvec,eigen.thres)
Umat=Umat[,1:Kthres]
Dvec=Dvec[1:Kthres]
LD2=matrixMultiply(Umat,t(Umat)*(Dvec^2))

Rinf=matrixMultiply(Umat,t(Umat)*(Dvec))

Bicvec=Bvar.inf=Bpv=c(1:length(Lvec))
Bbeta=matrix(0,m,length(Lvec))
Balpha=matrix(0,m,length(Lvec))
Bfit=list()
beta=beta1=0*z
for(i in 1:length(Lvec)){
error=1
iter=0
fit=NULL
alpha=z*0
while(error>max.eps&iter<max.iter){
beta1=beta
beta=beta*0
z1=z-matrixVectorMultiply(Rinf,alpha)
fit=susie_rss(z=z1,R=R,n=n,L=max(1,Lvec[i]),residual_variance=1,s_init=fit,estimate_prior_method="EM",max_iter=susie.iter)
beta=coef(fit)[-1]*sqrt(n)*(Lvec[i]!=0)
causal.cs=group.pip.filter(pip.summary=summary(fit)$var,pip.thres.cred=pip.cred.thres)
pip.alive=causal.cs$ind.keep
beta[-pip.alive]=0

res.upsilon=z-matrixVectorMultiply(R,beta)
outcome=matrixVectorMultiply(t(Umat),res.upsilon)
for(j in 1:reml.iter){
Hinv=1/(Dvec+1/var.inf)
alpha=matrixVectorMultiply(Umat,outcome*Hinv)
for(jj in 1:3){
df=sum(Hinv)
var.inf=min((sum(alpha^2)+df)/m,varinf.upper.boundary)
}
}

if(iter>3){error=max(abs(beta-beta1))}
iter=iter+1
if(score.test==T){
indbeta=top_K_pip(summary(fit)$vars,top_K=1,pip.min.thres=0.01,xQTL.pip.thres=pip.thres)
test=inf.test(res.inf=res.upsilon,LD=R,LD2=LD2,Theta=Theta,A=R[,which(fit$pip>=pip.thres)])
pv=test$pv
alpha=alpha*(pv<pv.thres)
}else{
pv=1
}
}
df1=Lvec[i]
df2=sum(Dvec*Hinv)*(pv<pv.thres)
res=z-matrixVectorMultiply(R,beta+alpha)
Bicvec[i]=log(sum(res*matrixVectorMultiply(Theta,res)))+log(n)/n*(df1+df2)
Bbeta[,i]=beta
Balpha[,i]=alpha
Bvar.inf[i]=var.inf
Bpv[i]=pv
Bfit[[i]]=fit
}
istar=which.min(Bicvec)
beta=Bbeta[,istar]
alpha=Balpha[,istar]
pv=Bpv[istar]
var.inf=Bvar.inf[istar]
fit=Bfit[[istar]]
return(list(eta=(alpha+beta)/sqrt(n),beta=beta/sqrt(n),alpha=alpha/sqrt(n),var.inf=var.inf,pv=pv,fit.susie=fit,L=Lvec[istar],Bicvec=Bicvec,df.inf=df2))
}
