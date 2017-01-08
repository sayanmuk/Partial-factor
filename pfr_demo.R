#source('PartialFactorRegression.R')
# uncomment the above line the first time you run in order to 
# compile the partial factor regression function


p <- 15
n <- 40
k <- 5


h <- .5
B <- matrix(rnorm(p*k,0,h),p,k)
W <- 1 + abs(rt(k,3))
B <- B%*%diag(sort(W,decreasing=FALSE))


u <- 1/6
Psi <- sqrt(diag(B%*%t(B)))/u


Fscores <- matrix(rnorm(k*n),k,n)
X <- B%*%Fscores + diag(Psi)%*%matrix(rnorm(p*n),p,n)

eta1 <- 5
eta <- c(eta1,rep(0,k-1)) 


gamtrue <- eta%*%t(B)%*%solve(B%*%t(B)+diag(Psi^2))

print(gamtrue)

M <- t(B)%*%solve(B%*%t(B)+diag(Psi^2))%*%B

sigma <- 1/5

s <- sqrt(M[1,1]+sigma^2)
y <- gamtrue%*%X + s*rnorm(n)




pfr_fit <- pfr(y,X,tot = 3000,burn = 1000,saveall=1)

fit_hs <- bhs(t(X),t(y),RJ=FALSE,T=4000)

pfr_est <- pfr_fit$betaest
hs_est <- colMeans(fit_hs$beta)


plot(pfr_est,gamtrue,pch=20)
points(hs_est,gamtrue,pch=20,col='red')

sdx <- apply(X,1,sd)
sdy <- sd(y)


D <- data.frame(t(y),t(X))
fit_bfa <- bfa_gauss(~.,data=D,num.factor=6,nsim=3000,burn=1000,print.status=0)

coef.bfa <- function(j){
B <- fit_bfa$post.loadings[1:p+1,,j]
A <- woodbury(B,fit_bfa$post.sigma2[j,1:p+1])

theta <- fit_bfa$post.loadings[1,,j]

coef <- theta%*%t(B)%*%A
return(coef)}

bfa_coeffs <- sapply(1:dim(fit_bfa$post.loadings)[3],coef.bfa)
bfa_est <- rowMeans(bfa_coeffs)
bfa_est <- bfa_est*sdy/sdx

points(bfa_est,gamtrue,pch=20,col='cyan')


abline(0,1)

cat('Factor model estimation error:',sqrt(s^2 + mean((bfa_est - gamtrue)^2))/s,'\n')
cat('Horseshoe regression estimation error:',sqrt(s^2 + mean((hs_est - gamtrue)^2))/s,'\n')
cat('Partial factor regression estimation error:',sqrt(s^2 + mean((pfr_est - gamtrue)^2))/s,'\n')


