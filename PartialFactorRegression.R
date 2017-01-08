library(inline)
library(Rcpp)
library(RcppArmadillo)
library(corpcor)
library(Matrix)
src='


using namespace arma;


imat klist =as<imat>(klistS);
imat plist = as<imat>(plistS);
arma::mat B = Rcpp::as<arma::mat>(BS);
arma::mat X = Rcpp::as<arma::mat>(XS);
arma::mat Y = Rcpp::as<arma::mat>(YS);
arma::mat Ff = Rcpp::as<arma::mat>(FfS);
int iter = as<int>(iterS);
arma::mat Lambda = Rcpp::as<arma::mat>(LambdaS);
arma::mat Psi = Rcpp::as<arma::mat>(PsiS);
arma::mat Theta = Rcpp::as<arma::mat>(ThetaS);
arma::mat t = Rcpp::as<arma::mat>(tS);
arma::mat v = Rcpp::as<arma::mat>(vS);
double tau = as<double>(tauS);
//arma::mat betasave = Rcpp::as<arma::mat>(betasaveS);
double sigma = as<double>(sigmaS);
double eta = as<double>(etaS);
int burn = as<int>(BurnInS);
int k = as<int>(kS);
int n = as<int>(nS);
int p = as<int>(pS);
int FlagLamSparse = as<int>(FlagLamSparseS);
int FlagThetSparse = as<int>(FlagThetSparseS);
int FlagBSparse = as<int>(FlagBSparseS);
double omega = as<double>(omegaS);
double alphab = as<double>(alphabS);
double alphaLam = as<double>(alphaLamS);
alphaLam=1-alphaLam;
double alphatheta = as<double>(alphathetaS);
Environment stats("package:stats");
Function pgamma = stats["pgamma"];
Function qgamma = stats["qgamma"];
Function dnorm = stats["dnorm"];
Function rbeta = stats["rbeta"];
//Function list = base["list"];
int acc = as<int>(accS);
int i;
int saveall = as<int>(saveallS);
ivec indLamZ = ones<ivec>(p);
mat Xvecind = zeros(p,n);
mat V0 = zeros(1,p);
mat RANDU;
mat RANDN;
mat betasave;
cube Fsave;
cube Bsave;
mat Thetasave;
mat Psisave;
mat Lambdasave; 
mat tausave; 
mat etasave;
mat omegasave;
mat sigmasave;
mat vsave;
cube tsave;

if(saveall==1){
betasave = zeros(iter,p);
Fsave = cube(k,n,iter);
Bsave = cube(p,k,iter);
Psisave = zeros(p,iter);
Thetasave = zeros(iter , k);
Lambdasave = zeros(iter , p);
tausave = zeros(1, iter);
etasave = zeros(1, iter);
omegasave = zeros(1, iter);
sigmasave = zeros(1, iter);
vsave = zeros(iter, p);
tsave = cube(p+1,k,iter);
}else{
//Rcpp::Rcout << "intialize betasave" << std::endl;
betasave = zeros(iter-burn,p);
}

Rcpp::Rcout << "The number of factors is ";
Rcpp::Rcout<< k <<std::endl;

for(i=0;i<iter;i++){

arma::mat Bold = B;
arma::mat Lambdaold=Lambda;
//Rcpp::Rcout << "Sample B" << std::endl;

//Sample B


arma::mat Ff2 = Ff%Ff; //find B^2
arma::mat ssB = sum(Ff2,1);
arma::mat resVarB = pow(omega, 2) * Psi % Psi;

arma::mat priorVarB = pow(tau, 2) * t.submat(span(0,p-1),span(0,k-1));
arma::mat cB = 1/((1/resVarB) * trans(ssB) + 1/priorVarB);
if(FlagBSparse==1){
arma::mat probsB = zeros(p,2);

if(k==1){
arma::mat mB = cB%(X*trans(Ff))/resVarB;
probsB.col(0)=1/(2*PI*cB.col(k-1))%exp(-pow(mB,2)/(2*cB.col(k-1)));
probsB.col(1)=1/(2*PI*priorVarB);
mat indprobs = alphab/(alphab+(1-alphab)*probsB.col(1)/probsB.col(0));
mat indB = floor(randu(p,1)+indprobs);
B.col(k-1)=indB%(mB+pow(cB.col(k-1),.5)%randn(p,1));
};
if(k>1){
imat kperm = shuffle(klist,1);
for(int j=0;j<k;j++){
int g = kperm(0,j);
mat AB = B;
AB.shed_cols(g,g);
mat FB = Ff;
FB.shed_rows(g,g);
mat YstarB = X - AB*FB;
mat mB = cB.col(g)%(YstarB*trans(Ff.row(g)))/resVarB;

for(int q=0;q<p;q++){
probsB(q,0)=as<double>(dnorm(0,mB(q,0),sqrt(cB(q,g))));
probsB(q,1)=as<double>(dnorm(0,0,sqrt(priorVarB(q,g))));
}
mat indprobs = alphab/(alphab+(1-alphab)*probsB.col(1)/probsB.col(0));
mat indB = floor(randu(p,1)+(1-indprobs)); // Generates bernoulli with probabilities indprobs
B.col(g)=indB%(mB+pow(cB.col(g),.5)%randn(p,1));
}
}
}else{//B has been flagged as not sparse
if(k==1){
arma::mat mB = cB%(X*trans(Ff))/resVarB;
B.col(k-1)=(mB+pow(cB.col(k-1),.5)%randn(p,1));
};
if(k>1){
imat kperm = shuffle(klist,1);
for(int j=0;j<k;j++){
int g = kperm(0,j);
mat AB = B;
AB.shed_cols(g,g);
mat FB = Ff;
FB.shed_rows(g,g);
mat YstarB = X - AB*FB;
mat mB = cB.col(g)%(YstarB*trans(Ff.row(g)))/resVarB;
B.col(g)=(mB+pow(cB.col(g),.5)%randn(p,1));

}

//Now Bprop is sampled
}
}
//Now Bprop is sampled

// compute residuals and diagonal Psis
mat residpropB = Y - Theta*Ff - Lambda*diagmat(1/(Psi*omega))*(X-B*Ff);
mat residB = Y - Theta*Ff - Lambda*diagmat(1/(Psi*omega))*(X-Bold*Ff);
mat proprat = -sum(log(1/pow((2*PI*pow(sigma,2)),.5)*exp(-pow(residB,2)/(2*pow(sigma,2)))),0) +
sum(log(1/pow((2*PI*pow(sigma,2)),.5)*exp(-pow(residpropB,2)/(2*pow(sigma,2)))),0);
RANDU=randu(1,1);
double u=RANDU(0,0);
if(exp(proprat(0,0))<u){
B=Bold;
}else{
acc++;
}

//End Sampling B
//Start Sampling Lambda
//Rcpp::Rcout << "Start Lambda" << std::endl;

mat Xvec = diagmat(1/(Psi*omega))*(X-B*Ff);
mat Yvec = Y - Theta*Ff;
if(FlagLamSparse==1){
RANDU=randu(1,1);
if(RANDU(0,0)<.5){
arma::mat ssB = sum(Xvec%Xvec,1);
double resVar = sigma*sigma;
priorVarB = eta*eta*v%v;
cB = 1/((1/resVar)*trans(ssB)+1/priorVarB);
mat probsL = zeros(1,2);
imat pperm = shuffle(plist,1);
for(int j=0;j<p;j++){
int g = pperm(0,j);
mat AL=Lambda;
AL.shed_cols(g,g);
mat Xmat = Xvec;
Xmat.shed_rows(g,g);
mat Ystar = Yvec - AL * Xmat;
mat m = cB.col(g)%(Ystar *trans(Xvec.row(g)))/resVar;
probsL(0,0)=as<double>(dnorm(0,m(0,0),sqrt(cB(0,g))));
probsL(0,1)=as<double>(dnorm(0,0,sqrt(priorVarB(0,g))));
mat indprobs = alphaLam/(alphaLam+(1-alphaLam)*probsL.col(1)/probsL.col(0));
RANDU=randu(1,1);
mat indB = floor(RANDU(0,0)+(1-indprobs)); // Generates bernoulli with probabilities 1-indprobs
RANDN=randn(1,1);
Lambda.col(g) = indB%(m+(pow(cB.col(g),.5)*RANDN(0,0)));
}
}else{
int sumL = 0;
for(int j=0;j<p;j++){
if(std::abs(Lambda(0,j))>0){
indLamZ(sumL)=j;
sumL++;
}
}
if(sumL>0){
V0.zeros(1,sumL);
Xvecind.zeros(sumL,n);
for(int l=0;l<sumL;l++){
Xvecind.row(l)=Xvec.row(indLamZ(l));
V0.col(l)=eta*eta*pow(v.col(indLamZ(l)),2);
}
mat XvecInd = Xvecind.submat(0, 0, sumL-1, n-1); 
mat V0Ind = V0.submat(0, 0, 0, sumL-1);
if(sumL==1){
mat Vmat = 1/((1/(sigma*sigma))*sum(pow(Xvecind.row(0),2),1)+1/V0.col(0)); 
RANDN=randn(1,1);
Lambda.col(indLamZ(0)) = (1/(sigma*sigma))*Vmat*sum(Xvecind.row(0)%Yvec,1)+pow(Vmat,.5)*RANDN(0,0);
}else{
mat Vmat = inv(sympd(((1/pow(sigma,2))*XvecInd*trans(XvecInd))+diagmat(1/V0Ind)));
mat L = trans(chol(Vmat));
//Rcpp::Rcout<<L<<std::endl;

mat Lambdaind=(1/pow(sigma,2))*Vmat*XvecInd*trans(Yvec)+L*randn(sumL);
for(int l=0;l<sumL;l++){
Lambda(0,indLamZ(l))=Lambdaind(l);
}
}
}
}
}else{//Lambda has benn flagged as not sparse

V0 = eta*eta*pow(v,2);
mat Vmat = inv(sympd(((1/pow(sigma,2))*Xvec*trans(Xvec))+diagmat(1/V0)));
mat L = trans(chol(Vmat));
mat Lambda=trans((1/pow(sigma,2))*Vmat*Xvec*trans(Yvec)+L*randn(p));

}

//End Lambda
//Sample PSi
//Rcpp::Rcout<<"Sample Psi"<<std::endl;
mat mPsi = sum(pow(X-B*Ff,2),1)/omega;

// sampCauchyVar=function(mu =mPsi, sigma=Psi Px1,num=n,alpha=1/2)

mat uPsi = randu(p,1)%(1/(1+1/(Psi%Psi)));

mat ubPsi = pow(uPsi,-1)-1;
double aPsi = (n+1)/2;
mat bPsi = mPsi/2;
for(int j=0;j<p;j++){
Psi(j,0) = as<double>(pgamma(ubPsi(j,0),aPsi,bPsi(j,0)));
}
mat tempPsi = randu(p)%Psi;
for(int j=0;j<p;j++){
Psi(j,0) = std::max(tempPsi(j,0),pow(10,-8));
}
for(int j=0; j<p;j++){
Psi(j,0) = as<double>(qgamma(Psi(j,0),aPsi,bPsi(j,0)));
}
Psi=1/sqrt(Psi);


//Rcpp::Rcout<<"Sample theta"<<std::endl;


Yvec = Y - Lambda*diagmat(1/(Psi*omega))*(X-B*Ff);
//Theta= sampTheta(Yvec = Ymat,Ff = Xmat,Theta = B,sigma^2 = resVar,t(tau^2*t[p+1,]^2) = priorVar,1=p,k=k,0,1-alphathet = alpha,RNORMT[iter]=RNORM);
mat ssFf = sum(Ff%Ff,1);
double resVar = pow(sigma,2); //resVar is already calculated as sigma^2
mat priorVar = trans(tau*tau*pow(t.row(p),2));
mat cT = 1/((1/resVar * ssFf) + 1/priorVar);
if(FlagThetSparse==1){
if(k==1){
mat probsT = zeros(1,2);

mat m = cT%(Yvec * trans(Ff))/resVar;
probsT.col(0) = 1/(2*PI*cT.col(0))%exp(-pow(m,2)/(2*cT.col(k-1)));
probsT.col(1) = 1/(2*PI*priorVar.col(0));
double probs = (1-alphatheta)/((1-alphatheta)+(alphatheta)*probsT(0,0)/probsT(0,1));


RANDU=randu(1,1);
if(floor(RANDU(0,0)+1-probs)>0){
RANDN=randn(1,1);


Theta.col(k-1) = (m+pow(cT,.5)*RANDN(0,0));
}else{


Theta(0,k-1) = 0;
}

}else{//k>1, 
imat kperm = shuffle(klist,1);

for(int j=0;j<k;j++){
int g = kperm(0,j);
mat AT = Theta;
AT.shed_cols(g,g);
mat XmatT = Ff;
XmatT.shed_rows(g,g);
mat YstarT = Yvec - AT*XmatT;
cT.row(g);


mat m = cT.row(g)*(YstarT * trans(Ff.row(g)))/resVar;
mat probsT = zeros(1,2);
probsT.col(0)=1/pow((2*PI*cT.row(g)),.5)%exp(-pow(m,2)/(2*cT.row(g)));
probsT.col(1)=1/pow((2*PI*priorVar.row(g)),.5);

double indprobsT = (1-alphatheta)/(1-alphatheta+alphatheta*probsT(0,1)/probsT(0,0));
RANDU=randu(1,1);
if(floor(RANDU(0,0)+1-indprobsT)>0){
RANDN=randn(1,1);
Theta.col(g) = m + pow(cT.row(g),.5)*RANDN(0,0);
}else{
Theta(0,g) = 0;
}
}
}
}else{//Theta flagged as not sparse
if(k==1){

mat m = cT%(Yvec * trans(Ff))/resVar;

RANDN=randn(1,1);


Theta.col(k-1) = (m+pow(cT,.5)*RANDN(0,0));


}else{//k>1, 
imat kperm = shuffle(klist,1);

for(int j=0;j<k;j++){
int g = kperm(0,j);
mat AT = Theta;
AT.shed_cols(g,g);
mat XmatT = Ff;
XmatT.shed_rows(g,g);
mat YstarT = Yvec - AT*XmatT;
cT.row(g);


mat m = cT.row(g)*(YstarT * trans(Ff.row(g)))/resVar;

RANDN=randn(1,1);
Theta.col(g) = m + pow(cT.row(g),.5)*RANDN(0,0);

}
}
}
//End Theta
//Begin Sigma
mat residSig = Y - Theta * Ff - Lambda * diagmat(1/(Psi*omega))*(X-B*Ff);
double mSig = accu(pow(residSig,2));
double etasig = 1/pow(sigma,2);
RANDU=randu(1,1);
double usig = RANDU(0,0)/(1+etasig);
double ubsig = 1/usig - 1;
double asig = (n+1)/2;
double bsig = mSig/2;

double ub2sig = as<double>(pgamma(ubsig, asig, bsig));
RANDU=randu(1,1);
double tempsig = RANDU(0,0)*ub2sig;
ub2sig = std::max(tempsig, pow(10,-10));
etasig = as<double>(qgamma(ub2sig, asig, bsig));
sigma = 1/pow(etasig, .5);

//Begin Eta
if(FlagLamSparse==1){
int sumL=0;
double meta = 0;
for(int j=0;j<p;j++){
if(std::abs(Lambda(0,j))>0){
sumL++;
meta = meta + pow(Lambda(0,j),2)/pow(v(0,j),2);
}
}
if(sumL<1){
eta = pow(1/as<double>(rbeta(1,.5,1)),.5);
}else{
double eat = 1/pow(eta,2);
RANDU=randu(1,1);
double ueta = RANDU(0,0)*pow(1+eat,-1.5);
double ubeta = pow(ueta,-2/3.0)-1;
double aeta = (sumL+1)/2;
double beta = meta/2;
double ub2eta = as<double>(pgamma(ubeta,aeta,beta));
RANDU=randu(1,1);
double tempeta = RANDU(0,0)*ub2eta;
double u2 = std::max(tempeta,pow(10,-10));
eat = as<double>(qgamma(u2,aeta,beta));
eta = pow(eat,-.5);
}
}else{//Lambda has been flagged as not sparse
double meta = 0;
for(int j=0;j<p;j++){
meta = meta + pow(Lambda(0,j),2)/pow(v(0,j),2);
}
double eat = 1/pow(eta,2);
RANDU=randu(1,1);
double ueta = RANDU(0,0)*pow(1+eat,-1.5);
double ubeta = pow(ueta,-2/3.0)-1;
double aeta = (p+1)/2;
double beta = meta/2;
double ub2eta = as<double>(pgamma(ubeta,aeta,beta));
RANDU=randu(1,1);
double tempeta = RANDU(0,0)*ub2eta;
double u2 = std::max(tempeta,pow(10,-10));
eat = as<double>(qgamma(u2,aeta,beta));
eta = pow(eat,-.5);
}
//End Eta
//Begin Samp V

if(FlagLamSparse==1){
for(int j=0;j<p;j++){
if(std::abs(Lambda(0,j)>0)){
double mv = pow(Lambda(0,j),2)/pow(eta,2);
double veta = 1/pow(v(0,j),2);
RANDU=randu(1,1);
double uv = RANDU(0,0)*pow(1+veta,-1.5);
double ubv = pow(uv,-2.0/3.0)-1;
double bv = mv/2;
double ub2v = as<double>(pgamma(ubv,1,bv));
RANDU=randu(1,1);
double tempv = RANDU(0,0)*ub2v;
veta = as<double>(qgamma(std::max(tempv, pow(10,-18)),1,bv));
v(0,j)=1/pow(veta,.5);
}
}
}else{//Lambda has been flagged as not sparse
for(int j=0;j<p;j++){
double mv = pow(Lambda(0,j),2)/pow(eta,2);
double veta = 1/pow(v(0,j),2);
RANDU=randu(1,1);
double uv = RANDU(0,0)*pow(1+veta,-1.5);
double ubv = pow(uv,-2.0/3.0)-1;
double bv = mv/2;
double ub2v = as<double>(pgamma(ubv,1,bv));
RANDU=randu(1,1);
double tempv = RANDU(0,0)*ub2v;
veta = as<double>(qgamma(std::max(tempv, pow(10,-18)),1,bv));
v(0,j)=1/pow(veta,.5);
}
}
//End Samp V
//Begin samp t
//Bmat=rbind(B/tau,Theta/tau);
// m=Bmat^2;
//ind=which(Bmat != 0,arr.ind=TRUE);
//t[ind]=sampCauchyVar(m[ind] = mu,t[ind] =sigma ,1 = num ,1 =alpha);, n=1 , alpha =1.5
mat Bmat = zeros(p+1,k);
Bmat.submat( 0, 0, p-1, k-1 ) = B;
Bmat.row(p)= Theta;
mat mt = pow(Bmat/tau,2);
for(int r=0;r<p+1;r++){
for(int c=0; c<k;c++){
if(mt(r,c)>0){
double veta = 1/pow(t(r,c),2);
RANDU=randu(1,1);
double uv = RANDU(0,0)*pow(1+veta,-1.5);
double ubv = pow(uv,-2/3.0)-1;
double bv = mt(r,c)/2;
double ub2v = as<double>(pgamma(ubv,1,bv));
RANDU=randu(1,1);
double tempv = RANDU(0,0)*ub2v;
veta = as<double>(qgamma(std::max(tempv, pow(10,-10)),1,bv));
t(r,c)=1/pow(veta,.5); 
}//close if
}//close for c
}//close for r

//Rcpp::Rcout<<"End Samp t"<<std::endl;
//Rcpp::Rcout<<"Begin Samp tau"<<std::endl;
//End Samp t
//Begin samp tau
if( (FlagBSparse==1) || (FlagThetSparse==1) ){
Bmat = Bmat/t;
//sum over the square of all nonzero values and count the number of nonzero values for the n of sampcauchy
double mtau = 0;
int sumL = 0;
for(int r=0;r<p+1;r++){
for(int c=0; c<k;c++){
if(std::abs(Bmat(r,c))>0){
mtau=mtau+pow(Bmat(r,c),2);
sumL++;
}//end if Bmat(r,c)>0
}//end row index
}//end col index
//Rcpp::Rcout<<"mTau is ";
//Rcpp::Rcout<<mtau<<std::endl;

if(sumL>0){
//sampCauchyVar(mtau = mu,tau =sigma,bn = num,1)   alpha=1.5, n=1;
double veta = 1/pow(tau,2);
RANDU=randu(1,1);
double uv = RANDU(0,0)*pow(1+veta,-1.5);
double ubv = pow(uv,-2/3.0)-1;
double bv = mtau/2;
double ub2v = as<double>(pgamma(ubv,(sumL)/2,bv));
RANDU=randu(1,1);
double tempv = RANDU(0,0)*ub2v;
veta = as<double>(qgamma(std::max(tempv, pow(10,-10)),(sumL)/2,bv));
tau=1/pow(veta,.5); 
}else{
tau = pow(1/as<double>(rbeta(1,.5,1))-1,.5);
}
}else{//Theta and B have been flagged as not sparse
Bmat = Bmat/t;
double mtau = 0;
for(int r=0;r<p+1;r++){
for(int c=0; c<k;c++){
mtau=mtau+pow(Bmat(r,c),2);
}//end row index
}//end col index
double veta = 1/pow(tau,2);
RANDU=randu(1,1);
double uv = RANDU(0,0)*pow(1+veta,-1.5);
double ubv = pow(uv,-2/3.0)-1;
double bv = mtau/2;
double ub2v = as<double>(pgamma(ubv,(p*k)/2,bv));
RANDU=randu(1,1);
double tempv = RANDU(0,0)*ub2v;
veta = as<double>(qgamma(std::max(tempv, pow(10,-10)),(p*k)/2,bv));
tau=1/pow(veta,.5); 

}


//Rcpp::Rcout<<"End Samp Tau"<<std::endl;
//Rcpp::Rcout<<"Begin Samp F"<<std::endl;
//End Samp tau
//Begin samp F
// Fsamp=sampF(B=B,X=X,omega^2*Psi^2 =Psi,sigma^2 = s,Theta =C,Lambda = Lambda,t(Y) =Y,k=k,n=n);
//  Ff=Fsamp[[1]];
//  invsigx=Fsamp[[2]];
//  betatemp=Fsamp[[3]];
mat PsiF = pow(Psi, 2)*pow(omega,2);
mat sigxF = B*trans(B) + diagmat(PsiF);
mat DF = diagmat(1/PsiF);
mat AF = trans(B)*DF*B+eye(k,k);
mat invsigx = DF - DF*B*solve(AF, trans(B))*DF;
mat VF = Lambda*diagmat(pow(PsiF,.5))+Theta*trans(B);
mat betal = VF*invsigx;
mat omegaF = (pow(sigma,2)) + Lambda*trans(Lambda) + Theta*trans(Theta);
mat gF = omegaF - betal*trans(VF);
mat invmat = invsigx + trans(betal)*betal/gF(0,0);
mat invmatF = zeros(p+1,p+1);
invmatF.submat(0,0,p-1,p-1) = invmat;
invmatF.submat(0,p,p-1,p) = -trans(betal)/gF(0,0);
invmatF.submat(p,0,p,p-1) = -betal/gF(0,0);
invmatF(p,p) = 1/gF(0,0);
mat XYF = zeros(p+1,n);
XYF.submat(0,0,p-1,n-1)=X;
XYF.row(p)=Y;
mat BCF = zeros(k,p+1);
BCF.submat(0,0,k-1,p-1)=trans(B);
BCF.col(p)=trans(Theta);
mat m2F = BCF*invmatF*XYF;
mat Smat = BCF*invmatF*trans(BCF);
Smat = eye(k,k) - Smat;
if(det(Smat)>0.001){
mat cLF = trans(chol(Smat));
Ff = m2F + cLF*randn(k,n);
}else{
mat cLF = trans(chol(Smat+.001*eye(k,k)));
Ff = m2F + cLF*randn(k,n);
}

//end samp F
//Rcpp::Rcout<<"Begin Samp F"<<std::endl;

//find betavec



mat betavec = Theta*trans(B)*invsigx+Lambda*diagmat(1/Psi)-Lambda*diagmat(1/Psi)*B*trans(B)*invsigx;


if(saveall == 1){
betasave.row(i) = betavec;
Fsave.slice(i)=Ff;
Thetasave.row(i) = Theta;
Psisave.col(i)=Psi;
Bsave.slice(i) = B;
Lambdasave.row(i)=Lambda;
tausave(0,i) = tau;
etasave(0,i) = eta;
omegasave(0,i)=omega;
sigmasave(0,i)=sigma;
tsave.slice(i) = t;
vsave.row(i) =v;

}else{
if(i>(burn-1)){
betasave.row(i-burn)=betavec;

}
}




}


if(saveall==0){

return wrap(betasave);
}else{
return List::create(Named("betasave") = betasave,
Named("Fsave")=Fsave,
Named("Thetasave")=Thetasave,
Named("Bsave") =Bsave,
Named("Psisave")=Psisave,
Named("Lambdasave") =Lambdasave,
Named("tausave")=tausave,
Named("etasave")=etasave,
Named("omegasave")=omegasave,
Named("sigmasave")=sigmasave,
Named("tsave")=tsave,
Named("vsave")=vsave);

}


'

PFRoption<-cxxfunction(signature(XS="double", YS="double", BS="double", 
                           FfS="double",
                           iterS="integer",LambdaS="double",PsiS="double",
                           ThetaS="double",sigmaS="double",
                           tS="double",vS="double",tauS="double",
                           BurnInS="integer",
                           kS="integer", nS="integer",pS="integer",
                           FlagLamSparseS="integer", FlagThetSparseS="integer", FlagBSparseS="integer",
                           omegaS="double", alphabS="double",
                           alphaLamS="double", alphathetaS="double", etaS="double",
                           accS="int",klistS="integer",plistS="integer",saveallS="integer"),
                 plugin="RcppArmadillo",
                 body=src)




pfr <-function(Y,X,tot=10000,burn=1000,k=-1,saveall=0,LamSparse=1,ThetaSparse=0,BSparse=0,sigma=-1){
  #WARNING: Require corpcor and Matrix packages
  # set number of factors to use 
  # if set to -1, the model will use the number of factors which has the
  # largest ratio between subsequent eigenvalues or p or n, whichever is
  # smaller
  if(is.null(dim(Y))){
    Y=as.matrix(Y)
  }
  
  
  if(dim(Y)[1]>1){
    Y=t(Y);
  }
  n=dim(Y)[2];
  
  
  
  #
  # Initialize PFR
  #
  #   this code picks initial values for parameters in the sparse partial
  #   factor regression model and preallocates memory for storing MCMC
  #   samples
  #
  #
  
  h=1;
  acc=0;
  if(dim(X)[2]!=n){
    if(dim(X)[1]==n){
      X=t(X);
    }else{
      print("X does not match Y")
    }
  }
  
  
  
  p=dim(X)[1];
  
  ds=apply(X,1,sd);
  X=diag(1/ds)%*%X;
  stdy=apply(Y,1,sd);
  Y=Y/stdy;
  
  tmp=svd(X); 
  S=tmp$u;
  D=diag(tmp$d);
  U=tmp$v;
  
  if(k==-1){
    p1=dim(D)[1];
    
    d=tmp$d/sum(tmp$d);
    d=d[1:(p1-1)]/d[2:p1];
    val=max(d);
    k=which.max(d);
    k=min(k,3)
  }
  
  
  
  
  
  alphab=.5;
  alphaLam=.1;
  alphathet=.9
  if (k > 1){
    B=S[,1:k]%*%diag(sqrt(diag(D[1:k,1:k])));
    tmp=Matrix::expand(lu(B));
    B=as.matrix(tmp$L);
    A=as.matrix(tmp$U);
    C=as.matrix(tmp$P);
    
    H1=diag(1/sqrt(abs(diag(A))));
    H2=diag(sqrt(abs(diag(A))));
    B=B%*%H2;
    Ff=H1%*%A%*%diag(sqrt(diag(D[1:k,1:k])))%*%t(U[,1:k]);
    
  }else{B = as.matrix(sqrt(D[1,1])*S[,1])
        H1 = sqrt(sum(B^2))
        A = t(B/H1)   
        Ff=as.matrix(H1*t(U[,1])*D[1,1]);
  }
  #print(dim(Ff))
  
  
  
  
  Ff=0.5*Ff+0.5*matrix(rnorm(dim(Ff)[1]*dim(Ff)[2],0,1),nrow=dim(Ff)[1],ncol=dim(Ff)[2]);
  tau = 1;
  if(sigma == -1){
    sigma = 0.1;
  }
  eta = 2;
  
  t = matrix(1,p+1,k);
  Psi = as.numeric(0.1*sqrt(t(apply(X-B%*%Ff,1,var))));
  
  omega=1;
  
  v=matrix(1,1,p);
  
  Theta=t(pseudoinverse(Ff%*%t(Ff))%*%Ff%*%t(Y))
  
  sigx=B%*%t(B)+diag(Psi);
  D=diag(1/Psi);
  A=t(B)%*%D%*%B+diag(k);
  invsigx=D-D%*%B%*%solve(A)%*%t(B)%*%D;
  
  
  Lambda=pseudoinverse(X%*%t(X))%*%X%*%t(Y)-t(Theta%*%t(B)%*%invsigx);
  Lambda=pseudoinverse(Lambda)%*%(diag(p)-B%*%t(B)%*%invsigx);
  Lambda=Lambda*.1
  #debug Lambda
  #print(Lambda)
  ##
  #   preallocate memory
  ##
  
  
  
  
  tup = matrix(1,p+1,p+1);
  tup[lower.tri(tup)]=0;
  tup=tup-diag(p+1);
  tup=tup[,1:k];
  
  
  
  klist=matrix(0:(k-1),1,k);
  plist=matrix(0:(p-1),1,p);
  
  
  
  
  if(saveall==0){
    betasave= PFRoption(X,Y,B,Ff,tot,Lambda,
                        matrix(Psi,p,1),Theta,sigma,t,v,tau, burn,k,n,p,LamSparse,ThetaSparse,BSparse,omega,alphab,alphaLam,alphathet,eta,acc,
                        klist,plist,saveall)
    #print('Exit PFR')
    betasave=t(betasave);
    betaest=(apply(betasave,1,mean)/ds)*stdy;
    piest=apply(betasave!=0,1,mean);
    
    output =list(betaest,piest,betasave) 
    names(output) = c('betaest','piest','betasave')
    return(output)
  }else{
    pfrmcmc = PFRoption(X,Y,B,Ff,tot,Lambda,
                        matrix(Psi,p,1),Theta,sigma,t,v,tau, burn,k,n,p,LamSparse, ThetaSparse, BSparse,omega,alphab,alphaLam,alphathet,eta,acc,
                        klist,plist,saveall)
    #print('Exit PFR')
    betasave=pfrmcmc$betasave;
    betasave=t(betasave);
    betaest=(apply(betasave[,(burn+1):tot],1,mean)/ds)*stdy;
    piest=apply(betasave!=0,1,mean);
    
    output = list(betaest,piest,betasave,pfrmcmc[2],pfrmcmc[3],pfrmcmc[4],pfrmcmc[5],pfrmcmc[6],pfrmcmc[7],
                  pfrmcmc[8],pfrmcmc[9],pfrmcmc[10],pfrmcmc[11],pfrmcmc[12]);
    names(output) = c('betaest','piest','betasave','Fsave','Thetasave','Bsave','Psisave','Lambdasave','tausave',
                      'etasave','omegasave','sigmasave','tsave','vsave')
    output$Fsave=output$Fsave$Fsave;
    output$Thetasave=output$Thetasave$Thetasave;
    output$Bsave=output$Bsave$Bsave;
    output$Psisave=output$Psisave$Psisave;
    output$Lambdasave=output$Lambdasave$Lambdasave;
    output$tausave=output$tausave$tausave;
    output$etasave=output$etasave$etasave;
    output$omegasave=output$omegasave$omegasave;
    output$sigmasave=output$sigmasave$sigmasave;
    output$tsave=output$tsave$tsave;
    output$vsave=output$vsave$vsave;
    return(output)
    
    
  }
  
  
}

