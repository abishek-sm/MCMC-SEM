library(msm)
library(MCMCpack)
library(mnormt)
m<-10
n<-300
p<-15

x<- matrix(runif(n*m,0,1),n,m) #n*m
x[x<0.3]<--1
x[x>0.66]<-1
x[x<=0.66 & x>0.3]<-0
mu<-sample(0:1,p,replace=T)#1*p
mu<-matrix(rep(mu,each=n),c(n,p))

b<-matrix(runif(m*p,0,1),m,p)#m*p
b<-matrix(0,m,p)
b<-matrix(sample(c(0,1),m*p,replace=T),m,p,byrow=T)

origb<-b
r<-diag(p)
e<-rmnorm(n,0,r)#n*p

lam<-matrix(sample(0:1,p*p,replace=T),p,p,byrow=T)#p*p matrix

for(i in 1:p)
for(j in 1:p)
if(i>=j)lam[i,j]=0
diag(lam)<-1

while(det(lam)==0) 
 { lam<-matrix(sample(0:1,p*p,replace=T),p,p,byrow=T)#p*p matrix 
   for(i in 1:p)
   for(j in 1:p)
   if(i>=j)lam[i,j]=0
   diag(lam)<-1
 }


origlam<-lam
y<- (mu + x%*%b + e)%*%solve(lam)
mujkp<-1
mujkm<-1
sigbkp<-sample(1,p,replace=T)#vector of p values
sigbkm<-sample(1,p,replace=T)#vector of p values
sigjkp<-2
sigjkm<-2
pjkp<-0.4
pjkm<-1
pbkp<-runif(p,0,1)#vector of p values
pbkm<-runif(p,0,1)#vector of p values
thetabk<-0.5
phibk<-2
psibk<-2
thetabkp<-0.5
thetabkm<-0.5
phibkp<-2
phibkm<-2
muklp<-1
siglap<-1
sigklp<-1
muklm <-1
siglam<-1
sigklm<-1
plap<-0.2
plam<-0.3
thetala<-1
phila<-1
psila<-1
thetalap<-0.5
thetalam<-0.5
philap<-2
philam<-2
nbkp<-numeric(p)
nbkm<-numeric(p)
b_agg<-0
lam_agg<-0
llig<-0
llikeys<-0
lam<-matrix(1,p,p)#initialising lambdas to 1
b<-matrix(0,m,p)#initialising b to 0s
for(i in 1:p)
for(j in 1:p)
if(i>j)lam[i,j]=0
lam<-diag(p)

for( wq in 1:50000)
{
#sample mu for the entire thing.
me<-colSums(y%*%lam -x%*%b)
mu<-rmnorm(1,me/n,r/n)

yilam<-y%*%lam

for( j in 1:m){
for( k in 1:p) {#for loop with j and k
temp<-0
for(i in 1:n)
temp<- temp+ (x[i,j]*(yilam[i,k]-mu[k] - (x[i,]%*%b[,k]-x[i,j]*b[j,k])))

mujkp <- sigbkp[k]*temp/(sigbkp[k]*sum(x[,j]^2)+ r[k,k]^-1)
sigjkp<-r[k,k]^-1*sigbkp[k]/(sigbkp[k]*sum(x[,j]^2 + r[k,k]^-1))

mujkm<- sigbkm[k]*temp/(sigbkm[k]*sum(x[,j]^2)+r[k,k]^-1)
sigjkm<-r[k,k]^-1*sigbkm[k]/(sigbkm[k]*sum(x[,j]^2 + r[k,k]^-1))

temp1<- 2*pbkp[k]*sqrt(sigjkp/sigbkp[k])*pnorm(mujkp/sqrt(sigjkp))
temp2<- 2*pbkm[k]*sqrt(sigjkm/sigbkm[k])*pnorm(-mujkm/sqrt(sigjkm))

pjkp <- temp1*exp((mujkp^2/(2*sigjkp))-(mujkm^2/(2*sigjkm)))/((1 -pbkp[k] -pbkm[k])*exp(-mujkm^2/(2*sigjkm)) +temp1*exp((mujkp^2/(2*sigjkp))-(mujkm^2/(2*sigjkm))) +temp2)
pjkm <- temp2*exp((mujkm^2/(2*sigjkm))-(mujkp^2/(2*sigjkp)))/((1 - pbkp[k] -pbkm[k])*exp(-mujkp^2/(2*sigjkp)) +temp1 +temp2*exp((mujkm^2/(2*sigjkm))-(mujkp^2/(2*sigjkp))))

urv<-runif(1, min=0, max=1)
  if(urv <(1-pjkp-pjkm))
     b[j,k]<-0
  else if( urv>= (1-pjkm-pjkp) & urv< (1-pjkp))
     b[j,k]<- rtnorm(1, mean=mujkm, sd=sqrt(sigjkm), lower=-Inf, upper=0)
  else
     b[j,k]<-rtnorm(1, mean=mujkp, sd=sqrt(sigjkp), lower=0, upper=Inf)

}
}

b_agg<-rbind(b_agg,b) 
#6.3 sampling beta hyperparameters

#for each k
for(k in 1:p)
{
nbkp[k]<-length(which(b[,k]>0))
nbkm[k]<-length(which(b[,k]<0))

x1<-rdirichlet(1,c(thetabk+nbkp,phibk+nbkm,psibk+m-nbkp-nbkm))
pbkp[k]<-x1[1]
pbkm[k]<-x1[2]

 temp1<-0
 temp2<-0
 for(j in 1:m)
 { temp1<- temp1 + (ifelse(b[j,k]>0,b[j,k],0)^2)
   temp2<- temp2 + (ifelse(b[j,k]<0,b[j,k],0)^2)
 }

sigbkp[k]<-1/rgamma(1,shape=(thetabkp+nbkp[k]/2),scale=(phibkp^-1 + temp1/2)^-1)
sigbkm[k]<-1/rgamma(1,shape=(thetabkm+nbkm[k]/2),scale=(phibkm^-1 + temp2/2)^-1)
}

#6.4 sampling r


for( k in 1:p)
{
  temp<-0
  for(i in 1:n)
  temp<-  temp + (yilam[i,k]-mu[k]-x[i,]%*%b[,k])^2
  r[k,k]<-rgamma(1,shape=n/2,scale=2/temp)
}


#6.5 sampling lambda
for (k in 1:p){
for(l in 1:p) {
if(k>=l) next
temp<-0
for(i in 1:n)
{ temp1<-0
for (r1 in 1:p)
 {
  temp1 <- temp1+ lam[k,r1]*y[i,r1] 
  }
  temp1<- temp1-lam[k,k]*y[i,k]-lam[k,l]*y[i,l]
  
 temp<-temp+(y[i,l]*(-y[i,k]-temp1 + mu[k] + x[i,]%*%b[,k]))
} 
 
muklp<- (siglap*temp)/(siglap*sum(y[,l]^2)+r[k,k]^-1)
sigklp<- r[k,k]*siglap/(siglap*sum(y[,l]^2)+r[k,k]^-1)
muklm<- (siglam*temp)/(siglam*sum(y[,l]^2)+r[k,k]^-1)
sigklm<- r[k,k]*siglam/(siglam*sum(y[,l]^2)+r[k,k]^-1)


temp3<- 2*plap*sqrt(sigklp/siglap)*pnorm(muklp/sqrt(sigklp))
temp4<- 2*plam*sqrt(sigklm/siglam)*pnorm(-muklm/sqrt(sigklm))

pklp<- temp3*exp(muklp^2/(2*sigklp)-(muklm^2/(2*sigklm)))/((1 - plam-plap)*exp(-muklm^2/(2*sigklm)) +temp3*exp(muklp^2/(2*sigklp)-(muklm^2/(2*sigklm))) + temp4)
pklm<- temp4*exp(muklm^2/(2*sigklm)-(muklp^2/(2*sigklp)))/((1 - plam-plap)*exp(-muklp^2/(2*sigklp)) +temp3 + temp4*exp(muklm^2/(2*sigklm)-(muklp^2/(2*sigklp))))

lamp<-lam
urv<-runif(1, min=0, max=1)

  if(urv <(1-pklp-pklm))
     lam[k,l]<-0
  else if( urv>= (1-pklm-pklp) & urv< (1-pklp))
     lam[k,l]<- rtnorm(1, mean=muklm, sd=sqrt(sigklm), lower=-Inf, upper=0)
  else
     lam[k,l]<-rtnorm(1, mean=muklp, sd=sqrt(sigklp), lower=0, upper=Inf)

if(det(lam)==0 | det(lamp)==0)
{}
else{     


mu1<-matrix(rep(mu,each=n),c(n,p))
prevval<-lamp[k,l]

logry<- (-1/2*(y %*%lam - mu1 - x%*%b)%*%solve(r)%*%t(y%*%lam-mu1-x%*%b))
logry1<- (-1/2*(y %*%lamp - mu1 - x%*%b)%*%solve(r)%*%t(y%*%lamp-mu1-x%*%b))

temp34<-sum(diag(logry))
temp45<-sum(diag(logry1))

bpl<-1
for(k11 in 1:p)
for(l11 in 1:p)
 if(k11==l11) next
 { if(lam[k11,l11]==0)
     bpl<-bpl*(1-plap-plam)
     else if(lam[k11,l11]>0)
     bpl<-bpl*(2*plap/sqrt(2*pi*siglap)*exp(-lam[k11,l11]^2/(2*siglap)))
     else
     bpl<-bpl*(2*plam/sqrt(2*pi*siglam)*exp(-lam[k11,l11]^2/(2*siglam)))
 }

bplp<-1
for(k11 in 1:p)
for(l11 in 1:p)
 if(k11==l11) next
 { if(lamp[k11,l11]==0)
     bplp<-bplp*(1-plap-plam)
     else if(lamp[k11,l11]>0)
     bplp<-bplp*(2*plap/sqrt(2*pi*siglap)*exp(-lamp[k11,l11]^2/(2*siglap)))
     else
     bplp<-bplp*(2*plam/sqrt(2*pi*siglam)*exp(-lamp[k11,l11]^2/(2*siglam)))
 }
 
logalpha1 <- (n*log(abs(det(lam))))+temp34+ log(bpl)-(n*log(abs(det(lamp))))-temp45-log(bplp)

if(exp(logalpha1)<1)
{ urv1<-runif(1, min=0, max=1)
  if(urv1>exp(logalpha1)) lam[k,l]<-lamp[k,l] 
}
}
}
}

#sampling lambda hyperparameters
nlap<-length(which(l>0))
nlam<-length(which(l<0))

x1<-rdirichlet(1,c(thetala+nlap,phila+nlam,psila+p-nlap-nlam))
plap<-x1[1]
plam<-x1[2]

 temp1<-0
 temp2<-0
 for(k in 1:p)
{ for( l in 1:p)
 { temp1<- temp1 + (ifelse(lam[k,l]>0,lam[k,l],0)^2)
   temp2<- temp2 + (ifelse(lam[k,l]<0,lam[k,l],0)^2)
 }
 }
 
siglap<-1/rgamma(1,shape=(thetalap+nlap/2),scale=(philap^-1 + temp1/2)^-1)
siglam<-1/rgamma(1,shape=(thetalam+nlam/2),scale=(philam^-1 + temp2/2)^-1)
lam_agg<-rbind(lam_agg,lam)

if(wq%%250 == 0 )
print(wq)

}

b_agg<-b_agg[-1,]
lam_agg<-lam_agg[-1,]
  
 beta_matrix<-matrix(nrow=m,ncol=p)
 for( j in 1:m) {
 for(k in 1:p) {
 if(j==m) j<-0
 sto<-0
 for(i in 1:nrow(b_agg)){
 if(i%%m ==j)
 sto<-rbind(sto,b_agg[i,k])
 } 
 sto<-sto[-1]
 if(j==0) j<-m
 beta_matrix[j,k]<-mean(sto[30000:50000],na.rm=T)}}
  
 
 lam_matrix<-matrix(nrow=p,ncol=p)
 for( j in 1:p) {
 for(k in 1:p) {
 if(j==p) j<-0
 sto1<-0
 for(i in 1:nrow(lam_agg)){
 if(i%%p ==j)
 sto1<-rbind(sto1,lam_agg[i,k])
 } 
 sto1<-sto1[-1]
 if(j==0) j<-p
 lam_matrix[j,k]<-mean(sto1[30000:50000],na.rm=T)
 }}
  

