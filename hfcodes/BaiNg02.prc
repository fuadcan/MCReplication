/********************************************************/
/* Gauss Procedure to estimate the number of factors    */
/* Bai and Ng (2002), DETERMINING THE NUMBER OF         */
/* FACTORS IN APPROXIMATE FACTOR MODELS, ECONOMETRICA   */
/* VOL70, 191-221.                                      */
/*                                                      */
/* {nbfactor}=bai02(xmat,kmax);                         */
/*                                                      */
/* <<INPUT>>                                            */
/* xmat: a (NxT) data matrix, where                     */
/*          N: cross section dimension                  */
/*          T: time series dimension                    */
/*                                                      */
/* kmax: maximum number of factors to set               */
/*                                                      */
/* <<OUTPUT>>                                           */
/* nbfactor: a (6 x 1) vector of estimated number of    */
/*          factors from 0,1,...,kmax                   */
/*          Each element is based on                    */
/*          Information Criteria defined in             */
/*          Bai and Ng (2002; p.201):                   */
/*          nbfactor=(PC1,PC2,PC3,IC1,IC2,IC3)'         */
/* penmat: a (maxk+1 x 6) matrix of the information     */
/*          criteria. Rows correspond to k=0,1,..., maxk*/
/*          and columns corresponds to the information  */
/*          criterion (PC1,PC2,PC3,IC1,IC2,IC3).        */
/*                                                      */
/*                                                      */
/* Takashi Yamagata March 2006                          */
/********************************************************/

proc(2)=bai02(xmat,kmax);
local Fmatmax,Lmax,vmax,f,l;
local x,n,t,k,penmat,res,v,c_nt,p1,p2,p3;
local pc1,pc2,pc3,ic1,ic2,ic3,khat;

{Fmatmax,Lmax,vmax}=factor(xmat,kmax);
x=(xmat-meanc(xmat)')./(stdc(xmat)');
n=cols(xmat);t=rows(xmat);

k=0;
penmat={};

do while k<=kmax;


if k==0;
res=x;
else;
L=x/Fmatmax[.,1:k];
res=x-Fmatmax[.,1:k]*L;
endif;
v=meanc(diag(res'res/T));/* V(k,Fhat)*/
c_nt=minc(sqrt(N)|sqrt(T));
p1=ln(n*t/(n+t))*k*(n+t)/(n*t);
p2=ln(c_nt^2)*k*(n+t)/(n*t);
p3=ln(c_nt^2)*k/(c_nt^2);

/* IC's in Bai Ng 2002 */
pc1=v+vmax*p1;
pc2=v+vmax*p2;
pc3=v+vmax*p3;
ic1=ln(v)+p1;
ic2=ln(v)+p2;
ic3=ln(v)+p3;

penmat=penmat|(pc1~pc2~pc3~ic1~ic2~ic3);


k=k+1;
endo;

khat=minindc(penmat)-1;
retp(khat,penmat);
endp;

/* Obtain Factor and 
factor loadings given k */
proc(3)=factor(xmat,k);
local x,n,t,ei,eivec,Fmat,L,res,v;

x=(xmat-meanc(xmat)')./(stdc(xmat)');/*standardization*/
n=cols(xmat);t=rows(xmat);/*N and T*/
if n>t;
{ei,eivec}=eighv(x*x');

eivec=rev(eivec')';

Fmat=sqrt(T)*eivec[.,1:k];/*Ftilde of Bai-Ng (2002)*/
L=(Fmat'x)'/T;/*Lambda-tilde of Bai-Ng (2002)*/

else;
{ei,eivec}=eighv(x'x);

eivec=rev(eivec')';

L=sqrt(N)*eivec[.,1:k];/*Lambda-bar of Bai-Ng (2002)*/
Fmat=x*L/N;/*Fbar of Bai-Ng (2002)*/

endif;

res=x-Fmat*L';

v=meanc(diag(res'res/T));

retp(Fmat,L,v);
endp;