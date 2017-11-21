/*************************************************************************/
/* THIS VERSION: 13/03/06                                               **/
/* COMPUTE tbar(P) (AND tbar*(P)) STATISTICS IN IM, PESARAN, SHIN (2003)**/
/* "Testing for unit roots in heterogeneous panels", Journal of         **/
/*  Econometrics 115, 53-74.                                            **/
/* For the truncated version, see Pesaran (2006),                       **/
/* " A SIMPLE PANEL UNIT ROOT TEST UNDER CROSS SECTION DEPENDENCE"      **/
/*  FEBRUARY 2006                                                       **/
/* NOTE:                                                                **/
/* Please report any problems you might have, to: ty228@econ.cam.ac.uk  **/
/*                                                                      **/
/* {tbarmat,adf_p_mat,adfs_p_mat,res_p_mat}=ips(var_mat,maxp,case);     **/
/*                                                                      **/
/* <<INPUTS>>                                                           **/
/* var_mat (T x N) matrix                                               **/
/* N            : number of cross section dimension                     **/
/* T		: number of observations for i's                            **/
/* maxp (more than or equal to zero) : maximum number of (ADF) lag order**/
/*                                    eg. maxp=0: ADF==DF               **/
/* case=1: no intercept nor trend                                       **/
/* case=2: with intercept                                               **/
/* case=3: with intercept and trend                                     **/
/* <<OUTPUTS>>                                                          **/
/* tbarmat: (maxp+1) x 2 matrix, first column reports                   **/
/*          tbar(p) and second column reports tbar*(p)                  **/
/*          for p=0,1,...,maxp in ascending order                       **/
/*     -tbar(p) statistic is defined as simple average of               **/
/*       adf(p)_i statistics                                            **/
/*     -tbar*(p) statistic is defined as simple average of              **/
/*       adf*(p)_i statistics, which is truncated version of adf(p)_i.  **/
/*    Notes: For a discussion about the truncation, See Pesaran (2006), **/
/*      " A Simple Panel Unit Root Test Under Cross Section Dependence" **/
/*                                                                      **/   
/* adf_p_mat: a N x (maxp+1) matrix, reports all ADF(p)_i statistics    **/   
/* adfs_p_mat: a N x (maxp+1) matrix, reports all ADF*(p)_i statistics  **/   
/* res_p_mat: a (maxp+1) set of						                    **/ 
/*		T-(maxp+1) x N matrix of residuals of  ADF(p)_i regressions	    **/
/*              concatenated vertically, from the top p=0,p=1,,,p=maxp  **/    
/*                                                                      **/   
/* IPS is sqrt(N)*(TBAR - E(TBAR))/sqrt(Var(TBAR)).                     **/
/* Appropriate E(TBAR) and Var(TBAR) are found in                       **/
/* Table 3 of Im,Pesaran,Shin (2003)                                    **/
/* NOTEs: Default setting is to print the results. It can be            **/
/*       suppressed by choosing "outpt=0" below (7th-line of the code   **/
/* Takashi Yamagata, 13 March 2006                                        **/
/*                                                                      **/
/*************************************************************************/
proc(4)=ips(var_mat,maxp,case);
local outpt,n,t,tlag,var_mat_1,Dvar_mat_temp,Dvar_mat;
local tt,y,z,x,xx,k1,k2,bvec,sevec,t_vec,t_bar,i,c_p,temp;
local trncl,truncu,trnc1,adf_p,adf_ps,tbar_p_mat,tbar_p,tbar_ps;
local adf_p_mat,adfs_p_mat,idvec,lgorder,res_p_mat,res_p;


/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
outpt=1;/* 1: reports ourput, 0:supress the output */
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
outwidth 256;

t=rows(var_mat);
n=cols(var_mat);
tlag=maxp+1;

var_mat_1=lag(var_mat);
Dvar_mat_temp=var_mat-var_mat_1;
    var_mat_1=var_mat_1[tlag+1:t,.];
    Dvar_mat=Dvar_mat_temp[tlag+1:t,.];

tt=t-tlag;

/*************/
y=vec(Dvar_mat);
x=vec(var_mat_1);

/****/
if case==1;k1=-6.12;k2=4.16;
elseif case==2;Z=ones(n*tt,1);k1=-6.19;k2=2.61;
elseif case==3;Z=ones(n*tt,1)~(ONES(N,1).*.seqa(1,1,tt));k1=-6.42;k2=1.70;
else;"please choose the 'case' only from 1,2, or 3";end;
endif;

/**************** ESTIMATION ****************/
c_p=0;
tbar_p_mat={};
adf_p_mat={};
adfs_p_mat={};
res_p_mat={};
do while c_p<=maxp;

if c_p>0;
    temp=(Dvar_mat_temp[tlag+1-c_p:t-c_p,.]);
    x=x~vec(temp);
endif;

/*************/
if case==1;xx=x;
else;xx=x~Z;
endif;

/*** MGCCE ****/
{bvec,sevec,t_vec,res_p}=mgsimple(y,xx,n,tt);

trncl=k1*(t_vec[.,1].<k1);truncu=k2*(t_vec[.,1].>k2);
trnc1=(t_vec[.,1].>k1).*(t_vec[.,1].<k2);

adf_p = t_vec[.,1];
adf_ps = (t_vec[.,1].*trnc1+trncl+truncu);

tbar_p=meanc(adf_p);
tbar_ps=meanc(adf_ps);

if outpt==1;
format /rd 2,0;
"@@@@@@@@@@@@@@@@@@@ ADF(";;C_P;;") @@@@@@@@@@@@@@@@@@@@";
format /rd 2,0;
"tbar STATISTIC FOR IPS(P): CASE";;case;;
if case==1;", NO INTERCEPT";
elseif case==2;", WITH INTERCEPT";
elseif case==3;", WITH INTERCEPT AND TREND";
endif;
format /rd 5,0;
"N=";;n;;", T=";;tt;;", (";;t;;"data points used)";
format /rd 8,3;
"   TBAR(P)            :";;tbar_p;
"   TBAR(P)*(truncated):";;tbar_ps;
format /rd 1,2;
"*Truncation is done for ADF_i in such a way that";
"when ADF_i<k1, ADF_i=k1 and when ADF_i>k2, ADF_i=k2,";
"where ";;"k1=";;k1;;" and k2=";;k2;;".";
"IPS is sqrt(N)*(TBAR - E(TBAR))/sqrt(Var(TBAR)).";
"Appropriate E(TBAR) and Var(TBAR) are found in Table 3 of Im,Pesaran,Shin (2003).";
"-------------------------------------------------------";
"";
endif;
tbar_p_mat=tbar_p_mat|(tbar_p~tbar_ps);
adf_p_mat=adf_p_mat~adf_p;
adfs_p_mat=adfs_p_mat~adf_ps;
res_p_mat=res_p_mat|res_p;

    c_p=c_p+1;
endo;

if outpt==1;
format /rd 8,0;
idvec=ftocv(seqa(1,1,n),1,0);
lgorder="id/p"~ftocv(seqa(0,1,maxp+1)',1,0);

format /rd 2,0;
"";
"*****************************************************";
"ADF_i(p) Statistics: Case";;case;;
if case==1;", NO INTERCEPT";
elseif case==2;", WITH INTERCEPT";
elseif case==3;", WITH INTERCEPT AND TREND";
endif;
format /rd 8,0;
$ (lgorder|(idvec~ftocv(adf_p_mat,1,3)));
"";
"*****************************************************";
format /rd 2,0;
"ADF*_i(p) Statistics (truncated): Case";;case;;
if case==1;", NO INTERCEPT";
elseif case==2;", WITH INTERCEPT";
elseif case==3;", WITH INTERCEPT and TREND";
endif;
format /rd 8,0;
$ (lgorder|(idvec~ftocv(adfs_p_mat,1,3)));
format /rd 1,2;
"*Truncation is done for ADF_i in such a way that";
"when ADF_i<k1, ADF_i=k1 and when ADF_i>k2, ADF_i=k2,";
"where ";;"k1=";;k1;;" and k2=";;k2;;".";
endif;

format /rd 16,8;
retp(tbar_p_mat,adf_p_mat,adfs_p_mat,res_p_mat);
endp;

/*** OLS regression for each cross section unit ****/
proc(4)=mgsimple(y,x,n,t);
local bvec,sevec,tvec,i,y_i,x_i,beta_i,e_i,zig2,se_i,t_i,res_mat;
bvec=zeros(n,cols(x));
sevec=zeros(n,cols(x));
tvec=zeros(n,cols(x));
res_mat={};
i=1;
do while i<=n;
    y_i=y[1+(i-1)*t:(i)*t,.];
    x_i=x[1+(i-1)*t:(i)*t,.];
    beta_i=y_i/x_i;
    e_i=y_i-x_i*beta_i;
    zig2=e_i'e_i/(T-cols(x));
    se_i=sqrt(diag(zig2*invpd(x_i'x_i)));
    t_i=beta_i./se_i;

    bvec[i,.]=beta_i';
    sevec[i,.]=se_i';
    tvec[i,.]=t_i';
    res_mat=res_mat~e_i;

i=i+1;
endo;

retp(bvec,sevec,tvec,res_mat);
endp;