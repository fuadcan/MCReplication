/*************************************************************************/
/* THIS VERSION: 13/03/06                                               **/
/* COMPUTE CIPS(P) AND CIPS*(P) STATISTICS IN M.H.PESARAN               **/
/* " A SIMPLE PANEL UNIT ROOT TEST UNDER CROSS SECTION DEPENDENCE"      **/
/*  FEBRUARY 2006                                                       **/
/* NOTE:                                                                **/
/* Please report any problems you might have, to: ty228@cam.ac.uk       **/
/*                                                                      **/
/* {CIPSmat,cadf_p_mat,cadfs_p_mat,res_p_mat}=cips(var_mat,maxp,case)   **/
/* <<INPUTS>>                                                           **/
/* var_mat (T x N) matrix                                               **/
/* N            : number of cross section dimension                     **/
/* T		: number of observations for i's                            **/
/* maxp (more than or equal to zero) : maximum number of (ADF) lag order**/
/*                                    eg. maxp=0: CADF==CDF             **/
/* case=1: no intercept nor trend                                       **/
/* case=2: with intercept                                               **/
/* case=3: with intercept and trend                                     **/
/* <<OUTPUTS>>                                                          **/
/* CIPSmat: (maxp+1) x 2 matrix, first column reports                   **/
/*          CIPS(p) and second column reports CIPS*(p)                  **/
/*          for p=0,1,...,maxp in ascending order                       **/
/*     -CIPS(p) statistic is defined as simple average of               **/
/*       CADF(p)_i statistics                                           **/
/*     -CIPS*(p) statistic is defined as simple average of              **/
/*       CADF*(p)_i statistics, which is truncated version of CADF_i.   **/   
/*                                                                      **/   
/* cadf_p_mat: a N x (maxp+1) matrix, reports all CADF(p)_i statistics  **/   
/* cadfs_p_mat: a N x (maxp+1) matrix, reports all CADF*(p)_i statistics**/
/* res_p_mat: a (maxp+1) set of						                    **/ 
/*		T-(maxp+1) x N matrix of residuals of CADF(p)_i regressions	    **/
/*              concatenated vertically, from the top p=0,p=1,,,p=maxp  **/    
/*                                                                      **/   
/* NB: Appropriate critical values for CIPS statistic are found         **/
/*  are in Table 2 in Pesaran (2006)                                    **/    
/*                                                                      **/
/* NOTEs: Default setting is to print the results. It can be            **/
/*       suppressed by choosing "outpt=0" below (7th-line of the code   **/
/* Takashi Yamagata, 13 March2006                                       **/
/*                                                                      **/
/*************************************************************************/
proc(4)=cips(var_mat,maxp,case);
local outpt,n,t,tlag,var_mat_1,Dvar_mat_temp,Dvar_mat,Dvar_CSM,var_1_CSM;
local tt,y,z,hi,x,h,k1,k2,bvec,sevec,t_vec_cce,t_bar,i,c_p,temp;
local trncl,truncu,trnc1,cadf_p,cadf_ps,cips_p_mat,cips_p,cips_ps;
local cadf_p_mat,cadfs_p_mat,idvec,lgorder,res_p_mat,res_p;


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

Dvar_CSM=meanc(Dvar_mat');
var_1_CSM=meanc(var_mat_1');

tt=t-tlag;

/*************/
y=vec(Dvar_mat);
Hi=Dvar_CSM~var_1_CSM;
x=vec(var_mat_1);

/****/
if case==1;k1=-6.12;k2=4.16;
elseif case==2;Z=ones(n*tt,1);k1=-6.19;k2=2.61;
elseif case==3;Z=ones(n*tt,1)~(ONES(N,1).*.seqa(1,1,tt));k1=-6.42;k2=1.70;
else;"please choose the 'case' only from 1,2, or 3";end;
endif;

/**************** ESTIMATION ****************/
c_p=0;
cips_p_mat={};
cadf_p_mat={};
cadfs_p_mat={};
res_p_mat={};
do while c_p<=maxp;

if c_p>0;
    temp=(Dvar_mat_temp[tlag+1-c_p:t-c_p,.]);
    Hi=Hi~meanc(temp'); 
    x=x~vec(temp);
endif;

/*************/
if case==1;H=(ONES(N,1).*.Hi);
else;H=Z~(ONES(N,1).*.Hi);
endif;

/*** MGCCE ****/
{bvec,sevec,t_vec_cce,res_p}=mgsimple(y,x~h,n,tt);

trncl=k1*(t_vec_cce[.,1].<k1);truncu=k2*(t_vec_cce[.,1].>k2);
trnc1=(t_vec_cce[.,1].>k1).*(t_vec_cce[.,1].<k2);

cadf_p = t_vec_cce[.,1];
cadf_ps = (t_vec_cce[.,1].*trnc1+trncl+truncu);

cips_p=meanc(cadf_p);
cips_ps=meanc(cadf_ps);

if outpt==1;
format /rd 2,0;
"@@@@@@@@@@@@@@@@@@@ CADF(";;C_P;;") @@@@@@@@@@@@@@@@@@@@";
format /rd 2,0;
"CIPS TEST: CASE";;case;;
if case==1;", NO INTERCEPT";
elseif case==2;", WITH INTERCEPT";
elseif case==3;", WITH INTERCEPT AND TREND";
endif;
format /rd 5,0;
"N=";;n;;", T=";;tt;;", (";;t;;"data points used)";
format /rd 8,3;
"   CIPS            :";;cips_p;
"   CIPS*(truncated):";;cips_ps;
"* Truncation is done for CADF_i in such a way that" ;
format /rd 1,2;
"when CADF_i<k1, CADF_i=k1 and when CADF_i>k2, CADF_i=k2,";
"where ";;"k1=";;k1;;" and k2=";;k2;;".";
"Appropriate critical values are in Table 2, Pesaran (2006).";
"-------------------------------------------------------";
"";
endif;
cips_p_mat=cips_p_mat|(cips_p~cips_ps);
cadf_p_mat=cadf_p_mat~cadf_p;
cadfs_p_mat=cadfs_p_mat~cadf_ps;

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
"CADF_i(p) Statistics: Case";;case;;
if case==1;", NO INTERCEPT";
elseif case==2;", WITH INTERCEPT";
elseif case==3;", WITH INTERCEPT AND TREND";
endif;
format /rd 8,0;
$ (lgorder|(idvec~ftocv(cadf_p_mat,1,3)));
"";
"Appropriate critical values are in Table 1, Pesaran (2006).";
"*****************************************************";
format /rd 2,0;
"CADF*_i(p) Statistics (truncated): Case";;case;;
if case==1;", NO INTERCEPT";
elseif case==2;", WITH INTERCEPT";
elseif case==3;", WITH INTERCEPT and TREND";
endif;
format /rd 8,0;
$ (lgorder|(idvec~ftocv(cadfs_p_mat,1,3)));
"* Truncation is done for CADF_i in such a way that" ;
format /rd 1,2;
"when CADF_i<k1, CADF_i=k1 and when CADF_i>k2, CADF_i=k2,";
"where ";;"k1=";;k1;;" and k2=";;k2;;".";
"Appropriate critical values are in Table 1, Pesaran (2006).";
endif;

format /rd 16,8;
retp(cips_p_mat,cadf_p_mat,cadfs_p_mat,res_p_mat);
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