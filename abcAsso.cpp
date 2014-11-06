 /*
  thorfinn thorfinn@binf.ku.dk dec17 2012
  
  has modified so no sites[].chromo etc are used
    
  anders albrecht@binf.ku.dk made this.

  part of angsd

  DRAGON positions are not offset correctly
*/
#include <algorithm>
#include <cmath>
#include <functional> 
#include <numeric>
#include <time.h>
#include <vector>
#include <zlib.h>
#include "kstring.h"
#include "shared.h"
#include "analysisFunction.h"

#include "abcFreq.h"
#include "abcAsso.h"



void abcAsso::printArg(FILE *argFile){
  fprintf(argFile,"-------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doAsso\t%d\n",doAsso);
  fprintf(argFile,"\t1: Frequency Test (Known Major and Minor)\n");
  fprintf(argFile,"\t2: Score Test\n");
  fprintf(argFile,"\t3: Frequency Test (Unknown Minor)\t\n");
  fprintf(argFile,"\t4: Robust Variance Score Test\t\n");  
  fprintf(argFile,"\t5: Robust Variance Score Burden Test (for rare variants)\t\n");    
  fprintf(argFile,"  Frequency Test Options:\n");
  fprintf(argFile,"\t-yBin\t\t%s\t(File containing disease status)\t\n\n",yfile);
  fprintf(argFile,"  Score Test Options:\n");
  fprintf(argFile,"\t-yBin\t\t%s\t(File containing disease status)\n",yfile);
  fprintf(argFile,"\t-yQuant\t\t%s\t(File containing phenotypes)\n",yfile);
  fprintf(argFile,"\t-minHigh\t%d\t(Require atleast minHigh number of high credible genotypes)\n",minHigh);
  fprintf(argFile,"\t-minCount\t%d\t(Require this number of minor alleles, estimated from MAF)\n",minCount);
  fprintf(argFile,"\t-numBootstraps\t%d\t(The number of bootstrap samples to generate when performing rare RVS burden testing)\n",numBootstraps);
  fprintf(argFile,"\t-cov\t\t%s\t(File containing additional covariates)\n",covfile);
  fprintf(argFile,"\t-model\t%d\n",model);
  fprintf(argFile,"\t1: Additive/Log-Additive (Default)\n");
  fprintf(argFile,"\t2: Dominant\n");
  fprintf(argFile,"\t3: Recessive\n\n");
  fprintf(argFile,"Examples:\n\tPerform Frequency Test\n\t  \'./angsd -yBin pheno.ybin -doAsso 1 -GL 1 -out out -doMajorMinor 1 -minLRT 24 -doMaf 2 -doSNP 1 -bam bam.filelist'\n");
  fprintf(argFile,"\tPerform Score Test\n\t  \'./angsd -yBin pheno.ybin -doAsso 2 -GL 1 -doPost 1 -out out -doMajorMinor 1 -minLRT 24 -doMaf 2 -doSNP 1 -bam bam.filelist'\n");
  fprintf(argFile,"\n");
}


void abcAsso::getOptions(argStruct *arguments){


  doAsso=angsd::getArg("-doAsso",doAsso,arguments);

  doMaf=angsd::getArg("-doMaf",doMaf,arguments);

  adjust=angsd::getArg("-adjust",adjust,arguments);
  model=angsd::getArg("-model",model,arguments);
  minCov=angsd::getArg("-minCov",minCov,arguments);
  dynCov=angsd::getArg("-dynCov",dynCov,arguments);
  minHigh=angsd::getArg("-minHigh",minHigh,arguments);
  doPrint=angsd::getArg("-doPrint",doPrint,arguments);
  minCount=angsd::getArg("-minCount",minCount,arguments);
  numBootstraps=angsd::getArg("-numBootstraps",numBootstraps,arguments);  
  sitePerm=angsd::getArg("-sitePerm",sitePerm,arguments);
  GL=angsd::getArg("-GL",GL,arguments);
  covfile=angsd::getArg("-cov",covfile,arguments);
  doPost=angsd::getArg("-doPost",doPost,arguments);
  yfile=angsd::getArg("-yBin",yfile,arguments);
  if(yfile!=NULL)
    isBinary=1;
  yfile=angsd::getArg("-yQuant",yfile,arguments);

  if(doPrint)
    fprintf(stderr,"finished [%s]\t[%s]\n",__FILE__,__FUNCTION__);


  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  if(doAsso && doMaf==0){
    fprintf(stderr,"Error: you must estimate the maf (-doMaf) in order to perform association \n");
    exit(0);

  }

  if(doAsso && yfile==NULL){
    fprintf(stderr,"Error: you must provide a phenotype file (-yBin or -yQuant) to perform association \n");
    exit(0);
   }
 if(doAsso && arguments->inputtype==INPUT_BEAGLE&&doAsso==1){
    fprintf(stderr,"Error: Only doAsso=2 can be performed on posterior input\n");
    exit(0);
  }
  if(doAsso && arguments->inputtype!=INPUT_BEAGLE&&(doAsso==2)&&doPost==0){
    fprintf(stderr,"Error: For doAsso=2 you must estimate the posterior probabilites for the genotypes (doPost!=0) \n");
    exit(0);
  }  

}
 
abcAsso::abcAsso(const char *outfiles,argStruct *arguments,int inputtype){
  bufstr.s=NULL;bufstr.l=bufstr.m=0;
  //default
  model=1;
  doPrint=0;
  doAsso=0;
  GL=0;
  doPost=0;
  isBinary=0;
  sitePerm=0;//not for users
  covfile=NULL;
  yfile=NULL;
  minHigh=10;
  minCount=10;
  numBootstraps=1000;
  dynCov=0;//not for users
  minCov=5;//not for users
  adjust=1;//not for users
  doMaf=0;
  //from command line


  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doAsso")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }
  

  getOptions(arguments);
  printArg(arguments->argumentFile);

  if(doAsso==0){
    shouldRun[index] =0;
    return;
  }

  //read phenotype
  if(isBinary)
    ymat = angsd::getMatrix(yfile,1,100000);
  else
    ymat = angsd::getMatrix(yfile,0,100000);

  //read covariates 
  if(covfile!=NULL)
    covmat = angsd::getMatrix(covfile,0,100000);
  else{
    covmat.x=0;
    covmat.y=0;
    covmat.matrix=NULL;
  }
  if(covfile!=NULL&&(covmat.x!=ymat.x)){
    fprintf(stderr,"The number of covariates (%d) does not match the number of phenotypes (%d)\n",covmat.x,ymat.x);
    exit(0);
  }


 
  if(!isBinary&&doAsso==1){
    fprintf(stderr,"Error: Only doAsso=2 can be performed on quantitative traits\n");
    exit(0);
  }

  //make output files
  MultiOutfile = new gzFile[ymat.y];
  const char* postfix;
  postfix=".lrt";

  if(covfile!=NULL&&minCov>0){
    int keepList[ymat.x];
    for(int i=0 ; i < ymat.x;i++) {
      keepList[i]=1;
      for(int yi=0;yi<ymat.y;yi++) {
	if(ymat.matrix[i][yi]==-999)
	  keepList[i]=0;
      }
      for(int ci=0;ci<covmat.y;ci++) {
	if(covmat.matrix[i][ci]==-999)
	  keepList[i]=0;
      }
    }
    int nCov=0;
    int count[covmat.y];
    for(int ci=0;ci<covmat.y;ci++) {
      count[ci]=0;
      for(int i=0 ; i < ymat.x;i++) {
	if(keepList[i]==0)
	  continue;
	if(covmat.matrix[i][ci]!=0){
	  count[ci]++;
	}
      }
  
      if(count[ci]<minCov){
	fprintf(stderr,"Error: Cov #%d only has %d non zero entries\n",ci,count[ci]);
      }
      else
	nCov++;

    }
    if(!dynCov&&covmat.y!=nCov){
      fprintf(stderr,"Error: Creating new covariant matrix with %d columns\n",nCov);
      exit(0);

    }
    else if(covmat.y!=nCov){
      //      angsd::printMatrix(covmat,stderr);
      fprintf(stderr,"Error: Creating new covariant matrix with %d columns\n",nCov);
      angsd::Matrix<double> newmat;
      newmat.x=covmat.x;
      newmat.y=nCov;
      newmat.matrix=new double*[covmat.x];
      for(int xi=0;xi<covmat.x;xi++){
	newmat.matrix[xi] = new double[nCov];
	int tempCount=0;
	for(int ci=0;ci<covmat.y;ci++){
	  if(count[ci]>minCov){
	    newmat.matrix[xi][tempCount]=covmat.matrix[xi][ci];
	    tempCount++;
	  }
	}
      }
      angsd::deleteMatrix(covmat);
      covmat=newmat;
      
    }
  }

  //open outfiles
  for(int i=0;i<ymat.y;i++){
    char ary[5000];
    snprintf(ary,5000,"%s%d.gz",postfix,i);
    MultiOutfile[i] = Z_NULL;
    MultiOutfile[i] = aio::openFileGz(outfiles,ary,GZOPT);
  }

  //print header
  for(int yi=0;yi<ymat.y;yi++){
    if(doAsso==2)
      gzprintf(MultiOutfile[yi],"Chromosome\tPosition\tMajor\tMinor\tFrequency\tN\tLRT\thigh_WT/HE/HO\n");
    else if(doAsso==4)
      gzprintf(MultiOutfile[yi],"Chromosome\tPosition\tMajor\tMinor\tFrequency\tN\tLRT_filtered\tLRT_rvs\tLRT_std\thigh_WT/HE/HO\n");
    else if(doAsso==5)
      gzprintf(MultiOutfile[yi],"Chromosome\tPosition\tMajor\tMinor\tFrequency\tLRT_rvs\tCAST\tnPerm\tCONT_high_WT/HE/HO\tCASE_high_WT/HE/HO\n");
    else
      gzprintf(MultiOutfile[yi],"Chromosome\tPosition\tMajor\tMinor\tFrequency\tLRT\n");
  }
}


abcAsso::~abcAsso(){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);


  if(doAsso==0)
    return;
  for(int i=0;i<ymat.y;i++)
    if(MultiOutfile[i]) gzclose(MultiOutfile[i]);
  delete [] MultiOutfile;

  if(covfile!=NULL)
    angsd::deleteMatrix(covmat);
  angsd::deleteMatrix(ymat);

}


void abcAsso::clean(funkyPars *pars){
  if(doAsso==0)
    return;

  assoStruct *assoc =(assoStruct*) pars->extras[index];

  if(assoc->stat!=NULL){
    for(int yi=0;yi<ymat.y;yi++){
      delete[]  assoc->stat[yi];
      delete[]  assoc->highWt[yi];
      delete[]  assoc->highHe[yi];
      delete[]  assoc->highHo[yi];
    }
  }

  delete[]  assoc->stat; 
  delete[]  assoc->std_LRT;
  delete[]  assoc->rvs_LRT;
  delete[]  assoc->burden;
  delete[]  assoc->numPerm;

  if(assoc->keepInd!=NULL)
    for( int yi =0;yi<ymat.y;yi++)
      delete[] assoc->keepInd[yi];
  delete[] assoc->keepInd;
  
  delete assoc;
}



void abcAsso::print(funkyPars *pars){
  if(doAsso==0)
    return;

  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  printDoAsso(pars);
}

assoStruct *allocAssoStruct(){
  assoStruct *assoc = new assoStruct;
  
  assoc->stat=NULL;
  assoc->keepInd=NULL;
  assoc->highWt = NULL;
  assoc->highHe = NULL;
  assoc->highHo = NULL;
  assoc->std_LRT = NULL;
  assoc->rvs_LRT = NULL;
  assoc->burden = NULL;
  assoc->numPerm = NULL;

  return assoc;
}


void abcAsso::run(funkyPars *pars){


  if(doAsso==0)
    return;
  
  assoStruct *assoc = allocAssoStruct();
  pars->extras[index] = assoc;

  if(doAsso==1||doAsso==3){
    frequencyAsso(pars,assoc);
  }
  else if(doAsso==2 || doAsso==4 || doAsso==5){
 
    assoc->highWt=new int**[ymat.y];
    assoc->highHe=new int**[ymat.y];
    assoc->highHo=new int**[ymat.y];

    for(int yi=0;yi<ymat.y;yi++){
      assoc->highWt[yi]=new int*[2];
      assoc->highHe[yi]=new int*[2];
      assoc->highHo[yi]=new int*[2];
      for(int i=0;i<=1;i++){
        assoc->highWt[yi][i]=new int[pars->numSites];
        assoc->highHe[yi][i]=new int[pars->numSites];
        assoc->highHo[yi][i]=new int[pars->numSites];     
      }
    }

    assoc->std_LRT=new double[pars->numSites];
    assoc->rvs_LRT=new double[pars->numSites];
    assoc->burden=new double[ymat.y];
    assoc->numPerm=new int[ymat.y];

    scoreAsso(pars,assoc);
  }

}

void abcAsso::frequencyAsso(funkyPars  *pars,assoStruct *assoc){

  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  if(pars->nInd!=ymat.x){
    fprintf(stderr,"The number of sequenced individuals (%d) does not match the number of phenotypes (%d)\n",pars->nInd,ymat.x);

    fflush(stderr);
    exit(0);
  }

  if(ymat.y!=1){
    fprintf(stderr,"Only one phenotype allowed for doAsso 1 or 3 \n");
    fflush(stderr);
    exit(0);
  }

  double **stat = new double*[ymat.y];
  for(int yi=0;yi<ymat.y;yi++)
    stat[yi] = new double[pars->numSites];
 
 
  
  int y[pars->nInd];
  for(int i=0;i<pars->nInd;i++)
    y[i]=ymat.matrix[i][0];

  double **like0;//genotype likelihood for controls
  double **like1;//genotype likelihood for cases
  double **likeAll;//genotype likelihood for cases and controls
  int Ncases=0; //number of cases
  int Ncontrols=0; //number of cases
  int Nall=0; //number of cases and controls
  int cases[pars->nInd];
  int controls[pars->nInd];
  int all[pars->nInd];

  for(int i=0;i<pars->nInd;i++){
    cases[i]=0;
    controls[i]=0;
    all[i]=0;
    if((int)y[i]==1){
      Ncases++;
      cases[i]=1;
      Nall++;
      all[i]=1;
    }
    if((int)y[i]==0){
      Ncontrols++;
      controls[i]=1;
      Nall++;
      all[i]=1;
    }
    if(doPrint)
      fprintf(stderr,"all, case, control: %d %d %d\n",all[i],cases[i],controls[i]);
  }
  if(doPrint)
    fprintf(stderr,"count complete [%s]\t[%s]\n",__FILE__,__FUNCTION__);


  if(doAsso==1){
    like0=angsd::get3likes(pars,controls);
    like1=angsd::get3likes(pars,cases);
    likeAll=angsd::get3likes(pars,all);
  }
  if(doAsso==3){//use all 10 genotype likes
    like0=angsd::getlikes(pars,controls);
    like1=angsd::getlikes(pars,cases);
    likeAll=angsd::getlikes(pars,all);
  }
 if(doPrint)
    fprintf(stderr,"like complete [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  for(int s=0;s<pars->numSites;s++){//loop overs sites
    stat[0][s]=-999;
    if(pars->keepSites[s]==0)
      continue;
    if(doAsso==1){

      if(doPrint)
	fprintf(stderr,"do freq [%s]\t[%s]\n",__FILE__,__FUNCTION__);


      double score0=abcFreq::likeFixedMinor(abcFreq::emFrequency_ext(like0[s],Ncontrols,NULL,Ncontrols),like0[s],Ncontrols);
      //likelihood for the cases
      double score1=abcFreq::likeFixedMinor(abcFreq::emFrequency_ext(like1[s],Ncases,NULL,Ncases),like1[s],Ncases);
      //likelihood for all individuals
      double scoreNull=abcFreq::likeFixedMinor(abcFreq::emFrequency_ext(likeAll[s],Nall,NULL,Nall),likeAll[s],Nall);
      //likelhood ratio statistics \sim chi^2
      double LRT=-2*(score0+score1-scoreNull);
      stat[0][s]=LRT;
    }
    if(doAsso==3){
      double score0=abcFreq::likeNoFixedMinor(abcFreq::emFrequencyNoFixed_ext(like0[s],Ncontrols,NULL,Ncontrols,pars->major[s],0),like0[s],Ncontrols,pars->major[s]);//last is posiiton in inner function, only used when printing when there are bugs
      //likelihood for the cases
      double score1=abcFreq::likeNoFixedMinor(abcFreq::emFrequencyNoFixed_ext(like1[s],Ncases,NULL,Ncases,pars->major[s],0),like1[s],Ncases,pars->major[s]);
      //likelihood for all individuals
      double scoreNull=abcFreq::likeNoFixedMinor(abcFreq::emFrequencyNoFixed_ext(likeAll[s],Nall,NULL,Nall,pars->major[s],0),likeAll[s],Nall,pars->major[s]);
      //likelhood ratio statistics \sim chi^2
      double LRT=-2*(score0+score1-scoreNull);
      stat[0][s]=LRT;
    }
   

  }
  assoc->stat=stat;
  for(int s=0;s<pars->numSites;s++){
    delete[] like0[s];
    delete[] like1[s];
    delete[] likeAll[s];
  }
   delete[] like0;
   delete[] like1;
   delete[] likeAll;


 if(doPrint)
    fprintf(stderr,"finish [%s]\t[%s]\n",__FILE__,__FUNCTION__);

}

void abcAsso::scoreAsso(funkyPars  *pars,assoStruct *assoc){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  if(pars->nInd!=ymat.x){
    fprintf(stderr,"The number of sequenced individuals (%d) does not match the number of phenotypes (%d)\n",pars->nInd,ymat.x);
    exit(0);
  }


  int **keepInd  = new int*[ymat.y];
  double **stat = new double*[ymat.y];
  for(int yi=0;yi<ymat.y;yi++){
    stat[yi] = new double[pars->numSites];
    keepInd[yi]= new int[pars->numSites];
  }
  
  if(doAsso == 5){
    stat = binomRVScoreEnvRare(pars,assoc);
  }
  else{
    for(int s=0;s<pars->numSites;s++){//loop overs sites
      if(pars->keepSites[s]==0)
        continue;
      
      
      int *keepListAll = new int[pars->nInd];
      for(int i=0 ; i<pars->nInd ;i++){
        keepListAll[i]=1;

      }

      for(int yi=0;yi<ymat.y;yi++) { //loop over phenotypes
        int *keepList = new int[pars->nInd];
        keepInd[yi][s]=0;
        for(int i=0 ; i<pars->nInd ;i++) {
  	keepList[i]=1;
  	if(keepListAll[i]==0||ymat.matrix[i][yi]==-999)
  	  keepList[i]=0;
  	if(covfile!=NULL)
  	  for(int ci=0;ci<covmat.y;ci++) {
  	    if(covmat.matrix[i][ci]==-999)
  	      keepList[i]=0;
  	  }


  	if(keepList[i]==1)
  	  keepInd[yi][s]++;
        }  
        double *y = new double[pars->nInd];
        for(int i=0 ; i<pars->nInd ;i++)
  	y[i]=ymat.matrix[i][yi]; 
   
        freqStruct *freq = (freqStruct *) pars->extras[6];
        stat[yi][s]=doAssociation(pars,pars->post[s],y,keepInd[yi][s],keepList,freq->freq[s],s,assoc);
        
        //cleanup
         delete [] y;
        delete [] keepList;

      } //phenotypes end
   
      delete [] keepListAll;
    } // sites end
  }

  assoc->stat=stat;
  assoc->keepInd=keepInd;
}


void abcAsso::getFit(double *res,double *Y,double *covMatrix,int nInd,int nEnv){

  /*
    linear regression. Fits a linear model. 
    res is the predicted values  (eta) = X%*%coef
    Y is the responce
    covMatrix is the design matrix (nInd x nEnv)
    nInd is the number of individuals
    nEnv is the number of predictors (including the intersept)
  */
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  double Xt_y[nEnv];
  double invXtX_Xt_y[nEnv];
  for(int i=0;i<nEnv;i++)
    Xt_y[i]=0;
  for(int i=0;i<nEnv;i++)
    invXtX_Xt_y[i]=0;


 
  //get t(X)%*%y
  for(int x=0;x<nEnv;x++)//col X
    for(int i=0;i<nInd;i++)
      Xt_y[x]+=covMatrix[x*nInd+i]*Y[i];

  //get inv(t(X)%*%X)
  double XtX[nEnv*nEnv];
  for(int i=0;i<nEnv*nEnv;i++)
    XtX[i]=0;

  for(int x=0;x<nEnv;x++)//col X
    for(int y=0;y<nEnv;y++)//row Xt
      for(int i=0;i<nInd;i++)
	XtX[x*nEnv+y]+=covMatrix[y*nInd+i]*covMatrix[x*nInd+i];

  double workspace[2*nEnv];
  angsd::matinv(XtX, nEnv, nEnv, workspace);


  //get (inv(t(X)%*%X))%*%(t(X)%*%y) //this is the coef!
 for(int x=0;x<nEnv;x++)//col X
   for(int y=0;y<nEnv;y++)//row Xt
      invXtX_Xt_y[x]+=XtX[y*nEnv+x]*Xt_y[y];

  //get X%*%(inv(t(X)%*%X))%*%(t(X)%*%y)
   for(int j=0;j<nInd;j++)//row Xt
     res[j]=0;
 
   for(int j=0;j<nInd;j++)//row Xt
     for(int x=0;x<nEnv;x++)
       res[j]+=covMatrix[x*nInd+j]*invXtX_Xt_y[x];
 
}

void abcAsso::getFitBin(double *res,double *Y,double *covMatrix,int nInd,int nEnv){

  double tol = 1e-6;
  /*
    logistic regression. Fits a logistic regression model. 
    res is the estimated coefficients
    Y is the responce
    covMatrix is the design matrix (nInd x nEnv)
    nInd is the number of individuals
    nEnv is the number of predictors (including the intersept)
    // R code
    getFitBin<-function(y,X){
        b<-rep(0,ncol(X))#estimates
	for(i in 1:20){ 
	    eta  <- as.vector(1/(1+exp(-(X%*%b))))
	    change<-solve(t(X*eta*(1-eta)) %*% X ) %*% t(X) %*% ( y - eta )
	    b<-b+ change
	    if(sum(abs(change))<1e-6)
	        break
	    }

	b
    }
    ///////
    coef<-getFitBin(y,X)
    yTilde <- sigm(X%*%coef)
  */
  
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);


  
  double coef[nEnv];
  double eta[nInd];
  double Xt_y[nEnv];
  double invXtX_Xt_y[nEnv];
  double XtX[nEnv*nEnv];


  for(int x=0;x<nEnv;x++)//col X
    coef[x]=0;

  for(int iter=0;iter<100;iter++){

    //set to zero
    for(int x=0;x<nEnv;x++){
      Xt_y[x]=0;
      invXtX_Xt_y[x]=0;
      XtX[x]=0;
    }
    for(int i=0;i<nInd;i++){
      eta[i]=0;
    }


    //eta <- 1/(1+exp(-(X%*%b)))
    for(int i=0;i<nInd;i++){
      for(int x=0;x<nEnv;x++)//col X
      	eta[i]+=coef[x]*covMatrix[x*nInd+i];
      eta[i] = 1.0/(1+exp(-eta[i])); 
    }

    

 
    //get t(X) %*% ( y - eta )
    for(int x=0;x<nEnv;x++)//col X
      for(int i=0;i<nInd;i++)
	Xt_y[x]+=covMatrix[x*nInd+i]*(Y[i]-eta[i]);

    //get solve(t(X*eta*(1-eta)) %*% X )
    for(int x=0;x<nEnv;x++)//col X
      for(int y=0;y<nEnv;y++)//row Xt
	for(int i=0;i<nInd;i++)
	  XtX[x*nEnv+y]+=covMatrix[y*nInd+i] * eta[i] * covMatrix[x*nInd+i];

    double workspace[2*nEnv];
    //    angsd::matinv(XtX, nEnv, nEnv, workspace);
    angsd::svd_inverse(XtX,nEnv,nEnv);
    //S = svd_inverse(S,flag);     
    //get (inv(t(X)%*%X))%*%(t(X)%*%y)
    for(int x=0;x<nEnv;x++)//col X
      for(int y=0;y<nEnv;y++)//row Xt
	invXtX_Xt_y[x]+=XtX[y*nEnv+x]*Xt_y[y];

    double diff = 0;
    for(int x=0;x<nEnv;x++)
      diff += fabs(invXtX_Xt_y[x]);
  
     for(int x=0;x<nEnv;x++)
       coef[x]+=invXtX_Xt_y[x];

     /*
     fprintf(stdout,"iter=%d\n",iter);
     for(int x=0;x<nEnv;x++)
       fprintf(stdout,"%f\t",invXtX_Xt_y[x]);
     fprintf(stdout,"diff=%f\n",diff);
  for(int x=0;x<nEnv;x++)
    fprintf(stdout,"%f\t",coef[x]);
  fprintf(stdout,"\n");
     */

     if(diff<tol)
       break;
  }

  
  //yTilde <- X%in%coef
  for(int j=0;j<nInd;j++)//row Xt
    res[j]=0;
  
  for(int j=0;j<nInd;j++)//row Xt
    for(int x=0;x<nEnv;x++)
      res[j]+=covMatrix[x*nInd+j]*coef[x];
  for(int j=0;j<nInd;j++)
    res[j]= angsd::sigm(res[j]);

  //   for(int x=0;x<nEnv;x++)
  //   fprintf(stdout,"%f\t",coef[x]);
  //fprintf(stdout,"\n");

  
}



double abcAsso::doAssociation(funkyPars *pars,double *postOrg,double *yOrg,int keepInd,int *keepList,double freq,int s,assoStruct *assoc){
  if(doPrint)
    fprintf(stderr,"Staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  double covMatrix[(covmat.y+1)*keepInd];
  double y[keepInd];
  double post[keepInd*3];
  int count=0;
  for(int i=0;i<pars->nInd;i++){
    if(keepList[i]){
      y[count]=yOrg[i];
      for(int g=0;g<3;g++)
	post[count*3+g]=postOrg[i*3+g];
      count++;
    }

  }

  if(count!=keepInd){
    fprintf(stderr,"[%s] wrong number of non missing\n",__FUNCTION__);
    fflush(stderr);
    exit(0);
  }

  int nEnv;
  if(adjust==1)
    nEnv=(covmat.y+1);
  else
    nEnv=1;



  int num;
  for(int j=0;j<keepInd;j++)
    covMatrix[j]=1;
  if(covmat.matrix!=NULL){
    num=0;
    for(int j=0;j<covmat.x;j++){
      if(keepList[j]==0)
	continue;
      for(int i=1;i<nEnv;i++)
	covMatrix[i*keepInd+num]=covmat.matrix[j][i-1];   
      num++;
    }
  }

  //fprintf(stderr,"number of inds %d\n",keepInd);

  // permutation
  if(sitePerm){
    if((covmat.y+1)==1){
      for(int i=0 ; i<keepInd ;i++){	
	int j = rand() % (keepInd);
	angsd::swapDouble(y[j],y[i]); 
      }
    }
    else{
      int col0=0; 
      for(int i=0 ; i<covmat.x ;i++) {
	if(keepList[i]==0)
	  continue;
	if(covmat.matrix[i][0]<0.5)
	  col0++;
      }
      for(int i=0 ; i<col0 ;i++) {
	int j = rand() % (col0);
	angsd::swapDouble(y[j],y[i]);
      }
      for(int i=0 ; i<keepInd-col0 ;i++) {
	int j = rand() % (keepInd-col0);
	angsd::swapDouble(y[j+col0],y[i+col0]);
      }
      if(col0<500||keepInd-col0<500){
	fprintf(stderr,"colTrouble %d %d\n",col0,keepInd-col0);
      }
    }

  }

  //


 
  double *yfit = new double[keepInd];
  if(nEnv==1){
    double mean=0;
    for(int i=0;i<keepInd;i++)
      mean+=y[i];
    mean=mean/keepInd;
    for(int i=0;i<keepInd;i++)
      yfit[i]=mean;
  }
  else{
    if(isBinary)
      getFitBin(yfit,y,covMatrix,keepInd,nEnv);
    else
      getFit(yfit,y,covMatrix,keepInd,nEnv);

  }
  
  //for(int i=0;i<keepInd;i++)
  //   fprintf(stdout,"%f\t",y[i]);
  // exit(0);


  if(model==2){
    for(int i=0 ; i<keepInd ;i++) {
      post[i*3+1]+=post[i*3+2];
      post[i*3+2]=0;
    }
  }
  if(model==3){
    for(int i=0 ; i<keepInd ;i++) {
      post[i*3+0]+=post[i*3+1];
      post[i*3+1]=post[i*3+2];
      post[i*3+2]=0;
    }
  }


  double stat;
  if(isBinary)
    if(doAsso==2)
      stat = binomScoreEnv(post,keepInd,y,yfit,covMatrix,nEnv,freq,assoc,s);
    else if(doAsso==4)
      stat = binomRVScoreEnv(post,keepInd,y,yfit,covMatrix,nEnv,freq,assoc,s);
  else
    stat = normScoreEnv(post,keepInd,y,yfit,covMatrix,nEnv,freq,assoc,s);

  delete[] yfit;
  return stat;

}




double abcAsso::normScoreEnv(double *post,int numInds, double *y, double *ytilde,double *cov,int nEnv,double freq,assoStruct *assoc,int s){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  int rankProb=0;
  double sum=0;
  double *Ex = angsd::allocArray<double>(numInds,0);
  double *Ex2 = angsd::allocArray<double>(numInds,0);
  double U=0;
  int highWT=0;
  int highHE=0;
  int highHO=0;
  
  double sumEx=0;
  
  for(int i=0;i<numInds;i++)
    sum+=pow(y[i]-ytilde[i],2);
  double var=sum/(numInds-nEnv);
  
  double Vaa[nEnv*nEnv];
  for(int x=0;x<nEnv*nEnv;x++)
    Vaa[x]=0;
  double Vab[nEnv];
  for(int x=0;x<nEnv;x++)
    Vab[x]=0;
  
  
  for(int i=0 ; i<numInds ;i++) {
    Ex[i]=post[i*3+1]+2*post[i*3+2];
    Ex2[i]=post[i*3+1]+4*post[i*3+2];
    U+=Ex[i]*(y[i]-ytilde[i])/var;
    
    //Vaa<-Vaa+1/var*Xe[tal,]%*%t(Xe[tal,])
    for(int Nx=0;Nx<nEnv;Nx++)
      for(int Ny=0;Ny<nEnv;Ny++){
	Vaa[Nx*nEnv+Ny]+= (1/var)*cov[Nx*numInds+i]*cov[Ny*numInds+i];
      }
    
    for(int x=0;x<nEnv;x++)
      Vab[x]+= (1/var)*Ex[i]*cov[x*numInds+i];
    
    //Vab<-Vab+1/var*Ex[tal]*cbind(Xe[tal,])
    
    if(post[i*3+0]>0.90)
      highWT++;
    if(post[i*3+1]>0.90)
      highHE++;
    if(post[i*3+2]>0.90)
      highHO++;
  }//recursion done
  
  assoc->highWt[0][0][s] = highWT;
  assoc->highHe[0][0][s] = highHE;
  assoc->highHo[0][0][s] = highHO;
  
  for(int i =0; i<numInds;i++)
    sumEx+=Ex[i];
  
  
  
  
  //  double Vab=sumEx/var;
  double Vbb=0;
  for(int i =0; i<numInds;i++)
    Vbb+=(1/var-pow(y[i]-ytilde[i],2)/pow(var,2))*Ex2[i]+pow(y[i]-ytilde[i],2)/pow(var,2)*pow(Ex[i],2);
  
  //I<-Vbb-t(Vab)%*%MASS::ginv(Vaa)%*%Vab
  double workspace[2*nEnv];
  rankProb=angsd::matinv(Vaa, nEnv, nEnv, workspace);
  
  double I =0;
  
  //inv(Vaa)%*%Vab
  double invVaa_Vab[nEnv];
  for(int x=0;x<nEnv;x++)
    invVaa_Vab[x]=0;
  
  //NB! Vaa is now the inverse 
  for(int Nx=0;Nx<nEnv;Nx++)
    for(int Ny=0;Ny<nEnv;Ny++)
      invVaa_Vab[Nx]+=Vaa[Nx*nEnv+Ny]*Vab[Ny];
  
  //I<-t(Vab)%*%MASS::ginv(Vaa)%*%Vab
  for(int x=0;x<nEnv;x++)
    I+=Vab[x]*invVaa_Vab[x];
  //I<-Vbb-t(Vab)%*%MASS::ginv(Vaa)%*%Vab
  //s  fprintf(stderr,"tVab_invVaa_Vab: %f\n",I);  
  I=Vbb-I;


  //the observed varians of the dispersion 
  double Vbs=0;
  for(int i =0; i<numInds;i++)
    Vbs+=Ex[i]*(y[i]-ytilde[i])/pow(var,2);
  
  double Vss=0;
  for(int i =0; i<numInds;i++)
    Vss+=pow(y[i],2)+2*y[i]*ytilde[i];
  Vss=Vss*pow(var,-3)-numInds/(4*M_PI*pow(var,2));

  //fprintf(stderr,"Vbs %f Vss %f\n",Vbs,Vss);
  I=I-pow(Vbs,2)/Vss;

  
  double lrt =pow(U,2)/I;

  int nGeno=0;
  if(highWT >= minHigh)
    nGeno++;
  if(highHE >= minHigh)
    nGeno++;
  if(highHO >= minHigh)
    nGeno++;

  if(nGeno<2)
    lrt=-999;//set_snan(lrt);
  if(freq*numInds*2 < minCount || (1-freq)*numInds*2 < minCount)
    lrt=-999;//set_snan(lrt);      set_snan(lrt);
  if(rankProb!=0)
    lrt=-99;
  //      fprintf(lrtfile,"\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",numInds,ytilde,var,U,Vaa,Vab,Vbb,I,lrt);
  //fprintf(lrtfile,"\t%d\t%f",numInds,lrt);
  /*
  if(verbose>0)
    fprintf(fp,"\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",numInds,ytilde[0],var,U,Vaa[0],Vab[0],Vbb,I,lrt);
    else
    fprintf(fp,"\t%d\t%f",numInds,lrt);
  */

  if((0||lrt>1000||I<-0.01||lrt<0)&&lrt!=-999&&lrt!=-99){//!std::isnan(lrt)){
    for(int i=0 ; i<numInds ;i++) {
      fprintf(stderr,"y: %f\t  yfit: %f \t post %f %f %f\tEx %f %f\tU %f cov: ",y[i],ytilde[i],post[i*3+0],post[i*3+1],post[i*3+2],Ex[i],Ex2[i],U);
      for(int j=0;j<nEnv;j++)
	fprintf(stderr,"%f\t",cov[j*numInds+i]);
      fprintf(stderr,"\n");
 
    }
    for(int j=0;j<pow(nEnv,2);j++)
      fprintf(stderr,"Vaa: %f\t",Vaa[j]); 
    fprintf(stderr,"rank %d\tlrt: %f\t",rankProb,lrt); 
    fprintf(stderr,"\n");                                                                                                            
 
  
    for(int j=0;j<nEnv;j++)
      fprintf(stderr,"Vab: %f\t",Vab[j]);
    fprintf(stderr,"\n");

    fflush(stderr);
    fflush(stdout);
    exit(0);
  }
  
  delete []  Ex;
  delete []  Ex2;
  return lrt;

}
  


double abcAsso::binomScoreEnv(double *post,int numInds, double *y, double *ytilde,double *cov,int nEnv,double freq,assoStruct *assoc,int s){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  double *Ex = angsd::allocArray<double>(numInds,0);
  double *Ex2 = angsd::allocArray<double>(numInds,0);
  int rankProb=0;

  double U=0;
  int highWT=0;
  int highHE=0;
  int highHO=0;
  double sumEx=0;
  double Vaa[nEnv*nEnv];
  for(int x=0;x<nEnv*nEnv;x++)
    Vaa[x]=0;
  double Vab[nEnv];
  for(int x=0;x<nEnv;x++)
    Vab[x]=0;


  for(int i=0 ; i<numInds ;i++) {
    Ex[i]=post[i*3+1]+2*post[i*3+2];
    Ex2[i]=post[i*3+1]+4*post[i*3+2];
    U+=Ex[i]*(y[i]-ytilde[i]);   

  // Vaa<-Vaa+yTilde[i]*(1-yTilde[i])*A[i,]%*%t(A[i,])
    for(int Nx=0;Nx<nEnv;Nx++)
      for(int Ny=0;Ny<nEnv;Ny++){
	Vaa[Nx*nEnv+Ny]+= ytilde[i]*(1-ytilde[i])*cov[Nx*numInds+i]*cov[Ny*numInds+i];
      }

    for(int x=0;x<nEnv;x++)
      Vab[x]+= ytilde[i]*(1-ytilde[i])*Ex[i]*cov[x*numInds+i];

    //Vba<-Vba+yTilde[i]*(1-yTilde[i])*A[i,]*Ex[i]

    if(post[i*3+0]>0.9)
      highWT++;
    if(post[i*3+1]>0.9)
      highHE++;
    if(post[i*3+2]>0.9)
      highHO++;
  }//recursion done
  assoc->highWt[0][0][s] = highWT;
  assoc->highHe[0][0][s] = highHE;
  assoc->highHo[0][0][s] = highHO;
 

    for(int i =0; i<numInds;i++)
      sumEx+=Ex[i];
    // double Vaa=ytilde[0]*(1-ytilde[0])*numInds;
    //    double Vab=ytilde[0]*(1-ytilde[0])*sumEx;
    double Vbb=0;
    for(int i =0; i<numInds;i++)
      Vbb+=(ytilde[i]*(1-ytilde[i])-pow(y[i]-ytilde[i],2))*Ex2[i]+pow(y[i]-ytilde[i],2)*pow(Ex[i],2);

    double workspace[2*nEnv];
    rankProb=angsd::matinv(Vaa, nEnv, nEnv, workspace);
 

    double I =0;

    double invVaa_Vab[nEnv];
    for(int x=0;x<nEnv;x++)
      invVaa_Vab[x]=0;

    for(int Nx=0;Nx<nEnv;Nx++)
      for(int Ny=0;Ny<nEnv;Ny++)
	invVaa_Vab[Nx]+=Vaa[Nx*nEnv+Ny]*Vab[Ny];

    for(int x=0;x<nEnv;x++)
      I+=Vab[x]*invVaa_Vab[x];
    I=Vbb-I;

    double lrt =pow(U,2)/I;

    int nGeno=0;
    if(highWT >= minHigh)
      nGeno++;
    if(highHE >= minHigh)
      nGeno++;
    if(highHO >= minHigh)
      nGeno++;
    if(nGeno<2)
      lrt=-999;
    //freq*numInds*2 is the expected number of minor alleles
    if(freq*numInds*2 < minCount || (1-freq)*numInds*2 < minCount)
      lrt=-999;
    if(rankProb!=0)
      lrt=-99;
    //dispersion matrix has zero corners

    if((lrt>1000||I<-0.01||lrt<0)&&lrt!=-999&&lrt!=-99){//!std::isnan(lrt)){
      for(int i=0 ; i<numInds ;i++) {
	fprintf(stderr,"y: %f\t  post %f %f %f\tEx %f %f\tU %f\n",y[i],post[i*3+0],post[i*3+1],post[i*3+2],Ex[i],Ex2[i],U);
      }
      exit(0);
    }

    delete []  Ex;
    delete []  Ex2;
   
    return lrt;
}



// This method implements the Robust Variance Score statistic described by Derkach et al (2014)
// DOI: 10.1093/bioinformatics/btu196
// Takes as input the posterior genotype probabilities (post), phenotypes (y), and MAF (freq),
// and returns a score statistic (chi-squared, 1df) for the association between phenotype and 
// observed data (through the unobserved genotype variable).
double abcAsso::binomRVScoreEnv(double *post,int numInds, double *y, double *ytilde,double *cov,int nEnv,double freq,assoStruct *assoc,int s){
  // For consistency with the rest of ANGSD.
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  // Array to store the conditional expected genotypes, E(Gij|Dij), for each individual.
  double *Ex = angsd::allocArray<double>(numInds,0);
  
  // Variable to hold the score.
  double U=0;

  // Keep track of the sum of E(Gij|Dij), for use when calculating the variance.
  // Seperately track the cases, controls, and combined values.
  double sumEx[3] = {0};

  // Track the high-confidence heterozygosity and homozygosity rates for filtering.
  int highHE=0;
  int highHO=0;
  int highWT=0;

  // Loop through each individual, retrieving E(Gij|Dij) and using it
  // to update the score statistic.
  for(int i=0 ; i<numInds ;i++) {
    
    // E(Gij|Dij), calculated as ∑gP(Gij=g|Dij), for g=0,1,2.
    Ex[i]=post[i*3+1]+2*post[i*3+2];

    // Score statistic. Y(i) represents the phenotype for this
    // individual, and Ytilde stores the mean (Ncase/N).
    U+=Ex[i]*(y[i]-ytilde[i]);

    // Sum E(Gij|Dij), both separately for cases and controls, and
    // in one overall term.
    sumEx[(int)y[i]]+=Ex[i];
    sumEx[2]+=Ex[i];

    // Track heterozygosity and homosygosity rates for filtering.
    if(post[i*3+0]>0.9)
      highWT++;
    if(post[i*3+1]>0.9)
      highHE++;
    if(post[i*3+2]>0.9)
      highHO++;
  }

  // Store the highHe and highHo rates in the assoc struct, so they can be printed later.
  assoc->highWt[0][0][s] = highWT;
  assoc->highHe[0][0][s] = highHE;
  assoc->highHo[0][0][s] = highHO;

  // Determine how many cases and controls there are.
  double ncase = ytilde[0]*numInds;
  double ncont = numInds - ncase;

  // Calculate the mean E(Gij|Dij) for the different groups (cases, controls, combined).
  double Extilde[3];
  Extilde[0] = sumEx[0]/ncont;
  Extilde[1] = sumEx[1]/ncase;
  Extilde[2] = sumEx[2]/numInds;

  // Array to hold the variance of E(Gij|Dij) for the different groups.
  double var[3] = {0};

  // Calculate Var[E(Gij|Dij)], [∑(Ex - Emean)^2]/N.
  for(int i =0; i<numInds;i++){
    // Update the case or control variance as appropriate.
    var[(int)y[i]]+=pow((Ex[i]-Extilde[(int)y[i]]),2);
    // Update the combined variance.
    var[2]+=pow((Ex[i]-Extilde[2]),2);
  }
  var[0]=var[0]/ncont;
  var[1]=var[1]/ncase;
  var[2]=var[2]/numInds;

  // Calculate the variance of the score, in order to construct the full test statistic.
    // STANDARD VARIANCE -> Var(S) = ∑cases[(1-Ybar)^2]Var(E(Gij|Dij)) + 
    //                               ∑controls[(Ybar)^2]Var(E(Gij|Dij))
    double I = ncase*pow((ncont/numInds),2)*var[2] + ncont*pow(ncase/numInds,2)*var[2];

    // ROBUST VARIANCE -> Var(S) = Ncase[(Ncontrol/N)^2]Var_case(E(Gij|Dij)) + 
    //                             Ncontrol[(Ncase/N)^2]Var_control(E(Gij|Dij))
    double I_RVS = ncase*pow((ncont/numInds),2)*var[1] + ncont*pow(ncase/numInds,2)*var[0];

  // Calculate the test statistic, Tj=(Sj^2)/var(Sj).
  double lrt =pow(U,2)/I;
  double lrt_RVS =pow(U,2)/I_RVS;

  // Keep track of the standard score statistic, to print later for testing/comparison.
  assoc->std_LRT[s] = lrt;
  assoc->rvs_LRT[s] = lrt_RVS;

  // If the MAF was too low, filter out this site.
  if(freq*numInds*2 < minCount || (1-freq)*numInds*2 < minCount){
    lrt_RVS=-999;
  }

  // If the number of confident heterozygous and homozygous calls is too small, filter out this site.
  int nGeno=0;
  if(highWT >= minHigh)
    nGeno++;
  if(highHE >= minHigh)
    nGeno++;
  if(highHO >= minHigh)
    nGeno++;
  if(nGeno<2)
    lrt_RVS=-999;

  // If the score statistic looks otherwise problematic, print out a warning to the user.
  if((lrt_RVS>1000||I<-0.01||lrt_RVS<0)&&lrt_RVS!=-999){
    fprintf(stderr,"WARNING - lrt: %f\tI: %f\tU %f\n",lrt_RVS,I_RVS,U);
  }

  // Tidy up the arrays.
  delete [] Ex;

  // Return the RVS lrt score 
  return lrt_RVS;
}

// This method implements the Robust Variance Score statistic for rare variants, as
// described by Derkach et al (2014). DOI: 10.1093/bioinformatics/btu196
// Takes as input the posterior genotype probabilities (post), phenotypes (y), and MAF (freq)
// for a series of rare variants, and returns a p-value representing the result
// of a CAST burden test for association over all the variants at once.
// NOTE: if the number of rare variants to be included exceeds 50, then it will be
// necessary to increase the chunk size using -nLines on the command line, so that
// all variants are processed together.
double** abcAsso::binomRVScoreEnvRare(funkyPars  *pars,assoStruct *assoc){

  // Seed the random number generator.
  srand (time(NULL));

  // Store the individual score statistics for each variant, to help work out how much each
  // variant may be contributing to the overall burden test.
  double **stat = new double*[ymat.y];

  // Repeat for each phenotype variable.
  for(int yi=0;yi<ymat.y;yi++){

    // Initialise the stat array for this phenotype variable.
    stat[yi] = new double[pars->numSites];
    for(int j=0 ; j<pars->numSites ;j++){
      stat[yi][j]=0;
    }

    // Set up a list of phenotypes and a table (nInd x nSites) of expected genotypes,
    // removing any individuals and sites that fail filtering. Also track the number
    // of cases and controls.
    std::vector<double> y;
    std::vector<std::vector<std::vector<double> > > expected_gt (2, std::vector<std::vector<double> >());

    // A position map to allow results to be saved back to the correct site.
    int *pos = new int[pars->numSites];
    for(int i=0 ; i<pars->nInd ;i++){
      
      // Check we have all the required covariate data for this individual.
      bool covfail=false;
      if(covfile!=NULL){
        for(int ci=0;ci<covmat.y;ci++) {
          if(covmat.matrix[i][ci]==-999)
            covfail=true;
        } 
      }

      // Add data for any individuals passing filtering.
      if(ymat.matrix[i][yi]!=-999 && covfail==false){

        // Add in the phenotype data.
        y.push_back((int)ymat.matrix[i][yi]);

        // Calculate expected genotypes for this individual at every site.
        // TODO: filter out bad sites as well - e.g. those with a lot of uncertainty.
        // Although, possibly this can be done during the first pass, when calculating
        // the baseline statistic - we have a dynamic vector now, so can delete bad sites
        // as we go.
        int keepSites=0;
        for(int j=0;j<pars->numSites;j++){
          
          if(pars->keepSites[j]!=0){

            // Initialise inner vector if necessary.
            if(y.size() == 1){
              std::vector<double> newSite;
              expected_gt.at(0).push_back(newSite);
              expected_gt.at(1).push_back(newSite);
              assoc->highWt[yi][0][j]=assoc->highWt[yi][1][j]=0;
              assoc->highHe[yi][0][j]=assoc->highHe[yi][1][j]=0;
              assoc->highHo[yi][0][j]=assoc->highHo[yi][1][j]=0;
            }

            // Model 1 (DEFAULT) - Additive / Log Additive
            if(model==1){
              expected_gt.at((int)ymat.matrix[i][yi]).at(keepSites).push_back(pars->post[j][i*3+1]+2*pars->post[j][i*3+2]);
            }
            // Model 2 - Dominant (het and non-ref hom contribute equally)
            else if(model==2){
              expected_gt.at((int)ymat.matrix[i][yi]).at(keepSites).push_back(pars->post[j][i*3+1]+pars->post[j][i*3+2]);   
            }
            // Model 3 - Recessive (het and ref hom contribute equally)
            else if(model==3){
              expected_gt.at((int)ymat.matrix[i][yi]).at(keepSites).push_back(pars->post[j][i*3+2]);   
            }

            // Add to the high confidence calls tally.
            if(pars->post[j][i*3+0]>0.9)
              assoc->highWt[yi][(int)ymat.matrix[i][yi]][j]++;
            if(pars->post[j][i*3+1]>0.9)
              assoc->highHe[yi][(int)ymat.matrix[i][yi]][j]++;
            if(pars->post[j][i*3+2]>0.9)
              assoc->highHo[yi][(int)ymat.matrix[i][yi]][j]++;

            pos[keepSites]=j;
            keepSites++;
          }
        }
      } 
    }

    // Calculate the average phenotype and subtract it from every element in the phenotype vector.
    double y_bar = std::accumulate(y.begin(),y.end(),0.0) / y.size();

    // Based on the number of cases and controls, compute a factor to use when combining variances.
    double *F = new double[2];
    F[1] = pow((expected_gt.at(0).at(0).size()/(expected_gt.at(1).at(0).size()+expected_gt.at(0).at(0).size()+0.0)),2); //* expected_gt.at(1).at(0).size()
    F[0] = pow((expected_gt.at(1).at(0).size()/(expected_gt.at(1).at(0).size()+expected_gt.at(0).at(0).size()+0.0)),2); //* expected_gt.at(0).at(0).size()


    // Get the scores for the single site statistics.
    for(int j=0;j<expected_gt.at(0).size();j++){
      stat[yi][pos[j]]=(0-y_bar)*std::accumulate(expected_gt.at(0).at(j).begin(),expected_gt.at(0).at(j).end(),0.0) + (1-y_bar)*std::accumulate(expected_gt.at(1).at(j).begin(),expected_gt.at(1).at(j).end(),0.0);
    }

    // Compute the baseline test statistic for the unpermuted sample.
    double baseline = calculateCAST(y_bar,F,expected_gt);

    // Get the variances for the single site statistics.
    double score = 0;
    for(int j=0;j<expected_gt.at(0).size();j++){
      double var[2] = {0};
      for(int n=0;n<=1;n++){
        var[n]=std::inner_product(expected_gt.at(n).at(j).begin(),expected_gt.at(n).at(j).end(),expected_gt.at(n).at(j).begin(),0.0);
      }
      stat[yi][pos[j]] = pow(stat[yi][pos[j]],2)/(F[1]*var[1] + F[0]*var[0]);
    }

    // Perform the permutation testing.
    int perm=0;
    int sig=0;
    int checkpoint=1000;   
    while(perm<1000000){

      std::vector<std::vector<std::vector<double> > > sample (2, std::vector<std::vector<double> >());

      // Randomly permute cases and controls (separately).
      for(int n=0;n<=1;n++){
        for(int i=0;i<expected_gt.at(n).at(0).size();i++){
            int x = rand() % expected_gt.at(n).at(0).size();
            for(int j=0;j<expected_gt.at(n).size();j++){
              if(i==0){
                std::vector<double> newSite;
                sample.at(n).push_back(newSite);
              }
              sample.at(n).at(j).push_back(expected_gt.at(n).at(j).at(x));
            }
        }
      }

      // Calculate a CAST statistic.
      double cast = calculateCAST(y_bar,F,sample);

      // Increment the sig counter if it is > baseline.
      if(cast>baseline){
        sig++;
      }

      // Increment the perm counter.
      perm++;

      // Periodically check the p-value, and exit if it is clearly not significant.
      if(perm==checkpoint){
        fprintf(stderr,"Yi: %d\tperm: %d\tp: %f\n",yi,perm,(sig+0.0)/perm);
        if(sig > 100)
          break;
        checkpoint*=10;
      }
    }

    // Store how many bootstrap samples were performed.
    assoc->numPerm[yi] = perm;

    // Compute the final p value, as the proportion of bootstrap samples with an absolute test
    // statistic greater than or equal to the original sample.
    assoc->burden[yi] = (sig+0.0)/perm;

  }

  // Return the individual LRT statistics.
  return stat;
}

// This method takes a vector of scores and a covariance matrix, and uses them to compute
// a CAST-style statistic for use in performing rare burden tests.
double abcAsso::calculateCAST(double y_bar, double *F, std::vector<std::vector<std::vector<double> > > &expected_gt){

  double score = 0;
  double cov[2] = {0};
  for(int n=0;n<=1;n++){
    for(int j=0;j<expected_gt.at(0).size();j++){
      
      // Sum up the expected genotypes.
      double sum = std::accumulate(expected_gt.at(n).at(j).begin(),expected_gt.at(n).at(j).end(),0.0);

      // Update the score
      score+=(n-y_bar)*sum;

      // Subtract the mean from each column
      double eg_bar = sum / expected_gt.at(n).at(j).size();
      std::transform(expected_gt.at(n).at(j).begin(),expected_gt.at(n).at(j).end(),expected_gt.at(n).at(j).begin(),std::bind2nd(std::minus<double>(),eg_bar));

      // Update the covariance.
      for(int k=0;k<=j;k++){
        cov[n]+=((k!=j)+1)*std::inner_product(expected_gt.at(n).at(j).begin(),expected_gt.at(n).at(j).end(),expected_gt.at(n).at(k).begin(),0.0);
      }
    }
  }

  // Return the overall variance for the burden test.
  double var = F[1]*cov[1] + F[0]*cov[0];

  // Return the score sum divided by the square root of the variance sum
  double cast=score/sqrt(var);
  if(cast<0){
    cast*=-1;
  }

  //fprintf(stderr,"SCORE: %f\tVAR: %f\tTEST: %f\n",score,var,cast);

  return cast;
}

void abcAsso::printDoAsso(funkyPars *pars){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  freqStruct *freq = (freqStruct *) pars->extras[6];
  assoStruct *assoc= (assoStruct *) pars->extras[index];
  for(int yi=0;yi<ymat.y;yi++){
    bufstr.l=0;
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0){//will skip sites that have been removed      
        continue;
     } 
      if(doAsso==2){
        ksprintf(&bufstr,"%s\t%d\t%c\t%c\t%f\t%d\t%f\t%d/%d/%d\n",header->name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],freq->freq[s],assoc->keepInd[yi][s],assoc->stat[yi][s],assoc->highWt[yi][0][s],assoc->highHe[yi][0][s],assoc->highHo[yi][0][s]);
      }else if(doAsso==4){
        ksprintf(&bufstr,"%s\t%d\t%c\t%c\t%f\t%d\t%f\t%f\t%f\t%d/%d/%d\n",header->name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],freq->freq[s],assoc->keepInd[yi][s],assoc->stat[yi][s],assoc->rvs_LRT[s],assoc->std_LRT[s],assoc->highWt[yi][0][s],assoc->highHe[yi][0][s],assoc->highHo[yi][0][s]);
      }else if(doAsso==5){
        ksprintf(&bufstr,"%s\t%d\t%c\t%c\t%f\t%f\t%f\t%d\t%d/%d/%d\t%d/%d/%d\n",header->name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],freq->freq[s],assoc->stat[yi][s],assoc->burden[yi],assoc->numPerm[yi],assoc->highWt[yi][0][s],assoc->highHe[yi][0][s],assoc->highHo[yi][0][s],assoc->highWt[yi][1][s],assoc->highHe[yi][1][s],assoc->highHo[yi][1][s]);
      }else{
        ksprintf(&bufstr,"%s\t%d\t%c\t%c\t%f\t%f\n",header->name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],freq->freq[s],assoc->stat[yi][s]);
      }
    }
    // gzwrite does not handle 0-sized writes very well (it messes up the checksum, and then the
    // resulting file appears corrupted to gzip, zcat etc). If -sites is being used we may encounter
    // an empty buffer: make sure none of these are passed to gzwrite.
    if(bufstr.s != NULL){
      gzwrite(MultiOutfile[yi],bufstr.s,bufstr.l);
    }
  }
}
