 /*
  thorfinn thorfinn@binf.ku.dk dec17 2012 
  has modified so no sites[].chromo etc are used
  
  anders albrecht@binf.ku.dk made this.


  logistic regression changed 612new/613

  part of angsd

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
  fprintf(argFile,"\t3: Adjusted Score Test\n");
  fprintf(argFile,"\t4: Adjusted Score Burden Test\n");
  fprintf(argFile,"  Frequency Test Options:\n");
  fprintf(argFile,"\t-yBin\t\t%s\t(File containing disease status)\t\n\n",yfile);
  fprintf(argFile,"  Score Test Options:\n");
  fprintf(argFile,"\t-yBin\t\t%s\t(File containing disease status)\n",yfile);
  fprintf(argFile,"\t-yQuant\t\t%s\t(File containing phenotypes)\n",yfile);
  fprintf(argFile,"\t-minHigh\t%d\t(Require atleast minHigh number of high credible genotypes)\n",minHigh);
  fprintf(argFile,"\t-minCount\t%d\t(Require this number of minor alleles, estimated from MAF)\n",minCount);
  fprintf(argFile,"\t-cov\t\t%s\t(File containing additional covariates)\n",covfile);
  fprintf(argFile,"\t-model\t%d\n",model);
  fprintf(argFile,"\t1: Additive/Log-Additive (Default)\n");
  fprintf(argFile,"\t2: Dominant\n");
  fprintf(argFile,"\t3: Recessive\n\n");
  fprintf(argFile,"\t-numBootstraps\t%d\t(The number of bootstrap samples to generate for burden testing)\n",numBootstraps);
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
  sitePerm=angsd::getArg("-sitePerm",sitePerm,arguments);
  GL=angsd::getArg("-GL",GL,arguments);
  covfile=angsd::getArg("-cov",covfile,arguments);
  doPost=angsd::getArg("-doPost",doPost,arguments);
  yfile=angsd::getArg("-yBin",yfile,arguments);
  if(yfile!=NULL)
    isBinary=1;
  yfile=angsd::getArg("-yQuant",yfile,arguments);
  numBootstraps=angsd::getArg("-numBootstraps",numBootstraps,arguments);
  numPermutations=numBootstraps;
  if(numBootstraps == -1){
    numBootstraps=1000000;
  }

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
  dynCov=0;//not for users
  minCov=5;//not for users
  adjust=1;//not for users
  doMaf=0;
  numBootstraps=0;
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
  int numOutfiles = ymat.y;
  if(doAsso >=3 && numBootstraps != 0)
    numOutfiles *= 3;
  MultiOutfile = new gzFile[numOutfiles];
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
  for(int i=0;i<numOutfiles;i++){
    int p=i;
    if(i>=ymat.y){
      postfix=".score";
      if(i>=ymat.y*2)
        postfix=".var";
      p=i%ymat.y;
    }
    char ary[5000];
    snprintf(ary,5000,"%s%d.gz",postfix,p);
    MultiOutfile[i] = Z_NULL;
    MultiOutfile[i] = aio::openFileGz(outfiles,ary,GZOPT);
  }

  //print header
  for(int yi=0;yi<numOutfiles;yi++){
    if(doAsso==2)
      gzprintf(MultiOutfile[yi],"Chromosome\tPosition\tMajor\tMinor\tFrequency\tN\tLRT\thigh_WT/HE/HO\n");
    else if(doAsso>=3){
      if(yi < ymat.y){
        gzprintf(MultiOutfile[yi],"Chromosome\tPosition\tMajor\tMinor\tFrequency\tN\tLRT\thigh_WT/HE/HO\tAF_case\tAF_ctrl\n");
      }
      else{
        gzprintf(MultiOutfile[yi],"Chromosome\tPosition\tMajor\tMinor\tFrequency\t");
        if(doAsso==4){
          gzprintf(MultiOutfile[yi],"ss0\t");
        }
        for(int b=0;b<=numBootstraps;b++){
          gzprintf(MultiOutfile[yi],"s%d\t",b);
        }
        gzprintf(MultiOutfile[yi],"\n");
      }
    }
    else
      gzprintf(MultiOutfile[yi],"Chromosome\tPosition\tMajor\tMinor\tFrequency\tLRT\n");
  }
}


abcAsso::~abcAsso(){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);


  if(doAsso==0)
    return;

  int numOutfiles = ymat.y;
  if(doAsso >=3 && numBootstraps != 0)
    numOutfiles *= 3;
  for(int i=0;i<numOutfiles;i++)
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

  if(assoc->stat!=NULL)
    for(int yi=0;yi<ymat.y;yi++)
      delete[] assoc->stat[yi];

  delete[] assoc->stat;
  
  delete[]  assoc->highWt;
  delete[]  assoc->highHe;
  delete[]  assoc->highHo;
  delete[]  assoc->afCase;
  delete[]  assoc->afCtrl;

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
  assoc->afCase = NULL;
  assoc->afCtrl = NULL;
  
  return assoc;
}


void abcAsso::run(funkyPars *pars){


  if(doAsso==0)
    return;
  
  assoStruct *assoc = allocAssoStruct();
  pars->extras[index] = assoc;

  if(doAsso==1){
    frequencyAsso(pars,assoc);
  }
  else if(doAsso>=2){
    assoc->highWt=new int[pars->numSites];
    assoc->highHe=new int[pars->numSites];
    assoc->highHo=new int[pars->numSites];
    assoc->afCase=new double[pars->numSites];
    assoc->afCtrl=new double[pars->numSites];

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
    fprintf(stderr,"Only one phenotype allowed for -doAsso 1\n");
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
  std::vector<std::vector<std::vector<scoreStruct> > > scores (ymat.y, std::vector<std::vector<scoreStruct> >());

  for(int yi=0;yi<ymat.y;yi++){
    stat[yi] = new double[pars->numSites];
    keepInd[yi]= new int[pars->numSites];
  }
  
  // Pull out the already-calculated frequency information.
  freqStruct *freq = (freqStruct *) pars->extras[6];

  // Loop over each phenotype.
  for(int yi=0;yi<ymat.y;yi++) { 

    fprintf(stderr,"Processing phenotype %d\n",yi);

    // Create a new phenotype status array, and then populate it 
    // using information from the phenotypes file.
    double *y = new double[pars->nInd];
    for(int i=0 ; i<pars->nInd ;i++)
      y[i]=ymat.matrix[i][yi]; 

    // Determine which individuals will be kept for the given
    // phenotype and covariates.
    int *keepList = new int[pars->nInd];
    int keptInd=0;
    for(int i=0 ; i<pars->nInd ;i++) {
      keepList[i]=1;
      // If the phenotype is unknown, remove this individual.
      if(ymat.matrix[i][yi]==-999)
        keepList[i]=0;
      // If we have any covariates supplied, check each in turn.
      if(covfile!=NULL)
        for(int ci=0;ci<covmat.y;ci++) {
          // Any missing covariate information leads to the
          // individual being excluded.
          if(covmat.matrix[i][ci]==-999)
            keepList[i]=0;
        }
      // If there were no grounds to exclude this individual,
      // add one to the count of kept individuals for this
      // site and phenotype.
      if(keepList[i]==1)
        keptInd++;
    }

    // If we are performing the original score test, perform the association test site by site.
    if(doAsso == 2){
      for(int s=0;s<pars->numSites;s++){
        
        // If this site was previously filtered out, skip over it.
        if(pars->keepSites[s]==0)
          continue;

        // Because we apparently care about storing this as a large array in memory, even though
        // the value will be the same at every site. I assume it is used elsewhere in ANGSD,
        // hence I keep it in the original format.
        keepInd[yi][s] = keptInd;

        // Do the actual association!
        stat[yi][s]=doAssociation(pars,pars->post[s],y,keepInd[yi][s],keepList,freq->freq[s],s,assoc);
      }
    }
    // If we are performing the adjusted score test (either as a single site or as a burden), send
    // the complete dataset off for processing.
    else if(doAsso == 3 || doAsso == 4){
      
      scores[yi]=doAdjustedAssociation(pars,y,keepList,assoc);
      
      // Compute the single-site test statistic, Tj = (Sj^2)/var(Sj), which is chi-squared with one degree. 
      // Add in an artificial direction for the association, to assist with analysis of results.
      for(int s=0;s<pars->numSites;s++){
        int direction = 1;
        if(scores[yi][s][0].score < 0)
          direction = -1;
        stat[yi][s]=(pow(scores[yi][s][0].score,2)/scores[yi][s][0].variance)*direction;
      }
    }

    //cleanup 
    delete [] y;
    delete [] keepList;
  }

  // Save the results for printing.
  assoc->stat=stat;
  assoc->scores=scores;
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
    stat = binomScoreEnv(post,keepInd,y,yfit,covMatrix,nEnv,freq,assoc,s);
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
  
  assoc->highWt[s] = highWT;
  assoc->highHe[s] = highHE;
  assoc->highHo[s] = highHO;
  
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
  assoc->highWt[s] = highWT;
  assoc->highHe[s] = highHE;
  assoc->highHo[s] = highHO;
 

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



// This method sets up the adjusted score test data, for both burden and single-site tests.
// It takes in a set of genotype probabilities, calculates adjusted expected genotypes from
// them, and performs bootstrap resampling when requested.
// This method implements the Robust Variance Score statistic described by Derkach et al (2014)
// DOI: 10.1093/bioinformatics/btu196, plus an additional AF estimate adjustment (unpublished).
// Takes as input the posterior genotype probabilities (post), phenotypes (y), and MAF (freq),
// and returns a score statistic (chi-squared, 1df) for the association between phenotype and 
// observed data (through the unobserved genotype variable).
std::vector<std::vector<scoreStruct> > abcAsso::doAdjustedAssociation(funkyPars *pars, double *y, int *keepList, assoStruct *assoc){

  // A matrix containing the scores and variances for each site, for the original sample
  // plus each of the bootstrap samples. For the burden test include an additional
  // list for the single site test results on the original sample.
  std::vector<std::vector<scoreStruct> > scores (pars->numSites,std::vector<scoreStruct>(doAsso-2,scoreStruct()));

  // A matrix to store the expected genotypes E(Gij|Dij) for each individual, at each site, split 
  // into cases and controls.
  std::vector<std::vector<std::vector<double> > > e_gij_dij (2, std::vector<std::vector<double> >());
  
  // A list to store Var(Gij=g|Dij) for each individual, at each site, split into cases
  // and controls.
  std::vector<std::vector<double> > v_gj_dj (2, std::vector<double> ());

  // A list of the alpha (IMPUTE2 info) scores x the number of samples, for each site.
  // This is calculated separately for cases and controls.
  std::vector<std::vector<double> > alphaN (2, std::vector<double> ());

  // Track how many sites we are actually keeping.
  int kept=0;

  // A position map to allow results to be saved back to the correct site.
  int *pos = new int[pars->numSites];

  // For each site j:
  for(int j=0;j<pars->numSites;j++){

    // If this site was previously filtered out, skip over it.
    if(pars->keepSites[j]==0)
      continue;

    // Add a new sites list, for both cases and controls.
    std::vector<double> newSite;
    e_gij_dij.at(0).push_back(newSite);
    e_gij_dij.at(1).push_back(newSite);

    // Initialise the Var(Gj|Dj) to 0 for this site (for both cases and controls).
    v_gj_dj.at(0).push_back(0.0);
    v_gj_dj.at(1).push_back(0.0);

    // Keep track of the number of high confidence genotype probabilities in the sample, at each site.
    int highWT=0;
    int highHE=0;
    int highHO=0;

    // Sum up the genotype probabilities to generate the expected genotype E(Gij|Dij)
    // for each individual: E(Gij|Dij) = ∑gP(Gij=g|Dij), for g=0,1,2.
    // Also compute the variance Var(Gij|Dij) for use in determing the amount of information
    // about the population allele frequency provided by the genotype probabilities:
    // E(Gij^2|Dij) = ∑(g^2)P(Gij=g|Dij), for g=0,1,2.
    // Var(Gij=g|Dij) = E(Gij^2|Dij) - E(Gij|Dij)^2. Sum this up for all individuals, to give the term v_gj_dj.
    for(int i=0 ; i<pars->nInd ;i++){
      if(keepList[i] == 1){
        e_gij_dij.at((int)y[i]).at(kept).push_back(pars->post[j][i*3+1]+2*pars->post[j][i*3+2]);
        v_gj_dj.at((int)y[i]).at(kept) += (pars->post[j][i*3+1]+4*pars->post[j][i*3+2]) - pow(pars->post[j][i*3+1]+2*pars->post[j][i*3+2],2);

        // Track the number of high confidence genotypes.
        if(pars->post[j][i*3+0]>0.9)
          highWT++;
        if(pars->post[j][i*3+1]>0.9)
          highHE++;
        if(pars->post[j][i*3+2]>0.9)
          highHO++;
      }
    }

    // Store the high confidence genotype counts.
    assoc->highWt[j] = highWT;
    assoc->highHe[j] = highHE;
    assoc->highHo[j] = highHO;

    // Use the two terms E(Gij|Dij) and Var(Gij|Dij) to compute alpha, separately for cases
    // and controls. Alpha is essentially the IMPUTE2 info score, and it describes the amount
    // of missing allele frequency information in the genotype probabilities, such that the
    // amount of available data is equivalent to a set of perfectly observed genotypes in a
    // sample of size aN. It is expected to be smaller for the lower coverage group (due to
    // reduced sensitivity in low coverage data): thus, do this separately for cases and controls.
    for(int n = 0; n <=1; n++){

      // Get the full sample size.
      int N = e_gij_dij.at(n).at(kept).size();

      // Compute the allele frequency estimate, ∑E(Gij|Dij)/2N.
      double af = std::accumulate(e_gij_dij.at(n).at(kept).begin(),e_gij_dij.at(n).at(kept).end(),0.0) / (2 * N);
      
      // Store the allele frequencies in the cases and controls for printing.
      if(n==0)
        assoc->afCtrl[j]=af;
      else
        assoc->afCase[j]=af;

      // This adjustment does not work if the allele frequency estimate is 0. Similarly,
      // for singletons there is no variance in the genotype probabiities, so a useful
      // alpha value will not be obtained.
      if(af > (1.0/(2 * N))){

        // The complete information, had we perfectly observed genotypes - 2N/(AF(1-AF)).
        double complete_info = (2 * N) / (af*(1-af));

        // The missing information, as we only have genotype probabilities - Var(Gij|Dij)/(AF^2((1-AF)^2)).
        double missing_info = v_gj_dj.at(n).at(kept) / (pow(af,2)*pow(1-af,2));

        // Alpha is the ratio of the observed information (complete - missing) to the complete information.
        // This is multiplied by N (the original sample size) to give an adjusted sample size, aN.
        alphaN[n].push_back(((complete_info - missing_info) / complete_info) * N);
      }
      else{
        alphaN[n].push_back(N);
      }
    }

    // Store the variant this particular keepSite position actually refers to.
    pos[kept]=j;

    // Update the 'kept' count.
    kept++;
  }

  // Compute the score Sj using the adjusted E(Gij|Dij) values.
  computeScore(0,pos,alphaN,e_gij_dij,scores);

  // Calculate the variance for each site.
  computeVariance(0,pos,alphaN,e_gij_dij,scores);

  // If doAsso == 4, calculate the variance-covariance matrix for the burden test using the
  // E(Gij|Dij) values across all sites. Sum up all the scores.
  if(doAsso == 4){
    computeScoreSums(0,pos,alphaN,e_gij_dij,scores);
    computeVarianceCovarianceMatrix(0,pos,alphaN,e_gij_dij,scores);
  }

  // Perform the requested number of bootstrap resampling steps: the default is to do no
  // resampling. If a positive value is given using the flag -numBootstraps, repeat that
  // many times. If the value -1 is passed, resample using adaptive permutation.
  // 
  if(numPermutations > 0 || numPermutations == -1){

    // Centre the E(Gij|Dij) values around their means (calculated separately for cases
    // and controls). This reduces the dimension of the difference between the groups
    // to variance only, allowing for resampling to be performed separately within cases
    // and controls, rather than across phenotype labels (which is not possible if read
    // depth is perfectly correlated with phenotype).
    // Note that we use the adjusted sample size to calculate the mean E(Gij|Dij) values.
    for(int n=0;n<=1;n++){
      for(int j=0;j<e_gij_dij.at(0).size();j++){
        double sum = std::accumulate(e_gij_dij.at(n).at(j).begin(),e_gij_dij.at(n).at(j).end(),0.0);
        double eg_bar = sum / alphaN[n][j];
        std::transform(e_gij_dij.at(n).at(j).begin(),e_gij_dij.at(n).at(j).end(),e_gij_dij.at(n).at(j).begin(),std::bind2nd(std::minus<double>(),eg_bar));
      }
    }

    // Seed the random number generator with a constant: I want the same sample number to
    // correspond to an identical set of permuted cases/controls every time this is run.
    srand (42);

    // Generate numBootstrap bootstrap samples from the set of centred E(Gij|Dij), and
    // compute the scores and variances using the appropriate method (doAsso == 3 or 4).
    int perm = 1;
    int checkpoint=1000;   
    
    // If we're doing adaptive permutation, we'll need to keep an eye on the significance.
    int sig = 0;
    double base = std::abs(scores[0][1].score/sqrt(scores[0][1].variance));

    // Currently, 1000000 is the maximum number of permutations that will be carried out.
    while(perm<=1000000){

      // Create a bootstrap sample by randomly permuting cases and controls (separately).
      std::vector<std::vector<std::vector<double> > > sample (2, std::vector<std::vector<double> >());
      for(int n=0;n<=1;n++){
        // Loop through nCase (or nControl) number of times.
        for(int i=0;i<e_gij_dij.at(n).at(0).size();i++){

            // Randomly select and individual (with replacement).
            int x = rand() % e_gij_dij.at(n).at(0).size();

            // Add the E(Gij|Dij) data for this individual (at all sites)
            // to the sample set of expected genotypes.
            for(int j=0;j<e_gij_dij.at(n).size();j++){

              // If this is the first selected individual, we'll need to set up the lists.
              if(i==0){
                std::vector<double> newSite;
                sample.at(n).push_back(newSite);
              }
              sample.at(n).at(j).push_back(e_gij_dij.at(n).at(j).at(x));
            }
        }
      }

      // Add a new sample scoreStruct to the list. This involves some unnecessary memory usage,
      // and an unfortunate extra loop through the number of sites, because I need to create one
      // for every site in order to fit in with the existing ANGSD structure - even though I don't
      // end up calculating values for every position.
      for(int j=0;j<pars->numSites;j++){
        scores[j].push_back(scoreStruct());
      }

      // Calculate the score and variance values (using the selected method).
      if(doAsso == 3){
        computeScore(perm,pos,alphaN,sample,scores);
        computeVariance(perm,pos,alphaN,sample,scores);
      }
      else if(doAsso == 4){
        computeScoreSums(perm,pos,alphaN,sample,scores);
        computeVarianceCovarianceMatrix(perm,pos,alphaN,sample,scores);
      }

      // If we are doing adaptive permutation, we'll need to keep an eye on the signficance.
      if(numPermutations  == -1){
        
        // Increment the sig counter if it is > baseline.
        double cast = std::abs(scores[0][perm+1].score/sqrt(scores[0][perm+1].variance));
        if(cast > base)
          sig++;

        // Periodically check the p-value, and exit if it is clearly not significant.
        if(perm==checkpoint){
          fprintf(stderr,"Perm: %d\tP-value: %f\n",perm,(sig*1.0/perm));
          if(sig > 100){
            numBootstraps=perm;
            break;
          }
          else
            checkpoint*=10;
        }
      }
      else if(perm == numBootstraps){
        break;
      }

      // Increment the perm counter.
      perm++;
    }
  }

  // Return a matrix of size num_sitesX(numBootstraps+1) for the single site test, or 
  // num_sitesX(numBootstraps+2) for the burden test (so that an additional single site column
  // for the unpermuted sample can be added). Each element of the matrix is a struct containing
  // two values: the score and the variance. 
  return scores;
}

// This method computes the score at all sites, using a matrix of expected genotypes.
void abcAsso::computeScore(int sample, int *pos, std::vector<std::vector<double> > alphaN, std::vector<std::vector<std::vector<double> > > e_gij_dij, std::vector<std::vector<scoreStruct> > &scores){

  // Process each site separately.
  for(int j=0;j<e_gij_dij.at(0).size();j++){

    // Compute the mean phenotype at this site using the adjusted sample sizes.
    // NOTE: Currently, it is a waste of processor time to compute this for every variant. Because
    // we do not exclude any individuals based on low MAX(genotype probability), this number will be
    // the same for every site, and every bootstrap permutation. However, we may wish to include some
    // filtering at a later stage, so this is left here to make updating the code simpler.
    double y_bar = alphaN[1][j] / (alphaN[0][j] + alphaN[1][j] + 0.0);

    // Calculate the score, ∑(Yi-Ybar)E(Gij|Dij).
    scores[pos[j]][sample].score = (0-y_bar)*std::accumulate(e_gij_dij.at(0).at(j).begin(),e_gij_dij.at(0).at(j).end(),0.0) + (1-y_bar)*std::accumulate(e_gij_dij.at(1).at(j).begin(),e_gij_dij.at(1).at(j).end(),0.0);
  }
}

// This method computes the variance of the adjusted score, var(Sj) using only a single site. It takes
// as input a matrix of expected genotypes, split into cases and controls, and a factor F
// describing the relative number of individuals in each group. 
void abcAsso::computeVariance(int sample, int *pos, std::vector<std::vector<double> > alphaN, std::vector<std::vector<std::vector<double> > > e_gij_dij, std::vector<std::vector<scoreStruct> > &scores){

  // Loop through each site.
  for(int j=0;j<e_gij_dij.at(0).size();j++){

    // It is expected that the variance will be greater in the lower coverage group
    // (as the genotypes are more uncertain in low coverage data), so the variance is
    // computed separately for each group and then combined, to avoid underestimation.
    scores[pos[j]][sample].variance = 0.0;
    for(int n=0;n<=1;n++){

      // Calculate the mean expected genotype at this site using the adjusted sample size.
      double sum = std::accumulate(e_gij_dij.at(n).at(j).begin(),e_gij_dij.at(n).at(j).end(),0.0);      
      double eg_bar = sum / alphaN[n][j];

      // Subtract this mean expected genotype from the E(Gij|Dij) for each individual.
      std::transform(e_gij_dij.at(n).at(j).begin(),e_gij_dij.at(n).at(j).end(),e_gij_dij.at(n).at(j).begin(),std::bind2nd(std::minus<double>(),eg_bar));

      // Add in the variance component (either cases or controls), using adjusted sample sizes.
      // var(Sj) = ∑(1-Y_bar)^2(Var_cases(E(Gij|Dij))) + ∑(Y_bar)^2(Var_controls(E(Gij|Dij))).
      double var = std::inner_product(e_gij_dij.at(n).at(j).begin(),e_gij_dij.at(n).at(j).end(),e_gij_dij.at(n).at(j).begin(),0.0);
      scores[pos[j]][sample].variance += pow(alphaN[std::abs(n-1)][j]/(alphaN[0][j]+alphaN[1][j]),2) * var;
    }
  }
}

// This method computes the sum of the scores at all sites, using a matrix of expected genotypes.
void abcAsso::computeScoreSums(int sample, int *pos, std::vector<std::vector<double> > alphaN, std::vector<std::vector<std::vector<double> > > e_gij_dij, std::vector<std::vector<scoreStruct> > &scores){

  // Get the scores for the single site statistics, removing any sites with no high-confidence variant calls in the process.
  double sum = 0.0;
  for(int j=0;j<e_gij_dij.at(0).size();j++){

    // Compute the mean phenotype at this site using the adjusted sample sizes.
    double y_bar = alphaN[1][j] / (alphaN[0][j] + alphaN[1][j] + 0.0);

    // Calculate the score for this site, and add it to the sum.
    sum += (0-y_bar)*std::accumulate(e_gij_dij.at(0).at(j).begin(),e_gij_dij.at(0).at(j).end(),0.0) + (1-y_bar)*std::accumulate(e_gij_dij.at(1).at(j).begin(),e_gij_dij.at(1).at(j).end(),0.0);
  }

  // NOTE: Currently this is unecessarily memory-intensive. It creates a full matrix of the same size as would be
  // required to do single-site permutations, but only fills in the first row and first column. The code should be
  // refactored so this doesn't happen!
  scores[0][sample+1].score = sum;
}

// This method takes expected genotypes from a number of sites, and computes the complete
// variance-covariance matrix for those sites. It takes as input a matrix of expected genotypes, 
// split into cases and controls, for a number of sites, and a factor F that describes the 
// relative number of individuals in each group. 
// Please note that this is the processing bottleneck, as calculating covariances is
// very computationally expensive. Therefore, this burden test should ideally only
// be performed after first evaluating the aggregation of single site results to filter
// out clearly insignificant burden regions.
void abcAsso::computeVarianceCovarianceMatrix(int sample, int *pos, std::vector<std::vector<double> > alphaN, std::vector<std::vector<std::vector<double> > > e_gij_dij, std::vector<std::vector<scoreStruct> > &scores){

  for(int j=0;j<e_gij_dij.at(0).size();j++){
    for(int n=0;n<=1;n++){

      // Calculate the mean expected genotype at this site using the adjusted sample size.
      double sum = std::accumulate(e_gij_dij.at(n).at(j).begin(),e_gij_dij.at(n).at(j).end(),0.0);      
      double eg_bar = sum / alphaN[n][j];

      // Subtract this mean expected genotype from the E(Gij|Dij) for each individual.
      std::transform(e_gij_dij.at(n).at(j).begin(),e_gij_dij.at(n).at(j).end(),e_gij_dij.at(n).at(j).begin(),std::bind2nd(std::minus<double>(),eg_bar));

      // Update the covariance.
      for(int k=0;k<=j;k++){

        // Add in the variance component (either cases or controls), using adjusted sample sizes.
        // var(Sj) = Ncase*(Nctrl/N)^2 * Vcase + Nctrl*(Ncase/N)^2 * Vctrl.
        double covar = std::inner_product(e_gij_dij.at(n).at(j).begin(),e_gij_dij.at(n).at(j).end(),e_gij_dij.at(n).at(k).begin(),0.0);
        double y_factor = pow((alphaN[std::abs(n-1)][j]+alphaN[std::abs(n-1)][k])/(alphaN[0][j]+alphaN[1][j]+alphaN[0][k]+alphaN[1][k]),2);

        // If this is one of the covariances, add it twice. The variances (diagonals on the matrix)
        // are only included once.
        scores[0][sample+1].variance += ((k!=j)+1) * y_factor * covar;
      }
    }
  }
}

void abcAsso::printDoAsso(funkyPars *pars){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  freqStruct *freq = (freqStruct *) pars->extras[6];
  assoStruct *assoc= (assoStruct *) pars->extras[index];

  int numOutfiles = ymat.y;
  if(doAsso >=3 && numBootstraps != 0)
    numOutfiles *=3;

  for(int yi=0;yi<numOutfiles;yi++){
    bufstr.l=0;

    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0){//will skip sites that have been removed      
	continue;
     } 
      if(doAsso==2){
	ksprintf(&bufstr,"%s\t%d\t%c\t%c\t%f\t%d\t%f\t%d/%d/%d\n",header->name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],freq->freq[s],assoc->keepInd[yi][s],assoc->stat[yi][s],assoc->highWt[s],assoc->highHe[s],assoc->highHo[s]);

      }
      else if(doAsso>=3){
        if(yi<ymat.y){
          ksprintf(&bufstr,"%s\t%d\t%c\t%c\t%f\t%d\t%f\t%d/%d/%d\t%f\t%f\n",header->name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],freq->freq[s],assoc->keepInd[yi][s],assoc->stat[yi][s],assoc->highWt[s],assoc->highHe[s],assoc->highHo[s],assoc->afCase[s],assoc->afCtrl[s]);
        }
        else if(yi<ymat.y*2){
          ksprintf(&bufstr,"%s\t%d\t%c\t%c\t%f\t%f\t%f\t",header->name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],freq->freq[s],assoc->afCase[s],assoc->afCtrl[s]);
          for(int b=0;b<numBootstraps+(doAsso-2);b++){
            if(s==0 || b == 0 || doAsso ==3)
              ksprintf(&bufstr,"%f\t",assoc->scores[yi%ymat.y][s][b].score);
          }
          ksprintf(&bufstr,"\n");  
        }
        else{
          ksprintf(&bufstr,"%s\t%d\t%c\t%c\t%f\t%f\t%f\t",header->name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],freq->freq[s],assoc->afCase[s],assoc->afCtrl[s]);
          for(int b=0;b<numBootstraps+(doAsso-2);b++){
            if(s==0 || b == 0 || doAsso ==3)
              ksprintf(&bufstr,"%f\t",assoc->scores[yi%ymat.y][s][b].variance);
          }
          ksprintf(&bufstr,"\n"); 
        }
      }
      else{
	ksprintf(&bufstr,"%s\t%d\t%c\t%c\t%f\t%f\n",header->name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],freq->freq[s],assoc->stat[yi][s]);

      }
    }

    // Add an extra information line to the end of the burden test results.
    if(doAsso == 4 && yi < ymat.y){

      // Generate a p-value for the burden test.
      double p_value = 0;
      double base = std::abs(assoc->scores[yi][0][1].score/sqrt(assoc->scores[yi][0][1].variance));
      for(int b=2;b<assoc->scores[yi][0].size();b++){
        double cast = std::abs(assoc->scores[yi][0][b].score/sqrt(assoc->scores[yi][0][b].variance));
        if(cast > base)
          p_value++;
      }
      p_value /= (assoc->scores[yi][0].size()-2);

      // Write it to file.
      ksprintf(&bufstr,"P-value for the complete burden test, after %d permutations: %f\n",(assoc->scores[yi][0].size()-2),p_value);
    }

    // gzwrite does not handle 0-sized writes very well (it messes up the checksum, and then the
    // resulting file appears corrupted to gzip, zcat etc). If -sites is being used we may encounter
    // an empty buffer: make sure none of these are passed to gzwrite.
    if(bufstr.s != NULL){
      gzwrite(MultiOutfile[yi],bufstr.s,bufstr.l);
    }
  }
}
