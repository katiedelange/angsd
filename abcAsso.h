#include "abc.h"

typedef struct{
  double score;
  double variance;
}scoreStruct;

typedef struct{

  double **stat;
  int **keepInd;
  int *highHe;
  int *highWt;
  int *highHo;
  double *afCase;
  double *afCtrl;
  std::vector<std::vector<std::vector<scoreStruct> > > scores;
}assoStruct;



class abcAsso:public abc{
private:
  kstring_t bufstr;
public:
  //none optional stuff
  gzFile *MultiOutfile;
  int doPrint;
  int minCov; //not for users
  int doMaf;
  int dynCov;//not for users
  int doAsso;
  int doPost;
  int GL;
  int sitePerm;  //not for users
  int isBinary;
  int minHigh;
  int minCount;
  int adjust;  //not for users
  int model;
  int numBootstraps;
  int numPermutations;
  void run(funkyPars  *pars);
  void print(funkyPars *pars);  
  void clean(funkyPars *pars);  
  void getOptions(argStruct *arguments);
  void printArg(FILE *argFile);

  abcAsso(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcAsso();
  //other stuff
  char *covfile;
  char *yfile;
  
  angsd::Matrix<double> ymat;
  angsd::Matrix<double> covmat;
  void scoreAsso(funkyPars  *pars,assoStruct *assoc);
  void frequencyAsso(funkyPars  *pars,assoStruct *assoc);
  double doAssociation(funkyPars *pars,double *post,double *y,int keepInd,int *keepList,double freq,int s,assoStruct *assoc);
  void getFit(double *res,double *Y,double *covMatrix,int nInd,int nEnv);
  void getFitBin(double *res,double *Y,double *covMatrix,int nInd,int nEnv);
  double normScoreEnv(double *post,int numInds, double *y, double *ytilde,double *cov,int nEnv,double freq,assoStruct *assoc,int s);
  double binomScoreEnv(double *post,int numInds, double *y, double *ytilde,double *cov,int nEnv,double freq,assoStruct *assoc,int s);
  std::vector<std::vector<scoreStruct> > doAdjustedAssociation(funkyPars *pars, double *phenotypes, int *keepList, assoStruct *assoc);
  void computeScore(int sample, int *pos, std::vector<std::vector<double> > alphaN, std::vector<std::vector<std::vector<double> > > e_gij_dij, std::vector<std::vector<scoreStruct> > &scores);
  void computeVariance(int sample, int *pos, std::vector<std::vector<double> > alphaN, std::vector<std::vector<std::vector<double> > > e_gij_dij, std::vector<std::vector<scoreStruct> > &scores);
  void computeScoreSums(int sample, int *pos, std::vector<std::vector<double> > alphaN, std::vector<std::vector<std::vector<double> > > e_gij_dij, std::vector<std::vector<scoreStruct> > &scores);
  void computeVarianceCovarianceMatrix(int sample, int *pos, std::vector<std::vector<double> > alphaN, std::vector<std::vector<std::vector<double> > > e_gij_dij, std::vector<std::vector<scoreStruct> > &scores);
  void printDoAsso(funkyPars *pars);
};
