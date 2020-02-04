#include <Rcpp.h>
using namespace Rcpp;

#define pi_cpp 3.141592653589793238462643383280

// [[Rcpp::export]]
double inorm(double x){
  double inormval=x*R::pnorm(x,0,1,1,0)+R::dnorm(x,0,1,0);
  return(inormval);
}

// [[Rcpp::export]]
double V(double a,double b,double c,double d,double s){
  double V_val=inorm((b-c)/s)-inorm((b-d)/s)-inorm((a-c)/s)+inorm((a-d)/s);
  return(sqrt(2*pi_cpp)*s*s*V_val/(b-a)/(d-c));
}

// [[Rcpp::export]]
NumericVector GCGK(NumericMatrix X,NumericVector parameters){
  Function rankC("rank");
  String strmax("max");
  String strmin("min");
  int n=X.nrow();
  int d=X.ncol();
  NumericMatrix FH(n,d);
  NumericMatrix FL(n,d);
  for(int i=0;i<d;i++){
    NumericVector vec1=X(_,i);
    NumericVector vec2=rankC(vec1,_["ties.method"]=strmax);
    FH(_,i)=vec2/n;
    NumericVector vec3=rankC(vec1,_["ties.method"]=strmin);
    FL(_,i)=(vec3-1)/n;
  }
  int pl=parameters.length();
  NumericVector Tstat(pl);
  for(int xpl=0;xpl<pl;xpl++){
    double s=parameters[xpl];
    double s1=0;
    for(int i=0;i<n;i++){
      for(int j=i;j<n;j++){
        double p1=1;
        for(int k=0;k<d;k++)
          p1=p1*V(FL(i,k),FH(i,k),FL(j,k),FH(j,k),s);
        if(i!=j)
          p1=2*p1;
        s1=s1+p1;
      }
    }
    s1=s1/(n*n);
    double s2=0;
    for(int i=0;i<n;i++){
      double p2=1;
      for(int j=0;j<d;j++)
        p2=p2*V(FL(i,j),FH(i,j),0,1,s);
      s2=s2+p2;
    }
    s2=s2/n;
    double s3=pow(V(0,1,0,1,s),d);
    double T=s1-2*s2+s3;
    if(T<0)
      T=0;
    Tstat[xpl]=sqrt(T);
  }
  return(Tstat);
}

// [[Rcpp::export]]
List GCGKTEST(NumericMatrix X,NumericVector parameters,int iteration){
  Function rankC("rank");
  Function sampleC("sample.int");
  String strmax("max");
  String strmin("min");
  int n=X.nrow();
  int d=X.ncol();
  NumericMatrix FH(n,d);
  NumericMatrix FL(n,d);
  for(int i=0;i<d;i++){
    NumericVector vec1=X(_,i);
    NumericVector vec2=rankC(vec1,_["ties.method"]=strmax);
    FH(_,i)=vec2/n;
    NumericVector vec3=rankC(vec1,_["ties.method"]=strmin);
    FL(_,i)=(vec3-1)/n;
  }
  int pl=parameters.length();
  NumericVector Tstat(pl);
  for(int xpl=0;xpl<pl;xpl++){
    double s=parameters[xpl];
    double s1=0;
    for(int i=0;i<n;i++){
      for(int j=i;j<n;j++){
        double p1=1;
        for(int k=0;k<d;k++)
          p1=p1*V(FL(i,k),FH(i,k),FL(j,k),FH(j,k),s);
        if(i!=j)
          p1=2*p1;
        s1=s1+p1;
      }
    }
    s1=s1/(n*n);
    double s2=0;
    for(int i=0;i<n;i++){
      double p2=1;
      for(int j=0;j<d;j++)
        p2=p2*V(FL(i,j),FH(i,j),0,1,s);
      s2=s2+p2;
    }
    s2=s2/n;
    double s3=pow(V(0,1,0,1,s),d);
    double T=s1-2*s2+s3;
    if(T<0)
      T=0;
    Tstat[xpl]=sqrt(T);
  }
  NumericMatrix ND(iteration,pl);
  for(int itr=0;itr<iteration;itr++){
    NumericMatrix PFH(n,d);
    NumericMatrix PFL(n,d);
    for(int k=0;k<d;k++){
      IntegerVector ind=sampleC(n);
      for(int i=0;i<n;i++){
        int j=ind[i]-1;
        PFH(i,k)=FH(j,k);
        PFL(i,k)=FL(j,k);
      }
    }
    for(int xpl=0;xpl<pl;xpl++){
      double s=parameters[xpl];
      double s1=0;
      for(int i=0;i<n;i++){
        for(int j=i;j<n;j++){
          double p1=1;
          for(int k=0;k<d;k++)
            p1=p1*V(PFL(i,k),PFH(i,k),PFL(j,k),PFH(j,k),s);
          if(i!=j)
            p1=2*p1;
          s1=s1+p1;
        }
      }
      s1=s1/(n*n);
      double s2=0;
      for(int i=0;i<n;i++){
        double p2=1;
        for(int j=0;j<d;j++)
          p2=p2*V(PFL(i,j),PFH(i,j),0,1,s);
        s2=s2+p2;
      }
      s2=s2/n;
      double s3=pow(V(0,1,0,1,s),d);
      double T=s1-2*s2+s3;
      if(T<0)
        T=0;
      ND(itr,xpl)=sqrt(T);
    }
  }
  List L=List::create(Tstat,ND);
  return(L);
}