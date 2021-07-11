#include <stdio.h>
#include <math.h>
#include <getopt.h>
#include <stdlib.h>


double logstirling(long n) {
  if (n==0) return 0;
  return lgamma(n+1);
}

double logglmultiplicity(int *a, int n, int N) {
  int i,j;
  double res=0;
  for (i=0;i<n;i++) {
    res+=logstirling(N+i)-logstirling(a[i])-logstirling(N+n-1-a[i]);
    for (j=0;j<i;j++) {
      res+=log(a[i]-a[j]);
    }
  }
  return res;
}

double loggldimension(int *a,int n,double q){
  int i,j;
  double res=0;
  for (i=0;i<n;i++) {
    for (j=0;j<i;j++) {
      res+=log(1-pow(q,a[i]-a[j]))-log(1-pow(q,i-j));
    }
  }
  return res;
}
int *maxa;
double maxlogmeasure=-10000;

double logmun(int* a,int* ap, int n, int N,double q,double t) {
  //  double res=logglmultiplicity(a,n,N)+loggldimension(a,n);/*-n*N*log(2);*/
  double res=loggldimension(ap,N,t)+loggldimension(a,n,q);/*-n*N*log(2);*/
  return res; 
}

void findmax(int n, int N,double q,double t) {
  long i,j;
  int *a=(int *)malloc(n * sizeof(int));
  int *ap=(int *)malloc(N * sizeof(int));  
  for (i=0;i<n;i++) {
    a[i]=i;
  }
  for (i=0;i<N;i++) {
    ap[i]=n+i;
  }
  printf("%lf %lf\n",loggldimension(ap,N,t),loggldimension(a,n,q));

  double logp;
  double newlogp;
  int cont=1;
  while (cont) {
    cont=0;
    for (i=0;i<n;i++) {
      logp=logmun(a,ap,n,N,q,t);
      //      printf("%lf\n",logp);
      if (a[i]-i>N-1) continue;
      ap[a[i]-i]-=1;
      a[i]+=1;
      newlogp=logmun(a,ap,n,N,q,t);
      //      printf("%lf\n",newlogp);      
      while(newlogp>logp) {
	cont=1;
	logp=newlogp;
	ap[a[i]-i]-=1;	
	a[i]+=1;
	if (a[i]-i>N-1) newlogp=-10000;
	else 	newlogp=logmun(a,ap,n,N,q,t);
      }
      a[i]-=1;
      ap[a[i]-i]+=1;	      
    }
  }
  for (i=n-1;i>=0;i--) {
    printf("%ld ",a[i]-i);
  }
}
  
int main(int argc, char **argv)
{
  long N=1000;
  long n=20;
  double q=1,t=1;
  int opt;  
  while ((opt = getopt(argc, argv, "N:n:q:t:")) != -1) {
    switch (opt) {
    case 'N':
      N=atol(optarg);
      break;
    case 'n':
      n = atol(optarg);
      break;
    case 'q':
      q=atof(optarg);
      break;
    case 't':
      t=atof(optarg);
      break;
    default: /* '?' */
      fprintf(stderr, "Usage: %s [-N tensor power ] [-n rank] [-q q-value] [-t t-value]\n",
	      argv[0]);
      exit(EXIT_FAILURE);
    }
  }
  printf("q:%lf, t:%lf\n",q,t);
  findmax(n,N,q,t);
  printf("\n");
  return 0;
}
  
