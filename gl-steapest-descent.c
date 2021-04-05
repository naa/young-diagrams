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

double loggldimension(int *a,int n){
  int i,j;
  double res=0;
  for (i=0;i<n;i++) {
    for (j=0;j<i;j++) {
      res+=log(a[i]-a[j])-log(i-j);
    }
  }
  return res;
}
int *maxa;
double maxlogmeasure=-10000;

double logmun(int* a, int n, int N) {
  double res=logglmultiplicity(a,n,N)+loggldimension(a,n)-n*N*log(2);
  return res; 
}

void findmax(int n, int N) {
  long i,j;
  int *a=(int *)malloc(n * sizeof(int));
  for (i=0;i<n;i++) {
    a[i]=i;
  }
  double logp;
  double newlogp;
  int cont=1;
  while (cont) {
    cont=0;
    for (i=0;i<n;i++) {
      logp=logmun(a,n,N);
      a[i]+=1;
      newlogp=logmun(a,n,N);
      while(newlogp>logp) {
	cont=1;
	logp=newlogp;
	a[i]+=1;
	newlogp=logmun(a,n,N);
      }
      a[i]-=1;
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
  int opt;  
  while ((opt = getopt(argc, argv, "N:n:")) != -1) {
    switch (opt) {
    case 'N':
      N=atol(optarg);
      break;
    case 'n':
      n = atol(optarg);
      break;
    default: /* '?' */
      fprintf(stderr, "Usage: %s [-N tensor power ] [-n rank] \n",
	      argv[0]);
      exit(EXIT_FAILURE);
    }
  }
  findmax(n,N);
  printf("\n");
  return 0;
}
  
