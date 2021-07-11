#include <stdio.h>
#include <math.h>
#include <getopt.h>
#include <stdlib.h>


double logstirling(long n) {
  if (n==0) return 0;
  return lgamma(n+1);
}

double logcnmultiplicity(int *a, int n, int N) {
  int i,j;
  double res=0;
  for (i=0;i<n;i++) {
    res+=logstirling(2*N-1+2*i)-logstirling(N+a[i]+n)-logstirling(N-a[i]+n-1)+log(2*a[i]);
    for (j=0;j<i;j++) {
      res+=log(a[i]*a[i]-a[j]*a[j]);
    }
  }
  return res;
}

double logcndimension(int *a,int n){
  int i,j;
  double res=0;
  for (i=0;i<n;i++) {
    res+=log(2*a[i]);
    for (j=0;j<i;j++) {
      res+=log(a[i]*a[i]-a[j]*a[j])-log(i-j)-log(2*n-i-j);
    }
  }
  return res;
}
int *maxa;
double maxlogmeasure=-10000;

//double logthetameasure(int *a, int k, int N, double theta) {
//  return -theta+k*log(theta)-logstirling(k)+logbnmultiplicity(a,k,N)+logbndimension(a,k);
//}

double logmun(int* a, int n, int N) {
  double res=logcnmultiplicity(a,n,N)+logcndimension(a,n)-2*n*N*log(2);
  return res; 
}

void findmax(int n, int N) {
  long i,j;
  int *a=(int *)malloc(n * sizeof(int));
  for (i=0;i<n;i++) {
    a[i]=i+1;
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
    printf("%d ",a[i]-i-1);
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
  
