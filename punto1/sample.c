#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#define PI 3.14159265358979323846264338327

double g(double x, float mu, float sigma) 
{ 
  return (double) ((1.0/(sqrt(2.0*PI*pow(sigma,2.0)))*(exp(-(pow(x-mu,2.0)))/(2.0*pow(sigma,2.0)))));
}

double* Metropolis(double*list, float mu, float sigma, int N)
{
  int i, j;
  srand48(N);
  for(i=0;i<N;i++)
  {    
      double r = 0.0;
      double p = g(r, mu, sigma);
      while(i < N) 
      {
            double rnext = r + 2*(drand48())-1;
            double pnext = g(rnext, mu, sigma);
            if (pnext >= p) 
            {
                     p = pnext;
                     r = rnext;
            } 
            else 
            {
                    double a = drand48();
                    if(a < (double) (pnext/p)) 
                    {
                            p = pnext;
                            r = rnext;
                     }
           }
                list[i] = r;
                i++;
       }
    }
}


int main(int argc, char **argv)
{
  double p, r;
  FILE *out;
  int N = 1000;
  char filename[128];
  float mu = 0.0;
  float sigma = 1.0;

int j;  
  for (j=0;j<8;j++)
  {
      FILE *out;
      char filename[128];
      double *list=malloc(sizeof(double)*N);
      Metropolis(list,mu, sigma, N);
      sprintf(filename, "sample_%d.dat", j);
    
      if(!(out = fopen(filename, "w"))){
                fprintf(stderr, "Problema abriendo el archivo\n");
                exit(1);
        }
      int i;
      for(i=0;i<N;i++){
                fprintf(out, "%f\n", list[i]);
        }

        fclose(out);
  }
    
    
    return 0;
}

