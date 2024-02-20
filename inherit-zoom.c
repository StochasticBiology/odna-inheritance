#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND drand48()    // random number between 0 and 1
#define NPOP 100         // population size
#define NGEN 2000        // number of generations in simulation
#define NSAMP 100        // number of samples to run

#define _MEANOUTPUT

// structure to store an individual's genetic makeup -- number of genomes of type a, b, c
typedef struct {
  int a, b, c;
} Ind;

// fitness of individual I in environment type env at time t
double fitness(Ind I, int env, int t, double scale, double penalty)
{
  double h, fitness;

  if(env == 0)
     fitness = (I.a+scale*I.b);
  else
     fitness = (t % env == 0 ? I.a+scale*I.b : I.b+scale*I.a);
  
  h = (I.b < I.a ? I.b : I.a);

  fitness -= penalty*h;
  return (fitness < 0 ? 0 : fitness);
}


// return a Poisson random variate with given mean
int poisson(double expectedValue)
{
  int n = 0; //counter of iteration
  double limit; 
  double x;  //pseudo random number
  limit = exp(-expectedValue);
  x = RND;
  while (x > limit) {
    n++;
    x *= RND;
  }
  return n;
}

int binomial(int n, double p)
{
  int i;
  int k = 0;
  
  for(i = 0; i < n; i++)
    {
      if(RND < p) k++;
    }
  return k;
}

int main(int argc, char *argv[])
{
  Ind I[NPOP], newI[NPOP];
  double f[NPOP], cs[NPOP];
  int i, j;
  int t;
  int mum, dad;
  int DUI;
  double LEAKAGE;
  double MU;
  int NDNA;
  FILE *fp, *fpm;
  int env = 0;
  double ball;
  int expt;
  double amu, bmu;
  int change;
  double pa, pb, r;
  double meanf, varf;
  char fstr[100];
  double scale, penalty;
  
  if(argc != 4) {
    printf("Please specify environmental change period, fitness scale, and heteroplasmy penalty\n");
    exit(0);
  }
  env = atoi(argv[1]);
  scale = atof(argv[2]);
  penalty = atof(argv[3]);
  
  // open file for output
#ifdef _FULLOUTPUT
  sprintf(fstr, "inherit-zoom-full-out-%i-%.3f-%.3f.csv", env, scale, penalty);
  fp = fopen(fstr, "w");
  fprintf(fp, "scale,penalty,env,nDNA,mu,DUI,leakage,expt,t,i,a,b,c,f\n");
#endif

#ifdef _MEANOUTPUT
  sprintf(fstr, "inherit-zoom-mean-out-%i-%.3f-%.3f.csv", env, scale, penalty);
  fpm = fopen(fstr, "w");
  fprintf(fpm, "scale,penalty,env,nDNA,mu,DUI,leakage,expt,t,mean.f,var.f\n");
#endif
  
  // loop over different environment types
  //  for(env = 0; env <= 6; env++)
    {
      // loop over DNA population size
      for(NDNA = 10; NDNA < 1000; NDNA *= 1.5)
	{
	  // loop over different mutation rates
	  // for(MU = 0; MU <= 0.01; MU *= 10)
	  MU = 0;
	    {
	      // loop over doubly-uniparental inheritance
	      DUI = 0;
	      // for(DUI = 0; DUI <= 1; DUI++)
		{
		  // loop over paternal leakage
		  for(LEAKAGE = 0; LEAKAGE <= 0.5; LEAKAGE *= 1.2)
		    {
		      // loop over instances for each parameterisation
		      for(expt = 0; expt < NSAMP; expt++)
			{
			  // initialise population. every even-index individual has NDNA a; every odd-index individual has NDNA b
			  // we picture the first half of the population as female and the second half as male
			  for(i = 0; i < NPOP; i++)
			    { I[i].a = (i%2 == 0 ? NDNA : 0); I[i].b = NDNA-I[i].a; I[i].c = 0; }

			  printf("%i %i %f %i %f %i\n", env, NDNA, MU, DUI, LEAKAGE, expt);
		      
			  // loop through generations
			  for(t = 0; t < NGEN; t++)
			    {
			      //// assign fitnesses to existing population

			      // assign fitnesses to first half of the population (mothers) and populate cumulative fitness sum
			      for(i = 0; i < (NPOP/2); i++)
				{
				  f[i] = fitness(I[i], env, t, scale, penalty);
				  cs[i] = (i == 0 ? f[i] : f[i]+cs[i-1]);
				}
			      // normalise
			      for(i = 0; i < (NPOP/2); i++)
				cs[i] /= cs[(NPOP/2)-1];

			      // assign fitnesses to second half of the population (fathers) and populate cumulative fitness sum
			      for(i = (NPOP/2); i < NPOP; i++)
				{
				  f[i] = fitness(I[i], env, t, scale, penalty);
				  cs[i] = (i == (NPOP/2) ? f[i] : f[i]+cs[i-1]);
				}
			      // normalise
			      for(i = (NPOP/2); i < NPOP; i++)
				cs[i] /= cs[NPOP-1];

#ifdef _MEANOUTPUT
			      if((expt < 4 && t % 50 == 0) || t == 1995)
				{
				  meanf = varf = 0;
				  for(i = 0; i < NPOP; i++)
				    meanf += f[i];
				  meanf /= NPOP;
				  for(i = 0; i < NPOP; i++)
				    varf += (f[i]-meanf)*(f[i]-meanf);
				  varf /= (NPOP-1);
				  fprintf(fpm, "%f,%f,%i,%i,%f,%i,%f,%i,%i,%f,%f\n", scale, penalty, env, NDNA, MU, DUI, LEAKAGE, expt, t, meanf, varf);
				}
#endif
			      
#ifdef _FULLOUTPUT
			      // output population
			      for(i = 0; i < NPOP; i++)
				{
				  if(t == 10 || t == 100 || t == 500 || t == 900)
				    {
				      fprintf(fp, "%f,%f,%i,%i,%f,%i,%f,%i,%i,%i,%i,%i,%i,%f\n", scale,penalty,env, NDNA, MU, DUI, LEAKAGE, expt, t, i, I[i].a, I[i].b, I[i].c, f[i]);
				    }
				}
#endif

			      //// build up new population

			      // loop over new individuals
			      for(i = 0; i < NPOP; i++)
				{
				  // select a mother from the first half of the population
				  ball = RND;
				  for(j = 0; cs[j] < ball; j++);
				  mum = j;
				  // select a father from the second half of the population
				  ball = RND;
				  for(j = NPOP/2; cs[j] < ball; j++);
				  dad = j;

				  // if we're allowing DUI, with prob 0.5, treat the father as the "mother" (hence inherit from father)
				  if(DUI && RND < 0.5)
				    {
				      j = mum; mum = dad; dad = j;
				    }

				  // construct new individual's genetic makeup from binomial draws with mean (1-leak)*mother + leak*father
				  newI[i].a = binomial(I[dad].a, LEAKAGE) + binomial(I[mum].a, 1.-LEAKAGE);
				  newI[i].b = binomial(I[dad].b, LEAKAGE) + binomial(I[mum].b, 1.-LEAKAGE);
				  newI[i].c = binomial(I[dad].c, LEAKAGE) + binomial(I[mum].c, 1.-LEAKAGE);

				  // apply mutations from binomial draws with mean a*MU, b*MU
				  if(MU > 0)
				    {
				      amu = binomial(newI[i].a, MU);
				      bmu = binomial(newI[i].b, MU);
				      newI[i].a -= amu; newI[i].c += amu;
				      newI[i].b -= bmu; newI[i].c += bmu;
				    }

				  // bring total cellular population to constant size (from above or below)
				  if(newI[i].a+newI[i].b+newI[i].c > NDNA) change = -1;
				  else change = 1;
			      
				  while(newI[i].a+newI[i].b+newI[i].c != NDNA)
				    {
				      pa = (double)newI[i].a/(newI[i].a+newI[i].b+newI[i].c);
				      pb = (double)newI[i].b/(newI[i].a+newI[i].b+newI[i].c);
				      r = RND;
				      if(r < pa) newI[i].a += change;
				      else if(r < pa+pb) newI[i].b += change;
				      else newI[i].c += change;
				    }
				}
			  
			      // update population
			      for(i = 0; i < NPOP; i++)
				{
				  I[i] = newI[i];
				}
			    }
			}
		      if(LEAKAGE == 0) LEAKAGE = 0.002;
		    }
		}
	      if(MU == 0) MU = 1e-6;
	    }
	}
    }
#ifdef _FULLOUTPUT  
  fclose(fp);
#endif
#ifdef _MEANOUTPUT
  fclose(fpm);
#endif
  
  return 0;
}
