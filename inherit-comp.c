#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND drand48()    // random number between 0 and 1
#define MAXPOP 2000         // population size
#define NGEN 500        // number of generations in simulation
#define NSAMP 10       // number of samples to run

#define _MEANOUTPUT
//#define _CHANGEOUTPUT

// structure to store an individual's genetic makeup -- number of genomes of type a, b, c
typedef struct {
  int a, b, c;
  double leakage;
  int NDNA;
} Ind;

// fitness of individual I in environment type env at time t
double fitness(Ind I, int env, int t, double scale, double penalty)
{
  double h, fitness;
  double ap, bp;

  ap = (double)I.a/I.NDNA;
  bp = (double)I.b/I.NDNA;
  
  if((t/env)%2 == 1) {
    fitness = (bp+scale*ap);
  } else {
    fitness = (ap+scale*bp);
  }
  
  h = (bp < ap ? bp : ap);

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

// produce gaussian random number
double gsl_ran_gaussian(const double sigma)
{
  double x, y, r2;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      x = -1 + 2 * RND;
      y = -1 + 2 * RND;

      /* see if it is in the unit circle */
      r2 = x * x + y * y;
    }
  while (r2 > 1.0 || r2 == 0);

  /* Box-Muller transform */
  return sigma * y * sqrt (-2.0 * log (r2) / r2);
}


int main(int argc, char *argv[])
{
  Ind I[MAXPOP], newI[MAXPOP];
  double f[MAXPOP], cs[MAXPOP];
  int NPOP;
  int i, j;
  int t;
  int mum, dad;
  int DUI;
  double MU;
  FILE *fp, *fpm, *fpc;
  int env = 0;
  double ball;
  int expt;
  double amu, bmu;
  int change;
  double pa, pb, r;
  double meanf, varf, meanh;
  char fstr[100];
  double scale, penalty;
  int reporttime;
  double meanmean, meanmeanh, meannorm;
  double meanw, meanwnorm;
  double TEMPLATE, FITDIFF;
  int arep, brep, crep;
  double leak, leak2, ndna, ndna2;
  
  if(argc != 6) {
    printf("Please specify Npop, fitness scale, heteroplasmy penalty, templating rate, rel int fitness of B allele\n");
    exit(0);
  }
  NPOP = atoi(argv[1]);
  scale = atof(argv[2]);
  penalty = atof(argv[3]);
  TEMPLATE = atof(argv[4]);
  FITDIFF = atof(argv[5]);
  
  // open file for output
#ifdef _FULLOUTPUT
  sprintf(fstr, "inherit-template-full-out-%i-%i-%i-%.3f-%.3f-%i-%i-%.5f-%.2f.csv", NPOP, ICs, env, scale, penalty, DET_REAMP, DET_LEAK, TEMPLATE, FITDIFF);
  fp = fopen(fstr, "w");
  fprintf(fp, "Npop,ICs,scale,penalty,det.reamp,det.leak,template,fitdiff,env,nDNA,mu,DUI,leakage,expt,t,i,a,b,c,f\n");
#endif

#ifdef _MEANOUTPUT
  sprintf(fstr, "inherit-comp-%i-%.3f-%.3f-%.5f-%.3f.csv", NPOP, scale, penalty, TEMPLATE, FITDIFF);
  fpm = fopen(fstr, "w");
  fprintf(fpm, "env,mu,Npop,scale,penalty,template,fitdiff,expt,t,curr.env,mean.f,var.f,mean.h,mean.leak,var.leak,mean.NDNA,var.NDNA\n");
#endif

#ifdef _CHANGEOUTPUT
  sprintf(fstr, "inherit-template-change-out-%i-%i-%i-%.3f-%.3f-%i-%i-%.5f-%.2f.csv", NPOP, ICs, env, scale, penalty, DET_REAMP, DET_LEAK, TEMPLATE, FITDIFF);
  fpc = fopen(fstr, "w");
  fprintf(fpc, "Npop,ICs,scale,penalty,det.reamp,det.leak,template,fitdiff,env,nDNA,mu,DUI,leakage,expt,end.mean.f,end.mean.h,window.mean.f\n");
#endif

  DUI = 0;
  for(env = 1; env <= 32; env *= 2)
    {
      for(MU = 0; MU <= 1e-2; MU *= 10)
	{
          printf("%i %f\n", env, MU);

	  // loop over instances for each parameterisation
	  for(expt = 0; expt < NSAMP; expt++)
	    {
	      // initialise population. every even-index individual has NDNA a; every odd-index individual has NDNA b
	      // we picture the first half of the population as female and the second half as male
	      for(i = 0; i < NPOP; i++)
		{
		  I[i].leakage = RND*0.5; I[i].NDNA = pow(10, RND*3);
		  if(I[i].NDNA < 5) I[i].NDNA = 5;
		  I[i].a = (i%2 == 0 ? I[i].NDNA : 0); I[i].b = I[i].NDNA-I[i].a; I[i].c = 0;
		}
			
	      meanmean = meannorm = meanmeanh = 0;
	      meanw = meanwnorm = 0;
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

		  meanf = varf = meanh = 0;
		  for(i = 0; i < NPOP; i++)
		    {
		      meanf += f[i];
		      meanh += (I[i].a < I[i].b ? I[i].a : I[i].b);
		    }
		  meanf /= NPOP;
		  meanh /= NPOP;
		  for(i = 0; i < NPOP; i++)
		    varf += (f[i]-meanf)*(f[i]-meanf);
		  varf /= (NPOP-1);

		  // summary stats of strategies
		  leak = leak2 = 0;
		  ndna = ndna2 = 0;
		  for(i = 0; i < NPOP; i++)
		    {
		      leak += newI[i].leakage;
		      leak2 += newI[i].leakage*newI[i].leakage;
		      ndna += newI[i].NDNA;
		      ndna2 += newI[i].NDNA*newI[i].NDNA;
		    }
		  leak /= NPOP; leak2 /= NPOP; ndna /= NPOP; ndna2 /= NPOP;
			    
#ifdef _MEANOUTPUT
		  fprintf(fpm, "%i,%f,   %i,%f,%f,  %f,%f,  %i,%i,%i,  %f,%f,%f,  %f,%f,%f,%f\n",
			  env, MU,
			  NPOP, scale, penalty,
			  TEMPLATE, FITDIFF,
			  expt, t, (t/env)%2, 
			  meanf, varf, meanh,
			  leak, leak2-leak*leak, ndna, ndna2-ndna*ndna);
#endif

		  // sample mean fitness after adaptive final section
		  if(t == NGEN-1) {
		    meanmean += meanf;
		    meannorm ++;
		    meanmeanh += meanh;
		    //	  }
		    // sample mean fitness over window of final env period before adaptive final section
		    //if((env==0 && t==NGEN-1) || (t>NGEN-1 - 2*env)) {
		    meanw += meanf;
		    meanwnorm ++;
		  }
			    
#ifdef _FULLOUTPUT
		  // output population
		  for(i = 0; i < NPOP; i++)
		    {
		      if(t == 10 || t == 100 || t == 500 || t == 900)
			{
			  fprintf(fp, "%i,%i,%f,%f,%i,%i,%f,%f,%i,%i,%f,%i,%f,%i,%i,%i,%i,%i,%i,%f\n", NPOP, ICs, scale,penalty, DET_REAMP, DET_LEAK, TEMPLATE, FITDIFF, env, NDNA, MU, DUI, LEAKAGE, expt, t, i, I[i].a, I[i].b, I[i].c, f[i]);
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

		      // if we're allowing DUI, for the males of the new population, treat the father as the "mother" (hence inherit from father)
		      if(DUI && i > NPOP/2)
			{
			  j = mum; mum = dad; dad = j;
			}

		      // construct new individual's strategy from mum and dad
		      newI[i].leakage = (I[mum].leakage + I[dad].leakage)/2. + gsl_ran_gaussian(0.02);
		      if(newI[i].leakage < 0) newI[i].leakage = 0;
		      if(newI[i].leakage > 0.5) newI[i].leakage = 0.5;
		      newI[i].NDNA = pow(10, log10(I[mum].NDNA) + log10(I[dad].NDNA) + gsl_ran_gaussian(0.05));
		      if(newI[i].NDNA < 5) newI[i].NDNA = 5;
		      if(newI[i].NDNA > 1000) newI[i].NDNA = 1000;
				
		      // construct new individual's genetic makeup from binomial draws with mean (1-leak)*mother + leak*father
		      do{
			newI[i].a = binomial(I[dad].a, newI[i].leakage*0.5) + binomial(I[mum].a, (1.-newI[i].leakage)*0.5);
			newI[i].b = binomial(I[dad].b, newI[i].leakage*0.5) + binomial(I[mum].b, (1.-newI[i].leakage)*0.5);
			newI[i].c = binomial(I[dad].c, newI[i].leakage*0.5) + binomial(I[mum].c, (1.-newI[i].leakage)*0.5);
		      }while(newI[i].a + newI[i].b + newI[i].c == 0);
				
		      // apply mutations from binomial draws with mean a*MU, b*MU
		      if(MU > 0)
			{
			  amu = binomial(newI[i].a, MU);
			  bmu = binomial(newI[i].b, MU);
			  newI[i].a -= amu; newI[i].c += amu;
			  newI[i].b -= bmu; newI[i].c += bmu;
			}

		      // templated repair
		      if(TEMPLATE > 0)
			{
			  // decide how many mutants get template-repaired overall
			  crep = binomial(newI[i].c, TEMPLATE*(newI[i].a+newI[i].b));
			  // decide how many were templated from type a
			  arep = binomial(crep, (double)(newI[i].a)/(newI[i].a+newI[i].b));
			  brep = crep-arep;
			  newI[i].a += arep; newI[i].b += brep;
			  newI[i].c -= crep;
			}
				
		      // bring total cellular population to constant size (from above or below)
		      if(newI[i].a+newI[i].b+newI[i].c > newI[i].NDNA) change = -1;
		      else change = 1;
			      
		      while(newI[i].a+newI[i].b+newI[i].c != newI[i].NDNA)
			{
			  pa = (double)newI[i].a/(newI[i].a+newI[i].b+newI[i].c);
			  pb = (double)newI[i].b*FITDIFF/(newI[i].a+newI[i].b+newI[i].c);
			  r = RND;
			  if(r < pa) newI[i].a += change;
			  else if(r < pa+pb*FITDIFF) newI[i].b += change;
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
	  if(MU == 0) MU = 1e-7;
	}
    }
#ifdef _FULLOUTPUT  
  fclose(fp);
#endif
#ifdef _MEANOUTPUT
  fclose(fpm);
#endif
#ifdef _CHANGEOUTPUT
  fclose(fpc);
#endif
  
  return 0;
}
