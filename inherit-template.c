#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND drand48()    // random number between 0 and 1
#define MAXPOP 2000         // population size
#define NGEN 500        // number of generations in simulation
#define NSAMP 100        // number of samples to run

//#define _MEANOUTPUT
#define _CHANGEOUTPUT

// structure to store an individual's genetic makeup -- number of genomes of type a, b, c
typedef struct {
  int a, b, c;
} Ind;

// fitness of individual I in environment type env at time t
double fitness(Ind I, int env, int t, double scale, double penalty)
{
  double h, fitness;

  if(t > env) {
    fitness = (I.b+scale*I.a);
  } else {
    fitness = (I.a+scale*I.b);
  }
  
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
  Ind I[MAXPOP], newI[MAXPOP];
  double f[MAXPOP], cs[MAXPOP];
  int NPOP;
  int i, j;
  int t;
  int mum, dad;
  int DUI;
  double LEAKAGE;
  double MU;
  int NDNA;
  FILE *fp, *fpm, *fpc;
  int env = 0;
  double ball;
  int expt;
  double amu, bmu;
  int change;
  double pa, pb, r, r1, r2;
  int da, db, dc, type1, type2;
  double meanf, varf, meanh;
  char fstr[100];
  double scale, penalty;
  int reporttime;
  double meanmean, meanmeanh, meannorm;
  double meanw, meanwnorm;
  int ICs, DET_REAMP, DET_LEAK;
  double TEMPLATE, FITNESSB, FITNESSC;
  int CLUSTER;
  int arep, brep, crep;
  
  if(argc != 12) {
    printf("Please specify Npop, initial conditions, environmental change period, fitness scale, heteroplasmy penalty, deterministic reamp, deterministic leakage, templating rate, rel int fitness of B allele, rel int fitness of C allele, cluster size\n");
    exit(0);
  }
  NPOP = atoi(argv[1]);
  ICs = atoi(argv[2]);
  env = atoi(argv[3]);
  scale = atof(argv[4]);
  penalty = atof(argv[5]);
  DET_REAMP = atoi(argv[6]);
  DET_LEAK = atoi(argv[7]);
  TEMPLATE = atof(argv[8]);
  FITNESSB = atof(argv[9]);
  FITNESSC = atof(argv[10]);
  CLUSTER = atoi(argv[11]);
  
  // open file for output
#ifdef _FULLOUTPUT
  sprintf(fstr, "inherit-template-full-out-%i-%i-%i-%.3f-%.3f-%i-%i-%.5f-%.2f-%.2f-%i.csv", NPOP, ICs, env, scale, penalty, DET_REAMP, DET_LEAK, TEMPLATE, FITNESSB, FITNESSC, CLUSTER);
  fp = fopen(fstr, "w");
  fprintf(fp, "Npop,ICs,scale,penalty,det.reamp,det.leak,template,fitness.b,fitness.c,env,nDNA,mu,DUI,leakage,cluster,expt,t,i,a,b,c,f\n");
#endif

#ifdef _MEANOUTPUT
  sprintf(fstr, "inherit-template-mean-out-%i-%i-%i-%.3f-%.3f-%i-%i-%.5f-%.2f-%.2f-%i.csv", NPOP, ICs, env, scale, penalty, DET_REAMP, DET_LEAK, TEMPLATE, FITNESSB, FITNESSC, CLUSTER);
  fpm = fopen(fstr, "w");
  fprintf(fpm, "Npop,ICs,scale,penalty,det.reamp,det.leak,template,fitness.b,fitness.c,env,nDNA,mu,DUI,leakage,cluster,expt,t,mean.f,var.f,mean.h\n");
#endif

#ifdef _CHANGEOUTPUT
  sprintf(fstr, "inherit-template-change-out-%i-%i-%i-%.3f-%.3f-%i-%i-%.5f-%.2f-%.2f-%i.csv", NPOP, ICs, env, scale, penalty, DET_REAMP, DET_LEAK, TEMPLATE, FITNESSB, FITNESSC, CLUSTER);
  fpc = fopen(fstr, "w");
  fprintf(fpc, "Npop,ICs,scale,penalty,det.reamp,det.leak,template,fitness.b,fitness.c,env,nDNA,mu,DUI,leakage,cluster,expt,end.mean.f,end.mean.h,window.mean.f\n");
#endif

  // loop over different environment types
  //  for(env = 0; env <= 6; env++)
  {
    // loop over DNA population size
    for(NDNA = 10; NDNA <= 1200; NDNA *= 2)
      {
	// loop over different mutation rates
	for(MU = 0; MU <= 1e-2; MU *= 10)
	  {
	    // loop over doubly-uniparental inheritance
	    for(DUI = 0; DUI <= 1; DUI++)
	    {
	      // loop over paternal leakage
	      for(LEAKAGE = 0; LEAKAGE <= 0.52; LEAKAGE *= 2)
		{
		  // loop over instances for each parameterisation
		    for(expt = 0; expt < NSAMP; expt++)
		      {
			// initialise population. every even-index individual has NDNA a; every odd-index individual has NDNA b
			// we picture the first half of the population as female and the second half as male
			for(i = 0; i < NPOP; i++)
			  {
			    if(ICs == 0)
			      {
				I[i].a = NDNA/2.; I[i].b = NDNA/2.; I[i].c = 0;
			      }
			    else
			      {
				I[i].a = (i%2 == 0 ? NDNA : 0); I[i].b = NDNA-I[i].a; I[i].c = 0;
			      }
			  }
			
			printf("%i %i %f %i %f %i\n", env, NDNA, MU, DUI, LEAKAGE, expt);

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
#ifdef _MEANOUTPUT
			    if(expt < 4)
			      fprintf(fpm, "%i,%i,%f,%f,%i,%i,%f,%f,%i,%i,%f,%i,%f,%i,%i,%f,%f,%f\n", NPOP, ICs, scale, penalty, DET_REAMP, DET_LEAK, TEMPLATE, FITNESS, env, NDNA, MU, DUI, LEAKAGE, expt, t, meanf, varf, meanh);
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
				    fprintf(fp, "%i,%i,%f,%f,%i,%i,%f,%f,%i,%i,%f,%i,%f,%i,%i,%i,%i,%i,%i,%i,%f\n", NPOP, ICs, scale,penalty, DET_REAMP, DET_LEAK, TEMPLATE, FITNESS, env, NDNA, MU, DUI, LEAKAGE, CLUSTER, expt, t, i, I[i].a, I[i].b, I[i].c, f[i]);
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

				// construct new individual's genetic makeup from binomial draws with mean (1-leak)*mother + leak*father
				do{
				  if(DET_LEAK) {
				    newI[i].a = CLUSTER*round(LEAKAGE*binomial(round((double)I[dad].a/CLUSTER), 0.5) + CLUSTER*(1.-LEAKAGE)*binomial(round((double)I[mum].a/CLUSTER), 0.5));
				    newI[i].a = CLUSTER*round(LEAKAGE*binomial(round((double)I[dad].b/CLUSTER), 0.5) + CLUSTER*(1.-LEAKAGE)*binomial(round((double)I[mum].b/CLUSTER), 0.5));
				    newI[i].a = CLUSTER*round(LEAKAGE*binomial(round((double)I[dad].c/CLUSTER), 0.5) + CLUSTER*(1.-LEAKAGE)*binomial(round((double)I[mum].c/CLUSTER), 0.5));
				  } else {
				    newI[i].a = CLUSTER*binomial(round((double)I[dad].a/CLUSTER), LEAKAGE*0.5) + CLUSTER*binomial(round((double)I[mum].a/CLUSTER), (1.-LEAKAGE)*0.5);
				    newI[i].b = CLUSTER*binomial(round((double)I[dad].b/CLUSTER), LEAKAGE*0.5) + CLUSTER*binomial(round((double)I[mum].b/CLUSTER), (1.-LEAKAGE)*0.5);
				    newI[i].c = CLUSTER*binomial(round((double)I[dad].c/CLUSTER), LEAKAGE*0.5) + CLUSTER*binomial(round((double)I[mum].c/CLUSTER), (1.-LEAKAGE)*0.5);
				  }
				}
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
				ntemplates = binomial(newI[i].a+newI[i].b+newI[i].c, TEMPLATE);
				for(i = 0; i < ntemplates; i++)
				  {
				    da = db = dc = 0;
				    pa = (double)newI[i].a / (newI[i].a + newI[i].b + newI[i].c);
				    pb = (double)newI[i].b / (newI[i].a + newI[i].b + newI[i].c);
				    r1 = RND; r2 = RND;
				    if(r1 < pa) type1 = 1;
				    else if(r2 < pa+pb) type1 = 2;
				    else type1 = 3;
				    if(r2 < pa) type2 = 1;
				    else if(r2 < pa+pb) type2 = 2;
				    else type2 = 3;
				    if(type1 == 1) { if(type2 == 2) { da++; db--; } if(type2 == 3) { da++; dc--; } }
				    if(type1 == 2) { if(type2 == 1) { db++; da--; } if(type2 == 3) { db++; dc--; } }
				    if(type1 == 3) { if(type2 == 1) { dc++; da--; } if(type2 == 2) { dc++; db--; } }
				    newI[i].a += da; newI[i].b += db; newI[i].c += dc;
				  }
			      }
				
			    // deterministic reamplification
				if(DET_REAMP)
				  {
				    newI[i].a *= 2; newI[i].b *= 2; newI[i].c *= 2;
				    int newsum = newI[i].a+newI[i].b+newI[i].c;
				    newI[i].a = round((double)newI[i].a*NDNA/newsum);
				    newI[i].b = round((double)newI[i].b*NDNA/newsum);
				    newI[i].c = round((double)newI[i].c*NDNA/newsum);
				  }
			      
				// bring total cellular population to constant size (from above or below)
				if(newI[i].a+newI[i].b+newI[i].c > NDNA) change = -1;
				else change = 1;
			      
				while(newI[i].a+newI[i].b+newI[i].c != NDNA)
				  {
				    pa = (double)newI[i].a          / (newI[i].a + FITNESSB*newI[i].b + FITNESSC*newI[i].c);
				    pb = (double)newI[i].b*FITNESSB / (newI[i].a + FITNESSB*newI[i].b + FITNESSC*newI[i].c);
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
#ifdef _CHANGEOUTPUT
		    fprintf(fpc, "%i,%i,%f,%f,%i,%i,%f,%f,%f,%i,%i,%f,%i,%f,%i,%i,%f,%f,%f\n", NPOP, ICs, scale,penalty,DET_REAMP, DET_LEAK, TEMPLATE, FITNESSB, FITNESSC, env, NDNA, MU, DUI, LEAKAGE, CLUSTER, expt,  meanmean/meannorm, meanmeanh/meannorm, meanw/meanwnorm);
#endif
		    
		      }
		    if(LEAKAGE == 0) LEAKAGE = 5e-4;
		  }
	      }
	    if(MU == 0) MU = 1e-7;
	  }
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
