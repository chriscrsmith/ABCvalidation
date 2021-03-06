#include "data_sumstat.h"
#include "assert.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define MISSING '2'

void set_ders(int** ders, char** list, int n, int segsites, int start, int end){
  assert(end <= n);
  int i,j,k;
  for(i=0; i<segsites; i++){
    
    k=0;
    for(j=start; j<end; j++){
      //printf("k:%d\n", k); //#p
      if(list[j][i] == '1'){
	ders[k][i] = j;
	k++;
      }
    }
  }
}
      
    


void set_maxder(int* maxder, char** list, int n, int segsites, int start, int end){
  if(end > n){
    fprintf(stderr, "end:%d, n:%d\n", end, n);
    assert(end <= n);
  }
  
  int i, j, max;
  for(i=0; i<segsites; i++){
    max=start;
    for(j=start; j<end; j++)
      if(list[j][i] == '1')  max=j;
    maxder[i] = max;
  }
}
    

/*set the allelic map */
void set_almap(int* almap, char** list, int n, int segsites, int start, int end){
  //printf("start: %d, end: %d, n: %d\n", start, end, n);
  if(end > n){
    fprintf(stderr, "end:%d, n:%d\n", end, n);
    assert(end <= n);
  }
  //printf("\nenn and n is %d and %d\n", end, n);
  
  int i,j;
  for(i=0; i<segsites; i++){
    int ones = 0;
    
    for( j=start; j<end; j++){
      ones += ( list[j][i] == '1' );
    }
    
    almap[i] = ones;
    //printf("i:%i\tsegsites:%i\talmap[]:%i\tstart:%d\tend:%d\n", i,segsites, almap[i], start, end);
  }
}


void set_missing(int* missing, char** list, int n, int segsites, int start, int end){
  if(end > n){
    fprintf(stderr, "end:%d, n:%d\n", end, n);
    assert(end <= n);
  }

  int i,j,mis=0;
  for(i=0; i<segsites; i++){
    mis=0;
    
    for( j=start; j<end; j++)
      mis += ( list[j][i] == MISSING );

    missing[i] = mis;
    //fprintf(stderr, "i:%i\tsegsites:%i\tmissing[]:%i\tstart:%d\tend:%d\n", i,segsites, missing[i], start, end);
  }
}


// number of segregating sites is based only on the '1's
int seg_sites(int* almap, int n, int segsites){
  int i;
  int s = 0;

  for(i=0; i<segsites; i++){
    if((almap[i] > 0) && (almap[i] < n))
      s++;
  }
  return s;
}

/* this is the Thomson estimator as it is described in Hudson et al 2007. 
   Notice that I have adapted that to missing data. Thus, for every segregating sites I divide by
   the size of the site (excluding missing data).
*/
int thomsonEst(double* thomson, int* almap, int n, int segsites, int* missing){
  
  double te = 0.;
  int i,nsam;
  assert(n>1);
  
  for(i=0; i<segsites; ++i){
    nsam = n - missing[i];
    
    if(nsam <= 0) return 0;
    
    if(almap[i] == 0 || almap[i] == nsam) continue;

    te += ((double)almap[i]/(double)nsam);
  }
  *thomson = te;

  return 1;
}

/* this is the estimate of variance of Thomson estimator, 
   see Hudson 2007: The variance of coalescent times estimates
*/
int thomsonVar(double* thomson, int* almap, int n, int segsites){
  double tv = 0.;
  int i,ones;
  assert(n > 1);
  
  int* sfs = malloc((n-1)*sizeof(int));
  for(i=0; i<n-1; ++i)
    sfs[i] = 0;

  for(i=0; i<segsites; ++i){
    
    if(almap[i] == 0 ||almap[i] == n) continue;

    ones = almap[i];
    if(ones > 0 && ones < n)
      sfs[ones-1] ++;
  }
  
  for(i=1; i<n; ++i)
    tv += i*i*sfs[i-1];
  tv /= ((double)n*n);
  *thomson = tv;
  free(sfs);
  return 1;
}

/*
  calculates the theta pi of simulated data stored in list
*/
int theta_pi(double* theta, int* almap, int n, int segsites, int* missing){
  double pi = 0.;
  int i,j, ones,  nsam;
  
  assert(n > 1);
  double denom;

  //denom = n*(n - 1);
  for(i=0; i<segsites; i++){
    nsam = n - missing[i];
    
    denom = nsam*(nsam - 1.);
    
    if(denom <= 0) continue;
      
    double ssh = 0.; // sum of site homozygosity;
    ones = almap[i];
    
    /* in cases of non-polymorphic sites add 0 to pi and continue the loop */
    if( ones == nsam || ones == 0)
      continue;
    /* for the 1's class of alleles */
    ssh += ((double)ones * (double)(ones - 1))/denom;
    /* for the 0's class of alleles */
    ssh += ((double)(nsam - ones) * (double)( nsam -ones - 1.))/ denom;
    
    pi += (1.0 - ssh);
  }
  *theta = pi;
  return 1;
}

// these denominators get the maximum n for their calculation. 
int denominators(int n, double* hn, double* sqhn, double* bn){
  int i;
  *hn=0.;
  *sqhn = 0.;
  for(i=1; i<n; i++){
    *hn += 1./(double)i;
    *sqhn += 1./(double)(i*i);
  }
  *bn = *sqhn + 1./((double)n*n);
  return 1;
}


  int theta_w(double* theta, int* almap, int n, int segsites, double denom, int* missing){
  double w = 0.;
  int i,j;
  int true_segsites = 0, nsam;
  
  for(i=0; i<segsites; i++){
    nsam = n - missing[i];
    if(nsam <= 0) continue;
    /* in cases of non-polymorphic sites add 0 to pi and continue the loop */
    if( almap[i] == nsam || almap[i] == 0)
      continue;
    
    /* increase the true segsites by one; */
    true_segsites++;
  }
  
  w = (true_segsites)/denom;
  *theta = w;
  //printf("*w:%e\n", w);
  return 1;
}

// this denominator works with the original sample size. i.e. does not know about missing data
double Dnominator(int n, int segsites, double hn, double sqhn){
  double epsilon = 1e-10;
  double a1 = 0., a2 = 0., b1 = 0., b2=0., c1=0., c2=0., e1=0., e2=0., denom=0.;
  int i;
  
  a1 = hn;
  a2 = sqhn;
  /* for(i=1; i<n; i++){ */
  /*   a1 += 1./(double)i; */
  /*   a2 += 1./pow( (double)i, 2.0 ); */
  /* } */
  b1 = ((double)n+1.)/(3.*((double)n-1.));
  b2 = (2.*((double)n*(double)n + (double)n + 3.))/(9.*(double)n*((double)n-1.));
  c1 = b1 - 1./a1;
  
  
  if(fabs(c1) < epsilon) c1 = 0.;
 
  c2 = b2 - (n+2)/(a1*n) + (a2/(a1*a1));
  if(fabs(c2) < epsilon) c2 = 0.;

  e1 = c1/a1;
  if(fabs(e1) < epsilon) e1 = 0.;

  e2 = c2/(a1*a1 + a2);
  if(fabs(e2) < epsilon) e2 = 0.;

    
  denom=(e1*segsites + e2*segsites*((double)segsites - 1.) );
  if(fabs(denom) < epsilon) denom = epsilon;

  if(denom <= 0){
    printf("sample size: %d\n", n);
    printf("segsites: %d\ne1: %e\ne2: %e\n", segsites, e1, e2);
    fprintf(stderr,
	    "a1: %e\na2: %e\nb1: %e\nb2: %e\nc1: %e\nc2: %e\ne1: %e\ne2: %e\n",
	    a1, a2, b1, b2, c1, c2, e1, e2);
  }
  assert(denom > 0);

  denom = sqrt(denom);
  return denom;
}

// missing data affect tajD through the numerator
int tajD(double* tajd, int segsites, int n, double thetaw, double thetap, double hn, double sqhn){
  if(n <= 0)
    return 0;
  
  double denom, td;
  denom = Dnominator(n, segsites, hn, sqhn);
  assert(denom > 0.);
  td = (thetap - thetaw)/denom;
  *tajd = td;
  return 1;
}



// missing data affect ZnS biasing it to higher values
int ZnS(double* ld, char** list, int n, int segsites, int filter, int* almap, int** ders, int* missing){
  
  double zns = 0.;
  int class = 0;
  int exc;
  
  //exc = (int*)malloc(segsites*sizeof(int));
  
  int i =0, j=0;
  int counter = 0;
  
  for(i=0; i<segsites-1; i++){
    //printf("i:%d\n", i);
    exc = ( (almap[i] == filter) || (almap[i] == n-filter) ) ? 1 : 0;
    
    if(exc) continue;
    for(j=i+1; j<segsites; j++){
      exc = ( (almap[j] == filter) || (almap[j] == n-filter) ) ? 1 : 0;

      if(!exc){
	//printf("*i,j: %d,%d - %d,%d\n", i,j,segsites -1, segsites);
	zns += r2(list, i, j, n, almap, ders, missing);
	//printf("i,j: %d,%d - %d,%d\n", i,j,segsites -1, segsites);
	counter ++;
      }
    }
  }
  
  if(counter == 0)
    return 0;
  
  zns /= counter;
  *ld = zns;
  //printf("zns %e\n", zns);
  //exit(-1);
  return 1;
}




int calc_wIbs(double* wAi,
	      double* wLi,
	      char** list,
	      int n,
	      int segsites,

	      int filter,
	      int* almap,
	      int** ders,
	      int* missing,
	      int start, 

	      int end,
	      double* segPositions,
	      int locusLength){

  int k,
    startPos,
    h1,
    h2;

  double af1,
    af2,
    currentTractLength,
    longestRun= 0.,
    numComparisons = 0.,
    totalSegments = 0,
    numIbsSegments;

  //char *currentSeq1, *currentSeq2;
  char *currentSeq1 = (char*)malloc((segsites+1)*(sizeof(char)));
  char *currentSeq2 = (char*)malloc((segsites+1)*(sizeof(char)));

  
  // find shared IBS tracts. NOTE: this approach- comparing every haplotype- is not ideal for short sequences
  if (segsites > 0){ // compare each pair of haplotypes
    for(h1=start; h1<(end-1); h1++){ // for each haplotype (in pop1)
      strncpy(currentSeq1, list[h1], segsites+1);
      for(h2=(h1+1); h2<end; h2++){ // compare with each haplotype (in pop2)
	strncpy(currentSeq2, list[h2], segsites+1);
	numComparisons = numComparisons + 1;
	numIbsSegments = 1; // minimum number of segments = 1
	if (strcmp(currentSeq1, currentSeq2) == 0){ // same haplotype
	  longestRun = 1; // proportion 1 
	}
	else{ // analyze differences between haplotypes
	  currentTractLength = 0;
	  startPos = -1; // initialize with "-1", meaning no shared position to begin the "run" with
	  for(k=0; k<segsites; k++){
	    if(currentSeq1[k] == currentSeq2[k]){ // current position is shared
	      if(startPos == -1){ // if no current start, set to current position. Otherwise, continue on.
		startPos = k;
	      }
	    }
	    else{ // reached a position that differs
	      numIbsSegments++;
	      if (startPos == -1){ // if the previous site was also different, or this is first seg site, then find segment length
		if(k == 0){
		  currentTractLength = segPositions[0]; // this is the first seg site, so just take distance from beginning of locus
		}
		else{
		  currentTractLength = segPositions[k]-segPositions[k-1]; // otherwise, take distance since last seg site
		}
	      }
	      else if(startPos == 0){ // if current start is the first seg site, then the whole locus before the first seg site is also shared
		currentTractLength = segPositions[k]; // from the start of the sequence all the way up until the current position is shared
	      }
	      else{ // current start of run is somewhere in the middle of the locus
		currentTractLength = segPositions[k] - segPositions[startPos-1]; // "-1" because all the way back to the pre-start seg site is shared
	      }
	      if( currentTractLength > longestRun ){ // is the current run the longest run encountered so far?
		longestRun = currentTractLength; 
	      }
	      startPos = -1; // start new run
	    }
	  }
	  if(startPos != -1){ // the next several lines evaluate the end-chunk of the locus, after looping through the seg sites.
	    currentTractLength = 1 - segPositions[startPos-1]; // if final segsite shared, the whole end of the locus is shared, past the last seg site
	  }
	  else{
	    currentTractLength = 1 - segPositions[segsites-1]; // if final seg site differs, just take that position to the end of the locus
	  }
	  if( currentTractLength > longestRun ){
	    longestRun = currentTractLength;
	  }
	}
	totalSegments += numIbsSegments;
      }
    }
    totalSegments = (double)1 / ( totalSegments / (double)numComparisons ); // 1 over the average number of segments
  }
  else{ // if 0 seg sites
    longestRun = 1;
    totalSegments = 1;
  }

  if(locusLength){ // if a locfile was supplied, use the locus lengths specified, otherwise output a (0,1) value
    totalSegments = (double)totalSegments * (double)locusLength;
    longestRun = longestRun * (double)locusLength;
  }
  *wAi = totalSegments;
  *wLi = longestRun;
  free(currentSeq1);
  free(currentSeq2);
  return 1;
}





int ZnA(double* ld, char** list, int n, int segsites, int filter){}

int FuLiD(){}

int  VarPi(){}


double r2(char** list, int x1, int x2, int n, int* almap, int** ders, int* missing){
  int i=0;
  int m1 = 0, der1=0;
  int m2 = 0, der2=0;
  int  m12 = 0;
  double  cor = 0.;
  
  double m1a, m2a, m12a;
  int use=x1;
  assert(n>0);
  
  double max_missing = (missing[x1] > missing[x2]) ? (double)missing[x1] : (double)missing[x2];
  double nsam = (double)n - max_missing;
   
  if(nsam <= 0) return 0.;
  
  //fprintf(stderr, "n: %d, x1: %d, x2: %d, nsam: %e, max_missing: %e\n", n, x1, x2, nsam, max_missing);

  m1 = almap[x1];
  m2 = almap[x2];
  // choose the min of the maxs
  int  min = m1;

  if(m2 < m1){
    min=m2;
    use = x2;
  }
    
  for(i=0; i<min; i++){ // notice <
          
    //
    m12 += ((list[ ders[i][use]  ][x1] == '1') && (list[ ders[i][use] ][x2] == '1'));
    /* if(m12 > 10 ){ */
/*       printf("m1:%d, m2:%d, m12:%d\n", m1, m2, m12); */
/*       printf("*: %d, %d, %d, %c, %c, %d, %d\n", i, ders[i][use], min, list[ ders[i][use]  ][x1], list[ ders[i][use] ][x2], x1, x2); */
/*     } */
  }
  
  if(max_missing > 0){
    m1 = m2 = 0;
    for(i=0; i<almap[x1]; i++)
      m1 += (list[ ders[i][x1] ][x1] == '1' && list[ ders[i][x1] ][x2] != MISSING);
    
    for(i=0; i<almap[x2]; i++)
      m2 += (list[ ders[i][x2] ][x2] == '1' && list[ ders[i][x2] ][x1] != MISSING);
    
  }
	
  
 
  //printf("max: %d\t%d\t%d\t%d\n", min, m1, m2, m12);
  if(m1 == 0 || m1 == nsam || m2 == nsam || m2 == 0){
    return 0.;
  }
  
  m1a = (double)m1/nsam;
  m2a = (double)m2/nsam;
  m12a = (double)m12/nsam;
  
  

  cor = (m12a - m1a*m2a)*(m12a - m1a*m2a)/(m1a*(1.-m1a)*m2a*(1.-m2a));
  //printf("cor: %f\n", cor);
  /* if(use == x1) */
/*     for(i=0; i<min; i++) */
/*       fprintf(stderr, "%c,%c at %d from %d \n", list[ ders[i][use]  ][x1], list[ ders[i][use]  ][x2], ders[i][use], use); */
/*   else */
/*     for(i=0; i<min; i++) */
/*       fprintf(stderr, "%c,%c at %d from %d \n", list[ ders[i][use]  ][x2], list[ ders[i][use]  ][x1], ders[i][use], use); */
 


  if(!(cor >=-0.0001 && cor <= 1.0001)){
    for(i=0; i<min; i++)
      fprintf(stderr, "%c,%c at %d from %d \n", list[ ders[i][use]  ][x1], list[ ders[i][use]  ][x2], ders[i][use], use);

    for(i=0; i<n; i++)
      fprintf(stderr, "%c,", list[i][x1]);
    fprintf(stderr, "\n");

    for(i=0; i<n; i++)
      fprintf(stderr, "%c,", list[i][x2]);
    fprintf(stderr, "\n");

    fprintf(stderr, "corr: %e\tm1: %d\tm2: %d\tm12: %d\n", cor, m1, m2, m12);
    exit(-1);
  }
  
  return cor;
}



    
/**
   Calculations of Fst
   It needs the weights for the populations because the total frequency is needed to be calculated
   The weights can be obtained by the initial command.   

   Code adapted from libsequence
*/
int calculations(double* weights, // larger pops have greater weight
		 char** list, // the polymorphic table
		 int* config, // the configuration of the sample
		 int npop,
		 int segsites,
		 int n,
		 double* segPositions,
		 int locusLength,

		 double* piT,
		 double* piS,
		 double* piB,
		 double* piD,

		 double** shared, // the fraction of share polymorphisms
		 double** private,
		 double** fixed_dif,
		 double** dxy,
		 double** shl,
		 double** lsh,
		 int derived){ 

  //printf("\n\nFst function %f\n\n", derived);
  //printf("\ncheck this out:  n = %i, segs = %i, list[0] = %s \n", n, segsites, list[0]);
  //printf("segPositions[0], segsites: %f %i\n", segPositions[0], segsites);
  //printf("\n\nlocus length %i\n\n", locusLength);
  //printf("\n");
  
  int i, 
    j,
    c,
    k=0, 
    start=0, 
    end=0,
    ni, nj,
    segs_ij=0,
    start_i,
    end_i,
    start_j,
    end_j,
    zeros_i,
    zeros_j,
    ind,
    g,
    h1,
    h2,
    newHaplotype,
    numHaps=0,
    popIndex,
    diffs,
    truesegs,
    startPos,
    numIbsSegments = 0;

  double w_ii_sq =0.,
    denom=0.,
    pi=0.,
    pi_ij=0.,
    sum_wi_wj=0.,
    weighted_pi_ij = 0.,
    weighted_pi_ii = 0.,
    af1,
    af2,
    currentTractLength,
    longestRun;

  int* almap_i = (int*)malloc(segsites*sizeof(int));
  int* almap_j = (int*)malloc(segsites*sizeof(int));
  int* missing_i = (int*)malloc(segsites*sizeof(int));
  int* missing_j = (int*)malloc(segsites*sizeof(int));
  
  int haploCounts[npop][n];
  char **haplos = (char**)malloc(n*sizeof(char*));
  if(haplos == NULL){
    printf("\ncannot allocate memory!\n");
    exit(-1);
  }

  static int* samples_b = NULL;
  static int* samples_e = NULL;
  static int first=1;

  if(first ==1){
    samples_b = (int*)malloc(npop * sizeof(int));
    samples_e = (int*)malloc(npop * sizeof(int));
    first = 0;
  }

  for(i=0; i<npop; i++){
    for(j=0; j<npop; j++){
      shared[i][j] = 0.;
      private[i][j] = 0.;
      fixed_dif[i][j] = 0.;
      dxy[i][j] = 0.;
      shl[i][j] = 0.;
      lsh[i][j] = 0.;
    }
  }
  
  end=0;
  start=0;
  samples_b[0] = 0;
  samples_e[0] = config[0];
  for(i=1; i<npop; i++){
    samples_b[i] = samples_b[i-1]+config[i-1];
    samples_e[i] = samples_b[i]+config[i];
    //printf("\nconfig %d is %d\n", i, config[i]);
  }

  // find distinct haplotypes
  for(i=0; i < n; i++){ // allocating memory to haplotype matrix
    haplos[i] = (char*)malloc(segsites*(sizeof(char)+1));
  }
  for(i=0; i<npop; i++){ // initializing haplotypeCounts matrix with zeros for each pop
    for(j=0; j<n; j++){
      haploCounts[i][j] = 0;
    }
  }
  if (segsites > 0){
    for (ind=0; ind<n; ind++){ // will want to upgrade this to a hash table one day
      newHaplotype = 1; // default setting. 1 means the haplotype is new to our list
      if (numHaps == 0){ // first haplotype. Including this block, because below there's a line comparing to empty (in this case) haplotype array
	strncpy(haplos[0], list[ind], segsites+1);
	numHaps = 1;
	haploCounts[0][0]++; // count this haplotype for the current sample (which is population 0 here, and haplotype 0)                           
      }
      else{
	for (i=0; i<npop; i++){ // determine which population this sample belongs to
	  start = samples_b[i];
	  end = samples_e[i];
	  if (ind >= start && ind < end){
	    popIndex = i;
	  }
	}
	for (g=0; g<numHaps; g++){ // compare with haplotypes already in list
	  if (strcmp(list[ind], haplos[g]) == 0){
	    newHaplotype = 0; // 0 means the haplotype is already in the list
	    haploCounts[popIndex][g]++;
	  }
	}
	if (newHaplotype == 1){
	  strncpy(haplos[numHaps], list[ind], segsites+1);
	  haploCounts[popIndex][numHaps]++;
	  numHaps += 1;
	}
      }
    }
  }
  
  for(i=0; i<npop; i++){ // over pops
    //printf("\ni %d, start %d, end %d\n", i, start, end); 
  
    start = samples_b[i];
    end = samples_e[i];
    //printf("\n%d, %d, %d, %d, %e\n", i, segsites, start, end, weights[i]);
    set_almap(almap_i, list, n, segsites, start, end);
    set_missing(missing_i, list, n, segsites, start, end);
    //printf("\n missing %i\n", missing_i[1]);
    
    pi=0.;
    w_ii_sq += weights[i]*weights[i];
    for(j=0; j<segsites; j++){ // over pol sites
      
      // set the sample size
      ni=config[i] - missing_i[j];
      if(ni < 2) continue; /* pi for a population of a sample size 0 or 1 is 0 */
      
      double ssh=0.;
      denom = (double)ni * ((double)ni-1.);

      ssh += (ni-almap_i[j] > 0) ? (double)(ni-almap_i[j])*(double)(ni-almap_i[j] - 1)/denom : 0.;
      ssh += (almap_i[j] > 0) ? (double)almap_i[j] * (double)(almap_i[j] - 1)/denom : 0.;
      pi += (1. - ssh);
      
    }//segsites
    weighted_pi_ii += weights[i]*weights[i]*pi;
    //    printf("%d, %d, %f, %f, %f\n", i, segsites, pi, weights[i], weighted_pi_ii);
  }//pops
  
  
  //between population divergence
  for(i=0; i<npop-1; i++){
    start_i=samples_b[ i];
    end_i=samples_e[i];
    
    set_almap(almap_i, list, n, segsites, start_i, end_i);
    set_missing(missing_i, list, n, segsites, start_i, end_i);
    
    for (j = i+1; j<npop; j++){
      start_j=samples_b[j];
      end_j=samples_e[j];
      //printf("BETWEEN i %d, start %d, end %d\n", j, start_j, end_j);
      set_almap(almap_j, list, n, segsites, start_j, end_j);
      set_missing(missing_j, list, n, segsites, start_j, end_j);

      pi_ij=0.;
      sum_wi_wj += weights[i] * weights[j];
      segs_ij = 0;

      for(k=0; k<segsites; k++){ // over segsites
	//printf("\n***k=%d\n", k);
	ni=config[i] - missing_i[k];
	nj=config[j] - missing_j[k];
	zeros_i = ni - almap_i[k];
	zeros_j = nj - almap_j[k];
	
	if((ni <=0) || (nj<=0) )
	  continue;
	
	//fprintf(stderr, "ni: %d, nj: %d, zi: %d, zj: %d, oni: %d, onj: %d\n", ni, nj, zeros_i, zeros_j, almap_i[k], almap_j[k]);

	pi_ij += (double)zeros_i*(double)almap_j[k]/(double)ni/(double)nj;
	pi_ij += (double)almap_i[k]*(double)zeros_j/(double)ni/(double)nj;
	
	if( ( (zeros_j + zeros_i) < (ni + nj) ) && 
	    ( (zeros_j + zeros_i) > 0) ){
	  segs_ij++; 
	}

	// private alleles
	//printf("site:%d\tzeros_i: %d, almap_i: %d,  zeros_j: %d, almap_j: %d\n", k, zeros_i, almap_i[k], zeros_j, almap_j[k]);
	if( (!zeros_i && zeros_j && almap_j[k]) || 
	    (!zeros_j && zeros_i && almap_i[k]) ||
	    (!almap_j[k] && almap_i[k] && zeros_i) ||
	    (!almap_i[k] && almap_j[k] && zeros_j) ){
	  private[i][j] = private[i][j] + 1.0;
	  //printf("private: %e\n", private[i][j]);
	}
	// shared alleles
	if( zeros_i && zeros_j && almap_i[k] && almap_j[k]){
	  shared[i][j] = shared[i][j] + 1.0;
	}
        // fixed alleles
	if( ( 
	     (derived == 0) && 
	     ( 
	      (!zeros_i && !almap_j[k]) || 
	       (!zeros_j && !almap_i[k]) 
	       )
	      )
	    ||
	    ( 
	     (derived == 1) && 
	      (!zeros_i && !almap_j[k]) 
	      )
	    ||
	    ( 
	     (derived == 2) && 
	      (!zeros_j && !almap_i[k]) 
	      )
	    ){
	  fixed_dif[i][j] = fixed_dif[i][j] + 1.0;
	}
	// dxy (Hahn, 2018; calculating dxy from unphased data)
	dxy[i][j] +=    ( ((double)zeros_i/(double)ni)*((double)almap_j[k]/(double)nj) )
	           +    ( ((double)almap_i[k]/(double)ni)*((double)zeros_j/(double)nj) );
      }

      // shl (average shared haplotype tract length) and lsh (longest shared haplotype tract length)
      longestRun = 0;
      if (segsites > 0){ // compare each haplotype from pop1 with each haplotype from pop2
	for(h1=0; h1<numHaps; h1++){ // for each haplotype
	  //printf("%i %i %s\n", haploCounts[i][h1], haploCounts[j][h1], haplos[h1]); 
	  if (haploCounts[i][h1] > 0){ // if it's present in pop1
	    for(h2=0; h2<numHaps; h2++){ // compare with each other haplotype
	      if (haploCounts[j][h2] > 0){ // if the current hap is in pop2, then compare
		//printf("%i %i %s %s\n", haploCounts[i][h1], haploCounts[j][h2], haplos[h1], haplos[h2]);
		numIbsSegments = 1; // minimum number of segments = 1
		if (h1 == h2){ // same haplotype.
		  longestRun = 1;
		}
		else{ // analyze differences between haplotypes
		  //printf("%s %s %i\n", haplos[h1], haplos[h2], segsites);
		  startPos = -1; // initialize with "-1", meaning no shared position to begin the "run" with
		  currentTractLength = 0;
		  for(k=0; k<segsites; k++){
		    if(haplos[h1][k] == haplos[h2][k]){ // current position is shared
		      if(startPos == -1){ // if no current start, set to current position. Otherwise, continue on. 
			startPos = k; 
		      }
		    }
		    else{ // reached a position that differs
		      numIbsSegments++;
		      if (startPos == -1){ // if the previous site was also different, then simply evaluate the distance since the last seg site
			if(k == 0){
			  currentTractLength = segPositions[0];
			}
			else{
			  currentTractLength = segPositions[k]-segPositions[k-1];
			}
		      }
		      else if(startPos == 0){ // if current start is position "0", then the whole locus before position "0" is also shared
			currentTractLength = segPositions[k]; // k, not k-1, because all the way up until the current position (not the previous) is shared
		      }
		      else{
			currentTractLength = segPositions[k] - segPositions[startPos-1]; // "-1" because all the way back to the pre-start position is shared
		      }
		      if( currentTractLength > longestRun ){ // is the current run the longest run encountered so far?
			longestRun = currentTractLength;
		      }
		      startPos = -1; // restart the run
		    }
		  }
		  if(startPos != -1){ // the next several lines evaluate the end-chunk of the locus, after looping through the seg sites.
		    currentTractLength = 1 - segPositions[startPos-1]; // if final segsite shared, the whole end of the locus is shared, past the last seg site 
		  }
		  else{
		    currentTractLength = 1 - segPositions[segsites-1]; // if final seg site differs, just take that position to the end of the locus
		  }
		  if( currentTractLength > longestRun ){
		    longestRun = currentTractLength;
		  }
		}
		af1 = (float)haploCounts[i][h1] / (float)config[i]; // many curly braces in this section; make sure this block is correctly inserted 
		af2 = (float)haploCounts[j][h2] / (float)config[j];
		shl[i][j] += (af1*af2) / (double)numIbsSegments;
	      }
	    }
	  }
	}
      }
      else{
	shl[i][j] = 1;
	longestRun = 1;
      }
      if(locusLength){ // if a locfile was supplied, use those locus lengths specified, otherwise output a (0,1) value
	shl[i][j] = shl[i][j] * (double)locusLength;
	longestRun = longestRun * (double)locusLength;
      }
      lsh[i][j] = longestRun;


      
      if(segs_ij > 0){
	shared[i][j] /= (double)segs_ij;
	fixed_dif[i][j] /= (double)segs_ij;
	private[i][j] /= (double)segs_ij;
      }
      else if(segs_ij == 0){
	shared[i][j] =	fixed_dif[i][j] = private[i][j] = dxy[i][j] = 0./0.;
      }
      weighted_pi_ij += weights[i]*weights[j]*pi_ij;
    }//pops_j
  }//pops_i

  *piT = weighted_pi_ii + 2.0*weighted_pi_ij;
  *piS = weighted_pi_ii / w_ii_sq;
  *piD = (*piT - *piS)/(2.0 * sum_wi_wj);
  *piB = weighted_pi_ij / sum_wi_wj;
  /* printf("\n"); */
/*   printf("piT: %e\n", *piT); */
/*   printf("piS: %e\n", *piS); */
/*   printf("piD: %e\n", *piD); */
/*   printf("piB: %e\n", *piB); */
/*   printf("sum_wi_wj: %e\n", sum_wi_wj); */
/*   printf("weighted_pi_ii: %e\n", weighted_pi_ii); */
/*   printf("weighted_pi_ij: %e\n", weighted_pi_ij); */


  free(almap_i);
  free(almap_j);
  free(missing_i);
  free(missing_j);
  for (i=0; i<n; i++){
    free(haplos[i]);
  }
  free(haplos);
  return 1;
}
	

double Fst_HSM(double piD, 
	       double piS){
  return piD/(piS + piD);
}

double Fst_Slatkin(double piD,
		   double piS){
  return piD/(2.0*piS + piD);
}

double Fst_HBK(double piS,
	       double piT){
  //  printf("piS:%e\tpiT:%e\n", piS, piT);
  return 1. - (piS/piT);
}

 


/**
   Calculations of pairwise Fst Values
   It needs the weights for the populations because the total frequency is needed to be calculated
   The weights can be obtained by the initial command.

   Code adapted from libsequence
*/
int pairwiseFstcalculations(int popi, int popj, // the two populations
			    double* weights, // larger pops have greater weight
			    char** list, // the polymorphic table
			    int* config, // the configuration of the sample
			    int npop,
			    int segsites,
			    int n,
			    
			    double* piT,
			    double* piS,
			    double* piB,
			    double* piD
			    ){ 

  

  //
  

  int i,
    j,
    k,
    start=0,
    end=0,
    ni, nj,
    segs_ij=0,
    start_i,
    end_i,
    start_j,
    end_j,
    zeros_i,
    zeros_j;

  double w_ii_sq =0.,
    denom=0.,
    pi=0.,
    pi_ij=0.,
    sum_wi_wj=0.,
    weighted_pi_ij = 0.,
    weighted_pi_ii = 0.;

  
  static int* samples_b = NULL;
  static int* samples_e = NULL;
  static int first = 1;
  
  int total_npop = npop; // save the total number of populations into npop;
  npop = 2;
  int popind[2]; // save the indices of populations
  popind[0] = popi;
  popind[1] = popj;

  end=0;
  start=0;

  if( first == 1){
    samples_b = (int*)malloc(total_npop * sizeof(int));
    samples_e = (int*)malloc(total_npop * sizeof(int));
    first = 0;
  }
  
  samples_b[0] = 0;
  samples_e[0] = samples_b[0] + config[0];
  for(i=1; i<total_npop; i++){
    samples_b[i] = samples_b[i-1] + config[i-1];
    samples_e[i] = samples_b[i] + config[i];
    //printf("config %d is %d\n", i, config[i]);
  }
  
  // printf("segsites: %d\n", segsites);

  int* almap_i = (int*)calloc(segsites,sizeof(int));
  int* almap_j = (int*)calloc(segsites,sizeof(int));

  int* missing_i = (int*)calloc(segsites, sizeof(int));
  int* missing_j = (int*)calloc(segsites, sizeof(int));
  
 

 
  //printf("Fst function\n");
  for(i=0; i<npop; i++){ // over pops
    

    start = samples_b[ popind[i] ];
    end = samples_e[ popind[i] ];
    /* printf("*i %d, start %d, end %d\n", popind[i], start, end); */
/*     printf("*%d, %d, %d, %d, %e\n", i, segsites, start, end, weights[i]); */
    set_almap(almap_i, list, n, segsites, start, end); // set the allelic map for the population
    set_missing(missing_i, list, n, segsites, start, end);
    pi=0.;
    w_ii_sq += ((double)total_npop/2.0 * weights[ popind[i] ] )* ((double)total_npop/2.0 *weights[ popind[i] ]);
    //w_ii_sq += .5 * .5;
    //fprintf(stderr, "weight %e, total: %d, weights: %e, %e\n", w_ii_sq, total_npop, weights[ popind[i] ], weights[ popind[i] ]);
    for(j=0; j<segsites; j++){ // over pol sites
      ni=config[ popind[i] ] - missing_i[j]; // configuration of the population
      if(ni < 2) continue;
      double ssh=0.;
      denom = (double)ni * ((double)ni-1.);


      ssh+= (ni-almap_i[j]>0) ?
	(double)(ni-almap_i[j])*(double)(ni-almap_i[j] - 1)/denom : 0.;
      ssh+=(almap_i[j] > 0)?(double)almap_i[j] * (double)(almap_i[j] - 1)/denom : 0.;
      pi += (1. - ssh);
      
    }//segsites
    weighted_pi_ii += weights[ popind[i] ]*weights[ popind[i] ]*pi;
    //printf("*%d, %d, %e\n", i, segsites, weights[i]);
  }//pops
  
  //between population divergence
  for(i=0; i<npop-1; i++){
    start_i=samples_b[ popind[i] ];
    end_i=samples_e[ popind[i] ];
    
    set_almap(almap_i, list, n, segsites, start_i, end_i);
    set_missing(missing_i, list, n, segsites, start_i, end_i);

    for (j = i+1; j<npop; j++){
      start_j=samples_b[ popind[j] ];
      end_j=samples_e[ popind[j] ];
      //printf("BETWEEN i %d, start %d, end %d\n", j, start_j, end_j);
      set_almap(almap_j, list, n, segsites, start_j, end_j);
      set_missing(missing_j, list, n, segsites, start_j, end_j);

      pi_ij=0.;
      sum_wi_wj += weights[ popind[i] ] * weights[popind[j] ];
      segs_ij = 0;

      for(k=0; k<segsites; k++){ // over segsites
	ni=config[ popind[i] ] - missing_i[k];
	nj=config[ popind[j] ] - missing_j[k];
	zeros_i = ni - almap_i[k];
	zeros_j = nj - almap_j[k];
	if((ni<=0) || (nj<=0) ) continue;

	pi_ij += (double)zeros_i*(double)almap_j[k]/(double)ni/(double)nj;
	pi_ij += (double)zeros_j*(double)almap_i[k]/(double)ni/(double)nj;
      }// sites
      weighted_pi_ij += weights[ popind[i] ] * weights[ popind[j] ] * pi_ij;
    }//pops_j
  }//pops_i

  *piT = weighted_pi_ii + 2.0*weighted_pi_ij;
  *piS = weighted_pi_ii / w_ii_sq;
  *piB = weighted_pi_ij / sum_wi_wj;
  *piD = (*piT - *piS)/(2.0 * sum_wi_wj);
  /* printf("*\n"); */
   /* printf("*piT: %e\n", *piT);  */
/*    printf("*piS: %e\n", *piS);  */
/*    printf("w_ii_sq: %e\n", w_ii_sq); */
/*   printf("*piD: %e\n", *piD); */
/*   printf("*piB: %e\n", *piB); */
/*   printf("*sum_wi_wj: %e\n", sum_wi_wj); */
/*   printf("*weighted_pi_ii: %e\n", weighted_pi_ii); */
/*   printf("*weighted_pi_ij: %e\n", weighted_pi_ij); */

  free(almap_i);
  free(almap_j);
  free(missing_i);
  free(missing_j);
  return 1;
}

	   
	
int thetaH(double *h, char** list, int n, int segs, double bm, int* almap, int start, int end, int* missing){
  int i, j, k;
  double f = 0.;
  int noc=0;

 
  //printf("******************\n");
        
  assert(start >=0 && end <= n);
  int nsam = end - start ; // start is zero and end the sample size


  if(bm <= 0.) // if no back mutation use the allelic map.
    for(i=0; i<segs; i++){
      nsam = end - start - missing[i];
      if(nsam < 2) continue;
      if(almap[i] == 0 || almap[i] == nsam) 
	continue;
      noc = almap[i];
      f += (2.*(double)noc*(double)noc)/((double)nsam*((double)nsam-1.));
      
      //fprintf(stderr, "noc: %d, nsam: %d, end : %d, start: %d\n", 
      //      noc, nsam, end, start);
      
      
    }
  /* if there is the possibility for back mutation then the allelic map
     is not exactly the one we observe.
     In this case the number of mis-inferred states is binomially distributed.
  */
  else if(bm > 0) {
    for(i=0; i<segs; i++){
      nsam = end - start - missing[i];
      noc = 0;
      for(j=start; j<end; j++){
	if(rand() > bm && (list[j][i] != MISSING) )
	  noc += (list[j][i] - '0');
	else if(list[j][i] != MISSING) 
	  noc += ('1' - list[j][i]);
      }
      if(noc == 0 || noc == nsam) continue;
      f += (2.*(double)noc*(double)noc)/((double)nsam*((double)nsam-1.));

      
    }
  }

  *h = f;
  return 1;
}


int hDenominator(double* hden, int n, int segs, double hn, double sqhn, double bn){
  double b1, b2, c1, c2, e1, e2, w;
  double den;
  b1 = (n+1.)/(3. * (n-1));
  b2 = (double)2*(n*n + n+3.)/(double)(9. * n * (n-1));
  c1 = b1 - 1./hn;
  c2 = b2 - (n + 2.)/(1. * n * hn) + sqhn/(hn*hn);
  e1 = c1/hn;
  e2 = c2/(hn*hn + sqhn);
  w = ((double) segs)/hn;
  
  *hden = sqrt((n - 2.)*w/(6*(n-1.)) + (18 * n*n*(n*3 + 2)*bn - (88*n*n*n + 9.*n*n - 13.*n+6)) * w*w/(9.*n*(n-1) * (n-1)));
  return 1;
}
  

int htest( double* h, int n, int segs, double thetaPi, double thetaH, double hn, double sqhn, double bn){
  double hden = 0;
  double fwh=0.;
    
  hDenominator(&hden, n, segs, hn, sqhn, bn);
  
  //fprintf(stderr, "thetapi: %e\tthetaH: %e\thden: %e\n", thetaPi, thetaH, hden);
  
  fwh = 0.5*(thetaPi - thetaH)/hden;
  
  *h = fwh;
  return 1;

}


int dvstat( char** list, int n, int segs, int start, int end, double* dvk, double* dvh){
  assert(end <= n);
  assert(start >= 0);
  int i,j, nsam=end-start, k=end-start;
  int* haplo = malloc((unsigned)n * sizeof(int));
  
  *dvk = k;
  *dvh = 1.0;

  for(i=0; i<n; i++)
    haplo[i] = 1;
  
  for(i=start; i<end-1; i++){
    if(!haplo[i]) continue;
    for(j=i+1; j<end; j++){
      if(!haplo[j]) continue;
      if(!mystrcmp(list[i], list[j], segs)){
	haplo[j] = 0;
	haplo[i]++;
	continue;
      } 
    }
  }
  
  
  for(i=start; i<end; i++){
    //printf("\n--i: %d, haplo: %d, nsam: %d, dvk: %e\n", i, haplo[i], nsam, *dvk);
    if(!haplo[i]){
      *dvk = *dvk - 1.;
      continue;
    }
    *dvh -= pow((double)haplo[i]/(double)nsam, 2);
    
  }
   
  free(haplo);
  return 1;
}


int mystrcmp(char* p1, char* p2, int length) 
{ 
  // A variable to hold the distance in ASCII charcters 
  // between the current two characters we're comparing. 
  int dist = 0, i=0;; 
  
  // Keep checking while the distance is 0, and we're 
  // not hitting NULL on one of the strings. 
  while (!dist && i<length){
    // Get the distance while incrementing the 
    // pointers to the next character. 
    dist = (*p2++) - (*p1++); 
    i++;
  }
  // Check the last distance and according to this 
  // return (1) if the first string is bigger, (-1) 
  // if the second string is bigger, or (0) if the 
  // strings are identical. 
  if (dist > 0) 
    return (-1); 
    else if (dist < 0) 
      return (1); 
  
  return (0); 
}


      
    
      
