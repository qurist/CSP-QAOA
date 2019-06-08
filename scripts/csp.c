/*---------------------------------------------------------------------------
  Title    : Simulator and data analysis tools for Quantum Approximate 
  Optimization Algorithm and Variational Quantum Eigensolvers.
  Authors  : Aniruddha Bapat (QuICS) and Stephen Jordan (Microsoft QuARC) 
  Year     : 2018
  Citation : If you use this code in your research, please cite it.
  -------------------------------------------------------------------------*/

#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include "csp.h"

//Print out the integer x in binary
void printbits(int x, int n) {
  int i;
  for(i = n-1; i >= 0; i--) printf("%i", (x>>i)&1);
}

//Print out the quantum state vector
void printvec(state psi) {
  int i;
  double p;
  for(i = 0; i < psi.N; i++) {
    p = psi.realcur[i]*psi.realcur[i] + psi.imagcur[i]*psi.imagcur[i];
    printf("%i(", i);
    printbits(i, psi.n);
    printf("):\t%e\t+i%e\t%e\n", psi.realcur[i], psi.imagcur[i], p);
  }
}

//Return the Euclidean norm of the state vector
double norm(state psi) {
  double val;
  int i;
  val = 0;
  for(i = 0; i < psi.N; i++) {
    val += psi.realcur[i]*psi.realcur[i];
    val += psi.imagcur[i]*psi.imagcur[i];
  }
  return val;
}

//Swap the current state with the buffered state
void swapbuf(state *psi) {
  double *tmp;
  tmp = psi->realcur;
  psi->realcur = psi->realbuf;
  psi->realbuf = tmp;
  tmp = psi->imagcur;
  psi->imagcur = psi->imagbuf;
  psi->imagbuf = tmp;
}  

//Perform a Hadamard transform on the n qubits in O(NlogN) time.
void Hadamard(state *psi) {
  int i,j;
  double root;
  root = 1.0/sqrt(2.0);
  for(i = 0; i < psi->n; i++) {
    //apply a Hadamard gate to the ith qubit and write result into buffer
    for(j = 0; j < psi->N; j++) {
      if((j&(1<<i)) == 0) {
        psi->realbuf[j] = root*(psi->realcur[j] + psi->realcur[j^(1<<i)]);
	psi->imagbuf[j] = root*(psi->imagcur[j] + psi->imagcur[j^(1<<i)]);
      }
      else {
        psi->realbuf[j] = root*(-psi->realcur[j] + psi->realcur[j^(1<<i)]);
	psi->imagbuf[j] = root*(-psi->imagcur[j] + psi->imagcur[j^(1<<i)]);
      }
    }
    swapbuf(psi); //swap cur with buf
  }
}

//Computes the Z2 inner product between bitstrings a and b of length n.
int Z2inner(int A, int B, int n) {
  int i;
  int telescope;
  int AandB;
  AandB = A&B;
  telescope = 0;
  for(i = 0; i < n; i++) telescope ^= (AandB>>i);
  return telescope&1;
}

//Convert an array of 1s and 0s to the corresponding number via binary place value.
//The zeroth element of the array contains the 1s place, the 1th element the 2s etc.
int b2i(int *bitstring, int n) {
  int i;
  int val;
  val = 0;
  for(i = 0; i < n; i++) val += bitstring[n-i-1]*(1<<i);
  return val;
}

//Convert a string of characters like "1100" to the corresponding number like 12.
//Note that in keeping with standard conventions, this is the reverse order for
//significant bits than in bin2int.
int cb2i(char *bitstring, int n) {
  int i;
  int val;
  val = 0;
  //the ascii code for '0' is 48, for '1' is 49, etc.
  for(i = 0; i < n; i++) val += ((int)bitstring[n-i-1]-48)*(1<<i);  
  return val;
}

//Returns the number of ones in the binary expansion of x
int HammingWeight(int x, int n) {
  int weight;
  int i;
  weight = 0;
  for(i = 0; i < n; i++) weight += (x>>i)&1;
  return weight;
}

//Returns b_i*x_i for x_i a bit in the binary expansion of x
double SkewWeight(int x, int n, double *b) {
  double weight;
  int i;
  weight = 0;
  for(i = 0; i < n; i++) weight += b[i]*((x>>i)&1);
  return weight;
}

//Unitarily evolve for time t according to H = - sum_j X_j
void evolveX(state *psi, double t) {
  int i;
  double angle;
  //Hadamard transform into the basis in which H is diagonal
  Hadamard(psi);
  //Do the time evolution in this basis
  for(i = 0; i < psi->N; i++) {
    angle = t*(double)(psi->n - 2*HammingWeight(i,psi->n));
    //multiply the amplitude by e^i(angle)
    psi->realbuf[i] = cos(angle)*psi->realcur[i] - sin(angle)*psi->imagcur[i];
    psi->imagbuf[i] = sin(angle)*psi->realcur[i] + cos(angle)*psi->imagcur[i];
  }
  swapbuf(psi);  //swap cur with buf
  Hadamard(psi); //Hadamard transform back to the computational basis
}

//The energy associated with a bit string according to a given CSP hamiltonian
//x stores the bit string
//n is the number of bits
double energy(int x, cspham *H) {
  int i;
  int j;
  int bitvar;
  int n, m;
  int var,sign,proj;
  double Esum;
  n   = H->n;
  m   = H->m;
  Esum = 0;
  int *clause;
  clause = (int *)malloc(n*sizeof(int));
  
  for(i=0;i<m;i++){//run through clause pointers
    clause = H->clauses[i];
    proj=1; //unsatisfied unless found otherwise
    
    for(j=0;j<H->lengths[i];j++){//run thru variables
      sign=0; // positive unless found otherwise
      var  = clause[j];

      if(var<0){
	sign = 1;
	var  = -var; //flip!
      }
      // now we have a positive var
      bitvar = (x>>(var-1))&1; //x bit at variable's index
      if(bitvar^sign){//if clause is satisfied
	proj=0;
	break;
      }
    }
    Esum+= H->weights[i]*(double)proj; //Add energy contribution with weights
  }
  return Esum;
}

//Unitarily evolve for time t according to H = \sum_i<j J_ij Z_i Z_j
void evolveZ(state *psi, cspham *H, double t) {
  int i;
  double angle;
  for(i = 0; i < psi->N; i++) {
    angle = (-t)*H->diag[i];
    //angle = (-t)*energy(i, H);
    psi->realbuf[i] = cos(angle)*psi->realcur[i] - sin(angle)*psi->imagcur[i];
    psi->imagbuf[i] = sin(angle)*psi->realcur[i] + cos(angle)*psi->imagcur[i];
  }
  swapbuf(psi);  //swap cur with buf
}

//Remove bruised zeros, i.e. things that are exactly zero but appear to
//be 10^-15 or something due to numerical noise.
void debruise(state *psi) {
  int i;            //counter variable
  double threshold; //amplitudes smaller than this will be zeroed
  threshold = 1.0E-12;
  for(i = 0; i < psi->N; i++) {
    if(fabs(psi->realcur[i]) < threshold) psi->realcur[i] = 0;
    if(fabs(psi->imagcur[i]) < threshold) psi->imagcur[i] = 0;
  }
}

//Allocate memory for the amplitudes and set N = 2^n
//Returns 1 on success, 0 on failure
int allocate(state *psi, int n) {
  psi->n = n;
  psi->N = 1<<n;
  psi->realcur = (double *)malloc(psi->N*sizeof(double));
  if(psi->realcur == NULL) return 0;
  psi->imagcur = (double *)malloc(psi->N*sizeof(double));
  if(psi->imagcur == NULL) return 0;  
  psi->realbuf = (double *)malloc(psi->N*sizeof(double));
  if(psi->realbuf == NULL) return 0;  
  psi->imagbuf = (double *)malloc(psi->N*sizeof(double));
  if(psi->imagbuf == NULL) return 0;
  return 1;
}

//Deallocate the memory for the amplitudes
void deallocate(state *psi) {
  free(psi->realcur);
  free(psi->imagcur);
  free(psi->realbuf);
  free(psi->imagbuf);
}

// A sketch for a function that reads in a MAXSAT instance and creates a generalized Z
// Hamiltonian
void readcnf(char *path, cspham *H){
  // Instance variables
  FILE *cnf;
  int    n; 
  int    m;
  // Counters, buffers, and other useful constants
  int    i, varno,clno; 
  int    lenbuff=128;
  char   ptype[lenbuff];
  char   buff[lenbuff];
  int    ibuf;
  double dbuf;
  int    *clause; 

  cnf = fopen(path, "r");
  clno = 0;
  while(!feof(cnf)){
    // FIRST CHARACTER DETERMINES IF LINE IS DATA OR NOT
    fscanf(cnf, "%s", buff);
    
    // PARSER FOR COMMENTS, META-DATA, AND TRASH
    if(buff[0]=='\n'){
      continue;
    }
    else if(buff[0]=='c'){
      fgets(buff, lenbuff, cnf);
      continue;
    }
    else if(buff[0]=='p'){
      fgets(buff, lenbuff, cnf);
      sscanf(buff, "%s %d %d\n", ptype, &n, &m);

      H->n       = n;
      H->N       = 1<<n;
      H->m       = m;
      H->clauses = (int   **)malloc(m*sizeof(int *));
      H->weights = (double *)malloc(m*sizeof(double));
      H->lengths = (int    *)malloc(m*sizeof(int));
      H->maxsat  = 0;
      H->diag = (double *)malloc(H->N*sizeof(double));


      fprintf(stderr, "Read a CSP instance of %d clauses on %d variables.\n", m, n);
      continue;
    }
    // CLAUSE PARSER
    
    // Unweighted cnf
    if(ptype[0]=='c'){
      clause  = (int *)malloc((n+1)*sizeof(int));

      sscanf(buff,"%d", &ibuf);
      clause[0]=ibuf;
      varno=1;

      while(ibuf!=0){
	fscanf(cnf, "%d", &ibuf);
	clause[varno]=ibuf;
	varno++;
      }

      fgets(buff, lenbuff, cnf); 
      H->clauses[clno] = clause; 
      H->weights[clno] = 1.0;
      H->lengths[clno] = varno-1;
      H->maxsat+=1;

      clno++;
      continue;  
    }
    
    // Weighted cnf instance
    else if(ptype[0]=='w'){
      clause  = (int    *)malloc((n+1)*sizeof(int));

      sscanf(buff, "%lf", &dbuf);
      H->weights[clno] = dbuf;
      H->maxsat+=dbuf;

      varno=0;
      ibuf =1;
      while(ibuf!=0){
	fscanf(cnf, "%d", &ibuf);
	clause[varno]=ibuf;
	varno++;
      }

      fgets(buff, lenbuff, cnf); // read the newline
      H->clauses[clno] = clause; // pass reference of clause to clauses 
      H->lengths[clno] = varno-1;

      clno++;
    }
  }
  fclose(cnf);
} 

//Deallocate the memory for the amplitudes
void deallocateH(cspham  *H) {
  int i;
  free(H->weights);
  free(H->lengths);
  for(i=0;i<H->m;i++){
    free(H->clauses[i]);
  }
  free(H->clauses);
  free(H->diag);
}

//Initialize psi to the uniform superposition
void uniform(state *psi) {
  int i;
  double amp;
  amp = 1.0/sqrt((double)psi->N);
  for(i = 0; i < psi->N; i++) {
    psi->realcur[i] = amp;
    psi->imagcur[i] = 0;
    psi->realbuf[i] = 0;
    psi->imagbuf[i] = 0;
  }
}  

//Initialize psi to the all zeros bit string
void zeros(state *psi) {
  int i;
  psi->realcur[0] = 1.0;
  psi->imagcur[0] = 0;
  psi->realbuf[0] = 0;
  psi->imagbuf[0] = 0;
  for(i = 1; i < psi->N; i++) {
    psi->realcur[i] = 0;
    psi->imagcur[i] = 0;
    psi->realbuf[i] = 0;
    psi->imagbuf[i] = 0;
  }
}

void energize(cspham *H){
  int i;
  for(i=0;i<H->N;i++){
    H->diag[i] = energy(i,H);
  }
}

//Return the expectation value of the Hamiltonian in state psi. The Hamiltonian is
//a CSP.
double expectH(state psi, cspham * H) {
  int x;
  double p;
  double expectZ;
  expectZ = 0;
  for(x = 0; x < psi.N; x++) {
    p = psi.realcur[x]*psi.realcur[x] + psi.imagcur[x]*psi.imagcur[x];
    //expectZ += p*energy(x,H);    
    expectZ += p*H->diag[x];
  }
  return expectZ;
}

// Evolve the state by a p-iteration QAOA protocol
void evolveQAOA(state *psi, cspham *H, double * beta, double * gamma, int p) {
  int i;
  for(i=0;i<p;i++) {
    evolveZ(psi, H, gamma[i]);
    evolveX(psi,    beta[i]);
  }    
}
