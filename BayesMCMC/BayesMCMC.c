#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MIN(x,y) x<y?x:y
#define kb 8.6173324E-5 /* in eV/K */
#define hbar 4.1356677E-15 /* in eV s */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


double ran1(long *idum);
double gasdev(long *idum);
double * init_array(int n_steps, double val);
double * load_matrix(char *filename, int *n, int *m);
void print_array(double *x, int n_steps);
void print_array2(int *x, int n_steps);
void print_matrix(double *A,  int n_rowA, int n_colsA);
double gammln(double xx);
double factrl(int n);
double *matrix_multiply(double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
double *matrix_multiply2(double *A,int *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
double sum_array(double *x, int N);
double std_array(double *x, int N);
double * std_array_along(double *A, int axis, int n_rowA, int n_colsA);
double * sum_array_along(double *A, int axis, int n_rowA, int n_colsA);
double * array_multiply(double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
double * array_multiply2(double *A, int *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
double * array_divide(double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
double * array_sum(double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
double * array_log(double *A,int n_rowA, int n_colsA);
double * array_exp(double *A,int n_rowA, int n_colsA);
double * array_factrl(int *A,int n_rowA, int n_colsA);
double * scalar_sum(double a, double *A,int n_rowA, int n_colsA);
double * scalar_multiply(double a, double *A,int n_rowA, int n_colsA);
double * array_normal(long *idum ,int n_rowA, int n_colsA);
double calculate_likelyhood(double *f,double *V,double *onesNbins , int *N ,double *Nfac,int MaxPop, int Nbins, int Tframes, double *histo);

int main(int argc, char **argv){


	int n_steps;
	int i;
	int j;
	int n_row, n_cols;
	
	
	int MaxPop, Nbins;
	int *N;
	double *Nfac;
	int Tframes;
	double *histo;
	double *f;
	double *V;
	double *fprime;
	double *Vprime;
	double *onesNbins;
	double *loglikely;
	double delta_l;
	double alpha;
	
	long *idum;
	long ed=-50;
	idum = &ed;
	
	
	MaxPop=atoi(argv[1]);
	Nbins=atoi(argv[2]);
	Tframes=atoi(argv[3]);
	n_steps=atoi(argv[4]);
	f = array_normal(idum ,MaxPop+1, 1);
	V = array_normal(idum ,Nbins, 1);
	
	N= malloc((MaxPop+1)*sizeof(int));
	for(i=0;i<MaxPop+1;i++){ N[i]=i; }
	
	Nfac=array_factrl(N,MaxPop+1, 1);
	
	onesNbins=init_array(Nbins,1.0);
	
	/*initialize */
	loglikely = init_array(n_steps,0.0);
	histo = array_normal(idum ,Nbins,MaxPop+1);
	

	/*first step*/
	loglikely[0]=calculate_likelyhood(f,V,onesNbins , N , Nfac,  MaxPop, Nbins,  Tframes, histo);
	
	/*monte carlo thermalization step*/
	for(i=1;i<n_steps;i++){
		loglikely[i]=loglikely[i-1];
	
		delta_l=calculate_likelyhood(f,V,onesNbins , N , Nfac,  MaxPop, Nbins,  Tframes, histo)-calculate_likelyhood(f,V,onesNbins , N , Nfac,  MaxPop, Nbins,  Tframes, histo);
		
		if (delta_l<=0.0) {
			loglikely[i]= calculate_likelyhood(fprime,Vprime,onesNbins , N , Nfac,  MaxPop, Nbins,  Tframes, histo);
			
			
		}
		else{
			alpha = drand48();
			if (log(alpha)<delta_l) {
			
				loglikely[i]= calculate_likelyhood(fprime,Vprime,onesNbins , N , Nfac,  MaxPop, Nbins,  Tframes, histo);
			}
		}
		
	


		
		
	}
	

	
	
	/*
	double *csv;
	double *csv2;
	csv=array_normal(idum, 3, 3);
	csv2=scalar_sum(100.0,scalar_multiply(10.0,array_normal(idum, 3, 3),3,3),3,3);
	print_array(csv2,9);
	print_matrix(csv2,3,3);
	print_matrix(sum_array_along(csv2,1,3,3),3,1);
	 */
	return(0);
	
}



double calculate_likelyhood(double *f,double *V,double *onesNbins , int *N ,double *Nfac,int MaxPop, int Nbins, int Tframes, double *histo)
{
	double loglike=0.0;
	double *logz;
	double *TeoHist;
	double *exponent;
	double *fmeanexp;
	double *VNexpav;
	int n_rowA, n_colsA;
	
	double *temp1,*temp2,*temp3,*temp4,*temp5,*temp6,*temp7,*temp8,*temp9;
	
	n_rowA=Nbins;
	n_colsA=MaxPop+1;
	
	temp1 = init_array(n_steps,0.0);
	temp2 = init_array(n_steps,0.0);
	temp3 = init_array(n_steps,0.0);
	temp4 = init_array(n_steps,0.0);
	temp5 = init_array(n_steps,0.0);
	temp6 = init_array(n_steps,0.0);
	temp7 = init_array(n_steps,0.0);
	temp8 = init_array(n_steps,0.0);
	temp9 = init_array(n_steps,0.0);
	
	exponent=scalar_multiply(-1.0, array_sum( matrix_multiply2( V, N, n_rowA , 1, 1, n_colsA),  matrix_multiply( onesNbins, f, n_rowA , 1, 1, n_colsA), n_rowA, n_colsA, n_rowA, n_colsA), n_rowA, n_colsA);
	TeoHist=array_divide(array_exp(exponent, n_rowA, n_colsA),matrix_multiply( onesNbins, Nfac, n_rowA , 1, 1, n_colsA), n_rowA, n_colsA, n_rowA, n_colsA);
	logz=array_log(sum_array_along(TeoHist, 1, n_rowA, n_colsA),n_rowA, n_colsA);
	fmeanexp=matrix_multiply( histo, f, n_rowA , n_colsA,n_colsA , 1);
	VNexpav=array_multiply(V,matrix_multiply2( histo, N, n_rowA , n_colsA,n_colsA , 1), n_rowA,1,n_rowA,1);
	
	loglike=-Tframes*(sum_array(logz, n_rowA )+sum_array(VNexpav, n_rowA )+sum_array(fmeanexp, n_rowA ));
	
	return loglike;
	
}




double ran1(long *idum)
/*
 “Minimal” random number generator of Park and Miller with Bays-Durham shuffle and added
 safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter idum between
 successive deviates in a sequence. RNMX should approximate the largest floating value that is less than 1.*/
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;
	if (*idum <= 0 || !iy) { /*Initialize.*/
		if (-(*idum) < 1) *idum=1; /*Be sure to prevent idum = 0.*/
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) { /*Load the shuffle table (after 8 warm-ups).*/
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ; /*Start here when not initializing.*/
	*idum=IA*(*idum-k*IQ)-IR*k; /*Compute idum=(IA*idum) % IM without overflows*/
	if(*idum < 0) *idum += IM; /* by Schrage’s method.*/
	j=iy/NDIV; /*Will be in the range 0..NTAB-1.*/
	iy=iv[j]; /*Output previously stored value and refill the shuffle table.*/
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX; /*Because users don’t expect endpoint values.*/
	else return temp;
}

double gasdev(long *idum)
/*Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum)
 as the source of uniform deviates.*/
{
	double ran1(long *idum);
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;
	if (*idum < 0) iset=0; /*Reinitialize.*/
	if (iset == 0) { /*We don’t have an extra deviate handy, so*/
		do {
			v1=2.0*ran1(idum)-1.0; /*pick two uniform numbers in the square extending from -1 to +1 in each direction,*/
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2; /*see if they are in the unit circle, and if they are not, try again.*/
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		/*Now make the Box-Muller transformation to get two normal deviates. Return one and
		 save the other for next time.*/
		gset=v1*fac;
		iset=1; /*Set flag.*/
		return v2*fac;
	} else { /*We have an extra deviate handy,*/
		iset=0; /*so unset the flag,*/
		return gset; /*and return it.*/
	}
}



double * init_array(int n_steps,double val){
	int i;
	double *x;
	if(!(x=malloc(sizeof(double) * n_steps))){
		fprintf(stderr,"problem with malloc\n");
		exit(1);
	}
	for(i=0;i<n_steps;i++){
		x[i] = val;
	}
	return x;
}

double *load_matrix(char *filename, int *n, int *m){
	double *matrix;
	FILE *in;
	int n_row, n_cols;
	int i;
	int j;
	
	if(!(in=fopen(filename, "r"))){
		printf("Problem opening file %s\n", filename);
		exit(1);
	}
	
	fscanf(in, "%d %d\n", &n_row, &n_cols);
	/*printf("%d %d\n", n_row, n_cols);*/
	
	matrix = malloc(n_row * n_cols * sizeof(double));
	
	for(i=0;i<n_row;i++){
		for(j=0;j<n_cols;j++){
			fscanf(in, "%f", &matrix[i*n_cols + j]);
		}
	}
	*n = n_row;
	*m = n_cols;
	return matrix;
}


void print_array(double *x, int n_steps){
	int i;
	for(i=0;i<n_steps-1;i++){
		fprintf(stdout, "%f,", x[i]);
	}
	fprintf(stdout, "%f\n", x[n_steps-1]);
	}

void print_array2(int *x, int n_steps){
	int i;
	for(i=0;i<n_steps-1;i++){
		fprintf(stdout, "%d,", x[i]);
	}
	fprintf(stdout, "%d\n", x[n_steps-1]);
}


void print_matrix(double *A,  int n_rowA, int n_colsA){
	int i,j;
	for(i=0;i<n_rowA;i++){
		for(j=0;j<n_colsA-1;j++){
			fprintf(stdout, "%f,", A[i*n_colsA + j]);
		}
	fprintf(stdout, "%f\n", A[i*n_colsA + (n_colsA -1)]);
	}
}



double gammln(double xx)
/*Returns the value ln[gamma(xx)] for xx > 0.*/
{
	/*Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
	 accuracy is good enough.*/
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
		return -tmp+log(2.5066282746310005*ser/x);
		}

double factrl(int n)
/*Returns the value n! as a floating-point number.*/
{
	double gammln(double xx);
	void nrerror(char error_text[]);
	static int ntop=4;
	static double a[33]={1.0,1.0,2.0,6.0,24.0}; /*Fill in table only as required.*/
	int j;
	if (n < 0) printf("Negative factorial in routine factrl");
		if (n > 32) return exp(gammln(n+1.0));
			/*Larger value than size of table is required. Actually, this big a value is going to overflow
			on many computers, but no harm in trying.*/
	while (ntop<n) { /*Fill in table up to desired value.*/
		j=ntop++;
		a[ntop]=a[j]*ntop;
	}
	return a[n];
}


double *matrix_multiply(double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB)
{
	int i, j,k;
	static double *matrix_C;
	matrix_C = malloc(n_rowA * n_colsB * sizeof(float));
	
	
	for(i=0;i<n_rowA;i++){
		for(j=0;j<n_colsB;j++){
			matrix_C[i*n_colsB + j]=0.0;
			for(k=0;k<n_colsA;k++){
				matrix_C[i*n_colsB + j] = matrix_C[i*n_colsB + j]+A[i*n_colsA + k]*B[k*n_colsB + j];
			}
			/*printf(" %f ", matrix_C[i*n_colsB + j]);*/
		}
		/*printf("\n");*/
	}

	return matrix_C;
}

double *matrix_multiply2(double *A,int *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB)
{
	int i, j,k;
	static double *matrix_C;
	matrix_C = malloc(n_rowA * n_colsB * sizeof(float));
	
	
	for(i=0;i<n_rowA;i++){
		for(j=0;j<n_colsB;j++){
			matrix_C[i*n_colsB + j]=0.0;
			for(k=0;k<n_colsA;k++){
				matrix_C[i*n_colsB + j] = matrix_C[i*n_colsB + j]+A[i*n_colsA + k]*B[k*n_colsB + j];
			}
			/*printf(" %f ", matrix_C[i*n_colsB + j]);*/
		}
		/*printf("\n");*/
	}
	
	return matrix_C;
}

double sum_array(double *x, int N){
	int i;
	double Sum=0.0;
	for (i=0; i<N; i++) {
		Sum+=x[i];
	}
	
	return Sum;
	
	
}

double std_array(double *x, int N)
{
	double sum = 0.0, mean, tmp = 0.0;
	
	int i;
	
	for(i=0; i<N; ++i)
	{
		sum += x[i];
	}
	mean = sum/N;
	
	for(i=0; i<N; ++i)
		tmp += pow(x[i] - mean, 2);
	
	return sqrt(tmp/N);
}


double * std_array_along(double *A, int axis, int n_rowA, int n_colsA)
{
	double sum = 0.0, mean, tmp = 0.0;
	static double *matrix_C;
	int i,j;
	
	if(axis==0){
		matrix_C = malloc(n_colsA* sizeof(double));
		for(j=0;j<n_colsA;j++){
			
			matrix_C[j]=0.0;
			mean=0.0;
			
			for(i=0;i<n_rowA;i++){
				mean=matrix_C[j]+A[i*n_colsA + j];
			}
			mean=mean/n_rowA;
			
			for(i=0; i<n_rowA; ++i){
				tmp += pow(A[i*n_colsA + j] - mean, 2);
			}
			matrix_C[j]=sqrt(tmp/n_rowA);
		}
		return matrix_C;
		
	}
	else if(axis==1){
		matrix_C = malloc(n_rowA * sizeof(double));
		for(i=0;i<n_rowA;i++){
			
			matrix_C[i]=0.0;
			mean=0.0;
			
			for(j=0;j<n_colsA;j++){
				matrix_C[i]=matrix_C[i]+A[i*n_colsA + j];
			}
			mean=mean/n_colsA;
			
			for(j=0; j<n_colsA; ++j){
				tmp += pow(A[i*n_colsA + j] - mean, 2);
			}
			
			matrix_C[i]=sqrt(tmp/n_colsA);
		}
		return matrix_C;
	}
	else{
		printf("axis not valid just 1 or 0 allowed \n");
		exit(1);
		
	}
}

double * sum_array_along(double *A, int axis, int n_rowA, int n_colsA){
	
	static double *matrix_C;
	int i,j;
	
	if(axis==0){
		matrix_C = malloc(n_colsA* sizeof(double));
		for(j=0;j<n_colsA;j++){
			matrix_C[j]=0.0;
			for(i=0;i<n_rowA;i++){
				matrix_C[j]=matrix_C[j]+A[i*n_colsA + j];
			}
		}
		return matrix_C;
		
	}
	else if(axis==1){
		matrix_C = malloc(n_rowA * sizeof(double));
		for(i=0;i<n_rowA;i++){
			matrix_C[i]=0.0;
			for(j=0;j<n_colsA;j++){
				matrix_C[i]=matrix_C[i]+A[i*n_colsA + j];
			}
		}
		return matrix_C;
	}
	else{
		printf("axis not valid just 1 or 0 allowed \n");
		exit(1);
		
	}
}

double *array_multiply(double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB)
{
	int i;
	int size_A,size_B,size_C;
	static double *matrix_C;
	
	size_A=n_rowA*n_colsA;
	size_B=n_rowB*n_colsB;
	
	if(size_A<=size_B){
		matrix_C = malloc(size_A * sizeof(double));
		size_C=size_A;
	}
	else{
		matrix_C = malloc(size_B * sizeof(double));
		size_C=size_B;
	}
	
	for(i=0;i<size_C;i++){
		matrix_C[i]=A[i]*B[i];
	}
	
	return matrix_C;
}

double *array_divide(double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB)
{
	int i;
	int size_A,size_B,size_C;
	static double *matrix_C;
	
	size_A=n_rowA*n_colsA;
	size_B=n_rowB*n_colsB;
	
	if(size_A<=size_B){
		matrix_C = malloc(size_A * sizeof(double));
		size_C=size_A;
	}
	else{
		matrix_C = malloc(size_B * sizeof(double));
		size_C=size_B;
	}
	
	for(i=0;i<size_C;i++){
		matrix_C[i]=A[i]/B[i];
	}
	
	return matrix_C;
}


double * array_sum(double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB){
	
	static double *matrix_C;
	int i,size;
	
	if(n_rowA!=n_rowB ||n_colsA!=n_colsB){
		printf("dimensions do not coincide \n");
		exit(1);
	}
	else{
		
		size=n_rowA*n_colsA;
		matrix_C = malloc(size * sizeof(double));
		for( i=0;i<size;i++){
			matrix_C[i]=A[i]+B[i];
		}
		
	}
	return matrix_C;
	
	
}
double * array_log(double *A, int n_rowA, int n_colsA){
	
	static double *matrix_C;
	int i,size;
	
	size=n_rowA*n_colsA;
	
	matrix_C = malloc(size * sizeof(double));
	for( i=0;i<size;i++){
		matrix_C[i]=log(A[i]);
	}
	return matrix_C;
}

double * array_exp(double *A,int n_rowA, int n_colsA)
{
	static double *matrix_C;
	int i,size;
	
	size=n_rowA*n_colsA;
	
	matrix_C = malloc(size * sizeof(double));
	for( i=0;i<size;i++){
		matrix_C[i]=exp(A[i]);
	}
	return matrix_C;
}

double * array_factrl(int *A, int n_rowA, int n_colsA){
	
	static double *matrix_C;
	int i,size;
	
	size=n_rowA*n_colsA;
	
	matrix_C = malloc(size * sizeof(double));
	for( i=0;i<size;i++ ){
		matrix_C[i]=factrl(A[i]);
	}
	return matrix_C;
}

double * scalar_sum(double a, double *A,int n_rowA, int n_colsA){
	
	static double *matrix_C;
	int i,size;
	
	size=n_rowA*n_colsA;
	
	matrix_C = malloc(size * sizeof(double));
	for( i=0;i<size;i++ ){
		matrix_C[i]=A[i]+a;
	}
	return matrix_C;
}

double * scalar_multiply(double a, double *A,int n_rowA, int n_colsA){
	
	static double *matrix_C;
	int i,size;
	
	size=n_rowA*n_colsA;
	
	matrix_C = malloc(size * sizeof(double));
	for( i=0;i<size;i++ ){
		matrix_C[i]=a*A[i];
	}
	return matrix_C;
}


double * array_normal(long *idum,int n_rowA, int n_colsA){
	static double *matrix_C;
	int i,size;
	double gasdev(long *idum);
	
	size=n_rowA*n_colsA;
	
	matrix_C = malloc(size * sizeof(double));
	for( i=0;i<size;i++ ){
		matrix_C[i]=gasdev(idum);
		
	}
	return matrix_C;
	
}

