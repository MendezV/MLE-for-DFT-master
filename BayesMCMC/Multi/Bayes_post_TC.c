#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MIN(x,y) x<y?x:y

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
void write_to_file(double *arr, int size);
double *init_array(int n_steps, double val);
float *load_matrix(char *filename, int *n, int *m);
void print_array(double *x, int n_steps);
void print_array2(int *x, int n_steps);
void print_matrix(double *A,  int n_rowA, int n_colsA);
double gammln(double xx);
double factrl(int n);
double factrl2(double n);
void matrix_multiply(double *matrix_C ,double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
void matrix_multiply2(double *matrix_C ,double *A,int *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
void matrix_multiply3(double *matrix_C ,int *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
double sum_array(double *x, int N);
double std_array(double *x, int N);
void std_array_along(double *matrix_C ,double *A, int axis, int n_rowA, int n_colsA);
void sum_array_along(double *matrix_C ,double *A, int axis, int n_rowA, int n_colsA);
void array_multiply(double *matrix_C ,double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
void array_multiply2(double *matrix_C ,double *A, int *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
void array_divide(double *matrix_C ,double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
void array_sum(double *matrix_C ,double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
void array_log(double *matrix_C ,double *A,int n_rowA, int n_colsA);
void array_exp(double *matrix_C ,double *A,int n_rowA, int n_colsA);
void array_factrl(double *matrix_C ,int *A,int n_rowA, int n_colsA);
void array_factrl2(double *matrix_C ,double *A,int n_rowA, int n_colsA);
void scalar_sum(double *matrix_C ,double a, double *A,int n_rowA, int n_colsA);
void scalar_multiply(double *matrix_C ,double a, double *A,int n_rowA, int n_colsA);
void array_normal(double *matrix_C ,long *idum ,int n_rowA, int n_colsA);
void array_assign(double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
void array_assign2(double *A,float *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
double calculate_likelyhood(double *f,double *V_M, double *V_F,double *onesNbins , double *N_M , double *N_F,double *N_Mfac, double *N_Ffac,int MaxPop, int Nbins, int Tframes, double *histo,double *temp0,double *temp1,double *temp2,double *temp3,double *temp4,double *temp5,double *temp6,double *temp7,double *temp8,double *temp9,double *temp10,double *temp11,double *temp12,double *temp13,double *temp14,double *temp15);


int main(int argc, char **argv){


	int n_steps;
	int i;
	int j;
	int n_row, n_cols;
	
	
	int MaxPop_F, MaxPop_M, MaxPop, Nbins;
	double *N_M;
	double *N_F;
	double *N_Mfac;
	double *N_Ffac;
	int Tframes;
	int size;
	double *histo;
	double *f;
	double *V_M;
	double *V_F;
	double *fprime;
	double *V_Mprime;
	double *V_Fprime;
	double *onesNbins;
	double *onesMaxPopM;
	double *onesMaxPopF;
	double *loglikely;
	double delta_l;
	double alpha;
	float *gv;
	int *temp00,*temp01; /*array settup*/
	double *temp02,*temp03; /*array settup*/
	double *temp0, *temp1,*temp2,*temp3,*temp4,*temp5,*temp6,*temp7,*temp8,*temp9,*temp10,*temp11,*temp12,*temp13,*temp14,*temp15; /*likelihood calculations*/
	double *temp16,*temp17,*temp18,*temp19, *temp20, *temp21,*temp22,*temp23,*temp24; /*monte carlo steps*/
	double *fmle;
	double *V_Mmle;
	double *V_Fmle;
	double maxlike;
	double likeprime, likenow;
	double ff,ffprime;
	double VV_M,VV_Mprime;
	double VV_F,VV_Fprime;
	double sigmaf,sigmaV;
/*variables for use of random number geneator*/
	long *idum;
	long ed=-50;
	idum = &ed;

	Tframes=atoi(argv[2]);
	n_steps=atoi(argv[3]);
	MaxPop_M=atoi(argv[4]);  /*rows in the reshaped histogram */
	MaxPop_F=atoi(argv[5]); /*columns in the reshaped histogram */
	
	/*reading experimental values*/
	gv = load_matrix(argv[1], &n_row, &n_cols);
	Nbins=n_row;
	MaxPop=n_cols;
	size=Nbins*MaxPop;
	
	/*initialize */
	
	loglikely = init_array(n_steps,0.0);
	histo=init_array(size,0.0);
	f=init_array(MaxPop,0.0);
	V_M=init_array(Nbins,0.0);
	V_F=init_array(Nbins,0.0);
	fprime=init_array(MaxPop,0.0);
	V_Mprime=init_array(Nbins,0.0);
	V_Fprime=init_array(Nbins,0.0);
	N_Mfac=init_array(size,0.0);
	N_Ffac=init_array(size,0.0);
	fmle=init_array(MaxPop,0.0);
	V_Mmle=init_array(Nbins,0.0);
	V_Fmle=init_array(Nbins,0.0);
	N_M=init_array(MaxPop,0.0);
	N_F=init_array(MaxPop,0.0);
	
	array_assign2(histo,gv, n_row,n_cols,n_row,n_cols);
	array_normal(f,idum ,MaxPop, 1);
	array_normal(V_M,idum ,Nbins, 1);
	array_normal(V_F,idum ,Nbins, 1);
	
	
	
	/*permanent arrays*/
	onesNbins=init_array(Nbins,1.0);
	onesMaxPopM=init_array(MaxPop_M+1,1.0);
	onesMaxPopF=init_array(MaxPop_F+1,1.0);
	
	temp00= malloc((MaxPop_M+1)*sizeof(int)); /*n_m not sized*/
	for(i=0;i<MaxPop_M+1;i++){ temp00[i]=i; }
	matrix_multiply2( N_M,onesMaxPopF, temp00 ,MaxPop_F+1,1,1,MaxPop_M+1);  /*n_m sized to maxpop, MaxPop_F+1 rows and MaxPop_M+1 cols*/
	
	temp01= malloc((MaxPop_F+1)*sizeof(int)); /*n_f not sized*/
	for(i=0;i<MaxPop_F+1;i++){ temp01[i]=i; }
	matrix_multiply3( N_F,temp01,onesMaxPopM,MaxPop_F+1,1,1,MaxPop_M+1); /*n_f sized to maxpop, MaxPop_F+1 rows and MaxPop_M+1 cols */
	
	
	temp02 = init_array(MaxPop,0.0); /*n_m! sized to Nbins*maxpop declaration*/
	temp03 = init_array(MaxPop,0.0); /*n_f! sized to Nbins*maxpop declaration*/
	array_factrl2(temp02,N_M,MaxPop, 1); /*n_m! sized to maxpop*/
	array_factrl2(temp03,N_F,MaxPop, 1); /*n_m! sized to maxpop */
	
	
	matrix_multiply( N_Mfac,onesNbins,temp02,Nbins,1,1,MaxPop); /*n_m! sized to Nbins*maxpop*/
	matrix_multiply( N_Ffac,onesNbins, temp03,Nbins,1,1,MaxPop); /*n_m! sized to Nbins*maxpop*/
	
	
	/*dummy arrays*/
	temp0 = init_array(size,0.0); /*v_mn_m*/
	temp1 = init_array(size,0.0); /*v_fn_f*/
	temp2 = init_array(size,0.0); /*f sized*/
	temp3 = init_array(size,0.0); /*f+v_mn_m*/
	temp4 = init_array(size,0.0); /*f+v_mn_m+v_fn_f*/
	temp5 = init_array(size,0.0); /*-f-v_mn_m-v_fn_f*/
	temp6 = init_array(size,0.0); /*e(-f-v_mn_m-v_fn_f)*/
	temp7 = init_array(size,0.0); /*e(-f-v_mn_m-v_fn_f))/n_m!*/
	temp8 = init_array(size,0.0); /*e(-f-v_mn_m-v_fn_f))/n_m!n_f!*/
	temp9 = init_array(n_row,0.0); /*z*/
	temp10 = init_array(n_row,0.0); /*logz*/
	temp11 = init_array(n_row,0.0); /*fmeanexp*/
	temp12 = init_array(n_row,0.0); /*N_Mav*/
	temp13 = init_array(n_row,0.0); /*N_Fav*/
	temp14 = init_array(n_row,0.0); /*V_MN_Mav*/
	temp15 = init_array(n_row,0.0); /*V_FN_Fav*/
	temp16 = init_array(n_cols,0.0); /*rand f*/
	temp17 = init_array(n_cols,0.0); /*scaling f*/
	temp18 = init_array(n_row,0.0); /*rand v_M*/
	temp19 = init_array(n_row,0.0); /*scaling V_M*/
	temp20 = init_array(n_row,0.0);  /*rand v_F*/
	temp21 = init_array(n_row,0.0); /*scaling V_F*/
	temp22 = init_array(n_cols,0.0); /* f'f*/
	temp23 = init_array(n_row,0.0); /*V_M'V_M*/
	temp24 = init_array(n_row,0.0); /*V_F'V_F*/
	
	/*first step*/
	loglikely[0]=calculate_likelyhood(f,V_M, V_F, onesNbins , N_M , N_F, N_Mfac, N_Ffac, MaxPop, Nbins, Tframes, histo, temp0, temp1,temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15);
	
	maxlike=loglikely[0];

	/*monte carlo thermalization step*/
	sigmaV=0.03;
	sigmaf=0.03;
	 for(i=1;i<n_steps;i++){
		 
		 f[0]=0.0;
		 f[1]=0.0;
		 f[1+MaxPop_M]=0.0;
		 
		 loglikely[i]=loglikely[i-1];
		 array_normal(temp16,idum ,MaxPop, 1);
		 array_normal(temp18,idum ,Nbins, 1);
		 array_normal(temp20,idum ,Nbins, 1);
		 
		 
		 scalar_multiply(temp17,sigmaf, temp16, MaxPop, 1);
		 scalar_multiply(temp19,sigmaV, temp18, Nbins, 1);
		 scalar_multiply(temp21,sigmaV, temp20, Nbins, 1);
		 
		 array_sum(fprime, temp17,  f, MaxPop, 1, MaxPop, 1);
		 array_sum(V_Mprime, temp19,  V_M, Nbins, 1, Nbins, 1);
		 array_sum(V_Fprime, temp21,  V_F, Nbins, 1, Nbins, 1);
		 
		 fprime[0]=0.0;
		 fprime[1]=0.0;
		 fprime[1+MaxPop_M]=0.0;
		 
		 
		
		  likenow= calculate_likelyhood(f,V_M, V_F, onesNbins , N_M , N_F, N_Mfac, N_Ffac, MaxPop, Nbins, Tframes, histo, temp0, temp1,temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15);
		  likeprime= calculate_likelyhood(fprime,V_Mprime, V_Fprime, onesNbins , N_M , N_F, N_Mfac, N_Ffac, MaxPop, Nbins, Tframes, histo, temp0, temp1,temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15);
		  delta_l=likeprime-likenow;
		 
		 /*
		 array_multiply(temp22, f,  f,  n_cols, 1,  n_cols,1);
		 array_multiply(temp23, V_M,  V_M, n_row, 1, n_row, 1);
		 array_multiply(temp24, V_F,  V_F, n_row, 1, n_row, 1);
		 
		 ff=sum_array(temp22,n_cols);
		 VV_M=sum_array(temp23,n_row);
		 VV_F=sum_array(temp24,n_row);
		 
		 array_multiply(temp22, fprime,  fprime,  n_cols,1,  n_cols,1);
		 array_multiply(temp23, V_Mprime,  V_Mprime, n_row, 1, n_row, 1);
		 array_multiply(temp24, V_Fprime,  V_Fprime, n_row, 1, n_row, 1);
		 
		 ffprime=sum_array(temp22,n_cols);
		 VV_Mprime=sum_array(temp23,n_row);
		 VV_Fprime=sum_array(temp24,n_row);
		 
		 likenow= calculate_likelyhood(f,V_M, V_F, onesNbins , N_M , N_F, N_Mfac, N_Ffac, MaxPop, Nbins, Tframes, histo, temp0, temp1,temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15)-0.5*ff/(sigmaf*sigmaf)-0.5*VV_M/(sigmaV*sigmaV)-0.5*VV_F/(sigmaV*sigmaV);
		 likeprime= calculate_likelyhood(fprime,V_Mprime, V_Fprime, onesNbins , N_M , N_F, N_Mfac, N_Ffac, MaxPop, Nbins, Tframes, histo, temp0, temp1,temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15)-0.5*ffprime/(sigmaf*sigmaf)-0.5*VV_Mprime/(sigmaV*sigmaV)-0.5*VV_Fprime/(sigmaV*sigmaV);
		 delta_l=likeprime-likenow;
		 */
		 
		 
		 if (delta_l>=0.0) {
			 	loglikely[i]= likeprime;
			 	array_assign(f,fprime,MaxPop,1,MaxPop,1);
			 	array_assign(V_M,V_Mprime, Nbins,1,Nbins,1);
			 	array_assign(V_F,V_Fprime, Nbins,1,Nbins,1);
		 	}
		 else{
			 	alpha = drand48();
			 	if (log(alpha)<=delta_l){
						loglikely[i]= likeprime;
						array_assign(f,fprime,MaxPop,1,MaxPop,1);
						array_assign(V_M,V_Mprime, Nbins,1,Nbins,1);
						array_assign(V_F,V_Fprime, Nbins,1,Nbins,1);
						}
		 	}
		 if(maxlike<loglikely[i]){
			 	maxlike=loglikely[i];
			 	array_assign(fmle,f,MaxPop,1,MaxPop,1);
			 	array_assign(V_Mmle,V_M, Nbins,1,Nbins,1);
			 	array_assign(V_Fmle,V_F, Nbins,1,Nbins,1);

		 }
	 }
	
	
	/*printing files with mle parameters*/
	
	FILE *in1;
	char filename1[100]="fmle.dat";
	/*printf("Writing to file: %s\n", filename1);*/
	/*opens the file, writes, closes the file*/
	in1 = fopen(filename1,"w");
	if(!in1){
		printf("problems opening the file %s\n", filename1);
		exit(1);
	}
	
	for(i=0;i<n_cols;i++){
		fprintf(in1, "%f\n", fmle[i]);
	}
	fclose(in1);
	
	
	FILE *in2;
	char filename2[100]="V_Mmle.dat";
	/*printf("Writing to file: %s\n", filename2);*/
	/*opens the file, writes, closes the file*/
	in2 = fopen(filename2,"w");
	if(!in2){
		printf("problems opening the file %s\n", filename2);
		exit(1);
	}
	
	for(i=0;i<n_row;i++){
		fprintf(in2, "%f\n", V_Mmle[i]);
	}
	fclose(in2);
	
	
	FILE *in3;
	char filename3[100]="V_Fmle.dat";
	/*printf("Writing to file: %s\n", filename2);*/
	/*opens the file, writes, closes the file*/
	in3 = fopen(filename3,"w");
	if(!in3){
		printf("problems opening the file %s\n", filename3);
		exit(1);
	}
	
	for(i=0;i<n_row;i++){
		fprintf(in3, "%f\n", V_Fmle[i]);
	}
	fclose(in3);
	
	FILE *in4;
	char filename4[100]="like.dat";
	/*printf("Writing to file: %s\n", filename2);*/
	/*opens the file, writes, closes the file*/
	in4 = fopen(filename4,"w");
	if(!in4){
		printf("problems opening the file %s\n", filename4);
		exit(1);
	}
	
	for(i=0;i<n_steps;i++){
		fprintf(in4, "%f\n", loglikely[i]);
	}
	fclose(in4);
	
	return(0);
	
}



double calculate_likelyhood(double *f,double *V_M, double *V_F,double *onesNbins , double *N_M , double *N_F,double *N_Mfac, double *N_Ffac,int MaxPop, int Nbins, int Tframes, double *histo,double *temp0,double *temp1,double *temp2,double *temp3,double *temp4,double *temp5,double *temp6,double *temp7,double *temp8,double *temp9,double *temp10,double *temp11,double *temp12,double *temp13,double *temp14,double *temp15)
{
	double loglike=0.0;
	int n_rowA, n_colsA, size;
	
	n_rowA=Nbins;
	n_colsA=MaxPop;
	size=n_rowA*n_colsA;
	
	
	matrix_multiply(temp0, V_M, N_M, n_rowA , 1, 1, n_colsA); /*v_mn_m*/
	matrix_multiply(temp1, V_F, N_F, n_rowA , 1, 1, n_colsA); /*v_fn_f*/
	matrix_multiply(temp2, onesNbins, f, n_rowA , 1, 1, n_colsA); /*f sized*/
	array_sum(temp3, temp2,  temp1, n_rowA, n_colsA, n_rowA, n_colsA); /*f+v_mn_m*/
	array_sum(temp4, temp3,  temp0, n_rowA, n_colsA, n_rowA, n_colsA); /*f+v_mn_m+v_fn_f*/
	scalar_multiply(temp5,-1.0, temp4, n_rowA, n_colsA); /*-f-v_mn_m-v_fn_f*/
	array_exp(temp6,temp5, n_rowA, n_colsA); /*e(-f-vn)*/
	array_divide(temp7,temp6 , N_Mfac , n_rowA, n_colsA, n_rowA, n_colsA); /*e(-f-v_mn_m-v_fn_f))/n_m!*/
	array_divide(temp8,temp7 , N_Ffac , n_rowA, n_colsA, n_rowA, n_colsA); /*e(-f-v_mn_m-v_fn_f))/n_f!n_m!*/
	sum_array_along(temp9,temp8, 1, n_rowA, n_colsA); /*z*/
	array_log(temp10,temp9,n_rowA, n_colsA); /*logz*/
	matrix_multiply(temp11, histo, f, n_rowA , n_colsA,n_colsA , 1); /*fmeanexp*/
	matrix_multiply(temp12, histo, N_M, n_rowA , n_colsA,n_colsA , 1); /*N_Mav*/
	matrix_multiply(temp13, histo, N_F, n_rowA , n_colsA,n_colsA , 1); /*N_Fav*/
	array_multiply(temp14,V_M,temp12, n_rowA,1,n_rowA,1); /*V_MN_Mav*/
	array_multiply(temp15,V_F,temp13, n_rowA,1,n_rowA,1); /*V_FN_Fav*/
	
	loglike=-Tframes*( sum_array(temp10, n_rowA )+sum_array(temp11, n_rowA )+sum_array(temp15, n_rowA )+sum_array(temp14, n_rowA ) );
	
	
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

float *load_matrix(char *filename, int *n, int *m){
	float *matrix;
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
	
	matrix = malloc(n_row * n_cols * sizeof(float));
	
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

double factrl2(double n)
/*Returns the value n! as a floating-point number.*/
{
	int m;
	m= (int)n;
	
	double gammln(double xx);
	void nrerror(char error_text[]);
	static int ntop=4;
	static double a[33]={1.0,1.0,2.0,6.0,24.0}; /*Fill in table only as required.*/
	int j;
	if (m < 0) printf("Negative factorial in routine factrl");
	if (m > 32) return exp(gammln(m+1.0));
	/*Larger value than size of table is required. Actually, this big a value is going to overflow
	 on many computers, but no harm in trying.*/
	while (ntop<m) { /*Fill in table up to desired value.*/
		j=ntop++;
		a[ntop]=a[j]*ntop;
	}
	return a[m];
}


void matrix_multiply(double *matrix_C ,double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB)
{
	int i, j,k;
	
	
	for(i=0;i<n_rowA;i++){
		for(j=0;j<n_colsB;j++){
			matrix_C[i*n_colsB + j]=0.0;
			for(k=0;k<n_colsA;k++){
				matrix_C[i*n_colsB + j] = matrix_C[i*n_colsB + j]+A[i*n_colsA + k]*B[k*n_colsB + j];
			}
		}

	}

}

void matrix_multiply2(double *matrix_C ,double *A,int *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB)
{
	int i, j,k;

	for(i=0;i<n_rowA;i++){
		for(j=0;j<n_colsB;j++){
			matrix_C[i*n_colsB + j]=0.0;
			for(k=0;k<n_colsA;k++){
				matrix_C[i*n_colsB + j] = matrix_C[i*n_colsB + j]+A[i*n_colsA + k]*B[k*n_colsB + j];
			}

		}
	
	}
	

}

void matrix_multiply3(double *matrix_C ,int *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB)
{
	int i, j,k;
	
	for(i=0;i<n_rowA;i++){
		for(j=0;j<n_colsB;j++){
			matrix_C[i*n_colsB + j]=0.0;
			for(k=0;k<n_colsA;k++){
				matrix_C[i*n_colsB + j] = matrix_C[i*n_colsB + j]+A[i*n_colsA + k]*B[k*n_colsB + j];
			}
		}
	}
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


void std_array_along(double *matrix_C ,double *A, int axis, int n_rowA, int n_colsA)
{
	double sum = 0.0, mean, tmp = 0.0;
	int i,j;
	
	if(axis==0){
		/*matrix_C = malloc(n_colsA* sizeof(double));*/
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
		
	}
	else if(axis==1){
		/*matrix_C = malloc(n_rowA * sizeof(double));*/
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
	
	}
	else{
		printf("axis not valid just 1 or 0 allowed \n");
		exit(1);
		
	}
}

void  sum_array_along(double *matrix_C ,double *A, int axis, int n_rowA, int n_colsA){
	
	int i,j;
	
	if(axis==0){
		/*matrix_C = malloc(n_colsA* sizeof(double));*/
		for(j=0;j<n_colsA;j++){
			matrix_C[j]=0.0;
			for(i=0;i<n_rowA;i++){
				matrix_C[j]=matrix_C[j]+A[i*n_colsA + j];
			}
		}
		
	}
	else if(axis==1){
		/*matrix_C = malloc(n_rowA * sizeof(double));*/
		for(i=0;i<n_rowA;i++){
			matrix_C[i]=0.0;
			for(j=0;j<n_colsA;j++){
				matrix_C[i]=matrix_C[i]+A[i*n_colsA + j];
			}
		}
	
	}
	else{
		printf("axis not valid just 1 or 0 allowed \n");
		exit(1);
		
	}
}

void array_multiply(double *matrix_C ,double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB)
{
	int i,size;
	
	
	if(n_rowA!=n_rowB ||n_colsA!=n_colsB){
		printf("dimensions do not coincide \n");
		exit(1);
	}
	else{
		
		size=n_rowA*n_colsA;
		for( i=0;i<size;i++){
			matrix_C[i]=A[i]*B[i];
		}
		
	}
	
}

void array_divide(double *matrix_C ,double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB)
{
	int i,size;
	
	
	if(n_rowA!=n_rowB ||n_colsA!=n_colsB){
		printf("dimensions do not coincide \n");
		exit(1);
	}
	else{
		
		size=n_rowA*n_colsA;
		for( i=0;i<size;i++){
			matrix_C[i]=A[i]/B[i];
		}
		
	}
	

}


void array_sum(double *matrix_C ,double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB){
	
	int i,size;
	
	if(n_rowA!=n_rowB || n_colsA!=n_colsB){
		printf("dimensions do not coincide \n");
		exit(1);
	}
	else{
		
		size=n_rowA*n_colsA;
		for( i=0;i<size;i++){
			matrix_C[i]=A[i]+B[i];
		}
		
	}
	
	
}
void array_log(double *matrix_C ,double *A, int n_rowA, int n_colsA){
	
	int i,size;
	
	size=n_rowA*n_colsA;
	
	for( i=0;i<size;i++){
		matrix_C[i]=log(A[i]);
	}

}

void array_exp(double *matrix_C,double *A,int n_rowA, int n_colsA)
{
	int i,size;
	
	size=n_rowA*n_colsA;
	
	for( i=0;i<size;i++){
		matrix_C[i]=exp(A[i]);
	}

}

void array_factrl(double *matrix_C ,int *A, int n_rowA, int n_colsA){
	
	int i,size;
	
	size=n_rowA*n_colsA;
	
	for( i=0;i<size;i++ ){
		matrix_C[i]=factrl(A[i]);
	}

}


void array_factrl2(double *matrix_C ,double *A, int n_rowA, int n_colsA){
	
	int i,size;
	
	size=n_rowA*n_colsA;
	
	for( i=0;i<size;i++ ){
		matrix_C[i]=factrl2(A[i]);
	}
	
}

void scalar_sum(double *matrix_C ,double a, double *A,int n_rowA, int n_colsA){
	
	int i,size;
	
	size=n_rowA*n_colsA;
	
	for( i=0;i<size;i++ ){
		matrix_C[i]=A[i]+a;
	}

}

void scalar_multiply(double *matrix_C ,double a, double *A,int n_rowA, int n_colsA){
	
	int i,size;
	
	size=n_rowA*n_colsA;
	
	for( i=0;i<size;i++ ){
		matrix_C[i]=a*A[i];
	}

}


void array_normal(double *matrix_C ,long *idum,int n_rowA, int n_colsA){
	
	int i,size;
	double gasdev(long *idum);
	
	size=n_rowA*n_colsA;
	
	for( i=0;i<size;i++ ){
		matrix_C[i]=gasdev(idum);
		
	}
	
}

void array_assign(double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB)
{
	int i,size;
	
	if(n_rowA!=n_rowB || n_colsA!=n_colsB){
		printf("dimensions do not coincide \n");
		exit(1);
	}
	else{
		
		size=n_rowA*n_colsA;
		for( i=0;i<size;i++){
			A[i]=B[i];
		}
		
	}
}

void array_assign2(double *A,float *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB)
{
	int i,size;
	
	if(n_rowA!=n_rowB || n_colsA!=n_colsB){
		printf("dimensions do not coincide \n");
		exit(1);
	}
	else{
		
		size=n_rowA*n_colsA;
		for( i=0;i<size;i++){
			A[i]=B[i];
		}
		
	}
}
