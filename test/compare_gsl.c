#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>

// From additional-lapack-func
#include <f2c.h>
#include <multiple_regression.h>

// Necessary for GSL
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multifit.h>

static long count_file_lines(FILE *f)
{
	int ch, lines=0;
	while(!feof(f))
	{
	  ch = fgetc(f);
	  if(ch == '\n')
	    {
	      lines++;
	    }
	}
	rewind(f);

	return lines;
}

static void read_input_file(double *mA, double *mB, long m, unsigned n, FILE *f)
{
	char buffer[1024];
	char *record,*line;
	int i=0,j=0;

	line=fgets(buffer,sizeof(buffer),f);//skipping first line
	while((line=fgets(buffer,sizeof(buffer),f))!=NULL && i<m)
	{
		record = strtok(line,",");
		mB[i] = atof(record);		
		record = strtok(NULL,",");
		j=0;
		while(record != NULL)
		{
			mA[i*n+j] = atof(record) ;
			++j;
			record = strtok(NULL,",");

		}
		i++;
	}

}

void lapack_multiple_reg_coeff(double *mA, double *mB, double *mX, long mm, int nn)
{

  /*  Arguments */
/*  ========= */

/*  TRANS   (input) CHARACTER*1 */
/*          = 'N': the linear system involves A; */
/*          = 'T': the linear system involves A**T. */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix A.  M >= 0. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix A.  N >= 0. */

/*  NRHS    (input) INTEGER */
/*          The number of right hand sides, i.e., the number of */
/*          columns of the matrices B and X. NRHS >=0. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          On entry, the M-by-N matrix A. */
/*          On exit, */
/*            if M >= N, A is overwritten by details of its QR */
/*                       factorization as returned by DGEQRF; */
/*            if M <  N, A is overwritten by details of its LQ */
/*                       factorization as returned by DGELQF. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(1,M). */

/*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS) */
/*          On entry, the matrix B of right hand side vectors, stored */
/*          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS */
/*          if TRANS = 'T'. */
/*          On exit, if INFO = 0, B is overwritten by the solution */
/*          vectors, stored columnwise: */
/*          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least */
/*          squares solution vectors; the residual sum of squares for the */
/*          solution in each column is given by the sum of squares of */
/*          elements N+1 to M in that column; */
/*          if TRANS = 'N' and m < n, rows 1 to N of B contain the */
/*          minimum norm solution vectors; */
/*          if TRANS = 'T' and m >= n, rows 1 to M of B contain the */
/*          minimum norm solution vectors; */
/*          if TRANS = 'T' and m < n, rows 1 to M of B contain the */
/*          least squares solution vectors; the residual sum of squares */
/*          for the solution in each column is given by the sum of */
/*          squares of elements M+1 to N in that column. */

/*  LDB     (input) INTEGER */
/*          The leading dimension of the array B. LDB >= MAX(1,M,N). */

/*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

/*  LWORK   (input) INTEGER */
/*          The dimension of the array WORK. */
/*          LWORK >= max( 1, MN + max( MN, NRHS ) ). */
/*          For optimal performance, */
/*          LWORK >= max( 1, MN + max( MN, NRHS )*NB ). */
/*          where MN = min(M,N) and NB is the optimum block size. */

/*          If LWORK = -1, then a workspace query is assumed; the routine */
/*          only calculates the optimal size of the WORK array, returns */
/*          this value as the first entry of the WORK array, and no error */
/*          message related to LWORK is issued by XERBLA. */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */
/*          > 0:  if INFO =  i, the i-th diagonal element of the */
/*                triangular factor of A is zero, so that A does not have */
/*                full rank; the least squares solution could not be */
/*                computed. */

/*  ===================================================================== */

  char trans = 'N';
  integer m = mm;
  integer n = nn;
  integer nrhs = 1; // number of columns of B and X (wich are vectors therefore nrhs=1)
  doublereal *A = malloc(sizeof(double)*n*m); // (/!\ modified at the output) contain the model and the different values of pararmters
  doublereal *B = malloc(sizeof(double)*m);
  
  /* for (int i=0; i < m; i++) */
  /* { */
  /*       B[i] = mB[i]; */
  /* 	for (int j=0; j < n; j++) */
  /* 	    A[i*n+j] = mA[i*n+j]; */
  /* } */

  for (int i=0; i < m; i++)
  {
        B[i] = mB[i];
	for (int j=0; j < n; j++)
	    A[j*m+i] = mA[i*n+j];
  }

  integer lda = m; 
  integer ldb = m; //
  integer info;
  
  integer lwork = n*2;
  doublereal *work = malloc(sizeof(double)*lwork); // (output)

  /* // Running CLAPACK */
  dgels_(&trans, &m, &n, &nrhs, A, &lda, B, &ldb, work, &lwork, &info);

  /* Check for the full rank */
  if( info != 0 ) {
    printf( "Problems with DGELS; info=%i\n");
    exit( 1 );
  }
	
  for (int j=0; j < n; j++)
	mX[j] = B[j];

  free(A);
  free(B);
  free(work);
}

void gsl_multiple_reg_coeff(double *mA, double *mB, double *mX, long m, int n)
{
	gsl_matrix *A = gsl_matrix_calloc(m, n);
	gsl_vector *B = gsl_vector_alloc(m);
	gsl_vector *X = gsl_vector_alloc(n);

	for (int i = 0; i < m; i++)
	{
	  gsl_vector_set(B, i, mB[i]);
          for (int j = 0; j < n; j++)
	    {
	      gsl_matrix_set(A, i, j, mA[i*n+j]);
	    }
	}
		
	double chisq;
	gsl_matrix *cov = gsl_matrix_alloc(n, n);
	gsl_multifit_linear_workspace * wspc = gsl_multifit_linear_alloc(m, n);
	gsl_multifit_linear(A, B, X, cov, &chisq, wspc);

	for(int i=0; i<n; i++)
	  mX[i] = gsl_vector_get(X, i);

	gsl_matrix_free(A);
	gsl_matrix_free(cov);
	gsl_vector_free(B);
	gsl_vector_free(X);
	gsl_multifit_linear_free(wspc);
}

void print_solution(double *mA, double *mB, double *mX, long m, int n)
{
   int i=0;
   int j=0;
   
   // Coefficients
   printf("\nCoefficients:");
   printf("  %g", mX[i]);
   for (i = 1; i < n; i++)
      printf(", %g", mX[i]);

   double sumR=0.;
   double sumB=0.;
   double ri=0.;
   
   for (int i = 0; i < m; i++)
   {
     sumB += mB[i]*mB[i];
     ri=0;
     for (int j = 0; j < n; j++)
       ri += mA[i*n+j]*mX[j];
     ri = ri - mB[i];
     sumR += ri*ri;
   }
   printf("\nError ||r||2 / ||b||2: %g\n", sqrt(sumR)/sqrt(sumB));
}

// Inspired from http://s-mat-pcs.oulu.fi/~mpa/matreng/eem5_5-1.htm
int first_example()
{
        unsigned m=3;
	unsigned n=2;

  	double mA[] = {1, -1, 1, 1, 2, 1};
  	double mB[] = {2, 4, 8};
	double *mX1 = (double *) malloc(n*sizeof(double));
	double *mX2 = (double *) malloc(n*sizeof(double));

	// LAPACK
	lapack_multiple_reg_coeff(mA, mB, mX1, m, n);
		  
	// GSL       
        gsl_multiple_reg_coeff(mA, mB, mX2, m, n); 
	  
	printf("##################\n");
        printf("\tLAPACK");
	print_solution(mA, mB, mX1, m, n);
	for (int i = 0; i < m; i++)
	{
	  printf("%g*%d ", mX1[0], (int)mA[i*n]);
	  for (int j = 1; j < n; j++)
              printf("+ %g*%d ", mX1[j], (int)mA[i*n+j]);
  	  printf("= %g\n", mB[i]);
        }
        printf("\tGSL");
	print_solution(mA, mB, mX2, m, n);
	for (int i = 0; i < m; i++)
	{
	  printf("%g*%d ", mX2[0],(int)mA[i*n]);
	  for (int j = 1; j < n; j++)
              printf("+ %g*%d ", mX2[j], (int)mA[i*n+j]);
  	  printf("= %g\n", mB[i]);
        }
	printf("##################\n");
	
	free(mX1);
	free(mX2);
	
	return 0;
}

double w[] = {	52.21, 53.12, 54.48, 55.84, 57.20,
		58.57, 59.93, 61.29, 63.11, 64.47,
		66.28, 68.10, 69.92, 72.19, 74.46 };
double h[] = {	1.47, 1.50, 1.52, 1.55, 1.57,
		1.60, 1.63, 1.65, 1.68, 1.70,
		1.73, 1.75, 1.78, 1.80, 1.83	};

int second_example()
{
        unsigned m=15;
	unsigned n=3;

  	double *mA = (double *) malloc(m*n*sizeof(double));
	double *mB = (double *) malloc(m*sizeof(double));
	double *mX1 = (double *) malloc(n*sizeof(double));
	double *mX2 = (double *) malloc(n*sizeof(double));

	for (int i = 0; i < m; i++)
	{
	  mB[i] = w[i];
  	  mA[i*n] = 1;
    	  mA[i*n+1] = h[i];
    	  mA[i*n+2] = h[i]*h[i];
	}
	

	// LAPACK
	lapack_multiple_reg_coeff(mA, mB, mX1, m, n);
		  
	// GSL       
        gsl_multiple_reg_coeff(mA, mB, mX2, m, n); 
	  
	printf("##################\n");
        printf("\tLAPACK");
	print_solution(mA, mB, mX1, m, n);
        printf("\tGSL");
	print_solution(mA, mB, mX2, m, n);
	printf("##################\n");

	free(mA);
	free(mB);	
	free(mX1);
	free(mX2);
	
	return 0;
}

int third_example()
{
        char *filepath = "input/input.csv";
	FILE *f;
	f = fopen(filepath, "a+");	
	
        unsigned m;
	m = count_file_lines(f);
	unsigned n=3;


  	double *mA = (double *) malloc(m*n*sizeof(double));
	double *mB = (double *) malloc(m*sizeof(double));
	double *mX1 = (double *) malloc(n*sizeof(double));
	double *mX2 = (double *) malloc(n*sizeof(double));
	
        read_input_file(mA, mB, m, n, f);

	// LAPACK
	lapack_multiple_reg_coeff(mA, mB, mX1, m, n);
		  
	// GSL       
        gsl_multiple_reg_coeff(mA, mB, mX2, m, n); 
	  
	printf("##################\n");
        printf("\tLAPACK");
	print_solution(mA, mB, mX1, m, n);
        printf("\tGSL");
	print_solution(mA, mB, mX2, m, n);
	printf("##################\n");

	free(mA);
	free(mB);	
	free(mX1);
	free(mX2);
	
	return 0;
}

void main()
{
  printf("\n\n\nMinimal example 3x2:\n");
  first_example();
  printf("\n\n\nSmall example 15x3:\n");
  second_example();
  printf("\n\n\nMedium example 5411x3 with input file:\n");
  third_example();
  return 0;
}
