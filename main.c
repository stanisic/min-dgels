#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>

// From src/additional
#include <f2c.h>
#include <multiple_regression.h>

void compute(double *mA, double *mB, double *mX, long mm, int nn)
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

  if(mm <= nn) {
    printf("\nERROR: This function is not intended for the use when number of parameters is larger than the number of observations. Check how your matrices A and B were allocated or simply add more benchmarks.\n");
    exit(1);
  }
  
  char trans = 'N'; // the linear system involves A
  integer m = mm; // the number of rows of the matrix A
  integer n = nn; // the number of columns of the matrix A
  integer nrhs = 1; // number of columns of B and X (which are vectors therefore nrhs=1)
  doublereal *A = malloc(sizeof(double)*n*m); // input/output array, dimension (LDA,N)
  doublereal *B = malloc(sizeof(double)*m); // input/output array, dimension (LDB,NRHS), on exit overwritten by the solution where rows 1 to n of B contain the least squares solution vector
  
  for(int i=0; i < m; i++)
  {
     B[i] = mB[i];
     for(int j=0; j < n; j++)
       A[j*m+i] = mA[i*n+j];
  }

  integer lda = m; // The leading dimension of the array A
  integer ldb = m; // The leading dimension of the array B
  integer info; // information on the output
  
  integer lwork = n+n; // min(M,N) + max(min(M,N), NRHS)
  doublereal *work = malloc(sizeof(double)*lwork); // (output)

  // Running CLAPACK
  dgels_(&trans, &m, &n, &nrhs, A, &lda, B, &ldb, work, &lwork, &info);

  // Check for the output
  if(info != 0) {
    printf("Problems with DGELS; info=%ld\n", info);
    exit(2);
  }
	
  for(int j=0; j < n; j++)
    mX[j] = B[j];

  free(A);
  free(B);
  free(work);
}

void validate()
{
  printf("\nValidation is not implemented yet!\n");
}

int main()
{
  // Simple example definition
  unsigned m=3;
  unsigned n=2;
  double mA[] = {1, -1, 1, 1, 2, 1};
  double mB[] = {2, 4, 8};
  double *mX = (double *) malloc(n*sizeof(double));

  // Printing input of the problem
  for (int i = 0; i < m; i++)
  {
    printf("x%d*%d ", 0, (int)mA[i*n]);
    for (int j = 1; j < n; j++)
      printf("+ x%d*%d ", j, (int)mA[i*n+j]);
    printf("= %g\n", mB[i]);
  }

  // Main function to compute coefficients
  compute(mA, mB, mX, m, n);

  // Printing output
  printf("\nCoefficients:");
  printf("  %g", mX[0]);
  for (int i = 1; i < n; i++)
      printf(", %g", mX[i]);
  
  for (int i = 0; i < m; i++)
  {
    printf("\n%g*%d ", mX[0], (int)mA[i*n]);
    for (int j = 1; j < n; j++)
      printf("+ %g*%d ", mX[j], (int)mA[i*n+j]);
    printf("= %g", mB[i]);
  }

  // Validating the model accuracy
  validate();
  
  return 0;
}
