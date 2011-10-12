//
//   Copyright (C) 2009 ViGIR Lab (http://vigir.missouri.edu)
//   Written by Guilherme N. DeSouza <desouzag@missouri.edu>
//
//
/*	
 *	Collection of useful and often used functions.
 *
 *	Started Jan 2001, Johnny Park
 */

#include "matrix.h"

/*
 * prints a matrix
 */
void PrintMatrix(matrix *M)
{
  register int i, j;

  for (i=0; i<M->r; i++) {
    for (j=0; j<M->c; j++) {
      printf("%15.4f", M->entry[i][j]);
    }
    printf("\n");
  }
  printf("\n");
  fflush(stdout);
}


/*
 * prints a vector
 */
void PrintVector(vector *p)
{
  register int i;

  for (i=0; i<p->n; i++) {
    printf("%15.4f", p->entry[i]);
  }
  printf("\n");
  fflush(stdout);
}


/*
 * returns determinant of 3x3 matrix
 */
double Determinant33(matrix *A)
{
  return(A->entry[0][0]*(A->entry[1][1]*A->entry[2][2]-A->entry[1][2]*A->entry[2][1]) -
         A->entry[0][1]*(A->entry[1][0]*A->entry[2][2]-A->entry[1][2]*A->entry[2][0]) +
         A->entry[0][2]*(A->entry[1][0]*A->entry[2][1]-A->entry[1][1]*A->entry[2][0])  );
}


/*
 * matrix transpose: M2 = M1'
 */
void Transpose(matrix *M1, matrix *M2)
{
  register int i, j;
  matrix *temp;

  if ((M1->r != M2->c) || (M1->c != M2->r)) {
    printf("ERROR in Transpose: Invalid matrix size\n");
	  exit(-1);
  }

  temp = AllocateMatrix(M1->c, M1->r);
  for (i=0; i<M1->r; i++) for (j=0; j<M1->c; j++) 
    temp->entry[j][i] = M1->entry[i][j];

  CopyMatrix(M2, temp);
  FreeMatrix(temp);
}

/*
 * replaces the entries of matrix with the values in the array v
 */
void SetMatrix(matrix *M, double *v)
{
  register int i, j;

  for (i=0; i<M->r; i++) {
    for (j=0; j<M->c; j++) {
      M->entry[i][j] = v[i*(M->c)+j];
    }
  }
}


/*
 * copy the entryies of the matrix to the array v
 */
void GetMatrix(matrix *M, double *v)
{
  register int i, j;

  for (i=0; i<M->r; i++) {
    for (j=0; j<M->c; j++) {
      v[i*(M->c)+j] = M->entry[i][j];
    }
  }
}


/*
 * replaces the entries of matrix with the values in the array v
 * indexed in row major
 */
void SetMatrixRowMajor(matrix *M, double *v)
{
  register int i, j;

  for (i=0; i<M->c; i++) {
    for (j=0; j<M->r; j++) {
      M->entry[i][j] = v[i*(M->r)+j];
    }
  }
}


/*
 * copy the entryies of the matrix to the array v
 * indexed in row major
 */
void GetMatrixRowMajor(matrix *M, double *v)
{
  register int i, j;

  for (i=0; i<M->c; i++) {
    for (j=0; j<M->r; j++) {
      v[i*(M->r)+j] = M->entry[i][j];
    }
  }
}


/*
 * allocates space for a row by col matrix
 * initial value of 0 assigned for all entries
 * returns pointer to new matrix
 */
matrix *AllocateMatrix(int row, int col)
{
  register int i;
  matrix *M;

  M = (matrix *) malloc (sizeof(matrix));

  M->r = row;
  M->c = col;

  M->entry = (double **) malloc (row * sizeof(double));
  for (i=0; i<row; i++) {
    M->entry[i] = (double *) calloc (col, sizeof(double));
  }
  return (M);
}


/*
 * allocates space for n-dimensional vector
 * returns pointer to new vector
 */
vector *AllocateVector(int num)
{
  vector *p;

  p = (vector *) malloc (sizeof(vector));

  p->n = num;
  p->entry = (double *) calloc (num, sizeof(double));
  return (p);
}


/*
 * frees memory of matrix M
 */
void FreeMatrix(matrix *M)
{
  register int i;

  for (i=0; i<M->r; i++) free(M->entry[i]);
  free(M->entry);
  free(M);
}


/*
 * frees memory of vector v
 */
void FreeVector(vector *v)
{
  free(v->entry);
  free(v);
}


/*
 * copies matrix: M1 = M2 where M2 is the source matrix
 */
void CopyMatrix(matrix *M1, matrix *M2)
{
  register int i, j;

  if ( (M1->r != M2->r) || (M1->c != M2->c) ) {
    printf("ERROR in CopyMatrix: Invalid matrix size\n");
    exit(-1);
  }
  
  for (i=0; i<M2->r; i++) {
    for (j=0; j<M2->c; j++) {
      M1->entry[i][j] = M2->entry[i][j];
    }
  }
}


/*
 * copies vector: v1 = v2 where v2 is the source vector
 */
void CopyVector(vector *v1, vector *v2)
{
  register int i;

  if (v1->n != v2->n) {
    printf("ERROR in CopyVector: Invalid vector size\n");
    exit(-1);
  }
  
  for (i=0; i<v2->n; i++) {
    v1->entry[i] = v2->entry[i];
  }
}

                 
/*
 *  matrix multiplication: M1 x M2 = M3
 */
void MultMatrix(matrix *M1, matrix *M2, matrix *M3)
{
  register int i, j, k;
  matrix *temp;

  if ( (M1->c != M2->r) || (M1->r != M3->r) || (M2->c != M3->c) ) {
    printf("ERROR in MultMatrix: Invalid matrix size\n");
	  exit(-1);
  }

  temp = AllocateMatrix(M3->r, M3->c);
  for (i=0; i<M1->r; i++) 
    for (j=0; j<M2->c; j++) 
      for (k=0; k<M1->c; k++) 
        temp->entry[i][j] += M1->entry[i][k] * M2->entry[k][j];

  CopyMatrix(M3, temp);
  FreeMatrix(temp);
} 


/*
 * matrix addition: M1 + M2 = M3
 */
void AddMatrix(matrix *M1, matrix *M2, matrix *M3)
{
  register int i, j;

  if ((M1->r != M2->r) || (M1->c != M2->c) || 
      (M1->r != M3->r) || (M1->c != M3->c)   ){
    printf("ERROR in AddMatrix: Invalid matrix size\n");
    exit(-1);
  }

  for (i=0; i<M1->r; i++) {
    for (j=0; j<M1->c; j++) {
      M3->entry[i][j] = M1->entry[i][j] + M2->entry[i][j];
    }
  }
}


/*
 * vector addition: p1 + p2 = p3
 */
void AddVector(vector *p1, vector *p2, vector *p3)
{
  register int i;

  if ((p1->n != p2->n) || (p1->n != p3->n)) {
    printf("ERROR in AddVector: Invalid vector size\n");
    exit(-1);
  }

  for (i=0; i<p1->n; i++) {
    p3->entry[i] = p1->entry[i] + p2->entry[i];
  }
}


/*
 * matrix subtration: M1 - M2 = M3
 */
void SubMatrix(matrix *M1, matrix *M2, matrix *M3)
{
  register int i, j;

  if ( (M1->r != M2->r) || (M1->c != M2->c) || 
       (M1->c != M3->c) || (M1->c != M3->c)   ){
    printf("ERROR in SubMatrix: Invalid matrix size\n");
    exit(-1);
  }

  for (i=0; i<M1->r; i++) {
    for (j=0; j<M1->c; j++) {
      M3->entry[i][j] = M1->entry[i][j] - M2->entry[i][j];
    }
  }
}


/*
 * vector subtration: p1 - p2 = p3
 */
void SubVector(vector *p1, vector *p2, vector *p3)
{
  register int i;

  if ((p1->n != p2->n) || (p1->n != p3->n)) {
    printf("ERROR in SubVector: Invalid vector size\n");
    exit(-1);
  }

  for (i=0; i<p1->n; i++) {
    p3->entry[i] = p1->entry[i] - p2->entry[i];
  }
}


/*
 * vector cross product: p1 x p2 = p3
 */
void CrossProduct(vector *p1, vector *p2, vector *p3)
{
  register int i;
  double copy[3];
  
  if ( (p1->n != 3) || (p2->n != 3) || (p3->n != 3) ) {
    printf("ERROR in CrossProduct: Invalid vector size\n");
    exit(-1);
  }
  
  copy[0] = (p1->entry[1]*p2->entry[2]) - (p1->entry[2]*p2->entry[1]);
  copy[1] = (p1->entry[2]*p2->entry[0]) - (p1->entry[0]*p2->entry[2]);
  copy[2] = (p1->entry[0]*p2->entry[1]) - (p1->entry[1]*p2->entry[0]);

  for (i=0; i<3; i++) {
    p3->entry[i] = copy[i];
  }
}


/*
 * vector dot product: p1 . p2 = X
 */
double DotProduct(vector *p1, vector *p2)
{
  register int i;
  double X;

  if (p1->n != p2->n) {
    printf("ERROR in DotProduct: Invalid vector size\n");
    exit(-1);
  }
  
  X = 0;
  for (i=0; i<p1->n; i++) {
    X += (p1->entry[i] * p2->entry[i]);
  }
  return (X);
}


/*
 * computes a unit vector of p1: p2 is a unit vector of p1
 */
void UnitVector(vector *p1, vector *p2)
{
  register int i;
  double mag = 0;

  if (p1->n != p2->n) {
    printf("ERROR in UnitVector: Invalid vector size\n");
    exit(-1);
  }

  for (i=0; i<p1->n; i++) {
    mag += SQ(p1->entry[i]);
  }
  mag = sqrt(mag);
  for (i=0; i<p1->n; i++) {
    p2->entry[i] = p1->entry[i] / mag;
  }
}

/* 
 * multiplies a constant number x to matrix
 * M2 = xM1
 */
void MultConstantToMatrix(double x, matrix *M1, matrix *M2)
{
  register int i, j;

  if ( (M1->r != M2->r) || (M1->c != M2->c) ) {
    printf("ERROR in MultConstantToMatrix: Invalid matrix size\n");
    exit(-1);
  }

  for (i=0; i<M1->r; i++) {
    for (j=0; j<M1->c; j++) {
      M2->entry[i][j] = M1->entry[i][j] * x;
    }
  }
}

/*
 * trace of matrix
 */
double Trace(matrix *M)
{
  register int i;
  double tr = 0;

  if (M->r != M->c) {
    printf("ERROR in Trace: Invalid matrix size\n");
    exit(-1);
  }

  for (i=0; i<M->r; i++) {
    tr += M->entry[i][i];
  }
  return (tr);
}


/*
 * computes rotation matrix R that describes a rotation of theta angle
 * about a unit vector rx, ry, rz
 */
void RotationAboutVector(double theta, double rx, double ry, double rz,
                         matrix *R)
{
  double V,    // 1 - cos(theta)
         C,    // cos(theta)
         S;    // sin(theta)
  double mag;

  if (R->r != 3 || R->c != 3) {
    printf("ERROR in RotationAboutVector: Invalid matrix size\n");
    exit(-1);
  }

  // make sure (rx, ry, rz) is a unit vector
  mag = sqrt(rx*rx + ry*ry + rz*rz);
  rx = rx / mag;
  ry = ry / mag;
  rz = rz / mag;

  // convert theta from degrees to radians
  theta = theta*PI / 180.0;
  C = cos(theta);
  S = sin(theta);
  V = 1.0 - C;

  // compute R
  R->entry[0][0] = rx*rx*V + C;
  R->entry[0][1] = rx*ry*V - rz*S;
  R->entry[0][2] = rx*rz*V + ry*S;
  R->entry[1][0] = rx*ry*V + rz*S;
  R->entry[1][1] = ry*ry*V + C;
  R->entry[1][2] = ry*rz*V - rx*S;
  R->entry[2][0] = rx*rz*V - ry*S;
  R->entry[2][1] = ry*rz*V + rx*S;
  R->entry[2][2] = rz*rz*V + C;
}


/*
 * anti-symmectric matrix
 * M2 = (M1 - M1')
 */
void AntiSymmetric(matrix *M1, matrix *M2)
{
  matrix *M1_t;
  
  if ( (M1->r != M1->c) || (M2->r != M2->c) || (M1->r != M2->c) ) {
    printf("ERROR in AntiSymmetric: Invalid matrix size\n");
    exit(-1);
  }

  M1_t = AllocateMatrix(M1->r, M1->c);
  Transpose(M1, M1_t);
  SubMatrix(M1, M1_t, M2);
  free(M1_t);
}


/*
 * computes rotation quaternion R
 */
void RotationQuaternion(vector *p, matrix *R)
{
  register int i;
  double q[4];
  
  if ( (R->r != 3) || (R->c != 3) ) {
    printf("ERROR in RotationQuaternion: Invalid matrix size\n");
    exit(-1);
  }
  if (p->n != 4) {
    printf("ERROR in RotationQuaternion: Invalid vector size\n");
    exit(-1);
  }

  // copy quaternion values
  for (i=0; i<4; i++) q[i] = p->entry[i];
 
  R->entry[0][0] = SQ(q[0]) + SQ(q[1]) - SQ(q[2]) - SQ(q[3]);
  R->entry[0][1] = 2 * (q[1]*q[2] - q[0]*q[3]);
  R->entry[0][2] = 2 * (q[1]*q[3] + q[0]*q[2]);
  R->entry[1][0] = 2 * (q[1]*q[2] + q[0]*q[3]);
  R->entry[1][1] = SQ(q[0]) + SQ(q[2]) - SQ(q[1]) - SQ(q[3]);
  R->entry[1][2] = 2 * (q[2]*q[3] - q[0]*q[1]);
  R->entry[2][0] = 2 * (q[1]*q[3] - q[0]*q[2]);
  R->entry[2][1] = 2 * (q[2]*q[3] + q[0]*q[1]);
  R->entry[2][2] = SQ(q[0]) + SQ(q[3]) - SQ(q[1]) - SQ(q[2]);
}

void QuaternionFromRotation(matrix *R, vector *p)
{
  register int i, j;
  double r[3][3];
  double t, s;
  
  if ( (R->r != 3) || (R->c != 3) ) {
    printf("ERROR in QuaternionFromRotation: Invalid matrix size\n");
    exit(-1);
  }
  if (p->n != 4) {
    printf("ERROR in QuaternionFromRotation: Invalid vector size\n");
    exit(-1);
  }

  // copy matrix
  for (i=0; i<3; i++) for (j=0; j<3; j++) {
    r[i][j] = R->entry[i][j];
  }

  // compute trace
  t = 1 + r[0][0] + r[1][1] + r[2][2];

  if (t > 0.00000001) {
    s = sqrt(t) * 2.0;
    p->entry[1] = (r[2][1] - r[1][2]) / s;
    p->entry[2] = (r[0][2] - r[2][0]) / s;
    p->entry[3] = (r[1][0] - r[0][1]) / s;
    p->entry[0] = 0.25 * s;
  } 
  else {
    if (r[0][0] > r[1][1] && r[0][0] > r[2][2]) {
      s = sqrt(1.0 + r[0][0] - r[1][1] - r[2][2]) * 2.0;
      p->entry[1] = 0.25 * s;
      p->entry[2] = (r[1][0] + r[0][1]) / s;
      p->entry[3] = (r[0][2] + r[2][0]) / s;
      p->entry[0] = (r[2][1] - r[1][2]) / s;
    }
    else if (r[1][1] > r[2][2]) {
      s = sqrt(1.0 - r[0][0] + r[1][1] - r[2][2]) * 2.0;
      p->entry[1] = (r[1][0] + r[0][1]) / s;
      p->entry[2] = 0.25 * s;
      p->entry[3] = (r[2][1] + r[1][2]) / s;
      p->entry[0] = (r[0][2] - r[2][0]) / s;
    }
    else {
      s = sqrt(1.0 - r[0][0] - r[1][1] + r[2][2]) * 2.0;
      p->entry[1] = (r[0][2] + r[2][0]) / s;
      p->entry[2] = (r[2][1] + r[1][2]) / s;
      p->entry[3] = 0.25 * s;
      p->entry[0] = (r[1][0] - r[0][1]) / s;
    }
  }
}


/*
 * computes eigenvalues and eigenvectors of a symmetric matrix A
 * using LAPACK dsyev_
 * 
 *    A: symmeric matrix
 *    Q: orthogonal matrix of eigenvectors
 *    v: eigenvalues
 *    
 *    
 * dsyev_(char *jobz, char *uplo, integer *n, doublereal *a,
 *        integer *lda, doublereal *w, doublereal *work, integer *lwork,
 *        integer *info);
 * 
 * 
 *  JOBZ    (input) CHARACTER*1
 *          = 'N':  Compute eigenvalues only;
 *          = 'V':  Compute eigenvalues and eigenvectors.
 *
 *  UPLO    (input) CHARACTER*1
 *          = 'U':  Upper triangle of A is stored;
 *          = 'L':  Lower triangle of A is stored.
 *
 *  N       (input) INTEGER
 *          The order of the matrix A.  N >= 0.
 *
 *  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
 *          On entry, the symmetric matrix A.  If UPLO = 'U', the
 *          leading N-by-N upper triangular part of A contains the
 *          upper triangular part of the matrix A.  If UPLO = 'L',
 *          the leading N-by-N lower triangular part of A contains
 *          the lower triangular part of the matrix A.
 *          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
 *          orthonormal eigenvectors of the matrix A.
 *          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
 *          or the upper triangle (if UPLO='U') of A, including the
 *          diagonal, is destroyed.
 *
 *  LDA     (input) INTEGER
 *          The leading dimension of the array A.  LDA >= max(1,N).
 *
 *  W       (output) DOUBLE PRECISION array, dimension (N)
 *          If INFO = 0, the eigenvalues in ascending order.
 *
 *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
 *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
 *
 *  LWORK   (input) INTEGER
 *          The length of the array WORK.  LWORK >= max(1,3*N-1).
 *          For optimal efficiency, LWORK >= (NB+2)*N,
 *          where NB is the blocksize for DSYTRD returned by ILAENV.
 *
 *          If LWORK = -1, then a workspace query is assumed; the routine
 *          only calculates the optimal size of the WORK array, returns
 *          this value as the first entry of the WORK array, and no error
 *          message related to LWORK is issued by XERBLA.
 *
 *  INFO    (output) INTEGER
 *          = 0:  successful exit
 *          < 0:  if INFO = -i, the i-th argument had an illegal value
 *          > 0:  if INFO = i, the algorithm failed to converge; i
 *                off-diagonal elements of an intermediate tridiagonal
 *                form did not converge to zero.
 */
int SymmetricEigens(matrix *A, matrix *Q, vector *v)
{
  register int i, j;
  char     jobz    = 'V';
  char     uplo    = 'U';
  integer  n;
  double   a[100];
  integer  lda;
  double   w[10];
  double   work[50];
  integer  lwork = 50;
  integer  info;
  
  
  // check input
  if ( (A->r != A->c) || (Q->r != Q->c) || 
       (A->r != Q->r) || (A->r != v->n) ) {
    fprintf(stderr, "ERROR in SymmetricEigens: Invalid matrix size\n");
    exit(1);
  }
  if (A->r > 10) {
    fprintf(stderr, "ERROR in SymmetricEigens: Above maximum size 10\n");
    exit(1);
  }

  // size of matrix
  n = A->r;
  lda = A->r;

  // copy A
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      a[i*n+j] = A->entry[j][i];  // switch to column major index
    }
  }

  // compute eigens
  dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);
    

  // copy results
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      Q->entry[j][i] = a[i*3+j];
    }
    v->entry[i] = w[i];
  }
  
  return info;
}
  
/*
 * Computes largest eigenvalue and corresponding eigenvector of A
 * using CLAPACK dsyevx_
 * 
 *    A: symmeric matrix
 *    evec: largest eigenvector
 *    *eval: largest eigenvalue 
 */
void SymmetricLargestEigens(matrix *A, vector *evec, double *eval)
{  
  register int i, j;
  char     jobz    = 'V';
  char     range   = 'I';
  char     uplo    = 'U';
  integer  n;
  double   a[100];
  integer  lda;
  double   vl, vu;
  integer  il, iu;
  double   abstol  = 0.000001;
  integer  m;
  double   w[2];
  double   z[10];
  integer  ldz;
  double   work[200];
  integer  lwork = 200;
  integer  iwork[50];
  integer  ifail[10];
  integer  info;

    // check input
  if ( (A->r != A->c) || (A->r != evec->n) ) {
    fprintf(stderr, "ERROR in SymmetricEigens: Invalid matrix size\n");
    exit(1);
  }
  if (A->r > 10) {
    fprintf(stderr, "ERROR in SymmetricEigens: Above maximum size 10\n");
    exit(1);
  }

  // size of matrix
  n = lda = ldz = il = iu = A->r;

  // copy A
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      a[i*n+j] = A->entry[j][i];  // switch to column major index
    }
  }

  // compute eigens
  dsyevx_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu,
          &abstol, &m, w, z, &ldz, work, &lwork, iwork, ifail, &info); 

  // copy results
  for (i=0; i<n; i++) {
    evec->entry[i] = z[i];
  }
  *eval = w[0];
}


/*
 * Computes smallest eigenvalue and corresponding eigenvector of A
 * using CLAPACK dsyevx_
 * 
 *    A: symmeric matrix
 *    evec: smallest eigenvector
 *    *eval: smallest eigenvalue 
 */
void SymmetricSmallestEigens(matrix *A, vector *evec, double *eval)
{
  register int i, j;
  char     jobz    = 'V';
  char     range   = 'I';
  char     uplo    = 'U';
  integer  n;
  double   a[100];
  integer  lda;
  double   vl, vu;
  integer  il, iu;
  double   abstol  = 0.000001;
  integer  m;
  double   w[2];
  double   z[10];
  integer  ldz;
  double   work[200];
  integer  lwork = 200;
  integer  iwork[50];
  integer  ifail[10];
  integer  info;

    // check input
  if ( (A->r != A->c) || (A->r != evec->n) ) {
    fprintf(stderr, "ERROR in SymmetricEigens: Invalid matrix size\n");
    exit(1);
  }
  if (A->r > 10) {
    fprintf(stderr, "ERROR in SymmetricEigens: Above maximum size 10\n");
    exit(1);
  }

  // size of matrix
  n = lda = ldz = A->r;

  il = iu = 1;
   
  // copy A
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      a[i*n+j] = A->entry[j][i];  // switch to column major index
    }
  }

  // compute eigens
  dsyevx_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu,
          &abstol, &m, w, z, &ldz, work, &lwork, iwork, ifail, &info); 

  // copy results
  for (i=0; i<n; i++) {
    evec->entry[i] = z[i];
  }
  *eval = w[0];
}


/*
*     SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGESV computes the solution to a real system of linear equations
*     A * X = B,
*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
*
*  The LU decomposition with partial pivoting and row interchanges is
*  used to factor A as
*     A = P * L * U,
*  where P is a permutation matrix, L is unit lower triangular, and U is
*  upper triangular.  The factored form of A is then used to solve the
*  system of equations A * X = B.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of linear equations, i.e., the order of the
*          matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the N-by-N coefficient matrix A.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (output) INTEGER array, dimension (N)
*          The pivot indices that define the permutation matrix P;
*          row i of the matrix was interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the N-by-NRHS matrix of right hand side matrix B.
*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
*                has been completed, but the factor U is exactly
*                singular, so the solution could not be computed.
*/


/*
 * Computes inverse of a square matrix M1
 */
void InverseMatrix(matrix *M1, matrix *M2)
{
  int     i, j;
  integer N;      // size of matrix
  double  *A,     // input matrix
          *B;     // N by N identity matrix
  integer *IPIV;  // pivot indices
  integer INFO;

  // must be a square matrix
  if ( M1->r != M1->c  || M2->r != M2->c || M1->r != M2->r){
    fprintf(stderr, "ERROR in InverseMatrix: Invalid matrix size\n");
    exit(-1);
  }

  N = M1->r;

  // allocate A, B, IPIV
  A = (double *) malloc (N * N * sizeof(double));
  B = (double *) calloc (N * N,  sizeof(double));
  IPIV = (integer *) malloc (N * sizeof(integer));


  // A = M, B = identity matrix
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      A[i*N+j] = M1->entry[j][i];  // switch to column major index
    }
    B[i*N+i] = 1;
  }
  
  // compute inverse
  dgesv_(&N, &N, A, &N, IPIV, B, &N, &INFO);

  if (INFO != 0) {
    fprintf(stderr, "ERROR in InverseMatrix: Inverse cannot be computed\n");
    exit(-1);
  }
/*  
  // copy inverse
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      M2->entry[i][j] = B[j*N+i]; // B is in column major index
    }
  }
*/
  // free memory
  free(A); free(B); free(IPIV);
}
 
