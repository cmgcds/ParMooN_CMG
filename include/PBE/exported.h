#ifndef MPIMIS_LAPACKEXPORTED
#define MPIMIS_LAPACKEXPORTED

extern "C" 
{

/********************************************************/
/************************* BLAS *************************/
/********************************************************/

void dgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
           const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
           const double *beta, double *c, const int *ldc);

void dgemv_(const char *trans, const int *m, const int *n, const double *alpha,
           const double *a, const int *lda, const double *x, const int *incx,
           const double *beta, double *y, const int *incy);

void dtrmv_(const char *uplo, const char *transa, const char *diag, const int *n,
           const double *a, const int *lda, double *b, const int *incx);

void   daxpy_(const int *n, const double *alpha, const double *x, const int *incx,
             double *y, const int *incy);

double dnrm2_(const int *n, const double *x, const int *incx);

void   dcopy_(const int *n, const double *x, const int *incx, double *y, const int *incy);

double ddot_(const int *n, const double *x, const int *incx, const double *y, const int *incy);

void   dscal_(const int *n, const double *a, double *x, const int *incx);

void dtrmm_(const char *side, const char *uplo, const char *transa, const char *diag,
           const int *m, const int *n, const double *alpha, const double *a, const int *lda,
           double *b, const int *ldb);

/********************************************************/
/************************ LAPACK ************************/
/********************************************************/

void    dgttrf_(int *n,double *dl,double *d,double *du,double *du2,int *ipiv,int *info);

void    dgttrs_(char* transposed,int *n,int* c,double *dl,double *d,double *du,double *du2,int *ipiv,double* x,int* m,int *info);

void    dgeqrf_(int *m,int *n,double *a,int *lda,double *tau,double *work,int *lwork,int *info);

void    dgesvd_(char *jobu,char *jobvt,int *m,int *n,double *a,int *lda,double *s,double *u,int *ldu,double *vt,int *ldvt,double *work,int *lwork,int *info);

void    dorgqr_(int *m,int *n,int *k,double *a,int *lda,double *tau,double *work,int *lwork,int *info);

void    dgetrf_(int *m,int *n,double *a,int *lda,int *ipiv,int *info);

void    dgetri_(int *n,double *a,int *lda,int *ipiv,double *work,int *lwork,int *info);
}

#endif // MPIMIS_LAPACKEXPORTED
