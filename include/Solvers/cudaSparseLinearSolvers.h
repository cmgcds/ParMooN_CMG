#ifndef __CUDA_REFACTOR_THIV__
#define __CUDA_REFACTOR_THIV__

#include <cuda_runtime.h>
#include "cusolverSp.h"
#include "cusolverRf.h"
#include "cusolverSp_LOWLEVEL_PREVIEW.h"




class cudaRefactor 
{
    public:
    
    //Default COnstructor 
    cudaRefactor()
    {

    };
    // ---------- Variable Declaration -------------- //

    // ----------- CuSparse  and Cusolver  handler Declararions
    cusolverRfHandle_t cusolverRfH = NULL; // refactorization
    cusolverSpHandle_t cusolverSpH = NULL; // reordering, permutation and 1st LU factorization
    cusparseHandle_t   cusparseH = NULL;   // residual evaluation
    cudaStream_t stream = NULL;
    cudaStream_t stream2 = NULL;

    cusparseMatDescr_t descrA = NULL; // A is a base-0 general matrix

    csrluInfoHost_t info = NULL; // opaque info structure for LU with parital pivoting

    int rowsA = 0; // number of rows of A
    int colsA = 0; // number of columns of A
    int nnzA  = 0; // number of nonzeros of A
    int baseA = 0; // base index in CSR format

    // ----- Permuted Matrix Arrays and vectors ----- //
    int* h_csrColIndB = NULL;
    int* h_csrRowPtrB = NULL;
    double* h_csrValB = NULL;
    double* h_x = NULL;
    double* h_xhat = NULL;
    double* h_bhat = NULL;
    int nRow;
    int nNNZ;

    // ------ Device Arrays ------------ //
    int* d_csrRowPtrA   = NULL;
    int* d_csrColIndA   = NULL;
    double* d_csrValA   = NULL;
    double* d_x         = NULL;
    double* d_b         = NULL;
    int* d_P            = NULL;
    int* d_Q            = NULL;
    double* d_T         = NULL;

    // --- Permutation Matrix Declaration ------ //
    int* h_Qreorder = NULL;
    int* h_mapBfromA = NULL;
    int *h_Plu = NULL; // <int> n
    int *h_Qlu = NULL; // <int> n

    int nnzL = 0;
    int *h_csrRowPtrL = NULL; // <int> n+1
    int *h_csrColIndL = NULL; // <int> nnzL
    double *h_csrValL = NULL; // <double> nnzL

    int nnzU = 0;
    int *h_csrRowPtrU = NULL; // <int> n+1
    int *h_csrColIndU = NULL; // <int> nnzU
    double *h_csrValU = NULL; // <double> nnzU

    int *h_P = NULL; // <int> n, P = Plu * Qreorder
    int *h_Q = NULL; // <int> n, Q = Qlu * Qreorder

    size_t size_perm = 0;
    size_t size_internal = 0; 
    size_t size_lu  = 0; 

    void *buffer_cpu = NULL;


    // TIME VALUES 
    double start, stop;
    double time_reorder;
    double time_perm;
    double time_sp_analysis;
    double time_sp_factor;
    double time_sp_solve;
    double time_sp_extract;
    double time_rf_assemble;
    double time_rf_reset;
    double time_rf_refactor;
    double time_rf_solve;

    //-- Values used in Solvers -----------//
        // the constants used in cusolverRf
    // nzero is the value below which zero pivot is flagged.
    // nboost is the value which is substitured for zero pivot.
    double nzero = 0.0;
    double nboost= 0.0;
    // the constant used in cusolverSp
    // singularity is -1 if A is invertible under tol
    // tol determines the condition of singularity
    // pivot_threshold decides pivoting strategy            
    int singularity = 0; 
    const double tol = 1.e-14;
    const double pivot_threshold = 1.0;
    // the constants used in cusolverRf
    const cusolverRfFactorization_t fact_alg = CUSOLVERRF_FACTORIZATION_ALG0; // default
    const cusolverRfTriangularSolve_t solve_alg = CUSOLVERRF_TRIANGULAR_SOLVE_ALG1; // default



    // Initialise the cuda handles for cusparse and cu solver  and cusolver RF routines
    void initialiseCudaHandles();


    // Performs Reordering and computes the LU Factorisation and stores the L and U Factors with their Permutations
    void LU_DecompositionHost(int* h_rowPtrA, int* h_colPtrB, 
                                        double* h_valuesA, int rowsA,int nnzA,
                                        double* h_b ,const char* reorder);

    // Intialise the methods needed methoed to initialise RF routines on the Device
    void cudaRefactorize(int* h_csrRowPtrA, int* h_csrColIndA, double* h_csrValA, int rowsA, int N);

    // Refactorise the system using cusolverRF
    void cudaRefactorSolve(double* h_b);

    //Solve the system using the L and U values obtained using refactorisation
    void cudaSolve();

    // Destroy Cuda handles
    void destroyCudahandles();

    //Reseting the reFactor SOlver Routines in order to do LU Factorisation Agian
    void resetCudaRF();



};





class cudaLowLevelQR 
{

    public :

        //Default COnstructor 
        cudaLowLevelQR()
        {

        };
        cusolverSpHandle_t cusolverSpH = NULL; // reordering, permutation and 1st LU factorization
        cusparseHandle_t   cusparseH = NULL;   // residual evaluation
        cudaStream_t stream = NULL;
        cusparseMatDescr_t descrA = NULL; // A is a base-0 general matrix

        csrqrInfoHost_t h_info = NULL; // opaque info structure for LU with parital pivoting
        csrqrInfo_t d_info = NULL; // opaque info structure for LU with parital pivoting

        int rowsA = 0; // number of rows of A
        int colsA = 0; // number of columns of A
        int nnzA  = 0; // number of nonzeros of A
        int baseA = 0; // base index in CSR format


        double *h_x = NULL; // <double> n,  x = A \ b
        
           // ------------------ Added for Permutation Details ----------------------------- //
   
    double *h_Qb = NULL; /* Q*b */
    double *h_z  = NULL; /* z = B \ (Q*b) */

    int *h_Q = NULL; /* <int> n */
                     /* reorder to reduce zero fill-in */
                     /* Q = symrcm(A) or Q = symamd(A) */
    /* B = Q*A*Q' or B = A(Q,Q) by MATLAB notation */
    int *h_csrRowPtrB = NULL; /* <int> n+1 */
    int *h_csrColIndB = NULL; /* <int> nnzA */
    double *h_csrValB = NULL; /* <double> nnzA */
    int *h_mapBfromA = NULL;  /* <int> nnzA */

    // ------------- END -- Added for Permutation Details ----------------------------- //




    size_t size_internal = 0; 
    size_t size_chol = 0; // size of working space for csrlu
    void *buffer_cpu = NULL; // working space for Cholesky
    void *buffer_gpu = NULL; // working space for Cholesky

    // ============= Allocatino in DEvice ================ //
        /* B = Q*A*Q' or B = A(Q,Q) by MATLAB notation */
     double *d_x = NULL; // <double> n, x = A \ b 
     int *d_csrRowPtrB = NULL;
     int *d_csrColIndB = NULL;
     double *d_csrValB = NULL;
     double *d_b =NULL;
     double *d_r =NULL;    // residual
     double *d_z = NULL; /* z = B \ Q*b */
     double *d_Qb = NULL; /* a copy of h_Qb */
     int *d_mapBfromA = NULL;  /* <int> nnzA */   // Mapping for values array
     int *d_Q = NULL; /* <int> n */               // Mapping for column pointer and row pointer


    // the constant used in cusolverSp
    // singularity is -1 if A is invertible under tol
    // tol determines the condition of singularity
    int singularity = 0; 
    const double tol = 1.e-14;

    void initialiseLowLevelCudaQR();

    void lowLevelQRFactorisation(int* h_rowPtrA, int* h_colPtrB, 
                                        double* h_valuesA, int rowsA,int nnzA,
                                        const char* reorder);

    void lowLevelQRSolve(double* h_b);

    void resetLowLevelQRRoutines();

};


class cudaLowLevelQR_Optimised
{

    public :

        //Default COnstructor 
        cudaLowLevelQR_Optimised()
        {

        };
        cusolverSpHandle_t cusolverSpH = NULL; // reordering, permutation and 1st LU factorization
        cusparseHandle_t   cusparseH = NULL;   // residual evaluation
        cudaStream_t stream = NULL;
        cusparseMatDescr_t descrA = NULL; // A is a base-0 general matrix

        csrqrInfoHost_t h_info = NULL; // opaque info structure for LU with parital pivoting
        csrqrInfo_t d_info = NULL; // opaque info structure for LU with parital pivoting

        int rowsA = 0; // number of rows of A
        int colsA = 0; // number of columns of A
        int nnzA  = 0; // number of nonzeros of A
        int baseA = 0; // base index in CSR format


        double *h_x = NULL; // <double> n,  x = A \ b
        
           // ------------------ Added for Permutation Details ----------------------------- //
   
    double *h_Qb = NULL; /* Q*b */
    double *h_z  = NULL; /* z = B \ (Q*b) */

    int *h_Q = NULL; /* <int> n */
                     /* reorder to reduce zero fill-in */
                     /* Q = symrcm(A) or Q = symamd(A) */
    /* B = Q*A*Q' or B = A(Q,Q) by MATLAB notation */
    int *h_csrRowPtrB = NULL; /* <int> n+1 */
    int *h_csrColIndB = NULL; /* <int> nnzA */
    double *h_csrValB = NULL; /* <double> nnzA */
    int *h_mapBfromA = NULL;  /* <int> nnzA */

    // ------------- END -- Added for Permutation Details ----------------------------- //

    // ----- Values Ne


    size_t size_internal = 0; 
    size_t size_chol = 0; // size of working space for csrlu
    void *buffer_cpu = NULL; // working space for Cholesky
    void *buffer_gpu = NULL; // working space for Cholesky

    // ============= Allocatino in DEvice ================ //
        /* B = Q*A*Q' or B = A(Q,Q) by MATLAB notation */
     double *d_x = NULL; // <double> n, x = A \ b 
     int *d_csrRowPtrB = NULL;
     int *d_csrColIndB = NULL;
     double *d_csrValB = NULL;
     double *d_b =NULL;
     double *d_r =NULL;    // residual
     double *d_z = NULL; /* z = B \ Q*b */
     double *d_Qb = NULL; /* a copy of h_Qb */
     int *d_mapBfromA = NULL;  /* <int> nnzA */   // Mapping for values array
     int *d_Q = NULL; /* <int> n */               // Mapping for column pointer and row pointer


    // the constant used in cusolverSp
    // singularity is -1 if A is invertible under tol
    // tol determines the condition of singularity
    int singularity = 0; 
    const double tol = 1.e-14;

    void initialiseLowLevelCudaQR_optimised();


    // Reorders and Sets up Cuda Values for QR Factorisation
    void lowLevelQR_Factorise_optimised(int* h_rowPtrA, int* h_colPtrB, 
                                        double* h_valuesA, int rowsA,int nnzA,
                                        const char* reorder);
    

    // Performs QR Factorisation on the Matrix based on Previously obtained values
    void lowLevelQR_ReFactorise_Optimised(double* h_valuesA);


    void lowLevelQRSolve_optimised(double* h_b);


    // Clears only the input matrix , and keeps  d_info ( QR in Device ) , and the permutation Matrices
    void resetLowLevelQR_Refactorise_optimised();


    // Master Resets QR code for fresh Computation
    void masterResetlowLevelQR_optimised();

};

#endif