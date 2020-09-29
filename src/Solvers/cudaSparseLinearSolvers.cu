
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
 #include <fstream>
 #include <iostream>

 #include "cusolverSp.h"
 #include "cusolverRf.h"
 
 #include "helper_string.h"
 #include "helper_cusolver.h"
 
 #include "cusolverSp_LOWLEVEL_PREVIEW.h"
 
 #include <cuda_runtime.h>
 #include "helper_cuda.h"
 
#include "cudaSparseLinearSolvers.h"


#include <memory>
#include <cstring>


using namespace std;

void cudaRefactor::initialiseCudaHandles()
{
    checkCudaErrors(cusolverSpCreate(&cusolverSpH));
    checkCudaErrors(cusparseCreate(&cusparseH));
    checkCudaErrors(cudaStreamCreate(&stream));
    checkCudaErrors(cudaStreamCreate(&stream2));


    checkCudaErrors(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    checkCudaErrors(cudaStreamCreateWithFlags(&stream2, cudaStreamNonBlocking));



    checkCudaErrors(cusolverSpSetStream(cusolverSpH, stream));   // Cusolver Stream and Handle link
    checkCudaErrors(cusolverSpSetStream(cusolverSpH, stream2));   // Cusolver Stream and Handle link

    // checkCudaErrors(cusparseSetStream(cusolverRfH, stream2));       // Cusparse stream and handle link


    checkCudaErrors(cusparseCreateMatDescr(&descrA));
    checkCudaErrors(cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL));

    checkCudaErrors(cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO));

}


void cudaRefactor::destroyCudahandles()
{
    if (cusolverRfH) { checkCudaErrors(cusolverRfDestroy(cusolverRfH)); }
    if (cusolverSpH) { checkCudaErrors(cusolverSpDestroy(cusolverSpH)); }
    if (cusparseH  ) { checkCudaErrors(cusparseDestroy(cusparseH)); }
    if (stream     ) { checkCudaErrors(cudaStreamDestroy(stream)); }
    if (descrA     ) { checkCudaErrors(cusparseDestroyMatDescr(descrA)); }
    if (info       ) { checkCudaErrors(cusolverSpDestroyCsrluInfoHost(info)); }
}


void cudaRefactor::LU_DecompositionHost(int* h_csrRowPtrA, int* h_csrColIndA, 
                                        double* h_csrValA, int _rowsA, int _nnzA,
                                         double* h_b , const char* reorder)
{
   rowsA = _rowsA;
   colsA = _rowsA;
   nnzA = _nnzA;
//    cout << " rows A : " << rowsA << " colsA : " << colsA << " nnz : " << nnzA <<endl;

   double* h_x;   // LOcal Declaration
    cout << " -------LU FACTOR ----------"  <<endl;

    /// -------- Memory Allocation for CPU arrays --------------- // 
    h_Qreorder   = (int*)malloc(sizeof(int)*colsA);
    h_csrColIndB = new int[nnzA]();
    h_csrRowPtrB = new int[rowsA+1]();
    h_csrValB    = (double*)malloc(sizeof(double)*nnzA);
    h_mapBfromA  = (int*   )malloc(sizeof(int)*nnzA);

    h_x    = (double*)malloc(sizeof(double)*colsA);
    // h_b    = (double*)malloc(sizeof(double)*rowsA);
    h_xhat = (double*)malloc(sizeof(double)*colsA);
    h_bhat = (double*)malloc(sizeof(double)*rowsA);

    assert(NULL != h_Qreorder); assert(NULL != h_csrRowPtrB); assert(NULL != h_csrColIndB);assert(NULL != h_csrValB);
    assert(NULL != h_mapBfromA);assert(NULL != h_x);assert(NULL != h_xhat);assert(NULL != h_bhat);

    // cout << " INFO :  Reordering The Matrix using " ;

    if ( 0 == strcmp(reorder, "symrcm") )
    {
        printf("  SYMRCM ==>  Q = symrcm(A)  \n" );
        checkCudaErrors(cusolverSpXcsrsymrcmHost(
            cusolverSpH, rowsA, nnzA,
            descrA, h_csrRowPtrA, h_csrColIndA, 
            h_Qreorder));
    }
    else if ( 0 == strcmp(reorder, "symamd") )
    {
        printf("   SYMAND ==>      Q = symamd(A)  \n" );
        checkCudaErrors(cusolverSpXcsrsymamdHost(
            cusolverSpH, rowsA, nnzA,
            descrA, h_csrRowPtrA, h_csrColIndA, 
            h_Qreorder));
    }

    else if ( 0 == strcmp(reorder, "metis") )
    {
        printf("    METIS ==>    Q = metis(A)  \n" );
        cusolverSpXcsrmetisndHost(
            cusolverSpH, rowsA, nnzA,
            descrA, h_csrRowPtrA, h_csrColIndA,NULL ,
            h_Qreorder);
    }
    else 
    {
        printf(" Error in Reordering \n" );
        // exit(0);
    }

    // printf("step 3: B = Q*A*Q^T\n");
    memcpy(h_csrRowPtrB, h_csrRowPtrA, sizeof(int)*(rowsA+1));
    memcpy(h_csrColIndB, h_csrColIndA, sizeof(int)*nnzA);

    // cout << " memcopy done " <<endl;
    // checkCudaErrors(cudaDeviceSynchronize());

    start = second();
    start = second();

    checkCudaErrors(cusolverSpXcsrperm_bufferSizeHost(
        cusolverSpH, rowsA, colsA, nnzA,
        descrA, h_csrRowPtrB, h_csrColIndB,
        h_Qreorder, h_Qreorder,
        &size_perm));

    if (buffer_cpu) {
        free(buffer_cpu); 
    }
    buffer_cpu = (void*)malloc(sizeof(char)*size_perm);
    assert(NULL != buffer_cpu);

    // h_mapBfromA = Identity 
    for(int j = 0 ; j < nnzA ; j++){
        h_mapBfromA[j] = j;
    }
    checkCudaErrors(cusolverSpXcsrpermHost(
        cusolverSpH, rowsA, colsA, nnzA,
        descrA, h_csrRowPtrB, h_csrColIndB,
        h_Qreorder, h_Qreorder,
        h_mapBfromA,
        buffer_cpu));

    // B = A( mapBfromA )
    for(int j = 0 ; j < nnzA ; j++){
        h_csrValB[j] = h_csrValA[ h_mapBfromA[j] ];
    }

    stop = second();
    time_perm = stop - start;

    // printf("step 4: solve A*x = b by LU(B) in cusolverSp\n");

    // printf("step 4.1: create opaque info structure\n");
    checkCudaErrors(cusolverSpCreateCsrluInfoHost(&info));

    // printf("step 4.2: analyze LU(B) to know structure of Q and R, and upper bound for nnz(L+U)\n");
    start = second();
    start = second();

    checkCudaErrors(cusolverSpXcsrluAnalysisHost(
        cusolverSpH, rowsA, nnzA,
        descrA, h_csrRowPtrB, h_csrColIndB,
        info));

    stop = second();
    time_sp_analysis = stop - start;

    // printf("step 4.3: workspace for LU(B)\n");
    checkCudaErrors(cusolverSpDcsrluBufferInfoHost(
        cusolverSpH, rowsA, nnzA,
        descrA, h_csrValB, h_csrRowPtrB, h_csrColIndB,
        info,
        &size_internal,
        &size_lu));

    if (buffer_cpu) { 
        free(buffer_cpu); 
    }
    buffer_cpu = (void*)malloc(sizeof(char)*size_lu);
    assert(NULL != buffer_cpu);

    // printf("step 4.4: compute Ppivot*B = L*U \n");
    start = second();
    start = second();
    checkCudaErrors(cusolverSpDcsrluFactorHost(
        cusolverSpH, rowsA, nnzA,
        descrA, h_csrValB, h_csrRowPtrB, h_csrColIndB,
        info, pivot_threshold,
        buffer_cpu));

    stop = second();
    time_sp_factor = stop - start;

    // TODO: check singularity by tol
    // printf("step 4.5: check if the matrix is singular \n");
    checkCudaErrors(cusolverSpDcsrluZeroPivotHost(
        cusolverSpH, info, tol, &singularity));

    if ( 0 <= singularity)
    {
        fprintf(stderr, "Error: A is not invertible, singularity=%d\n", singularity);
        exit(0);
    }

    // printf("step 4.6: solve A*x = b \n");
    // printf("    i.e.  solve B*(Qx) = Q*b \n");
    start = second();
    start = second();

    // b_hat = Q*b
    for(int j = 0 ; j < rowsA ; j++){
        h_bhat[j] = h_b[h_Qreorder[j]];
    }
    // B*x_hat = b_hat
    checkCudaErrors(cusolverSpDcsrluSolveHost(
        cusolverSpH, rowsA, h_bhat, h_xhat, info, buffer_cpu));

    // x = Q^T * x_hat    // Solution Replaces B Vector
    for(int j = 0 ; j < rowsA ; j++){
        h_b[h_Qreorder[j]] = h_xhat[j];
    }

    stop = second();
    time_sp_solve = stop - start;

    // printf("step 5: extract P, Q, L and U from P*B*Q^T = L*U \n");
    // printf("        L has implicit unit diagonal\n");
    start = second();
    start = second();

    checkCudaErrors(cusolverSpXcsrluNnzHost(
        cusolverSpH,
        &nnzL,
        &nnzU,
        info));

    h_Plu = (int*)malloc(sizeof(int)*rowsA);
    h_Qlu = (int*)malloc(sizeof(int)*colsA);

    h_csrValL    = (double*)malloc(sizeof(double)*nnzL);
    h_csrRowPtrL = (int*)malloc(sizeof(int)*(rowsA+1)); 
    h_csrColIndL = (int*)malloc(sizeof(int)*nnzL);

    h_csrValU    = (double*)malloc(sizeof(double)*nnzU);
    h_csrRowPtrU = (int*)malloc(sizeof(int)*(rowsA+1)); 
    h_csrColIndU = (int*)malloc(sizeof(int)*nnzU);

    assert(NULL != h_Plu);
    assert(NULL != h_Qlu);

    assert(NULL != h_csrValL);
    assert(NULL != h_csrRowPtrL);
    assert(NULL != h_csrColIndL);

    assert(NULL != h_csrValU);
    assert(NULL != h_csrRowPtrU);
    assert(NULL != h_csrColIndU);

    checkCudaErrors(cusolverSpDcsrluExtractHost(
        cusolverSpH,
        h_Plu,
        h_Qlu,
        descrA,
        h_csrValL, 
        h_csrRowPtrL,
        h_csrColIndL,
        descrA,
        h_csrValU,
        h_csrRowPtrU,
        h_csrColIndU,
        info,
        buffer_cpu));

    stop = second();
    time_sp_extract = stop - start;

    // printf("nnzL = %d, nnzU = %d\n", nnzL, nnzU);
    

    printf("step 6: form P*A*Q^T = L*U\n");

    
    h_P = (int*)malloc(sizeof(int)*rowsA);
    h_Q = (int*)malloc(sizeof(int)*colsA);
    cudaHostAlloc( (void**)&h_P, sizeof(int)*rowsA, cudaHostAllocDefault) ;
    cudaHostAlloc((void**)&h_Q, sizeof(int)*colsA, cudaHostAllocDefault);
    // checkCudaErrors( cudaHostRegister(h_P, sizeof(int)*rowsA, cudaHostRegisterPortable) );
    // checkCudaErrors( cudaHostRegister(h_Q, sizeof(int)*colsA, cudaHostRegisterPortable) );


    assert(NULL != h_P);
    assert(NULL != h_Q);

    // printf("step 6.1: P = Plu*Qreroder\n");
    // gather operation, P = Qreorder(Plu)
    for(int j = 0 ; j < rowsA ; j++){
        h_P[j] = h_Qreorder[h_Plu[j]];
    }

    // printf("step 6.2: Q = Qlu*Qreorder \n");
    // gather operation, Q = Qreorder(Qlu)
    for(int j = 0 ; j < colsA ; j++){
        h_Q[j] = h_Qreorder[h_Qlu[j]];
    }

    delete[] h_x;


        printf("--- REFACTORIZE ");
    checkCudaErrors(cusolverRfCreate(&cusolverRfH));



    printf("step 8: set parameters for cusolverRf \n");
    // numerical values for checking "zeros" and for boosting.
    checkCudaErrors(cusolverRfSetNumericProperties(cusolverRfH, nzero, nboost));

    // choose algorithm for refactorization and solve
    checkCudaErrors(cusolverRfSetAlgs(cusolverRfH, fact_alg, solve_alg));

    // matrix mode: L and U are CSR format, and L has implicit unit diagonal
    checkCudaErrors(cusolverRfSetMatrixFormat(
        cusolverRfH, CUSOLVERRF_MATRIX_FORMAT_CSR, CUSOLVERRF_UNIT_DIAGONAL_ASSUMED_L));

    // fast mode for matrix assembling
    checkCudaErrors(cusolverRfSetResetValuesFastMode(
        cusolverRfH, CUSOLVERRF_RESET_VALUES_FAST_MODE_ON));
    // cout << " Rows a : "<< rowsA << " colsa : " << colsA << " nnzA : "<< nnzA <<endl; 

    checkCudaErrors(cudaMalloc((void **)&d_csrRowPtrA, sizeof(int)*(rowsA+1)));
    checkCudaErrors(cudaMalloc((void **)&d_csrColIndA, sizeof(int)*nnzA));
    checkCudaErrors(cudaMalloc((void **)&d_csrValA   , sizeof(double)*nnzA));
    checkCudaErrors(cudaMalloc((void **)&d_x, sizeof(double)*colsA));
    checkCudaErrors(cudaMalloc((void **)&d_b, sizeof(double)*rowsA));
    checkCudaErrors(cudaMalloc((void **)&d_P, sizeof(int)*rowsA));
    checkCudaErrors(cudaMalloc((void **)&d_Q, sizeof(int)*colsA));
    checkCudaErrors(cudaMalloc((void **)&d_T, sizeof(double)*rowsA*1));
    // cout<< " Finished LU Host " <<endl;
}



void cudaRefactor::cudaRefactorize(int* h_csrRowPtrA, int* h_csrColIndA, double* h_csrValA, int rowsA, int N_iteration)
{


    checkCudaErrors(cudaMemcpyAsync(d_csrRowPtrA, h_csrRowPtrA, sizeof(int)*(rowsA+1), cudaMemcpyHostToDevice,stream));
    checkCudaErrors(cudaMemcpyAsync(d_csrColIndA, h_csrColIndA, sizeof(int)*nnzA     , cudaMemcpyHostToDevice,stream));
    checkCudaErrors(cudaMemcpyAsync(d_csrValA   , h_csrValA   , sizeof(double)*nnzA  , cudaMemcpyHostToDevice,stream));

            // Allocate Memory on Device
    if(N_iteration == 1)
    {
        checkCudaErrors(cudaMemcpyAsync(d_P, h_P, sizeof(int)*rowsA, cudaMemcpyHostToDevice,stream));
        checkCudaErrors(cudaMemcpyAsync(d_Q, h_Q, sizeof(int)*colsA, cudaMemcpyHostToDevice,stream));
        checkCudaErrors(cusolverRfSetupHost(
                rowsA, nnzA, h_csrRowPtrA, h_csrColIndA, h_csrValA,nnzL, 
                h_csrRowPtrL, h_csrColIndL, h_csrValL, nnzU, h_csrRowPtrU, h_csrColIndU, h_csrValU, h_P, h_Q, 
                cusolverRfH));
        checkCudaErrors(cusolverRfAnalyze(cusolverRfH));
    }

    start = second();
    // checkCudaErrors(cudaStreamSynchronize(stream));
    
    checkCudaErrors(cusolverRfResetValues(
        rowsA,nnzA,
        d_csrRowPtrA, d_csrColIndA, d_csrValA,
        d_P,
        d_Q,
        cusolverRfH));
    

    checkCudaErrors(cusolverRfRefactor(cusolverRfH));

    // checkCudaErrors(cudaDeviceSynchronize());
    stop = second();
    time_rf_refactor = stop - start;
    printf("time_rf_refactor :  %f", time_rf_refactor);
    // checkCudaErrors(cudaDeviceSynchronize());


    // Delete the Non used Device arrays

}


void cudaRefactor::cudaRefactorSolve(double* h_b)
{
    // cout << "----- Refactor and solve --- " <<endl;
    
    checkCudaErrors(cudaMemcpyAsync(d_x, h_b, sizeof(double)*rowsA, cudaMemcpyHostToDevice,stream));
    start = second();
    start = second();

    checkCudaErrors(cudaStreamSynchronize(stream));
    
    checkCudaErrors(cusolverRfSolve(cusolverRfH, d_P, d_Q, 1, d_T, rowsA, d_x, rowsA));

    // checkCudaErrors(cudaDeviceSynchronize());
    stop = second();
    time_rf_solve = stop - start;

    checkCudaErrors(cudaMemcpy(h_b, d_x, sizeof(double)*colsA, cudaMemcpyDeviceToHost));

    // printf("===== timing profile \n");
    // // printf(" reorder A   : %f sec\n", time_reorder);
    // printf(" B = Q*A*Q^T : %f sec\n", time_perm);s
    // printf("\n");
    // printf(" cusolverSp LU analysis: %f sec\n", time_sp_analysis);
    // printf(" cusolverSp LU factor  : %f sec\n", time_sp_factor);
    // printf(" cusolverSp LU solve   : %f sec\n", time_sp_solve);
    // printf(" cusolverSp LU extract : %f sec\n", time_sp_extract);
    // printf("\n");
    // printf(" cusolverRf assemble : %f sec\n", time_rf_assemble);
    // printf(" cusolverRf reset    : %f sec\n", time_rf_reset);
    // printf(" cusolverRf refactor : %f sec\n", time_rf_refactor);
    // printf(" cusolverRf solve    : %f sec\n", time_rf_solve);



   
}


// Deletes all the LU factor Entities that was created as part of the Routines and calls back the Destroy handles routines
void cudaRefactor::resetCudaRF()
{
    // cout << " ------ REST CUDA RF " <<endl;

    // Delete all array created 
    if (h_Qreorder  ) { free(h_Qreorder); }
    
    if (h_csrRowPtrB) { free(h_csrRowPtrB); }
    if (h_csrColIndB) { free(h_csrColIndB); }
    if (h_csrValB   ) { free(h_csrValB   ); }
    if (h_mapBfromA ) { free(h_mapBfromA ); }

    if (h_xhat) { free(h_xhat); }
    if (h_bhat) { free(h_bhat); }

    // if (buffer_cpu) { free(buffer_cpu); }

    if (h_Plu) { free(h_Plu); }
    if (h_Qlu) { free(h_Qlu); }
    if (h_csrRowPtrL) { free(h_csrRowPtrL); }
    if (h_csrColIndL) { free(h_csrColIndL); }
    if (h_csrValL   ) { free(h_csrValL   ); }
    if (h_csrRowPtrU) { free(h_csrRowPtrU); }
    if (h_csrColIndU) { free(h_csrColIndU); }
    if (h_csrValU   ) { free(h_csrValU   ); }

    if (h_P) { checkCudaErrors(cudaFreeHost(h_P)); }
    if (h_Q) { checkCudaErrors(cudaFreeHost(h_Q)); }

    
    // if (h_P) { free(h_P); }
    // if (h_Q) { free(h_Q); }

    if (d_x) { checkCudaErrors(cudaFree(d_x)); }
    if (d_P) { checkCudaErrors(cudaFree(d_P)); }
    if (d_Q) { checkCudaErrors(cudaFree(d_Q)); }
    if (d_T) { checkCudaErrors(cudaFree(d_T)); }

    if (d_csrValA   ) { checkCudaErrors(cudaFree(d_csrValA)); }
    if (d_csrRowPtrA) { checkCudaErrors(cudaFree(d_csrRowPtrA)); }
    if (d_csrColIndA) { checkCudaErrors(cudaFree(d_csrColIndA)); }
    if (d_b) { checkCudaErrors(cudaFree(d_b)); }

    // Destroy ALl the handles
    destroyCudahandles();
}























































void cudaLowLevelQR::initialiseLowLevelCudaQR()
{

    checkCudaErrors(cusolverSpCreate(&cusolverSpH));
    checkCudaErrors(cusparseCreate(&cusparseH));
    checkCudaErrors(cudaStreamCreate(&stream));
    checkCudaErrors(cusolverSpSetStream(cusolverSpH, stream));
    checkCudaErrors(cusparseSetStream(cusparseH, stream));

    checkCudaErrors(cusparseCreateMatDescr(&descrA));

    checkCudaErrors(cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL));

    if (baseA) 
    {
        checkCudaErrors(cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ONE));
    }
    else
    {
        checkCudaErrors(cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO));
    }
}


void cudaLowLevelQR::lowLevelQRFactorisation(int* h_csrRowPtrA, int* h_csrColIndA, 
                                double* h_csrValA, int _rowsA,int _nnzA
                                 ,const char* reorder)

{

    rowsA = _rowsA;
    colsA = _rowsA;
    nnzA  = _nnzA;

    double tt1 , tt2;

    tt1 = second();
    // Allocate Memory on the Host
    h_x        = (double*)malloc(sizeof(double)*colsA);
    h_Qb       = (double*)malloc(sizeof(double)*rowsA);
   //  h_r        = (double*)malloc(sizeof(double)*rowsA);
    h_Q          = (int*   )malloc(sizeof(int)*colsA);
    h_csrRowPtrB = (int*   )malloc(sizeof(int)*(rowsA+1));
    h_csrColIndB = (int*   )malloc(sizeof(int)*nnzA);
    h_csrValB    = (double*)malloc(sizeof(double)*nnzA);
    h_mapBfromA  = (int*   )malloc(sizeof(int)*nnzA);

    tt2 = second();

    cout << " Malloc : " << tt2 - tt1  <<endl;


    assert(NULL != h_x);
    assert(NULL != h_csrRowPtrB );
    assert(NULL != h_csrColIndB );
    assert(NULL != h_csrValB    );
    assert(NULL != h_mapBfromA  );
    assert(NULL != h_Q  );
    



    //Allocate Memory on the Device 
     checkCudaErrors(cudaMalloc((void **)&d_csrRowPtrB, sizeof(int)*(rowsA+1)));
     checkCudaErrors(cudaMalloc((void **)&d_csrColIndB, sizeof(int)*nnzA));
     checkCudaErrors(cudaMalloc((void **)&d_csrValB   , sizeof(double)*nnzA));
     checkCudaErrors(cudaMalloc((void **)&d_Q, sizeof(int)*colsA));
     checkCudaErrors(cudaMalloc((void **)&d_z, sizeof(double)*colsA));
     checkCudaErrors(cudaMalloc((void **)&d_x, sizeof(double)*colsA));
     checkCudaErrors(cudaMalloc((void **)&d_Qb, sizeof(double)*rowsA));

        // tt1 = second();
     // Reorder the Host Matrix
     if ( 0 == strcmp(reorder, "symrcm") )
     {
        printf(" SYMRCM Reordering Performed ");
         printf("step 2.1: Q = symrcm(A) \n");
         checkCudaErrors(cusolverSpXcsrsymrcmHost(
            cusolverSpH, rowsA, nnzA,
             descrA, h_csrRowPtrA,h_csrColIndA,
             h_Q));
     }
     else if ( 0 == strcmp(reorder, "symamd") )
     {
        printf(" SYMAND Reordering Performed ");
         printf("step 2.1: Q = symamd(A) \n");
         checkCudaErrors(cusolverSpXcsrsymamdHost(
            cusolverSpH, rowsA, nnzA,
             descrA,h_csrRowPtrA,h_csrColIndA,
             h_Q));
     }
     else if ( 0 == strcmp(reorder, "metis") )
     {
        //  printf("step 2.1: Q = metis(A) \n");
        checkCudaErrors(cusolverSpXcsrmetisndHost(
           cusolverSpH, rowsA, nnzA,
          descrA, h_csrRowPtrA, h_csrColIndA,
            NULL, /* default setting. */
            h_Q)); 

     }
     else 
     {
         fprintf(stderr, "Error: %s is unknown reordering\n", reorder);
         exit(0);
     }


    // Get a current copy of the matrix to array B
    memcpy(h_csrRowPtrB, h_csrRowPtrA, sizeof(int)*(rowsA+1));
    memcpy(h_csrColIndB, h_csrColIndA, sizeof(int)*nnzA);


     /* h_mapBfromA = Identity */
     for(int j = 0 ; j < nnzA ; j++)
     h_mapBfromA[j] = j;
 
    size_t size_perm;
    checkCudaErrors(cusolverSpXcsrperm_bufferSizeHost(
         cusolverSpH, rowsA, colsA, nnzA,
         descrA, h_csrRowPtrB, h_csrColIndB,
         h_Q, h_Q,
         &size_perm));
 
    //  if (buffer_cpu)
    //  {
    //      free(buffer_cpu);
    //  }
     buffer_cpu = (void*)malloc(sizeof(char)*size_perm);
     assert(NULL != buffer_cpu);
 
     // Compute Permutation 
     checkCudaErrors(cusolverSpXcsrpermHost(
         cusolverSpH, rowsA, colsA, nnzA,
         descrA, h_csrRowPtrB, h_csrColIndB,
         h_Q, h_Q,
         h_mapBfromA,
         buffer_cpu));
     
     // Calculate B = Q*A*Q'
     for(int j = 0 ; j < nnzA ; j++)
     {
         h_csrValB[j] = h_csrValA[ h_mapBfromA[j] ];
     }

     

    //  cudaStreamSynchronize(stream);
    //  tt2 = second();
    //  cout << " PRe QR TIME TAKEN " << tt2 - tt1 <<endl;

    // Get a current copy of the matrix to array B
    checkCudaErrors(cudaMemcpyAsync(d_csrRowPtrB, h_csrRowPtrB, sizeof(int)*(rowsA+1), cudaMemcpyHostToDevice, stream));
    checkCudaErrors(cudaMemcpyAsync(d_csrColIndB, h_csrColIndB, sizeof(int)*nnzA     , cudaMemcpyHostToDevice, stream));
    checkCudaErrors(cudaMemcpyAsync(d_csrValB   , h_csrValB   , sizeof(double)*nnzA  , cudaMemcpyHostToDevice, stream));
    checkCudaErrors(cudaMemcpyAsync(d_Q         , h_Q         , sizeof(int)*rowsA    , cudaMemcpyHostToDevice, stream));



    // printf("step 2: create opaque info structure\n");
    checkCudaErrors(cusolverSpCreateCsrqrInfo(&d_info));

    // cudaStreamSynchronize(stream); tt1 = second();
    // printf("step 3: analyze qr(A) to know structure of L\n");
    checkCudaErrors(cusolverSpXcsrqrAnalysis(
    cusolverSpH, rowsA, colsA, nnzA,
    descrA, d_csrRowPtrB, d_csrColIndB,
    d_info));
    // cudaStreamSynchronize(stream); tt2 = second();
    // cout << " QR Analysis TIME TAKEN " << tt2 - tt1 <<endl;
    
    // cudaStreamSynchronize(stream); tt1 = second();

    // printf("step 4: workspace for qr(A)\n");
    checkCudaErrors(cusolverSpDcsrqrBufferInfo(
           cusolverSpH, rowsA, colsA, nnzA,
           descrA, d_csrValB, d_csrRowPtrB, d_csrColIndB,
           d_info,
           &size_internal,
           &size_chol));
    
    //    printf("GPU buffer size = %lld bytes\n", (signed long long)size_chol);
    //    if (buffer_gpu) {
    //        checkCudaErrors(cudaFree(buffer_gpu));
    //    }
   checkCudaErrors(cudaMalloc(&buffer_gpu, sizeof(char)*size_chol));

   cudaStreamSynchronize(stream); tt2 = second(); 
   cout << " QR Buffer setup TIME TAKEN " << tt2 - tt1 <<endl;
   
   const double zero = 0.0;

    cudaStreamSynchronize(stream);      tt1 = second(); 

    //    printf("step 5: compute A = L*L^T \n");
   checkCudaErrors(cusolverSpDcsrqrSetup(
           cusolverSpH, rowsA, colsA, nnzA,
           descrA, d_csrValB, d_csrRowPtrB, d_csrColIndB,
           zero,
           d_info));

    cudaStreamSynchronize(stream); tt2 = second();
    cout << " QR SETUP TIME TAKEN " << tt2 - tt1 <<endl;

    double t1  =  second();
    checkCudaErrors(cusolverSpDcsrqrFactor(
        cusolverSpH, rowsA, colsA, nnzA,
        NULL, NULL,
        d_info,
        buffer_gpu));
    
    // cudaStreamSynchronize(stream);
    // double t2  =  second();
    // cout << " QR FACTORISATION in kernel : " << t2 - t1 <<endl;


    //    cudaStreamSynchronize(stream);


    //    printf("step 6: check if the matrix is singular \n");
    //    checkCudaErrors(cusolverSpDcsrqrZeroPivot(
        //    cusolverSpH, d_info, tol, &singularity));

    //    if ( 0 <= singularity){
    //        fprintf(stderr, "Error: A is not invertible, singularity=%d\n", singularity);
    //         exit(0);
    //    }
        

}


void cudaLowLevelQR::lowLevelQRSolve(double* h_b)
{
    
    for(int row = 0 ; row < rowsA ; row++)
        h_Qb[row] = h_b[h_Q[row]];
        
    checkCudaErrors(cudaMemcpyAsync(d_Qb , h_Qb  , sizeof(double)*rowsA , cudaMemcpyHostToDevice, stream));


    // Solve the System 
    printf("step 7: solve A*x = b \n");
    checkCudaErrors(cusolverSpDcsrqrSolve(
            cusolverSpH, rowsA, colsA, d_Qb, d_z, d_info, buffer_gpu));
     
    //    cudaStreamSynchronize(stream);

    // SCatter the Solution Vector based on the Permutation Vector
     
    checkCudaErrors(cusparseDsctr (cusparseH,
                                    rowsA,d_z,d_Q,d_x,
                                    CUSPARSE_INDEX_BASE_ZERO));
    
    checkCudaErrors(cudaMemcpyAsync(h_b, d_x, sizeof(double)*rowsA, cudaMemcpyDeviceToHost,stream));
    
    cudaStreamSynchronize(stream);
}


void cudaLowLevelQR::resetLowLevelQRRoutines()
{
    if (cusolverSpH) { checkCudaErrors(cusolverSpDestroy(cusolverSpH)); }
    if (cusparseH  ) { checkCudaErrors(cusparseDestroy(cusparseH)); }
    if (stream     ) { checkCudaErrors(cudaStreamDestroy(stream)); }
    if (descrA     ) { checkCudaErrors(cusparseDestroyMatDescr(descrA)); }
   //  if (h_info     ) { checkCudaErrors(cusolverSpDestroyCsrqrInfoHost(h_info)); }
    if (d_info     ) { checkCudaErrors(cusolverSpDestroyCsrqrInfo(d_info)); }


    // Delete HOst Allocated Arrays
    if (h_csrValB  ) { free(h_csrValB); }
    if (h_csrRowPtrB) { free(h_csrRowPtrB); }
    if (h_csrColIndB) { free(h_csrColIndB); }
    if (h_mapBfromA) { free(h_mapBfromA); }
    if (h_x   ) { free(h_x); }
    if (h_Q  ) { free(h_Q); }
    if (h_Qb){ free(h_Qb); }

    // Delete Temp Buffers
    if (buffer_cpu) { free(buffer_cpu); }
    if (buffer_gpu) { checkCudaErrors(cudaFree(buffer_gpu)); }

    // Delete all Allocated Arrays - Device
    if (d_csrValB  ) { checkCudaErrors(cudaFree(d_csrValB)); }
     if (d_csrRowPtrB) { checkCudaErrors(cudaFree(d_csrRowPtrB)); }
     if (d_csrColIndB) { checkCudaErrors(cudaFree(d_csrColIndB)); }
     if (d_x) { checkCudaErrors(cudaFree(d_x)); }
     if (d_b) { checkCudaErrors(cudaFree(d_b)); }
     if (d_Q) { checkCudaErrors(cudaFree(d_Q)); }
     if (d_Qb) { checkCudaErrors(cudaFree(d_Qb)); }

     if (d_z) { checkCudaErrors(cudaFree(d_z)); }

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void cudaLowLevelQR_Optimised::initialiseLowLevelCudaQR_optimised()
{

    checkCudaErrors(cusolverSpCreate(&cusolverSpH));
    checkCudaErrors(cusparseCreate(&cusparseH));
    checkCudaErrors(cudaStreamCreate(&stream));
    checkCudaErrors(cusolverSpSetStream(cusolverSpH, stream));
    checkCudaErrors(cusparseSetStream(cusparseH, stream));

    checkCudaErrors(cusparseCreateMatDescr(&descrA));

    checkCudaErrors(cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL));


    // ------------- Opaque Data Structure for QR Factorisation -------------------------- //
    checkCudaErrors(cusolverSpCreateCsrqrInfo(&d_info));


    if (baseA) 
    {
        checkCudaErrors(cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ONE));
    }
    else
    {
        checkCudaErrors(cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO));
    }


}

void cudaLowLevelQR_Optimised::lowLevelQR_Factorise_optimised(int* h_csrRowPtrA, int* h_csrColIndA, 
                                                    double* h_csrValA, int _rowsA,int _nnzA
                                                    ,const char* reorder)
{

    rowsA = _rowsA;
    colsA = _rowsA;
    nnzA  = _nnzA;

    // double tt1 , tt2;

    // tt1 = second();
    // Allocate Memory on the Host
    h_Qb       = (double*)malloc(sizeof(double)*rowsA);
   //  h_r        = (double*)malloc(sizeof(double)*rowsA);
    h_Q          = (int*   )malloc(sizeof(int)*colsA);
    h_csrRowPtrB = (int*   )malloc(sizeof(int)*(rowsA+1));
    h_csrColIndB = (int*   )malloc(sizeof(int)*nnzA);
    h_csrValB    = (double*)malloc(sizeof(double)*nnzA);
    h_mapBfromA  = (int*   )malloc(sizeof(int)*nnzA);

    // tt2 = second();

    // cout << " Malloc : " << tt2 - tt1  <<endl;


    assert(NULL != h_csrRowPtrB );
    assert(NULL != h_csrColIndB );
    assert(NULL != h_csrValB    );
    assert(NULL != h_mapBfromA  );
    assert(NULL != h_Q  );
    assert(NULL != h_Qb  );

    



    //Allocate Memory on the Device 
     checkCudaErrors(cudaMalloc((void **)&d_csrRowPtrB, sizeof(int)*(rowsA+1)));
     checkCudaErrors(cudaMalloc((void **)&d_csrColIndB, sizeof(int)*nnzA));
     checkCudaErrors(cudaMalloc((void **)&d_csrValB   , sizeof(double)*nnzA));
     checkCudaErrors(cudaMalloc((void **)&d_Q, sizeof(int)*colsA));
     checkCudaErrors(cudaMalloc((void **)&d_z, sizeof(double)*colsA));
     checkCudaErrors(cudaMalloc((void **)&d_x, sizeof(double)*colsA));
     checkCudaErrors(cudaMalloc((void **)&d_Qb, sizeof(double)*rowsA));

        // tt1 = second();
     // Reorder the Host Matrix
     if ( 0 == strcmp(reorder, "symrcm") )
     {
        printf(" SYMRCM Reordering Performed ");
         printf("step 2.1: Q = symrcm(A) \n");
         checkCudaErrors(cusolverSpXcsrsymrcmHost(
            cusolverSpH, rowsA, nnzA,
             descrA, h_csrRowPtrA,h_csrColIndA,
             h_Q));
     }
     else if ( 0 == strcmp(reorder, "symamd") )
     {
        printf(" SYMAND Reordering Performed ");
         printf("step 2.1: Q = symamd(A) \n");
         checkCudaErrors(cusolverSpXcsrsymamdHost(
            cusolverSpH, rowsA, nnzA,
             descrA,h_csrRowPtrA,h_csrColIndA,
             h_Q));
     }
     else if ( 0 == strcmp(reorder, "metis") )
     {
         printf("step 2.1: Q = metis(A) \n");
        checkCudaErrors(cusolverSpXcsrmetisndHost(
           cusolverSpH, rowsA, nnzA,
          descrA, h_csrRowPtrA, h_csrColIndA,
            NULL, /* default setting. */
            h_Q)); 

     }
     else 
     {
         fprintf(stderr, "Error: %s is unknown reordering\n", reorder);
         exit(0);
     }


    // Get a current copy of the matrix to array B
    memcpy(h_csrRowPtrB, h_csrRowPtrA, sizeof(int)*(rowsA+1));
    memcpy(h_csrColIndB, h_csrColIndA, sizeof(int)*nnzA);


     /* h_mapBfromA = Identity */
     for(int j = 0 ; j < nnzA ; j++)
     h_mapBfromA[j] = j;
 
    size_t size_perm;
    checkCudaErrors(cusolverSpXcsrperm_bufferSizeHost(
         cusolverSpH, rowsA, colsA, nnzA,
         descrA, h_csrRowPtrB, h_csrColIndB,
         h_Q, h_Q,
         &size_perm));
 
    //  if (buffer_cpu)
    //  {
    //      free(buffer_cpu);
    //  }
     buffer_cpu = (void*)malloc(sizeof(char)*size_perm);
     assert(NULL != buffer_cpu);
 
     // Compute Permutation 
     checkCudaErrors(cusolverSpXcsrpermHost(
         cusolverSpH, rowsA, colsA, nnzA,
         descrA, h_csrRowPtrB, h_csrColIndB,
         h_Q, h_Q,
         h_mapBfromA,
         buffer_cpu));

        // Calculate B = Q*A*Q'
     for(int j = 0 ; j < nnzA ; j++)
     {
         h_csrValB[j] = h_csrValA[ h_mapBfromA[j] ];
     }



    // Get a current copy of the matrix to array B
    checkCudaErrors(cudaMemcpyAsync(d_csrRowPtrB, h_csrRowPtrB, sizeof(int)*(rowsA+1), cudaMemcpyHostToDevice, stream));
    checkCudaErrors(cudaMemcpyAsync(d_csrColIndB, h_csrColIndB, sizeof(int)*nnzA     , cudaMemcpyHostToDevice, stream));
    checkCudaErrors(cudaMemcpyAsync(d_csrValB   , h_csrValB   , sizeof(double)*nnzA  , cudaMemcpyHostToDevice, stream));
    checkCudaErrors(cudaMemcpyAsync(d_Q         , h_Q         , sizeof(int)*rowsA    , cudaMemcpyHostToDevice, stream));



    // printf("step 2: create opaque info structure\n");
    checkCudaErrors(cusolverSpCreateCsrqrInfo(&d_info));

    // cudaStreamSynchronize(stream); tt1 = second();
    // printf("step 3: analyze qr(A) to know structure of L\n");
    checkCudaErrors(cusolverSpXcsrqrAnalysis(
                    cusolverSpH, rowsA, colsA, nnzA,
                    descrA, d_csrRowPtrB, d_csrColIndB,
                    d_info));
    
    // printf("step 4: workspace for qr(A)\n");
    checkCudaErrors(cusolverSpDcsrqrBufferInfo(
        cusolverSpH, rowsA, colsA, nnzA,
        descrA, d_csrValB, d_csrRowPtrB, d_csrColIndB,
        d_info,
        &size_internal,
        &size_chol));
      
        checkCudaErrors(cudaMalloc(&buffer_gpu, sizeof(char)*size_chol));

        const double zero = 0.0;

        // printf("step 5: compute A = L*L^T \n");
        checkCudaErrors(cusolverSpDcsrqrSetup(
                cusolverSpH, rowsA, colsA, nnzA,
                descrA, d_csrValB, d_csrRowPtrB, d_csrColIndB,
                zero,
                d_info));

        checkCudaErrors(cusolverSpDcsrqrFactor(
            cusolverSpH, rowsA, colsA, nnzA,
            NULL, NULL,
            d_info,
            buffer_gpu));
}

void cudaLowLevelQR_Optimised::lowLevelQR_ReFactorise_Optimised(double* h_valA)
{

    // Allocate  Memory For values Array only
    h_csrValB    = (double*)malloc(sizeof(double)*nnzA);
    assert(NULL != h_csrValB);

    // Map the Values array based on permutation matrix
    for(int j = 0 ; j < nnzA ; j++)
    {
        h_csrValB[j] = h_valA[ h_mapBfromA[j] ];
    }

    // Allocate Memory on GPU
    checkCudaErrors(cudaMalloc((void **)&d_csrValB   , sizeof(double)*nnzA));
    checkCudaErrors(cudaMalloc((void **)&d_z, sizeof(double)*colsA));
    checkCudaErrors(cudaMalloc((void **)&d_x, sizeof(double)*colsA));
    checkCudaErrors(cudaMalloc((void **)&d_Qb, sizeof(double)*rowsA));

    checkCudaErrors(cudaMemcpyAsync(d_csrValB   , h_csrValB   , sizeof(double)*nnzA  , cudaMemcpyHostToDevice, stream));


    checkCudaErrors(cudaMalloc(&buffer_gpu, sizeof(char)*size_chol));
    
    double zero = 0.0;
    checkCudaErrors(cusolverSpDcsrqrSetup(
        cusolverSpH, rowsA, colsA, nnzA,
        descrA, d_csrValB, d_csrRowPtrB, d_csrColIndB,
        zero,
        d_info));

    checkCudaErrors(cusolverSpDcsrqrFactor(
        cusolverSpH, rowsA, colsA, nnzA,
        NULL, NULL,
        d_info,
        buffer_gpu));

}



void cudaLowLevelQR_Optimised::lowLevelQRSolve_optimised(double* h_b)
{
    
    for(int row = 0 ; row < rowsA ; row++)
        h_Qb[row] = h_b[h_Q[row]];
        
    checkCudaErrors(cudaMemcpyAsync(d_Qb , h_Qb  , sizeof(double)*rowsA , cudaMemcpyHostToDevice, stream));


    // Solve the System 
    printf("step 7: solve A*x = b \n");
    checkCudaErrors(cusolverSpDcsrqrSolve(
            cusolverSpH, rowsA, colsA, d_Qb, d_z, d_info, buffer_gpu));
     
    //    cudaStreamSynchronize(stream);

    // SCatter the Solution Vector based on the Permutation Vector
     
    checkCudaErrors(cusparseDsctr (cusparseH,
                                    rowsA,d_z,d_Q,d_x,
                                    CUSPARSE_INDEX_BASE_ZERO));
    
    checkCudaErrors(cudaMemcpyAsync(h_b, d_x, sizeof(double)*rowsA, cudaMemcpyDeviceToHost,stream));
    
    cudaStreamSynchronize(stream);
}


void cudaLowLevelQR_Optimised::resetLowLevelQR_Refactorise_optimised()
{

    // Delete HOst Allocated Arrays
    if (h_csrValB  ) { free(h_csrValB); }



    // Delete Temp Buffers
    if (buffer_cpu) { free(buffer_cpu); }
    if (buffer_gpu) { checkCudaErrors(cudaFree(buffer_gpu)); }

    // Delete all Allocated Arrays - Device
    if (d_csrValB  ) { checkCudaErrors(cudaFree(d_csrValB)); }
     if (d_x) { checkCudaErrors(cudaFree(d_x)); }
     if (d_Qb) { checkCudaErrors(cudaFree(d_Qb)); }
     if (d_z) { checkCudaErrors(cudaFree(d_z)); }

}


void cudaLowLevelQR_Optimised::masterResetlowLevelQR_optimised()
{
    if (cusolverSpH) { checkCudaErrors(cusolverSpDestroy(cusolverSpH)); }
    if (cusparseH  ) { checkCudaErrors(cusparseDestroy(cusparseH)); }
    if (stream     ) { checkCudaErrors(cudaStreamDestroy(stream)); }
    if (descrA     ) { checkCudaErrors(cusparseDestroyMatDescr(descrA)); }
    if (d_info     ) { checkCudaErrors(cusolverSpDestroyCsrqrInfo(d_info)); }

    // Delete HOst Allocated Arrays
    if (h_csrValB  ) { free(h_csrValB); }
    if (h_csrRowPtrB) { free(h_csrRowPtrB); }
    if (h_csrColIndB) { free(h_csrColIndB); }
    if (h_mapBfromA) { free(h_mapBfromA); }
    // if (h_x   ) { free(h_x); }
    if (h_Q  ) { free(h_Q); }
    if (h_Qb){ free(h_Qb); }

    // Delete Temp Buffers
    // if (buffer_cpu) { free(buffer_cpu); }
    if (buffer_gpu) { checkCudaErrors(cudaFree(buffer_gpu)); }



    // Delete all Allocated Arrays - Device
    if (d_csrValB  ) { checkCudaErrors(cudaFree(d_csrValB)); }
     if (d_csrRowPtrB) { checkCudaErrors(cudaFree(d_csrRowPtrB)); }
     if (d_csrColIndB) { checkCudaErrors(cudaFree(d_csrColIndB)); }
     if (d_x) { checkCudaErrors(cudaFree(d_x)); }
     if (d_b) { checkCudaErrors(cudaFree(d_b)); }
     if (d_Q) { checkCudaErrors(cudaFree(d_Q)); }
     if (d_Qb) { checkCudaErrors(cudaFree(d_Qb)); }

     if (d_z) { checkCudaErrors(cudaFree(d_z)); }

}
