/// THIS FILE IS TO BE INCLUDED IN MIDDLE OF THE TNSE2D ROUTINES 
// TO FACILITATE THE SETUP FOR DO PARAMETERS

// STAND ALONE THIS FILE HAS NO USE ON ITS OWN -- tHE FILE IS JUST A COPY PASTE CONTENT OF THE CODE FROM 
// MAIN CODE BLOCK.


	///////////////////////////////////////////////////////////////////////////////////////////////
////////// -------- REALISATION DATA GENERATION ----------------------------------------- //////
///////////////////////////////////////////////////////////////////////////////////////////////
int N_Realisations  = 5;
double LengthScale  = 10;
double EigenPercent = 0.9999;

    double *org_x_coord     = new double[N_U];
    double *org_y_coord     = new double[N_U];
    double *x_coord         = new double[N_U];
    double *y_coord         = new double[N_U];
    int *mappingArray       = new int[N_U];

    int i=0;
    int N = (2*pow(2,TDatabase::ParamDB->UNIFORM_STEPS  )) + 1;
    for ( int i = 0 ; i < N_U; i++)
    {
        int local_i = i/N;
        int local_j = i%N;
       
        x_coord[i] =  double(1.0/(N-1)) * local_i;
        y_coord[i] =  double(1.0/(N-1)) * local_j;
    }
cout << " End File Read" <<endl;
Velocity_FeSpace->GetDOFPosition(org_x_coord,org_y_coord);

for ( int i=0 ; i < N_U; i++)   // Generated Values
    {  
        // get the generated Value
        double xx = x_coord[i];
        double yy = y_coord[i];
        bool foundFlag = false;

        for ( int j=0 ; j<N_U;j++)  // Actual parmooN Co-ordinates
        {  
            if(abs(xx - org_x_coord[j]) < 1e-10 &&  abs(yy - org_y_coord[j]) < 1e-10 )
            {
                mappingArray[i] = j;
                foundFlag = true;
            }
        }

        if(!foundFlag) cerr<< " DOF NOT FOUND FOR " << i << " position : " << setw(8) << org_x_coord[i]<<setw(8) <<org_y_coord[i] <<endl;
    }

 double* x  =  new double[N_U];
 double* y  =  new double[N_U];

    for ( int i = 0 ; i < N_U; i++ )
    {
        int local_i = i/N;
        int local_j = i%N;

        x[i] =  double(1.0/(N-1)) * local_j;
        y[i] =  double(1.0/(N-1)) * local_i;
    }

double *C = new double[N_U*N_U];  //MATRIX
double *C1 = new double[N_U*N_U];  //MATRIX  - Corelation Matrix 
 double norm = 0;
    for( int i =0  ; i < N_U ; i++ )
    {
        double actual_x = x[i];
        double actual_y = y[i];

        for ( int j=0 ; j < N_U; j++)
        {
            double local_x = x[j];
            double local_y = y[j];

            double r = sqrt( pow((actual_x - local_x),2 ) + pow((actual_y - local_y),2 ));
            
            // CO -Relation
            C[j*N_U + i] = exp ( (- 1.0 * r )/ (LengthScale) );
            C1[j*N_U  + i] = exp ( (- 1.0 * r )/ (LengthScale) );


            //if(TDatabase::ParamDB->WRITE_PS == 0)
            //{
            //    double sig_r1 = exp (-1.0/(1.0 - pow(( 2*actual_x - 1),4) ) )  * exp ( -1.0/ ( 1 - pow(( 2*actual_y - 1),4) ) ) ;
            //    double sig_r2 = exp (-1.0/(1.0 - pow(( 2*local_x - 1),4) ) )  * exp ( -1.0/ ( 1 - pow(( 2*local_y - 1),4) ) ) ; 
            //
            //    // Co Variance
            //    C[j*N_DOF + i] *= sig_r1 * sig_r2 * 5.0;
            //}

           // else if(TDatabase::ParamDB->WRITE_PS == 1)
           // {
                double E = 0.03;
                double disp = 0.0;
                double power = 2;
                double sig_r1 = exp ( - pow( ( 2*actual_x - 1 - disp),power)  / (E) )  / (2*3.14159265359 * sqrt(E))  * exp ( - pow(( 2*actual_x - 1-disp),power)  / (E) )  / (2*3.14159265359 * sqrt(E)) ;
                double sig_r2 = exp ( - pow(( 2*local_x - 1 -disp),power)  / (E) )  / (2*3.14159265359 * sqrt(E))  * exp ( - pow(( 2*local_y - 1-disp),power)  / (E) ) / (2*3.14159265359 * sqrt(E)); 
                // Co Variance
                C[j*N_U+ i] *= sig_r1 * sig_r2 ;
            //}

            //else{
            //    cout << "Error " <<endl;
            //    exit(0);
            //}

            norm += C[j*N + i]*C[j*N + i];
        }

    }

std::ofstream fileo;
    fileo.open("Corelation.txt");

    for ( int i=0 ; i < N_U ; i++)
    {
        for ( int j=0 ; j < N_U ; j++)
        {
            fileo << C1[i*N_U + j] ;
            if(j!= N_U-1 ) fileo<<",";
        }
        fileo<<endl;
    }

    fileo.close();

std::ofstream fileo_r;
    fileo.open("Covarriance.txt");

    for ( int i=0 ; i < N_U ; i++)
    {
        for ( int j=0 ; j < N_U ; j++)
        {
            fileo_r << C[i*N_U + j] ;
            if(j!= N_U-1 ) fileo_r<<",";
        }
        fileo_r<<endl;
    }

    fileo_r.close();

/////////////////////////////SVD//////////////////////////////////////////////
 // Declare SVD parameters
    MKL_INT m1 = N_U, n = N_U, lda = N_U, ldu = N_U, ldvt = N_U, info;
    double superb[std::min(N_U,N_U)-1];

    double* S = new double[N_U];
    double* U = new double[N_U*N_U];
    double* Vt = new double[N_U*N_U];
    cout << " REALISATIONS COMPUTED " <<endl;
    info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'A', 'A', m1, n, C, lda,
                        S, U, ldu, Vt, ldvt, superb );

    cout << endl <<endl;

    if( info > 0 ) {
      printf( "The algorithm computing SVD failed to converge.\n" );
      exit( 1 );
    }
    cout << " REALISATIONS COMPUTED " <<endl;
    int energyVal = 0;
  
    double sumSingularVal = 0;
    for( int i=0;i<N_U;i++) sumSingularVal += S[i];
    double val = 0;
    for( energyVal =0 ; energyVal< N_U; energyVal++)
    {
        val += S[energyVal];
        if(val/sumSingularVal > 0.8) break;
    }

 cout << " MODES : "  << energyVal+1 <<endl;

  int modDim = energyVal+1;
       
    double* Ut = new double[N_U*modDim]();
    double* Z  = new double[N_Realisations*modDim]();

    double* RealizationVector = new double[N_U * N_Realisations]();
   // -------------- Generate Random Number Based on Normal Distribution -------------------------//
    int k=0;
    int skip = N_U - modDim;
    int count =0;

     for ( int i = 0 ; i < N_U*N_U ; i++ )
    {  
        // cout << "i val " << i <<endl;
        if(count < modDim )
        {
            Ut[k] =  U[i];

            count++;
            k++;
        }
        else
        {
            i += skip;
            count = 0;
            i--;
        }
       
    }

 ///////////
 for( int k = 0 ; k < modDim ; k++)
    {  
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution<> d{0,1};

        double* norm1 = new double[N_Realisations];

        for(int n=0; n<N_Realisations; ++n) {
            Z[k*N_Realisations + n] =  S[k] * d(gen);
        }
    }

    cout << " N_Realisations : " << N_Realisations <<endl;
    cout << " MULT START "<<endl;
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,N_U,N_Realisations, modDim , 1.0, Ut,modDim,Z,N_Realisations,0.0,RealizationVector,N_Realisations);
    cout << " MULT DONE "<<endl;
    // printMatrix(RealizationVector, N_DOF,N_Realisations);

    // mkl_dimatcopy('R','T', N_DOF,N_Realisations,1.0,RealizationVector,N_DOF,N_Realisations);
    //cout << " COPY DONE "<<endl;

    cout << " REALISATIONS COMPUTED " <<endl;   


 //////////////////////////////////End of Realization/////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////
/////////// DO - Initialization /////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

double* MeanVector = new double[N_U * 1]();
for(int i=0; i<N_U; ++i) {
  for(int j = 0; j < N_Realisations; ++j){
            MeanVector[i] +=  (RealizationVector[j*N_U+i]/N_Realisations);
        }
}

double* PerturbationVector = new double[N_U * N_Realisations]();
for(int i = 0; i < N_U; ++i){
  for(int j = 0; j < N_Realisations; ++j){
    PerturbationVector[j*N_U+i] = RealizationVector[j*N_U+i] - MeanVector[i];
  }
}
//================================================================================================
/////////////////////////////DO - Initialization SVD//////////////////////////////////////////////
//================================================================================================
 // Declare SVD parameters
    MKL_INT mDO = N_U, nDO = N_Realisations, ldaDO = N_U, lduDO = N_U, ldvtDO = N_Realisations, infoDO;
    double superbDO[std::min(N_U,N_U)-1];

    double* a = new double[N_U * N_Realisations]();
    memcpy(a,PerturbationVector,N_U*N_Realisations*SizeOfDouble);

    double* Sg = new double[N_Realisations];
    double* L = new double[N_U*N_U];
    double* Rt = new double[N_Realisations*N_Realisations];
    cout << " REALISATIONS COMPUTED " <<endl;
    infoDO = LAPACKE_dgesvd( LAPACK_COL_MAJOR, 'S', 'A', mDO, nDO, a, ldaDO,
                        Sg, L, lduDO, Rt, ldvtDO, superbDO );

    cout << endl <<endl;

    if( infoDO > 0 ) {
      printf( "The algorithm computing SVD for DO failed to converge.\n" );
      exit( 1 );
    }
    cout << " DO SVD COMPUTED " <<endl;
//////DO - SVD End///////////////////////////////

///////DO - Subspace dimension calculation //////
    int s = 0;
    double valDO = 0.0;
    double sumSingularValDO = 0;
    for( int i=0;i<N_Realisations;i++) sumSingularValDO += Sg[i];
    while( valDO/sumSingularValDO < 0.99999999)
    {
        valDO += Sg[s];
        s++;
    }

    cout << " SUBSPACE DIMENSION : "  << s+1 <<endl;

  int subDim = s+1;
  ////////Subspace dimension calculated//////////////////

/////Projection Matrix///////////
////////////////////////////////
double* ProjectionVector = new double[N_Realisations * N_Realisations]();
cout << "PROJ VECTOR MULT START "<<endl;
    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,N_Realisations,N_Realisations, N_U , 1.0, PerturbationVector,N_U,L,N_U,0.0,ProjectionVector,N_Realisations);
    cout << "PROJ VECTOR MULT DONE "<<endl;

///// Initialize Coefficient Matrix - First subDim columns of Projection Matrix ////////////////////
double* CoeffVector = new double[N_Realisations * subDim]();
memcpy(CoeffVector, ProjectionVector, N_Realisations*subDim*SizeOfDouble);

////////////Initialize Mode Vector - First subDim columns of Left Singular Vector//////////////////
double* ModeVector = new double[N_U* subDim]();
memcpy(ModeVector, L, N_U*subDim*SizeOfDouble);

////////////////////////////////////////////DO - Initialization Ends//////////////////////////////////////
///////================================================================================//////////////////

double* Cov = new double[subDim* subDim]();
cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,N_Realisations,subDim, subDim, (1.0/(N_Realisations-1)), CoeffVector,subDim,CoeffVector,subDim,0.0,Cov,subDim);

// Assign the Cov Array to the global pointer - Tdatabase
TDatabase::ParamDB->COVARIANCE_MATRIX_DO = Cov;



///cblas_dgemm Usage///////////////////////////////////////////////////////////////////////////////////////
//aplha(A*B) + beta(C)
//1 - CblasColMajor/CblasRowMajor
//2 - CblasTrans/CblasNoTrans - for A
//3 - CblasTrans/CblasNoTrans - for B
// cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
//            m, n, k, alpha, A, k, B, n, beta, C, n);

// The arguments provide options for how Intel MKL performs the operation. In this case:

// CblasRowMajor
//     Indicates that the matrices are stored in row major order, with the elements of each row of the matrix stored contiguously as shown in the figure above.
// CblasNoTrans
//     Enumeration type indicating that the matrices A and B should not be transposed or conjugate transposed before multiplication.
// m, n, k
//     Integers indicating the size of the matrices:
//         A: m rows by k columns
//         B: k rows by n columns
//         C: m rows by n columns
// alpha
//     Real value used to scale the product of matrices A and B
// A
//     Array used to store matrix A
// k
//     Leading dimension of array A, or the number of elements between successive rows (for row major storage) in memory.
//     In the case of this exercise the leading dimension is the same as the number of columns
// B
//     Array used to store matrix B
// n
//     Leading dimension of array B, or the number of elements between successive rows (for row major storage)
//     in memory. In the case of this exercise the leading dimension is the same as the number of columns
// beta
//     Real value used to scale matrix C
// C
//     Array used to store matrix C
// n
//     Leading dimension of array C, or the number of elements between successive rows (for row major storage)
//     in memory. In the case of this exercise the leading dimension is the same as the number of columns

////////////////////////////////////////////////////////////
///////Co-Skewness Matrix

double* M = new double[subDim*subDim*subDim]();
for(int i = 0; i < subDim; i++){
  for(int j = 0; j < subDim; j++){
    for(int k = 0; k < subDim; k++){
      for(int p = 0; p < subDim; p++){

        M[subDim*subDim*i + subDim*j + k] += ((CoeffVector[subDim*i+p]*CoeffVector[subDim*j+p]*CoeffVector[subDim*k+p])/(N_Realisations-1));

      }

    }
  }
}


	
	
	int N_S  = subDim;
	TDatabase::ParamDB->N_ENERGY_MODES = N_S;

	// Now Combine these two arrays To form a Vect Function Array for the same. 
	// U_tilde_t, V_tilde_t
	// u1 , u2
	// U   [ u1 u2 ]
	// Vect Function  = (U,N_U,2);
	// Here , the second component of velocity is not not being initialised, we will make it to zero
	// The second compinent of modeVector -- modevector_2 

	// Initialize mode vector_2 with zeros  ( The second component of velocity )
	double* ModeVector_2 = new double[N_U* subDim]();


	double* FluctuatingVelocityArray = new double[N_S * N_U * 2]();
	// Assemble Fluctuating Velocity Array
	// Fluctuating Velocity Array =    [ u_tilde_T_1 v_tilde_T_1 u_tilde_T_2 v_tilde_T_2  ........ u_tilde_T_S v_tilde_T_S  ]
	
	for(int kk = 0 ; kk < subDim ; kk++)
	{
		for ( int jj=0 ; jj < N_U ; jj++)
		{
			FluctuatingVelocityArray[kk*subDim + jj] = ModeVector[kk*subDim + jj];
			FluctuatingVelocityArray[(kk+1)*subDim  + jj] = ModeVector_2[kk*subDim + jj];
		}
	}

    double* FluctuatingPressureArray = new double[N_P]();