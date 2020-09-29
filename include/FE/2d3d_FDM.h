void ComputeCoordinates(int i, int N, int Nz, double *x, double *y, double *z);

void Build_3D_FDM_Matrix(double h, TFEFunction2D *velocity1, TFEFunction2D *velocity2, 
                         TFEFunction2D *concent_C, double *f_old, 
			 double *velo1, double *velo2);

void Evalute_f_at_outflow(int n_dof_FDM, int Nx, int Nz, double *f);

void Integral_For_Particle_Increase_Term(TFESpace2D *fespace, TFEFunction2D *fefct,
					 int n_dof_FDM, int Nx, int Nz, double *f);
