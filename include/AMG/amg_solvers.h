/****************************************************************************/
/*                                                                            */
/* File:          amg_solve_solver.h                                                  */
/*                                                                            */
/* Purpose:   header file for amg_solve_solver.c                                  */
/*                                                                            */
/* Author:    Volker John                                                   */
/*            Otto--von--Guericke Universitaet Magdeburg                    */
/*            Institut fuer Analysis und Numerik                            */
/*            Postfach 4120                                                 */
/*            39016 Magdeburg                                               */
/*            email : volker.john@mathematik.uni-magdeburg.de               */
/*                                                                            */
/* History:   1998/08/06 start                                              */
/*                                                                            */
/* Remarks:                                                                       */
/*                                                                            */
/****************************************************************************/

int prepare_ls_solve(AMG_SolverContext *sc,int *A_length, int *B_length, int depth);
int clear_ls_solve(AMG_SolverContext *sc,int depth);
int ls_solve (AMG_VECTOR *x_in, AMG_VECTOR *b_in);

int prepare_cg_solve(AMG_SolverContext *sc,int *A_length, int *B_length, int depth);
int clear_cg_solve(AMG_SolverContext *sc,int depth);
int cg_solve (AMG_VECTOR *x, AMG_VECTOR *b);

int prepare_bcgs_solve(AMG_SolverContext *sc,int *A_length, int *B_length, int depth);
int clear_bcgs_solve(AMG_SolverContext *sc,int depth);
int bcgs_solve (AMG_VECTOR *x, AMG_VECTOR *b);

int prepare_mixed_bcgs_cgs_solve(AMG_SolverContext *sc,int *A_length, int *B_length, int depth);
int clear_mixed_bcgs_cgs_solve(AMG_SolverContext *sc,int depth);
int mixed_bcgs_cgs_solve (AMG_VECTOR *x, AMG_VECTOR *b);

int prepare_gmres_solve(AMG_SolverContext *sc, int *A_length, int *B_length, int depth);
int clear_gmres_solve(AMG_SolverContext *sc, int depth);
int gmres_left_solve (AMG_VECTOR *x, AMG_VECTOR *b);
int gmres_right_solve (AMG_VECTOR *x, AMG_VECTOR *b);

int prepare_gmres_flexible_solve(AMG_SolverContext *sc, int *A_length, int *B_length, int depth);
int clear_gmres_flexible_solve(AMG_SolverContext *sc, int depth);
int gmres_flexible_solve (AMG_VECTOR *x, AMG_VECTOR *b);

int prepare_exact_solve(AMG_SolverContext *sc,int *A_length, int *B_length, int depth);
int clear_exact_solve(AMG_SolverContext *sc);
int exact_solve (AMG_VECTOR *x, AMG_VECTOR *b);

void AMG_GeneratePlaneRotation(double dx, double dy, double *cs, double *sn);
void AMG_ApplyPlaneRotation(double *dx, double *dy, double cs, double sn);
void AMG_Update(AMG_VECTOR *x, int Len_x, int k, AMG_VECTOR **h, 
                AMG_VECTOR *s, AMG_VECTOR **v);

int prepare_lcd_solve(AMG_SolverContext *sc, int *A_length, int *B_length, int depth);
int clear_lcd_solve(AMG_SolverContext *sc, int depth);
int lcd_solve (AMG_VECTOR *x, AMG_VECTOR *b);

