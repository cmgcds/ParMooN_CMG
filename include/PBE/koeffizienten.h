#ifndef __KOEFFIZIENTEN_H
#define __KOEFFIZIENTEN_H

//#include "funktion.h"
#include "config.h"


typedef struct {
	int  zeilen;									/* Anzahl der Zeilen */
	int  spalten;									/* Anzahl der Spalten */
	fepc_real_t  **a;							/* Matrixglieder a[i][j]; i-te Zeile, j-te Spalte */
} matrix_t;

typedef matrix_t *  matrix_p;

/* Datenstruktur entspicht einer d1 x d2 x d3 Matrix mit den Eintraegen a[i][j][k]
mit 0<=i<d1 , 0<=j<d2 , 0<=k<d3 */
typedef struct {
	int  d1;
	int  d2;
	int  d3;
	fepc_real_t  ***a;				/* Eintraege a[i][j][k] der Matrix */
} matrix3_t;

typedef matrix3_t *  matrix3_p;


/* Initialisieren einer Matrix mit m-Zeilen und n-Spalten */
matrix_p
matrix_new(int m,int n);

void
matrix_del(matrix_p matrix);

/*Initialisieren einer d1 x d2 x d3 Matrix*/
matrix3_p
matrix3_new(int d1,int d2,int d3);

void
matrix3_del(matrix3_p matrix);

/* Algorithmus zur Berechnnung der eindimensionalen Xi-Koeffizienten */
matrix_p
koeffizienten_xi_1dim(int grad);

/* Berechnung der eindimensionalen Gammakoeffizienten mit h = 1.0 und l = 0 */
matrix3_p
koeffizienten_gamma_1dim(int grad);

/* Berechnung der mehrdimensionalen Gamma-Koeffizienten bezueglich des Level level, der Vektoren r, alpha, mu, kappa
 und des Gittervektors mesh. Es gilt 0<=alpha, mu, kappa und -1<=r<=0. */
fepc_real_t
koeffizienten_gamma(int level, int r, fepc_real_t mesh, int alpha, int mu, int kappa, matrix3_p gamma_koef);

#endif



