#ifndef __FALTUNG_HILFE_H
#define __FALTUNG_HILFE_H

#include "koeffizienten.h"
#include "funktion.h"

/* Berechnung der Restriktion siehe Dokumentation */
void
faltung_hilfe_restrict(folgen_vektor_p0_p f, matrix_p xi_koef, folgen_vektor_p0_p back);

/* Funktion zur Berechnung von Gamma[l][l], siehe Dokumentation */
void
faltung_hilfe_lambda(folgen_vektor_p0_p g, int level, fepc_real_t mesh, matrix3_p gamma_koef, folgen_matrix_p0_p back);

/* Funktion zur Berechnung von Gamma[l], siehe Dokumentation */
void
faltung_hilfe_LAMBDA(folgen_matrix_p0_p Gamma, matrix_p xi_koef, folgen_matrix_p0_p back);

/* Funktion zur Berechnung der Gamma Daten fuer den Algorithmus A1 und B1 (Fall l'<=l siehe Dokumentation) */
void
faltung_hilfe_Gamma_bauen_1( func_p0_p g, fepc_real_t mesh, matrix3_p gamma_koef, matrix_p xi_koef, func2_p0_p back );

/* Funktion zur Berechnung der Gamma Daten fuer den Algorithmus A2 (Fall l'>l siehe Dokumentation) */
void
faltung_hilfe_Gamma_bauen_2( func_p0_p g, fepc_real_t mesh, matrix3_p gamma_koef, matrix_p xi_koef, func2_p0_p back );

/* Berechnung der Restriktion siehe Dokumentation */
void
faltung_hilfe_restrict(folgen_vektor_p1_p f, matrix_p xi_koef, folgen_vektor_p1_p back);

/* Funktion zur Berechnung von Gamma[l][l], siehe Dokumentation */
void
faltung_hilfe_lambda(folgen_vektor_p1_p g, int level, fepc_real_t mesh, matrix3_p gamma_koef, folgen_matrix_p1_p back);

/* Funktion zur Berechnung von Gamma[l], siehe Dokumentation */
void
faltung_hilfe_LAMBDA(folgen_matrix_p1_p Gamma, matrix_p xi_koef, folgen_matrix_p1_p back);

/* Funktion zur Berechnung der Gamma Daten fuer den Algorithmus A1 und B1 (Fall l'<=l siehe Dokumentation) */
void
faltung_hilfe_Gamma_bauen_1( func_p1_p g, fepc_real_t mesh, matrix3_p gamma_koef, matrix_p xi_koef, func2_p1_p back );

/* Funktion zur Berechnung der Gamma Daten fuer den Algorithmus A2 (Fall l'>l siehe Dokumentation) */
void
faltung_hilfe_Gamma_bauen_2( func_p1_p g, fepc_real_t mesh, matrix3_p gamma_koef, matrix_p xi_koef, func2_p1_p back );

#endif