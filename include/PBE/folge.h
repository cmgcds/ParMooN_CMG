#ifndef __FOLGE_H
#define __FOLGE_H

#include "config.h"

void
faltung(fepc_real_t* a, int n_a, fepc_real_t* b, int n_b, fepc_real_t* c);


/************************************* Constant ansatz *******************************************/

typedef struct {
	int lang;					/* Laenge der Folge */
	fepc_real_t *glied;	/* Speicherung d. Folgenwerte durch multidimensionales Array dessen Dimension durch Vektor size repraesentiert wird */
} folge_p0_t;

typedef folge_p0_t*  folge_p0_p;

/* Folge wird initialisiert. Beachte: alle Gliedwerte werden gleich 0.0 gesetzt */
folge_p0_p
folge_new(int lang);

/* Folge wird geloescht */
void
folge_del(folge_p0_p f);

/* Faltung der Folgen f und g wird berechnet ueber die Fouriertransformation */
void
folge_faltung(folge_p0_p f, folge_p0_p g, folge_p0_p back);

/* Gibt die einzelnen Details einer Folge auf dem Bildschirm wieder. Je groesser der Wert info ist,umso
 mehr Informationen werden wiedergegeben. 0, 1, sind moegliche Infowerte. 0 wenig Info, 1 viele Infos */
void
folge_print(folge_p0_p f, int info);

void
folge_print_aram(folge_p0_p f);

/* Ermittelt ob der Vektor r im Traeger der Folge f liegt. Falls ja, dann Rueckgabe des Folgengliedes
 von f an der Stelle r, falls nein Rueckgabe von 0.0 . */
fepc_real_t
folge_glied( int r, folge_p0_p f );

/* Addition von Folgen f und g */
void
folge_add(folge_p0_p f, folge_p0_p g, folge_p0_p back);

/* Rueckgabe einer einer Kopie der Folge f */
void
folge_copy( folge_p0_p f, folge_p0_p back);

/* Uebertragen der gesamten Folgenwerte der Folge f auf eine Folge mit der gleichen Gitterstruktur wie
 die Folge g */
void
folge_projekt(folge_p0_p f, folge_p0_p g, folge_p0_p back);


/************************************* Linear ansatz *******************************************/

typedef struct {
	int lang;					/* Laenge der Folge */
	fepc_real_t *glied;	/* Speicherung d. Folgenwerte durch multidimensionales Array dessen Dimension durch Vektor size repraesentiert wird */
} folge_p1_t;

typedef folge_p1_t *  folge_p1_p;

/* Folge wird initialisiert. Beachte: alle Gliedwerte werden gleich 0.0 gesetzt */
folge_p1_p
folge_p1_new(int lang);

/* Folge wird geloescht */
void
folge_p1_del(folge_p1_p f);

/* Faltung der Folgen f und g wird berechnet ueber die Fouriertransformation */
void
folge_p1_faltung(folge_p1_p f, folge_p1_p g, folge_p1_p back);

/* Gibt die einzelnen Details einer Folge auf dem Bildschirm wieder. Je groesser der Wert info ist,umso
 mehr Informationen werden wiedergegeben. 0, 1, sind moegliche Infowerte. 0 wenig Info, 1 viele Infos */
void
folge_p1_print(folge_p1_p f, int info);

void
folge_p1_print_aram(folge_p1_p f);

/* Ermittelt ob der Vektor r im Traeger der Folge f liegt. Falls ja, dann Rueckgabe des Folgengliedes
 von f an der Stelle r, falls nein Rueckgabe von 0.0 . */
fepc_real_t
folge_p1_glied( int r, folge_p1_p f );

/* Addition von Folgen f und g */
void
folge_p1_add(folge_p1_p f, folge_p1_p g, folge_p1_p back);

/* Rueckgabe einer einer Kopie der Folge f */
void
folge_p1_copy( folge_p1_p f, folge_p1_p back);

/* Uebertragen der gesamten Folgenwerte der Folge f auf eine Folge mit der gleichen Gitterstruktur wie
 die Folge g */
void
folge_p1_projekt(folge_p1_p f, folge_p1_p g, folge_p1_p back);

#endif