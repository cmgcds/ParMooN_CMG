#ifndef __FALTUNG_H
#define __FALTUNG_H

#include "faltung_hilfe.h"

/* Berechnet die Projektion der Faltung der Funktionen f und g auf die Gitterstruktur von w.
 Von f und g sind die Gitterstrukturen sowie die dazugehoerigen Eintraege gegeben. Von w ist
 nur die Gitterstruktur gegeben.

 Eingabewerte: Funktionen f,g mit Gitterstruktur und Eintraegen (dh. zugehoerigen Folgen)
							Funktion w mit Gitterstruktur ohne Eintraege
 Ausgabewert:	Faltung von f und g, projeziert auf die Gitterstruktur von w */


/* Berechnung der FEPC-Faltung */
void
faltung_fepc(func_p0_p f, func_p0_p g, func_p0_p w, fepc_real_t h, func_p0_p back);

void
faltung_fepc(func_p1_p f, func_p1_p g, func_p1_p w, fepc_real_t h, func_p1_p back);

#endif
