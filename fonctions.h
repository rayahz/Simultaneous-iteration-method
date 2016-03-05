/* AUTHORS : Yohan CHATELAIN && Rayhana ZIARA
*  Students in M2 IHPS 2015 / 2016
*  Project regarding a Simultaneous Iteration Method
*  LANGUAGE USED : C
*/

#ifndef FONCTIONS_H
#define FONCTIONS_H

#include <stdlib.h>
#include <stdio.h>
#include <cblas.h>
#include <lapacke.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

void affichage(int, int, double*);
void simultaneous_iteration(int, int, int, double*);
void matrice_test(int, int, double*);
void comparaison(int, int, double*);
void copy(int, int, double *, double *);
double norme_Frobeinius(int, int, double *, double *);
double mediane(int, double *);
void saisie_matrice(int, int, double *);
void saisie_proprietes(int *, int *, int *);
void matMat(int, int, int, double *, double *, double *);

#endif
