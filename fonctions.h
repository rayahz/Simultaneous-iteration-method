#ifndef FONCTIONS_H
#define FONCTIONS_H

#include <stdlib.h>
#include <stdio.h>
#include <cblas.h>
#include <lapacke.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
//#include <assert.h>

#define M_PI 3.14159265358979323846
#define DBL_EPSILON 0.00001

void affichage(int, int, double*);
void simultaneous_iteration(int, int, int, double*);
void matrice_test(int, int, double*);
void comparaison(int, int, double*);
void mat_AMn(int, int, double*);
void copy(int, int, double *, double *);
double norme_Frobeinius(int, int, double *, double *);
double mediane(int, double *);
void saisie_matrice(int, int, double *);
void saisie_proprietes(int *, int *, int *);

#endif
