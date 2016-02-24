#ifndef FONCTIONS_H
#define FONCTIONS_H

//#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <cblas.h>
#include <lapacke.h>
#include <math.h>
#include <string.h>

#define M_PI 3.14159265358979323846

struct info_t 
{
	int nproc;
	int rang;
	int nloc;
	int ideb;
	int ifin;
	double temps;
};

// permet d'initialiser l'environnement MPI
void MPI_initialize(int, char **, struct info_t *);

// permet d'afficher le temps d'execution de chaque processus
void print_time(struct info_t *, double);

void affichage(int, int, double*);
double norme(int, int, int, double*);
void matMat(int, int, double*, double*, double*);
void matVec(int, int, int, double*, double*, double*);
void simultaneous_iteration(int, int, double*);
void matrice_test(int, int, double*);
void comparaison(int);

#endif
