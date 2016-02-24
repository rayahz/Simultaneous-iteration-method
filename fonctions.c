#include "fonctions.h"

/*
  PROCEDURE : MPI_initialize
  DESCRIPTION : permet d'initialiser l'environnement MPI
  IN : nombre d'arguments, vecteur d'arguments et structure info_t
  OUT : /
*/
/*void MPI_initialize(int argc, char **argv, struct info_t *info)*/
/*{*/
/*	int Q, R;*/

/*	MPI_Init(&argc, &argv);*/
/*	MPI_Comm_rank(MPI_COMM_WORLD, &(info->rang));*/
/*	MPI_Comm_size(MPI_COMM_WORLD, &(info->nproc));*/

/*	info->temps = 0.0;*/
/*	Q = n / info->nproc;*/
/*	R = n % info->nproc;*/

/*	if(info->rang < R)*/
/*	{*/
/*		info->nloc = Q + 1;*/
/*		info->ideb = info->rang * (Q + 1);*/
/*		info->ifin = info->ideb + info->nloc;*/
/*	}else{*/
/*		info->nloc = Q;*/
/*		info->ideb = R * (Q+1) + (info->rang - R) * Q;*/
/*		info->ifin = info->ideb + info->nloc;*/
/*	}*/
/*}*/


/*  PROCEDURE : print_time
  DESCRIPTION : permet d'afficher le temps d'execution de chaque processus
  IN : structure info_t et temps d'execution
  OUT : /
*/
/*void print_time(struct info_t *info, double runtime)*/
/*{*/
/*	int i;*/
/*	double runtime_in_seconds = runtime / 1000000;*/
/*	double runtime_t[info->nproc];*/

/*	MPI_Gather(&runtime_in_seconds, 1, MPI_DOUBLE, &runtime_t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);*/

/*	if(info->rang == 0)*/
/*	{*/
/*		for(i = 0; i < info->nproc; i++)*/
/*		{*/
/*			if(info->temps < runtime_t[i])*/
/*				info->temps = runtime_t[i];*/
/*				*/
/*			fprintf(stdout, "[%d] Temps de l'execution : %e s\n", i, runtime_t[i]);*/
/*		}*/
/*	}*/
/*}*/

void affichage(int M, int N, double *A)
{
	for(int i = 0; i < M; i++)
	{
		for(int j = 0; j < N; j++)
			fprintf(stdout, " %.2lf\t", A[i * N + j]);
		fprintf(stdout, "\n");
	}
}

double norme(int ligne, int colonne, int k, double *X)
{
	double resultat = 0.0;

	for(int i = 0; i < ligne; i++)
		resultat += X[i * colonne + k] * X[i * colonne + k];

	return sqrt(resultat);
}

double norme2(int ligne, int colonne, int k, double *Q1, double *Q2)
{
	double resultat = 0.0;
	double *test = calloc(ligne * colonne, sizeof(double));

	for(int i = 0; i < ligne; i++)
		test[i * colonne + i] = 1.;

	for(int i = 0; i < ligne; i++)
	{
			resultat += Q1[i * colonne + k-1] * Q2[i * colonne + k];
	}

	free(test);
	return sqrt(resultat);
}

void matMat(int ligne, int colonne, double *A, double *B, double *Z)
{
	// Z = A * B
	for(int i = 0; i < ligne; i++)
	{
		Z[i] = 0.;
		for(int j = 0; j < colonne; j++)
			Z[i] += A[i * colonne + j] * B[i * colonne + j];
	}
}

void matVec(int ligne, int colonne, int k, double *A, double *B, double *Z)
{
	for(int i = 0; i < ligne; i++)
	{
		for(int j = 0; j < colonne; j++)
			Z[i * colonne + j] = A[i * colonne + j] * B[j * ligne + k];
	}
}

void simultaneous_iteration(int ligne, int colonne, double *A)
{
	double *W = calloc(ligne * colonne, sizeof(double));
	double *Q = calloc(ligne * colonne, sizeof(double));
	double tau[ligne];

	for(int i = 0; i < ligne; i++)
		Q[i * colonne] = 1;

	// DGEQRF computes a QR factorization of a real M-by-N matrix A = Q * R
	LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, ligne, colonne, Q, ligne, tau);

	// DORGQR generates an M-by-N real matrix Q with orthonormal columns
	LAPACKE_dorgqr(LAPACK_ROW_MAJOR, ligne, colonne, 1, Q, ligne, tau);

	int k = 0;
	while(norme(ligne, colonne, k, Q) > 0.001 || k < ligne)
	{
		// W = A * Q^k-1
		matVec(ligne, colonne, k, A, Q, W);

		// W = Q * R
		LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, ligne, colonne, W, ligne, tau); 

#ifdef DEBUG
		if(k == ligne - 1)
		{
			fprintf(stdout, "\nMatrice R\n");
			affichage(ligne, colonne, W);
		}
#endif

		// W = Q
		LAPACKE_dorgqr(LAPACK_ROW_MAJOR, ligne, colonne, k, W, ligne, tau);

		for(int i = 0; i < ligne * colonne; i++)
			Q[i] = W[i];

		k++;
	}

	fprintf(stdout, "Iterations %d - Norme = %lf \n", k - 1, norme(ligne, colonne, k - 1, Q));

#ifdef DEBUG
	fprintf(stdout, "\nMatrice Q\n");
	affichage(ligne, colonne, Q);
#endif

	free(W);
	free(Q);
}

void matrice_test(int ligne, int colonne, double *A)
{
	for(int i = 0; i < ligne; i++)
	{
		for(int j = 0; j < colonne; j++)
		{
			if(j <= i)
				A[i * colonne + j] = ligne + 1 - j;
			else
				A[i * colonne + j] = ligne + 1 - i;
		}
	}
}

void comparaison(int ligne)
{
	double *lambda = malloc(ligne * sizeof(double));

	for(int i = 0; i < ligne; i++)
		lambda[i] = 1 / (2 - 2 * cos(((2 * i - 1) * M_PI) / (2 * ligne + 1)));

	fprintf(stdout, "\nComparaison\n");
	affichage(ligne, 1, lambda);
	free(lambda);
}
