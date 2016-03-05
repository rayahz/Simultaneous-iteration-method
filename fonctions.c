/* AUTHORS : Yohan CHATELAIN && Rayhana ZIARA
*  Students in M2 IHPS 2015 / 2016
*  Project regarding a Simultaneous Iteration Method
*  LANGUAGE USED : C
*/

#include "fonctions.h"

/*
*  NAME        : affichage
*  DESCRIPTION : permet d'afficher une matrice ou vecteur
*  IN          : nombre de ligne, nombre de colonne, matrice ou vecteur à afficher
*  OUT         : /
*  DEBUG       : /
*/
void affichage(int M, int N, double *A)
{
	for(int i = 0; i < M; i++)
	{
		for(int j = 0; j < N; j++)
			fprintf(stdout, " %.2lf\t", A[i * N + j]);
		fprintf(stdout, "\n");
	}
}

/*
*  NAME        : copy
*  DESCRIPTION : permet de copier une matrice dans une autre matrice
*  IN          : nombre de ligne, nombre de colonne, matrice destination, matrice source
*  OUT         : /
*  DEBUG       : /
*/
void copy(int ligne, int colonne, double *dest, double *src)
{
	memcpy(dest, src, ligne * colonne * sizeof(double));
}

/*
*  NAME        : norme_Frobeinius
*  DESCRIPTION : permet de calculer la norme de Frobeinus
*  IN          : nombre de ligne, nombre de colonne, matrice à k-1, matrice k
*  OUT         : norme
*  DEBUG       : /
*/
double norme_Frobeinius(int ligne, int colonne, double *Q_old, double *Q)
{
	double resultat = 0.0;

	for(int i = 0; i < ligne * colonne; i++)
		resultat += (Q[i] - Q_old[i]) * (Q[i] - Q_old[i]);

	return sqrt(resultat);
}

/*
*  NAME        : mediane
*  DESCRIPTION : permet de calculer la médiane d'un vecteur
*  IN          : nombre de ligne, vecteur
*  OUT         : mediane
*  DEBUG       : /
*/
double mediane(int n, double *T)
{
	int j;
	double x;

	for(int i = 0; i < n - 1; i++)
	{
		x = T[i];
		j = i;
		while(j > 0 && T[j - 1] > x)
		{
			T[j] = T[j - 1];
			j--;
		}
		T[j] = x;
	}

	return T[n / 2];
}

/*
*  NAME        : simultaneous_iteration
*  DESCRIPTION : permet d'effectuer la methode des iterations simultanees
*  IN          : nombre de ligne, nombre de colonne, nombre de valeurs propres à calculer, matrice
*  OUT         : /
*  DEBUG       : affichage de la matrice Q à chaque itération
*/
void simultaneous_iteration(int ligne, int colonne, int nb_eigen, double *A)
{
	double *W = calloc(ligne * colonne, sizeof(double));
	double *Q = calloc(ligne * colonne, sizeof(double));
	double *Q_old = calloc(ligne * colonne, sizeof(double));
	double *lambda = calloc(nb_eigen, sizeof(double));
	double *med = calloc(ligne, sizeof(double));
	double *err = calloc(nb_eigen, sizeof(double));
	double tau[ligne];

	for(int i = 0; i < nb_eigen; i++)
		Q[i * nb_eigen + colonne] = 1.;

	// DGEQRF computes a QR factorization of a real M-by-N matrix A = Q * R
	LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, ligne, nb_eigen, Q, nb_eigen, tau);

	// DORGQR generates an M-by-N real matrix Q with orthonormal columns
	LAPACKE_dorgqr(LAPACK_ROW_MAJOR, ligne, nb_eigen, nb_eigen, Q, nb_eigen, tau);

	int k = 0;
	do
	{
		// copie de Q dans Q_old
		copy(ligne, nb_eigen, Q_old, Q);

		// W = A * Q^k-1
		matMat(ligne, nb_eigen, colonne, A, Q, W);

		// Q * R = W
		LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, ligne, nb_eigen, W, nb_eigen, tau);

		// W = Q
		LAPACKE_dorgqr(LAPACK_ROW_MAJOR, ligne, nb_eigen, nb_eigen, W, nb_eigen, tau);;

#ifdef DEBUG
		fprintf(stdout, "\nMatrice Q\n");
		affichage(ligne, nb_eigen, W);
#endif

		// copie de W dans Q
		copy(ligne, nb_eigen, Q, W);
		k++;
	} while(norme_Frobeinius(ligne, nb_eigen, Q_old, Q) > 1E-6);

	fprintf(stdout, "\nIterations %d - Norme = %e \n", k - 1, norme_Frobeinius(ligne, nb_eigen, Q_old, Q));

	// W = AQ
	matMat(ligne, nb_eigen, colonne, A, Q, W);

	// calcul des valeurs propres
	for(int j = 0; j < nb_eigen; j++)
	{
		for(int i = 0; i < ligne; i++)
			med[i] = fabs(W[i * nb_eigen + j] / Q[i * nb_eigen + j]);
		lambda[j] = mediane(ligne, med);
	}

	fprintf(stdout, "\nValeurs propres issues de la méthode\n");
	affichage(nb_eigen, 1, lambda);

#ifdef DEBUG
	fprintf(stdout, "\nMatrice Q\n");
	affichage(ligne, nb_eigen, Q);
#endif

	free(err);
	free(lambda);
	free(med);
	free(W);
	free(Q);
	free(Q_old);
}

/*
*  NAME        : matrice_test
*  DESCRIPTION : permet de creer une matrice test tridiagonale
*  IN          : nombre de ligne, nombre de colonne, matrice
*  OUT         : /
*  DEBUG       : /
*/
void matrice_test(int ligne, int colonne, double *A)
{
	for(int i = 1; i <= ligne; i++)
	{
		for(int j = 1; j <= colonne; j++)
		{
			if(j >= i)
				A[(i - 1) * colonne + (j - 1)] = ligne + 1 - j;
			else
				A[(i - 1) * colonne + (j - 1)] = ligne + 1 - i;
		}
	}
}

/*
*  NAME        : matMat
*  DESCRIPTION : permet d'effectuer le produit matriciel
*  IN          : nombre de ligne de A, nombre de colonne de A, nombre de ligne de B, matrice A, matrice B et matrice resultat C
*  OUT         : /
*  DGEMM       : utilise blas au lieu de notre implementation
*/
void matMat(int m, int n, int l, double *A, double *B, double *C)
{
#ifdef DGEMM
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, l, 1.0, A, l, B, n, 0.0, C, n);
#else
	int i, j;
	#pragma omp parallel for private(i,j)
	for(int k = 0; k < n; k++)
	{
		for(i = 0; i < m; i++)
		{
			C[i * n + k] = 0.0;
			for(j = 0; j < l; j++)
				C[i * n + k] += A[i * l + j] * B[j * n + k];
		}
	}
#endif
}

/*
*  NAME        : comparaison
*  DESCRIPTION : permet de verifier les valeurs propres obtenus avec notre implementation à ceux de lapacke
*  IN          : nombre de ligne, nombre de colonne, matrice
*  OUT         : /
*  DEBUG       : /
*/
void comparaison(int ligne, int colonne, double *a)
{
	double *eigenvalues = malloc(ligne * sizeof(double));
	double tmp;

	fprintf(stdout, "\nValeurs propres issues de LAPACk\n");
	LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'L', ligne, a, ligne, eigenvalues);

	for(int i = 0; i < ligne / 2; i++)
	{
		tmp = eigenvalues[i];
		eigenvalues[i] = eigenvalues[ligne - i - 1];
		eigenvalues[ligne - i - 1] = tmp;
	}

	affichage(ligne, 1, eigenvalues);
	free(eigenvalues);
}

/*
*  NAME        : saisie_matrice
*  DESCRIPTION : permet de saisir manuellement une matrice
*  IN          : nombre de ligne, nombre de colonne, matrice
*  OUT         : /
*  DEBUG       : /
*/
void saisie_matrice(int ligne, int colonne, double *a)
{
	fprintf(stdout, "\nSaisie de la matrice\n");

	for(int i = 0; i < ligne; i++)
	{
		for(int j = 0; j < colonne; j++)
		{
			fprintf(stdout, "(%d, %d) ?\t", i, j);
			double valeur;
			int erreur = scanf("%lf", &valeur);

			if(erreur == 0)
				a[i * colonne + j] = 0.;
			else 
				a[i * colonne + j] = valeur;
		}
	}

	fprintf(stdout, "\nVoici la matrice que vous avez saisie\n");
	affichage(ligne, colonne, a);
}

/*
*  NAME        : saisie_proprietes
*  DESCRIPTION : permet de saisir le nombre de lignes, colonnes et valeurs propres à calculer
*  IN          : nombre de ligne, nombre de colonne, nombres de valeurs propres
*  OUT         : /
*  DEBUG       : /
*/
void saisie_proprietes(int *ligne, int *colonne, int *nb_eigen)
{
	fprintf(stdout, "Nombre de ligne ?\t");
	if(scanf("%d", ligne) == 0)
		fprintf(stderr, "Veuillez entrer au moins un nombre\n");

	fprintf(stdout, "Nombre de colonne ?\t");
	if(scanf("%d", colonne) == 0)
		fprintf(stderr, "Veuillez entrer au moins un nombre\n");

	// uniquement les matrices carrees
	if(*ligne != *colonne)
	{
		fprintf(stderr, "La matrice doit être carrée\n");
		abort();
		exit(EXIT_FAILURE);
	}

	fprintf(stdout, "Nombre de valeurs propres ?\t");
	if(scanf("%d", nb_eigen) == 0)
		fprintf(stderr, "Veuillez entrer au moins un nombre\n");

	// nombre de valeurs propres <= nombres de colonnes
	if(*nb_eigen - *colonne > 0)
	{
		fprintf(stderr, "Il ne peut pas y avoir plus de valeurs propres que de colonnes\n");
		abort();
		exit(EXIT_FAILURE);
	}
	
}
