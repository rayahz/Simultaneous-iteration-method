#include "fonctions.h"

void affichage(int M, int N, double *A)
{
	for(int i = 0; i < M; i++)
	{
		for(int j = 0; j < N; j++)
			fprintf(stdout, " %.2lf\t", A[i * N + j]);
		fprintf(stdout, "\n");
	}
}

void copy(int ligne, int colonne, double *dest, double *src)
{
	memcpy(dest, src, ligne * colonne * sizeof(double));
}

double norme_Frobeinius(int ligne, int colonne, double *Q_old, double *Q)
{
	double resultat = 0.0;

	for(int i = 0; i < ligne * colonne; i++)
		resultat += (Q[i] - Q_old[i]) * (Q[i] - Q_old[i]);

	return sqrt(resultat);
}

void matMat(int m, int n, int k, double *A, double *B, double *C)
{
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, A, k, B, n, 0.0, C, n);
}

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
		memcpy(Q_old, Q, ligne * nb_eigen * sizeof(double));

		// W = A * Q^k-1
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, ligne, nb_eigen, colonne, 1.0, A, colonne, Q, nb_eigen, 0.0, W, nb_eigen);

		// Q * R = W
		LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, ligne, nb_eigen, W, nb_eigen, tau);

		// W = Q
		LAPACKE_dorgqr(LAPACK_ROW_MAJOR, ligne, nb_eigen, nb_eigen, W, nb_eigen, tau);;

#ifdef DEBUG
		fprintf(stdout, "\nMatrice Q\n");
		affichage(ligne, nb_eigen, W);
#endif

		copy(ligne, nb_eigen, Q, W);
		k++;
	} while(norme_Frobeinius(ligne, nb_eigen, Q_old, Q) > 1E-6);

	fprintf(stdout, "Iterations %d - Norme = %lf \n", k - 1, norme_Frobeinius(ligne, nb_eigen, Q_old, Q));

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, ligne, nb_eigen, colonne, 1.0, A, colonne, Q, nb_eigen, 0.0, W, nb_eigen);

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

void matrice_test(int ligne, int colonne, double *A)
{
	for(int i = 1; i <= ligne; i++)
	{
		for(int j = 1; j <= colonne; j++)
		{
			if(j >= i)
				A[(i-1) * colonne + (j-1)] = ligne + 1 - j;
			else
				A[(i-1) * colonne + (j-1)] = ligne + 1 - i;
		}
	}
}

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

void mat_AMn(int ligne, int colonne, double *M)
{
	M[0] =  1.0;
	M[1] = -0.1;
	
	for(int i = 1; i < colonne - 1; i++)
	{
		M[i * colonne + i-1] =  0.1;
		M[i * colonne + i  ] =  1.0 + i;
		M[i * colonne + i+1] = -0.1;
	}

	M[colonne * colonne - 2] = 0.1;
	M[colonne * colonne - 1] = colonne;
}
