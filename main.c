#include "fonctions.h"

int main(int argc, char **argv)
{
	int ligne = 10, colonne = ligne;
	double *matrice = malloc(ligne * colonne * sizeof(double));

	double m_a = 1.0, m_b = 2.0;
	for(int i = 0; i < ligne; i++) 
	{
		for(int j = 0; j < colonne; j++)
		{
			if(i == j)
				matrice[i * ligne + j] = m_a;
			else if (i == j + 1 || i == j - 1)
				matrice[i * ligne + j] = m_b;
			else
				matrice[i * ligne + j] = 0.0;
		}
	}

	fprintf(stdout, "Matrice de base\n");
/*	matrice_test(ligne, colonne, matrice);*/
/*	affichage(ligne, colonne, matrice);*/

	simultaneous_iteration(ligne, colonne, matrice);

/*	comparaison(ligne);*/
	free(matrice);

	return 0;
}
