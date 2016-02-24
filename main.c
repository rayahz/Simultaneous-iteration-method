#include "fonctions.h"

int main(int argc, char **argv)
{
	int ligne = 100, colonne = ligne;
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

/*	matrice_test(ligne, colonne, matrice);*/
/*	comparaison(ligne);*/

#ifdef DEBUG
	fprintf(stdout, "Matrice de base\n");
	affichage(ligne, colonne, matrice);
#endif

	simultaneous_iteration(ligne, colonne, matrice);

	free(matrice);

	return 0;
}
