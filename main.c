#include "fonctions.h"

#ifdef __i386
extern __inline__ uint64_t rdtsc(void)
{
	uint64_t x;
	__asm__ volatile ("rdtsc" : "=A" (x));
	return x;
}
#elif defined __amd64
extern __inline__ uint64_t rdtsc(void)
{
	uint64_t a, d;
	__asm__ volatile ("rdtsc" : "=a" (a), "=d" (d));
	return (d<<32) | a;
}
#endif

int main(int argc, char **argv)
{
	/* *** DECLARATION DES VARIABLES *** */
	int ligne = 11, colonne = ligne, nb_eigen = 10;
	double *matrice = calloc(ligne * colonne, sizeof(double));
	uint64_t timer_start, timer_end;

	/* *** INITIALISATION DES VARIABLES *** */
	matrice_test(ligne, colonne, matrice);

#ifdef DEBUG
	fprintf(stdout, "Matrice de base\n");
	affichage(ligne, colonne, matrice);
#endif

	/* *** CORPS DU PROGRAMME PRINCIPAL *** */
	timer_start = rdtsc();
	simultaneous_iteration(ligne, colonne, nb_eigen, matrice);
	timer_end = rdtsc();

	comparaison(ligne, colonne, matrice);
	fprintf(stdout, "Temps de l'execution : %e s\n", (timer_end - timer_start) / 2.9E9);

	/* *** LIBERATION DES RESSOURCES *** */
	free(matrice);

	return 0;
}
