#include "bases.h"
#include "matrice.h"
#include "dct.h"

/*
 * La fonction calculant les coefficients de la DCT (et donc de l'inverse)
 * car la matrice de l'inverse DCT est la transposée de la matrice DCT
 *
 * Cette fonction prend beaucoup de temps.
 * il faut que VOUS l'utilisiez le moins possible (UNE SEULE FOIS)
 *
 * FAITES LES CALCULS EN "double"
 *
 * La valeur de Pi est : M_PI
 *
 * Pour ne pas avoir de problèmes dans la suite du TP, indice vos tableau
 * avec [j][i] et non [i][j].
 */

void coef_dct(Matrice *table)
{
	double n = table->height;
	double sqrtN = sqrt(n);
	double sqrtIJ = sqrt(2) / sqrtN;
	for (int j = 0; j < table->height; ++j) {
		for (int i = 0; i < table->width; ++i) {
			if (j == 0)
				table->t[j][i] = 1 / sqrtN;
			else
				table->t[j][i] = sqrtIJ * cos(((2 * i + 1) * j * M_PI / (2 * n)));
		}
	}
}

/*
 * La fonction calculant la DCT ou son inverse.
 *
 * Cette fonction va être appelée très souvent pour faire
 * la DCT du son ou de l'image (nombreux paquets).
 */

void dct(int   inverse,		/* ==0: DCT, !=0 DCT inverse */
	 int nbe,		/* Nombre d'échantillons  */
	 const float *entree,	/* Le son avant transformation (DCT/INVDCT) */
	 float *sortie		/* Le son après transformation */
	 )
{
	static Matrice *dct = NULL;
	static Matrice * tmp = NULL;
	if(dct == NULL){
		dct = allocation_matrice_float(nbe, nbe);
		tmp = allocation_matrice_float(nbe, nbe);
		coef_dct(dct);
		transposition_matrice(dct, tmp);
	}	

	if(inverse) {
		produit_matrice_vecteur(tmp, entree, sortie);
	}
	else {
		produit_matrice_vecteur(dct, entree, sortie);
	}
	
}
