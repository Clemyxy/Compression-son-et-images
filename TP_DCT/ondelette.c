#include "bases.h"
#include "bitstream.h"
#include "sf.h"
#include "intstream.h"
#include "image.h"
#include "rle.h"
#include "exception.h"
#include "matrice.h"
#include "ondelette.h"

/*
 * Cette fonction effectue UNE SEULE itération d'une ondelette 1D
 * Voici quelques exemples de calculs
 *
 * Entree            Sortie
 * A                   A
 * A B               (A+B)/2 (A-B)/2
 * A B C             (A+B)/2    C    (A-B)/2
 * A B C D           (A+B)/2 (C+D)/2 (A-B)/2 (C-D)/2
 * A B C D E         (A+B)/2 (C+D)/2   E     (A-B)/2 (C-D)/2
 * A B C D E F       (A+B)/2 (C+D)/2 (E+F)/2 (A-B)/2 (C-D)/2 (E-F)/2
 */

void ondelette_1d(const float *entree, float *sortie, int nbe)
{
    for(int i = 0; i < nbe / 2; ++i) {
        *sortie++ = (entree[2 * i] + entree[2 * i + 1]) / 2.f;
    }
    if(nbe % 2 == 1) {
        *sortie++ = entree[nbe - 1];
    }
    for(int i = 0; i < nbe / 2; ++i) {
        *sortie++ = (entree[2 * i] - entree[2 * i + 1]) / 2.f;
    }
}

/*
 * Comme pour la DCT, on applique dans un sens puis dans l'autre.
 *
 * La fonction reçoit "image" et la modifie directement.
 *
 * Vous pouvez allouer plusieurs images intermédiaires pour
 * simplifier la fonction.
 *
 * Par exemple pour une image  3x6
 *    3x6 ondelette horizontale
 *    On transpose, on a donc une image 6x3
 *    6x3 ondelette horizontale
 *    On transpose à nouveau, on a une image 3x6
 *    On ne travaille plus que sur les basses fréquences (moyennes)
 *    On ne travaille donc que sur le haut gauche de l'image de taille 2x3
 *
 * On recommence :
 *    2x3 horizontal   
 *    transposee => 3x2
 *    3x2 horizontal
 *    transposee => 2x3
 *    basse fréquences => 1x2
 *
 * On recommence :
 *    1x2 horizontal
 *    transposee => 2x1
 *    2x1 horizontal (ne fait rien)
 *    transposee => 1x2
 *    basse fréquences => 1x1
 *
 * L'image finale ne contient qu'un seul pixel de basse fréquence.
 * Les autres sont des blocs de plus ou moins haute fréquence.
 * Sur une image 8x8 :
 *
 * M   	F1H  F2H  F2H  F3H  F3H  F3H  F3H
 * F1V 	F1HV F2H  F2H  F3H  F3H  F3H  F3H
 * F2V 	F2V  F2HV F2HV F3H  F3H  F3H  F3H
 * F2V 	F2V  F2HV F2HV F3H  F3H  F3H  F3H
 * F3V 	F3V  F3V  F3V  F3HV F3HV F3HV F3HV
 * F3V 	F3V  F3V  F3V  F3HV F3HV F3HV F3HV
 * F3V 	F3V  F3V  F3V  F3HV F3HV F3HV F3HV
 * F3V 	F3V  F3V  F3V  F3HV F3HV F3HV F3HV
 *
 * La fréquence F2 est plus petite (moins haute) que la fréquence F3
 * F1H  Indique que c'est une fréquence horizontale
 * F1V  Indique que c'est une fréquence verticale
 * F1HV Indique que c'est une fréquence calculée dans les 2 directions
 * 
 */

void ondelette_2d(Matrice *image)
{
    int h = image->height, w = image->width;
    Matrice* tmp0 = allocation_matrice_float(h, w);
    Matrice* tmp1 = allocation_matrice_float(w, h);
    Matrice* tmp2 = allocation_matrice_float(w, h);
    while(h > 1 || w > 1) {
        for(int i = 0; i < h; ++i) {
            ondelette_1d(image->t[i], tmp0->t[i], w);
        }
        transposition_matrice_partielle(tmp0, tmp1, h, w);
        for(int i = 0; i < w; ++i) {
            ondelette_1d(tmp1->t[i], tmp2->t[i], h);
        }
        transposition_matrice_partielle(tmp2, image, w, h);
        h = (h + 1) / 2;
        w = (w + 1) / 2;
    }
    liberation_matrice_float(tmp0);
    liberation_matrice_float(tmp1);
    liberation_matrice_float(tmp2);
}

/*
 * Quantification de l'ondelette.
 * La facteur de qualité initial s'applique à la fréquence la plus haute.
 * Quand on divise la fréquence par 2 on divise qualité par 8
 * tout en restant supérieur à 1.
 * Une qualité de 1 indique que l'on a pas de pertes.
 */

void quantif_ondelette(Matrice *image, float qualite)
{
    int h = image->height, w = image->width;
    while(qualite > 1 && (h > 1 || w > 1)){
        int halfh = h / 2 + 1, halfw = w / 2 + 1;
        for(int i = 0; i < h; ++i){
            for(int j = 0; j < w; ++j){
                if (i > halfh || j > halfw)
                    image->t[i][j] /= qualite;
            }
        }
        h = halfh;
        w = halfw;
        qualite /= 8;
    }
}

/*
 * Sortie des coefficients dans le bonne ordre afin
 * d'être bien compressé par la RLE.
 * Cette fonction n'est pas optimale, elle devrait faire
 * un parcours de Péano sur chacun des blocs.
 */

void codage_ondelette(Matrice *image, FILE *f)
 {
  int j, i ;
  float *t, *pt ;
  struct intstream *entier, *entier_signe ;
  struct bitstream *bs ;
  struct shannon_fano *sf ;
  int hau, lar ;

  /*
   * Conversion de la matrice en une table linéaire
   * Pour pouvoir utiliser la fonction "compresse"
   */
  hau = image->height ;
  lar = image->width ;
  ALLOUER(t, hau*lar) ;
  pt = t ;

  while( hau != 1 || lar != 1 )
    {
      for(j=0; j<hau; j++)
	for(i=0; i<lar; i++)
	  if ( j>=(hau+1)/2 || i>=(lar+1)/2 )
	    *pt++ = image->t[j][i] ;

      hau = (hau+1)/2 ;
      lar = (lar+1)/2 ;
    }
  *pt = image->t[0][0] ;
  /*
   * Compression RLE avec Shannon-Fano
   */
  bs = open_bitstream("-", "w") ;
  sf = open_shannon_fano() ;
  entier = open_intstream(bs, Shannon_fano, sf) ;
  entier_signe = open_intstream(bs, Shannon_fano, sf) ;

  compresse(entier, entier_signe, image->height*image->width, t) ;

  close_intstream(entier) ;
  close_intstream(entier_signe) ;
  close_bitstream(bs) ;
  free(t) ;
 }
  


/*
*******************************************************************************
* Fonctions inverses
*******************************************************************************
*/

void ondelette_1d_inverse(const float *entree, float *sortie, int nbe)
{
    for(int i = 0; i < nbe / 2; ++i) {
        float sum = entree[i];
        float diff;
        if (nbe % 2 == 1)
            diff = entree[nbe / 2 + i + 1]; 
        else
            diff = entree[nbe / 2 + i];
        *sortie++ = sum + diff;
        *sortie++ = sum - diff;
    }
    if(nbe % 2 == 1) {
        *sortie++ = entree[nbe / 2];
    }
}

void ondelette_2d_inverse_recursive(Matrice *image, int h, int w){
    if(h > 1 || w > 1) {
        ondelette_2d_inverse_recursive(image, (h + 1) / 2, (w + 1) / 2);
    }
    Matrice* tmp0 = allocation_matrice_float(image->height, image->width);
    Matrice* tmp1 = allocation_matrice_float(image->width, image->height);
    Matrice* tmp2 = allocation_matrice_float(image->width, image->height);
    for(int i = 0; i < h; ++i) {
        ondelette_1d_inverse(image->t[i], tmp0->t[i], w);
    }
    transposition_matrice_partielle(tmp0, tmp1, h, w);
    for(int i = 0; i < w; ++i) {
        ondelette_1d_inverse(tmp1->t[i], tmp2->t[i], h);
    }
    transposition_matrice_partielle(tmp2, image, w, h);
    liberation_matrice_float(tmp0);
    liberation_matrice_float(tmp1);
    liberation_matrice_float(tmp2);
}
void ondelette_2d_inverse(Matrice *image)
{
   ondelette_2d_inverse_recursive(image, image->height, image->width);
}


void dequantif_ondelette(Matrice *image, float qualite)
{
int h = image->height, w = image->width;
    while(qualite > 1 && (h > 1 || w > 1)){
        int halfh = h / 2 + 1, halfw = w / 2 + 1;
        for(int i = 0; i < h; ++i){
            for(int j = 0; j < w; ++j){
                if (i > halfh || j > halfw)
                    image->t[i][j] *= qualite;
            }
        }
        h = halfh;
        w = halfw;
        qualite /= 8;
    }
}

void decodage_ondelette(Matrice *image, FILE *f)
 {
  int j, i ;
  float *t, *pt ;
  struct intstream *entier, *entier_signe ;
  struct bitstream *bs ;
  struct shannon_fano *sf ;
  int largeur = image->width, hauteur = image->height ;

  /*
   * Decompression RLE avec Shannon-Fano
   */
  ALLOUER(t, hauteur*largeur) ;
  bs = open_bitstream("-", "r") ;
  sf = open_shannon_fano() ;
  entier = open_intstream(bs, Shannon_fano, sf) ;
  entier_signe = open_intstream(bs, Shannon_fano, sf) ;

  decompresse(entier, entier_signe, hauteur*largeur, t) ;

  close_intstream(entier) ;
  close_intstream(entier_signe) ;
  close_bitstream(bs) ;

  /*
   * Met dans la matrice
   */
  pt = t ;
  while( hauteur != 1 || largeur != 1 )
    {
      for(j=0; j<hauteur; j++)
	for(i=0; i<largeur; i++)
	  if ( j>=(hauteur+1)/2 || i>=(largeur+1)/2 )
	      image->t[j][i] = *pt++ ;

      hauteur = (hauteur+1)/2 ;
      largeur = (largeur+1)/2 ;
    }
  image->t[0][0] = *pt++ ;

  free(t) ;
 }
  
/*
 * Programme de test.
 * La ligne suivante compile, compresse et décompresse l'image
 * et affiche la taille compressée.

export QUALITE=1  # Qualité de "quantification"
export SHANNON=1  # Si 1, utilise shannon-fano dynamique
ondelette <DONNEES/bat710.pgm 1 >xxx && ls -ls xxx && ondelette_inv <xxx | xv -

 */

void ondelette_encode_image(float qualite)
 {
  struct image *image ;
  Matrice *im ;
  int i, j ;

  image = lecture_image(stdin) ;
  assert(fwrite(&image->hauteur, 1, sizeof(image->hauteur), stdout)
	 == sizeof(image->hauteur)) ;
  assert(fwrite(&image->largeur, 1, sizeof(image->largeur), stdout)
	 == sizeof(image->largeur));
  assert(fwrite(&qualite       , 1, sizeof(qualite)       , stdout)
	 == sizeof(qualite));

  im = allocation_matrice_float(image->hauteur, image->largeur) ;
  for(j=0; j<image->hauteur; j++)
    for(i=0; i<image->largeur; i++)
      im->t[j][i] = image->pixels[j][i] ;

  fprintf(stderr, "Compression ondelette, image %dx%d\n"
	  , image->largeur, image->hauteur) ;
  ondelette_2d     (im) ;
  fprintf(stderr, "Quantification qualité = %g\n", qualite) ;
  quantif_ondelette(im, qualite) ;
  fprintf(stderr, "Codage\n") ;
  codage_ondelette(im, stdout) ;

  //  affiche_matrice_float(im, image->hauteur, image->largeur) ;
 }

void ondelette_decode_image()
 {
  int hauteur, largeur ;
  float qualite ;
  struct image *image ;
  Matrice *im ;

  assert(fread(&hauteur, 1, sizeof(hauteur), stdin) == sizeof(hauteur)) ;
  assert(fread(&largeur, 1, sizeof(largeur), stdin) == sizeof(largeur)) ;
  assert(fread(&qualite, 1, sizeof(qualite), stdin) == sizeof(qualite)) ;

  im = allocation_matrice_float(hauteur, largeur) ;

  fprintf(stderr, "Décodage\n") ;
  decodage_ondelette(im, stdin ) ;

  fprintf(stderr, "Déquantification qualité = %g\n", qualite) ;
  dequantif_ondelette(im, qualite) ;

  fprintf(stderr, "Décompression ondelette, image %dx%d\n", largeur, hauteur) ;
  ondelette_2d_inverse (im) ;

  //  affiche_matrice_float(im, hauteur, largeur) ;
  image = creation_image_a_partir_de_matrice_float(im) ;
  ecriture_image(stdout, image) ;
 }


