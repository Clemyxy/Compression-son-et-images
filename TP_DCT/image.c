#include "image.h"
#include "bases.h"
#include "bit.h"
#include <stddef.h>
#include <stdio.h>
#include <string.h>


/*
 * Lecture d'une ligne du fichier.
 * On saute les lignes commençant par un "#" (commentaire)
 * On simplifie en considérant que les lignes ne dépassent pas MAXLIGNE
 */

void lire_ligne(FILE *f, char *ligne) {
    size_t len = MAXLIGNE;
    while (getline(&ligne, &len, f) && ligne[0] == '#');
}

/*
 * Allocation d'une image
 */

struct image* allocation_image(int hauteur, int largeur) {
    struct image* im;
    ALLOUER(im, 1);

    im->hauteur = hauteur;
    im->largeur = largeur;
    ALLOUER(im->pixels, hauteur);
    for(int i = 0; i<hauteur; ++i) {
        ALLOUER(im->pixels[i], largeur);
    }

    return im;
}

/*
 * Libération image
 */

void liberation_image(struct image* image) {
    for (int i = 0; i < image->hauteur; i++)
        free(image->pixels[i]);
    free(image->pixels);
    free(image);
}

/*
 * Allocation et lecture d'un image au format PGM.
 * (L'entête commence par "P5\nLargeur Hauteur\n255\n"
 * Avec des lignes de commentaire possibles avant la dernière.
 */

struct image* lecture_image(FILE *f) {
    int hauteur, largeur;
    char buffer[MAXLIGNE];

    lire_ligne(f, buffer);
    lire_ligne(f, buffer);
    sscanf(buffer, "%d %d\n", &largeur, &hauteur);
    lire_ligne(f, buffer);

    struct image* im = allocation_image(hauteur, largeur);

    for(int i = 0; i < hauteur; i++)
        fread(im->pixels[i], sizeof(unsigned char), im->largeur, f);

    return im;
}

/*
 * Écriture de l'image (toujours au format PGM)
 */

void ecriture_image(FILE *f, const struct image *image) {
    char buffer[MAXLIGNE];
    sprintf(buffer, "P5\n%d %d\n255\n", image->largeur, image->hauteur);
    fwrite(buffer, sizeof(char), strlen(buffer), f);
    
    for(int i = 0; i < image->hauteur; i++)
        for(int j = 0; j < image->largeur; j++)
            fputc(image->pixels[i][j], f);
}
