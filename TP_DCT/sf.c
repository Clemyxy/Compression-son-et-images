/*
 * Le but du shannon-fano dynamique est de ne pas transmettre la table
 * des occurrences. Pour ce faire, on ajoute à la table un symbole ESCAPE
 * qui permet l'ajout d'éléments à la table.
 * Le décompresseur sait qu'après un événement ESCAPE il trouvera
 * la valeur (et non le code) d'un événement à ajouter à la table.
 */


#include "bits.h"
#include "sf.h"

#define VALEUR_ESCAPE 0x7fffffff /* Plus grand entier positif */

struct evenement
 {
  int valeur ;
  int nb_occurrences ;
 } ;

struct shannon_fano
 {
  int nb_evenements ;
  struct evenement evenements[200000] ;
 } ;

/*
 * Allocation des la structure et remplissage des champs pour initialiser
 * le tableau des événements avec l'événement ESCAPE (avec une occurrence).
 */
struct shannon_fano* open_shannon_fano()
{
  struct shannon_fano* tmp;
  ALLOUER(tmp, 1); 
  tmp->nb_evenements = 1;
  tmp->evenements[0].valeur = VALEUR_ESCAPE;
  tmp->evenements[0].nb_occurrences = 1;

  return tmp;
}

/*
 * Fermeture (libération mémoire)
 */
void close_shannon_fano(struct shannon_fano *sf)
{
  free(sf);
}

/*
 * En entrée l'événement (sa valeur, pas son code shannon-fano).
 * En sortie la position de l'événement dans le tableau "evenements"
 * Si l'événement n'est pas trouvé, on retourne la position
 * de l'événement ESCAPE.
 */

static int trouve_position(const struct shannon_fano *sf, int evenement)
{
  for (int i = 0; i < sf->nb_evenements; i++){
    if(evenement == sf->evenements[i].valeur){
      return i;
    }
  }
  for (int i = 0; i < sf->nb_evenements; i++){
    if(VALEUR_ESCAPE == sf->evenements[i].valeur){
      return i;
    }
  }

  return -1;
}

/*
 * Soit le sous-tableau "evenements[position_min..position_max]"
 * Les bornes sont incluses.
 *
 * LES FORTES OCCURRENCES SONT LES PETITS INDICES DU TABLEAU
 *
 * Il faut trouver la séparation.
 * La séparation sera codée comme l'indice le plus fort
 * des fortes occurrences.
 *
 * La séparation minimise la valeur absolue de la différence
 * de la somme des occurrences supérieures et inférieures.
 *
 * L'algorithme (trivial) n'est pas facile à trouver, réfléchissez bien.
 */
static int trouve_separation(const struct shannon_fano *sf
			     , int position_min
			     , int position_max)
{
  int min_tot = 0;
  int max_tot = 0;
  int tot_occ = 0;
  for(int i = 0; i < sf->nb_evenements; ++i) {
    tot_occ += sf->evenements[i].nb_occurrences;
  }

  for(int i = position_min; i < position_max; ++i) {
    min_tot += sf->evenements[i].nb_occurrences;
    max_tot = tot_occ - min_tot;

    if(min_tot - max_tot > 0) {
      return i;
    }
  }
  return position_min; /* pour enlever un warning du compilateur */

}

/*
 * Cette fonction (simplement itérative)
 * utilise "trouve_separation" pour générer les bons bit dans "bs"
 * le code de l'événement "sf->evenements[position]".
 * Coupe tableau en deux, si a gauche 0, si a droite 1 -> refaire jusqu'a obtenir la position a coder -> sh dans bs
 */

static void encode_position(struct bitstream *bs,struct shannon_fano *sf,
		     int position)
{
  int posmin = 0, posmax = sf->nb_evenements - 1;
  int sep;
  while(posmin != posmax){
    sep = trouve_separation(sf,posmin,posmax);
    if(position <= sep){
      put_bit(bs, 0);
      posmax = sep;
    }
    else{
      put_bit(bs, 1);
      posmin = sep + 1;
    }
  }
}

/*
 * Cette fonction incrémente le nombre d'occurrence de
 * "sf->evenements[position]"
 * Puis elle modifie le tableau pour qu'il reste trié par nombre
 * d'occurrence (un simple échange d'événement suffit)
 *
 * Les faibles indices correspondent aux grand nombres d'occurrences
 */

static void incremente_et_ordonne(struct shannon_fano *sf, int position)
{
  int nb = sf->evenements[position].nb_occurrences;
  int i = position - 1;

  while(i >= 0 && nb == sf->evenements[i].nb_occurrences) {
      i--;
  }
  
  sf->evenements[position].nb_occurrences++;
  struct evenement temp = sf->evenements[i+1];
  sf->evenements[i+1] = sf->evenements[position];
  sf->evenements[position] = temp;
}

/*
 * Cette fonction trouve la position de l'événement puis l'encode.
 * Si la position envoyée est celle de ESCAPE, elle fait un "put_bits"
 * de "evenement" pour envoyer le code du nouvel l'événement.
 * Elle termine en appelant "incremente_et_ordonne" pour l'événement envoyé.
 * 
 */
void put_entier_shannon_fano(struct bitstream *bs
			     ,struct shannon_fano *sf, int evenement)
{
  int pos = trouve_position(sf, evenement);
  encode_position(bs, sf, pos);
  if(sf->evenements[pos].valeur == VALEUR_ESCAPE){
    put_bits(bs, sizeof(evenement) * 8, evenement);
    sf->evenements[sf->nb_evenements].valeur = evenement;
    sf->evenements[sf->nb_evenements].nb_occurrences = 1;
    sf->nb_evenements++;
  }
  incremente_et_ordonne(sf,pos);
}

/*
 * Fonction inverse de "encode_position"
 */
static int decode_position(struct bitstream *bs,struct shannon_fano *sf)
{
  int min = 0, max = sf->nb_evenements - 1;
  int s;
 
  while(min != max) {
      s = trouve_separation(sf, min, max);
      if(get_bit(bs)) {
        min = s + 1;
      } else {
        max = s;
      }
  }
  
  return min;
}

/*
 * Fonction inverse de "put_entier_shannon_fano"
 *
 * Attention au piège : "incremente_et_ordonne" change le tableau
 * donc l'événement trouvé peut changer de position.
 */
int get_entier_shannon_fano(struct bitstream *bs, struct shannon_fano *sf)
{
int p = decode_position(bs, sf);

  int evenement = sf->evenements[p].valeur;
  if(evenement == VALEUR_ESCAPE) {
    evenement = get_bits(bs, 8 * sizeof(int));
    sf->evenements[sf->nb_evenements].valeur = evenement;
    sf->evenements[sf->nb_evenements].nb_occurrences = 1;
    sf->nb_evenements++;
  }

  incremente_et_ordonne(sf, p);
  return evenement;
}

/*
 * Fonctions pour les tests, NE PAS MODIFIER, NE PAS UTILISER.
 */
int sf_get_nb_evenements(struct shannon_fano *sf)
 {
   return sf->nb_evenements ;
 }
void sf_get_evenement(struct shannon_fano *sf, int i, int *valeur, int *nb_occ)
 {
   *valeur = sf->evenements[i].valeur ;
   *nb_occ = sf->evenements[i].nb_occurrences ;
 }
int sf_table_ok(const struct shannon_fano *sf)
 {
  int i, escape ;

  escape = 0 ;
  for(i=0;i<sf->nb_evenements;i++)
    {
    if ( i != 0
        && sf->evenements[i-1].nb_occurrences<sf->evenements[i].nb_occurrences)
	{
	   fprintf(stderr, "La table des événements n'est pas triée\n") ;
	   return(0) ;
	}
    if ( sf->evenements[i].valeur == VALEUR_ESCAPE )
	escape = 1 ;
    }
 if ( escape == 0 )
	{
	   fprintf(stderr, "Pas de ESCAPE dans la table !\n") ;
	   return(0) ;
	}
 return 1 ;
 }