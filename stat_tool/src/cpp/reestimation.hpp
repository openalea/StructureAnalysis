/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2015 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: reestimation.hpp 18015 2015-04-23 07:04:17Z guedon $
 *
 *       Forum for V-Plants developers:
 *
 *  ----------------------------------------------------------------------------
 *
 *                      GNU General Public Licence
 *
 *       This program is free software; you can redistribute it and/or
 *       modify it under the terms of the GNU General Public License as
 *       published by the Free Software Foundation; either version 2 of
 *       the License, or (at your option) any later version.
 *
 *       This program is distributed in the hope that it will be useful,
 *       but WITHOUT ANY WARRANTY; without even the implied warranty of
 *       MERCHANTABILITY or FITNESS For A PARTICULAR PURPOSE. See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public
 *       License along with this program; see the file COPYING. If not,
 *       write to the Free Software Foundation, Inc., 59
 *       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  ----------------------------------------------------------------------------
 */



#ifndef REESTIMATION_HPP
#define REESTIMATION_HPP



#include <math.h>

#include <boost/math/special_functions/digamma.hpp>

#include "tool/util_math.h"

#include "stat_label.h"

using namespace std;
using namespace boost::math;


namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Reestimation.
 *
 *  argument : nombre de valeurs.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::init(int inb_value)

{
  nb_element = 0;
  nb_value = inb_value;
  alloc_nb_value = nb_value;
  offset = 0;
  max = 0;
  mean = D_DEFAULT;
  variance = D_DEFAULT;

  if (nb_value == 0) {
    frequency = NULL;
  }

  else {
    register int i;

    frequency = new Type[nb_value];

    for (i = 0;i < nb_value;i++) {
      frequency[i] = 0;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Reestimation.
 *
 *  argument : reference sur un objet Reestimation.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::copy(const Reestimation<Type> &histo)

{
  register int i;


  nb_value = histo.nb_value;
  alloc_nb_value = nb_value;
  nb_element = histo.nb_element;
  offset = histo.offset;
  max = histo.max;
  mean = histo.mean;
  variance = histo.variance;

  frequency = new Type[nb_value];

  for (i = 0;i < nb_value;i++) {
    frequency[i] = histo.frequency[i];
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe Reestimation.
 *
 *  argument : reference sur un objet Reestimation.
 *
 *--------------------------------------------------------------*/

template <typename Type>
Reestimation<Type>::Reestimation(const Reestimation<Type> &histo)

{
  copy(histo);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par fusion de la classe Reestimation.
 *
 *  arguments : nombre d'objets Reestimation, pointeur sur les objets Reestimation.
 *
 *--------------------------------------------------------------*/

template <typename Type>
Reestimation<Type>::Reestimation(int nb_histo , const Reestimation<Type> **histo)

{
  register int i , j;


  nb_value = 0;
  nb_element = 0;
  for (i = 0;i < nb_histo;i++) {
    if (histo[i]->nb_value > nb_value) {
      nb_value = histo[i]->nb_value;
    }
    nb_element += histo[i]->nb_element;
  }
  alloc_nb_value = nb_value;

  frequency = new Type[nb_value];

  for (i = 0;i < nb_value;i++) {
    frequency[i] = 0;
  }
  for (i = 0;i < nb_histo;i++) {
    for (j = histo[i]->offset;j < histo[i]->nb_value;j++) {
      frequency[j] += histo[i]->frequency[j];
    }
  }

  // calcul des caracteristiques de la loi empirique

  offset_computation();
  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Reestimation.
 *
 *--------------------------------------------------------------*/

template <typename Type>
Reestimation<Type>::~Reestimation()

{
  delete [] frequency;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Reestimation.
 *
 *  argument : reference sur un objet Reestimation.
 *
 *--------------------------------------------------------------*/

template <typename Type>
Reestimation<Type>& Reestimation<Type>::operator=(const Reestimation<Type> &histo)

{
  if (&histo != this) {
    delete [] frequency;

    copy(histo);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des caracteristiques d'un objet Reestimation.
 *
 *  arguments : stream, flag ecriture des parametres de forme, flag commentaire.
 *
 *--------------------------------------------------------------*/

template <typename Type>
ostream& Reestimation<Type>::ascii_characteristic_print(ostream &os , bool shape ,
                                                        bool comment_flag) const

{
  os << STAT_label[STATL_SAMPLE_SIZE] << ": " << nb_element << endl;

  if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
    if (comment_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_MEAN] << ": " << mean << "   "
       << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
       << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance) << endl;

    if ((shape) && (variance > 0.)) {
      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << skewness_computation() << "   "
         << STAT_label[STATL_KURTOSIS_COEFF] << ": " << kurtosis_computation() << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des caracteristiques d'un objet Reestimation pour une variable circulaire.
 *
 *  arguments : stream, flag commentaire.
 *
 *--------------------------------------------------------------*/

template <typename Type>
ostream& Reestimation<Type>::ascii_circular_characteristic_print(ostream &os ,
                                                                 bool comment_flag) const

{
  double mean_direction[4];


  os << STAT_label[STATL_SAMPLE_SIZE] << ": " << nb_element << endl;

  mean_direction_computation(mean_direction);

  if (comment_flag) {
    os << "# ";
  }
  os << STAT_label[STATL_MEAN_DIRECTION] << ": " << mean_direction[3];
  if (mean_direction[2] > 0.) {
    os << "   " << STAT_label[STATL_MEAN_RESULTANT_LENGTH] << ": " << mean_direction[2]
       << "   " << STAT_label[STATL_CIRCULAR_STANDARD_DEVIATION] << ": "
       << 180 * sqrt(-2 * log(mean_direction[2])) / M_PI;
  }
  os << endl;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Visualisation d'un objet Reestimation.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

template <typename Type>
ostream& Reestimation<Type>::print(ostream &os) const

{
  register int i;


  os << endl;
  ascii_characteristic_print(os);

  os << "offset : " << offset;
  if (max > 0) {
    os << "   maximum : " << max;
  }
  os << endl;

  os << "frequencies (" << nb_value << " values) : ";
  for (i = 0;i < nb_value;i++) {
    os << frequency[i] << " ";
  }
  os << endl;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de valeurs d'une loi empirique.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::nb_value_computation()

{
  Type *pfrequency;


  pfrequency = frequency + alloc_nb_value;
  nb_value = alloc_nb_value;

  while ((*--pfrequency == 0) && (nb_value > 0)) {
    nb_value--;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de valeurs de frequence nulle a partir de 0.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::offset_computation()

{
  Type *pfrequency;


  pfrequency = frequency;
  offset = 0;

  while ((*pfrequency++ == 0) && (offset < nb_value - 1)) {
    offset++;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'effectif total d'une loi empirique.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::nb_element_computation()

{
  register int i;


  nb_element = 0;
  for (i = offset;i < nb_value;i++) {
    nb_element += frequency[i];
  }
}


/*--------------------------------------------------------------*
 *
 *  Recherche de la frequence maximum d'une loi empirique.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::max_computation()

{
  register int i;


  max = 0;
  for (i = offset;i < nb_value;i++) {
    if (frequency[i] > max) {
      max = frequency[i];
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la moyenne d'une loi empirique.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::mean_computation()

{
  if (nb_element > 0) {
    register int i;


    mean = 0.;
    for (i = offset;i < nb_value;i++) {
      mean += frequency[i] * i;
    }
    mean /= nb_element;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la variance d'une loi empirique.
 *
 *  argument : flag biais.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::variance_computation(bool bias)

{
  if (mean != D_DEFAULT) {
    variance = 0.;

    if (nb_element > 1) {
      register int i;
      double diff;


      for (i = offset;i < nb_value;i++) {
        diff = i - mean;
        variance += frequency[i] * diff * diff;
      }

      switch (bias) {
      case false :
        variance /= (nb_element - 1);
        break;
      case true :
        variance /= nb_element;
        break;
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'ecart absolu moyen d'une loi empirique.
 *
 *--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::mean_absolute_deviation_computation() const

{
  register int i;
  double mean_absolute_deviation = D_DEFAULT;


  if ((mean != D_DEFAULT) && (nb_element > 0)) {
    mean_absolute_deviation = 0.;
    for (i = offset;i < nb_value;i++) {
      mean_absolute_deviation += frequency[i] * fabs(i - mean);
    }
    mean_absolute_deviation /= nb_element;
  }

  return mean_absolute_deviation;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du log de la moyenne geometrique.
 *
 *--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::log_geometric_mean_computation() const

{
  register int i;
  double log_geometric_mean = D_DEFAULT;


  if (nb_element - frequency[0] > 0) {
    log_geometric_mean = 0.;
    for (i = MAX(offset , 2);i < nb_value;i++) {
      log_geometric_mean += frequency[i] * log(i);
    }
    log_geometric_mean /= (nb_element - frequency[0]);
  }

  return log_geometric_mean;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du coefficient d'asymetrie d'une loi empirique.
 *
 *--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::skewness_computation() const

{
  register int i;
  double skewness = D_INF , diff;


  if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
    skewness = 0.;

    if ((nb_element > 2) && (variance > 0.)) {
      for (i = offset;i < nb_value;i++) {
        diff = i - mean;
        skewness += frequency[i] * diff * diff * diff;
      }

      skewness = skewness * nb_element /
                 ((nb_element - 1) * (double)(nb_element - 2) * pow(variance , 1.5));
    }
  }

  return skewness;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'exces d'applatissement d'une loi empirique :
 *  exces d'applatissement = coefficient d'applatissement - 3.
 *
 *--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::kurtosis_computation() const

{
  register int i;
  double kurtosis = D_INF , diff;


  if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
    if (variance == 0.) {
      kurtosis = -2.;
    }

    else {
      kurtosis = 0.;
      for (i = offset;i < nb_value;i++) {
        diff = i - mean;
        kurtosis += frequency[i] * diff * diff * diff * diff;
      }
      kurtosis = kurtosis / ((nb_element - 1) * variance * variance) - 3.;
    }
  }

  return kurtosis;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la direction moyenne d'une variable circulaire.
 *
 *  argument : pointeur sur la direction moyenne.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::mean_direction_computation(double *mean_direction) const

{
  register int i , j;


  mean_direction[0] = 0.;
  mean_direction[1] = 0.;

  for (i = offset;i < nb_value;i++) {
    mean_direction[0] += frequency[i] * cos(i * M_PI / 180);
    mean_direction[1] += frequency[i] * sin(i * M_PI / 180);
  }

  mean_direction[0] /= nb_element;
  mean_direction[1] /= nb_element;

  mean_direction[2] = sqrt(mean_direction[0] * mean_direction[0] +
                           mean_direction[1] * mean_direction[1]);

  if (mean_direction[2] > 0.) {
    mean_direction[3] = atan(mean_direction[1] / mean_direction[0]);

    if (mean_direction[0] < 0.) {
      mean_direction[3] += M_PI;
    }
    else if (mean_direction[1] < 0.) {
      mean_direction[3] += 2 * M_PI;
    }

    mean_direction[3] *= (180 / M_PI);
  }

  else {
    mean_direction[3] = D_DEFAULT;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la quantite d'information d'une loi empirique.
 *
 *--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::information_computation() const

{
  register int i;
  double information = D_INF;


  if (nb_element > 0) {
    information = 0.;
    for (i = offset;i < nb_value;i++) {
      if (frequency[i] > 0) {
        information += frequency[i] * log((double)frequency[i] / (double)nb_element);
      }
    }
  }

  return information;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance d'une loi discrete pour un echantillon.
 *
 *  argument : reference sur un objet Distribution.
 *
 *--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::likelihood_computation(const Distribution &dist) const

{
  register int i;
  double likelihood = 0.;


  if (nb_element > 0) {
    if ((offset < dist.offset) || (nb_value > dist.nb_value)) {
      likelihood = D_INF;
    }

    else {
      for (i = offset;i < nb_value;i++) {
        if (frequency[i] > 0) {
          if (dist.mass[i] > 0.) {
            likelihood += frequency[i] * log(dist.mass[i]);
          }
          else {
            likelihood = D_INF;
            break;
          }
        }
      }

      if ((likelihood != D_INF) && (dist.complement > 0.)) {
        likelihood -= nb_element * log(1. - dist.complement);
      }
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Estimation d'une loi discrete a partir d'un echantillon.
 *
 *  argument : pointeur sur une loi discrete.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::distribution_estimation(Distribution *dist) const

{
  if (nb_element > 0) {
    register int i;


    dist->offset = offset;
    dist->nb_value = nb_value;

    for (i = 0;i < offset;i++) {
      dist->mass[i] = 0.;
    }
    for (i = offset;i < nb_value;i++) {
      dist->mass[i] = (double)frequency[i] / (double)nb_element;
    }

    dist->cumul_computation();

    dist->max = (double)max / (double)nb_element;
    dist->mean = mean;
    dist->variance_computation();
  }
}


/*--------------------------------------------------------------*
 *
 *  Estimation d'une loi discrete a partir d'un echantillon au sens
 *  d'une vraisemblance penalisee.
 *
 *  arguments : pointeur sur une loi discrete, poids de la penalisation,
 *              type de penalisation (difference 1ere, seconde ou entropie),
 *              type de gestion des effets de bord (zero a l'exterieur du support ou
 *              prolongation de la loi).
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::penalized_likelihood_estimation(Distribution *dist , double weight ,
                                                         int type , double *penalty , int outside) const

{
  if (nb_element > 0) {
    register int i;
    int iter;
    double ratio , inf_ratio , sup_ratio , norm , inf_norm , sup_norm;


    dist->penalty_computation(weight , type , penalty , outside);

#   ifdef DEBUG
    {
      switch (type) {
      case FIRST_DIFFERENCE :
        norm = dist->first_difference_norm_computation();
        break;
      case SECOND_DIFFERENCE :
        norm = dist->second_difference_norm_computation();
        break;
      case ENTROPY :
        norm = dist->information_computation();
        break;
      }

      norm = -norm * weight;
      for (i = dist->offset;i < dist->nb_value;i++) {
        norm += frequency[i] * log(dist->mass[i]);
      }

      cout << "\nPenalizedLikelihood: " << norm << endl;
    }
#   endif

    // calcul de la constante de normalisation

    inf_norm = 0.;

    do {
      inf_norm += 0.05 * nb_element;

      for (i = offset;i < nb_value;i++) {
        if (inf_norm + penalty[i] <= 0.) {
          break;
        }
      }
    }
    while ((i < nb_value) && (inf_norm < nb_element / 2.));

    sup_norm = 2. * nb_element;

    if ((type != ENTROPY) && (nb_element + weight > sup_norm)) {
      sup_norm = nb_element + weight;
    }

#   ifdef DEBUG
    cout << "Initialization : " << inf_norm << " " << sup_norm;
#   endif

    inf_ratio = 0.;
    sup_ratio = 0.;

    for (i = offset;i < nb_value;i++) {
      if (sup_norm + penalty[i] > 0.) {
        inf_ratio += frequency[i] / (sup_norm + penalty[i]);
      }
      if (inf_norm + penalty[i] > 0.) {
        sup_ratio += frequency[i] / (inf_norm + penalty[i]);
      }
    }

    iter = 0;
    do {
      iter++;

      ratio = 0.;
      norm = (inf_norm + sup_norm) / 2.;

      for (i = offset;i < nb_value;i++) {
        if (norm + penalty[i] > 0.) {
          ratio += frequency[i] / (norm + penalty[i]);
        }

        else {

#         ifdef MESSAGE
          cout << "\nRATIO ERROR" << "   " << STAT_label[STATL_ITERATION] << " " << iter
               << "   norm: " << norm << "   ratio: " << ratio << " | " << i << " ("
               << inf_ratio << ", " << sup_ratio << ")" << endl;
#         endif

          break;
        }
      }

      if (i < nb_value) {
        break;
      }

      else {
        if (ratio < 1.) {
          inf_ratio = ratio;
          sup_norm = norm;
        }
        else {
          sup_ratio = ratio;
          inf_norm = norm;
        }
      }
    }
    while ((sup_ratio - inf_ratio > BISECTION_RATIO_THRESHOLD) && (iter < BISECTION_NB_ITER));

    if ((sup_ratio - inf_ratio <= BISECTION_RATIO_THRESHOLD) || (iter >= BISECTION_NB_ITER)) {
      norm = (inf_norm + sup_norm) / 2.;

#     ifdef DEBUG
      cout << "   norm : " << norm << " | " << nb_element << "  (" << iter << ")  ";
#     endif

      // reestimation de la loi

      for (i = 0;i < offset;i++) {
        dist->mass[i] = 0.;
      }
      for (i = offset;i < nb_value;i++) {
        dist->mass[i] = frequency[i] / (norm + penalty[i]);
      }
      for (i = nb_value;i < dist->nb_value;i++) {
        dist->mass[i] = 0.;
      }

      dist->offset_computation();
      dist->nb_value_computation();
      dist->cumul_computation();
      dist->max_computation();
      dist->mean_computation();
      dist->variance_computation();
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une loi binomiale a partir
 *  d'un echantillon.
 *
 *  arguments : loi discrete parametrique, borne inferieure minimum,
 *              flag sur la borne inferieure.
 *
 *--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::binomial_estimation(DiscreteParametric *dist , int min_inf_bound ,
                                               bool min_inf_bound_flag) const

{
  register int i , j;
  int max_inf_bound , inf_bound , sup_bound , sup_sup_bound , est_sup_bound ,
      min_sup_bound , max_sup_bound;
  double shift_mean , probability , likelihood , inf_bound_likelihood , max_likelihood = D_INF;


  if (mean - min_inf_bound >= variance) {

    // calcul de l'intervalle teste sur la borne inferieure

    if (!min_inf_bound_flag) {
      max_inf_bound = min_inf_bound;
    }
    else {
      max_inf_bound = MIN(offset , (int)(mean - variance));
    }

    if (variance == 0.) {
      if (mean <= max_inf_bound) {
        max_likelihood = 0.;
        dist->init((int)mean , (int)mean + 1 , D_DEFAULT , 0.);
      }
      else {
        if (mean < dist->alloc_nb_value) {
          max_likelihood = 0.;
          dist->init(min_inf_bound , (int)mean , D_DEFAULT , 1.);
        }
      }
    }

    else {
      sup_sup_bound = MIN(dist->alloc_nb_value ,
                          (int)(nb_value * SAMPLE_NB_VALUE_COEFF)) - 1;

      // estimation pour chaque bornes inferieure et superieure possibles
      // de la probabilite a partir de la moyenne de l'echantillon

#     ifdef DEBUG
//      cout << "\nplage de recherche : " << min_inf_bound
//           << " | " << max_inf_bound << " (" << sup_sup_bound << ")" << endl;

      register int k = 0;
#     endif

      for (i = min_inf_bound;i <= max_inf_bound;i++) {
        shift_mean = mean - i;

        if (shift_mean > variance) {
          est_sup_bound = (int)round(i + shift_mean * shift_mean / (shift_mean - variance));

          min_sup_bound = MAX(nb_value - 1 , i + 1);
          if (min_sup_bound < MIN(sup_sup_bound , est_sup_bound - SUP_BOUND_MARGIN)) {
            min_sup_bound = MIN(sup_sup_bound , est_sup_bound - SUP_BOUND_MARGIN);
          }
          max_sup_bound = MIN(sup_sup_bound , est_sup_bound + SUP_BOUND_MARGIN);
          if (max_sup_bound < min_sup_bound) {
            max_sup_bound = min_sup_bound;
          }

#         ifdef DEBUG
//          cout << min_sup_bound << " " << max_sup_bound << " | " << est_sup_bound << endl;

          k += max_sup_bound - min_sup_bound + 1;
#         endif

          inf_bound_likelihood = D_INF;
          dist->inf_bound = i;

          for (j = min_sup_bound;j <= max_sup_bound;j++) {
            dist->sup_bound = j;
            dist->probability = (mean - i) / (j - i);

            dist->binomial_computation(1 , 's');
            likelihood = dist->likelihood_computation(*this);

#           ifdef DEBUG
//            cout << i << " " << j << " " << dist->probability << " " << likelihood << endl;
#           endif

            if (likelihood > max_likelihood) {
              max_likelihood = likelihood;
              inf_bound = dist->inf_bound;
              sup_bound = dist->sup_bound;
              probability = dist->probability;
            }

            if (likelihood > inf_bound_likelihood) {
              inf_bound_likelihood = likelihood;
            }
          }

          if (inf_bound_likelihood < max_likelihood) {
            break;
          }
        }

        else {
          break;
        }
      }

      // mise a jour de la loi estimee

      if (max_likelihood != D_INF) {
        dist->init(inf_bound , sup_bound , D_DEFAULT , probability);
      }

#     ifdef DEBUG
//      cout << "\nnombre de cas : "
//           << (max_inf_bound - min_inf_bound + 1) * (2 * SUP_BOUND_MARGIN + 1)
//           << " | nombre de calculs : " << k << endl;
#     endif
    }
  }

  return max_likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une loi de Poisson a partir
 *  d'un echantillon.
 *
 *  arguments : loi discrete parametrique, borne inferieure minimum,
 *              flag sur la borne inferieure, seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::poisson_estimation(DiscreteParametric *dist , int min_inf_bound ,
                                              bool min_inf_bound_flag , double cumul_threshold) const

{
  register int i;
  int max_inf_bound , inf_bound;
  double diff , parameter , likelihood , max_likelihood = D_INF;


  // calcul de l'intervalle teste sur la borne inferieure

  if (!min_inf_bound_flag) {
    max_inf_bound = min_inf_bound;
  }

  else {
    diff = mean - variance;
    if (diff < 0.) {
      diff--;
    }

    min_inf_bound = MAX(min_inf_bound , (int)round(diff) - INF_BOUND_MARGIN);
    max_inf_bound = MIN(offset , (int)round(diff) + INF_BOUND_MARGIN);
  }

  if (((mean - max_inf_bound) * POISSON_RATIO < variance) &&
      ((mean - min_inf_bound) / POISSON_RATIO > variance)) {

    // estimation pour chaque borne inferieure possible du parametre
    // a partir de la moyenne de l'echantillon

#   ifdef DEBUG
//    cout << "\nplage de recherche : " << min_inf_bound
//         << " | " << max_inf_bound << endl;
#   endif

    for (i = max_inf_bound;i >= min_inf_bound;i--) {
      dist->inf_bound = i;
      dist->parameter = mean - i;

      dist->poisson_computation(nb_value , cumul_threshold , 's');
      likelihood = dist->likelihood_computation(*this);

      if (likelihood > max_likelihood) {
        max_likelihood = likelihood;
        inf_bound = dist->inf_bound;
        parameter = dist->parameter;
      }

      else {
        if (likelihood < max_likelihood) {
          break;
        }
      }
    }

    // mise a jour de la loi estimee

    if (max_likelihood != D_INF) {
      dist->init(inf_bound , I_DEFAULT , parameter , D_DEFAULT);
    }

#   ifdef DEBUG
//    cout << "\nnombre de cas : " << (max_inf_bound - min_inf_bound + 1)
//         << " | nombre de calculs : " << (max_inf_bound - i + 1) << endl;
#   endif
  }

  return max_likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une loi binomiale negative
 *  a partir d'un echantillon.
 *
 *  arguments : loi discrete parametrique, borne inferieure minimum,
 *              flag sur la borne inferieure, seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::negative_binomial_estimation(DiscreteParametric *dist , int min_inf_bound ,
                                                        bool min_inf_bound_flag , double cumul_threshold) const

{
  register int i;
  int max_inf_bound , inf_bound;
  double diff , shift_mean , parameter , probability , likelihood , max_likelihood = D_INF;


  // calcul de l'intervalle teste sur la borne inferieure

  if (!min_inf_bound_flag) {
    max_inf_bound = min_inf_bound;
  }

  else {
    diff = mean - variance;
    if (diff < 0.) {
      diff--;
    }

    min_inf_bound = MAX(min_inf_bound , (int)diff + 1);
    max_inf_bound = offset;
  }

  if (mean - max_inf_bound < variance) {

    // estimation pour chaque borne inferieure possible du parametre et
    // de la probabilite a partir de la moyenne et de la variance de l'echantillon

#   ifdef DEBUG
//    cout << "\nplage de recherche : " << min_inf_bound
//         << " | " << max_inf_bound << endl;
#   endif

    for (i = max_inf_bound;i >= min_inf_bound;i--) {
      dist->inf_bound = i;

      shift_mean = mean - i;
      dist->parameter = shift_mean * shift_mean / (variance - shift_mean);
      dist->probability = shift_mean / variance;

#     ifdef DEBUG
//      cout << i << " : " dist->parameter << " | " << dist->probability << endl;
#     endif

      dist->negative_binomial_computation(nb_value , cumul_threshold , 's');
      likelihood = dist->likelihood_computation(*this);

      if (likelihood > max_likelihood) {
        max_likelihood = likelihood;
        inf_bound = dist->inf_bound;
        parameter = dist->parameter;
        probability = dist->probability;
      }

      else {
        if (likelihood < max_likelihood) {
          break;
        }
      }
    }

    // mise a jour de la loi estimee

    if (max_likelihood != D_INF) {
      dist->init(inf_bound , I_DEFAULT , parameter , probability);
    }

#   ifdef DEBUG
//    cout << "\nnombre de cas : " << max_inf_bound - min_inf_bound + 1
//         << " | nombre de calculs : " << (max_inf_bound - i + 1) << endl;
#   endif
  }

  return max_likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une loi discrete elementaire
 *  (binomiale negative, binomiale, Poisson) a partir d'un echantillon.
 *
 *  arguments : loi discrete parametrique, borne inferieure minimum,
 *              flag sur la borne inferieure, seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::parametric_estimation(DiscreteParametric *dist , int min_inf_bound ,
                                                 bool min_inf_bound_flag , double cumul_threshold) const

{
  double likelihood;


  switch (dist->ident) {
  case BINOMIAL :
    likelihood = binomial_estimation(dist , min_inf_bound , min_inf_bound_flag);
    break;
  case POISSON :
    likelihood = poisson_estimation(dist , min_inf_bound , min_inf_bound_flag , cumul_threshold);
    break;
  case NEGATIVE_BINOMIAL :
    likelihood = negative_binomial_estimation(dist , min_inf_bound , min_inf_bound_flag , cumul_threshold);
    break;
  }

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Estimation du type et des parametres d'une loi discrete elementaire
 *  (binomiale negative, binomiale, Poisson) a partir d'un echantillon.
 *
 *  arguments : loi discrete parametrique, borne inferieure minimum,
 *              flag sur la borne inferieure, seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::type_parametric_estimation(DiscreteParametric *dist , int min_inf_bound ,
                                                      bool min_inf_bound_flag , double cumul_threshold) const

{
  double likelihood , max_likelihood;
  DiscreteParametric *bdist;


  bdist = new DiscreteParametric(dist->alloc_nb_value);

  max_likelihood = binomial_estimation(dist , min_inf_bound , min_inf_bound_flag);
  if (max_likelihood != D_INF) {
    dist->ident = BINOMIAL;
  }

# ifdef DEBUG
  max_likelihood = D_INF;  // pour les donnees de tremblements de terre
  likelihood = poisson_estimation(bdist , 0 , false , cumul_threshold);
# endif

  likelihood = poisson_estimation(bdist , min_inf_bound , min_inf_bound_flag , cumul_threshold);
  if (likelihood > max_likelihood) {
    bdist->ident = POISSON;
    max_likelihood = likelihood;
    dist->equal_size_copy(*bdist);
    dist->copy(*bdist);
  }

  likelihood = negative_binomial_estimation(bdist , min_inf_bound , min_inf_bound_flag , cumul_threshold);
  if (likelihood > max_likelihood) {
    bdist->ident = NEGATIVE_BINOMIAL;
    max_likelihood = likelihood;
    dist->equal_size_copy(*bdist);
    dist->copy(*bdist);
  }

  delete bdist;

  return max_likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Estimation du type et des parametres d'une loi discrete elementaire
 *  (binomiale negative, binomiale, Poisson) a partir d'un echantillon.
 *
 *  arguments : borne inferieure minimum, flag sur la borne inferieure,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

template <typename Type>
DiscreteParametric* Reestimation<Type>::type_parametric_estimation(int min_inf_bound , bool flag ,
                                                                   double cumul_threshold) const

{
  double likelihood;
  DiscreteParametric *dist;


  // creation d'un objet DiscreteParametric

  dist = new DiscreteParametric((int)(nb_value * SAMPLE_NB_VALUE_COEFF));

  // estimation des parametres de la loi

  likelihood = type_parametric_estimation(dist , min_inf_bound , flag , cumul_threshold);

  // mise a jour de la loi estimee

  if (likelihood != D_INF) {
    dist->computation(nb_value , cumul_threshold);
  }
  else {
    delete dist;
    dist = NULL;
  }

  return dist;
}


/*--------------------------------------------------------------*
 *
 *  Combinaison des quantites de reestimation de la loi et
 *  de la loi biaisee par la longueur (processus en equilibre).
 *
 *  arguments : pointeurs sur les quantites de reestimation de la loi biaisee
 *              par la longueur, moyenne de la loi.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::equilibrium_process_combination(const Reestimation<Type> *length_bias_reestim ,
                                                         double imean)

{
  if (nb_element + length_bias_reestim->nb_element > 0) {
    register int i;


    for (i = offset;i < nb_value;i++) {
      frequency[i] = (frequency[i] + length_bias_reestim->frequency[i]) *
                     (nb_element + length_bias_reestim->nb_element) /
                     (nb_element + length_bias_reestim->nb_element * i / imean);
    }

    nb_element_computation();
    max_computation();
    mean_computation();
    variance_computation();
  }
}


/*--------------------------------------------------------------*
 *
 *  Estimation d'une loi a partir de quantites de reestimation de la loi et
 *  de la loi biaisee par la longueur (processus en equilibre).
 *
 *  arguments : pointeurs sur les quantites de reestimation de la loi biaisee
 *              par la longueur et sur la loi, moyenne de la loi.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::equilibrium_process_estimation(const Reestimation<Type> *length_bias_reestim ,
                                                        Distribution *dist , double imean) const

{
  if (nb_element + length_bias_reestim->nb_element > 0) {
    register int i;


    for (i = 0;i < offset;i++) {
      dist->mass[i] = 0.;
    }
    for (i = offset;i < nb_value;i++) {
      dist->mass[i] = (frequency[i] + length_bias_reestim->frequency[i]) /
                      (nb_element + length_bias_reestim->nb_element * i / imean);
    }
    for (i = nb_value;i < dist->nb_value;i++) {
      dist->mass[i] = 0.;
    }

    dist->offset_computation();
    dist->nb_value_computation();
    dist->cumul_computation();

    // renormalisation de la loi

    for (i = dist->offset;i < dist->nb_value;i++) {
      dist->mass[i] /= dist->cumul[dist->nb_value - 1];  //BRICE OLIVIER : cumul[nb_value - 1] BECOMES cumul[dist->nb_value - 1]
    }

    dist->cumul_computation();
    dist->max_computation();
    dist->mean_computation();
    dist->variance_computation();
  }
}


/*--------------------------------------------------------------*
 *
 *  Estimation d'une loi a partir de quantites de reestimation de la loi et
 *  de la loi biaisee par la longueur au sens d'une vraisemblance penalisee
 *  (processus en equilibre).
 *
 *  arguments : pointeurs sur les quantites de reestimation de la loi biaisee
 *              par la longueur et sur la loi, moyenne de la loi, poids de la penalisation,
 *              type de penalisation (difference 1ere, seconde ou entropie),
 *              type de gestion des effets de bord (zero a l'exterieur du support ou
 *              prolongation de la loi).
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::penalized_likelihood_equilibrium_process_estimation(const Reestimation<Type> *length_bias_reestim ,
                                                                             Distribution *dist , double imean ,
                                                                             double weight , int type ,
                                                                             double *penalty , int outside) const

{
  if (nb_element + length_bias_reestim->nb_element > 0) {
    register int i;
    int iter;
    double ratio , inf_ratio , sup_ratio , norm , inf_norm , sup_norm;


    dist->penalty_computation(weight , type , penalty , outside);

    // calcul de la constante de normalisation

    inf_norm = 0.;

    do {
      inf_norm += 0.05 * nb_element;

      for (i = offset;i < nb_value;i++) {
        if (inf_norm + length_bias_reestim->nb_element * i / imean + penalty[i] <= 0.) {
          break;
        }
      }
    }
    while ((i < nb_value) && (inf_norm < nb_element / 2.));

    sup_norm = 2. * nb_element;

    if ((type != ENTROPY) && (nb_element + weight > sup_norm)) {
      sup_norm = nb_element + weight;
    }

    inf_ratio = 0.;
    sup_ratio = 0.;

    for (i = offset;i < nb_value;i++) {
      if (sup_norm + length_bias_reestim->nb_element * i / imean + penalty[i] > 0.) {
        inf_ratio += (frequency[i] + length_bias_reestim->frequency[i]) /
                     (sup_norm + length_bias_reestim->nb_element * i / imean + penalty[i]);
      }
      if (inf_norm + length_bias_reestim->nb_element * i / imean + penalty[i] > 0.) {
        sup_ratio += (frequency[i] + length_bias_reestim->frequency[i]) /
                     (inf_norm + length_bias_reestim->nb_element * i / imean + penalty[i]);
      }
    }

    iter = 0;
    do {
      iter++;

      ratio = 0.;
      norm = (inf_norm + sup_norm) / 2.;

      for (i = offset;i < nb_value;i++) {
        if (norm + length_bias_reestim->nb_element * i / imean + penalty[i] > 0.) {
          ratio += (frequency[i] + length_bias_reestim->frequency[i]) /
                   (norm + length_bias_reestim->nb_element * i / imean + penalty[i]);
        }

        else {

#         ifdef MESSAGE
          cout << "\nRATIO ERROR" << "   " << STAT_label[STATL_ITERATION] << " " << iter
               << "   norm: " << norm << "   ratio: " << ratio << " | " << i << " ("
               << inf_ratio << ", " << sup_ratio << ")" << endl;
#         endif

          break;
        }
      }

      if (i < nb_value) {
        break;
      }

      else {
        if (ratio < 1.) {
          inf_ratio = ratio;
          sup_norm = norm;
        }
        else {
          sup_ratio = ratio;
          inf_norm = norm;
        }
      }
    }
    while ((sup_ratio - inf_ratio > BISECTION_RATIO_THRESHOLD) && (iter < BISECTION_NB_ITER));

    if ((sup_ratio - inf_ratio <= BISECTION_RATIO_THRESHOLD) || (iter >= BISECTION_NB_ITER)) {
      norm = (inf_norm + sup_norm) / 2.;

      // reestimation de la loi

      for (i = 0;i < offset;i++) {
        dist->mass[i] = 0.;
      }
      for (i = offset;i < nb_value;i++) {
        dist->mass[i] = (frequency[i] + length_bias_reestim->frequency[i]) /
                        (norm + length_bias_reestim->nb_element * i / imean + penalty[i]);
      }
      for (i = nb_value;i < dist->nb_value;i++) {
        dist->mass[i] = 0.;
      }

      dist->offset_computation();
      dist->nb_value_computation();
      dist->cumul_computation();
      dist->max_computation();
      dist->mean_computation();
      dist->variance_computation();
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Estimation d'une loi d'occupation d'un etat par l'estimateur de Kaplan-Meier.
 *
 *  arguments : pointeurs sur les temps de sejour censures a droite,
 *              sur les quantites de reestimation, sur les fonctions de survie
 *              correspondant aux temps de sejour complets et censures a droite,
 *              flag pour le calcul des caracteristiques des quantites de reestimation.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::state_occupancy_estimation(const Reestimation<Type> *final_run ,
                                                    Reestimation<double> *occupancy_reestim ,
                                                    Type *occupancy_survivor ,
                                                    Type *censored_occupancy_survivor ,
                                                    bool characteristic_computation)

{
  register int i;
  int max_nb_value;
  double hazard_rate , hazard_product;


  if (nb_value > 0) {
    occupancy_survivor[nb_value - 1] = frequency[nb_value - 1];
    for (i = nb_value - 2;i >= 1;i--) {
      occupancy_survivor[i] = occupancy_survivor[i + 1] + frequency[i];
    }
  }

  max_nb_value = MAX(nb_value + 1 , final_run->nb_value);
  for (i = max_nb_value - 1;i >= final_run->nb_value;i--) {
    censored_occupancy_survivor[i] = 0;
  }
  censored_occupancy_survivor[final_run->nb_value - 1] = final_run->frequency[final_run->nb_value - 1];
  for (i = final_run->nb_value - 2;i >= 2;i--) {
    censored_occupancy_survivor[i] = censored_occupancy_survivor[i + 1] + final_run->frequency[i];
  }

  hazard_product = 1.;
  for (i = 1;i < nb_value;i++) {
    if (occupancy_survivor[i] + censored_occupancy_survivor[i + 1] > 0) {
      hazard_rate = (double)frequency[i] / (double)(occupancy_survivor[i] +
                     censored_occupancy_survivor[i + 1]);
      occupancy_reestim->frequency[i] = hazard_rate * hazard_product * (nb_element +
                                        final_run->nb_element);
      hazard_product *= (1. - hazard_rate);
    }

    else {
      break;
    }
  }

  if (characteristic_computation) {
    occupancy_reestim->nb_value = i;
  }

  if (i == nb_value) {
    for (i = nb_value;i < final_run->nb_value - 1;i++) {
      occupancy_reestim->frequency[i] = 0.;
    }
    occupancy_reestim->frequency[final_run->nb_value - 1] = hazard_product *
                                                            (nb_element + final_run->nb_element);

    if ((characteristic_computation) &&
        (occupancy_reestim->frequency[final_run->nb_value - 1] > 0)) {
      occupancy_reestim->nb_value = final_run->nb_value;
    }
  }

  if (characteristic_computation) {
    occupancy_reestim->offset = offset;
    occupancy_reestim->nb_element_computation();
    occupancy_reestim->max_computation();
    occupancy_reestim->mean_computation();
    occupancy_reestim->variance_computation();
  }
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une loi gamma a partir
 *  d'une distribution de frequences empiriques.
 *
 *  arguments : loi continue parametrique, iteration EM.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::gamma_estimation(ContinuousParametric *dist , int iter) const

{
  register int i;
  double log_geometric_mean , diff;
//  double bvariance;


  if (frequency[0] / nb_element > GAMMA_ZERO_FREQUENCY_THRESHOLD) {
    dist->shape = 0;
    dist->scale = D_DEFAULT;
  }

  else {
    if (variance > 0.) {
/*      if (sqrt(variance) < mean * GAMMA_VARIATION_COEFF_THRESHOLD) {
        bvariance = mean * mean * GAMMA_VARIATION_COEFF_THRESHOLD * GAMMA_VARIATION_COEFF_THRESHOLD;
      }
      else {
        bvariance = variance;
      } */
//      dist->shape = mean * mean / variance;
//      dist->scale = variance / mean;

      // Hawang & Huang (2012), Ann. Inst. Statist. Math. 54(4), 840-847

      dist->shape = mean * mean / variance - 1. / (double)nb_element;
      dist->scale = mean / dist->shape;

      if ((dist->shape >= GAMMA_SHAPE_PARAMETER_THRESHOLD) &&
          (nb_element < GAMMA_FREQUENCY_THRESHOLD)) {
        log_geometric_mean = log_geometric_mean_computation();
/*        i = 0;   a revoir

#       ifdef DEBUG
        cout << "\n" << STAT_word[STATW_SHAPE] << " : " << dist->shape << "   "
             << STAT_word[STATW_SCALE] << " : " << dist->scale << endl;
#       endif

        do {
          dist->scale = exp(log_geometric_mean - digamma(dist->shape));
          dist->shape = mean / dist->scale;
          i++;

#         ifdef DEBUG
          cout << STAT_word[STATW_SHAPE] << " : " << dist->shape  << "   "
               << STAT_word[STATW_SCALE] << " : " << dist->scale << endl;
#         endif

        }
        while (i < MIN(GAMMA_ITERATION_FACTOR * iter , GAMMA_MAX_NB_ITERATION)); */

        // approximations Johnson, Kotz & Balakrishnan, Continuous Univariate Distributions, vol. 1, 2nd ed., pp. 361-362

//        dist->shape = mean / (2 * (mean - exp(log_geometric_mean))) - 1./12.;
        diff = log(mean) - log_geometric_mean;
        dist->shape = (1 + sqrt(1 + 4 * diff / 3)) / (4 * diff);
        dist->scale = mean / dist->shape;
      }
    }

    else {
      dist->shape = GAMMA_MIN_SHAPE_PARAMETER;
      dist->scale = GAMMA_DEFAULT_SCALE_PARAMETER;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une loi zero-inflated gamma a partir
 *  d'une distribution de frequences empiriques.
 *
 *  arguments : loi continue parametrique, iteration EM.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::zero_inflated_gamma_estimation(ContinuousParametric *dist , int iter) const

{
  if (frequency[0] / nb_element > GAMMA_ZERO_FREQUENCY_THRESHOLD) {
    dist->zero_probability = 1.;
    dist->shape = D_DEFAULT;
    dist->scale = D_DEFAULT;
  }

  else {
    register int i;
    double bmean , bvariance , diff , log_geometric_mean;


    dist->zero_probability = frequency[0] / nb_element;

    bmean = 0.;
    for (i = MIN(offset , 1);i < nb_value;i++) {
      bmean += frequency[i] * i;
    }
    bmean /= (nb_element - frequency[0]);

    bvariance = 0.;

    if (nb_element - frequency[0] > 0) {
      for (i = MIN(offset , 1);i < nb_value;i++) {
        diff = i - bmean;
        bvariance += frequency[i] * diff * diff;
      }
      bvariance /= (nb_element - frequency[0]);
    }

    if (bvariance > 0.) {
/*      if (sqrt(bvariance) < bmean * GAMMA_VARIATION_COEFF_THRESHOLD) {
        bvariance = bmean * bmean * GAMMA_VARIATION_COEFF_THRESHOLD * GAMMA_VARIATION_COEFF_THRESHOLD;
      } */
//      dist->shape = bmean * bmean / bvariance;
//      dist->scale = bvariance / bmean;

      // Hawang & Huang (2012), Ann. Inst. Statist. Math. 54(4), 840-847

      dist->shape = bmean * bmean / bvariance - 1. / (double)(nb_element - frequency[0]);
      dist->scale = bmean / dist->shape;

      if ((dist->shape >= GAMMA_SHAPE_PARAMETER_THRESHOLD) &&
          (nb_element - frequency[0] < GAMMA_FREQUENCY_THRESHOLD)) {
        log_geometric_mean = log_geometric_mean_computation();
/*        i = 0;   a revoir

#       ifdef DEBUG
        cout << "\n" << STAT_word[STATW_SHAPE] << " : " << dist->shape << "   "
             << STAT_word[STATW_SCALE] << " : " << dist->scale << endl;
#       endif

        do {
          dist->scale = exp(log_geometric_mean - digamma(dist->shape));
          dist->shape = bmean / dist->scale;
          i++;

#         ifdef DEBUG
          cout << STAT_word[STATW_SHAPE] << " : " << dist->shape  << "   "
               << STAT_word[STATW_SCALE] << " : " << dist->scale << endl;
#         endif

        }
        while (i < MIN(GAMMA_ITERATION_FACTOR * iter , GAMMA_MAX_NB_ITERATION)); */

        // approximations Johnson, Kotz & Balakrishnan, Continuous Univariate Distributions, vol. 1, 2nd ed., pp. 361-362

//        dist->shape = bmean / (2 * (bmean - exp(log_geometric_mean))) - 1./12.;
        diff = log(bmean) - log_geometric_mean;
        dist->shape = (1 + sqrt(1 + 4 * diff / 3)) / (4 * diff);
        dist->scale = bmean / dist->shape;
      }
    }

    else {
      dist->shape = GAMMA_MIN_SHAPE_PARAMETER;
      dist->scale = GAMMA_DEFAULT_SCALE_PARAMETER;
    }
  }
}


};  // namespace stat_tool



#endif
