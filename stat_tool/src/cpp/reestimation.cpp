/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2010 CIRAD/INRIA Virtual Plants
 *
 *       File author(s): Y. Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
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



#ifndef REESTIMATION_C
#define REESTIMATION_C



#include "stat_label.h"
#include "tool/util_math.h"

using namespace std;



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
    Type *pfrequency;

    frequency = new Type[nb_value];

    pfrequency = frequency;
    for (i = 0;i < nb_value;i++) {
      *pfrequency++ = 0;
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
  Type *pfrequency , *cfrequency;


  nb_value = histo.nb_value;
  alloc_nb_value = nb_value;
  nb_element = histo.nb_element;
  offset = histo.offset;
  max = histo.max;
  mean = histo.mean;
  variance = histo.variance;

  // copie des frequences

  frequency = new Type[nb_value];

  pfrequency = frequency;
  cfrequency = histo.frequency;
  for (i = 0;i < nb_value;i++) {
    *pfrequency++ = *cfrequency++;
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
  Type *pfrequency , *cfrequency;


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

  pfrequency = frequency;
  for (i = 0;i < nb_value;i++) {
    *pfrequency++ = 0;
  }

  for (i = 0;i < nb_histo;i++) {
    pfrequency = frequency + histo[i]->offset;
    cfrequency = histo[i]->frequency + histo[i]->offset;
    for (j = histo[i]->offset;j < histo[i]->nb_value;j++) {
      *pfrequency++ += *cfrequency++;
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
  Type *pfrequency;


  pfrequency = frequency + offset;
  nb_element = 0;

  for (i = offset;i < nb_value;i++) {
    nb_element += *pfrequency++;
  }
}


/*--------------------------------------------------------------*
 *
 *  Recherche de la valeur de frequence maximum d'une loi empirique.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::max_computation()

{
  register int i;
  Type *pfrequency;


  pfrequency = frequency + offset;
  max = 0;

  for (i = offset;i < nb_value;i++) {
    if (*pfrequency > max) {
      max = *pfrequency;
    }
    pfrequency++;
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
    Type *pfrequency;


    pfrequency = frequency + offset;
    mean = 0.;

    for (i = offset;i < nb_value;i++) {
      mean += *pfrequency++ * i;
    }

    mean /= nb_element;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la variance d'une loi empirique.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::variance_computation()

{
  if (mean != D_DEFAULT) {
    variance = 0.;

    if (nb_element > 1) {
      register int i;
      Type *pfrequency;
      double diff;


      pfrequency = frequency + offset;
      for (i = offset;i < nb_value;i++) {
        diff = i - mean;
        variance += *pfrequency++ * diff * diff;
      }

      variance /= (nb_element - 1);
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
  Type *pfrequency;
  double mean_absolute_deviation = D_DEFAULT;


  if ((mean != D_DEFAULT) && (nb_element > 0)) {
    pfrequency = frequency + offset;
    mean_absolute_deviation = 0.;

    for (i = offset;i < nb_value;i++) {
      mean_absolute_deviation += *pfrequency++ * fabs(i - mean);
    }

    mean_absolute_deviation /= nb_element;
  }

  return mean_absolute_deviation;
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
  Type *pfrequency;
  double skewness = D_INF , diff;


  if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
    skewness = 0.;

    if ((nb_element > 2) && (variance > 0.)) {
      pfrequency = frequency + offset;
      for (i = offset;i < nb_value;i++) {
        diff = i - mean;
        skewness += *pfrequency++ * diff * diff * diff;
      }

      skewness = skewness * nb_element / ((nb_element - 1) * (double)(nb_element - 2) * pow(variance , 1.5));
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
  Type *pfrequency;
  double kurtosis = D_INF , diff;


  if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
    if (variance == 0.) {
      kurtosis = -2.;
    }

    else {
      pfrequency = frequency + offset;
      kurtosis = 0.;

      for (i = offset;i < nb_value;i++) {
        diff = i - mean;
        kurtosis += *pfrequency++ * diff * diff * diff * diff;
      }

      kurtosis = kurtosis / ((nb_element - 1) * variance * variance) - 3.;
    }
  }

  return kurtosis;
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
  Type *pfrequency;
  double information = D_INF;


  if (nb_element > 0) {
    pfrequency = frequency + offset;
    information = 0.;

    for (i = offset;i < nb_value;i++) {
      if (*pfrequency > 0) {
        information += *pfrequency * log((double)*pfrequency / (double)nb_element);
      }
      pfrequency++;
    }
  }

  return information;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance d'une loi donnee pour un echantillon.
 *
 *  argument : reference sur un objet Distribution.
 *
 *--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::likelihood_computation(const Distribution &dist) const

{
  register int i;
  Type *pfrequency;
  double likelihood = 0. , *pmass;


  if (nb_element > 0) {
    if ((offset < dist.offset) || (nb_value > dist.nb_value)) {
      likelihood = D_INF;
    }

    else {
      pfrequency = frequency + offset;
      pmass = dist.mass + offset;

      for (i = offset;i < nb_value;i++) {
        if (*pfrequency > 0) {
          if (*pmass > 0.) {
            likelihood += *pfrequency * log(*pmass);
          }
          else {
            likelihood = D_INF;
            break;
          }
        }

        pfrequency++;
        pmass++;
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
    Type *pfrequency;
    double *pmass;


    dist->offset = offset;
    dist->nb_value = nb_value;

    pmass = dist->mass;
    for (i = 0;i < offset;i++) {
      *pmass++ = 0.;
    }

    pfrequency = frequency + offset;
    for (i = offset;i < nb_value;i++) {
      *pmass++ = (double)*pfrequency++ / (double)nb_element;
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
    Type *pfrequency;
    double ratio , inf_ratio , sup_ratio , norm , inf_norm , sup_norm , *ppenalty , *pmass;


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

      ppenalty = penalty + offset;
      for (i = offset;i < nb_value;i++) {
        if (inf_norm + *ppenalty++ <= 0.) {
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

    pfrequency = frequency + offset;
    ppenalty = penalty + offset;
    inf_ratio = 0.;
    sup_ratio = 0.;

    for (i = offset;i < nb_value;i++) {
      if (sup_norm + *ppenalty > 0.) {
        inf_ratio += *pfrequency / (sup_norm + *ppenalty);
      }
      if (inf_norm + *ppenalty > 0.) {
        sup_ratio += *pfrequency / (inf_norm + *ppenalty);
      }
      pfrequency++;
      ppenalty++;
    }

    iter = 0;
    do {
      iter++;

      pfrequency = frequency + offset;
      ppenalty = penalty + offset;
      ratio = 0.;
      norm = (inf_norm + sup_norm) / 2.;

      for (i = offset;i < nb_value;i++) {
        if (norm + *ppenalty > 0.) {
          ratio += *pfrequency / (norm + *ppenalty);
        }

        else {

#         ifdef MESSAGE
          cout << "\nRATIO ERROR" << "   " << STAT_label[STATL_ITERATION] << " " << iter
               << "   norm: " << norm << "   ratio: " << ratio << " | " << i << " ("
               << inf_ratio << ", " << sup_ratio << ")" << endl;
#         endif

          break;
        }

        pfrequency++;
        ppenalty++;
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

      pmass = dist->mass;
      for (i = 0;i < offset;i++) {
        *pmass++ = 0.;
      }

      pfrequency = frequency + offset;
      ppenalty = penalty + offset;
      for (i = offset;i < nb_value;i++) {
        *pmass++ = *pfrequency++ / (norm + *ppenalty++);
      }

      for (i = nb_value;i < dist->nb_value;i++) {
        *pmass++ = 0.;
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
 *  arguments : loi parametrique, borne inferieure minimum,
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
 *  arguments : loi parametrique, borne inferieure minimum,
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
 *  arguments : loi parametrique, borne inferieure minimum,
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
 *  arguments : loi parametrique, borne inferieure minimum,
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
 *  arguments : loi parametrique, borne inferieure minimum,
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
    Type *pfrequency , *lfrequency;


    pfrequency = frequency + offset;
    lfrequency = length_bias_reestim->frequency + offset;

    for (i = offset;i < nb_value;i++) {
      *pfrequency = (*pfrequency + *lfrequency++) * (nb_element + length_bias_reestim->nb_element) /
                    (nb_element + length_bias_reestim->nb_element * i / imean);
      *pfrequency++;
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
    Type *pfrequency , *lfrequency;
    double *pmass;


    pmass = dist->mass;
    for (i = 0;i < offset;i++) {
      *pmass++ = 0.;
    }

    pfrequency = frequency + offset;
    lfrequency = length_bias_reestim->frequency + offset;
    for (i = offset;i < nb_value;i++) {
      *pmass++ = (*pfrequency++ + *lfrequency++) /
                 (nb_element + length_bias_reestim->nb_element * i / imean);
    }

    for (i = nb_value;i < dist->nb_value;i++) {
      *pmass++ = 0.;
    }

    dist->offset_computation();
    dist->nb_value_computation();
    dist->cumul_computation();

    // renormalisation de la loi

    pmass = dist->mass + dist->offset;
    for (i = dist->offset;i < dist->nb_value;i++) {
      *pmass++ /= dist->cumul[nb_value - 1];
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
    Type *pfrequency , *lfrequency;
    double ratio , inf_ratio , sup_ratio , norm , inf_norm , sup_norm , *ppenalty , *pmass;


    dist->penalty_computation(weight , type , penalty , outside);

    // calcul de la constante de normalisation

    inf_norm = 0.;

    do {
      inf_norm += 0.05 * nb_element;

      ppenalty = penalty + offset;
      for (i = offset;i < nb_value;i++) {
        if (inf_norm + length_bias_reestim->nb_element * i / imean + *ppenalty++ <= 0.) {
          break;
        }
      }
    }
    while ((i < nb_value) && (inf_norm < nb_element / 2.));

    sup_norm = 2. * nb_element;

    if ((type != ENTROPY) && (nb_element + weight > sup_norm)) {
      sup_norm = nb_element + weight;
    }

    pfrequency = frequency + offset;
    lfrequency = length_bias_reestim->frequency + offset;
    ppenalty = penalty + offset;
    inf_ratio = 0.;
    sup_ratio = 0.;

    for (i = offset;i < nb_value;i++) {
      if (sup_norm + length_bias_reestim->nb_element * i / imean + *ppenalty > 0.) {
        inf_ratio += (*pfrequency + *lfrequency) /
                     (sup_norm + length_bias_reestim->nb_element * i / imean + *ppenalty);
      }
      if (inf_norm + length_bias_reestim->nb_element * i / imean + *ppenalty > 0.) {
        sup_ratio += (*pfrequency + *lfrequency) /
                     (inf_norm + length_bias_reestim->nb_element * i / imean + *ppenalty);
      }
      pfrequency++;
      lfrequency++;
      ppenalty++;
    }

    iter = 0;
    do {
      iter++;

      pfrequency = frequency + offset;
      lfrequency = length_bias_reestim->frequency + offset;
      ppenalty = penalty + offset;
      ratio = 0.;
      norm = (inf_norm + sup_norm) / 2.;

      for (i = offset;i < nb_value;i++) {
        if (norm + length_bias_reestim->nb_element * i / imean + *ppenalty > 0.) {
          ratio += (*pfrequency + *lfrequency) /
                   (norm + length_bias_reestim->nb_element * i / imean + *ppenalty);
        }

        else {

#         ifdef MESSAGE
          cout << "\nRATIO ERROR" << "   " << STAT_label[STATL_ITERATION] << " " << iter
               << "   norm: " << norm << "   ratio: " << ratio << " | " << i << " ("
               << inf_ratio << ", " << sup_ratio << ")" << endl;
#         endif

          break;
        }

        pfrequency++;
        lfrequency++;
        ppenalty++;
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

      pmass = dist->mass;
      for (i = 0;i < offset;i++) {
        *pmass++ = 0.;
      }

      pfrequency = frequency + offset;
      lfrequency = length_bias_reestim->frequency + offset;
      ppenalty = penalty + offset;
      for (i = offset;i < nb_value;i++) {
        *pmass++ = (*pfrequency++ + *lfrequency++) /
                   (norm + length_bias_reestim->nb_element * i / imean + *ppenalty++);
      }

      for (i = nb_value;i < dist->nb_value;i++) {
        *pmass++ = 0.;
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
  double hazard_rate , hazard_product , *preestim;
  Type *pfrequency , *psurvivor;


  if (nb_value > 0) {
    psurvivor = occupancy_survivor + nb_value - 1;
    pfrequency = frequency + nb_value - 1;
    *psurvivor = *pfrequency;
    for (i = nb_value - 2;i >= 1;i--) {
      psurvivor--;
      *psurvivor = *(psurvivor + 1) + *--pfrequency;
    }
  }

  max_nb_value = MAX(nb_value + 1 , final_run->nb_value);
  psurvivor = censored_occupancy_survivor + max_nb_value - 1;
  for (i = max_nb_value - 1;i >= final_run->nb_value;i--) {
    *psurvivor-- = 0;
  }
  pfrequency = final_run->frequency + final_run->nb_value - 1;
  *psurvivor = *pfrequency;
  for (i = final_run->nb_value - 2;i >= 2;i--) {
    psurvivor--;
    *psurvivor = *(psurvivor + 1) + *--pfrequency;
  }

  preestim = occupancy_reestim->frequency + 1;
  pfrequency = frequency + 1;
  hazard_product = 1.;
  for (i = 1;i < nb_value;i++) {
    if (occupancy_survivor[i] + censored_occupancy_survivor[i + 1] > 0) {
      hazard_rate = (double)*pfrequency++ / (double)(occupancy_survivor[i] +
                     censored_occupancy_survivor[i + 1]);
      *preestim++ = hazard_rate * hazard_product * (nb_element +
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
      *preestim++ = 0.;
    }
    *preestim = hazard_product * (nb_element + final_run->nb_element);

    if ((characteristic_computation) && (*preestim > 0)) {
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



#endif
