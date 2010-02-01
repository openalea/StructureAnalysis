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



#include <math.h>
#include "stat_tools.h"
#include "distribution.h"
#include "mixture.h"
#include "stat_label.h"

using namespace std;



/*--------------------------------------------------------------*
 *
 *  Calcul de la quantite d'information d'un objet Mixture_data.
 *
 *--------------------------------------------------------------*/

double Mixture_data::information_computation() const

{
  register int i;
  double information , buff;


  information = weight->information_computation();

  if (information != D_INF) {
    for (i = 0;i < nb_component;i++) {
      if (weight->frequency[i] > 0) {
        buff = component[i]->information_computation();

        if (buff != D_INF) {
          information += buff;
        }
        else {
          information = D_INF;
          break;
        }
      }
    }
  }

  return information;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance d'un objet Mixture_data pour un melange de lois.
 *
 *  argument : reference sur un objet Mixture_data.
 *
 *--------------------------------------------------------------*/

double Mixture::likelihood_computation(const Mixture_data &mixt_histo) const

{
  register int i;
  double likelihood , buff;


  likelihood = weight->likelihood_computation(*(mixt_histo.weight));

  if (likelihood != D_INF) {
    for (i = 0;i < mixt_histo.nb_component;i++) {
      buff = component[i]->likelihood_computation(*(mixt_histo.component[i]));

      if (buff != D_INF) {
        likelihood += buff;
      }
      else {
        likelihood = D_INF;
        break;
      }
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'un melange de lois.
 *
 *  arguments : nombre minimum de valeurs de chaque composante,
 *              seuil sur la fonction de repartition,
 *              flag pour calculer les composantes.
 *
 *--------------------------------------------------------------*/

void Mixture::computation(int min_nb_value , double cumul_threshold , bool component_flag)

{
  register int i , j;
  double *pweight , *pmass;


  // calcul de la loi des ponderations

  if (weight->ident != NONPARAMETRIC) {
    weight->computation(1 , cumul_threshold);
  }

  // calcul des composantes

  if (component_flag) {
    for (i = 0;i < nb_component;i++) {
      component[i]->computation(min_nb_value , cumul_threshold);
    }
  }

  // calcul de la loi resultante

  nb_value = 0;
  for (i = 0;i < nb_component;i++) {
    if (component[i]->nb_value > nb_value) {
      nb_value = component[i]->nb_value;
    }
  }

  offset = nb_value;
  for (i = 0;i < nb_component;i++) {
    if (component[i]->offset < offset) {
      offset = component[i]->offset;
    }
  }

  pmass = mass - 1;
  for (i = 0;i < nb_value;i++) {
    pweight = weight->mass;
    *++pmass = 0.;
    for (j = 0;j < nb_component;j++) {
      if (i < component[j]->nb_value) {
        *pmass += *pweight * component[j]->mass[i];
      }
      pweight++;
    }
  }

  cumul_computation();

  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Initialisation d'un melange de lois a partir d'un histogramme.
 *
 *  arguments : reference sur l'histogramme, flags sur les composantes connues,
 *              borne inferieure minimum, flag sur le decalage des composantes.
 *
 *--------------------------------------------------------------*/

void Mixture::init(const Histogram &histo , bool *estimate ,
                   int min_inf_bound , bool component_flag)

{
  register int i , j = -1;
  int nb_element = 0 , threshold , *pfrequency;
  double cumul_weight = 0. , shift_mean;


  pfrequency = histo.frequency;
  for (i = 0;i < nb_component;i++) {

    if (estimate[i]) {
      if ((i > 0) && (component_flag) && (component[i]->ident != BINOMIAL)) {
        component[i]->inf_bound = j;
      }
      else {
        component[i]->inf_bound = min_inf_bound;
      }
    }

    threshold = (int)round(histo.nb_element * (cumul_weight + 0.5 * weight->mass[i]));
    while (nb_element < threshold) {
      nb_element += *pfrequency++;
      j++;
    }

    shift_mean = j - component[i]->inf_bound;
    if (shift_mean < 1.) {
      shift_mean = 1.;
    }
 
    threshold = (int)round(histo.nb_element * (cumul_weight + 0.75 * weight->mass[i]));
    while (nb_element < threshold) {
      nb_element += *pfrequency++;
      j++;
    }

    if (estimate[i]) {
      switch (component[i]->ident) {

      case BINOMIAL : {
        component[i]->sup_bound = histo.nb_value - 1;
        component[i]->probability = shift_mean / (component[i]->sup_bound - component[i]->inf_bound);
        break;
      }

      case POISSON : {
        component[i]->parameter = shift_mean;
        break;
      }

      case NEGATIVE_BINOMIAL : {
        component[i]->parameter = MIXTURE_PARAMETER;
        component[i]->probability = component[i]->parameter / (shift_mean + component[i]->parameter);
        break;
      }
      }
    }
  }

# ifdef DEBUG
  cout << endl;
  for (i = 0;i < nb_component;i++) {
    cout << "weights : " << weight->mass[i] << "  ";
    component[i]->ascii_print(cout);
  }
# endif

}


/*--------------------------------------------------------------*
 *
 *  Calcul des histogrammes correspondant a chacune des composantes
 *  (estimateur EM d'un melange de lois).
 *
 *  arguments : pointeur sur un objet Mixture_data,
 *              effectif theorique de l'histogramme.
 *
 *--------------------------------------------------------------*/

void Mixture::expectation_step(Mixture_data *mixt_histo , int nb_element) const

{
  register int i , j , k;
  int component_index , value_index , *pfrequency , *mfrequency;
  double scale , sum , max_frequency , *pweight , *mmass , **rfrequency;


  scale = (double)nb_element / (double)mixt_histo->nb_element;

  rfrequency = new double*[nb_component];
  for (i = 0;i < nb_component;i++) {
    rfrequency[i] = new double[mixt_histo->nb_value];
  }

  mfrequency = mixt_histo->frequency + mixt_histo->offset;
  mmass = mass + mixt_histo->offset;
  sum = 0.;

  for (i = mixt_histo->offset;i < mixt_histo->nb_value;i++) {
    if ((*mfrequency > 0) && (*mmass > 0.)) {
      pweight = weight->mass;

      // repartition de l'effectif d'une classe entre les histogrammes
      // correspondant a chacune des composantes

      for (j = 0;j < nb_component;j++) {
        pfrequency = mixt_histo->component[j]->frequency + i;

        if ((i >= component[j]->inf_bound) && (i < component[j]->nb_value)) {
          rfrequency[j][i] = scale * *mfrequency * *pweight * component[j]->mass[i] / *mmass;
          *pfrequency = (int)rfrequency[j][i];
          rfrequency[j][i] -= *pfrequency;
          if (rfrequency[j][i] > 0.) {
            sum += rfrequency[j][i];
          }
        }

        else {
          rfrequency[j][i] = 0.;
          *pfrequency = 0;
        }

        pweight++;
      }
    }

    else {
      for (j = 0;j < nb_component;j++) {
        mixt_histo->component[j]->frequency[i] = 0;
      }
    }

    mfrequency++;
    mmass++;
  }

  // prise en compte des arrondis

  for (i = 0;i < (int)round(sum);i++) {
    max_frequency = 0.;
    mfrequency = mixt_histo->frequency + mixt_histo->offset;
    mmass = mass + mixt_histo->offset;

    for (j = mixt_histo->offset;j < mixt_histo->nb_value;j++) {
      if ((*mfrequency > 0) && (*mmass > 0.)) {
        for (k = 0;k < nb_component;k++) {
          if (rfrequency[k][j] > max_frequency) {
            max_frequency = rfrequency[k][j];
            component_index = k;
            value_index = j;
          }
        }
      }

      mfrequency++;
      mmass++;
    }

    rfrequency[component_index][value_index] = 0.;
    (mixt_histo->component[component_index]->frequency[value_index])++;
  }

  for (i = 0;i < nb_component;i++) {
    delete [] rfrequency[i];
  }
  delete [] rfrequency;

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < nb_component;i++) {
    mixt_histo->component[i]->nb_value_computation();
    mixt_histo->component[i]->offset_computation();
    mixt_histo->component[i]->nb_element_computation();
    mixt_histo->component[i]->max_computation();
    mixt_histo->component[i]->mean_computation();
    mixt_histo->component[i]->variance_computation();

#   ifdef DEBUG
    cout << *mixt_histo->component[i];
#   endif

  }

  // mise a jour de l'histogramme des poids

  for (i = 0;i < nb_component;i++) {
    mixt_histo->weight->frequency[i] = mixt_histo->component[i]->nb_element;
  }
  mixt_histo->weight->nb_element = nb_element;
  mixt_histo->weight->max_computation();
}


/*--------------------------------------------------------------*
 *
 *  Correction des histogrammes dans le cas d'une loi binomiale ou
 *  d'une loi de Poisson
 *
 *  arguments : pointeur sur un objet Mixture_data,
 *              flags sur les composantes connues, borne inferieure minimum.
 *
 *--------------------------------------------------------------*/

void Mixture::variance_correction(Mixture_data *mixt_histo ,
                                  bool *estimate , int min_inf_bound) const

{
  register int i;
  double skewness;
  Histogram *pcomponent , *ncomponent;


  for (i = 0;i < nb_component;i++) {
    pcomponent = mixt_histo->component[i];
    if ((estimate[i]) &&
        (((component[i]->ident == BINOMIAL) && (pcomponent->variance > pcomponent->mean - min_inf_bound)) ||
         ((component[i]->ident == POISSON) && (pcomponent->variance > (pcomponent->mean - min_inf_bound) / POISSON_RATIO)))) {

      skewness = pcomponent->skewness_computation();

      if ((i == 0) || ((i < nb_component - 1) && (skewness > 0.))) {
        ncomponent = mixt_histo->component[i + 1];
        while (((component[i]->ident == BINOMIAL) && (pcomponent->variance > pcomponent->mean - min_inf_bound)) ||
               ((component[i]->ident == POISSON) && (pcomponent->variance > (pcomponent->mean - min_inf_bound) / POISSON_RATIO))) {
          (pcomponent->frequency[pcomponent->nb_value - 1])--;
          (ncomponent->frequency[pcomponent->nb_value - 1])++;

          pcomponent->nb_value_computation();
          (pcomponent->nb_element)--;
          pcomponent->mean_computation();
          pcomponent->variance_computation();
        }
      }

      else {
        ncomponent = mixt_histo->component[i - 1];
        while (((component[i]->ident == BINOMIAL) && (pcomponent->variance > pcomponent->mean - min_inf_bound)) ||
               ((component[i]->ident == POISSON) && (pcomponent->variance > (pcomponent->mean - min_inf_bound) / POISSON_RATIO))) {
          (pcomponent->frequency[pcomponent->offset])--;
          (ncomponent->frequency[pcomponent->offset])++;

          pcomponent->offset_computation();
          (pcomponent->nb_element)--;
          pcomponent->mean_computation();
          pcomponent->variance_computation();
        }
      }

      // extraction des caracteristiques des histogrammes

      pcomponent->max_computation();

      ncomponent->nb_value_computation();
      ncomponent->offset_computation();
      ncomponent->nb_element_computation();
      ncomponent->max_computation();
      ncomponent->mean_computation();
      ncomponent->variance_computation();
    }
  }

  // mise a jour de l'histogramme des poids

  for (i = 0;i < nb_component;i++) {
    mixt_histo->weight->frequency[i] = mixt_histo->component[i]->nb_element;
  }
  mixt_histo->weight->max_computation();
}


/*--------------------------------------------------------------*
 *
 *  Test de l'ordre des composantes du melange dans le cas de melanges "heterogenes".
 *
 *--------------------------------------------------------------*/

bool Mixture::component_order_test() const

{
  bool order = true;
  register int i;


  for (i = 1;i < nb_component;i++) {
    if (component[i]->ident != component[0]->ident) {
      break;
    }
  }

  if (i < nb_component) {
    for (i = 1;i < nb_component;i++) {
      if (component[i]->mean < component[i - 1]->mean) {
        order = false;
        break;
      }
    }
  }

  return order;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un melange de lois par l'algorithme EM.
 *
 *  arguments : reference sur un objet Format_error, pointeur sur les composantes connues,
 *              flags sur les composantes connues, borne inferieure minimum du melange,
 *              flag sur la borne inferieure du melange, flag sur les bornes
 *              inferieures des composantes, pas pour l'initialisation des ponderations.
 *
 *--------------------------------------------------------------*/

Mixture* Histogram::mixture_estimation(Format_error &error , const Mixture &imixt ,
                                       bool *estimate , int min_inf_bound , bool mixt_flag ,
                                       bool component_flag , double weight_step) const

{
  bool status = true;
  register int i , j , k;
  int nb_component = imixt.nb_component , inf_bound[MIXTURE_NB_COMPONENT] , sup_bound[MIXTURE_NB_COMPONENT];
  double step , likelihood , previous_likelihood , max_likelihood = D_INF , weight[MIXTURE_NB_COMPONENT] ,
         parameter[MIXTURE_NB_COMPONENT] , probability[MIXTURE_NB_COMPONENT];
  Mixture *mixt;
  Mixture_data *mixt_histo;


  mixt = NULL;
  error.init();

  if ((min_inf_bound < 0) || (min_inf_bound > 1) || (min_inf_bound > offset)) {
    status = false;
    error.update(STAT_error[STATR_MIN_INF_BOUND]);
  }
  if ((weight_step < MIN_WEIGHT_STEP) || (weight_step > MAX_WEIGHT_STEP)) {
    status = false;
    error.update(STAT_error[STATR_WEIGHT_STEP]);
  }
  if ((!component_flag) && (mixt_flag)) {
    status = false;
    error.update(STAT_error[STATR_DISTRIBUTION_FLAG]);
  }

  if (status) {

    // creation d'un objet Mixture

    mixt = new Mixture(imixt , estimate , (int)(nb_value * SAMPLE_NB_VALUE_COEFF));
    mixt->mixture_data = new Mixture_data(*this , nb_component);
    mixt_histo = mixt->mixture_data;

    // estimation pour chaque ponderation initiale possible des parametres
    // des composantes inconnues du melange au sens du maximum de vraisemblance

    for (i = 0;i < nb_component;i++) {
      for (step = weight_step;step < 1. - weight_step + 1.e-2;step += weight_step) {
        likelihood = D_INF;

        // initialisation des ponderations

        mixt->weight->mass[i] = step;
        for (j = 0;j < nb_component;j++) {
          if (j != i) {
            mixt->weight->mass[j] = (1. - step) / (nb_component - 1);
          }
        }

        // initialisation des parametres des composantes du melange 
        // a partir des ponderations initiales

        mixt->init(*this , estimate , min_inf_bound , component_flag);
        mixt->computation(nb_value);

        j = 0;
        do {
          j++;

          // calcul des histogrammes correspondant aux composantes

          mixt->expectation_step(mixt_histo , (int)round(nb_element * MAX(sqrt(mixt->variance) , 1.) * MIXTURE_COEFF));
          mixt->variance_correction(mixt_histo , estimate , min_inf_bound);

          // reestimation des ponderations

          for (k = 0;k < nb_component;k++) {
            mixt->weight->mass[k] = (double)mixt_histo->weight->frequency[k] /
                                    (double)mixt_histo->weight->nb_element;
          }

          // reestimation des parametres des composantes inconnues

          for (k = 0;k < nb_component;k++) {
            if (estimate[k]) {
              if (((k == 0) && (!mixt_flag)) || (!component_flag)) {
                mixt_histo->component[k]->Reestimation<int>::parametric_estimation(mixt->component[k] ,
                                                                                   min_inf_bound , false);
              }

              else {
                mixt_histo->component[k]->Reestimation<int>::parametric_estimation(mixt->component[k] ,
                                                                                   min_inf_bound , true);
              }
            }
          }

          // calcul du melange estime et de la log-vraisemblance correspondante,
          // test de l'ordre des composantes dans le cas de melanges heterogenes

          mixt->computation(nb_value);
          previous_likelihood = likelihood;
          likelihood = mixt->Distribution::likelihood_computation(*this);

          if (!mixt->component_order_test()) {
            likelihood = D_INF;
          }
        }
        while (((likelihood - previous_likelihood) / -likelihood > MIXTURE_LIKELIHOOD_DIFF) &&
               (j < MIXTURE_NB_ITER) && (likelihood != D_INF));

#       ifdef DEBUG
        cout << "\nnumber of iterations : " << j << "  "
             << STAT_label[STATL_LIKELIHOOD] << " : " << likelihood << endl;
        for (j = 0;j < nb_component;j++) {
          cout << "weights : " << mixt->weight->mass[j] << "  ";
          mixt->component[j]->ascii_print(cout);
        }
#       endif

        // mise a jour des parametres optimaux au sens du maximum de vraisemblance

        if (likelihood > max_likelihood) {
          max_likelihood = likelihood;
          for (j = 0;j < nb_component;j++) {
            weight[j] = mixt->weight->mass[j];
            if (estimate[j]) {
              inf_bound[j] = mixt->component[j]->inf_bound;
              sup_bound[j] = mixt->component[j]->sup_bound;
              parameter[j] = mixt->component[j]->parameter;
              probability[j] = mixt->component[j]->probability;
            }
          }
        }
      }

      // cas 2 composantes

      if ((mixt->nb_component == 2) && (i == 0)) {
        i++;
      }
    }

    likelihood = max_likelihood;

    if (max_likelihood != D_INF) {

      // mise a jour du melange estime

      for (i = 0;i < nb_component;i++) {
        mixt->weight->mass[i] = weight[i];
        if (estimate[i]) {
          mixt->component[i]->init(inf_bound[i] , sup_bound[i] ,
                                   parameter[i] , probability[i]);
        }
      }

      mixt->weight->cumul_computation();
      mixt->weight->max_computation();

      mixt->computation(nb_value);

      // mise a jour du nombre de parametres inconnus

      mixt->nb_parameter = nb_component - 1;
      for (i = 0;i < nb_component;i++) {
        if (estimate[i]) {
          mixt->component[i]->nb_parameter_update();
          mixt->nb_parameter += mixt->component[i]->nb_parameter;
          if (((i == 0) && (!mixt_flag)) || ((i > 0) && (!component_flag))) {
            (mixt->component[i]->nb_parameter)--;
            (mixt->nb_parameter)--;
          }
        }
      }

      mixt->expectation_step(mixt_histo , nb_element);
    }

    else {
      delete mixt;
      mixt = NULL;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }
  }

  return mixt;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un melange de lois par l'algorithme EM.
 *
 *  arguments : reference sur un objet Format_error, pointeur sur les composantes connues,
 *              borne inferieure minimum du melange, flag sur la borne inferieure
 *              du melange, flag sur les bornes inferieures des composantes,
 *              pas pour l'initialisation des ponderations.
 *
 *--------------------------------------------------------------*/

Mixture* Histogram::mixture_estimation(Format_error &error , const Mixture &imixt ,
                                       int min_inf_bound , bool mixt_flag ,
                                       bool component_flag , double weight_step) const

{
  bool estimate[MIXTURE_NB_COMPONENT];
  register int i;
  Mixture *mixt;


  for (i = 0;i < imixt.nb_component;i++) {
    estimate[i] = true;
  }

  mixt = mixture_estimation(error , imixt , estimate , min_inf_bound ,
                            mixt_flag , component_flag , weight_step);

  return mixt;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un melange de lois par l'algorithme EM.
 *
 *  arguments : reference sur un objet Format_error, nombre de composantes,
 *              identificateur des composantes, borne inferieure minimum du melange,
 *              flag sur la borne inferieure du melange, flag sur les bornes
 *              inferieures des composantes, pas pour l'initialisation des ponderations.
 *
 *--------------------------------------------------------------*/

Mixture* Histogram::mixture_estimation(Format_error &error , int nb_component , int *ident ,
                                       int min_inf_bound , bool mixt_flag ,
                                       bool component_flag , double weight_step) const

{
  bool estimate[MIXTURE_NB_COMPONENT];
  register int i;
  const Parametric *pcomponent[MIXTURE_NB_COMPONENT];
  Mixture *imixt , *mixt;


  error.init();

  if ((nb_component < 2) || (nb_component > MIXTURE_NB_COMPONENT)) {
    mixt = NULL;
    error.update(STAT_error[STATR_NB_DISTRIBUTION]);
  }

  else {
    for (i = 0;i < nb_component;i++) {
      pcomponent[i] = new Parametric(0 , ident[i]);
      estimate[i] = true;
    }

    imixt = new Mixture(nb_component , pcomponent);

    for (i = 0;i < nb_component;i++) {
      delete pcomponent[i];
    }

    mixt = mixture_estimation(error , *imixt , estimate , min_inf_bound ,
                              mixt_flag , component_flag , weight_step);

    delete imixt;
  }

  return mixt;
}


/*--------------------------------------------------------------*
 *
 *  Estimation du nombre de composantes et des parametres
 *  d'un melange de lois par l'algorithme EM.
 *
 *  arguments : reference sur un objet Format_error, stream, nombres minimum et
 *              maximum de composantes, identificateur des composantes, borne inferieure minimum
 *              du melange, flag sur la borne inferieure du melange, flag sur
 *              les bornes inferieures des composantes, type de penalisation (AIC(c)/BIC(c)),
 *              pas pour l'initialisation des ponderations.
 *
 *--------------------------------------------------------------*/

Mixture* Histogram::mixture_estimation(Format_error &error , std::ostream &os , int min_nb_component ,
                                       int max_nb_component , int *ident , int min_inf_bound ,
                                       bool mixt_flag , bool component_flag , int penalty_type ,
                                       double weight_step) const

{
  bool status = true , estimate[MIXTURE_NB_COMPONENT];
  register int i;
  int nb_parameter[MIXTURE_NB_COMPONENT + 1];
  double penalty , max_likelihood , likelihood[MIXTURE_NB_COMPONENT + 1] ,
         penalized_likelihood[MIXTURE_NB_COMPONENT + 1];
  const Parametric *pcomponent[MIXTURE_NB_COMPONENT];
  Parametric_model *dist;
  Mixture *imixt , *mixt , *pmixt;


  mixt = NULL;
  error.init();

  if ((min_nb_component < 1) || (min_nb_component >= max_nb_component)) {
    status = false;
    error.update(STAT_error[STATR_NB_DISTRIBUTION]);
  }
  if ((max_nb_component > MIXTURE_NB_COMPONENT) || (max_nb_component <= min_nb_component)) {
    status = false;
    error.update(STAT_error[STATR_NB_DISTRIBUTION]);
  }

  if (status) {
    for (i = 0;i < max_nb_component;i++) {
      pcomponent[i] = new Parametric(0 , ident[i]);
      estimate[i] = true;
    }

    // calcul du terme de penalisation

    switch (penalty_type) {
    case AIC :
      penalty = 1.;
      break;
    case BIC :
      penalty = 0.5 * log((double)nb_element);
      break;
    }

    likelihood[min_nb_component] = D_INF;

    if (min_nb_component == 1) {
//      dist = parametric_estimation(error , ident[0] , min_inf_bound , mixt_flag);
      dist = type_parametric_estimation(error , min_inf_bound , mixt_flag);

      if (dist) {
        likelihood[min_nb_component] = dist->likelihood_computation(*this);
        nb_parameter[min_nb_component] = dist->nb_parameter_computation();

        if ((likelihood[min_nb_component] != D_INF) && (penalty_type == BICc)) {
          penalized_likelihood[min_nb_component] = likelihood[min_nb_component] -
                                                   0.5 * nb_parameter[min_nb_component] * log((double)nb_element);
        }

#       ifdef MESSAGE
        os << "\n";
        dist->ascii_print(os);
#       endif

        delete dist;
      }
    }

    else {
      imixt = new Mixture(min_nb_component , pcomponent);

      mixt = mixture_estimation(error , *imixt , estimate , min_inf_bound ,
                                mixt_flag , component_flag , weight_step);
      delete imixt;

      if (mixt) {
        likelihood[min_nb_component] = mixt->Distribution::likelihood_computation(*this);
        nb_parameter[min_nb_component] = mixt->nb_parameter_computation();

        if ((likelihood[min_nb_component] != D_INF) && (penalty_type == BICc)) {
          penalized_likelihood[min_nb_component] = likelihood[min_nb_component] -
                                                   0.5 * mixt->penalty_computation();
        }
      }
    }

    if (likelihood[min_nb_component] != D_INF) {
      if (penalty_type == AICc) {
        if (nb_parameter[min_nb_component] < nb_element - 1) {
          penalized_likelihood[min_nb_component] = likelihood[min_nb_component] -
                                                   (double)(nb_parameter[min_nb_component] * nb_element) /
                                                   (double)(nb_element - nb_parameter[min_nb_component] - 1);
        }
        else {
          penalized_likelihood[min_nb_component] = D_INF;
        }
      }

      else if (penalty_type != BICc) {
        penalized_likelihood[min_nb_component] = likelihood[min_nb_component] -
                                                 nb_parameter[min_nb_component] * penalty;
      }
    }

    else {
      penalized_likelihood[min_nb_component] = D_INF;
    }

    max_likelihood = penalized_likelihood[min_nb_component];

    for (i = min_nb_component + 1;i <= max_nb_component;i++) {
      imixt = new Mixture(i , pcomponent);

      pmixt = mixture_estimation(error , *imixt , estimate , min_inf_bound ,
                                 mixt_flag , component_flag , weight_step);
      delete imixt;

      if (pmixt) {
        likelihood[i] = pmixt->Distribution::likelihood_computation(*this);
        nb_parameter[i] = pmixt->nb_parameter_computation();

        if (penalty_type == AICc) {
          if (nb_parameter[i] < nb_element - 1) {
            penalized_likelihood[i] = likelihood[i] - (double)(nb_parameter[i] * nb_element) /
                                      (double)(nb_element - nb_parameter[i] - 1);
          }
          else {
            penalized_likelihood[i] = D_INF;
          }
        }

        else if (penalty_type == BICc) {
          penalized_likelihood[i] = likelihood[i] - 0.5 * pmixt->penalty_computation();
        }

        else {
          penalized_likelihood[i] = likelihood[i] - nb_parameter[i] * penalty;
        }

        if (penalized_likelihood[i] > max_likelihood) {
          max_likelihood = penalized_likelihood[i];
          delete mixt;
          mixt = pmixt;
        }
        else {
          delete pmixt;
        }
      }

      else {
        likelihood[i] = D_INF;
      }
    }

#   ifdef MESSAGE
    {
      double norm = 0. , weight[MIXTURE_NB_COMPONENT + 1];

      for (i = min_nb_component;i <= max_nb_component;i++) {
        if (likelihood[i] != D_INF) {
          weight[i] = exp(penalized_likelihood[i] - max_likelihood);
          norm += weight[i];
        }
      }

      for (i = min_nb_component;i <= max_nb_component;i++) {
        if (likelihood[i] != D_INF) {
          os << "\n" << i << " " << STAT_label[i == 1 ? STATL_DISTRIBUTION : STATL_DISTRIBUTIONS]
             << "   2 * " << STAT_label[STATL_LIKELIHOOD] << ": " << 2 * likelihood[i] << "   "
             << nb_parameter[i] << " " << STAT_label[STATL_FREE_PARAMETERS]
             << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("
             << STAT_criterion_word[penalty_type] << "): " << 2 * penalized_likelihood[i] << "   "
             << STAT_label[STATL_WEIGHT] << ": " << weight[i] / norm << endl;
        }
      }
    }
#   endif

    for (i = 0;i < max_nb_component;i++) {
      delete pcomponent[i];
    }
  }

  return mixt;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par un melange de lois.
 *
 *  arguments : reference sur un objet Format_error, effectif.
 *
 *--------------------------------------------------------------*/

Mixture_data* Mixture::simulation(Format_error &error , int nb_element) const

{
  register int i , j;
  int value;
  Mixture_data *mixt_histo;


  error.init();

  if ((nb_element < 1) || (nb_element > DIST_NB_ELEMENT)) {
    mixt_histo = NULL;
    error.update(STAT_error[STATR_SAMPLE_SIZE]);
  }

  else {

    // creation d'un objet Mixture_data

    mixt_histo = new Mixture_data(*this);
    mixt_histo->mixture = new Mixture(*this , false);

    for (i = 0;i < nb_element;i++) {

      // ponderation

      j = weight->simulation();
      (mixt_histo->weight->frequency[j])++;

      // composante

      value = component[j]->simulation();
      (mixt_histo->component[j]->frequency[value])++;
      (mixt_histo->frequency[value])++;
    }

    // extraction des caracteristiques des histogrammes

    mixt_histo->nb_value_computation();
    mixt_histo->offset_computation();
    mixt_histo->nb_element = nb_element;
    mixt_histo->max_computation();
    mixt_histo->mean_computation();
    mixt_histo->variance_computation();

    mixt_histo->weight->nb_value_computation();
    mixt_histo->weight->offset_computation();
    mixt_histo->weight->nb_element = nb_element;
    mixt_histo->weight->max_computation();

    for (i = 0;i < mixt_histo->nb_component;i++) {
      mixt_histo->component[i]->nb_value_computation();
      mixt_histo->component[i]->offset_computation();
      mixt_histo->component[i]->nb_element_computation();
      mixt_histo->component[i]->max_computation();
      mixt_histo->component[i]->mean_computation();
      mixt_histo->component[i]->variance_computation();
    }
  }

  return mixt_histo;
}
