/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2002 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): J.-B. Durand and Y. Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: mixture_algorithms.cpp 3256 2007-06-06 12:32:54Z dufourko $
 *
 *       Forum for AMAPmod developers: amldevlp@cirad.fr
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
#include "mv_mixture.h"
#include "stat_label.h"
#include "vectors.h"
#include "markovian.h"

using namespace std;



/*--------------------------------------------------------------*
 *
 *  Calcul de la quantite d'information d'un objet Mv_Mixture_data.
 *
 *--------------------------------------------------------------*/

double Mv_Mixture_data::information_computation() const

{
  register int i, var;
  double information , buff;

  information = weight->information_computation();

  if (information != D_INF) {
      for (i = 0;i < nb_component;i++) {
	if (weight->frequency[i] > 0) {
	  buff = 0.;
	  for (var = 0;var < nb_variable;var++) {
	    buff += component[var][i]->information_computation();
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
  }

  return information;
}

/*--------------------------------------------------------------*
 *
 *  Calcul de la densite conditionnelle d'un objet Mv_Mixture_data pour un melange de lois.
 *
 *  argument : reference sur un objet Mv_Mixture_data.
 *
 *--------------------------------------------------------------*/

void Mv_Mixture::get_output_conditional_distribution(const Mv_Mixture_data &mixt_data, 
						     double** &output_cond,
						     bool log_computation) const
{
  int n, i, var;

  if (output_cond == NULL) {
    output_cond = new double*[mixt_data.nb_vector];
    for (n = 0; n < mixt_data.nb_vector; n++)
      output_cond[n] = NULL;
  }

  for (n = 0; n < mixt_data.nb_vector; n++)
    if (output_cond[n] == NULL)
      output_cond[n] = new double[nb_component];

  for (n = 0; n < mixt_data.nb_vector; n++) {
    for (i = 0; i < nb_component; i++) {
      if (log_computation)
	output_cond[n][i] = 0.;
      else
	output_cond[n][i] = 1.;
      for (var = 0; var < mixt_data.nb_variable; var++) {
	if (pcomponent[var] != NULL) {
	  if (log_computation) {
	    if (pcomponent[var]->observation[i]->mass[mixt_data.int_vector[n][var]] > 0)
	      output_cond[n][i] += log(pcomponent[var]->observation[i]->mass[mixt_data.int_vector[n][var]]);
	    else {
	      output_cond[n][i] = D_INF;
	      break;
	    }
	  }
	  else
	    output_cond[n][i]*= pcomponent[var]->observation[i]->mass[mixt_data.int_vector[n][var]];
	}
	else {
	  if (log_computation) {
	    if (npcomponent[var]->get_observation(i)->mass[mixt_data.int_vector[n][var]] > 0)
	      output_cond[n][i] += log(npcomponent[var]->get_observation(i)->mass[mixt_data.int_vector[n][var]]);
	    else {
	      output_cond[n][i] = D_INF;
	      break;
	    }
	  }
	  else
	    output_cond[n][i]*= npcomponent[var]->get_observation(i)->mass[mixt_data.int_vector[n][var]];
	}
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance d'un objet Mv_Mixture_data pour un melange de lois.
 *
 *  argument : reference sur un objet Mv_Mixture_data.
 *
 *--------------------------------------------------------------*/

double Mv_Mixture::likelihood_computation(const Mv_Mixture_data &mixt_data, 
					  bool log_computation) const

{
  register int n, i, var;
  double likelihood = 0. , buff , **output_cond = NULL;
  double *pmass = NULL;
 
  get_output_conditional_distribution(mixt_data, output_cond, false);

  for (n = 0; n < mixt_data.nb_vector; n++) {
    buff = 0.;
    for (i = 0;i < mixt_data.nb_component;i++)
      buff += weight->mass[i] * output_cond[n][i];

    if (buff >= D_INF)
      likelihood += buff;
    else {
      likelihood = D_INF;
      break;
    }
  }

  if (log_computation) {
    if (likelihood > 0)
      likelihood = log(likelihood);
    else
      likelihood = D_INF;
  }

  for (n = 0; n < mixt_data.nb_vector; n++) {
    delete [] output_cond[n];
    output_cond[n] = NULL;
    }

  delete [] output_cond;
  output_cond = NULL;
  
  return likelihood;
}

/*--------------------------------------------------------------*
 *
 *  Initialisation d'un melange de lois a partir d'un histogramme.
 *
 *  arguments : flag sur le decalage des composantes.
 *
 *--------------------------------------------------------------*/

void Mv_Mixture::init(bool component_flag)

{
  /*  register int i , j = -1;
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
  */
}


/*--------------------------------------------------------------*
 *
 *  Calcul des histogrammes correspondant a chacune des composantes
 *  (estimateur EM d'un melange de lois).
 *
 *  arguments : pointeur sur un objet Mv_Mixture_data,
 *              effectif theorique de l'histogramme.
 *
 *--------------------------------------------------------------*/

void Mv_Mixture::expectation_step(Mv_Mixture_data *mixt_data , int nb_element) const

{
  /* register int i , j , k;
  int component_index , value_index , *pfrequency , *mfrequency;
  double scale , sum , max_frequency , *pweight , *mmass , **rfrequency;


  scale = (double)nb_element / (double)mixt_data->nb_element;

  rfrequency = new double*[nb_component];
  for (i = 0;i < nb_component;i++) {
    rfrequency[i] = new double[mixt_data->nb_value];
  }

  mfrequency = mixt_data->frequency + mixt_data->offset;
  mmass = mass + mixt_data->offset;
  sum = 0.;

  for (i = mixt_data->offset;i < mixt_data->nb_value;i++) {
    if ((*mfrequency > 0) && (*mmass > 0.)) {
      pweight = weight->mass;

      // repartition de l'effectif d'une classe entre les histogrammes
      // correspondant a chacune des composantes

      for (j = 0;j < nb_component;j++) {
        pfrequency = mixt_data->component[j]->frequency + i;

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
        mixt_data->component[j]->frequency[i] = 0;
      }
    }

    mfrequency++;
    mmass++;
  }

  // prise en compte des arrondis

  for (i = 0;i < (int)round(sum);i++) {
    max_frequency = 0.;
    mfrequency = mixt_data->frequency + mixt_data->offset;
    mmass = mass + mixt_data->offset;

    for (j = mixt_data->offset;j < mixt_data->nb_value;j++) {
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
    (mixt_data->component[component_index]->frequency[value_index])++;
  }

  for (i = 0;i < nb_component;i++) {
    delete [] rfrequency[i];
  }
  delete [] rfrequency;

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < nb_component;i++) {
    mixt_data->component[i]->nb_value_computation();
    mixt_data->component[i]->offset_computation();
    mixt_data->component[i]->nb_element_computation();
    mixt_data->component[i]->max_computation();
    mixt_data->component[i]->mean_computation();
    mixt_data->component[i]->variance_computation();

#   ifdef DEBUG
    cout << *mixt_data->component[i];
#   endif

  }

  // mise a jour de l'histogramme des poids

  for (i = 0;i < nb_component;i++) {
    mixt_data->weight->frequency[i] = mixt_data->component[i]->nb_element;
  }
  mixt_data->weight->nb_element = nb_element;
  mixt_data->weight->max_computation();*/
}



/*--------------------------------------------------------------*
 *
 *  Estimation du nombre de composantes et des parametres
 *  d'un melange de lois par l'algorithme EM.
 *
 *  arguments : reference sur un objet Format_error, modele initial,
 *  nombre d'iteration maximal, flag sur l'utilisation forcee ou
 *  non de lois d'observation de meme nature que le modele initial
 *  (parametriques ou non)
 *
 *--------------------------------------------------------------*/

Mv_Mixture* Vectors::mixture_estimation(Format_error &error, const Mv_Mixture &imixture, 
					int nb_iter, bool force_param) const
{
  bool status = true , estimate[MIXTURE_NB_COMPONENT];
  register int i;
  int nb_parameter[MIXTURE_NB_COMPONENT + 1];
  double penalty , max_likelihood , likelihood[MIXTURE_NB_COMPONENT + 1] ,
         penalized_likelihood[MIXTURE_NB_COMPONENT + 1];
  const Parametric *pcomponent[MIXTURE_NB_COMPONENT];
  Parametric_model *dist;
  Mv_Mixture *imixt , *mixt , *pmixt;

  return mixt;
}

/*--------------------------------------------------------------*
 *
 *  Estimation du nombre de composantes et des parametres
 *  d'un melange de lois par l'algorithme EM.
 *
 *  arguments : reference sur un objet Format_error, nombre de composantes,
 *  nombre d'iterations maximal, flag sur l'utilisation forcee ou
 *  non de lois d'observation parametriques (pour chaque variable)
 *
 *--------------------------------------------------------------*/

Mv_Mixture*  Vectors::mixture_estimation(Format_error &error, int nb_component, 
					 int nb_iter, bool *force_param) const
{
  bool status = true , estimate[MIXTURE_NB_COMPONENT];
  register int i;
  int nb_parameter[MIXTURE_NB_COMPONENT + 1];
  double penalty , max_likelihood , likelihood[MIXTURE_NB_COMPONENT + 1] ,
         penalized_likelihood[MIXTURE_NB_COMPONENT + 1];
  const Parametric *pcomponent[MIXTURE_NB_COMPONENT];
  Parametric_model *dist;
  Mv_Mixture *imixt , *mixt , *pmixt;

  return mixt;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par un melange de lois.
 *
 *  arguments : reference sur un objet Format_error, effectif.
 *
 *--------------------------------------------------------------*/

Mv_Mixture_data* Mv_Mixture::simulation(Format_error &error , int nb_element) const

{
  int k , n, var;
  int value;
  int *iidentifier, **iint_vector;
  Mv_Mixture_data *mixt_data = NULL;
  Vectors *vec = NULL;
  Histogram *hweight = NULL;
  Histogram ***hcomponent = NULL;


  error.init();

  if ((nb_element < 1) || (nb_element > DIST_NB_ELEMENT))
    error.update(STAT_error[STATR_SAMPLE_SIZE]);
  else {

    // creation d'un objet Mv_Mixture_data
    iidentifier = new int[nb_element];
    iint_vector = new int*[nb_element];
    hweight = new Histogram(nb_component);
    hcomponent = new Histogram**[nb_var];

    for (var = 0; var < nb_var; var++) {
      hcomponent[var] = new Histogram*[nb_component];
      for (k = 0; k < nb_component; k++)
	hcomponent[var][k] = new Histogram(nb_element);
    }

    for (n = 0; n < nb_element; n++) {

      iint_vector[n] = new int[nb_var];
      iidentifier[n] = n+1;

      // ponderation (simulation de l'etat)

      k = weight->simulation();
      (hweight->frequency[k])++;

      // composante
      for (var = 0; var < nb_var; var++) {
	if (pcomponent[var] != NULL)
	  value = pcomponent[var]->observation[k]->simulation();
	else
	  value = npcomponent[var]->get_observation(k)->simulation();
	(hcomponent[var][k]->frequency[value])++;	
	iint_vector[n][var] = value;
      }
    } // end for (n)

    vec = new Vectors(nb_element, iidentifier, nb_var, iint_vector);

    for (n = 0; n < nb_element; n++) {
      delete [] iint_vector[n];
      iint_vector[n] = NULL;
    }

    // calcul des caracteristiques des histogrammes
    for (var = 0; var < nb_var; var++) {
      for (k = 0; k < nb_component; k++) {
	hcomponent[var][k]->nb_value_computation();
	hcomponent[var][k]->offset_computation();
	hcomponent[var][k]->nb_element_computation();
	hcomponent[var][k]->max_computation();
	hcomponent[var][k]->mean_computation();
	hcomponent[var][k]->variance_computation();
      }
    }

    hweight->nb_value_computation();
    hweight->offset_computation();
    hweight->nb_element_computation();
    hweight->max_computation();
    hweight->mean_computation();
    hweight->variance_computation();
    
    delete [] iint_vector;
    iint_vector = NULL;
    delete [] iidentifier;
    iidentifier = NULL;

    // extraction des caracteristiques des histogrammes
    mixt_data = new Mv_Mixture_data(*vec, nb_component);
    mixt_data->mixture = new Mv_Mixture(*this , false);

    mixt_data->weight = hweight;
    mixt_data->component = hcomponent;
  }
  return mixt_data;
}
