/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2002 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): Y. Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
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
#include <sstream>
#include "stat_tool/stat_tools.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequences.h"
#include "semi_markov.h"
#include "sequence_label.h"

using namespace std;


extern double interval_bisection(Reestimation<double> *distribution_reestim ,
                                 Reestimation<double> *length_bias_reestim);
extern int cumul_method(int nb_value , const double *cumul , double scale = 1.);



/*--------------------------------------------------------------*
 *
 *  Calcul de la loi stationnaire pour une semi-chaine de Markov en equilibre.
 *
 *--------------------------------------------------------------*/

void Semi_markov::initial_probability_computation()

{
  register int i , j , k;
  double sum , *state , *state_out , **state_in;
  Parametric *occupancy;


  state = new double[nb_state];
  state_out = new double[nb_state];

  state_in = new double*[STATIONARY_PROBABILITY_LENGTH];
  for (i = 0;i < STATIONARY_PROBABILITY_LENGTH;i++) {
    state_in[i] = new double[nb_state];
  }

  i = 0;

  do {
    if (i > 0) {
      sum = 0.;
    }

    for (j = 0;j < nb_state;j++) {
      if (i > 0) {
        sum += fabs(state_in[i - 1][j] - state_out[j]);
      }

      switch (state_subtype[j]) {

      // cas etat semi-markovien

      case SEMI_MARKOVIAN : {
        if (i == 0) {
          state[j] = initial[j];
        }
        else {
          state[j] += state_in[i - 1][j] - state_out[j];
        }

        occupancy = nonparametric_process[0]->sojourn_time[j];
        state_out[j] = 0.;

        for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
          if (k < i + 1) {
            state_out[j] += occupancy->mass[k] * state_in[i - k][j];
          }
          else {
            state_out[j] += forward[j]->mass[k] * initial[j];
          }
        }
        break;
      }

      // cas etat markovien

      case MARKOVIAN : {
        if (i == 0) {
          state_out[j] = initial[j];
        }
        else {
          state_out[j] = state_in[i - 1][j];
        }
        break;
      }
      }
    }

    for (j = 0;j < nb_state;j++) {
      state_in[i][j] = 0.;
      for (k = 0;k < nb_state;k++) {
        state_in[i][j] += transition[k][j] * state_out[k];
      }
    }

#   ifdef DEBUG
//    if ((i > 0) && (i % 100 == 0)) {
      cout << i << "  ";
      for (j = 0;j < nb_state;j++) {
        cout << state[j] << " ";
      }
      cout << " | " << sum / nb_state << endl;
//    }
#   endif

    i++;
  }
  while (((i == 1) || (sum / nb_state > STATIONARY_PROBABILITY_THRESHOLD)) &&
         (i < STATIONARY_PROBABILITY_LENGTH));

# ifdef DEBUG
  cout << "\n" << SEQ_label[SEQL_LENGTH] << ": "  << i << endl;
# endif

  for (j = 0;j < nb_state;j++) {
    switch (state_subtype[j]) {

    // cas etat semi-markovien

    case SEMI_MARKOVIAN :
      initial[j] = state_in[i - 1][j] - state_out[j] + state[j];
//      initial[j] = state[j];
      break;

    // cas etat markovien

    case MARKOVIAN :
      initial[j] = state_in[i - 1][j];
//      initial[j] = state_out[j];
      break;
    }
  }

  // renormalisation pour tenir compte des seuils appliques sur
  // les fonctions de repartition des lois d'occupation des etats

  sum = 0.;
  for (i = 0;i < nb_state;i++) {
    sum += initial[i];
  }

  if (sum < 1.) {
    for (i = 0;i < nb_state;i++) {
      initial[i] /= sum;
    }
  }

  delete [] state;
  delete [] state_out;

  for (i = 0;i < STATIONARY_PROBABILITY_LENGTH;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance de sequences pour une semi-chaine de Markov.
 *
 *  arguments : reference sur un objet Markovian_sequences, indice de la sequence.
 *
 *--------------------------------------------------------------*/

double Semi_markov::likelihood_computation(const Markovian_sequences &seq , int index) const

{
  register int i , j , k , m;
  int nb_value , occupancy , *pstate , **poutput;
  double likelihood = 0. , proba;


  // verification de la compatibilite entre le modele et les donnees

  if (nb_output_process + 1 == seq.nb_variable) {
    for (i = 0;i <= nb_output_process;i++) {
      if (nonparametric_process[i]) {
        nb_value = nonparametric_process[i]->nb_value;
      }
      else {
        nb_value = parametric_process[i]->nb_value;
      }

      if ((seq.marginal[i]) && (nb_value < seq.marginal[i]->nb_value)) {
        likelihood = D_INF;
        break;
      }
    }
  }

  else {
    likelihood = D_INF;
  }

  if (likelihood != D_INF) {
    if (nb_output_process > 0) {
      poutput = new int*[nb_output_process];
    }

    for (i = 0;i < seq.nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
        pstate = seq.int_sequence[i][0];

        proba = initial[*pstate];
        if (proba > 0.) {
          likelihood += log(proba);
        }
        else {
          likelihood = D_INF;
          break;
        }

        if (nb_output_process > 0) {
          for (j = 0;j < nb_output_process;j++) {
            poutput[j] = seq.int_sequence[i][j + 1];
          }
        }

        j = 0;
        do {
          if (j > 0) {
            pstate++;

            proba = transition[*(pstate - 1)][*pstate];
            if (proba > 0.) {
              likelihood += log(proba);
            }
            else {
              likelihood = D_INF;
              break;
            }
          }

          if (transition[*pstate][*pstate] < 1.) {
            j++;
            occupancy = 1;

            if (state_subtype[*pstate] == SEMI_MARKOVIAN) {
              while ((j < seq.length[i]) && (*(pstate + 1) == *pstate)) {
                j++;
                occupancy++;
                pstate++;
              }

              proba = 0.;
              if ((type == 'e') && (j == occupancy)) {
                if (occupancy < forward[*pstate]->nb_value) {
                  if (j < seq.length[i]) {
                    proba = forward[*pstate]->mass[occupancy];
                  }
                  else {
                    proba = (1. - forward[*pstate]->cumul[occupancy - 1]);
                  }
                }
              }

              else {
                if (occupancy < nonparametric_process[0]->sojourn_time[*pstate]->nb_value) {
                  if (j < seq.length[i]) {
                    proba = nonparametric_process[0]->sojourn_time[*pstate]->mass[occupancy];
                  }
                  else {
                    proba = (1. - nonparametric_process[0]->sojourn_time[*pstate]->cumul[occupancy - 1]);
                  }
                }
              }

              if (proba > 0.) {
                likelihood += log(proba);
              }
              else {
                likelihood = D_INF;
                break;
              }
            }
          }

          else {
            occupancy = seq.length[i] - j;
            j += occupancy;
          }

          if (nb_output_process > 0) {
            for (k = 0;k < occupancy;k++) {
              for (m = 0;m < nb_output_process;m++) {
                if (nonparametric_process[m + 1]) {
                  proba = nonparametric_process[m + 1]->observation[*pstate]->mass[*poutput[m]++];
                }
                else {
                  proba = parametric_process[m + 1]->observation[*pstate]->mass[*poutput[m]++];
                }

                if (proba > 0.) {
                  likelihood += log(proba);
                }
                else {
                  likelihood = D_INF;
                  break;
                }
              }

              if (likelihood == D_INF) {
                break;
              }
            }

            if (likelihood == D_INF) {
              break;
            }
          }
        }
        while (j < seq.length[i]);

        if (likelihood == D_INF) {
          break;
        }
      }
    }

    if (nb_output_process > 0) {
      delete [] poutput;
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance de sequences pour une semi-chaine de Markov.
 *
 *  argument : reference sur un objet Semi_markov_data.
 *
 *--------------------------------------------------------------*/

double Semi_markov::likelihood_computation(const Semi_markov_data &seq) const

{
  register int i , j;
  int nb_value;
  double buff , likelihood;
  Histogram **initial_run , **final_run , **single_run;


  // verification de la compatibilite entre le modele et les donnees

  if (nb_output_process + 1 == seq.nb_variable) {
    for (i = 0;i <= nb_output_process;i++) {
      if (nonparametric_process[i]) {
        nb_value = nonparametric_process[i]->nb_value;
      }
      else {
        nb_value = parametric_process[i]->nb_value;
      }

      if (nb_value < seq.marginal[i]->nb_value) {
        likelihood = D_INF;
        break;
      }
    }
  }

  else {
    likelihood = D_INF;
  }

  if (likelihood != D_INF) {
    likelihood = Chain::likelihood_computation(*(seq.chain_data));

    if (likelihood != D_INF) {
      if (type == 'e') {

        // creation des histogrammes des temps de sejour censures

        initial_run = new Histogram*[seq.marginal[0]->nb_value];
        for (i = 0;i < seq.marginal[0]->nb_value;i++) {
          initial_run[i] = new Histogram(seq.max_length);
        }

        final_run = new Histogram*[seq.marginal[0]->nb_value];
        for (i = 0;i < seq.marginal[0]->nb_value;i++) {
          final_run[i] = new Histogram(seq.max_length);
        }

        single_run = new Histogram*[seq.marginal[0]->nb_value];
        for (i = 0;i < seq.marginal[0]->nb_value;i++) {
          single_run[i] = new Histogram(seq.max_length + 1);
        }

        // mise a jour des histogrammes des temps de sejour censures

        seq.censored_sojourn_time_histogram_computation(initial_run , final_run , single_run);
      }

      for (i = 0;i < nb_state;i++) {
        if (state_subtype[i] == SEMI_MARKOVIAN) {
          buff = nonparametric_process[0]->sojourn_time[i]->likelihood_computation(*(seq.characteristics[0]->sojourn_time[i]));

          if (buff != D_INF) {
            likelihood += buff;
          }
          else {
            likelihood = D_INF;
            break;
          }

          switch (type) {

          case 'o' : {
            buff = nonparametric_process[0]->sojourn_time[i]->survivor_likelihood_computation(*(seq.characteristics[0]->final_run[i]));

            if (buff != D_INF) {
              likelihood += buff;
            }
            else {
              likelihood = D_INF;
            }
            break;
          }

          case 'e' : {
            buff = nonparametric_process[0]->sojourn_time[i]->survivor_likelihood_computation(*(final_run[i]));

            if (buff != D_INF) {
              likelihood += buff;
              buff = forward[i]->likelihood_computation(*(initial_run[i]));

              if (buff != D_INF) {
                likelihood += buff;
                buff = forward[i]->survivor_likelihood_computation(*(single_run[i]));

                if (buff != D_INF) {
                  likelihood += buff;
                }
                else {
                  likelihood = D_INF;
                }
              }

              else {
                likelihood = D_INF;
              }
            }

            else {
              likelihood = D_INF;
            }
            break;
          }
          }

          if (likelihood == D_INF) {
            break;
          }
        }
      }

      if (type == 'e') {
        for (i = 0;i < seq.marginal[0]->nb_value;i++) {
          delete initial_run[i];
        }
        delete [] initial_run;

        for (i = 0;i < seq.marginal[0]->nb_value;i++) {
          delete final_run[i];
        }
        delete [] final_run;

        for (i = 0;i < seq.marginal[0]->nb_value;i++) {
          delete single_run[i];
        }
        delete [] single_run;
      }
    }

    if (likelihood != D_INF) {
      for (i = 1;i <= nb_output_process;i++) {
        if (nonparametric_process[i]) {
          for (j = 0;j < nb_state;j++) {
            buff = nonparametric_process[i]->observation[j]->likelihood_computation(*(seq.observation[i][j]));

            if (buff != D_INF) {
              likelihood += buff;
            }
            else {
              likelihood = D_INF;
              break;
            }
          }
        }

        else {
          for (j = 0;j < nb_state;j++) {
            buff = parametric_process[i]->observation[j]->likelihood_computation(*(seq.observation[i][j]));

            if (buff != D_INF) {
              likelihood += buff;
            }
            else {
              likelihood = D_INF;
              break;
            }
          }
        }

        if (likelihood == D_INF) {
          break;
        }
      }
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Comptage des etats initiaux et des transitions.
 *
 *  arguments : reference sur un objet Chain_data,
 *              flags sur les probabilites de rester dans un etat.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::transition_count_computation(const Chain_data &chain_data ,
                                                       const Semi_markov *smarkov) const

{
  register int i , j;
  int *pstate;


  for (i = 0;i < chain_data.nb_state;i++) {
    chain_data.initial[i] = 0;
  }

  for (i = 0;i < chain_data.nb_row;i++) {
    for (j = 0;j < chain_data.nb_state;j++) {
      chain_data.transition[i][j] = 0;
    }
  }

  // extraction des etats initiaux et des transitions

  for (i = 0;i < nb_sequence;i++) {
    pstate = int_sequence[i][0];
    (chain_data.initial[*pstate])++;

    for (j = 1;j < length[i];j++) {
      pstate++;
      (chain_data.transition[*(pstate - 1)][*pstate])++;
    }
  }

  if (smarkov) {
    for (i = 0;i < chain_data.nb_state;i++) {
      if (smarkov->state_subtype[i] == SEMI_MARKOVIAN) {
        chain_data.transition[i][i] = 0;
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction des comptages des etats initiaux et des transitions.
 *
 *  argument : flags sur les probabilites de rester dans un etat.
 *
 *--------------------------------------------------------------*/

void Semi_markov_data::build_transition_count(const Semi_markov *smarkov)

{
  chain_data = new Chain_data('o' , marginal[0]->nb_value , marginal[0]->nb_value);
  transition_count_computation(*chain_data , smarkov);
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance des temps de sejour dans un etat pour
 *  une semi-chaine de Markov ordinaire.
 *
 *  arguments : histogrammes des temps de sejour complets et censures a droite.
 *
 *--------------------------------------------------------------*/

double Parametric::state_occupancy_likelihood_computation(const Histogram &sojourn_time ,
                                                          const Histogram &final_run) const

{
  double likelihood , buff;


  likelihood = likelihood_computation(sojourn_time);

  if (likelihood != D_INF) {
    buff = survivor_likelihood_computation(final_run);

    if (buff != D_INF) {
      likelihood += buff;
    }
    else {
      likelihood = D_INF;
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance des temps de sejour dans un etat pour
 *  une semi-chaine de Markov en equilibre.
 *
 *  arguments : loi de l'intervalle de temps residuel, histogrammes des temps
 *              de sejour complets, censures a droite, a gauche et
 *              de la longueur des sequences dans le cas d'un seul etat visite.
 *
 *--------------------------------------------------------------*/

double Parametric::state_occupancy_likelihood_computation(const Forward &forward , const Histogram &sojourn_time ,
                                                          const Histogram &final_run , const Histogram &initial_run ,
                                                          const Histogram &single_run) const

{
  double likelihood , buff;


  likelihood = likelihood_computation(sojourn_time);

  if (likelihood != D_INF) {
    buff = survivor_likelihood_computation(final_run);

    if (buff != D_INF) {
      likelihood += buff;
      buff = forward.likelihood_computation(initial_run);

      if (buff != D_INF) {
        likelihood += buff;
        buff = forward.survivor_likelihood_computation(single_run);

        if (buff != D_INF) {
          likelihood += buff;
        }
        else {
          likelihood = D_INF;
        }
      }

      else {
        likelihood = D_INF;
      }
    }

    else {
      likelihood = D_INF;
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des quantites de reestimation correspondant a la loi d'occupation d'un etat
 *  (estimateur EM d'une semi-chaine de Markov ordinaire).
 *
 *  arguments : histogrammes des temps de sejour complets et censures a droite,
 *              pointeurs sur les quantites de reestimation de la loi d'occupation de l'etat.
 *
 *--------------------------------------------------------------*/

void Parametric::expectation_step(const Histogram &sojourn_time , const Histogram &final_run ,
                                  Reestimation<double> *occupancy_reestim , int iter) const

{
  register int i;
  int *pfrequency;
  double sum , *ofrequency , *pmass , *pcumul;


  // initialisations

  ofrequency = occupancy_reestim->frequency;
  for (i = 0;i < occupancy_reestim->alloc_nb_value;i++) {
    *ofrequency++ = 0.;
  }

  // calcul des quantites de reestimation de la loi d'occupation de l'etat

  ofrequency = occupancy_reestim->frequency + sojourn_time.offset;
  pfrequency = sojourn_time.frequency + sojourn_time.offset;
  for (i = sojourn_time.offset;i < sojourn_time.nb_value;i++) {
    *ofrequency++ += *pfrequency++;
  }

  pfrequency = final_run.frequency + final_run.offset;
  pcumul = cumul + final_run.offset - 1;
  ofrequency = occupancy_reestim->frequency + final_run.offset;
  pmass = mass + final_run.offset;
  sum = 0.;

  for (i = final_run.offset;i < final_run.nb_value;i++) {
    sum += *pfrequency++ / (1. - *pcumul++);
    *ofrequency++ += *pmass++ * sum;
  }
  for (i = final_run.nb_value;i < nb_value;i++) {
    *ofrequency++ += *pmass++ * sum;
  }

  occupancy_reestim->nb_value_computation();
  occupancy_reestim->offset_computation();
  occupancy_reestim->nb_element_computation();

# ifdef DEBUG
  if ((iter < 10) || ((iter < 100) && (iter % 10 == 0)) ||
      ((iter < 1000) && (iter % 100 == 0)) || (iter % 1000 == 0)) {
    occupancy_reestim->max_computation();
    occupancy_reestim->mean_computation();
    occupancy_reestim->variance_computation();

    cout << "\nquantites de reestimation loi d'occupation de l'etat :" << *occupancy_reestim << endl;
  }
# endif

}


/*--------------------------------------------------------------*
 *
 *  Calcul des quantites de reestimation correspondant a la loi d'occupation d'un etat
 *  (estimateur EM d'une semi-chaine de Markov en equilibre).
 *
 *  arguments : histogrammes des temps de sejour complets, censures a droite, a gauche et
 *              de la longueur des sequences dans le cas d'un seul etat visite,
 *              pointeurs sur les quantites de reestimation de la loi d'occupation de l'etat et
 *              de la loi biaisee par la longueur, combinaison ou non des quantites de reestimation,
 *              methode de calcul de la moyenne de la loi d'occupation de l'etat.
 *
 *--------------------------------------------------------------*/

void Parametric::expectation_step(const Histogram &sojourn_time , const Histogram &final_run ,
                                  const Histogram &initial_run , const Histogram &single_run ,
                                  Reestimation<double> *occupancy_reestim ,
                                  Reestimation<double> *length_bias_reestim , int iter ,
                                  bool combination , int mean_computation) const

{
  register int i , j;
  int reestim_offset , reestim_nb_value , *pfrequency;
  double sum , occupancy_mean , *ofrequency , *lfrequency , *pmass , *pcumul , *norm;


  // initialisations

  ofrequency = occupancy_reestim->frequency;
  lfrequency = length_bias_reestim->frequency;
  for (i = 0;i < occupancy_reestim->alloc_nb_value;i++) {
    *ofrequency++ = 0.;
    *lfrequency++ = 0.;
  }

  // calcul des quantites de reestimation de la loi d'occupation de l'etat

  ofrequency = occupancy_reestim->frequency + sojourn_time.offset;
  pfrequency = sojourn_time.frequency + sojourn_time.offset;
  for (i = sojourn_time.offset;i < sojourn_time.nb_value;i++) {
    *ofrequency++ += *pfrequency++;
  }

  if (final_run.nb_element > 0) {
    pfrequency = final_run.frequency + final_run.offset;
    pcumul = cumul + final_run.offset - 1;
    ofrequency = occupancy_reestim->frequency + final_run.offset;
    pmass = mass + final_run.offset;
    sum = 0.;

    for (i = final_run.offset;i < final_run.nb_value;i++) {
      sum += *pfrequency++ / (1. - *pcumul++);
      *ofrequency++ += *pmass++ * sum;
    }
    for (i = final_run.nb_value;i < nb_value;i++) {
      *ofrequency++ += *pmass++ * sum;
    }
  }

  // calcul des quantites de reestimation de la loi biaisee par la longueur

  if (initial_run.nb_element > 0) {
    pfrequency = initial_run.frequency + initial_run.offset;
    pcumul = cumul + initial_run.offset - 1;
    lfrequency = length_bias_reestim->frequency + initial_run.offset;
    pmass = mass + initial_run.offset;
    sum = 0.;

    for (i = initial_run.offset;i < initial_run.nb_value;i++) {
      sum += *pfrequency++ / (1. - *pcumul++);
      *lfrequency++ += *pmass++ * sum;
    }
    for (i = initial_run.nb_value;i < nb_value;i++) {
      *lfrequency++ += *pmass++ * sum;
    }
  }

  if (single_run.nb_element > 0) {
    norm = new double[single_run.nb_value];

    for (i = single_run.offset;i < single_run.nb_value;i++) {
      if (single_run.frequency[i] > 0) {
        pmass = mass + i;
        norm[i] = 0.;
        for (j = i;j < nb_value;j++) {
          norm[i] += (j + 1 - i) * *pmass++;
        }
      }
    }

    lfrequency = length_bias_reestim->frequency + single_run.offset;
    pmass = mass + single_run.offset;
    for (i = single_run.offset;i < nb_value;i++) {
      pfrequency = single_run.frequency + single_run.offset;
      sum = 0.;
      for (j = single_run.offset;j <= MIN(i , single_run.nb_value - 1);j++) {
        if ((*pfrequency > 0) && (norm[j] > 0.)) {
          sum += *pfrequency * (i + 1 - j) / norm[j];
               }
        pfrequency++;
      }

      *lfrequency++ += *pmass++ * sum;
    }

    delete [] norm;
  }

  reestim_offset = 1;
  reestim_nb_value = occupancy_reestim->alloc_nb_value;

  ofrequency = occupancy_reestim->frequency + occupancy_reestim->alloc_nb_value;
  lfrequency = length_bias_reestim->frequency + occupancy_reestim->alloc_nb_value;
  while ((*--ofrequency == 0) && (*--lfrequency == 0) && (reestim_nb_value > 2)) {
    reestim_nb_value--;
  }
  occupancy_reestim->nb_value = reestim_nb_value;
  length_bias_reestim->nb_value = reestim_nb_value;

  ofrequency = occupancy_reestim->frequency + reestim_offset;
  lfrequency = length_bias_reestim->frequency + reestim_offset;
  while ((*ofrequency++ == 0) && (*lfrequency++ == 0) && (reestim_offset < reestim_nb_value - 1)) {
    reestim_offset++;
  }
  occupancy_reestim->offset = reestim_offset;
  length_bias_reestim->offset = reestim_offset;

  occupancy_reestim->nb_element_computation();
  length_bias_reestim->nb_element_computation();

# ifdef DEBUG
  if ((iter < 10) || ((iter < 100) && (iter % 10 == 0)) ||
      ((iter < 1000) && (iter % 100 == 0)) || (iter % 1000 == 0)) {
    occupancy_reestim->max_computation();
    occupancy_reestim->mean_computation();
    occupancy_reestim->variance_computation();

    length_bias_reestim->max_computation();
    length_bias_reestim->mean_computation();
    length_bias_reestim->variance_computation();

    cout << "\nquantites de reestimation loi d'occupation de l'etat :" << *occupancy_reestim << endl;
    cout << "\nquantites de reestimation loi biaisee par la longueur :" << *length_bias_reestim << endl;
  }
# endif

  if (combination) {
    switch (mean_computation) {
    case COMPUTED :
      occupancy_mean = interval_bisection(occupancy_reestim , length_bias_reestim);
      break;
    case ONE_STEP_LATE :
      occupancy_mean = mean;
      break;
    }

#   ifdef DEBUG
    cout << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_MEAN] << ": " << occupancy_mean << endl;
#   endif

    occupancy_reestim->equilibrium_process_combination(length_bias_reestim , occupancy_mean);

#   ifdef DEBUG
    if ((iter < 10) || ((iter < 100) && (iter % 10 == 0)) ||
        ((iter < 1000) && (iter % 100 == 0)) || (iter % 1000 == 0)) {
      cout << "\nquantites de reestimation loi d'occupation de l'etat :" << *occupancy_reestim << endl;
    }
#   endif

  }
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une semi-chaine de Markov
 *  a partir d'un echantillon de sequences.
 *
 *  arguments : reference sur un objet Format_error, stream,
 *              type de processus ('o' : ordinaire, 'e' : en equilibre),
 *              type d'estimateur pour la reestimation des lois d'occupation des etats,
 *              flag sur le calcul des lois de comptage, nombre d'iterations,
 *              methode de calcul de la moyenne des lois d'occupation des etats
 *              (semi-chaine de Markov en equilibre).
 *
 *--------------------------------------------------------------*/

Semi_markov* Markovian_sequences::semi_markov_estimation(Format_error &error , ostream &os ,
                                                         char model_type , int estimator , bool counting_flag ,
                                                         int nb_iter , int mean_computation) const

{
  bool status = true;
  register int i , j;
  int nb_likelihood_decrease , *occupancy_survivor , *censored_occupancy_survivor , nb_value[1];
  double likelihood , previous_likelihood , hlikelihood , occupancy_mean;
  Parametric *occupancy;
  Forward *forward;
  Reestimation<double> *occupancy_reestim , *length_bias_reestim;
  Semi_markov *smarkov;
  Semi_markov_data *seq;
  Histogram *complete_run , *censored_run , *pfinal_run , *hreestim , **initial_run ,
            **final_run , **single_run;
  const Histogram *prun[3];


  smarkov = 0;
  error.init();

  if ((type[0] != INT_VALUE) && (type[0] != STATE)) {
    status = false;
    ostringstream correction_message;
    correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
    error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
  }

  else {
    if ((marginal[0]->nb_value < 2) || (marginal[0]->nb_value > NB_STATE)) {
      status = false;
      error.update(SEQ_error[SEQR_NB_STATE]);
    }

    else if (!characteristics[0]) {
      for (i = 0;i < marginal[0]->nb_value;i++) {
        if (marginal[0]->frequency[i] == 0) {
          status = false;
          ostringstream error_message;
          error_message << SEQ_error[SEQR_MISSING_STATE] << " " << i;
          error.update((error_message.str()).c_str());
        }
      }
    }
  }

  if (nb_variable > 1) {
    if (nb_variable > 2) {
      status = false;
      error.correction_update(STAT_error[STATR_NB_VARIABLE] , "1 or 2");
    }

    if ((type[1] != INT_VALUE) && (type[1] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << 2 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }

    else {
      if (test_hidden(1)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << 2 << ": "
                      << SEQ_error[SEQR_OVERLAP];
        error.update((error_message.str()).c_str());
      }

      if (marginal[1]->nb_value > NB_STATE) {
        status = false;
        error.update(SEQ_error[SEQR_NB_OUTPUT]);
      }

/*      if (!characteristics[1]) {
        for (i = 0;i < marginal[1]->nb_value;i++) {
          if (marginal[1]->frequency[i] == 0) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << 2 << ": "
                          << STAT_error[STATR_MISSING_VALUE] << " " << i;
            error.update((error_message.str()).c_str());
          }
        }
      } */
    }
  }

  if (status) {
    if (max_length > COUNTING_MAX_LENGTH) {
      counting_flag = false;
    }

    if (model_type == 'e') {

      // creation des histogrammes des temps de sejour censures

      initial_run = new Histogram*[marginal[0]->nb_value];
      for (i = 0;i < marginal[0]->nb_value;i++) {
        initial_run[i] = new Histogram(max_length);
      }

      final_run = new Histogram*[marginal[0]->nb_value];
      for (i = 0;i < marginal[0]->nb_value;i++) {
        final_run[i] = new Histogram(max_length);
      }

      single_run = new Histogram*[marginal[0]->nb_value];
      for (i = 0;i < marginal[0]->nb_value;i++) {
        single_run[i] = new Histogram(max_length + 1);
      }

      // mise a jour des histogrammes des temps de sejour censures

      censored_sojourn_time_histogram_computation(initial_run , final_run , single_run);
    }

    if (nb_variable == 2) {
      nb_value[0] = marginal[1]->nb_value;
    }

    smarkov = new Semi_markov(model_type , marginal[0]->nb_value , nb_variable - 1 , nb_value);
    smarkov->semi_markov_data = new Semi_markov_data(*this , 'c' , (model_type == 'e' ? true : false));

    seq = smarkov->semi_markov_data;
    seq->state_variable_init();
    seq->build_transition_count();

    for (i = 0;i < smarkov->nb_state;i++) {
      if ((seq->characteristics[0]->sojourn_time[i]->nb_element > 0) ||
          ((model_type == 'e') && (initial_run[i]->nb_element > 0))) {
        seq->chain_data->transition[i][i] = 0;
      }
    }

    // estimation des parametres de la chaine de Markov sous-jacente

    seq->chain_data->estimation(*smarkov);
    smarkov->component_computation();

    if ((model_type == 'e') && (smarkov->nb_component > 1)) {
      delete smarkov;
      smarkov = 0;
      error.correction_update(STAT_parsing[STATP_CHAIN_STRUCTURE] , STAT_parsing[STATP_IRREDUCIBLE]);
    }

    else {

      // estimation des lois d'occupation des etats

      if (estimator != PARTIAL_LIKELIHOOD) {
        occupancy_survivor = new int[max_length];
        censored_occupancy_survivor = new int[max_length + 1];
        occupancy_reestim = new Reestimation<double>(max_length + 1);
      }

      smarkov->state_subtype = new int[smarkov->nb_state];
      smarkov->nonparametric_process[0]->absorption = new double[smarkov->nb_state];
      smarkov->nonparametric_process[0]->sojourn_time = new Parametric*[smarkov->nb_state];
      smarkov->forward = new Forward*[smarkov->nb_state];
      for (i = 0;i < smarkov->nb_state;i++) {
        smarkov->forward[i] = 0;
        smarkov->nonparametric_process[0]->sojourn_time[i] = 0;
      }

      for (i = 0;i < smarkov->nb_state;i++) {
        if ((seq->characteristics[0]->sojourn_time[i]->nb_element == 0) &&
            ((model_type == 'o') || (initial_run[i]->nb_element == 0))) {
          smarkov->state_subtype[i] = MARKOVIAN;
          smarkov->nonparametric_process[0]->absorption[i] = 1.;
        }

        else {
          smarkov->nonparametric_process[0]->absorption[i] = 0.;

          if ((model_type == 'e') && (seq->characteristics[0]->sojourn_time[i]->nb_element == 0)) {
            occupancy = 0;
          }

          else if ((estimator != PARTIAL_LIKELIHOOD) && (model_type == 'e') && 
                   ((initial_run[i]->nb_element > 0) || (single_run[i]->nb_element > 0))) {

            //  initialisation de la loi d'occupation de l'etat

            prun[0] = seq->characteristics[0]->sojourn_time[i];
            prun[1] = seq->characteristics[0]->sojourn_time[i];
            complete_run = new Histogram(2 , prun);

//            prun[0] = seq->initial_run[0][i];
//            prun[1] = final_run[i];

//            prun[0] = initial_run[i];
//            prun[1] = seq->characteristics[0]->final_run[i];
//            censored_run = new Histogram(2 , prun);

            prun[0] = initial_run[i];
            prun[1] = final_run[i];
            prun[2] = single_run[i];
            censored_run = new Histogram(3 , prun);

            complete_run->state_occupancy_estimation(censored_run , occupancy_reestim ,
                                                     occupancy_survivor ,
                                                     censored_occupancy_survivor);
            delete complete_run;
            delete censored_run;

            if ((estimator == KAPLAN_MEIER) && (single_run[i]->nb_element == 0)) {
              occupancy = occupancy_reestim->type_parametric_estimation(1 , true , OCCUPANCY_THRESHOLD);
            }

            else {
              occupancy = new Parametric(NEGATIVE_BINOMIAL , 1 , I_DEFAULT , 1. ,
                                         (occupancy_reestim->mean > 1. ? 1. / occupancy_reestim->mean : 0.99) ,
                                         OCCUPANCY_THRESHOLD);
              occupancy->init(NONPARAMETRIC , I_DEFAULT , I_DEFAULT , D_DEFAULT , D_DEFAULT);
              forward = new Forward(*occupancy , occupancy->alloc_nb_value);

              delete occupancy_reestim;
              occupancy_reestim = new Reestimation<double>(MAX(occupancy->alloc_nb_value , max_length + 1));
              length_bias_reestim = new Reestimation<double>(MAX(occupancy->alloc_nb_value , max_length + 1));

              likelihood = D_INF;
              j = 0;

              do {
                j++;

                // calcul des quantites de reestimation de la loi d'occupation de l'etat

                occupancy->expectation_step(*(seq->characteristics[0]->sojourn_time[i]) ,
                                            *(final_run[i]) , *(initial_run[i]) ,
                                            *(single_run[i]) , occupancy_reestim ,
                                            length_bias_reestim , j);

                switch (mean_computation) {
                case COMPUTED :
                  occupancy_mean = interval_bisection(occupancy_reestim , length_bias_reestim);
                  break;
                case ONE_STEP_LATE :
                  occupancy_mean = occupancy->mean;
                  break;
                }

                occupancy_reestim->equilibrium_process_estimation(length_bias_reestim , occupancy ,
                                                                  occupancy_mean);
                forward->computation(*occupancy);

                previous_likelihood = likelihood;
                likelihood = occupancy->state_occupancy_likelihood_computation(*forward , *(seq->characteristics[0]->sojourn_time[i]) ,
                                                                               *(final_run[i]) , *(initial_run[i]) ,
                                                                               *(single_run[i]));

#               ifdef MESSAGE
                if ((j < 10) || ((j < 100) && (j % 10 == 0)) || ((j < 1000) && (j % 100 == 0)) || (j % 1000 == 0)) {
                  os << STAT_label[STATL_ITERATION] << " " << j << "   "
                     << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
                     << STAT_label[STATL_SMOOTHNESS] << ": " << occupancy->second_difference_norm_computation() << endl;
                }
#               endif

              }
              while ((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (j < OCCUPANCY_NB_ITER) &&
                       ((likelihood - previous_likelihood) / -likelihood > OCCUPANCY_LIKELIHOOD_DIFF)) ||
                      ((nb_iter != I_DEFAULT) && (j < nb_iter))));

              if (likelihood != D_INF) {

#               ifdef MESSAGE
                os << "\n" << STAT_label[STATL_STATE] << " " << i << " - "
                   << j << " " << STAT_label[STATL_ITERATIONS] << "   "
                   << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
                   << STAT_label[STATL_SMOOTHNESS] << ": " << occupancy->second_difference_norm_computation() << "\n" << endl;
#               endif

                hreestim = new Histogram(MAX(occupancy->alloc_nb_value , max_length + 1));

                likelihood = D_INF;
                nb_likelihood_decrease = 0;

                j = 0;
                do {
                  j++;

                  // calcul des quantites de reestimation de la loi d'occupation de l'etat

                  occupancy->expectation_step(*(seq->characteristics[0]->sojourn_time[i]) ,
                                              *(final_run[i]) , *(initial_run[i]) ,
                                              *(single_run[i]) , occupancy_reestim ,
                                              length_bias_reestim , j , true , mean_computation);

                  hreestim->update(occupancy_reestim , (int)(occupancy_reestim->nb_element *
                                   MAX(sqrt(occupancy_reestim->variance) , 1.) * OCCUPANCY_COEFF));
                  hlikelihood = hreestim->Reestimation<int>::type_parametric_estimation(occupancy , 1 , true ,
                                                                                        OCCUPANCY_THRESHOLD);

                  if (hlikelihood == D_INF) {
                    likelihood = D_INF;
                  }

                  else {
                    occupancy->computation(hreestim->nb_value , OCCUPANCY_THRESHOLD);
                    forward->copy(*occupancy);
                    forward->computation(*occupancy);

                    previous_likelihood = likelihood;
                    likelihood = occupancy->state_occupancy_likelihood_computation(*forward , *(seq->characteristics[0]->sojourn_time[i]) ,
                                                                                   *(final_run[i]) , *(initial_run[i]) ,
                                                                                   *(single_run[i]));

                    if (likelihood < previous_likelihood) {
                      nb_likelihood_decrease++;
                    }
                    else {
                      nb_likelihood_decrease = 0;
                    }

#                   ifdef MESSAGE
                    if ((j < 10) || ((j < 100) && (j % 10 == 0)) || ((j < 1000) && (j % 100 == 0)) || (j % 1000 == 0)) {
                      os << STAT_label[STATL_ITERATION] << " " << j << "   "
                         << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
                         << STAT_label[STATL_SMOOTHNESS] << ": " << occupancy->second_difference_norm_computation() << endl;
                    }
#                   endif

                  }
                }
                while ((likelihood != D_INF) && (j < OCCUPANCY_NB_ITER) &&
                       (((likelihood - previous_likelihood) / -likelihood > OCCUPANCY_LIKELIHOOD_DIFF) ||
                        (hlikelihood == D_INF) || (nb_likelihood_decrease == 1)));

                delete hreestim;

                if (likelihood != D_INF) {

#                 ifdef MESSAGE
                  os << "\n" << STAT_label[STATL_STATE] << " " << i << " - "
                     << j << " " << STAT_label[STATL_ITERATIONS] << "   "
                     << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
                     << STAT_label[STATL_SMOOTHNESS] << ": " << occupancy->second_difference_norm_computation() << "\n" << endl;
#                 endif

                }

                else {
                  delete occupancy;
                  occupancy = 0;
                }
              }

              else {
                delete occupancy;
                occupancy = 0;
              }

              delete forward;
              delete length_bias_reestim;
            }
          }

          else if ((estimator != PARTIAL_LIKELIHOOD) && (((model_type == 'o') && (seq->characteristics[0]->final_run[i]->nb_element > 0) &&
                     (seq->characteristics[0]->final_run[i]->nb_value > seq->characteristics[0]->sojourn_time[i]->nb_value)) || ((model_type == 'e') &&
                     (final_run[i]->nb_element > 0) && (final_run[i]->nb_value > seq->characteristics[0]->sojourn_time[i]->nb_value)))) {
            switch (model_type) {
            case 'o' :
              pfinal_run = seq->characteristics[0]->final_run[i];
              break;
            case 'e' :
              pfinal_run = final_run[i];
              break;
            }

            //  initialisation de la loi d'occupation de l'etat

            seq->characteristics[0]->sojourn_time[i]->state_occupancy_estimation(pfinal_run , occupancy_reestim ,
                                                                                 occupancy_survivor ,
                                                                                 censored_occupancy_survivor);
            occupancy = new Parametric(NEGATIVE_BINOMIAL , 1 , I_DEFAULT , 1. ,
                                       1. / occupancy_reestim->mean , OCCUPANCY_THRESHOLD);
            occupancy->init(NONPARAMETRIC , I_DEFAULT , I_DEFAULT , D_DEFAULT , D_DEFAULT);

            delete occupancy_reestim;
            occupancy_reestim = new Reestimation<double>(MAX(occupancy->alloc_nb_value , max_length + 1));

            likelihood = D_INF;
            j = 0;

            do {
              j++;

              // calcul des quantites de reestimation de la loi d'occupation de l'etat

              occupancy->expectation_step(*(seq->characteristics[0]->sojourn_time[i]) ,
                                          *pfinal_run , occupancy_reestim , j);
              occupancy_reestim->distribution_estimation(occupancy);

              previous_likelihood = likelihood;
              likelihood = occupancy->state_occupancy_likelihood_computation(*(seq->characteristics[0]->sojourn_time[i]) ,
                                                                             *pfinal_run);

#             ifdef MESSAGE
              if ((j < 10) || ((j < 100) && (j % 10 == 0)) || ((j < 1000) && (j % 100 == 0)) || (j % 1000 == 0)) {
                os << STAT_label[STATL_ITERATION] << " " << j << "   "
                   << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
                   << STAT_label[STATL_SMOOTHNESS] << ": " << occupancy->second_difference_norm_computation() << endl;
              }
#             endif

            }
            while ((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (j < OCCUPANCY_NB_ITER) &&
                     ((likelihood - previous_likelihood) / -likelihood > OCCUPANCY_LIKELIHOOD_DIFF)) ||
                    ((nb_iter != I_DEFAULT) && (j < nb_iter))));

            if (likelihood != D_INF) {

#             ifdef MESSAGE
              os << "\n" << STAT_label[STATL_STATE] << " " << i << " - "
                 << j << " " << STAT_label[STATL_ITERATIONS] << "   "
                 << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
                 << STAT_label[STATL_SMOOTHNESS] << ": " << occupancy->second_difference_norm_computation() << "\n" << endl;
#             endif

              hreestim = new Histogram(MAX(occupancy->alloc_nb_value , max_length + 1));

              likelihood = D_INF;
              nb_likelihood_decrease = 0;

              j = 0;
              do {
                j++;

                // calcul des quantites de reestimation de la loi d'occupation de l'etat

                occupancy->expectation_step(*(seq->characteristics[0]->sojourn_time[i]) ,
                                            *pfinal_run , occupancy_reestim , j);

                hreestim->update(occupancy_reestim , (int)(occupancy_reestim->nb_element *
                                 MAX(sqrt(occupancy_reestim->variance) , 1.) * OCCUPANCY_COEFF));
                hlikelihood = hreestim->Reestimation<int>::type_parametric_estimation(occupancy , 1 , true ,
                                                                                      OCCUPANCY_THRESHOLD);

                if (hlikelihood == D_INF) {
                  likelihood = D_INF;
                }

                else {
                  occupancy->computation(hreestim->nb_value , OCCUPANCY_THRESHOLD);

                  previous_likelihood = likelihood;
                  likelihood = occupancy->state_occupancy_likelihood_computation(*(seq->characteristics[0]->sojourn_time[i]) ,
                                                                                 *pfinal_run);

                  if (likelihood < previous_likelihood) {
                    nb_likelihood_decrease++;
                  }
                  else {
                    nb_likelihood_decrease = 0;
                  }

#                 ifdef MESSAGE
                  if ((j < 10) || ((j < 100) && (j % 10 == 0)) || ((j < 1000) && (j % 100 == 0)) || (j % 1000 == 0)) {
                    os << STAT_label[STATL_ITERATION] << " " << j << "   "
                       << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
                       << STAT_label[STATL_SMOOTHNESS] << ": " << occupancy->second_difference_norm_computation() << endl;
                  }
#                 endif

                }
              }
              while ((likelihood != D_INF) && (j < OCCUPANCY_NB_ITER) &&
                     (((likelihood - previous_likelihood) / -likelihood > OCCUPANCY_LIKELIHOOD_DIFF) ||
                      (hlikelihood == D_INF) || (nb_likelihood_decrease == 1)));

              delete hreestim;

              if (likelihood != D_INF) {

#               ifdef MESSAGE
                os << "\n" << STAT_label[STATL_STATE] << " " << i << " - "
                   << j << " " << STAT_label[STATL_ITERATIONS] << "   "
                   << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
                   << STAT_label[STATL_SMOOTHNESS] << ": " << occupancy->second_difference_norm_computation() << "\n" << endl;
#               endif

              }

              else {
                delete occupancy;
                occupancy = 0;
              }
            }

            else {
              delete occupancy;
              occupancy = 0;
            }
          }

          else if ((estimator != PARTIAL_LIKELIHOOD) && (((model_type == 'o') && (seq->characteristics[0]->final_run[i]->nb_element > 0)) ||
                    ((model_type == 'e') && (final_run[i]->nb_element > 0)))) {
            seq->characteristics[0]->sojourn_time[i]->state_occupancy_estimation((model_type == 'o' ? seq->characteristics[0]->final_run[i] : final_run[i]) ,
                                                                                 occupancy_reestim , occupancy_survivor ,
                                                                                 censored_occupancy_survivor);

            occupancy = occupancy_reestim->type_parametric_estimation(1 , true , OCCUPANCY_THRESHOLD);

/*          occupancy = new Parametric(occupancy_reestim->nb_value);
            occupancy_reestim->distribution_estimation(occupancy); */
          }

          else {
            occupancy = seq->characteristics[0]->sojourn_time[i]->Reestimation<int>::type_parametric_estimation(1 , true ,
                                                                                                                OCCUPANCY_THRESHOLD);
          }

          if (occupancy) {
            if (occupancy->mean == 1.) {
              smarkov->state_subtype[i] = MARKOVIAN;
            }

            else {
              smarkov->state_subtype[i] = SEMI_MARKOVIAN;
              smarkov->nonparametric_process[0]->sojourn_time[i] = new Parametric(*occupancy);
              if (smarkov->state_type[i] == 'r') {
                smarkov->forward[i] = new Forward(*(smarkov->nonparametric_process[0]->sojourn_time[i]));
              }
            }

            delete occupancy;
          }

          else {
            delete smarkov;
            smarkov = 0;
            error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
            break;
          }
        }
      }

      if (estimator != PARTIAL_LIKELIHOOD) {
        delete [] occupancy_survivor;
        delete [] censored_occupancy_survivor;
        delete occupancy_reestim;
      }
    }

    if (model_type == 'e') {
      for (i = 0;i < marginal[0]->nb_value;i++) {
        delete initial_run[i];
      }
      delete [] initial_run;

      for (i = 0;i < marginal[0]->nb_value;i++) {
        delete final_run[i];
      }
      delete [] final_run;

      for (i = 0;i < marginal[0]->nb_value;i++) {
        delete single_run[i];
      }
      delete [] single_run;
    }

    if (smarkov) {
      if (model_type == 'e') {
        for (i = 0;i < smarkov->nb_state;i++) {
          smarkov->initial[i] = 1. / (double)smarkov->nb_state;
        }
        smarkov->initial_probability_computation();
      }

      // estimation des lois d'observation

      if (smarkov->nb_output_process == 1) {
        seq->build_observation_histogram();

        for (i = 0;i < smarkov->nb_state;i++) {
          seq->observation[1][i]->distribution_estimation(smarkov->nonparametric_process[1]->observation[i]);
        }
      }

      // calcul de la vraisemblance et des lois caracteristiques du modele

      seq->likelihood = smarkov->likelihood_computation(*seq);

#     ifdef DEBUG
      cout << "\n" << STAT_label[STATL_LIKELIHOOD] << ": " << seq->likelihood << " | "
           << smarkov->likelihood_computation(*seq , I_DEFAULT) << endl;
#     endif

      if (seq->likelihood == D_INF) {
        delete smarkov;
        smarkov = 0;
        error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
      }

      else {
        smarkov->characteristic_computation(*seq , counting_flag , I_DEFAULT , false);
      }
    }
  }

  return(smarkov);
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de differentes semi-chaines de Markov pour un ensemble
 *  de sequences.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre de semi-chaines
 *              de Markov, pointeur sur les semi-chaines de Markov, path.
 *
 *--------------------------------------------------------------*/

bool Markovian_sequences::comparison(Format_error &error , ostream &os , int nb_model ,
                                     const Semi_markov **ismarkov , const char *path) const

{
  bool status = true;
  register int i , j;
  double **likelihood;


  error.init();

  if ((type[0] != INT_VALUE) && (type[0] != STATE)) {
    status = false;
    ostringstream correction_message;
    correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
    error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
  }

  else if (!characteristics[0]) {
    for (i = 0;i < marginal[0]->nb_value;i++) {
      if (marginal[0]->frequency[i] == 0) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_error[SEQR_MISSING_STATE] << " " << i;
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (nb_variable > 1) {
    if (nb_variable > 2) {
      status = false;
      error.correction_update(STAT_error[STATR_NB_VARIABLE] , "1 or 2");
    }

    if ((type[1] != INT_VALUE) && (type[1] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << 2 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }

    else {
      if (test_hidden(1)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << 2 << ": "
                      << SEQ_error[SEQR_OVERLAP];
        error.update((error_message.str()).c_str());
      }

      if (!characteristics[1]) {
        for (i = 0;i < marginal[1]->nb_value;i++) {
          if (marginal[1]->frequency[i] == 0) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << 2 << ": "
                          << STAT_error[STATR_MISSING_VALUE] << " " << i;
            error.update((error_message.str()).c_str());
          }
        }
      }
    }
  }

  for (i = 0;i < nb_model;i++) {
    if (ismarkov[i]->nb_output_process + 1 != nb_variable) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_SEMI_MARKOV_CHAIN] << " " << i + 1 << ": "
                    << SEQ_error[SEQR_NB_OUTPUT_PROCESS];
      error.update((error_message.str()).c_str());
    }

    else {
      if (ismarkov[i]->nonparametric_process[0]->nb_value < marginal[0]->nb_value) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEMI_MARKOV_CHAIN] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_NB_STATE];
        error.update((error_message.str()).c_str());
      }

      if (nb_variable == 2) {
        if (ismarkov[i]->nonparametric_process[1]->nb_value < marginal[1]->nb_value) {
          status = false;
          ostringstream error_message;
          error_message << SEQ_label[SEQL_SEMI_MARKOV_CHAIN] << " " << i + 1 << ": "
                        << SEQ_error[SEQR_NB_OUTPUT];
          error.update((error_message.str()).c_str());
        }
      }
    }
  }

  if (status) {
    likelihood = new double*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      likelihood[i] = new double[nb_model];
    }

    // pour chaque sequence, calcul de la vraisemblance pour chaque modele possible

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < nb_model;j++) {
        likelihood[i][j] = ismarkov[j]->likelihood_computation(*this , i);
      }
    }

#   ifdef MESSAGE
    likelihood_write(os , nb_model , likelihood , SEQ_label[SEQL_SEMI_MARKOV_CHAIN] , true);
#   endif

    if (path) {
      status = likelihood_write(error , path , nb_model , likelihood , SEQ_label[SEQL_SEMI_MARKOV_CHAIN]);
    }

    for (i = 0;i < nb_sequence;i++) {
      delete [] likelihood[i];
    }
    delete [] likelihood;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une semi-chaine de Markov.
 *
 *  arguments : reference sur un objet Format_error,
 *              histogramme des longueurs des sequences,
 *              flag sur le calcul des lois de comptage,
 *              flag calcul d'une divergence de Kullback-Leibler.
 *
 *--------------------------------------------------------------*/

Semi_markov_data* Semi_markov::simulation(Format_error &error , const Histogram &hlength ,
                                          bool counting_flag , bool divergence_flag) const

{
  bool status = true;
  register int i , j , k , m;
  int cumul_length , occupancy , *pstate , **poutput;
  Semi_markov *smarkov;
  Semi_markov_data *seq;


  seq = 0;
  error.init();

  if ((hlength.nb_element < 1) || (hlength.nb_element > NB_SEQUENCE)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }
  if (hlength.offset < 2) {
    status = false;
    error.update(SEQ_error[SEQR_SHORT_SEQUENCE_LENGTH]);
  }
  if (hlength.nb_value - 1 > MAX_LENGTH) {
    status = false;
    error.update(SEQ_error[SEQR_LONG_SEQUENCE_LENGTH]);
  }

  if (status) {
    cumul_length = 0;
    for (i = hlength.offset;i < hlength.nb_value;i++) {
      cumul_length += i * hlength.frequency[i];
    }

    if (cumul_length > CUMUL_LENGTH) {
      status = false;
      error.update(SEQ_error[SEQR_CUMUL_SEQUENCE_LENGTH]);
    }
  }

  if (status) {
    if (hlength.nb_value - 1 > COUNTING_MAX_LENGTH) {
      counting_flag = false;
    }

    // initialisations

    seq = new Semi_markov_data(hlength , nb_output_process + 1);
    seq->type[0] = STATE;

    seq->semi_markov = new Semi_markov(*this , false);

    smarkov = seq->semi_markov;
    smarkov->create_cumul();
    smarkov->cumul_computation();

    if (smarkov->nb_output_process > 0) {
      poutput = new int*[smarkov->nb_output_process];
    }

    for (i = 0;i < seq->nb_sequence;i++) {
      pstate = seq->int_sequence[i][0];
      *pstate = cumul_method(smarkov->nb_state , smarkov->cumul_initial);

      for (j = 0;j < smarkov->nb_output_process;j++) {
        poutput[j] = seq->int_sequence[i][j + 1];
      }

      j = 0;
      do {
        if (j > 0) {
          pstate++;
          *pstate = cumul_method(smarkov->nb_state , smarkov->cumul_transition[*(pstate - 1)]);
        }

        switch (smarkov->state_subtype[*pstate]) {

        case SEMI_MARKOVIAN : {
          if ((smarkov->type == 'e') && (j == 0)) {
            occupancy = smarkov->forward[*pstate]->simulation();
          }
          else {
            occupancy = smarkov->nonparametric_process[0]->sojourn_time[*pstate]->simulation();
          }

          if (j + occupancy > seq->length[i]) {
            occupancy = seq->length[i] - j;
          }
          break;
        }

        case MARKOVIAN : {
          if (smarkov->transition[*pstate][*pstate] < 1.) {
            occupancy = 1;
          }
          else {
            occupancy = seq->length[i] - j;
          }
          break;
        }
        }

        j += occupancy;

        for (k = 1;k < occupancy;k++) {
          pstate++;
          *pstate = *(pstate - 1);
        }

        for (k = 0;k < occupancy;k++) {
          for (m = 0;m < smarkov->nb_output_process;m++) {
            if (smarkov->nonparametric_process[m + 1]) {
              *poutput[m]++ = smarkov->nonparametric_process[m + 1]->observation[*pstate]->simulation();
            }
            else {
              *poutput[m]++ = smarkov->parametric_process[m + 1]->observation[*pstate]->simulation();
            }
          }
        }
      }
      while (j < seq->length[i]);
    }

    smarkov->remove_cumul();

    if (smarkov->nb_output_process > 0) {
      delete [] poutput;
    }

    // extraction des caracteristiques des sequences simulees

    seq->min_value[0] = 0;
    seq->max_value[0] = nb_state - 1;
    seq->build_marginal_histogram(0);

    for (i = 1;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);
      seq->build_marginal_histogram(i);
    }

    seq->build_transition_count(smarkov);
    seq->build_observation_histogram();
    seq->build_characteristic(I_DEFAULT , true , (type == 'e' ? true : false));

/*    if ((seq->max_value[0] < nb_state - 1) || (!(seq->characteristics[0]))) {
      delete seq;
      seq = 0;
      error.update(SEQ_error[SEQR_STATES_NOT_REPRESENTED]);
    }

    else if (!divergence_flag) { */
    if (!divergence_flag) {
      smarkov->characteristic_computation(*seq , counting_flag);

      // calcul de la vraisemblance

      seq->likelihood = smarkov->likelihood_computation(*seq);
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une semi-chaine de Markov.
 *
 *  arguments : reference sur un objet Format_error,
 *              nombre et longueur des sequences,
 *              flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

Semi_markov_data* Semi_markov::simulation(Format_error &error , int nb_sequence ,
                                          int length , bool counting_flag) const

{
  bool status = true;
  Semi_markov_data *seq;


  seq = 0;
  error.init();

  if ((nb_sequence < 1) || (nb_sequence > NB_SEQUENCE)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }
  if (length < 2) {
    status = false;
    error.update(SEQ_error[SEQR_SHORT_SEQUENCE_LENGTH]);
  }
  if (length > MAX_LENGTH) {
    status = false;
    error.update(SEQ_error[SEQR_LONG_SEQUENCE_LENGTH]);
  }

  if (status) {
    Histogram hlength(length + 1);

    hlength.nb_element = nb_sequence;
    hlength.offset = length;
    hlength.max = nb_sequence;
    hlength.mean = length;
    hlength.variance = 0.;
    hlength.frequency[length] = nb_sequence;

    seq = simulation(error , hlength , counting_flag);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une semi-chaine de Markov.
 *
 *  arguments : reference sur un objet Format_error, nombre de sequences,
 *              reference sur un objet Markovian_sequences,
 *              flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

Semi_markov_data* Semi_markov::simulation(Format_error &error , int nb_sequence ,
                                          const Markovian_sequences &iseq , bool counting_flag) const

{
  Histogram *hlength;
  Semi_markov_data *seq;


  error.init();

  if ((nb_sequence < 1) || (nb_sequence > NB_SEQUENCE)) {
    seq = 0;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }

  else {
    hlength = iseq.hlength->frequency_scale(nb_sequence);

    seq = simulation(error , *hlength , counting_flag);
    delete hlength;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de semi-chaines de Markov par calcul de divergences de Kullback-Leibler.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre de semi-chaines
 *              de Markov, pointeur sur les semi-chaines de Markov,
 *              histogramme des longueurs des sequences, path.
 *
 *--------------------------------------------------------------*/

Distance_matrix* Semi_markov::divergence_computation(Format_error &error , ostream &os ,
                                                     int nb_model , const Semi_markov **ismarkov ,
                                                     Histogram **hlength , const char *path) const

{
  bool status = true , lstatus;
  register int i , j , k;
  int cumul_length;
  double ref_likelihood , target_likelihood , **likelihood;
  const Semi_markov **smarkov;
  Markovian_sequences *iseq , *seq;
  Semi_markov_data *simul_seq;
  Distance_matrix *dist_matrix;
  ofstream *out_file;


  dist_matrix = 0;
  error.init();

  for (i = 0;i < nb_model - 1;i++) {
    if (ismarkov[i]->type != type) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_SEMI_MARKOV_CHAIN] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_MODEL_TYPE];
      error.update((error_message.str()).c_str());
    }

    if (ismarkov[i]->nb_output_process == nb_output_process) {
      if (ismarkov[i]->nb_state != nb_state) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEMI_MARKOV_CHAIN] << " " << i + 2 << ": "
                      << SEQ_error[SEQR_NB_STATE];
        error.update((error_message.str()).c_str());
      }

      if (nb_output_process == 1) {
        if (ismarkov[i]->nonparametric_process[1]->nb_value != nonparametric_process[1]->nb_value) {
          status = false;
          ostringstream error_message;
          error_message << SEQ_label[SEQL_SEMI_MARKOV_CHAIN] << " " << i + 2 << ": "
                        << SEQ_error[SEQR_NB_OUTPUT];
          error.update((error_message.str()).c_str());
        }
      }
    }

    else if ((nb_output_process == 0) && (ismarkov[i]->nb_output_process == 1)) {
      if (ismarkov[i]->nonparametric_process[1]->nb_value != nb_state) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEMI_MARKOV_CHAIN] << " " << i + 2 << ": "
                      << SEQ_error[SEQR_NB_OUTPUT];
        error.update((error_message.str()).c_str());
      }
    }

    else {  //  if ((nb_output_process == 1) && (ismarkov[i]->nb_output_process == 0))
      if (ismarkov[i]->nb_state != nonparametric_process[1]->nb_value) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEMI_MARKOV_CHAIN] << " " << i + 2 << ": "
                      << SEQ_error[SEQR_NB_STATE];
        error.update((error_message.str()).c_str());
      }
    }
  }

  for (i = 0;i < nb_model;i++) {
    lstatus = true;

    if ((hlength[i]->nb_element < 1) || (hlength[i]->nb_element > NB_SEQUENCE)) {
      lstatus = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << " "
                    << i + 1 << ": "  << SEQ_error[SEQR_NB_SEQUENCE];
      error.update((error_message.str()).c_str());
    }
    if (hlength[i]->offset < 2) {
      lstatus = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << " "
                    << i + 1 << ": "  << SEQ_error[SEQR_SHORT_SEQUENCE_LENGTH];
      error.update((error_message.str()).c_str());
    }
    if (hlength[i]->nb_value - 1 > MAX_LENGTH) {
      lstatus = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << " "
                    << i + 1 << ": "  << SEQ_error[SEQR_LONG_SEQUENCE_LENGTH];
      error.update((error_message.str()).c_str());
    }

    if (!lstatus) {
      status = false;
    }

    else {
      cumul_length = 0;
      for (j = hlength[i]->offset;j < hlength[i]->nb_value;j++) {
        cumul_length += j * hlength[i]->frequency[j];
      }

      if (cumul_length > CUMUL_LENGTH) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << " "
                      << i + 1 << ": "  << SEQ_error[SEQR_CUMUL_SEQUENCE_LENGTH];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {
    out_file = 0;

    if (path) {
      out_file = new ofstream(path);

      if (!out_file) {
        error.update(STAT_error[STATR_FILE_NAME]);

#       ifdef MESSAGE
        os << error;
#       endif

      }
    }

    smarkov = new const Semi_markov*[nb_model];

    smarkov[0] = this;
    for (i = 1;i < nb_model;i++) {
      smarkov[i] = ismarkov[i - 1];
    }

    dist_matrix = new Distance_matrix(nb_model , SEQ_label[SEQL_SEMI_MARKOV_CHAIN]);

    for (i = 0;i < nb_model;i++) {

      // simulation d'un echantillon de sequences a partir d'une semi-chaine de Markov

      simul_seq = smarkov[i]->simulation(error , *hlength[i] , false , true);

      likelihood = new double*[simul_seq->nb_sequence];
      for (j = 0;j < simul_seq->nb_sequence;j++) {
        likelihood[j] = new double[nb_model];
      }

      ref_likelihood = 0.;
      for (j = 0;j < simul_seq->nb_sequence;j++) {
        likelihood[j][i] = smarkov[i]->likelihood_computation(*simul_seq , j);
        ref_likelihood += likelihood[j][i];
      }

      if (smarkov[i]->nb_output_process == 1) {
        iseq = simul_seq->remove_variable_1();
      }
      else {
        iseq = simul_seq;
      }

      // calcul des vraisemblances de l'echantillon pour chacune des semi-chaines de Markov

      for (j = 0;j < nb_model;j++) {
        if (j != i) {
          if (smarkov[j]->nb_output_process == 1) {
            seq = iseq->transcode(error , smarkov[j]->nonparametric_process[1]);
          }
          else {
            seq = iseq;
          }

          target_likelihood = 0.;
          for (k = 0;k < seq->nb_sequence;k++) {
            likelihood[k][j] = smarkov[j]->likelihood_computation(*seq , k);
            if (target_likelihood != D_INF) {
              if (likelihood[k][j] != D_INF) {
                target_likelihood += likelihood[k][j];
              }
              else {
                target_likelihood = D_INF;
              }
            }
          }

          if (target_likelihood != D_INF) {
            dist_matrix->update(i + 1 , j + 1 , ref_likelihood - target_likelihood , seq->cumul_length);
          }

          if (smarkov[j]->nb_output_process == 1) {
            delete seq;
          }
        }
      }

#     ifdef MESSAGE
      os << SEQ_label[SEQL_SEMI_MARKOV_CHAIN] << " " << i + 1 << ": " << simul_seq->nb_sequence << " "
         << SEQ_label[SEQL_SIMULATED] << " " << SEQ_label[simul_seq->nb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << endl;
      simul_seq->likelihood_write(os , nb_model , likelihood , SEQ_label[SEQL_SEMI_MARKOV_CHAIN]);
#     endif

      if (out_file) {
        *out_file << SEQ_label[SEQL_SEMI_MARKOV_CHAIN] << " " << i + 1 << ": " << simul_seq->nb_sequence << " "
                  << SEQ_label[SEQL_SIMULATED] << " " << SEQ_label[simul_seq->nb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << endl;
        simul_seq->likelihood_write(*out_file , nb_model , likelihood , SEQ_label[SEQL_SEMI_MARKOV_CHAIN]);
      }

      for (j = 0;j < simul_seq->nb_sequence;j++) {
        delete [] likelihood[j];
      }
      delete [] likelihood;

      if (smarkov[i]->nb_output_process == 1) {
        delete iseq;
      }
      delete simul_seq;
    }

    if (out_file) {
      out_file->close();
      delete out_file;
    }

    delete smarkov;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de semi-chaines de Markov par calcul de divergences de Kullback-Leibler.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre de semi-chaines
 *              de Markov, pointeur sur les semi-chaines de Markov,
 *              nombre et longueur des sequences, path.
 *
 *--------------------------------------------------------------*/

Distance_matrix* Semi_markov::divergence_computation(Format_error &error , ostream &os ,
                                                     int nb_model , const Semi_markov **smarkov ,
                                                     int nb_sequence , int length , const char *path) const

{
  bool status = true;
  register int i;
  Histogram **hlength;
  Distance_matrix *dist_matrix;


  dist_matrix = 0;
  error.init();

  if ((nb_sequence < 1) || (nb_sequence > NB_SEQUENCE)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }
  if (length < 2) {
    status = false;
    error.update(SEQ_error[SEQR_SHORT_SEQUENCE_LENGTH]);
  }
  if (length > MAX_LENGTH) {
    status = false;
    error.update(SEQ_error[SEQR_LONG_SEQUENCE_LENGTH]);
  }

  if (status) {
    hlength = new Histogram*[nb_model];

    hlength[0] = new Histogram(length + 1);

    hlength[0]->nb_element = nb_sequence;
    hlength[0]->offset = length;
    hlength[0]->max = nb_sequence;
    hlength[0]->mean = length;
    hlength[0]->variance = 0.;
    hlength[0]->frequency[length] = nb_sequence;

    for (i = 1;i < nb_model;i++) {
      hlength[i] = new Histogram(*hlength[0]);
    }

    dist_matrix = divergence_computation(error , os , nb_model , smarkov , hlength , path);

    for (i = 0;i < nb_model;i++) {
      delete hlength[i];
    }
    delete [] hlength;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de semi-chaines de Markov par calcul de divergences de Kullback-Leibler.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre de semi-chaines
 *              de Markov, pointeur sur les semi-chaines de Markov,
 *              pointeurs sur des objets Markovian_sequences, path.
 *
 *--------------------------------------------------------------*/

Distance_matrix* Semi_markov::divergence_computation(Format_error &error , ostream &os ,
                                                     int nb_model , const Semi_markov **smarkov ,
                                                     int nb_sequence , const Markovian_sequences **seq ,
                                                     const char *path) const

{
  register int i;
  Histogram **hlength;
  Distance_matrix *dist_matrix;


  error.init();

  if ((nb_sequence < 1) || (nb_sequence > NB_SEQUENCE)) {
    dist_matrix = 0;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }

  else {
    hlength = new Histogram*[nb_model];
    for (i = 0;i < nb_model;i++) {
      hlength[i] = seq[i]->hlength->frequency_scale(nb_sequence);
    }

    dist_matrix = divergence_computation(error , os , nb_model , smarkov , hlength , path);

    for (i = 0;i < nb_model;i++) {
      delete hlength[i];
    }
    delete [] hlength;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Semi_markov_iterator.
 *
 *  argument : pointeur sur un objet Semi_markov.
 *
 *--------------------------------------------------------------*/

Semi_markov_iterator::Semi_markov_iterator(Semi_markov *ismarkov)

{
  semi_markov = ismarkov;
  (semi_markov->nb_iterator)++;

  if ((!(semi_markov->cumul_initial)) || (!(semi_markov->cumul_transition))) {
    semi_markov->create_cumul();
    semi_markov->cumul_computation();
  }

  state = I_DEFAULT;
  occupancy = 0;
  counter = I_DEFAULT;
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Semi_markov_iterator.
 *
 *  argument : reference sur un objet Semi_markov_iterator.
 *
 *--------------------------------------------------------------*/

void Semi_markov_iterator::copy(const Semi_markov_iterator &it)

{
  semi_markov = it.semi_markov;
  (semi_markov->nb_iterator)++;

  state = it.state;
  occupancy = it.occupancy;
  counter = it.counter;
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Semi_markov_iterator.
 *
 *--------------------------------------------------------------*/

Semi_markov_iterator::~Semi_markov_iterator()

{
  (semi_markov->nb_iterator)--;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Semi_markov_iterator.
 *
 *  argument : reference sur un objet Semi_markov_iterator.
 *
 *--------------------------------------------------------------*/

Semi_markov_iterator& Semi_markov_iterator::operator=(const Semi_markov_iterator &it)

{
  if (&it != this) {
    (semi_markov->nb_iterator)--;
    copy(it);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une semi-chaine de Markov.
 *
 *  arguments : sequence, longueur de la sequence, flag initialisation.
 *
 *--------------------------------------------------------------*/

bool Semi_markov_iterator::simulation(int **int_seq , int length , bool initialization)

{
  bool status;


  if ((state == I_DEFAULT) && (!initialization)) {
    status = false;
  }

  else {
    register int i , j;
    int *pstate , **poutput;


    status = true;

    if (semi_markov->nb_output_process > 0) {
      poutput = new int*[semi_markov->nb_output_process];
    }

    if (initialization) {
      state = cumul_method(semi_markov->nb_state , semi_markov->cumul_initial);

      switch (semi_markov->state_subtype[state]) {

      case SEMI_MARKOVIAN : {
        switch (semi_markov->type) {
        case 'o' :
          occupancy = semi_markov->nonparametric_process[0]->sojourn_time[state]->simulation();
          break;
        case 'e' :
          occupancy = semi_markov->forward[state]->simulation();
          break;
        }
        break;
      }

      case MARKOVIAN : {
        if (semi_markov->transition[state][state] < 1.) {
          occupancy = 1;
        }
        break;
      }
      }

      counter = 0;
    }

    pstate = int_seq[0];
    for (i = 0;i < semi_markov->nb_output_process;i++) {
      poutput[i] = int_seq[i + 1];
    }

    for (i = 0;i < length;i++) {
      counter++;
      *pstate++ = state;
      for (j = 0;j < semi_markov->nb_output_process;j++) {
        if (semi_markov->nonparametric_process[j + 1]) {
          *poutput[j]++ = semi_markov->nonparametric_process[j + 1]->observation[state]->simulation();
        }
        else {
          *poutput[j]++ = semi_markov->parametric_process[j + 1]->observation[state]->simulation();
        }
      }

      if ((semi_markov->transition[state][state] < 1.) && (counter == occupancy)) {
        state = cumul_method(semi_markov->nb_state , semi_markov->cumul_transition[state]);

        switch (semi_markov->state_subtype[state]) {
        case SEMI_MARKOVIAN :
          occupancy = semi_markov->nonparametric_process[0]->sojourn_time[state]->simulation();
          break;
        case MARKOVIAN :
          occupancy = 1;
          break;
        }

        counter = 0;
      }
    }

    if (semi_markov->nb_output_process > 0) {
      delete [] poutput;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une semi-chaine de Markov.
 *
 *  arguments : longueur de la sequence, flag initialisation.
 *
 *--------------------------------------------------------------*/

int** Semi_markov_iterator::simulation(int length , bool initialization)

{
  register int i;
  int **int_seq;


  if ((state == I_DEFAULT) && (!initialization)) {
    int_seq = 0;
  }

  else {
    int_seq = new int*[semi_markov->nb_output_process + 1];
    for (i = 0;i <= semi_markov->nb_output_process;i++) {
      int_seq[i] = new int[length];
    }

    simulation(int_seq , length , initialization);
  }

  return int_seq;
}
