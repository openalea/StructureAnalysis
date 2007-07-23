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
#include <iomanip>
#include "stat_tool/stat_tools.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/regression.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequences.h"
#include "markov.h"
#include "sequence_label.h"
#include "tool/config.h"

using namespace std;


extern int column_width(int value);
extern int column_width(int nb_value , const double *value , double scale = 1.);
extern int cumul_method(int nb_value , const double *cumul , double scale = 1.);



/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance de sequences pour une chaine de Markov.
 *
 *  arguments : reference sur un objet Markovian_sequences, indice de la sequence.
 *
 *--------------------------------------------------------------*/

double Markov::likelihood_computation(const Markovian_sequences &seq , int index) const

{
  register int i , j , k;
  int *pstate , *mstate , **poutput , power[ORDER];
  double likelihood = 0. , proba , **ptransition;


  // verification de la compatibilite entre le modele et les donnees

  if (nb_output_process + 1 == seq.nb_variable) {
    for (i = 0;i <= nb_output_process;i++) {
      if ((seq.marginal[i]) && (process[i]->nb_value < seq.marginal[i]->nb_value)) {
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

    i = 1;
    for (j = 0;j < order;j++) {
      power[j] = i;
      i *= nb_state;
    }

    for (i = 0;i < seq.nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
        pstate = seq.sequence[i][0];

        proba = initial[*pstate];
        if (proba > 0.) {
          likelihood += log(proba);
        }
        else {
          likelihood = D_INF;
          break;
        }

        for (j = 0;j < nb_output_process;j++) {
          poutput[j] = seq.sequence[i][j + 1];

          proba = process[j + 1]->observation[*pstate]->mass[*poutput[j]];
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

        for (j = 1;j < seq.length[i];j++) {
          ptransition = transition;
          mstate = pstate;

          for (k = 0;k < MIN(j - 1 , order);k++) {
            ptransition += *mstate-- * power[order - 1 - k];
          }
          for (k = MIN(j - 1 , order);k < order;k++) {
            ptransition += *mstate * power[order - 1 - k];
          }

          proba = *(*ptransition + *++pstate);
          if (proba > 0.) {
            likelihood += log(proba);
          }
          else {
            likelihood = D_INF;
            break;
          }

          for (k = 0;k < nb_output_process;k++) {
            proba = process[k + 1]->observation[*pstate]->mass[*++poutput[k]];
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

    if (nb_output_process > 0) {
      delete [] poutput;
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance de sequences pour une chaine de Markov.
 *
 *  argument : reference sur un objet Markov_data.
 *
 *--------------------------------------------------------------*/

double Markov::likelihood_computation(const Markov_data &seq) const

{
  register int i , j;
  double buff , likelihood = 0.;


  // verification de la compatibilite entre le modele et les donnees

  if (nb_output_process + 1 == seq.nb_variable) {
    for (i = 0;i <= nb_output_process;i++) {
      if (process[i]->nb_value < seq.marginal[i]->nb_value) {
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
      for (i = 1;i <= nb_output_process;i++) {
        for (j = 0;j < nb_state;j++) {
          buff = process[i]->observation[j]->likelihood_computation(*(seq.observation[i][j]));

          if (buff != D_INF) {
            likelihood += buff;
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
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Comptage des etats initiaux et des transitions.
 *
 *  arguments : reference sur un objet Chain_data, ordre,
 *              flag pour prendre en compte le debut de la sequence.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::transition_count_computation(const Chain_data &chain_data ,
                                                       int order , bool begin) const

{
  register int i , j , k;
  int **ttransition , *pstate , *mstate , power[ORDER];


  for (i = 0;i < chain_data.nb_state;i++) {
    chain_data.initial[i] = 0;
  }

  for (i = 0;i < chain_data.nb_row;i++) {
    for (j = 0;j < chain_data.nb_state;j++) {
      chain_data.transition[i][j] = 0;
    }
  }

  // extraction des etats initiaux et des transitions

  i = 1;
  for (j = 0;j < order;j++) {
    power[j] = i;
    i *= chain_data.nb_state;
  }

  for (i = 0;i < nb_sequence;i++) {
    pstate = sequence[i][0];
    (chain_data.initial[*pstate])++;

    if (!begin) {
       pstate += order - 1;
    }
    for (j = (begin ? 1 : order);j < length[i];j++) {
      ttransition = chain_data.transition;
      mstate = pstate;
      for (k = 0;k < MIN(j - 1 , order);k++) {
        ttransition += *mstate-- * power[order - 1 - k];
      }
      for (k = j - 1;k < order;k++) {
        ttransition += *mstate * power[order - 1 - k];
      }

      (*(*ttransition + *++pstate))++;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction des comptages des etats initiaux et des transitions.
 *
 *  arguments : ordre, flag pour prendre en compte le debut de la sequence.
 *
 *--------------------------------------------------------------*/

void Markov_data::build_transition_count(int order , bool begin)

{
  chain_data = new Chain_data('o' , marginal[0]->nb_value ,
                              (int)pow((double)marginal[0]->nb_value , order));
  transition_count_computation(*chain_data , order , begin);
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une chaine de Markov a partir
 *  d'un echantillon de sequences.
 *
 *  arguments : reference sur un objet Format_error,
 *              ordre de la chaine de Markov, flags sur le calcul
 *              des lois de comptage et des lois caracteristiques.
 *
 *--------------------------------------------------------------*/

Markov* Markovian_sequences::markov_estimation(Format_error &error , int order ,
                                               bool counting_flag , bool characteristic_flag) const

{
  bool status = true;
  register int i;
  int nb_value[1];
  Markov *markov;
  Markov_data *seq;


  markov = 0;
  error.init();

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

  if ((order < 1) || (order > ORDER)) {
    status = false;
    error.update(SEQ_error[SEQR_ORDER]);
  }
  else {
    if ((int)pow((double)marginal[0]->nb_value , order + 1) > NB_PARAMETER) {
      status = false;
      error.update(SEQ_error[SEQR_NB_PARAMETER]);
    }
  }

  if (nb_variable > 1) {
    if (nb_variable > 2) {
      status = false;
      error.correction_update(STAT_error[STATR_NB_VARIABLE] , "1 or 2");
    }

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

  if (status) {
    if (nb_variable == 2) {
      nb_value[0] = marginal[1]->nb_value;
    }

    markov = new Markov(marginal[0]->nb_value , order , nb_variable - 1 , nb_value);
    markov->markov_data = new Markov_data(*this);

    seq = markov->markov_data;
    seq->state_variable_init();
    seq->build_transition_count(order);

    // estimation des parametres de la chaine de Markov

    seq->chain_data->estimation(*markov);

    // estimation des lois d'observation

    if (markov->nb_output_process == 1) {
      seq->build_observation_histogram();

      for (i = 0;i < markov->nb_state;i++) {
        seq->observation[1][i]->distribution_estimation(markov->process[1]->observation[i]);
      }
    }

    // calcul de la vraisemblance et des lois caracteristiques du modele

    seq->likelihood = markov->likelihood_computation(*seq);

    if (seq->likelihood == D_INF) {
      delete markov;
      markov = 0;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    else {
      if (characteristic_flag) {
        markov->component_computation();
        markov->characteristic_computation(*seq , counting_flag , I_DEFAULT , false);
      }
    }
  }

  return markov;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une chaine de Markov d'ordre 0
 *  a partir d'un echantillon de sequences.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Markov* Markovian_sequences::markov_order0_estimation(Format_error &error) const

{
  bool status = true;
  register int i , j;
  int nb_value[1];
  Markov *markov;
  Markov_data *seq;


  markov = 0;
  error.init();

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

  if (nb_variable > 1) {
    if (nb_variable > 2) {
      status = false;
      error.correction_update(STAT_error[STATR_NB_VARIABLE] , "1 or 2");
    }

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

  if (status) {
    if (nb_variable == 2) {
      nb_value[0] = marginal[1]->nb_value;
    }

    markov = new Markov(marginal[0]->nb_value , 1 , nb_variable - 1 , nb_value);
    markov->markov_data = new Markov_data(*this);

    seq = markov->markov_data;
    seq->state_variable_init();
    seq->build_transition_count(1);

    // estimation des parametres de la chaine de Markov

    for (i = 0;i < markov->nb_state;i++) {
/*      sum = seq->chain_data->initial[i];
      for (j = 0;j < markov->nb_state;j++) {
        sum += seq->chain_data->transition[j][i];
      }
      markov->initial[i] = (double)sum / (double)seq->cumul_length; */

      markov->initial[i] = (double)marginal[0]->frequency[i] / (double)marginal[0]->nb_element;
      for (j = 0;j < markov->nb_state;j++) {
        markov->transition[j][i] = markov->initial[i];
      }
    }

    // estimation des lois d'observation

    if (markov->nb_output_process == 1) {
      seq->build_observation_histogram();

      for (i = 0;i < markov->nb_state;i++) {
        seq->observation[1][i]->distribution_estimation(markov->process[1]->observation[i]);
      }
    }

    // calcul de la vraisemblance

    seq->likelihood = markov->likelihood_computation(*seq);

    if (seq->likelihood == D_INF) {
      delete markov;
      markov = 0;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }
  }

  return markov;
}


/*--------------------------------------------------------------*
 *
 *  Correction de la vraisemblance de sequences pour une chaine de Markov.
 *
 *  argument : reference sur un objet Markov_data.
 *
 *--------------------------------------------------------------*/

double Markov::likelihood_correction(const Markov_data &seq) const

{
  register int i , j;
  int *cinitial;
  double correction , *pinitial;


  correction = 0.;

  cinitial = seq.chain_data->initial;
  pinitial = initial;
  for (i = 0;i < nb_state;i++) {
    if (*cinitial > 0) {
      correction += *cinitial * log(*pinitial);
    }
    cinitial++;
    pinitial++;
  }

  if (nb_output_process > 0) {
    for (i = 0;i < seq.nb_sequence;i++) {
      for (j = 0;j < nb_output_process;j++) {
        correction += log(process[j + 1]->observation[seq.sequence[i][0][0]]->mass[seq.sequence[i][j + 1][0]]);
      }
    }
  }

  return correction;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une chaine de Markov a partir
 *  d'un echantillon de sequences.
 *
 *  arguments : reference sur un objet Format_error, stream,
 *              table de transcodage des symboles, type de penalisation (AIC(c)/BIC),
 *              ordre de la chaine de Markov, flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

Markov* Markovian_sequences::markov_lumpability_estimation(Format_error &error , ostream &os ,
                                                           int *symbol , int penalty_type ,
                                                           int order , bool counting_flag) const

{
  bool status = true , *presence;
  register int i;
  int max_symbol , nb_state[2] , nb_parameter[2];
  double penalty , max_likelihood , likelihood[2] , penalized_likelihood[2];
  Markov *markov , *lumped_markov;
  Markovian_sequences *seq;


  markov = 0;
  error.init();

  if (nb_variable > 1) {
    status = false;
    error.correction_update(STAT_error[STATR_NB_VARIABLE] , 1);
  }

  else {
    max_symbol = 0;
    for (i = 0;i < marginal[0]->nb_value;i++) {
      if ((symbol[i] < 0) || (symbol[i] >= marginal[0]->nb_value - 1)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_SYMBOL] << " " << symbol[i] << " "
                      << STAT_error[STATR_NOT_ALLOWED];
        error.update((error_message.str()).c_str());
      }
      else if (symbol[i] > max_symbol) {
        max_symbol = symbol[i];
      }
    }

    if (max_symbol == 0) {
      status = false;
      error.update(STAT_error[STATR_NB_SYMBOL]);
    }

    if (status) {
      presence = new bool[max_symbol + 1];
      for (i = 0;i <= max_symbol;i++) {
        presence[i] = false;
      }

      for (i = 0;i < marginal[0]->nb_value;i++) {
        presence[symbol[i]] = true;
      }

      for (i = 0;i <= max_symbol;i++) {
        if (!presence[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_error[STATR_MISSING_SYMBOL] << " " << i;
          error.update((error_message.str()).c_str());
        }
      }

      delete [] presence;
    }
  }

  if ((order < 1) || (order > ORDER)) {
    status = false;
    error.update(SEQ_error[SEQR_ORDER]);
  }

  if (status) {
    markov = markov_estimation(error , order , false , false);

    if (markov) {

      // calcul du terme de compensation

      switch (penalty_type) {
      case AIC :
        penalty = 1.;
        break;
      case BIC :
        penalty = 0.5 * log((double)cumul_length);
        break;
      }

      nb_state[1] = markov->nb_state;
      nb_parameter[1] = markov->nb_parameter_computation();
      likelihood[1] = markov->markov_data->likelihood -
                      markov->likelihood_correction(*(markov->markov_data));

      if (penalty_type == AICc) {
        if (nb_parameter[1] < cumul_length - 1) {
          penalized_likelihood[1] = likelihood[1] - (double)(nb_parameter[1] * cumul_length) /
                                    (double)(cumul_length - nb_parameter[1] - 1);
        }
        else {
          penalized_likelihood[1] = D_INF;
        }
      }

      else {
        penalized_likelihood[1] = likelihood[1] - nb_parameter[1] * penalty;
      }

      max_likelihood = penalized_likelihood[1];

      seq = transcode(error , 1 , symbol , true);

      lumped_markov = seq->markov_estimation(error , order , false , false);

      if (lumped_markov) {

#       ifdef MESSAGE
        {
          register int j , k;
          int nb_output , sum , lumped_nb_parameter , *pstate , *poutput , *pfrequency ,
              ***observation_data;
          double lumped_likelihood , lumped_penalized_likelihood;


          // 2eme propriete d'aggregation (probabilites d'observation attachees aux transitions)

          observation_data = new int**[seq->marginal[0]->nb_value];
          for (i = 0;i < seq->marginal[0]->nb_value;i++) {
            observation_data[i] = new int*[seq->marginal[1]->nb_value];
            for (j = 0;j < seq->marginal[1]->nb_value;j++) {
              observation_data[i][j] = new int[seq->marginal[1]->nb_value];
            }
            for (j = 0;j < seq->marginal[0]->nb_value;j++) {
              pfrequency = observation_data[i][j];
              for (k = 0;k < seq->marginal[1]->nb_value;k++) {
                *pfrequency++ = 0;
              }
            }
          }

          // accumulation des frequences d'observation

          for (i = 0;i < seq->nb_sequence;i++) {
            pstate = seq->sequence[i][0] + 1;
            poutput = seq->sequence[i][1] + 1;
            for (j = 1;j < seq->length[i];j++) {
              observation_data[*(pstate - 1)][*pstate][*poutput]++;
              pstate++;
              poutput++;
            }
          }

          // estimation des probabilites d'observation, calcul de la vraisemblance et
          // du nombre de parametres independants

          lumped_nb_parameter = lumped_markov->Chain::nb_parameter_computation();
          lumped_likelihood = lumped_markov->Chain::likelihood_computation(*(lumped_markov->markov_data->chain_data) , false);

          if (penalty_type == AICc) {
            if (lumped_nb_parameter < cumul_length - 1) {
              lumped_penalized_likelihood = lumped_likelihood - (double)(lumped_nb_parameter * cumul_length) /
                                            (double)(cumul_length - lumped_nb_parameter - 1);
            }
            else {
              lumped_penalized_likelihood = D_INF;
            }
          }

          else {
            lumped_penalized_likelihood = lumped_likelihood - lumped_nb_parameter * penalty;
          }

          os << "\n" << lumped_markov->nb_state << " " << STAT_label[STATL_STATES]
             << "   2 * " << SEQ_label[SEQL_MARKOV_CHAIN] << " " << STAT_label[STATL_LIKELIHOOD] << ": "
             << 2 * lumped_likelihood << "   " << lumped_nb_parameter << " "
             << STAT_label[lumped_nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
             << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("
             << STAT_criterion_word[penalty_type] << "): " << 2 * lumped_penalized_likelihood << endl;

          os << "\n" << STAT_word[STATW_OBSERVATION_PROBABILITIES] << endl;

          for (i = 0;i < seq->marginal[0]->nb_value;i++) {
            for (j = 0;j < seq->marginal[0]->nb_value;j++) {
              nb_output = 0;
              sum = 0;
              pfrequency = observation_data[i][j];
              for (k = 0;k < seq->marginal[1]->nb_value;k++) {
                if (*pfrequency > 0) {
                  nb_output++;
                  sum += *pfrequency;
                }
                pfrequency++;
              }

              if (nb_output > 1) {
                os << i << " -> " << j << " : ";

                lumped_nb_parameter += (nb_output - 1);

                pfrequency = observation_data[i][j];
                for (k = 0;k < seq->marginal[1]->nb_value;k++) {
                  if (*pfrequency > 0) {
                    os << k << " (" << (double)*pfrequency / (double)sum << ") | ";

                    lumped_likelihood += *pfrequency * log((double)*pfrequency / (double)sum);
                  }
                  pfrequency++;
                }

                os << endl;
              }
            }
          }

          if (penalty_type == AICc) {
            if (lumped_nb_parameter < cumul_length - 1) {
              lumped_penalized_likelihood = lumped_likelihood - (double)(lumped_nb_parameter * cumul_length) /
                                            (double)(cumul_length - lumped_nb_parameter - 1);
            }
            else {
              lumped_penalized_likelihood = D_INF;
            }
          }

          else {
            lumped_penalized_likelihood = lumped_likelihood - lumped_nb_parameter * penalty;
          }

          os << "\n" << lumped_markov->nb_state << " " << STAT_label[STATL_STATES]
             << "   2 * " << STAT_label[STATL_LIKELIHOOD] << ": " << 2 * lumped_likelihood << "   "
             << lumped_nb_parameter << " " << STAT_label[lumped_nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
             << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("
             << STAT_criterion_word[penalty_type] << "): " << 2 * lumped_penalized_likelihood << endl;

          // 3eme propriete d'aggregation (classique)

          for (i = 0;i < seq->marginal[0]->nb_value;i++) {
            for (j = 0;j < seq->marginal[1]->nb_value;j++) {
              pfrequency = observation_data[i][j];
              for (k = 0;k < seq->marginal[1]->nb_value;k++) {
                *pfrequency++ = 0;
              }
            }
          }

          // accumulation des frequences d'observation

          for (i = 0;i < seq->nb_sequence;i++) {
            pstate = seq->sequence[i][0] + 1;
            poutput = seq->sequence[i][1] + 1;
            for (j = 1;j < seq->length[i];j++) {
              observation_data[*pstate][*(poutput - 1)][*poutput]++;
              pstate++;
              poutput++;
            }
          }

          // estimation des probabilites d'observation, calcul de la vraisemblance et
          // du nombre de parametres independants

          lumped_nb_parameter = lumped_markov->Chain::nb_parameter_computation();
          lumped_likelihood = lumped_markov->Chain::likelihood_computation(*(lumped_markov->markov_data->chain_data) , false);

          os << "\n" << STAT_word[STATW_OBSERVATION_PROBABILITIES] << endl;

          for (i = 0;i < seq->marginal[0]->nb_value;i++) {
            for (j = 0;j < seq->marginal[1]->nb_value;j++) {
              nb_output = 0;
              sum = 0;
              pfrequency = observation_data[i][j];
              for (k = 0;k < seq->marginal[1]->nb_value;k++) {
                if (*pfrequency > 0) {
                  nb_output++;
                  sum += *pfrequency;
                }
                pfrequency++;
              }

              if (nb_output > 1) {
                os << j << ", " << i << " : ";

                lumped_nb_parameter += (nb_output - 1);

                pfrequency = observation_data[i][j];
                for (k = 0;k < seq->marginal[1]->nb_value;k++) {
                  if (*pfrequency > 0) {
                    os << k << " (" << (double)*pfrequency / (double)sum << ") | ";

                    lumped_likelihood += *pfrequency * log((double)*pfrequency / (double)sum);
                  }
                  pfrequency++;
                }

                os << endl;
              }
            }
          }

          if (penalty_type == AICc) {
            if (lumped_nb_parameter < cumul_length - 1) {
              lumped_penalized_likelihood = lumped_likelihood - (double)(lumped_nb_parameter * cumul_length) /
                                            (double)(cumul_length - lumped_nb_parameter - 1);
            }
            else {
              lumped_penalized_likelihood = D_INF;
            }
          }

          else {
            lumped_penalized_likelihood = lumped_likelihood - lumped_nb_parameter * penalty;
          }

          os << "\n" << lumped_markov->nb_state << " " << STAT_label[STATL_STATES]
             << "   2 * " << STAT_label[STATL_LIKELIHOOD] << ": " << 2 * lumped_likelihood << "   "
             << lumped_nb_parameter << " " << STAT_label[lumped_nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
             << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("
             << STAT_criterion_word[penalty_type] << "): " << 2 * lumped_penalized_likelihood << endl;

          for (i = 0;i < seq->marginal[0]->nb_value;i++) {
            for (j = 0;j < seq->marginal[1]->nb_value;j++) {
              delete [] observation_data[i][j];
            }
            delete [] observation_data[i];
          }
          delete [] observation_data;
        }
#       endif

        nb_state[0] = lumped_markov->nb_state;
        nb_parameter[0] = lumped_markov->nb_parameter_computation();
        likelihood[0] = lumped_markov->markov_data->likelihood -
                        lumped_markov->likelihood_correction(*(lumped_markov->markov_data));

        if (penalty_type == AICc) {
          if (nb_parameter[0] < cumul_length - 1) {
            penalized_likelihood[0] = likelihood[0] - (double)(nb_parameter[0] * cumul_length) /
                                      (double)(cumul_length - nb_parameter[0] - 1);
          }
          else {
            penalized_likelihood[0] = D_INF;
          }
        }

        else {
          penalized_likelihood[0] = likelihood[0] - nb_parameter[0] * penalty;
        }

#       ifdef DEBUG
        if (penalized_likelihood[0] > max_likelihood) {
          markov->ascii_write(os);
        }
        else {
          lumped_markov->ascii_write(os);
        }
#       endif

        if (penalized_likelihood[0] > max_likelihood) {
          max_likelihood = penalized_likelihood[0];
          delete markov;
          markov = lumped_markov;
        }
        else {
          delete lumped_markov;
        }

#       ifdef DEBUG
        lumpability_test(error , symbol , os , order);
#       endif

#       ifdef MESSAGE
        {
/*          double norm = 0. , weight[2];

          for (i = 0;i < 2;i++) {
            weight[i] = exp(penalized_likelihood[i] - max_likelihood);
            norm += weight[i];
          } */

          for (i = 0;i < 2;i++) {
            os << "\n" << nb_state[i] << " " << STAT_label[STATL_STATES]
               << "   2 * " << STAT_label[STATL_LIKELIHOOD] << ": " << 2 * likelihood[i] << "   "
               << nb_parameter[i] << " " << STAT_label[nb_parameter[i] == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
               << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("
               << STAT_criterion_word[penalty_type] << "): " << 2 * penalized_likelihood[i] << endl;
//               << "   " << STAT_label[STATL_WEIGHT] << ": " << weight[i] / norm << endl;
          }
        }
#       endif

      }

      delete seq;
    }

    // calcul des lois caracteristiques du modele

    markov->component_computation();
    markov->characteristic_computation(*(markov->markov_data) , counting_flag , I_DEFAULT , false);
  }

  return markov;
}


/*--------------------------------------------------------------*
 *
 *  Test d'aggregation d'etats pour une chaine de Markov.
 *
 *  arguments : reference sur un objet Format_error, stream,
 *              table de transcodage des symboles, ordre de la chaine de Markov.
 *
 *--------------------------------------------------------------*/

bool Markovian_sequences::lumpability_test(Format_error &error , ostream &os , int *symbol ,
                                           int order) const

{
  bool status = true , *presence;
  register int i , j;
  int max_symbol , df , sum , *ftransition , power[ORDER] , state_index[ORDER];
  double value , var1 , var2 , **ptransition;
  Test *test;
  Markov *markov , *lumped_markov;
  Markovian_sequences *seq;


  error.init();

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

  if (nb_variable > 1) {
    status = false;
    error.correction_update(STAT_error[STATR_NB_VARIABLE] , 1);
  }

  else {
    max_symbol = 0;
    for (i = 0;i < marginal[0]->nb_value;i++) {
      if ((symbol[i] < 0) || (symbol[i] >= marginal[0]->nb_value - 1)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_SYMBOL] << " " << symbol[i] << " "
                      << STAT_error[STATR_NOT_ALLOWED];
        error.update((error_message.str()).c_str());
      }
      else if (symbol[i] > max_symbol) {
        max_symbol = symbol[i];
      }
    }

    if (max_symbol == 0) {
      status = false;
      error.update(STAT_error[STATR_NB_SYMBOL]);
    }

    if (status) {
      presence = new bool[max_symbol + 1];
      for (i = 0;i <= max_symbol;i++) {
        presence[i] = false;
      }

      for (i = 0;i < marginal[0]->nb_value;i++) {
        presence[symbol[i]] = true;
      }

      for (i = 0;i <= max_symbol;i++) {
        if (!presence[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_error[STATR_MISSING_SYMBOL] << " " << i;
          error.update((error_message.str()).c_str());
        }
      }

      delete [] presence;
    }
  }

  if ((order < 1) || (order > ORDER)) {
    status = false;
    error.update(SEQ_error[SEQR_ORDER]);
  }

  if (status) {
    markov = markov_estimation(error , order , false , false);

    seq = transcode(error , 1 , symbol , true);

    lumped_markov = seq->markov_estimation(error , order , false , false);

    df = markov->nb_parameter_computation() - lumped_markov->nb_parameter_computation();

    i = 1;
    for (j = 0;j < markov->order;j++) {
      power[j] = i;
      i *= lumped_markov->nb_state;
    }

    for (i = 0;i < markov->order;i++) {
      state_index[i] = 0;
    }

    value = 0.;

    for (i = 0;i < markov->nb_row;i++) {
      ftransition = markov->markov_data->chain_data->transition[i];
      sum = 0;
      for (j = 0;j < markov->nb_state;j++) {
        sum += *ftransition++;
      }

      if (sum > 0) {
        ptransition = lumped_markov->transition;
        for (j = 0;j < markov->order;j++) {
          ptransition += symbol[state_index[j]] * power[j];
        }

        ftransition = markov->markov_data->chain_data->transition[i];
        for (j = 0;j < markov->nb_state;j++) {
          var1 = (double)sum * lumped_markov->process[1]->observation[symbol[j]]->mass[j] *
                 *(*ptransition + symbol[j]);
          if (var1 > 0.) {
            var2 = *ftransition - var1;
            value += var2 * var2 / var1;
          }
          ftransition++;
        }
      }

      // mise a jour des indices des etats

      for (j = 0;j < markov->order;j++) {
        if (state_index[j] < markov->nb_state - 1) {
          state_index[j]++;
          break;
        }
        else {
          state_index[j] = 0;
        }
      }
    }

    test = new Test(CHI2 , true , df , I_DEFAULT , value);

    test->chi2_critical_probability_computation();

#   ifdef MESSAGE
    os << *test;
#   endif

    delete test;

    value = 2 * (markov->markov_data->likelihood - markov->likelihood_correction(*(markov->markov_data)) -
            (lumped_markov->markov_data->likelihood - lumped_markov->likelihood_correction(*(lumped_markov->markov_data))));

    test = new Test(CHI2 , true , df , I_DEFAULT , value);

    test->chi2_critical_probability_computation();

#   ifdef MESSAGE
    os << "\n" << STAT_label[STATL_LIKELIHOOD_RATIO_TEST] << "\n" << *test;
#   endif

    delete test;

#   ifdef MESSAGE
    {
      register int k;
      int nb_output , sum , lumped_nb_parameter , *pstate , *poutput , *pfrequency ,
          ***observation_data;
      double lumped_likelihood , *pproba , ***observation_proba;


      // 2eme propriete d'aggregation (probabilites d'observation attachees aux transitions)

      observation_data = new int**[seq->marginal[0]->nb_value];
      observation_proba = new double**[seq->marginal[0]->nb_value];
      for (i = 0;i < seq->marginal[0]->nb_value;i++) {
        observation_data[i] = new int*[seq->marginal[1]->nb_value];
        observation_proba[i] = new double*[seq->marginal[1]->nb_value];
        for (j = 0;j < seq->marginal[1]->nb_value;j++) {
          observation_data[i][j] = new int[seq->marginal[1]->nb_value];
          observation_proba[i][j] = new double[seq->marginal[1]->nb_value];
        }
        for (j = 0;j < seq->marginal[0]->nb_value;j++) {
          pfrequency = observation_data[i][j];
          for (k = 0;k < seq->marginal[1]->nb_value;k++) {
            *pfrequency++ = 0;
          }
        }
      }

      // accumulation des frequences d'observation

      for (i = 0;i < seq->nb_sequence;i++) {
        pstate = seq->sequence[i][0] + 1;
        poutput = seq->sequence[i][1] + 1;
        for (j = 1;j < seq->length[i];j++) {
          observation_data[*(pstate - 1)][*pstate][*poutput]++;
          pstate++;
          poutput++;
        }
      }

      // estimation des probabilites d'observation, calcul de la vraisemblance et
      // du nombre de parametres independants

      lumped_nb_parameter = lumped_markov->Chain::nb_parameter_computation();
      lumped_likelihood = lumped_markov->Chain::likelihood_computation(*(lumped_markov->markov_data->chain_data) , false);

      for (i = 0;i < seq->marginal[0]->nb_value;i++) {
        for (j = 0;j < seq->marginal[0]->nb_value;j++) {
          nb_output = 0;
          sum = 0;
          pfrequency = observation_data[i][j];
          for (k = 0;k < seq->marginal[1]->nb_value;k++) {
            if (*pfrequency > 0) {
              nb_output++;
              sum += *pfrequency;
            }
            pfrequency++;
          }

          if (nb_output > 1) {
            lumped_nb_parameter += (nb_output - 1);

            pfrequency = observation_data[i][j];
            for (k = 0;k < seq->marginal[1]->nb_value;k++) {
              if (*pfrequency > 0) {
                lumped_likelihood += *pfrequency * log((double)*pfrequency / (double)sum);
              }
              pfrequency++;
            }
          }

          if (sum > 0) {
            pproba = observation_proba[i][j];
            pfrequency = observation_data[i][j];
            for (k = 0;k < seq->marginal[1]->nb_value;k++) {
              *pproba++ = (double)*pfrequency++ / (double)sum;
            }
          }
        }
      }

      df = markov->nb_parameter_computation() - lumped_nb_parameter;

      i = 1;
      for (j = 0;j < markov->order;j++) {
        power[j] = i;
        i *= lumped_markov->nb_state;
      }

      for (i = 0;i < markov->order;i++) {
        state_index[i] = 0;
      }

      value = 0.;

      for (i = 0;i < markov->nb_row;i++) {
        ftransition = markov->markov_data->chain_data->transition[i];
        sum = 0;
        for (j = 0;j < markov->nb_state;j++) {
          sum += *ftransition++;
        }

        if (sum > 0) {
          ptransition = lumped_markov->transition;
          for (j = 0;j < markov->order;j++) {
            ptransition += symbol[state_index[j]] * power[j];
          }

          ftransition = markov->markov_data->chain_data->transition[i];
          for (j = 0;j < markov->nb_state;j++) {
            var1 = (double)sum * observation_proba[symbol[state_index[0]]][symbol[j]][j] *
                   *(*ptransition + symbol[j]);
            if (var1 > 0.) {
              var2 = *ftransition - var1;
              value += var2 * var2 / var1;
            }
            ftransition++;
          }
        }

        // mise a jour des indices des etats

        for (j = 0;j < markov->order;j++) {
          if (state_index[j] < markov->nb_state - 1) {
            state_index[j]++;
            break;
          }
          else {
            state_index[j] = 0;
          }
        }
      }

      test = new Test(CHI2 , true , df , I_DEFAULT , value);

      test->chi2_critical_probability_computation();

      os << "\n" << *test;

      delete test;

      value = 2 * (markov->markov_data->likelihood - markov->likelihood_correction(*(markov->markov_data)) - lumped_likelihood);

      test = new Test(CHI2 , true , df , I_DEFAULT , value);

      test->chi2_critical_probability_computation();

      os << "\n" << STAT_label[STATL_LIKELIHOOD_RATIO_TEST] << "\n" << *test;

      delete test;

      // 3eme propriete d'aggregation (classique)

      for (i = 0;i < seq->marginal[0]->nb_value;i++) {
        for (j = 0;j < seq->marginal[1]->nb_value;j++) {
          pfrequency = observation_data[i][j];
          for (k = 0;k < seq->marginal[1]->nb_value;k++) {
            *pfrequency++ = 0;
          }
        }
      }

      // accumulation des frequences d'observation

      for (i = 0;i < seq->nb_sequence;i++) {
        pstate = seq->sequence[i][0] + 1;
        poutput = seq->sequence[i][1] + 1;
        for (j = 1;j < seq->length[i];j++) {
          observation_data[*pstate][*(poutput - 1)][*poutput]++;
          pstate++;
          poutput++;
        }
      }

      // estimation des probabilites d'observation, calcul de la vraisemblance et
      // du nombre de parametres independants

      lumped_nb_parameter = lumped_markov->Chain::nb_parameter_computation();
      lumped_likelihood = lumped_markov->Chain::likelihood_computation(*(lumped_markov->markov_data->chain_data) , false);

      for (i = 0;i < seq->marginal[0]->nb_value;i++) {
        for (j = 0;j < seq->marginal[1]->nb_value;j++) {
          nb_output = 0;
          sum = 0;
          pfrequency = observation_data[i][j];
          for (k = 0;k < seq->marginal[1]->nb_value;k++) {
            if (*pfrequency > 0) {
              nb_output++;
              sum += *pfrequency;
            }
            pfrequency++;
          }

          if (nb_output > 1) {
            lumped_nb_parameter += (nb_output - 1);

            pfrequency = observation_data[i][j];
            for (k = 0;k < seq->marginal[1]->nb_value;k++) {
              if (*pfrequency > 0) {
                lumped_likelihood += *pfrequency * log((double)*pfrequency / (double)sum);
              }
              pfrequency++;
            }
          }

          if (sum > 0) {
            pproba = observation_proba[i][j];
            pfrequency = observation_data[i][j];
            for (k = 0;k < seq->marginal[1]->nb_value;k++) {
              *pproba++ = (double)*pfrequency++ / (double)sum;
            }
          }
        }
      }

      df = markov->nb_parameter_computation() - lumped_nb_parameter;

      i = 1;
      for (j = 0;j < markov->order;j++) {
        power[j] = i;
        i *= lumped_markov->nb_state;
      }

      for (i = 0;i < markov->order;i++) {
        state_index[i] = 0;
      }

      value = 0.;

      for (i = 0;i < markov->nb_row;i++) {
        ftransition = markov->markov_data->chain_data->transition[i];
        sum = 0;
        for (j = 0;j < markov->nb_state;j++) {
          sum += *ftransition++;
        }

        if (sum > 0) {
          ptransition = lumped_markov->transition;
          for (j = 0;j < markov->order;j++) {
            ptransition += symbol[state_index[j]] * power[j];
          }

          ftransition = markov->markov_data->chain_data->transition[i];
          for (j = 0;j < markov->nb_state;j++) {
            var1 = (double)sum * observation_proba[symbol[j]][state_index[0]][j] *
                   *(*ptransition + symbol[j]);
            if (var1 > 0.) {
              var2 = *ftransition - var1;
              value += var2 * var2 / var1;
            }
            ftransition++;
          }
        }

        // mise a jour des indices des etats

        for (j = 0;j < markov->order;j++) {
          if (state_index[j] < markov->nb_state - 1) {
            state_index[j]++;
            break;
          }
          else {
            state_index[j] = 0;
          }
        }
      }

      test = new Test(CHI2 , true , df , I_DEFAULT , value);

      test->chi2_critical_probability_computation();

      os << "\n" << *test;

      delete test;

      value = 2 * (markov->markov_data->likelihood - markov->likelihood_correction(*(markov->markov_data)) - lumped_likelihood);

      test = new Test(CHI2 , true , df , I_DEFAULT , value);

      test->chi2_critical_probability_computation();

      os << "\n" << STAT_label[STATL_LIKELIHOOD_RATIO_TEST] << "\n" << *test;

      delete test;

      for (i = 0;i < seq->marginal[0]->nb_value;i++) {
        for (j = 0;j < seq->marginal[1]->nb_value;j++) {
          delete [] observation_data[i][j];
          delete [] observation_proba[i][j];
        }
        delete [] observation_data[i];
        delete [] observation_proba[i];
      }
      delete [] observation_data;
      delete [] observation_proba;
    }
#   endif

    delete markov;
    delete seq;
    delete lumped_markov;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Estimation de l'ordre et des parametres d'une chaine de Markov
 *  a partir d'un echantillon de sequences.
 *
 *  arguments : reference sur un objet Format_error, stream,
 *              ordre maximum de la chaine de Markov, type de penalisation (AIC(c)/BIC),
 *              flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

Markov* Markovian_sequences::markov_order_estimation(Format_error &error , ostream &os ,
                                                     int max_order , int penalty_type ,
                                                     bool counting_flag) const

{
  register int i;
  int nb_parameter[ORDER + 1];
  double penalty , max_likelihood , likelihood[ORDER + 1] , penalized_likelihood[ORDER + 1];
  Markov *markov , *pmarkov;


  error.init();

  if ((max_order < 1) || (max_order > ORDER)) {
    markov = 0;
    error.update(SEQ_error[SEQR_ORDER]);
  }

  else {

    // estimation des parametres d'une chaine de Markov d'ordre 0

    markov = markov_order0_estimation(error);

    if (markov) {

      // calcul du terme de compensation

      switch (penalty_type) {
      case AIC :
        penalty = 1.;
        break;
      case BIC :
        penalty = 0.5 * log((double)cumul_length);
        break;
      }

      likelihood[0] = markov->markov_data->likelihood;
      nb_parameter[0] = markov->nb_state - 1;

      if (penalty_type == AICc) {
        if (nb_parameter[0] < cumul_length - 1) {
          penalized_likelihood[0] = likelihood[0] - (double)(nb_parameter[0] * cumul_length) /
                                    (double)(cumul_length - nb_parameter[0] - 1);
        }
        else {
          penalized_likelihood[0] = D_INF;
        }
      }

      else {
        penalized_likelihood[0] = likelihood[0] - nb_parameter[0] * penalty;
      }

      max_likelihood = penalized_likelihood[0];

      for (i = 1;i <= max_order;i++) {
        pmarkov = markov_estimation(error , i , false , false);

        if (pmarkov) {
          likelihood[i] = pmarkov->markov_data->likelihood;
          nb_parameter[i] = pmarkov->nb_parameter_computation();

          if (penalty_type == AICc) {
            if (nb_parameter[i] < cumul_length - 1) {
              penalized_likelihood[i] = likelihood[i] - (double)(nb_parameter[i] * cumul_length) /
                                        (double)(cumul_length - nb_parameter[i] - 1);
            }
            else {
              penalized_likelihood[i] = D_INF;
            }
          }

          else {
            penalized_likelihood[i] = likelihood[i] - nb_parameter[i] * penalty;
          }

#         ifdef DEBUG
          order_test(error , i , i - 1 , os);
#         endif

          if (penalized_likelihood[i] > max_likelihood) {
            max_likelihood = penalized_likelihood[i];
            delete markov;
            markov = pmarkov;
          }
          else {
            delete pmarkov;
          }
        }

        else {
          likelihood[i] = D_INF;
        }
      }

#     ifdef MESSAGE
      {
        double norm = 0. , weight[ORDER + 1];

        for (i = 0;i <= max_order;i++) {
          if (likelihood[i] != D_INF) {
            weight[i] = exp(penalized_likelihood[i] - max_likelihood);
            norm += weight[i];
          }
        }

        for (i = 0;i <= max_order;i++) {
          if (likelihood[i] != D_INF) {
            os << "\n" << STAT_label[STATL_ORDER] << " " << i
               << "   2 * " << STAT_label[STATL_LIKELIHOOD] << ": " << 2 * likelihood[i] << "   "
               << nb_parameter[i] << " " << STAT_label[nb_parameter[i] == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
               << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("
               << STAT_criterion_word[penalty_type] << "): " << 2 * penalized_likelihood[i] << "   "
               << STAT_label[STATL_WEIGHT] << ": " << weight[i] / norm << endl;
          }
        }
      }
#     endif

      // calcul des lois caracteristiques du modele

      markov->component_computation();
      markov->characteristic_computation(*(markov->markov_data) , counting_flag , I_DEFAULT , false);
    }
  }

  return markov;
}


/*--------------------------------------------------------------*
 *
 *  Test de l'ordre d'une chaine de Markov.
 *
 *  arguments : reference sur un objet Format_error, stream, ordres.
 *
 *--------------------------------------------------------------*/

bool Markovian_sequences::order_test(Format_error &error , ostream &os , int order1 ,
                                     int order2) const

{
  bool status = true;
  register int i , j;
  int df , power , sum1 , sum2 , *ptransition1 , *ptransition2;
  double value , var1 , var2;
  Test *test;
  Chain_data *chain_data1 , *chain_data2;


  error.init();

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

  if (nb_variable > 1) {
    if (nb_variable > 2) {
      status = false;
      error.correction_update(STAT_error[STATR_NB_VARIABLE] , "1 or 2");
    }

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

  if ((order1 < 1) || (order1 > ORDER)) {
    status = false;
    error.update(SEQ_error[SEQR_ORDER]);
  }
  if ((order2 < 0) || (order2 > ORDER - 1)) {
    status = false;
    error.update(SEQ_error[SEQR_ORDER]);
  }
  if (order1 <= order2) {
    status = false;
    error.update(SEQ_error[SEQR_ORDER]);
  }

  if (status) {
    chain_data1 = new Chain_data('o' , marginal[0]->nb_value ,
                                 (int)pow((double)marginal[0]->nb_value , order1));
    transition_count_computation(*chain_data1 , order1);
    chain_data2 = new Chain_data(*chain_data1 , order1 , order2);

    df = chain_data1->nb_parameter_computation() - chain_data2->nb_parameter_computation();

    power = 1;
    for (i = order2;i < order1;i++) {
      power *= chain_data2->nb_state;
    }

    value = 0.;
    sum2 = I_DEFAULT;

    for (i = 0;i < chain_data1->nb_row;i++) {
      ptransition1 = chain_data1->transition[i];
      sum1 = 0;
      for (j = 0;j < chain_data1->nb_state;j++) {
        sum1 += *ptransition1++;
      }

      if (sum1 > 0) {
        if (sum2 == I_DEFAULT) {
          ptransition2 = chain_data2->transition[i / power];
          sum2 = 0;
          for (j = 0;j < chain_data2->nb_state;j++) {
            sum2 += *ptransition2++;
          }
        }

        ptransition2 = chain_data2->transition[i / power];
        ptransition1 = chain_data1->transition[i];
        for (j = 0;j < chain_data2->nb_state;j++) {
          var1 = (double)sum1 * *ptransition2++ / (double)sum2;
          if (var1 > 0.) {
            var2 = *ptransition1 - var1;
            value += var2 * var2 / var1;
          }
          ptransition1++;
        }
      }

      if ((i + 1) % power == 0) {
        sum2 = I_DEFAULT;
      }
    }

    delete chain_data1;
    delete chain_data2;

    test = new Test(CHI2 , true , df , I_DEFAULT , value);

    test->chi2_critical_probability_computation();

#   ifdef MESSAGE
    os << *test;
#   endif

    delete test;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des comptage des transitions pour differents ordres.
 *
 *  arguments : stream, ordre maximum, pointeurs sur les objets Chain_data,
 *              flag pour prendre en compte le debut de la sequence.
 *
 *--------------------------------------------------------------*/

ostream& Chain_data::transition_count_ascii_write(ostream &os , int max_order ,
                                                  const Chain_data **chain_data , bool begin) const

{
  register int i , j , k , m;
  int buff , width , initial_count , child_nb_parameter , *pfrequency , **memory_count ,
      **nb_parameter , state_index[ORDER];
  long old_adjust;
  double initial_likelihood , child_likelihood , **memory_likelihood;
  Distribution *marginal_dist;
  Histogram *marginal_histo;
  Chain **chain;
  const Chain_data **pchain_data;


  old_adjust = os.setf(ios::right , ios::adjustfield);

  pchain_data = new const Chain_data*[max_order + 1];
  pchain_data[1] = this;
  for (i = 2;i <= max_order;i++) {
    pchain_data[i] = chain_data[i - 2];
  }

  marginal_histo = new Histogram(pchain_data[1]->nb_state);

  memory_count = new int*[max_order + 1];
  memory_count[0] = new int[1];
  for (i = 1;i <= max_order;i++) {
    memory_count[i] = new int[pchain_data[i]->nb_row];
  }

  memory_likelihood = new double*[max_order + 1];
  memory_likelihood[0] = new double[1];
  for (i = 1;i <= max_order;i++) {
    memory_likelihood[i] = new double[pchain_data[i]->nb_row];
  }

  nb_parameter = new int*[max_order + 1];
  nb_parameter[0] = new int[1];
  for (i = 1;i <= max_order;i++) {
    nb_parameter[i] = new int[pchain_data[i]->nb_row];
  }

  pfrequency = marginal_histo->frequency;
  for (i = 0;i < pchain_data[1]->nb_state;i++) {
    switch (begin) {
    case false :
      *pfrequency = 0;
      break;
    case true :
      *pfrequency = pchain_data[1]->initial[i];
      break;
    }

    for (j = 0;j < pchain_data[1]->nb_state;j++) {
      *pfrequency += pchain_data[1]->transition[j][i];
    }
    pfrequency++;
  }

  marginal_histo->nb_element_computation();
  marginal_histo->max_computation();

  width = column_width(marginal_histo->nb_element) + ASCII_SPACE;

  if (begin) {
    os << "\n" << SEQ_label[SEQL_INITIAL_COUNTS] << endl;

    for (i = 0;i < max_order;i++) {
      os << "  ";
    }
    os << "  ";

    for (i = 0;i < pchain_data[1]->nb_state;i++) {
      os << setw(width) << i;
    }
    os << endl;

    for (i = 0;i < max_order;i++) {
      os << "  ";
    }
    os << "  ";

    initial_count = 0;
    for (i = 0;i < pchain_data[1]->nb_state;i++) {
      initial_count += pchain_data[1]->initial[i];
      os << setw(width) << pchain_data[1]->initial[i];
    }
    os << "  " << setw(width) << initial_count << endl;
  }

  os << "\n" << SEQ_label[SEQL_TRANSITION_COUNTS] << endl;

  for (i = 0;i < max_order;i++) {
    os << "  ";
  }
  os << "  ";

  for (i = 0;i < pchain_data[1]->nb_state;i++) {
    os << setw(width) << i;
  }
  os << endl;

  for (i = 0;i < max_order;i++) {
    os << "  ";
  }
  os << "  ";

  pfrequency = marginal_histo->frequency;
  for (i = 0;i < marginal_histo->nb_value;i++) {
    os << setw(width) << *pfrequency++;
  }
  os << "  " << setw(width) << marginal_histo->nb_element << endl;

  memory_count[0][0] = marginal_histo->nb_element;

  for (i = 1;i <= max_order;i++) {
    os << "\n";
    for (j = 0;j < i;j++) {
      state_index[j] = 0;
    }

    for (j = 0;j < pchain_data[i]->nb_row;j++) {
      for (k = i;k < max_order;k++) {
        os << "  ";
      }
      for (k = 0;k < i;k++) {
        os << state_index[k] << " ";
      }
      os << "  ";

      memory_count[i][j] = 0;
      for (k = 0;k < pchain_data[i]->nb_state;k++) {
        memory_count[i][j] += pchain_data[i]->transition[j][k];
        os << setw(width) << pchain_data[i]->transition[j][k];
      }
      os << "  " << setw(width) << memory_count[i][j] << endl;

      // mise a jour des indices des etats

      for (k = 0;k < i;k++) {
        if (state_index[k] < pchain_data[i]->nb_state - 1) {
          state_index[k]++;
          break;
        }
        else {
          state_index[k] = 0;
        }
      }
    }
  }
  os << endl;

  marginal_dist = new Distribution(*marginal_histo);

  chain = new Chain*[max_order + 1];
  for (i = 1;i <= max_order;i++) {
    chain[i] = new Chain('o' , pchain_data[i]->nb_state ,
                         (int)pow((double)pchain_data[i]->nb_state , i) , true);
    pchain_data[i]->estimation(*chain[i]);
  }

  // calcul des largeurs des colonnes

  width = column_width(marginal_dist->nb_value , marginal_dist->mass);
  for (i = 1;i <= max_order;i++) {
    for (j = 0;j < chain[i]->nb_row;j++) {
      buff = column_width(chain[i]->nb_state , chain[i]->transition[j]);
      if (buff > width) {
        width = buff;
      }
    }
  }
  width += ASCII_SPACE;

  if (begin) {
    os << "\n" << STAT_word[STATW_INITIAL_PROBABILITIES] << endl;

    for (i = 0;i < max_order;i++) {
      os << "  ";
    }
    os << "  ";

    for (i = 0;i < chain[1]->nb_state;i++) {
      os << setw(width) << i;
    }
    os << endl;

    for (i = 0;i < max_order;i++) {
      os << "  ";
    }
    os << "  ";

    for (i = 0;i < chain[1]->nb_state;i++) {
      os << setw(width) << chain[1]->initial[i];
    }
    os << endl;
  }

  os << "\n" << STAT_word[STATW_TRANSITION_PROBABILITIES] << endl;

  for (i = 0;i < max_order;i++) {
    os << "  ";
  }
  os << "  ";

  for (i = 0;i < chain[1]->nb_state;i++) {
    os << setw(width) << i;
  }
  os << endl;

  for (i = 0;i < max_order;i++) {
    os << "  ";
  }
  os << "  ";

  for (i = 0;i < marginal_dist->nb_value;i++) {
    os << setw(width) << marginal_dist->mass[i];

    if (marginal_histo->frequency[i] > 0) {
      marginal_dist->mass[i] = marginal_histo->frequency[i] * log(marginal_dist->mass[i]);
    }
    else {
      marginal_dist->mass[i] = 0.;
    }
  }
  os << endl;

  for (i = 1;i <= max_order;i++) {
    os << "\n";
    for (j = 0;j < i;j++) {
      state_index[j] = 0;
    }

    for (j = 0;j < chain[i]->nb_row;j++) {
      for (k = i;k < max_order;k++) {
        os << "  ";
      }
      for (k = 0;k < i;k++) {
        os << state_index[k] << " ";
      }
      os << "  ";

      for (k = 0;k < chain[i]->nb_state;k++) {
        os << setw(width) << chain[i]->transition[j][k];
      }
      os << endl;

      // mise a jour des indices des etats

      for (k = 0;k < i;k++) {
        if (state_index[k] < chain[i]->nb_state - 1) {
          state_index[k]++;
          break;
        }
        else {
          state_index[k] = 0;
        }
      }
    }
  }

  // calcul des largeurs des colonnes

  marginal_dist->cumul_computation();
  width = column_width(marginal_dist->nb_value , marginal_dist->cumul) + ASCII_SPACE;

  os << "\n\n" << SEQ_label[SEQL_LIKELIHOODS] << endl;

  if (begin) {
    for (i = 0;i < max_order;i++) {
      os << "  ";
    }
    os << "  ";

    for (i = 0;i < chain[1]->nb_state;i++) {
      os << setw(width) << i;
    }
    os << endl;

    for (i = 0;i < max_order;i++) {
      os << "  ";
    }
    os << "  ";

    initial_likelihood = 0.;
    for (i = 0;i < chain[1]->nb_state;i++) {
      if (pchain_data[1]->initial[i] > 0) {
        initial_likelihood += pchain_data[1]->initial[i] * log(chain[1]->initial[i]);
        os << setw(width) << pchain_data[1]->initial[i] * log(chain[1]->initial[i]);
      }
      else {
        os << setw(width) << 0.;
      }
    }
    os << "  " << setw(width) << initial_likelihood << "\n\n";
  }

  for (i = 0;i < max_order;i++) {
    os << "  ";
  }
  os << "  ";

  for (i = 0;i < chain[1]->nb_state;i++) {
    os << setw(width) << i;
  }
  os << setw(width) << " " << STAT_label[STATL_FREE_PARAMETERS]
     << "    " << SEQ_label[SEQL_DELTA] << " " << STAT_label[STATL_FREE_PARAMETERS]
     << "    " << SEQ_label[SEQL_DELTA] << " " << STAT_criterion_word[BIC]
     << "    " << SEQ_label[SEQL_DELTA] << " " << STAT_criterion_word[BICc] << endl;

  for (i = 0;i < max_order;i++) {
    os << "  ";
  }
  os << "  ";

  memory_likelihood[0][0] = 0.;
  nb_parameter[0][0] = 0;
  for (i = 0;i < marginal_dist->nb_value;i++) {
    if (marginal_histo->frequency[i] > 0) {
      nb_parameter[0][0]++;
      memory_likelihood[0][0] += marginal_dist->mass[i];
    }
    os << setw(width) << marginal_dist->mass[i];
  }
  if (nb_parameter[0][0] > 0) {
    nb_parameter[0][0]--;
  }
  os << "  " << setw(width) << memory_likelihood[0][0]
     << "  " << setw(width) << nb_parameter[0][0] << endl;

  for (i = 1;i <= max_order;i++) {
    os << "\n";
    for (j = 0;j < i;j++) {
      state_index[j] = 0;
    }

    m = 0;
    child_likelihood = 0.;
    child_nb_parameter = 0;

    for (j = 0;j < chain[i]->nb_row;j++) {
      for (k = i;k <max_order;k++) {
        os << "  ";
      }
      for (k = 0;k < i;k++) {
        os << state_index[k] << " ";
      }
      os << "  ";

      memory_likelihood[i][j] = 0.;
      nb_parameter[i][j] = 0;
      for (k = 0;k < chain[i]->nb_state;k++) {
        if (pchain_data[i]->transition[j][k] > 0) {
          nb_parameter[i][j]++;
          memory_likelihood[i][j] += pchain_data[i]->transition[j][k] * log(chain[i]->transition[j][k]);
          os << setw(width) << pchain_data[i]->transition[j][k] * log(chain[i]->transition[j][k]);
        }
        else {
          os << setw(width) << 0.;
        }
      }
      if (nb_parameter[i][j] > 0) {
        nb_parameter[i][j]--;
      }

      child_likelihood += memory_likelihood[i][j];
      child_nb_parameter += nb_parameter[i][j];

      os << "  " << setw(width) << memory_likelihood[i][j]
         << "  " << setw(width) << nb_parameter[i][j];

      // mise a jour des indices des etats

      for (k = 0;k < i;k++) {
        if (state_index[k] < chain[i]->nb_state - 1) {
          state_index[k]++;
          break;
        }
        else {
          state_index[k] = 0;

          if (k == 0) {
            os << "  " << setw(width) << child_nb_parameter - nb_parameter[i - 1][m]
               << "  " << setw(width) << 2 * (child_likelihood - memory_likelihood[i - 1][m]) -
               (child_nb_parameter - nb_parameter[i - 1][m]) * log((double)marginal_histo->nb_element);

            if (memory_count[i - 1][m] > 0) {
              os << "  " << setw(width) << 2 * (child_likelihood - memory_likelihood[i - 1][m]) -
                 (child_nb_parameter - nb_parameter[i - 1][m]) * log((double)memory_count[i - 1][m]);
            }
            else {
              os << "  " << setw(width) << 0;
            }

            m++;
            child_likelihood = 0.;
            child_nb_parameter = 0;
          }
        }
      }

      os << endl;
    }
  }

  delete marginal_dist;

  for (i = 1;i <= max_order;i++) {
    delete chain[i];
  }
  delete [] chain;

  delete [] pchain_data;
  delete marginal_histo;

  for (i = 0;i <= max_order;i++) {
    delete [] memory_count[i];
  }
  delete [] memory_count;

  for (i = 0;i <= max_order;i++) {
    delete [] memory_likelihood[i];
  }
  delete [] memory_likelihood;

  for (i = 0;i <= max_order;i++) {
    delete [] nb_parameter[i];
  }
  delete [] nb_parameter;

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des comptage des transitions pour differents ordres dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path, ordre maximum,
 *              pointeurs sur les objets Chain_data,
 *              flag pour prendre en compte le debut de la sequence.
 *
 *--------------------------------------------------------------*/

bool Chain_data::transition_count_ascii_write(Format_error &error , const char *path , int max_order ,
                                              const Chain_data **chain_data , bool begin) const

{
  bool status;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    transition_count_ascii_write(out_file , max_order , chain_data , begin);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des comptage des transitions pour differents ordres
 *  dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path, ordre maximum,
 *              pointeurs sur les objets Chain_data,
 *              flag pour prendre en compte le debut de la sequence.
 *
 *--------------------------------------------------------------*/

bool Chain_data::transition_count_spreadsheet_write(Format_error &error , const char *path , int max_order ,
                                                    const Chain_data **chain_data , bool begin) const

{
  register int i , j , k , m;
  bool status;
  int initial_count , child_nb_parameter , *pfrequency , **memory_count , **nb_parameter ,
      state_index[ORDER];
  double initial_likelihood , child_likelihood , **memory_likelihood;
  Distribution *marginal_dist;
  Histogram *marginal_histo;
  Chain **chain;
  const Chain_data **pchain_data;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    pchain_data = new const Chain_data*[max_order + 1];
    pchain_data[1] = this;
    for (i = 2;i <= max_order;i++) {
      pchain_data[i] = chain_data[i - 2];
    }

    marginal_histo = new Histogram(pchain_data[1]->nb_state);

    memory_count = new int*[max_order + 1];
    memory_count[0] = new int[1];
    for (i = 1;i <= max_order;i++) {
      memory_count[i] = new int[pchain_data[i]->nb_row];
    }

    nb_parameter = new int*[max_order + 1];
    nb_parameter[0] = new int[1];
    for (i = 1;i <= max_order;i++) {
      nb_parameter[i] = new int[pchain_data[i]->nb_row];
    }

    memory_likelihood = new double*[max_order + 1];
    memory_likelihood[0] = new double[1];
    for (i = 1;i <= max_order;i++) {
      memory_likelihood[i] = new double[pchain_data[i]->nb_row];
    }

    pfrequency = marginal_histo->frequency;
    for (i = 0;i < pchain_data[1]->nb_state;i++) {
      switch (begin) {
      case false :
        *pfrequency = 0;
        break;
      case true :
        *pfrequency = pchain_data[1]->initial[i];
        break;
      }

      for (j = 0;j < pchain_data[1]->nb_state;j++) {
        *pfrequency += pchain_data[1]->transition[j][i];
      }
      pfrequency++;
    }

    marginal_histo->nb_element_computation();

    if (begin) {
      out_file << "\n" << SEQ_label[SEQL_INITIAL_COUNTS] << endl;

      for (i = 0;i < pchain_data[1]->nb_state;i++) {
        out_file << "\t" << i;
      }
      out_file << endl;

      initial_count = 0;
      for (i = 0;i < pchain_data[1]->nb_state;i++) {
        initial_count += pchain_data[1]->initial[i];
        out_file << "\t" << pchain_data[1]->initial[i];
      }
      out_file << "\t" << initial_count << endl;
    }

    out_file << "\n" << SEQ_label[SEQL_TRANSITION_COUNTS] << endl;

    for (i = 0;i < pchain_data[1]->nb_state;i++) {
      out_file << "\t" << i;
    }
    out_file << endl;

    pfrequency = marginal_histo->frequency;
    for (i = 0;i < marginal_histo->nb_value;i++) {
      out_file << "\t" << *pfrequency++;
    }
    out_file << "\t" << marginal_histo->nb_element << endl;

    memory_count[0][0] = marginal_histo->nb_element;

    for (i = 1;i <= max_order;i++) {
      out_file << "\n";
      for (j = 0;j < i;j++) {
        state_index[j] = 0;
      }

      for (j = 0;j < pchain_data[i]->nb_row;j++) {
        for (k = 0;k < i;k++) {
          out_file << state_index[k] << " ";
        }

        memory_count[i][j] = 0;
        for (k = 0;k < pchain_data[i]->nb_state;k++) {
          memory_count[i][j] += pchain_data[i]->transition[j][k];
          out_file << "\t" << pchain_data[i]->transition[j][k];
        }
        out_file << "\t" << memory_count[i][j] << endl;

        // mise a jour des indices des etats

        for (k = 0;k < i;k++) {
          if (state_index[k] < pchain_data[i]->nb_state - 1) {
            state_index[k]++;
            break;
          }
          else {
            state_index[k] = 0;
          }
        }
      }
    }
    out_file << endl;

    marginal_dist = new Distribution(*marginal_histo);

    chain = new Chain*[max_order + 1];
    for (i = 1;i <= max_order;i++) {
      chain[i] = new Chain('o' , pchain_data[i]->nb_state ,
                           (int)pow((double)pchain_data[i]->nb_state , i) , true);
      pchain_data[i]->estimation(*chain[i]);
    }

    if (begin) {
      out_file << "\n" << STAT_word[STATW_INITIAL_PROBABILITIES] << endl;

      for (i = 0;i < chain[1]->nb_state;i++) {
        out_file << "\t" << i;
      }
      out_file << endl;

      for (i = 0;i < chain[1]->nb_state;i++) {
        out_file << "\t" << chain[1]->initial[i];
      }
      out_file << endl;
    }

    out_file << "\n" << STAT_word[STATW_TRANSITION_PROBABILITIES] << endl;

    for (i = 0;i < chain[1]->nb_state;i++) {
      out_file << "\t" << i;
    }
    out_file << endl;

    for (i = 0;i < marginal_dist->nb_value;i++) {
      out_file << "\t" << marginal_dist->mass[i];
    }
    out_file << endl;

    for (i = 1;i <= max_order;i++) {
      out_file << "\n";
      for (j = 0;j < i;j++) {
        state_index[j] = 0;
      }

      for (j = 0;j < chain[i]->nb_row;j++) {
        for (k = 0;k < i;k++) {
          out_file << state_index[k] << " ";
        }

        for (k = 0;k < chain[i]->nb_state;k++) {
          out_file << "\t" << chain[i]->transition[j][k];
        }
        out_file << endl;

        // mise a jour des indices des etats

        for (k = 0;k < i;k++) {
          if (state_index[k] < chain[i]->nb_state - 1) {
            state_index[k]++;
            break;
          }
          else {
            state_index[k] = 0;
          }
        }
      }
    }

    out_file << "\n\n" << SEQ_label[SEQL_LIKELIHOODS] << endl;

    if (begin) {
      for (i = 0;i < chain[1]->nb_state;i++) {
        out_file << "\t" << i;
      }
      out_file << endl;

      initial_likelihood = 0.;
      for (i = 0;i < chain[1]->nb_state;i++) {
        if (pchain_data[1]->initial[i] > 0) {
          initial_likelihood += pchain_data[1]->initial[i] * log((double)chain[1]->initial[i]);
          out_file << "\t" << pchain_data[1]->initial[i] * log((double)chain[1]->initial[i]);
        }
        else {
          out_file << "\t" << 0.;
        }
      }
      out_file << "\t" << initial_likelihood << "\n\n";
    }

    for (i = 0;i < chain[1]->nb_state;i++) {
      out_file << "\t" << i;
    }
    out_file << "\t\t" << STAT_label[STATL_FREE_PARAMETERS]
             << "\t" << SEQ_label[SEQL_DELTA] << " " << STAT_label[STATL_FREE_PARAMETERS]
             << "\t" << SEQ_label[SEQL_DELTA] << " " << STAT_criterion_word[BIC]
             << "\t" << SEQ_label[SEQL_DELTA] << " " << STAT_criterion_word[BICc] << endl;

    memory_likelihood[0][0] = 0.;
    nb_parameter[0][0] = 0;
    for (i = 0;i < marginal_dist->nb_value;i++) {
      if (marginal_histo->frequency[i] > 0) {
        nb_parameter[0][0]++;
        memory_likelihood[0][0] += marginal_histo->frequency[i] * log((double)marginal_dist->mass[i]);
        out_file << "\t" << marginal_histo->frequency[i] * log((double)marginal_dist->mass[i]);
      }
      else {
        out_file << "\t" << 0;
      }
    }
    if (nb_parameter[0][0] > 0) {
      nb_parameter[0][0]--;
    }
    out_file << "\t" << memory_likelihood[0][0] << "\t" << nb_parameter[0][0] << endl;

    for (i = 1;i <= max_order;i++) {
      out_file << "\n";
      for (j = 0;j < i;j++) {
        state_index[j] = 0;
      }

      m = 0;
      child_likelihood = 0.;
      child_nb_parameter = 0;

      for (j = 0;j < chain[i]->nb_row;j++) {
        for (k = 0;k < i;k++) {
          out_file << state_index[k] << " ";
        }

        memory_likelihood[i][j] = 0.;
        nb_parameter[i][j] = 0;
        for (k = 0;k < chain[i]->nb_state;k++) {
          if (pchain_data[i]->transition[j][k] > 0) {
            nb_parameter[i][j]++;
            memory_likelihood[i][j] += pchain_data[i]->transition[j][k] * log((double)chain[i]->transition[j][k]);
            out_file << "\t" << pchain_data[i]->transition[j][k] * log((double)chain[i]->transition[j][k]);
          }
          else {
            out_file << "\t" << 0;
          }
        }
        if (nb_parameter[i][j] > 0) {
          nb_parameter[i][j]--;
        }

        child_likelihood += memory_likelihood[i][j];
        child_nb_parameter += nb_parameter[i][j];

        out_file << "\t" << memory_likelihood[i][j] << "\t" << nb_parameter[i][j];

        // mise a jour des indices des etats

        for (k = 0;k < i;k++) {
          if (state_index[k] < chain[i]->nb_state - 1) {
            state_index[k]++;
            break;
          }
          else {
            state_index[k] = 0;

            if (k == 0) {
              out_file << "\t" << child_nb_parameter - nb_parameter[i - 1][m]
                       << "\t" << 2 * (child_likelihood - memory_likelihood[i - 1][m]) -
                       (child_nb_parameter - nb_parameter[i - 1][m]) *
                       log((double)marginal_histo->nb_element);

              if (memory_count[i - 1][m] > 0) {
                out_file << "\t" << 2 * (child_likelihood - memory_likelihood[i - 1][m]) -
                         (child_nb_parameter - nb_parameter[i - 1][m]) *
                         log((double)memory_count[i - 1][m]);
              }
              else {
                out_file << "\t" << 0;
              }

              m++;
              child_likelihood = 0.;
              child_nb_parameter = 0;
            }
          }
        }

        out_file << endl;
      }
    }

    delete marginal_dist;

    for (i = 1;i <= max_order;i++) {
      delete chain[i];
    }
    delete [] chain;

    delete [] pchain_data;
    delete marginal_histo;

    for (i = 0;i <= max_order;i++) {
      delete [] memory_count[i];
    }
    delete [] memory_count;

    for (i = 0;i <= max_order;i++) {
      delete [] memory_likelihood[i];
    }
    delete [] memory_likelihood;

    for (i = 0;i <= max_order;i++) {
      delete [] nb_parameter[i];
    }
    delete [] nb_parameter;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Comptage des transitions pour differents ordres.
 *
 *  arguments : indice de la variable, ordre maximum,
 *              flag pour prendre en compte le debut de la sequence,
 *              path, format de fichier ('a' : ASCII, 's' : Spreadsheet).
 *
 *--------------------------------------------------------------*/

bool Markovian_sequences::transition_count_0(Format_error &error , ostream &os , int max_order ,
                                             bool begin , const char *path , char format) const

{
  bool status = true;
  register int i;
  const Chain_data **chain_data;


  error.init();

  if (nb_variable > 1) {
    status = false;
    error.correction_update(STAT_error[STATR_NB_VARIABLE] , 1);
  }
  if ((max_order < 1) || (max_order > ORDER)) {
    status = false;
    error.update(SEQ_error[SEQR_ORDER]);
  }

  if (status) {
    chain_data = new const Chain_data*[max_order + 1];

    switch (begin) {

    case false : {
      chain_data[max_order] = new Chain_data('o' , marginal[0]->nb_value ,
                                             (int)pow((double)marginal[0]->nb_value , max_order));
      transition_count_computation(*chain_data[max_order] , max_order , begin);
      for (i = max_order - 1;i >= 1;i--) {
        chain_data[i] = new Chain_data(*chain_data[i + 1] , i + 1 , i);
      }
      break;
    }

    case true : {
      for (i = 1;i <= max_order;i++) {
        chain_data[i] = new Chain_data('o' , marginal[0]->nb_value ,
                                       (int)pow((double)marginal[0]->nb_value , i));
        transition_count_computation(*chain_data[i] , i , begin);
      }
      break;
    }
    }

#   ifdef MESSAGE
    chain_data[1]->transition_count_ascii_write(os , max_order , chain_data + 2 , begin);
#   endif

    if (path) {
      switch (format) {
      case 'a' :
        status = chain_data[1]->transition_count_ascii_write(error , path , max_order , chain_data + 2 , begin);
        break;
      case 's' :
        status = chain_data[1]->transition_count_spreadsheet_write(error , path , max_order , chain_data + 2 , begin);
        break;
      }
    }

    for (i = 1;i <= max_order;i++) {
      delete chain_data[i];
    }
    delete [] chain_data;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de differentes chaines de Markov pour un ensemble de sequences.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre de chaines de Markov,
 *              pointeur sur les chaines de Markov, path.
 *
 *--------------------------------------------------------------*/

bool Markovian_sequences::comparison(Format_error &error , ostream &os , int nb_model ,
                                     const Markov **imarkov , const char *path) const

{
  bool status = true;
  register int i , j;
  double **likelihood;


  error.init();

  if (!characteristics[0]) {
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

    if (test_hidden(1)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << 2 << ": "
                    << SEQ_error[SEQR_OVERLAP];
      error.update((error_message.str()).c_str());
    }

    else if (!characteristics[1]) {
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

  for (i = 0;i < nb_model;i++) {
    if (imarkov[i]->nb_output_process + 1 != nb_variable) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_MARKOV_CHAIN] << " " << i + 1 << ": "
                    << SEQ_error[SEQR_NB_OUTPUT_PROCESS];
      error.update((error_message.str()).c_str());
    }

    else {
      if (imarkov[i]->process[0]->nb_value < marginal[0]->nb_value) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_MARKOV_CHAIN] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_NB_STATE];
        error.update((error_message.str()).c_str());
      }

      if (nb_variable == 2) {
        if (imarkov[i]->process[1]->nb_value < marginal[1]->nb_value) {
          status = false;
          ostringstream error_message;
          error_message << SEQ_label[SEQL_MARKOV_CHAIN] << " " << i + 1 << ": "
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
        likelihood[i][j] = imarkov[j]->likelihood_computation(*this , i);
      }
    }

#   ifdef MESSAGE
    likelihood_write(os , nb_model , likelihood , SEQ_label[SEQL_MARKOV_CHAIN] , true);
#   endif

    if (path) {
      status = likelihood_write(error , path , nb_model , likelihood , SEQ_label[SEQL_MARKOV_CHAIN]);
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
 *  Simulation par une chaine de Markov.
 *
 *  arguments : reference sur un objet Format_error,
 *              histogramme des longueurs des sequences,
 *              flag sur le calcul des lois de comptage,
 *              flag calcul d'une divergence de Kullback-Leibler.
 *
 *--------------------------------------------------------------*/

Markov_data* Markov::simulation(Format_error &error , const Histogram &hlength ,
                                bool counting_flag , bool divergence_flag) const

{
  bool status = true;
  register int i , j , k;
  int cumul_length , *pstate , *mstate , **poutput , power[ORDER];
  double **pcumul;
  Markov *markov;
  Markov_data *seq;


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

    // initialisations

    seq = new Markov_data(nb_output_process + 1 , hlength , false);
    seq->type[0] = STATE;

    seq->markov = new Markov(*this , false);

    markov = seq->markov;
    markov->create_cumul();
    markov->cumul_computation();

    if (markov->nb_output_process > 0) {
      poutput = new int*[markov->nb_output_process];
    }

    i = 1;
    for (j = 0;j < markov->order;j++) {
      power[j] = i;
      i *= markov->nb_state;
    }

    for (i = 0;i < seq->nb_sequence;i++) {
      pstate = seq->sequence[i][0];

      *pstate = cumul_method(markov->nb_state , markov->cumul_initial);
      for (j = 0;j < markov->nb_output_process;j++) {
        poutput[j] = seq->sequence[i][j + 1];
        *poutput[j] = markov->process[j + 1]->observation[*pstate]->simulation();
      }

      for (j = 1;j < seq->length[i];j++) {
        pcumul = markov->cumul_transition;
        mstate = pstate;

        for (k = 0;k < MIN(j - 1 , markov->order);k++) {
          pcumul += *mstate-- * power[markov->order - 1 - k];
        }
        for (k = MIN(j - 1 , markov->order);k < markov->order;k++) {
          pcumul += *mstate * power[markov->order - 1 - k];
        }

        *++pstate = cumul_method(markov->nb_state , *pcumul);
        for (k = 0;k < markov->nb_output_process;k++) {
          *++poutput[k] = markov->process[k + 1]->observation[*pstate]->simulation();
        }
      }
    }

    markov->remove_cumul();

    if (markov->nb_output_process > 0) {
      delete [] poutput;
    }

    // extraction des caracteristiques des sequences simulees

    for (i = 0;i < seq->nb_variable;i++) {
      seq->max_value_computation(i);
      seq->build_marginal_histogram(i);
    }

    seq->build_transition_count(markov->order);
    seq->build_observation_histogram();
    seq->build_characteristic();

    if (!divergence_flag) {
      markov->characteristic_computation(*seq , counting_flag);

      // calcul de la vraisemblance

      seq->likelihood = markov->likelihood_computation(*seq);
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une chaine de Markov.
 *
 *  arguments : reference sur un objet Format_error,
 *              nombre et longueur des sequences,
 *              flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

Markov_data* Markov::simulation(Format_error &error , int nb_sequence ,
                                int length , bool counting_flag) const

{
  bool status = true;
  Markov_data *seq;


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
 *  Simulation par une chaine de Markov.
 *
 *  arguments : reference sur un objet Format_error, nombre de sequences,
 *              reference sur un objet Markovian_sequences,
 *              flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

Markov_data* Markov::simulation(Format_error &error , int nb_sequence ,
                                const Markovian_sequences &iseq , bool counting_flag) const

{
  Histogram *hlength;
  Markov_data *seq;


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
 *  Comparaison de chaines de Markov par calcul de divergences de Kullback-Leibler.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre de chaines
 *              de Markov, pointeur sur les chaines de Markov,
 *              histogrammes des longueurs des sequences, path.
 *
 *--------------------------------------------------------------*/

Distance_matrix* Markov::divergence_computation(Format_error &error , ostream &os ,
                                                int nb_model , const Markov **imarkov ,
                                                Histogram **hlength , const char *path) const

{
  bool status = true , lstatus;
  register int i , j , k;
  int cumul_length;
  double ref_likelihood , target_likelihood , **likelihood;
  const Markov **markov;
  Markovian_sequences *iseq , *seq;
  Markov_data *simul_seq;
  Distance_matrix *dist_matrix;
  ofstream *out_file;


  dist_matrix = 0;
  error.init();

  for (i = 0;i < nb_model - 1;i++) {
    if (imarkov[i]->nb_output_process == nb_output_process) {
      if (imarkov[i]->nb_state != nb_state) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_MARKOV_CHAIN] << " " << i + 2 << ": "
                      << SEQ_error[SEQR_NB_STATE];
        error.update((error_message.str()).c_str());
      }

      if (nb_output_process == 1) {
        if (imarkov[i]->process[1]->nb_value != process[1]->nb_value) {
          status = false;
          ostringstream error_message;
          error_message << SEQ_label[SEQL_MARKOV_CHAIN] << " " << i + 2 << ": "
                        << SEQ_error[SEQR_NB_OUTPUT];
          error.update((error_message.str()).c_str());
        }
      }
    }

    else if ((nb_output_process == 0) && (imarkov[i]->nb_output_process == 1)) {
      if (imarkov[i]->process[1]->nb_value != nb_state) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_MARKOV_CHAIN] << " " << i + 2 << ": "
                      << SEQ_error[SEQR_NB_OUTPUT];
        error.update((error_message.str()).c_str());
      }
    }

    else {  // if ((nb_output_process == 1) && (imarkov[i]->nb_output_process == 0))
      if (imarkov[i]->nb_state != process[1]->nb_value) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_MARKOV_CHAIN] << " " << i + 2 << ": "
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

    markov = new const Markov*[nb_model];

    markov[0] = this;
    for (i = 1;i < nb_model;i++) {
      markov[i] = imarkov[i - 1];
    }

    dist_matrix = new Distance_matrix(nb_model , SEQ_label[SEQL_MARKOV_CHAIN]);

    for (i = 0;i < nb_model;i++) {

      // simulation d'un echantillon de sequences a partir d'une chaine de Markov

      simul_seq = markov[i]->simulation(error , *hlength[i] , true);

      likelihood = new double*[simul_seq->nb_sequence];
      for (j = 0;j < simul_seq->nb_sequence;j++) {
        likelihood[j] = new double[nb_model];
      }

      ref_likelihood = 0.;
      for (j = 0;j < simul_seq->nb_sequence;j++) {
        likelihood[j][i] = markov[i]->likelihood_computation(*simul_seq , j);
        ref_likelihood += likelihood[j][i];
      }

      if (markov[i]->nb_output_process == 1) {
        iseq = simul_seq->remove_variable_1();
      }
      else {
        iseq = simul_seq;
      }

      // calcul des vraisemblances de l'echantillon pour chacune des chaines de Markov

      for (j = 0;j < nb_model;j++) {
        if (j != i) {
          if (markov[j]->nb_output_process == 1) {
            seq = iseq->transcode(error , markov[j]->process[1]);
          }
          else {
            seq = iseq;
          }

          target_likelihood = 0.;
          for (k = 0;k < seq->nb_sequence;k++) {
            likelihood[k][j] = markov[j]->likelihood_computation(*seq , k);
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

          if (markov[j]->nb_output_process == 1) {
            delete seq;
          }
        }
      }

#     ifdef MESSAGE
      os << SEQ_label[SEQL_MARKOV_CHAIN] << " " << i + 1 << ": " << simul_seq->nb_sequence << " "
         << SEQ_label[SEQL_SIMULATED] << " " << SEQ_label[simul_seq->nb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << endl;
      simul_seq->likelihood_write(os , nb_model , likelihood , SEQ_label[SEQL_MARKOV_CHAIN]);
#     endif

      if (out_file) {
        *out_file << SEQ_label[SEQL_MARKOV_CHAIN] << " " << i + 1 << ": " << simul_seq->nb_sequence << " "
                  << SEQ_label[SEQL_SIMULATED] << " " << SEQ_label[simul_seq->nb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << endl;
        simul_seq->likelihood_write(*out_file , nb_model , likelihood , SEQ_label[SEQL_MARKOV_CHAIN]);
      }

      for (j = 0;j < simul_seq->nb_sequence;j++) {
        delete [] likelihood[j];
      }
      delete [] likelihood;

      if (markov[i]->nb_output_process == 1) {
        delete iseq;
      }
      delete simul_seq;
    }

    if (out_file) {
      out_file->close();
      delete out_file;
    }

    delete markov;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de chaines de Markov par calcul de divergences de Kullback-Leibler.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre de chaines
 *              de Markov, pointeur sur les chaines de Markov,
 *              nombre et longueur des sequences, path.
 *
 *--------------------------------------------------------------*/

Distance_matrix* Markov::divergence_computation(Format_error &error , ostream &os ,
                                                int nb_model , const Markov **markov ,
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

    dist_matrix = divergence_computation(error , os , nb_model , markov , hlength , path);

    for (i = 0;i < nb_model;i++) {
      delete hlength[i];
    }
    delete [] hlength;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de chaines de Markov par calcul de divergences de Kullback-Leibler.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre de chaines
 *              de Markov, pointeur sur les chaines de Markov,
 *              pointeurs sur des objets Markovian_sequences, path.
 *
 *--------------------------------------------------------------*/

Distance_matrix* Markov::divergence_computation(Format_error &error , ostream &os ,
                                                int nb_model , const Markov **markov ,
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

    dist_matrix = divergence_computation(error , os , nb_model , markov , hlength , path);

    for (i = 0;i < nb_model;i++) {
      delete hlength[i];
    }
    delete [] hlength;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Markov_iterator.
 *
 *  argument : pointeur sur un objet Markov.
 *
 *--------------------------------------------------------------*/

Markov_iterator::Markov_iterator(Markov *imarkov)

{
  markov = imarkov;
  (markov->nb_iterator)++;

  if ((!(markov->cumul_initial)) || (!(markov->cumul_transition))) {
    markov->create_cumul();
    markov->cumul_computation();
  }

  chain = new Chain(*markov);

  state = new int[markov->order];
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Markov_iterator.
 *
 *  argument : reference sur un objet Markov_iterator.
 *
 *--------------------------------------------------------------*/

void Markov_iterator::copy(const Markov_iterator &it)

{
  register int i;


  markov = it.markov;
  (markov->nb_iterator)++;

  chain = new Chain(*(it.chain));

  state = new int[markov->order];
  for (i = 0;i < markov->order;i++) {
    state[i] = it.state[i];
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Markov_iterator.
 *
 *--------------------------------------------------------------*/

Markov_iterator::~Markov_iterator()

{
  (markov->nb_iterator)--;

  delete chain;

  delete [] state;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Markov_iterator.
 *
 *  argument : reference sur un objet Markov_iterator.
 *
 *--------------------------------------------------------------*/

Markov_iterator& Markov_iterator::operator=(const Markov_iterator &it)

{
  if (&it != this) {
    (markov->nb_iterator)--;

    delete chain;

    delete [] state;

    copy(it);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une chaine de Markov.
 *
 *  arguments : sequence, longueur de la sequence, flag initialisation.
 *
 *--------------------------------------------------------------*/

void Markov_iterator::simulation(int **seq , int length , bool initialization)

{
  register int i , j;
  int offset = 0 , *pstate , **poutput , power[ORDER];
  double **cumul;


  if (markov->nb_output_process > 0) {
    poutput = new int*[markov->nb_output_process];
  }

  pstate = seq[0];
  for (i = 0;i < markov->nb_output_process;i++) {
    poutput[i] = seq[i + 1];
  }

  i = 1;
  for (j = 0;j < markov->order;j++) {
    power[j] = i;
    i *= markov->nb_state;
  }

  if (initialization) {
    *pstate = cumul_method(markov->nb_state , markov->cumul_initial);

    state[0] = *pstate++;
    for (i = 1;i < markov->order;i++) {
      state[i] = state[0];
    }

    for (i = 0;i < markov->nb_output_process;i++) {
      *poutput[i]++ = markov->process[i + 1]->observation[state[0]]->simulation();
    }

    offset++;
  }

  for (i = offset;i < length;i++) {
    cumul = chain->cumul_transition;
    for (j = 0;j < markov->order;j++) {
      cumul += state[j] * power[markov->order - 1 - j];
    }

    *pstate = cumul_method(markov->nb_state , *cumul);

    for (j = markov->order - 1;j > 0;j--) {
      state[j] = state[j - 1];
    }
    state[0] = *pstate++;

    for (j = 0;j < markov->nb_output_process;j++) {
      *poutput[j]++ = markov->process[j + 1]->observation[state[0]]->simulation();
    }
  }

  if (markov->nb_output_process > 0) {
    delete [] poutput;
  }
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une chaine de Markov.
 *
 *  arguments : longueur de la sequence, flag initialisation.
 *
 *--------------------------------------------------------------*/

int** Markov_iterator::simulation(int length , bool initialization)

{
  register int i;
  int **seq;


  seq = new int*[markov->nb_output_process + 1];
  for (i = 0;i <= markov->nb_output_process;i++) {
    seq[i] = new int[length];
  }

  simulation(seq , length , initialization);

  return seq;
}
