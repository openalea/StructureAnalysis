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
 *       $Id: categorical_sequence_process2.cpp 18044 2015-04-23 09:32:23Z guedon $
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



#include <sstream>
#include <iomanip>

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

#include "sequences.h"
#include "sequence_label.h"
#include "tool/config.h"

#include <math.h>

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet CategoricalSequenceProcess.
 *
 *  arguments : stream, indice du processus d'observation,
 *              pointeurs sur les lois d'observation empiriques,
 *              la loi marginale empiriques et
 *              les caracteristiques des sequences observees,
 *              flag niveau de detail, flag fichier ,
 *              pointeurs sur les lois de l'intervalle de temps residuel.
 *
 *--------------------------------------------------------------*/

ostream& CategoricalSequenceProcess::ascii_print(ostream &os , int process ,
                                                 FrequencyDistribution **empirical_observation ,
                                                 FrequencyDistribution *marginal_distribution ,
                                                 const SequenceCharacteristics *characteristics ,
                                                 bool exhaustive , bool file_flag ,
                                                 Forward **forward) const

{
  register int i , j;
  int buff , width[2];
  long old_adjust;
  double *pmass , scale[NB_STATE];
  const Distribution *pobservation[NB_STATE];
  double sum_probabilities, current_probability;

  old_adjust = os.setf(ios::left , ios::adjustfield);

  if (observation) {
    for (i = 0;i < nb_state;i++) {
      os << "\n" << STAT_word[STATW_STATE] << " " << i << " "
         << STAT_word[STATW_OBSERVATION_DISTRIBUTION] << endl;
      pmass = observation[i]->mass + observation[i]->offset;
      sum_probabilities = 0.0;
      current_probability = 0.0;
      for (j = observation[i]->offset;j < observation[i]->nb_value - 1;j++) {
        current_probability = (int)*pmass + floor((*pmass - (int)*pmass)*10000)/10000;
        sum_probabilities += current_probability;
        if (*pmass > 0.) {
          os << STAT_word[STATW_OUTPUT] << " " << j << " : " << current_probability << endl;
        }
        pmass++;
      }
      j = observation[i]->nb_value-1;
      os << STAT_word[STATW_OUTPUT] << " " << j << " : " << 1 - sum_probabilities << endl;

      if ((empirical_observation) && (empirical_observation[i]->nb_element > 0) && (exhaustive)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "   | " << STAT_label[STATL_STATE] << " " << i << " "
           << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
           << " | " << STAT_label[STATL_STATE] << " " << i << " "
           << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION] << endl;

        observation[i]->ascii_print(os , file_flag , false , false , empirical_observation[i]);
      }
    }

    // calcul des largeurs des colonnes

    width[0] = column_width(nb_state - 1) + ASCII_SPACE;

    width[1] = 0;
    for (i = 0;i < nb_state;i++) {
      buff = column_width(observation[i]->nb_value , observation[i]->mass);
      if (buff > width[1]) {
        width[1] = buff;
      }
    }
    width[1] += ASCII_SPACE;

    // ecriture de la matrice des probabilites d'observation

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_OBSERVATION_PROBABILITIY_MATRIX] << endl;

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << setw(width[0]) << " ";
    for (i = 0;i < nb_value;i++) {
      os << setw(width[1]) << i;
    }

    for (i = 0;i < nb_state;i++) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << setw(width[0]) << i;
      for (j = 0;j < nb_value;j++) {
        os << setw(width[1]) << observation[i]->mass[j];
      }
    }
    os << endl;

    if (marginal_distribution) {
      double likelihood , information;
      Test test(CHI2);


      if ((weight) && (mixture)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_THEORETICAL] << " " << STAT_label[STATL_WEIGHTS] << ":";

        for (i = 0;i < nb_state;i++) {
          os << " " << weight->mass[i];
        }
        os << endl;

        likelihood = mixture->likelihood_computation(*marginal_distribution);
        information = marginal_distribution->information_computation();

        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << likelihood / marginal_distribution->nb_element << ")" << endl;

        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_MAX_LIKELIHOOD] << ": " << information << "   ("
           << STAT_label[STATL_INFORMATION] << ": " << information / marginal_distribution->nb_element << ")" << endl;

        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_DEVIANCE] << ": " << 2 * (information - likelihood) << endl;

        mixture->chi2_fit(*marginal_distribution , test);
        os << "\n";
        test.ascii_print(os , file_flag);

        if (exhaustive) {
          for (i = 0;i < nb_state;i++) {
            pobservation[i] = observation[i];
            scale[i] = weight->mass[i] * marginal_distribution->nb_element;
          }

          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << "   | " << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
          for (i = 0;i < nb_state;i++) {
            os << " | " << STAT_label[STATL_STATE] << " " << i << " " << STAT_label[STATL_OBSERVATION]
               << " " << STAT_label[STATL_DISTRIBUTION];
          }
          os << " | " << STAT_label[STATL_MIXTURE] << " | " << STAT_label[STATL_CUMULATIVE]
             << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION]
             << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE]
             << " " << STAT_label[STATL_FUNCTION] << endl;

          mixture->ascii_print(os , nb_state , pobservation , scale ,
                               file_flag , true , marginal_distribution);
        }
      }

      if ((restoration_weight) && (restoration_mixture)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_RESTORATION] << " " << STAT_label[STATL_WEIGHTS] << ":";

        for (i = 0;i < nb_state;i++) {
          os << " " << restoration_weight->mass[i];
        }
        os << endl;

        likelihood = restoration_mixture->likelihood_computation(*marginal_distribution);
        information = marginal_distribution->information_computation();

        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << likelihood / marginal_distribution->nb_element << ")" << endl;

        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_MAX_LIKELIHOOD] << ": " << information << "   ("
           << STAT_label[STATL_INFORMATION] << ": " << information / marginal_distribution->nb_element << ")" << endl;

        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_DEVIANCE] << ": " << 2 * (information - likelihood) << endl;

        restoration_mixture->chi2_fit(*marginal_distribution , test);
        os << "\n";
        test.ascii_print(os , file_flag);

        if (exhaustive) {
          for (i = 0;i < nb_state;i++) {
            pobservation[i] = observation[i];
            scale[i] = restoration_weight->mass[i] * marginal_distribution->nb_element;
          }

          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << "   | " << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
          for (i = 0;i < nb_state;i++) {
            os << " | " << STAT_label[STATL_STATE] << " " << i << " " << STAT_label[STATL_OBSERVATION]
               << " " << STAT_label[STATL_DISTRIBUTION];
          }
          os << " | " << STAT_label[STATL_MIXTURE] << " | " << STAT_label[STATL_CUMULATIVE]
             << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION]
             << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE]
             << " " << STAT_label[STATL_FUNCTION] << endl;

          restoration_mixture->ascii_print(os , nb_state , pobservation , scale ,
                                           file_flag , true , marginal_distribution);
        }
      }
    }
  }

  if (((index_value) || (characteristics)) && (exhaustive)) {
    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "  ";

    for (i = 0;i < nb_value;i++) {
      if ((characteristics) && (i < characteristics->nb_value)) {
        os << " | " << SEQ_label[SEQL_OBSERVED] << " "
           << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i;
      }
      if (index_value) {
        os << " | " << SEQ_label[SEQL_THEORETICAL] << " "
           << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i;
      }
    }
    if (characteristics) {
      os << " | " << STAT_label[STATL_FREQUENCY];
    }
    os << endl;

    if (index_value) {
      index_value->ascii_print(os , file_flag ,
                               (characteristics ? characteristics->index_value : NULL));
    }
    else {
      characteristics->index_value->ascii_print(os , file_flag);
    }
  }

  if ((first_occurrence) || (characteristics)) {
    for (i = 0;i < nb_value;i++) {
      if (first_occurrence) {
        if (no_occurrence[i] > 0.) {
          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << SEQ_label[SEQL_NO_OCCURRENCE] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << ": " << no_occurrence[i] << endl;
        }

        if (first_occurrence[i]) {
          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << SEQ_label[SEQL_FIRST_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << " " << STAT_label[STATL_DISTRIBUTION] << endl;
          first_occurrence[i]->ascii_characteristic_print(os , false , file_flag);
        }
      }

      if ((characteristics) && (i < characteristics->nb_value)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_FIRST_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
           << " " << i << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
        characteristics->first_occurrence[i]->ascii_characteristic_print(os , false , file_flag);
      }

      if ((((first_occurrence) && (first_occurrence[i])) ||
           ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->first_occurrence[i]->nb_element > 0))) && (exhaustive)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "  ";
        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->first_occurrence[i]->nb_element > 0)) {
          os << " | " << SEQ_label[SEQL_FIRST_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        }

        if ((first_occurrence) && (first_occurrence[i])) {
          os << " | " << SEQ_label[SEQL_FIRST_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << " " << STAT_label[STATL_DISTRIBUTION];
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->first_occurrence[i]->nb_element > 0)) {
            os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
               << STAT_label[STATL_FUNCTION];
          }
          os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION] << endl;

          first_occurrence[i]->ascii_print(os , file_flag , true , true ,
                                           (((characteristics) && (i < characteristics->nb_value) && (characteristics->first_occurrence[i]->nb_element > 0)) ? characteristics->first_occurrence[i] : NULL));
        }

        else {
          os << endl;
          characteristics->first_occurrence[i]->ascii_print(os , file_flag);
        }
      }
    }
  }

  if ((recurrence_time) || (characteristics)) {
    for (i = 0;i < nb_value;i++) {
      if (recurrence_time) {
        if (leave[i] > 0.) {
          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << SEQ_label[SEQL_LEAVING] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << ": " << leave[i] << endl;
        }

        if (recurrence_time[i]) {
          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
          recurrence_time[i]->ascii_characteristic_print(os , false , file_flag);
        }
      }

      if ((characteristics) && (i < characteristics->nb_value)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
           << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
        characteristics->recurrence_time[i]->ascii_characteristic_print(os , false , file_flag);
      }

      if ((((recurrence_time) && (recurrence_time[i])) ||
           ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->recurrence_time[i]->nb_element > 0))) && (exhaustive)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "  ";
        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->recurrence_time[i]->nb_element > 0)) {
          os << " | " << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        }

        if ((recurrence_time) && (recurrence_time[i])) {
          os << " | " << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_DISTRIBUTION];
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->recurrence_time[i]->nb_element > 0)) {
            os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
               << STAT_label[STATL_FUNCTION];
          }
          os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION] << endl;

          recurrence_time[i]->ascii_print(os , file_flag , true , true ,
                                          (((characteristics) && (i < characteristics->nb_value) && (characteristics->recurrence_time[i]->nb_element > 0)) ?
                                           characteristics->recurrence_time[i] : NULL));
        }

        else {
          os << endl;
          characteristics->recurrence_time[i]->ascii_print(os , file_flag);
        }
      }
    }
  }

  if ((sojourn_time) || (characteristics)) {
    for (i = 0;i < nb_value;i++) {
      if (sojourn_time) {
        if (absorption[i] > 0.) {
          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << SEQ_label[SEQL_ABSORPTION] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << ": " << absorption[i] << endl;
        }

        if (sojourn_time[i]) {
          if (sojourn_time[i]->ident != CATEGORICAL) {
            os << "\n" << STAT_word[STATW_STATE] << " " << i << " "
               << SEQ_word[SEQW_OCCUPANCY_DISTRIBUTION] << endl;
            sojourn_time[i]->ascii_print(os);
          }

          else {
            os << "\n";
            if (file_flag) {
              os << "# ";
            }
            os << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
               << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
          }

#         ifdef DEBUG
          sojourn_time[i]->ascii_characteristic_print(os , (sojourn_time[i]->ident == CATEGORICAL ? false : true) , file_flag);
#         endif

//          sojourn_time[i]->ascii_parametric_characteristic_print(os , true , file_flag);
          sojourn_time[i]->ascii_parametric_characteristic_print(os , (sojourn_time[i]->ident == CATEGORICAL ? false : true) , file_flag);

#         ifdef MESSAGE
          if (file_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_VARIATION_COEFF] << ": " << (sojourn_time[i]->ident == CATEGORICAL ? sqrt(sojourn_time[i]->variance) / sojourn_time[i]->mean : sqrt(sojourn_time[i]->parametric_variance_computation()) / sojourn_time[i]->parametric_mean_computation()) << endl;
#         endif

        }
      }

      if ((characteristics) && (i < characteristics->nb_value)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
           << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
        characteristics->sojourn_time[i]->ascii_characteristic_print(os , false , file_flag);
      }

      if ((((sojourn_time) && (sojourn_time[i])) ||
           ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->sojourn_time[i]->nb_element > 0))) && (exhaustive)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "  ";
        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->sojourn_time[i]->nb_element > 0)) {
          os << " | " << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        }

        if ((sojourn_time) && (sojourn_time[i])) {
          os << " | " << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION];
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->sojourn_time[i]->nb_element > 0)) {
            os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
               << STAT_label[STATL_FUNCTION];
          }
          os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION] << endl;

          sojourn_time[i]->Distribution::ascii_print(os , file_flag , true , false ,
                                                     (((characteristics) && (i < characteristics->nb_value) && (characteristics->sojourn_time[i]->nb_element > 0)) ?
                                                      characteristics->sojourn_time[i] : NULL));
        }

        else {
          os << endl;
          characteristics->sojourn_time[i]->ascii_print(os , file_flag);
        }
      }

      if ((forward) && (forward[i])) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_STATE] << " " << i << " " << SEQ_label[SEQL_FORWARD] << " "
           << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
        forward[i]->ascii_characteristic_print(os , false , file_flag);
      }

      if ((characteristics) && (i < characteristics->nb_value) &&
          (characteristics->initial_run)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_INITIAL_RUN] << " - "
           << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
           << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
        characteristics->initial_run[i]->ascii_characteristic_print(os , false , file_flag);

        if ((characteristics->initial_run[i]->nb_element > 0) && (exhaustive)) {
          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << "   | " << SEQ_label[SEQL_INITIAL_RUN] << " - "
             << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];

          if ((forward) && (forward[i])) {
            os << " | " << STAT_label[STATL_STATE] << " " << i << " " << SEQ_label[SEQL_FORWARD] << " "
               << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION]
               << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
               << STAT_label[STATL_FUNCTION] << " | " << STAT_label[STATL_CUMULATIVE]
               << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

            forward[i]->Distribution::ascii_print(os , file_flag , true , false ,
                                                  characteristics->initial_run[i]);
          }

          else {
            os << endl;
            characteristics->initial_run[i]->ascii_print(os , file_flag);
          }
        }
      }

      if ((characteristics) && (i < characteristics->nb_value)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_FINAL_RUN] << " - "
           << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
           << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
        characteristics->final_run[i]->ascii_characteristic_print(os , false , file_flag);
      }

      if ((((forward) && (forward[i])) ||
           ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->final_run[i]->nb_element > 0))) && (exhaustive)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "  ";
        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->final_run[i]->nb_element > 0)) {
          os << " | " << SEQ_label[SEQL_FINAL_RUN] << " - "
             << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        }

        if ((forward) && (forward[i])) {
          os << " | " << STAT_label[STATL_STATE] << " " << i << " " << SEQ_label[SEQL_FORWARD] << " "
             << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION];
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->final_run[i]->nb_element > 0)) {
            os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
               << STAT_label[STATL_FUNCTION];
          }
          os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION] << endl;

          forward[i]->Distribution::ascii_print(os , file_flag , true , false ,
                                                (((characteristics) && (i < characteristics->nb_value) && (characteristics->final_run[i]->nb_element > 0)) ?
                                                 characteristics->final_run[i] : NULL));
        }

        else {
          os << endl;
          characteristics->final_run[i]->ascii_print(os , file_flag);
        }
      }
    }
  }

  if ((nb_run) || ((characteristics) && (characteristics->nb_run))) {
    for (i = 0;i < nb_value;i++) {
      if (nb_run) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        if (length->variance == 0.) {
          os << SEQ_label[SEQL_NB_RUN_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
             << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
        }
        else {
          os << SEQ_label[SEQL_MIXTURE_OF] << SEQ_label[SEQL_NB_RUN_OF]
             << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTIONS] << endl;
        }
        nb_run[i]->ascii_characteristic_print(os , (length->variance > 0. ? false : true) , file_flag);
      }

      if ((characteristics) && (i < characteristics->nb_value) &&
          (characteristics->nb_run)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_NB_RUN_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
           << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
        characteristics->nb_run[i]->ascii_characteristic_print(os , (length->variance > 0. ? false : true) , file_flag);
      }

      if (exhaustive) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "  ";
        if ((characteristics) && (i < characteristics->nb_value) && (characteristics->nb_run) &&
            (characteristics->nb_run[i]->nb_element > 0)) {
          os << " | " << SEQ_label[SEQL_NB_RUN_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        }

        if (nb_run) {
          if (length->variance == 0.) {
            os << " | " << SEQ_label[SEQL_NB_RUN_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
               << " " << i << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
               << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION];
          }
          else {
            os << " | " << SEQ_label[SEQL_MIXTURE_OF] << SEQ_label[SEQL_NB_RUN_OF]
               << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
               << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTIONS];
          }
          if ((characteristics) && (i < characteristics->nb_value) && (characteristics->nb_run) &&
              (characteristics->nb_run[i]->nb_element > 0)) {
            os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
               << STAT_label[STATL_FUNCTION];
          }
          if (length->variance == 0.) {
            os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
               << STAT_label[STATL_FUNCTION] << endl;
          }
          else {
            os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE] << " "
               << STAT_label[STATL_FUNCTION] << endl;
          }

          nb_run[i]->ascii_print(os , file_flag , true , false ,
                                 (((characteristics) && (i < characteristics->nb_value) && (characteristics->nb_run) && (characteristics->nb_run[i]->nb_element > 0)) ? characteristics->nb_run[i] : NULL));
        }

        else {
          os << endl;
          characteristics->nb_run[i]->ascii_print(os , file_flag);
        }
      }
    }
  }

  if ((nb_occurrence) || ((characteristics) && (characteristics->nb_occurrence))) {
    for (i = 0;i < nb_value;i++) {
      if (nb_occurrence) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        if (length->variance == 0.) {
          os << SEQ_label[SEQL_NB_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
             << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
        }
        else {
          os << SEQ_label[SEQL_MIXTURE_OF] << SEQ_label[SEQL_NB_OCCURRENCE_OF]
             << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTIONS] << endl;
        }
        nb_occurrence[i]->ascii_characteristic_print(os , (length->variance > 0. ? false : true) , file_flag);
      }

      if ((characteristics) && (i < characteristics->nb_value) &&
          (characteristics->nb_occurrence)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_NB_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
           << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
        characteristics->nb_occurrence[i]->ascii_characteristic_print(os , (length->variance > 0. ? false : true) , file_flag);
      }

      if (exhaustive) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "  ";
        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->nb_occurrence) &&
            (characteristics->nb_occurrence[i]->nb_element > 0)) {
          os << " | " << SEQ_label[SEQL_NB_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        }

        if (nb_occurrence) {
          if (length->variance == 0.) {
            os << " | " << SEQ_label[SEQL_NB_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
               << " " << i << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
               << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION];
          }
          else {
            os << " | " << SEQ_label[SEQL_MIXTURE_OF] << SEQ_label[SEQL_NB_OCCURRENCE_OF]
               << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
               << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTIONS];
          }
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->nb_occurrence) &&
              (characteristics->nb_occurrence[i]->nb_element > 0)) {
            os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
               << STAT_label[STATL_FUNCTION];
          }
          if (length->variance == 0.) {
            os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
               << STAT_label[STATL_FUNCTION] << endl;
          }
          else {
            os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE] << " "
               << STAT_label[STATL_FUNCTION] << endl;
          }

          nb_occurrence[i]->ascii_print(os , file_flag , true , false ,
                                        (((characteristics) && (i < characteristics->nb_value) && (characteristics->nb_occurrence) && (characteristics->nb_occurrence[i]->nb_element > 0)) ? characteristics->nb_occurrence[i] : NULL));
        }

        else {
          os << endl;
          characteristics->nb_occurrence[i]->ascii_print(os , file_flag);
        }
      }
    }
  }

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet CategoricalSequenceProcess au format tableur.
 *
 *  arguments : stream, indice du processus d'observation,
 *              pointeurs sur les lois d'observation empiriques,
 *              la loi marginale empiriques empiriques,
 *              sur les caracteristiques des sequences observees et
 *              sur les lois de l'intervalle de temps residuel.
 *
 *--------------------------------------------------------------*/

ostream& CategoricalSequenceProcess::spreadsheet_print(ostream &os , int process ,
                                                       FrequencyDistribution **empirical_observation ,
                                                       FrequencyDistribution *marginal_distribution ,
                                                       const SequenceCharacteristics *characteristics ,
                                                       Forward **forward) const

{
  register int i , j;
  double *pmass , scale[NB_STATE];
  const Distribution *pobservation[NB_STATE];
  Curves *smoothed_curves;


  if (observation) {
    for (i = 0;i < nb_state;i++) {
      os << "\n" << STAT_word[STATW_STATE] << " " << i << "\t"
         << STAT_word[STATW_OBSERVATION_DISTRIBUTION] << endl;
      pmass = observation[i]->mass + observation[i]->offset;
      for (j = observation[i]->offset;j < observation[i]->nb_value;j++) {
        if (*pmass > 0.) {
          os << STAT_word[STATW_OUTPUT] << "\t" << j << "\t" << *pmass << endl;
        }
        pmass++;
      }

      if ((empirical_observation) && (empirical_observation[i]->nb_element > 0)) {
        os << "\n\t" << STAT_label[STATL_STATE] << " " << i << " "
           << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
           << "\t" << STAT_label[STATL_STATE] << " " << i << " "
           << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION] << endl;

        observation[i]->spreadsheet_print(os , false , false , false , empirical_observation[i]);
      }
    }

    if (marginal_distribution) {
      double likelihood , information;
      Test test(CHI2);


      if ((weight) && (mixture)) {
        os << "\n" << STAT_label[STATL_THEORETICAL] << " " << STAT_label[STATL_WEIGHTS];
        for (i = 0;i < nb_state;i++) {
          os << "\t" << weight->mass[i];
        }
        os << endl;

        likelihood = mixture->likelihood_computation(*marginal_distribution);
        information = marginal_distribution->information_computation();

        os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << likelihood << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / marginal_distribution->nb_element << endl;
        os << STAT_label[STATL_MAX_LIKELIHOOD] << "\t" << information << "\t"
           << STAT_label[STATL_INFORMATION] << "\t" << information / marginal_distribution->nb_element << endl;
        os << STAT_label[STATL_DEVIANCE] << "\t" << 2 * (information - likelihood) << endl;

        mixture->chi2_fit(*marginal_distribution , test);
        os << "\n";
        test.spreadsheet_print(os);

        for (i = 0;i < nb_state;i++) {
          pobservation[i] = observation[i];
          scale[i] = weight->mass[i] * marginal_distribution->nb_element;
        }

        os << "\n\t" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        for (i = 0;i < nb_state;i++) {
          os << "\t" << STAT_label[STATL_STATE] << " " << i << " " << STAT_label[STATL_OBSERVATION]
             << " " << STAT_label[STATL_DISTRIBUTION];
        }
        os << "\t" << STAT_label[STATL_MIXTURE] << "\t" << STAT_label[STATL_CUMULATIVE]
           << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION]
           << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE]
           << " " << STAT_label[STATL_FUNCTION] << endl;

        mixture->spreadsheet_print(os , nb_state , pobservation , scale , true ,
                                   marginal_distribution);
      }

      if ((restoration_weight) && (restoration_mixture)) {
        os << "\n" << STAT_label[STATL_RESTORATION] << " " << STAT_label[STATL_WEIGHTS];
        for (i = 0;i < nb_state;i++) {
          os << "\t" << restoration_weight->mass[i];
        }
        os << endl;

        likelihood = restoration_mixture->likelihood_computation(*marginal_distribution);
        information = marginal_distribution->information_computation();

        os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << likelihood << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / marginal_distribution->nb_element << endl;
        os << STAT_label[STATL_MAX_LIKELIHOOD] << "\t" << information << "\t"
           << STAT_label[STATL_INFORMATION] << "\t" << information / marginal_distribution->nb_element << endl;
        os << STAT_label[STATL_DEVIANCE] << "\t" << 2 * (information - likelihood) << endl;

        restoration_mixture->chi2_fit(*marginal_distribution , test);
        os << "\n";
        test.spreadsheet_print(os);

        for (i = 0;i < nb_state;i++) {
          pobservation[i] = observation[i];
          scale[i] = restoration_weight->mass[i] * marginal_distribution->nb_element;
        }

        os << "\n\t" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        for (i = 0;i < nb_state;i++) {
          os << "\t" << STAT_label[STATL_STATE] << " " << i << " " << STAT_label[STATL_OBSERVATION]
             << " " << STAT_label[STATL_DISTRIBUTION];
        }
        os << "\t" << STAT_label[STATL_MIXTURE] << "\t" << STAT_label[STATL_CUMULATIVE]
           << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION]
           << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE]
           << " " << STAT_label[STATL_FUNCTION] << endl;

        restoration_mixture->spreadsheet_print(os , nb_state , pobservation , scale , true ,
                                               marginal_distribution);
      }
    }
  }

  if ((index_value) || (characteristics)) {
    os << "\n";
    for (i = 0;i < nb_value;i++) {
      if ((characteristics) && (i < characteristics->nb_value)) {
        os << "\t" << SEQ_label[SEQL_OBSERVED] << " "
           << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i;
      }
      if (index_value) {
        os << "\t" << SEQ_label[SEQL_THEORETICAL] << " "
           << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i;
      }
    }
    if (characteristics) {
      os << "\t" << STAT_label[STATL_FREQUENCY];
    }
    os << endl;

    if (index_value) {
      index_value->spreadsheet_print(os , (characteristics ? characteristics->index_value : NULL));
    }
    else {
      characteristics->index_value->spreadsheet_print(os);
    }

    if (characteristics) {
      smoothed_curves = new Curves(*(characteristics->index_value) , 's');

      os << "\n" << SEQ_label[SEQL_SMOOTHED_OBSERVED_PROBABILITIES] << endl;
      for (i = 0;i < nb_value;i++) {
        if (i < characteristics->nb_value) {
          os << "\t" << SEQ_label[SEQL_OBSERVED] << " "
             << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i;
        }
        if (index_value) {
          os << "\t" << SEQ_label[SEQL_THEORETICAL] << " "
             << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i;
        }
      }
      os << "\t" << STAT_label[STATL_FREQUENCY] << endl;

      if (index_value) {
        index_value->spreadsheet_print(os , smoothed_curves);
      }
      else {
        smoothed_curves->spreadsheet_print(os);
      }

      delete smoothed_curves;
    }
  }

  if ((first_occurrence) || (characteristics)) {
    for (i = 0;i < nb_value;i++) {
      if (first_occurrence) {
        if (no_occurrence[i] > 0.) {
          os << "\n" << SEQ_label[SEQL_NO_OCCURRENCE] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << ": " << no_occurrence[i] << endl;
        }

        if (first_occurrence[i]) {
          os << "\n" << SEQ_label[SEQL_FIRST_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << " " << STAT_label[STATL_DISTRIBUTION] << endl;
          first_occurrence[i]->spreadsheet_characteristic_print(os);
        }
      }

      if ((characteristics) && (i < characteristics->nb_value)) {
        os << "\n" << SEQ_label[SEQL_FIRST_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
           << " " << i << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
        characteristics->first_occurrence[i]->spreadsheet_characteristic_print(os);
      }

      if (((first_occurrence) && (first_occurrence[i])) ||
          ((characteristics) && (i < characteristics->nb_value) &&
           (characteristics->first_occurrence[i]->nb_element > 0))) {
        os << "\n";
        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->first_occurrence[i]->nb_element > 0)) {
          os << "\t" << SEQ_label[SEQL_FIRST_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        }

        if ((first_occurrence) && (first_occurrence[i])) {
          os << "\t" << SEQ_label[SEQL_FIRST_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << " " << STAT_label[STATL_DISTRIBUTION];
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->first_occurrence[i]->nb_element > 0)) {
            os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
               << STAT_label[STATL_FUNCTION];
          }
          os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION] << endl;

          first_occurrence[i]->spreadsheet_print(os , true , false , true ,
                                                 (((characteristics) && (i < characteristics->nb_value) && (characteristics->first_occurrence[i]->nb_element > 0)) ? characteristics->first_occurrence[i] : NULL));
        }

        else {
          os << endl;
          characteristics->first_occurrence[i]->spreadsheet_print(os);
        }
      }
    }
  }

  if ((recurrence_time) || (characteristics)) {
    for (i = 0;i < nb_value;i++) {
      if (recurrence_time) {
        if (leave[i] > 0.) {
          os << "\n" << SEQ_label[SEQL_LEAVING] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << ": " << leave[i] << endl;
        }

        if (recurrence_time[i]) {
          os << "\n" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
          recurrence_time[i]->spreadsheet_characteristic_print(os);
        }
      }

      if ((characteristics) && (i < characteristics->nb_value)) {
        os << "\n" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
           << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
        characteristics->recurrence_time[i]->spreadsheet_characteristic_print(os);
      }

      if (((recurrence_time) && (recurrence_time[i])) ||
          ((characteristics) && (i < characteristics->nb_value) &&
           (characteristics->recurrence_time[i]->nb_element > 0))) {
        os << "\n";
        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->recurrence_time[i]->nb_element > 0)) {
          os << "\t" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        }

        if ((recurrence_time) && (recurrence_time[i])) {
          os << "\t" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_DISTRIBUTION];
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->recurrence_time[i]->nb_element > 0)) {
            os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
               << STAT_label[STATL_FUNCTION];
          }
          os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION] << endl;

          recurrence_time[i]->spreadsheet_print(os , true , false , true ,
                                                (((characteristics) && (i < characteristics->nb_value) && (characteristics->recurrence_time[i]->nb_element > 0)) ?
                                                 characteristics->recurrence_time[i] : NULL));
        }

        else {
          os << endl;
          characteristics->recurrence_time[i]->spreadsheet_print(os);
        }
      }
    }
  }

  if ((sojourn_time) || (characteristics)) {
    for (i = 0;i < nb_value;i++) {
      if (sojourn_time) {
        if (absorption[i] > 0.) {
          os << "\n" << SEQ_label[SEQL_ABSORPTION] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << ": " << absorption[i] << endl;
        }

        if (sojourn_time[i]) {
          if (sojourn_time[i]->ident != CATEGORICAL) {
            os << "\n" << STAT_word[STATW_STATE] << " " << i << "\t"
               << SEQ_word[SEQW_OCCUPANCY_DISTRIBUTION] << endl;
            sojourn_time[i]->spreadsheet_print(os);
          }
          else {
            os << "\n" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
               << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
          }

          sojourn_time[i]->spreadsheet_parametric_characteristic_print(os , (sojourn_time[i]->ident == CATEGORICAL ? false : true));
        }
      }

      if ((characteristics) && (i < characteristics->nb_value)) {
        os << "\n" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
           << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
        characteristics->sojourn_time[i]->spreadsheet_characteristic_print(os);
      }

      if (((sojourn_time) && (sojourn_time[i])) ||
          ((characteristics) && (i < characteristics->nb_value) &&
           (characteristics->sojourn_time[i]->nb_element > 0))) {
        os << "\n";
        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->sojourn_time[i]->nb_element > 0)) {
          os << "\t" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        }

        if ((sojourn_time) && (sojourn_time[i])) {
          os << "\t" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION];
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->sojourn_time[i]->nb_element > 0)) {
            os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
               << STAT_label[STATL_FUNCTION];
          }
          os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION] << endl;

          sojourn_time[i]->Distribution::spreadsheet_print(os , true , false , false ,
                                                           (((characteristics) && (i < characteristics->nb_value) && (characteristics->sojourn_time[i]->nb_element > 0)) ?
                                                            characteristics->sojourn_time[i] : NULL));
        }

        else {
          os << endl;
          characteristics->sojourn_time[i]->spreadsheet_print(os);
        }
      }

      if ((forward) && (forward[i])) {
        os << "\n" << STAT_label[STATL_STATE] << " " << i << " " << SEQ_label[SEQL_FORWARD] << " "
           << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
        forward[i]->spreadsheet_characteristic_print(os);
      }

      if ((characteristics) && (i < characteristics->nb_value) &&
          (characteristics->initial_run)) {
        os << "\n" << SEQ_label[SEQL_INITIAL_RUN] << " - "
           << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
           << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
        characteristics->initial_run[i]->spreadsheet_characteristic_print(os);

        if (characteristics->initial_run[i]->nb_element > 0) {
          os << "\n\t" << SEQ_label[SEQL_INITIAL_RUN] << " - "
             << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];

          if ((forward) && (forward[i])) {
            os << "\t" << STAT_label[STATL_STATE] << " " << i << " " << SEQ_label[SEQL_FORWARD] << " "
               << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION]
               << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
               << STAT_label[STATL_FUNCTION] << "\t" << STAT_label[STATL_CUMULATIVE] << " "
               << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

            forward[i]->Distribution::spreadsheet_print(os , true , false , false ,
                                                        characteristics->initial_run[i]);
          }

          else {
            os << endl;
            characteristics->initial_run[i]->spreadsheet_print(os);
          }
        }
      }

      if ((characteristics) && (i < characteristics->nb_value)) {
        os << "\n" << SEQ_label[SEQL_FINAL_RUN] << " - "
           << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
           << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
        characteristics->final_run[i]->spreadsheet_characteristic_print(os);
      }

      if (((forward) && (forward[i])) || ((characteristics) && (i < characteristics->nb_value) &&
           (characteristics->final_run[i]->nb_element > 0))) {
        os << "\n";
        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->final_run[i]->nb_element > 0)) {
          os << "\t" << SEQ_label[SEQL_FINAL_RUN] << " - "
             << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        }

        if ((forward) && (forward[i])) {
          os << "\t" << STAT_label[STATL_STATE] << " " << i << " " << SEQ_label[SEQL_FORWARD] << " "
             << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION];
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->final_run[i]->nb_element > 0)) {
            os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
               << STAT_label[STATL_FUNCTION];
          }
          os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION] << endl;

          forward[i]->Distribution::spreadsheet_print(os , true , false , false ,
                                                      (((characteristics) && (i < characteristics->nb_value) && (characteristics->final_run[i]->nb_element > 0)) ?
                                                       characteristics->final_run[i] : NULL));

        }

        else {
          os << endl;
          characteristics->final_run[i]->spreadsheet_print(os);
        }
      }
    }
  }

  if ((nb_run) || ((characteristics) && (characteristics->nb_run))) {
    for (i = 0;i < nb_value;i++) {
      if (nb_run) {
        if (length->variance == 0.) {
          os << "\n" << SEQ_label[SEQL_NB_RUN_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
             << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
        }
        else {
          os << "\n" << SEQ_label[SEQL_MIXTURE_OF] << SEQ_label[SEQL_NB_RUN_OF]
             << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTIONS] << endl;
        }
        nb_run[i]->spreadsheet_characteristic_print(os , (length->variance > 0. ? false : true));
      }

      if ((characteristics) && (i < characteristics->nb_value) && (characteristics->nb_run)) {
        os << "\n" << SEQ_label[SEQL_NB_RUN_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
           << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
        characteristics->nb_run[i]->spreadsheet_characteristic_print(os , (length->variance > 0. ? false : true));
      }

      os << "\n";
      if ((characteristics) && (i < characteristics->nb_value) && (characteristics->nb_run) &&
          (characteristics->nb_run[i]->nb_element > 0)) {
        os << "\t" << SEQ_label[SEQL_NB_RUN_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
           << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      }

      if (nb_run) {
        if (length->variance == 0.) {
          os << "\t" << SEQ_label[SEQL_NB_RUN_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
             << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION];
        }
        else {
          os << "\t" << SEQ_label[SEQL_MIXTURE_OF] << SEQ_label[SEQL_NB_RUN_OF]
             << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTIONS];
        }
        if ((characteristics) && (i < characteristics->nb_value) && (characteristics->nb_run) &&
            (characteristics->nb_run[i]->nb_element > 0)) {
          os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION];
        }
        if (length->variance == 0.) {
          os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION] << endl;
        }
        else {
          os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE] << " "
             << STAT_label[STATL_FUNCTION] << endl;
        }

        nb_run[i]->spreadsheet_print(os , true , false , false ,
                                     (((characteristics) && (i < characteristics->nb_value) && (characteristics->nb_run) && (characteristics->nb_run[i]->nb_element > 0)) ? characteristics->nb_run[i] : NULL));
      }

      else {
        os << endl;
        characteristics->nb_run[i]->spreadsheet_print(os);
      }
    }
  }

  if ((nb_occurrence) || ((characteristics) && (characteristics->nb_occurrence))) {
    for (i = 0;i < nb_value;i++) {
      if (nb_occurrence) {
        if (length->variance == 0.) {
          os << "\n" << SEQ_label[SEQL_NB_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
             << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
        }
        else {
          os << "\n" << SEQ_label[SEQL_MIXTURE_OF] << SEQ_label[SEQL_NB_OCCURRENCE_OF]
             << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTIONS] << endl;
        }
        nb_occurrence[i]->spreadsheet_characteristic_print(os , (length->variance > 0. ? false : true));
      }

      if ((characteristics) && (i < characteristics->nb_value) &&
          (characteristics->nb_occurrence)) {
        os << "\n" << SEQ_label[SEQL_NB_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
           << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
        characteristics->nb_occurrence[i]->spreadsheet_characteristic_print(os , (length->variance > 0. ? false : true));
      }

      os << "\n";
      if ((characteristics) && (i < characteristics->nb_value) &&
          (characteristics->nb_occurrence) &&
          (characteristics->nb_occurrence[i]->nb_element > 0)) {
        os << "\t" << SEQ_label[SEQL_NB_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
           << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      }

      if (nb_occurrence) {
        if (length->variance == 0.) {
          os << "\t" << SEQ_label[SEQL_NB_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
             << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION];
        }
        else {
          os << "\t" << SEQ_label[SEQL_MIXTURE_OF] << SEQ_label[SEQL_NB_OCCURRENCE_OF]
             << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
             << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTIONS];
        }
        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->nb_occurrence) &&
            (characteristics->nb_occurrence[i]->nb_element > 0)) {
          os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION];
        }
        if (length->variance == 0.) {
          os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION] << endl;
        }
        else {
          os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE] << " "
             << STAT_label[STATL_FUNCTION] << endl;
        }

        nb_occurrence[i]->spreadsheet_print(os , true , false , false ,
                                            (((characteristics) && (i < characteristics->nb_value) && (characteristics->nb_occurrence) && (characteristics->nb_occurrence[i]->nb_element > 0)) ? characteristics->nb_occurrence[i] : NULL));
      }

      else {
        os << endl;
        characteristics->nb_occurrence[i]->spreadsheet_print(os);
      }
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet CategoricalSequenceProcess.
 *
 *  arguments : prefixe des fichiers, titre des figures,
 *              indice du processus d'observation,
 *              pointeurs sur les lois d'observation empiriques,
 *              la loi marginale empiriques,
 *              les caracteristiques des sequences observees,
 *              la loi empirique des longueurs des sequences et
 *              les lois de l'intervalle de temps residuel.
 *
 *--------------------------------------------------------------*/

bool CategoricalSequenceProcess::plot_print(const char *prefix , const char *title , int process ,
                                            FrequencyDistribution **empirical_observation ,
                                            FrequencyDistribution *marginal_distribution ,
                                            const SequenceCharacteristics *characteristics ,
                                            const FrequencyDistribution *length_distribution ,
                                            Forward **forward) const

{
  bool status = false , start;
  register int i , j , k , m;
  int index_length , nb_histo , nb_dist , histo_index , dist_index , *dist_nb_value;
  double *scale;
  Curves *smoothed_curves;
  const Distribution **pdist;
  const FrequencyDistribution **phisto;
  ostringstream data_file_name[2];


  // ecriture des fichiers de donnees

  if ((index_value) || (characteristics)) {
    if (characteristics) {
      index_length = characteristics->index_value->plot_length_computation();
      if (characteristics->index_value->frequency[index_length - 1] < MAX_FREQUENCY) {
        smoothed_curves = new Curves(*(characteristics->index_value) , 's');
      }
      else {
        smoothed_curves = NULL;
      }
    }

    data_file_name[0] << prefix << process << 0 << ".dat";

    if (index_value) {
      if (characteristics) {
        status = index_value->plot_print((data_file_name[0].str()).c_str() , index_length ,
                                         characteristics->index_value , smoothed_curves);
      }
      else {
        status = index_value->plot_print((data_file_name[0].str()).c_str());
      }
    }

    else {
      status = characteristics->index_value->plot_print((data_file_name[0].str()).c_str() ,
                                                        index_length , smoothed_curves);
    }

    if (characteristics) {
      delete smoothed_curves;
    }
  }

  if (status) {
    pdist = new const Distribution*[6 * nb_value + 3 * nb_state + 2];
    dist_nb_value = new int[6 * nb_value + 3 * nb_state + 2];
    scale = new double[6 * nb_value + 3 * nb_state + 2];
    phisto = new const FrequencyDistribution*[7 * nb_value + nb_state + 3];

    nb_histo = 0;
    nb_dist = 0;

    if (length_distribution) {
      phisto[nb_histo++] = length_distribution;
    }

    if ((first_occurrence) || (characteristics)) {
      for (i = 0;i < nb_value;i++) {
        if ((first_occurrence) && (first_occurrence[i])) {
          pdist[nb_dist] = first_occurrence[i];

          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->first_occurrence[i]->nb_element > 0)) {
            phisto[nb_histo] = characteristics->first_occurrence[i];
            dist_nb_value[nb_dist] = MIN(first_occurrence[i]->nb_value , phisto[nb_histo]->nb_value * 3);
            scale[nb_dist++] = phisto[nb_histo++]->nb_element /
                               (1. - first_occurrence[i]->complement);
//                               first_occurrence[i]->cumul[first_occurrence[i]->nb_value - 1];
          }
          else {
            dist_nb_value[nb_dist] = first_occurrence[i]->nb_value;
            scale[nb_dist++] = 1.;
          }
        }

        else if ((characteristics) && (i < characteristics->nb_value) &&
                 (characteristics->first_occurrence[i]->nb_element > 0)) {
          phisto[nb_histo++] = characteristics->first_occurrence[i];
        }
      }
    }

    if ((recurrence_time) || (characteristics)) {
      for (i = 0;i < nb_value;i++) {
        if ((recurrence_time) && (recurrence_time[i])) {
          pdist[nb_dist] = recurrence_time[i];

          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->recurrence_time[i]->nb_element > 0)) {
            phisto[nb_histo] = characteristics->recurrence_time[i];
            dist_nb_value[nb_dist] = MIN(recurrence_time[i]->nb_value , phisto[nb_histo]->nb_value * 3);
            scale[nb_dist++] = phisto[nb_histo++]->nb_element /
                               (1. - recurrence_time[i]->complement);
          }
          else {
            dist_nb_value[nb_dist] = recurrence_time[i]->nb_value;
            scale[nb_dist++] = 1.;
          }
        }

        else if ((characteristics) && (i < characteristics->nb_value) &&
                 (characteristics->recurrence_time[i]->nb_element > 0)) {
          phisto[nb_histo++] = characteristics->recurrence_time[i];
        }
      }
    }

    if ((sojourn_time) || (characteristics)) {
      for (i = 0;i < nb_value;i++) {
        if ((sojourn_time) && (sojourn_time[i])) {
          pdist[nb_dist] = sojourn_time[i];
          dist_nb_value[nb_dist] = sojourn_time[i]->nb_value;

          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->sojourn_time[i]->nb_element > 0)) {
            phisto[nb_histo] = characteristics->sojourn_time[i];
            if (sojourn_time[i]->cumul[sojourn_time[i]->nb_value - 1] < CUMUL_THRESHOLD) {
              scale[nb_dist++] = phisto[nb_histo++]->nb_element /
                                 (1. - sojourn_time[i]->complement);
//                                 sojourn_time[i]->cumul[sojourn_time[i]->nb_value - 1];
            }
            else {
              scale[nb_dist++] = phisto[nb_histo++]->nb_element;
            }
          }
          else {
            scale[nb_dist++] = 1.;
          }
        }

        else if ((characteristics) && (i < characteristics->nb_value) &&
                 (characteristics->sojourn_time[i]->nb_element > 0)) {
          phisto[nb_histo++] = characteristics->sojourn_time[i];
        }

        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->initial_run) &&
            (characteristics->initial_run[i]->nb_element > 0)) {
          if ((forward) && (forward[i])) {
            pdist[nb_dist] = forward[i];
            dist_nb_value[nb_dist] = forward[i]->nb_value;
            phisto[nb_histo] = characteristics->initial_run[i];
            scale[nb_dist++] = phisto[nb_histo++]->nb_element;
          }

          else {
            phisto[nb_histo++] = characteristics->initial_run[i];
          }
        }

        if ((forward) && (forward[i])) {
          pdist[nb_dist] = forward[i];
          dist_nb_value[nb_dist] = forward[i]->nb_value;

          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->final_run[i]->nb_element > 0)) {
            phisto[nb_histo] = characteristics->final_run[i];
            scale[nb_dist++] = phisto[nb_histo++]->nb_element;
          }
          else {
            scale[nb_dist++] = 1.;
          }
        }

        else if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->final_run[i]->nb_element > 0)) {
          phisto[nb_histo++] = characteristics->final_run[i];
        }
      }
    }

    if ((nb_run) || (nb_occurrence) ||
        ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence))) {
      for (i = 0;i < nb_value;i++) {
        if (nb_run) {
          pdist[nb_dist] = nb_run[i];

          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->nb_run) && (characteristics->nb_run[i]->nb_element > 0)) {
            phisto[nb_histo] = characteristics->nb_run[i];
            dist_nb_value[nb_dist] = nb_run[i]->plot_nb_value_computation(phisto[nb_histo]);
            scale[nb_dist++] = phisto[nb_histo++]->nb_element;
          }
          else {
            dist_nb_value[nb_dist] = nb_run[i]->plot_nb_value_computation();
            scale[nb_dist++] = 1.;
          }
        }

        else if ((characteristics) && (i < characteristics->nb_value) &&
                 (characteristics->nb_run) && (characteristics->nb_run[i]->nb_element > 0)) {
          phisto[nb_histo++] = characteristics->nb_run[i];
        }

        if (nb_occurrence) {
          pdist[nb_dist] = nb_occurrence[i];

          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->nb_occurrence) &&
              (characteristics->nb_occurrence[i]->nb_element > 0)) {
            phisto[nb_histo] = characteristics->nb_occurrence[i];
            dist_nb_value[nb_dist] = nb_occurrence[i]->plot_nb_value_computation(phisto[nb_histo]);
            scale[nb_dist++] = phisto[nb_histo++]->nb_element;
          }
          else {
            dist_nb_value[nb_dist] = nb_occurrence[i]->plot_nb_value_computation();
            scale[nb_dist++] = 1.;
          }
        }

        else if ((characteristics) && (i < characteristics->nb_value) &&
                 (characteristics->nb_occurrence) &&
                 (characteristics->nb_occurrence[i]->nb_element > 0)) {
          phisto[nb_histo++] = characteristics->nb_occurrence[i];
        }
      }
    }

    if (observation) {
      for (i = 0;i < nb_state;i++) {
        pdist[nb_dist] = observation[i];
        dist_nb_value[nb_dist] = observation[i]->nb_value;

        if ((empirical_observation) && (empirical_observation[i]->nb_element > 0)) {
          phisto[nb_histo++] = empirical_observation[i];
          scale[nb_dist++] = empirical_observation[i]->nb_element;
        }
        else {
          scale[nb_dist++] = 1.;
        }
      }

      if (marginal_distribution) {
        if ((weight) && (mixture)) {
          for (i = 0;i < nb_state;i++) {
            pdist[nb_dist] = observation[i];
            dist_nb_value[nb_dist] = observation[i]->nb_value;
            scale[nb_dist++] = weight->mass[i] * marginal_distribution->nb_element;
          }

          pdist[nb_dist] = mixture;
          dist_nb_value[nb_dist] = mixture->nb_value;
          phisto[nb_histo++] = marginal_distribution;
          scale[nb_dist++] = marginal_distribution->nb_element;
        }

        if ((restoration_weight) && (restoration_mixture)) {
          for (i = 0;i < nb_state;i++) {
            pdist[nb_dist] = observation[i];
            dist_nb_value[nb_dist] = observation[i]->nb_value;
            scale[nb_dist++] = restoration_weight->mass[i] * marginal_distribution->nb_element;
          }

          pdist[nb_dist] = restoration_mixture;
          dist_nb_value[nb_dist] = restoration_mixture->nb_value;
          phisto[nb_histo++] = marginal_distribution;
          scale[nb_dist++] = marginal_distribution->nb_element;
        }
      }
    }

    data_file_name[1] << prefix << process << 1 << ".dat";
    status = ::plot_print((data_file_name[1].str()).c_str() , nb_dist , pdist , scale ,
                          dist_nb_value , nb_histo , phisto);

    if (status) {

      // ecriture des fichiers de commandes et des fichiers d'impression

      for (i = 0;i < 2;i++) {
        ostringstream file_name[2];

        switch (i) {
        case 0 :
          file_name[0] << prefix << process << 1 << ".plot";
          break;
        case 1 :
          file_name[0] << prefix << process << 1 << ".print";
          break;
        }

        ofstream out_file((file_name[0].str()).c_str());

        if (i == 1) {
          out_file << "set terminal postscript" << endl;
          file_name[1] << label(prefix) << process << 1 << ".ps";
          out_file << "set output \"" << file_name[1].str() << "\"\n\n";
        }

        out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n";

        if (characteristics) {
          if (characteristics->index_value->frequency[index_length - 1] < MAX_FREQUENCY) {
            out_file << "set title \"";
            if (title) {
              out_file << title << " - ";
            }
            if (process > 0) {
              out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process << " - ";
            }
            out_file << SEQ_label[SEQL_SMOOTHED_OBSERVED_PROBABILITIES] << "\"\n\n";

            if (index_length - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }

            out_file << "plot [0:" << index_length - 1 << "] [0:1] ";

            j = 0;
            for (k = 0;k < nb_value;k++) {
              if (k < characteristics->nb_value) {
                j++;
                out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                         << (index_value ? nb_value : 0) + characteristics->nb_value + k + 1
                         << " title \"" << SEQ_label[SEQL_OBSERVED] << " "
                         << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                         << k << "\" with linespoints";
              }
              if (index_value) {
                j++;
                if (k < characteristics->nb_value) {
                  out_file << ",\\" << endl;
                }
                out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                         << k + 1 << " title \"" << SEQ_label[SEQL_THEORETICAL] << " "
                         << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                         << k << "\" with linespoints";
              }

              if ((j == PLOT_NB_CURVE) && (k < nb_value - 1)) {
                out_file << endl;
                if (i == 0) {
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                }
                out_file << "\nplot [0:" << index_length - 1 << "] [0:1] ";
              }

              else {
                if (k < nb_value - 1) {
                  out_file << ",\\";
                }
                out_file << endl;
              }
            }

            if (index_length - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;
          }

          out_file << "set title";
          if ((title) || (process > 0)) {
            out_file << " \"";
            if (title) {
              out_file << title;
              if (process > 0) {
                out_file << " - ";
              }
            }
            if (process > 0) {
              out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
            }
            out_file << "\"";
          }
          out_file << "\n\n";

          if (index_length - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }

          out_file << "plot [0:" << index_length - 1 << "] [0:1] ";

          j = 0;
          for (k = 0;k < nb_value;k++) {
            if (k < characteristics->nb_value) {
              j++;
              out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                       << (index_value ? nb_value : 0) + k + 1 << " title \"" << SEQ_label[SEQL_OBSERVED] << " "
                       << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                       << k << "\" with linespoints";
            }
            if (index_value) {
              j++;
              if (k < characteristics->nb_value) {
                out_file << ",\\" << endl;
              }
              out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                       << k + 1 << " title \"" << SEQ_label[SEQL_THEORETICAL] << " "
                       << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                       << k << "\" with linespoints";
            }

            if ((j == PLOT_NB_CURVE) && (k < nb_value - 1)) {
              out_file << endl;
              if (i == 0) {
                out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
              }
              out_file << "\nplot [0:" << index_length - 1 << "] [0:1] ";
            }

            else {
              if (k < nb_value - 1) {
                out_file << ",\\";
              }
              out_file << endl;
            }
          }

          if (index_length - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;

          if (length_distribution->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }
          if ((int)(length_distribution->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [0:" << length_distribution->nb_value - 1 << "] [0:"
                   << (int)(length_distribution->max * YSCALE) + 1 << "] \""
                   << label((data_file_name[1].str()).c_str()) << "\" using 1 title \""
                   << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                   << "\" with impulses" << endl;

          if (length_distribution->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
          if ((int)(length_distribution->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics autofreq" << endl;
          }
        }

        else {
          out_file << "set title";
          if ((title) || (process > 0)) {
            out_file << " \"";
            if (title) {
              out_file << title;
              if (process > 0) {
                out_file << " - ";
              }
            }
            if (process > 0) {
              out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
            }
            out_file << "\"";
          }
          out_file << "\n\n";

          if (index_value->length - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }

          out_file << "plot [0:" << index_value->length - 1 << "] [0:1] ";

          for (j = 0;j < nb_value;j++) {
            out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                     << j + 1 << " title \"" << SEQ_label[SEQL_THEORETICAL] << " "
                     << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                     << j << "\" with linespoints";
            if (j < nb_value - 1) {
              out_file << ",\\";
            }
            out_file << endl;
          }

          if (index_value->length - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
        }

        if (i == 1) {
          out_file << "\nset terminal x11" << endl;
        }

        out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }

      histo_index = 1;
      dist_index = 0;

      if ((first_occurrence) || (characteristics)) {
        for (i = 0;i < 2;i++) {
          ostringstream file_name[2];

          switch (i) {
          case 0 :
            file_name[0] << prefix << process << 2 << ".plot";
            break;
          case 1 :
            file_name[0] << prefix << process << 2 << ".print";
            break;
          }

          ofstream out_file((file_name[0].str()).c_str());

          if (i == 1) {
            out_file << "set terminal postscript" << endl;
            file_name[1] << label(prefix) << process << 2 << ".ps";
            out_file << "set output \"" << file_name[1].str() << "\"\n\n";
          }

          out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                   << "set title";
          if ((title) || (process > 0)) {
            out_file << " \"";
            if (title) {
              out_file << title;
              if (process > 0) {
                out_file << " - ";
              }
            }
            if (process > 0) {
              out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
            }
            out_file << "\"";
          }
          out_file << "\n\n";

          j = histo_index;
          k = dist_index;

          start = true;
          for (m = 0;m < nb_value;m++) {
            if ((first_occurrence) && (first_occurrence[m])) {
              if (!start) {
                if (i == 0) {
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                }
                out_file << endl;
              }
              else {
                start = false;
              }

              if (MAX(dist_nb_value[k] , 2) - 1 < TIC_THRESHOLD) {
                out_file << "set xtics 0,1" << endl;
              }

              if ((characteristics) && (m < characteristics->nb_value) &&
                  (characteristics->first_occurrence[m]->nb_element > 0)) {
                out_file << "plot [0:" << MAX(dist_nb_value[k] , 2) - 1 << "] [0:"
                         << (int)(MAX(phisto[j]->max , pdist[k]->max * scale[k]) * YSCALE) + 1
                         << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using " << j + 1
                         << " title \"" << SEQ_label[SEQL_FIRST_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                         << " " << m << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses,\\" << endl;
                out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1
                         << " title \"" << SEQ_label[SEQL_FIRST_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                         << " " << m << " " << STAT_label[STATL_DISTRIBUTION] << "\" with linespoints" << endl;
                j++;
              }

              else {
                out_file << "plot [0:" << MAX(dist_nb_value[k] , 2) - 1 << "] [0:"
                         << MIN(pdist[k]->max * YSCALE , 1.) << "] \""
                         << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1
                         << " title \"" << SEQ_label[SEQL_FIRST_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                         << " " << m << " " << STAT_label[STATL_DISTRIBUTION] << "\" with linespoints" << endl;
              }

              if (MAX(dist_nb_value[k] , 2) - 1 < TIC_THRESHOLD) {
                out_file << "set xtics autofreq" << endl;
              }
              k++;
            }

            else if ((characteristics) && (m < characteristics->nb_value) &&
                     (characteristics->first_occurrence[m]->nb_element > 0)) {
              if (!start) {
                if (i == 0) {
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                }
                out_file << endl;
              }
              else {
                start = false;
              }

              if (MAX(phisto[j]->nb_value , 2) - 1 < TIC_THRESHOLD) {
                out_file << "set xtics 0,1" << endl;
              }
              if ((int)(phisto[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
                out_file << "set ytics 0,1" << endl;
              }

              out_file << "plot [0:" << MAX(phisto[j]->nb_value , 2) - 1 << "] [0:"
                       << (int)(phisto[j]->max * YSCALE) + 1 << "] \""
                       << label((data_file_name[1].str()).c_str()) << "\" using " << j + 1
                       << " title \"" << SEQ_label[SEQL_FIRST_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                       << " " << m << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses" << endl;

              if (MAX(phisto[j]->nb_value , 2) - 1 < TIC_THRESHOLD) {
                out_file << "set xtics autofreq" << endl;
              }
              if ((int)(phisto[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
                out_file << "set ytics autofreq" << endl;
              }
              j++;
            }
          }

          if (i == 1) {
            out_file << "\nset terminal x11" << endl;
          }

          out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
        }

        histo_index = j;
        dist_index = k;
      }

      if ((recurrence_time) || (characteristics)) {
        for (i = 0;i < 2;i++) {
          ostringstream file_name[2];

          switch (i) {
          case 0 :
            file_name[0] << prefix << process << 3 << ".plot";
            break;
          case 1 :
            file_name[0] << prefix << process << 3 << ".print";
            break;
          }

          ofstream out_file((file_name[0].str()).c_str());

          if (i == 1) {
            out_file << "set terminal postscript" << endl;
            file_name[1] << label(prefix) << process << 3 << ".ps";
            out_file << "set output \"" << file_name[1].str() << "\"\n\n";
          }

          out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                   << "set title";
          if ((title) || (process > 0)) {
            out_file << " \"";
            if (title) {
              out_file << title;
              if (process > 0) {
                out_file << " - ";
              }
            }
            if (process > 0) {
              out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
            }
            out_file << "\"";
          }
          out_file << "\n\n";

          j = histo_index;
          k = dist_index;

          start = true;
          for (m = 0;m < nb_value;m++) {
            if ((recurrence_time) && (recurrence_time[m])) {
              if (!start) {
                if (i == 0) {
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                }
                out_file << endl;
              }
              else {
                start = false;
              }

              if (dist_nb_value[k] - 1 < TIC_THRESHOLD) {
                out_file << "set xtics 0,1" << endl;
              }

              if ((characteristics) && (m < characteristics->nb_value) &&
                  (characteristics->recurrence_time[m]->nb_element > 0)) {
                out_file << "plot [0:" << dist_nb_value[k] - 1 << "] [0:"
                         << (int)(MAX(phisto[j]->max , pdist[k]->max * scale[k]) * YSCALE) + 1
                         << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using " << j + 1
                         << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << m
                         << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                         << "\" with impulses,\\" << endl;
                out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1
                         << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << m
                         << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_DISTRIBUTION]
                         << "\" with linespoints" << endl;
                j++;
              }

              else {
                out_file << "plot [0:" << dist_nb_value[k] - 1 << "] [0:"
                         << MIN(pdist[k]->max * YSCALE , 1.) << "] \""
                         << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1
                         << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << m
                         << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_DISTRIBUTION]
                         << "\" with linespoints" << endl;
              }

              if (dist_nb_value[k] - 1 < TIC_THRESHOLD) {
                out_file << "set xtics autofreq" << endl;
              }
              k++;
            }

            else if ((characteristics) && (m < characteristics->nb_value) &&
                     (characteristics->recurrence_time[m]->nb_element > 0)) {
              if (!start) {
                if (i == 0) {
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                }
                out_file << endl;
              }
              else {
                start = false;
              }

              if (phisto[j]->nb_value - 1 < TIC_THRESHOLD) {
                out_file << "set xtics 0,1" << endl;
              }
              if ((int)(phisto[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
                out_file << "set ytics 0,1" << endl;
              }

              out_file << "plot [0:" << phisto[j]->nb_value - 1 << "] [0:"
                       << (int)(phisto[j]->max * YSCALE) + 1 << "] \""
                       << label((data_file_name[1].str()).c_str()) << "\" using " << j + 1
                       << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << m
                       << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                       << "\" with impulses" << endl;

              if (phisto[j]->nb_value - 1 < TIC_THRESHOLD) {
                out_file << "set xtics autofreq" << endl;
              }
              if ((int)(phisto[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
                out_file << "set ytics autofreq" << endl;
              }
              j++;
            }
          }

          if (i == 1) {
            out_file << "\nset terminal x11" << endl;
          }

          out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
        }

        histo_index = j;
        dist_index = k;
      }

      if ((sojourn_time) || (characteristics)) {
        for (i = 0;i < 2;i++) {
          ostringstream file_name[2];

          switch (i) {
          case 0 :
            file_name[0] << prefix << process << 4 << ".plot";
            break;
          case 1 :
            file_name[0] << prefix << process << 4 << ".print";
            break;
          }

          ofstream out_file((file_name[0].str()).c_str());

          if (i == 1) {
            out_file << "set terminal postscript" << endl;
            file_name[1] << label(prefix) << process << 4 << ".ps";
            out_file << "set output \"" << file_name[1].str() << "\"\n\n";
          }

          out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                   << "set title";
          if ((title) || (process > 0)) {
            out_file << " \"";
            if (title) {
              out_file << title;
              if (process > 0) {
                out_file << " - ";
              }
            }
            if (process > 0) {
              out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
            }
            out_file << "\"";
          }
          out_file << "\n\n";

          j = histo_index;
          k = dist_index;

          start = true;
          for (m = 0;m < nb_value;m++) {
            if ((sojourn_time) && (sojourn_time[m])) {
              if (!start) {
                if (i == 0) {
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                }
                out_file << endl;
              }
              else {
                start = false;
              }

              if (dist_nb_value[k] - 1 < TIC_THRESHOLD) {
                out_file << "set xtics 0,1" << endl;
              }

              if ((characteristics) && (m < characteristics->nb_value) &&
                  (characteristics->sojourn_time[m]->nb_element > 0)) {
                out_file << "plot [0:" << dist_nb_value[k] - 1 << "] [0:"
                         << (int)(MAX(phisto[j]->max , pdist[k]->max * scale[k]) * YSCALE) + 1
                         << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using " << j + 1
                         << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << m
                         << " " << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                         << "\" with impulses,\\" << endl;
                out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1
                         << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << m
                         << " " << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION];
                sojourn_time[m]->plot_title_print(out_file);
                out_file << "\" with linespoints" << endl;
                j++;
              }

              else {
                out_file << "plot [0:" << dist_nb_value[k] - 1 << "] [0:"
                         << MIN(pdist[k]->max * YSCALE , 1.) << "] \""
                         << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1
                         << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << m
                         << " " << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION];
                sojourn_time[m]->plot_title_print(out_file);
                out_file << "\" with linespoints" << endl;
              }

              if (dist_nb_value[k] - 1 < TIC_THRESHOLD) {
                out_file << "set xtics autofreq" << endl;
              }
              k++;
            }

            else if ((characteristics) && (m < characteristics->nb_value) &&
                     (characteristics->sojourn_time[m]->nb_element > 0)) {
              if (!start) {
                if (i == 0) {
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                }
                out_file << endl;
              }
              else {
                start = false;
              }

              if (phisto[j]->nb_value - 1 < TIC_THRESHOLD) {
                out_file << "set xtics 0,1" << endl;
              }
              if ((int)(phisto[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
                out_file << "set ytics 0,1" << endl;
              }

              out_file << "plot [0:" << phisto[j]->nb_value - 1 << "] [0:"
                       << (int)(phisto[j]->max * YSCALE) + 1 << "] \""
                       << label((data_file_name[1].str()).c_str()) << "\" using " << j + 1
                       << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << m
                       << " " << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                       << "\" with impulses" << endl;

              if (phisto[j]->nb_value - 1 < TIC_THRESHOLD) {
                out_file << "set xtics autofreq" << endl;
              }
              if ((int)(phisto[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
                out_file << "set ytics autofreq" << endl;
              }
              j++;
            }

            if ((characteristics) && (m < characteristics->nb_value) &&
                (characteristics->initial_run) &&
                (characteristics->initial_run[m]->nb_element > 0)) {
              if (!start) {
                if (i == 0) {
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                }
                out_file << endl;
              }
              else {
                start = false;
              }

              if ((forward) && (forward[m])) {
                if (dist_nb_value[k] - 1 < TIC_THRESHOLD) {
                  out_file << "set xtics 0,1" << endl;
                }

                out_file << "plot [0:" << dist_nb_value[k] - 1 << "] [0:"
                         << (int)(MAX(phisto[j]->max , pdist[k]->max * scale[k]) * YSCALE) + 1
                         << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using " << j + 1
                         << " title \"" << SEQ_label[SEQL_INITIAL_RUN] << " - " << STAT_label[STATL_STATE]
                         << " " << m << " " << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                         << "\" with impulses,\\" << endl;
                out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1
                         << " title \"" << STAT_label[STATL_STATE] << " " << m << " " << SEQ_label[SEQL_FORWARD]
                         << " " << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION]
                         << "\" with linespoints" << endl;

                if (dist_nb_value[k] - 1 < TIC_THRESHOLD) {
                  out_file << "set xtics autofreq" << endl;
                }
                k++;
              }

              else {
                if (phisto[j]->nb_value - 1 < TIC_THRESHOLD) {
                  out_file << "set xtics 0,1" << endl;
                }
                if ((int)(phisto[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
                  out_file << "set ytics 0,1" << endl;
                }

                out_file << "plot [0:" << phisto[j]->nb_value - 1 << "] [0:"
                         << (int)(phisto[j]->max * YSCALE) + 1 << "] \""
                         << label((data_file_name[1].str()).c_str()) << "\" using " << j + 1
                         << " title \"" << SEQ_label[SEQL_INITIAL_RUN] << " - "
                         << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << m
                         << " " << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                         << "\" with impulses" << endl;

                if (phisto[j]->nb_value - 1 < TIC_THRESHOLD) {
                  out_file << "set xtics autofreq" << endl;
                }
                if ((int)(phisto[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
                  out_file << "set ytics autofreq" << endl;
                }
              }

              j++;
            }

            if ((forward) && (forward[m])) {
              if (!start) {
                if (i == 0) {
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                }
                out_file << endl;
              }
              else {
                start = false;
              }

              if (dist_nb_value[k] - 1 < TIC_THRESHOLD) {
                out_file << "set xtics 0,1" << endl;
              }

              if ((characteristics) && (m < characteristics->nb_value) &&
                  (characteristics->final_run[m]->nb_element > 0)) {
                out_file << "plot [0:" << dist_nb_value[k] - 1 << "] [0:"
                         << (int)(MAX(phisto[j]->max , pdist[k]->max * scale[k]) * YSCALE) + 1
                         << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using " << j + 1
                         << " title \"" << SEQ_label[SEQL_FINAL_RUN] << " - " << STAT_label[STATL_STATE]
                         << " " << m << " " << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                         << "\" with impulses,\\" << endl;
                out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1
                         << " title \"" << STAT_label[STATL_STATE] << " " << m << " " << SEQ_label[SEQL_FORWARD]
                         << " " << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION]
                         << "\" with linespoints" << endl;
                j++;
              }

              else {
                out_file << "plot [0:" << dist_nb_value[k] - 1 << "] [0:"
                         << MIN(pdist[k]->max * YSCALE , 1.) << "] \""
                         << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1
                         << " title \"" << STAT_label[STATL_STATE] << " " << m << " " << SEQ_label[SEQL_FORWARD]
                         << " " << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION]
                         << "\" with linespoints" << endl;
              }

              if (dist_nb_value[k] - 1 < TIC_THRESHOLD) {
                out_file << "set xtics autofreq" << endl;
              }
              k++;
            }

            else if ((characteristics) && (m < characteristics->nb_value) &&
                     (characteristics->final_run[m]->nb_element > 0)) {
              if (!start) {
                if (i == 0) {
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                }
                out_file << endl;
              }
              else {
                start = false;
              }

              if (phisto[j]->nb_value - 1 < TIC_THRESHOLD) {
                out_file << "set xtics 0,1" << endl;
              }
              if ((int)(phisto[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
                out_file << "set ytics 0,1" << endl;
              }

              out_file << "plot [0:" << phisto[j]->nb_value - 1 << "] [0:"
                       << (int)(phisto[j]->max * YSCALE) + 1 << "] \""
                       << label((data_file_name[1].str()).c_str()) << "\" using " << j + 1
                       << " title \"" << SEQ_label[SEQL_FINAL_RUN] << " - "
                       << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << m
                       << " " << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                       << "\" with impulses" << endl;

              if (phisto[j]->nb_value - 1 < TIC_THRESHOLD) {
                out_file << "set xtics autofreq" << endl;
              }
              if ((int)(phisto[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
                out_file << "set ytics autofreq" << endl;
              }
              j++;
            }
          }

          if (i == 1) {
            out_file << "\nset terminal x11" << endl;
          }

          out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
        }

        histo_index = j;
        dist_index = k;
      }

      if ((nb_run) || (nb_occurrence) ||
          ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence))) {
        for (i = 0;i < 2;i++) {
          ostringstream file_name[2];

          switch (i) {
          case 0 :
            file_name[0] << prefix << process << 5 << ".plot";
            break;
          case 1 :
            file_name[0] << prefix << process << 5 << ".print";
            break;
          }

          ofstream out_file((file_name[0].str()).c_str());

          if (i == 1) {
            out_file << "set terminal postscript" << endl;
            file_name[1] << label(prefix) << process << 5 << ".ps";
            out_file << "set output \"" << file_name[1].str() << "\"\n\n";
          }

          out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                   << "set title";
          if ((title) || (process > 0)) {
            out_file << " \"";
            if (title) {
              out_file << title;
              if (process > 0) {
                out_file << " - ";
              }
            }
            if (process > 0) {
              out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
            }
            out_file << "\"";
          }
          out_file << "\n\n";

          j = histo_index;
          k = dist_index;

          start = true;
          for (m = 0;m < nb_value;m++) {
            if (nb_run) {
              if (!start) {
                if (i == 0) {
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                }
                out_file << endl;
              }
              else {
                start = false;
              }

              if (dist_nb_value[k] - 1 < TIC_THRESHOLD) {
                out_file << "set xtics 0,1" << endl;
              }

              if ((characteristics) && (m < characteristics->nb_value) &&
                  (characteristics->nb_run) &&
                  (characteristics->nb_run[m]->nb_element > 0)) {
                out_file << "plot [0:" << dist_nb_value[k] - 1 << "] [0:"
                         << (int)(MAX(phisto[j]->max , pdist[k]->max * scale[k]) * YSCALE) + 1
                         << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using " << j + 1
                         << " title \"" << SEQ_label[SEQL_NB_RUN_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                         << " " << m << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                         << "\" with impulses,\\" << endl;
                out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1;
                if (length->variance == 0.) {
                  out_file << " title \"" << SEQ_label[SEQL_NB_RUN_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                           << " " << m << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
                           << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION];
                }
                else {
                  out_file << " title \"" << SEQ_label[SEQL_MIXTURE_OF] << SEQ_label[SEQL_NB_RUN_OF]
                           << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                           << " " << m << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION];
                }
                out_file << "\" with linespoints" << endl;
                j++;
              }

              else {
                out_file << "plot [0:" << dist_nb_value[k] - 1 << "] [0:"
                         << MIN(pdist[k]->max * YSCALE , 1.) << "] \""
                         << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1
                         << " title \"" << SEQ_label[SEQL_NB_RUN_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                         << " " << m << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
                         << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION]
                         << "\" with linespoints" << endl;
              }

              if (dist_nb_value[k] - 1 < TIC_THRESHOLD) {
                out_file << "set xtics autofreq" << endl;
              }
              k++;
            }

            else if ((characteristics) && (m < characteristics->nb_value) &&
                     (characteristics->nb_run) &&
                     (characteristics->nb_run[m]->nb_element > 0)) {
              if (!start) {
                if (i == 0) {
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                }
                out_file << endl;
              }
              else {
                start = false;
              }

              if (phisto[j]->nb_value - 1 < TIC_THRESHOLD) {
                out_file << "set xtics 0,1" << endl;
              }
              if ((int)(phisto[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
                out_file << "set ytics 0,1" << endl;
              }

              out_file << "plot [0:" << phisto[j]->nb_value - 1 << "] [0:"
                       << (int)(phisto[j]->max * YSCALE) + 1 << "] \""
                       << label((data_file_name[1].str()).c_str()) << "\" using " << j + 1
                       << " title \"" << SEQ_label[SEQL_NB_RUN_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                       << " " << m << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                       << "\" with impulses" << endl;

              if (phisto[j]->nb_value - 1 < TIC_THRESHOLD) {
                out_file << "set xtics autofreq" << endl;
              }
              if ((int)(phisto[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
                out_file << "set ytics autofreq" << endl;
              }
              j++;
            }

            if (nb_occurrence) {
              if (!start) {
                if (i == 0) {
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                }
                out_file << endl;
              }
              else {
                start = false;
              }

              if (dist_nb_value[k] - 1 < TIC_THRESHOLD) {
                out_file << "set xtics 0,1" << endl;
              }

              if ((characteristics) && (m < characteristics->nb_value) &&
                  (characteristics->nb_occurrence) &&
                  (characteristics->nb_occurrence[m]->nb_element > 0)) {
                out_file << "plot [0:" << dist_nb_value[k] - 1 << "] [0:"
                         << (int)(MAX(phisto[j]->max , pdist[k]->max * scale[k]) * YSCALE) + 1
                         << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using " << j + 1
                         << " title \"" << SEQ_label[SEQL_NB_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                         << " " << m << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                         << "\" with impulses,\\" << endl;
                out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1;
                if (length->variance == 0.) {
                  out_file << " title \"" << SEQ_label[SEQL_NB_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                           << " " << m << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
                           << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION];
                }
                else {
                  out_file << " title \"" << SEQ_label[SEQL_MIXTURE_OF] << SEQ_label[SEQL_NB_OCCURRENCE_OF]
                           << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                           << " " << m << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION];
                }
                out_file << "\" with linespoints" << endl;
                j++;
              }

              else {
                out_file << "plot [0:" << dist_nb_value[k] - 1 << "] [0:"
                         << MIN(pdist[k]->max * YSCALE , 1.) << "] \""
                         << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1
                         << " title \"" << SEQ_label[SEQL_NB_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                         << " " << m << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
                         << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION]
                         << "\" with linespoints" << endl;
              }

              if (dist_nb_value[k] - 1 < TIC_THRESHOLD) {
                out_file << "set xtics autofreq" << endl;
              }
              k++;
            }

            else if ((characteristics) && (m < characteristics->nb_value) &&
                     (characteristics->nb_occurrence) &&
                     (characteristics->nb_occurrence[m]->nb_element > 0)) {
              if (i == 0) {
                out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
              }
              out_file << endl;

              if (phisto[j]->nb_value - 1 < TIC_THRESHOLD) {
                out_file << "set xtics 0,1" << endl;
              }
              if ((int)(phisto[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
                out_file << "set ytics 0,1" << endl;
              }

              out_file << "plot [0:" << phisto[j]->nb_value - 1 << "] [0:"
                       << (int)(phisto[j]->max * YSCALE) + 1 << "] \""
                       << label((data_file_name[1].str()).c_str()) << "\" using " << j + 1
                       << " title \"" << SEQ_label[SEQL_NB_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                       << " " << m << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                       << "\" with impulses" << endl;

              if (phisto[j]->nb_value - 1 < TIC_THRESHOLD) {
                out_file << "set xtics autofreq" << endl;
              }
              if ((int)(phisto[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
                out_file << "set ytics autofreq" << endl;
              }
              j++;
            }
          }

          if ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence)) {
            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            if (length_distribution->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
            }
            if ((int)(length_distribution->max * YSCALE) + 1 < TIC_THRESHOLD) {
              out_file << "set ytics 0,1" << endl;
            }

            out_file << "plot [0:" << length_distribution->nb_value - 1 << "] [0:"
                     << (int)(length_distribution->max * YSCALE) + 1 << "] \""
                     << label((data_file_name[1].str()).c_str()) << "\" using 1 title \""
                     << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                     << "\" with impulses" << endl;

            if (length_distribution->nb_value - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
            if ((int)(length_distribution->max * YSCALE) + 1 < TIC_THRESHOLD) {
              out_file << "set ytics autofreq" << endl;
            }
          }

          if (i == 1) {
            out_file << "\nset terminal x11" << endl;
          }

          out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
        }

        histo_index = j;
        dist_index = k;
      }

      if (observation) {
        for (i = 0;i < 2;i++) {
          ostringstream file_name[2];

          switch (i) {
          case 0 :
            file_name[0] << prefix << process << 0 << ".plot";
            break;
          case 1 :
            file_name[0] << prefix << process << 0 << ".print";
            break;
          }

          ofstream out_file((file_name[0].str()).c_str());

          if (i == 1) {
            out_file << "set terminal postscript" << endl;
            file_name[1] << label(prefix) << process << 0 << ".ps";
            out_file << "set output \"" << file_name[1].str() << "\"\n\n";
          }

          out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                   << "set title";
          if (title) {
            out_file << " \"" << title << " - " << STAT_label[STATL_OUTPUT_PROCESS]
                     << " " << process << "\"";
          }
          out_file << "\n\n";

          j = histo_index;
          k = dist_index;

          for (m = 0;m < nb_state;m++) {
            if (dist_nb_value[k] - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }

            if ((empirical_observation) && (empirical_observation[m]->nb_element > 0)) {
              out_file << "plot [0:" << dist_nb_value[k] - 1 << "] [0:"
                       << (int)(MAX(phisto[j]->max , pdist[k]->max * scale[k]) * YSCALE) + 1
                       << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using " << j + 1
                       << " title \"" << STAT_label[STATL_STATE] << " " << m << " "
                       << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                       << "\" with impulses,\\" << endl;
              out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1
                       << " title \"" << STAT_label[STATL_STATE] << " " << m << " "
                       << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION]
                       << "\" with linespoints" << endl;
              j++;
            }

            else {
              out_file << "plot [0:" << dist_nb_value[k] - 1 << "] [0:"
                       << MIN(pdist[k]->max * YSCALE , 1.) << "] \""
                       << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1
                       << " title \"" << STAT_label[STATL_STATE] << " " << m << " "
                       << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION]
                       << "\" with linespoints" << endl;
            }

            if (dist_nb_value[k] - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
            k++;

            if ((i == 0) && (m < nb_state - 1)) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
          }

          if (marginal_distribution) {
            if ((weight) && (mixture)) {
              if (i == 0) {
                out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
              }
              out_file << endl;

              out_file << "set title \"";
              if (title) {
                out_file << title << " - ";
              }
              out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process << " - "
                       << STAT_label[STATL_THEORETICAL] << " " << STAT_label[STATL_WEIGHTS] << "\"\n\n";

              if (nb_value - 1 < TIC_THRESHOLD) {
                out_file << "set xtics 0,1" << endl;
              }

              out_file << "plot [0:" << nb_value - 1 << "] [0:"
                       << (int)(MAX(marginal_distribution->max , mixture->max * marginal_distribution->nb_element) * YSCALE) + 1
                       << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using " << j + 1
                       << " title \"" << STAT_label[STATL_MARGINAL] << " "
                       << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses,\\" << endl;
              j++;

              for (m = 0;m < nb_state;m++) {
                out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1
                         << " title \"" << STAT_label[STATL_STATE] << " " << m << " " << STAT_label[STATL_OBSERVATION]
                         << " " << STAT_label[STATL_DISTRIBUTION] << "\" with linespoints,\\" << endl;
                k++;
              }

              out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1
                       << " title \"" << STAT_label[STATL_MIXTURE] << "\" with linespoints" << endl;
              k++;

              if (nb_value - 1 < TIC_THRESHOLD) {
                out_file << "set xtics autofreq" << endl;
              }
            }

            if ((restoration_weight) && (restoration_mixture)) {
              if (i == 0) {
                out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
              }
              out_file << endl;

              out_file << "set title \"";
              if (title) {
                out_file << title << " - ";
              }
              out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process << " - "
                       << STAT_label[STATL_RESTORATION] << " " << STAT_label[STATL_WEIGHTS] << "\"\n\n";

              if (nb_value - 1 < TIC_THRESHOLD) {
                out_file << "set xtics 0,1" << endl;
              }

              out_file << "plot [0:" << nb_value - 1 << "] [0:"
                       << (int)(MAX(marginal_distribution->max , restoration_mixture->max * marginal_distribution->nb_element) * YSCALE) + 1
                       << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using " << j + 1
                       << " title \"" << STAT_label[STATL_MARGINAL] << " "
                       << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses,\\" << endl;
              j++;

              for (m = 0;m < nb_state;m++) {
                out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1
                         << " title \"" << STAT_label[STATL_STATE] << " " << m << " " << STAT_label[STATL_OBSERVATION]
                         << " " << STAT_label[STATL_DISTRIBUTION] << "\" with linespoints,\\" << endl;
                k++;
              }

              out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1
                       << " title \"" << STAT_label[STATL_MIXTURE] << "\" with linespoints" << endl;
              k++;

              if (nb_value - 1 < TIC_THRESHOLD) {
                out_file << "set xtics autofreq" << endl;
              }
            }
          }

          if (i == 1) {
            out_file << "\nset terminal x11" << endl;
          }

          out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
        }
      }

      delete [] pdist;
      delete [] dist_nb_value;
      delete [] scale;
      delete [] phisto;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet CategoricalSequenceProcess.
 *
 *  arguments : reference sur un objet MultiPlotSet, indice du MultiPlot,
 *              indice du processus d'observation,
 *              pointeurs sur les lois d'observation empiriques,
 *              la loi marginale empiriques,
 *              les caracteristiques des sequences observees,
 *              la loi empirique des longueurs des sequences et
 *              les lois de l'intervalle de temps residuel.
 *
 *--------------------------------------------------------------*/

void CategoricalSequenceProcess::plotable_write(MultiPlotSet &plot , int &index , int process ,
                                                FrequencyDistribution **empirical_observation ,
                                                FrequencyDistribution *marginal_distribution ,
                                                const SequenceCharacteristics *characteristics ,
                                                const FrequencyDistribution *length_distribution ,
                                                Forward **forward) const

{
  register int i , j;
  int index_length , dist_nb_value;
  double scale , max;
  Curves *smoothed_curves;
  ostringstream title , legend;


  // calcul du nombre de vues

/*  nb_plot_set = 0;

  if ((index_value) || (characteristics)) {
    nb_plot_set++;

    if (characteristics) {
      index_length = characteristics->index_value->plot_length_computation();

      if (characteristics->index_value->frequency[index_length - 1] < MAX_FREQUENCY) {
        nb_plot_set++;
      }
      nb_plot_set++;
    }
  }

  if ((first_occurrence) || (characteristics)) {
    for (i = 0;i < nb_value;i++) {
      if ((first_occurrence) && (first_occurrence[i])) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->first_occurrence[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }
  }

  if ((recurrence_time) || (characteristics)) {
    for (i = 0;i < nb_value;i++) {
      if ((recurrence_time) && (recurrence_time[i])) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->recurrence_time[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }
  }

  if ((sojourn_time) || (characteristics)) {
    for (i = 0;i < nb_value;i++) {
      if ((sojourn_time) && (sojourn_time[i])) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->sojourn_time[i]->nb_element > 0)) {
        nb_plot_set++;
      }

      if ((characteristics) && (i < characteristics->nb_value) &&
          (characteristics->initial_run) &&
          (characteristics->initial_run[i]->nb_element > 0)) {
        nb_plot_set++;
      }

      if ((forward) && (forward[i])) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->final_run[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }
  }

  if ((nb_run) || (nb_occurrence) ||
      ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence))) {
    for (i = 0;i < nb_value;i++) {
      if (nb_run) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->nb_run) && (characteristics->nb_run[i]->nb_element > 0)) {
        nb_plot_set++;
      }

      if (nb_occurrence) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->nb_occurrence) &&
               (characteristics->nb_occurrence[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }

    if ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence)) {
      nb_plot_set++;
    }
  }

  if (observation) {
    if (empirical_observation) {
      nb_plot_set += nb_state;
    }
    else {
      nb_plot_set++;
    }

    if (marginal_distribution) {
      if ((weight) && (mixture)) {
        nb_plot_set++;
      }
      if ((restoration_weight) && (restoration_mixture)) {
        nb_plot_set++;
      }
    }
  } */

  if ((index_value) || (characteristics)) {
    plot.variable_nb_viewpoint[process]++;
  }
  if ((first_occurrence) || (characteristics)) {
    plot.variable_nb_viewpoint[process]++;
  }
  if ((recurrence_time) || (characteristics)) {
    plot.variable_nb_viewpoint[process]++;
  }
  if ((sojourn_time) || (characteristics)) {
    plot.variable_nb_viewpoint[process]++;
  }
  if ((nb_run) || (nb_occurrence) ||
      ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence))) {
    plot.variable_nb_viewpoint[process]++;
  }

  if (characteristics) {
    index_length = characteristics->index_value->plot_length_computation();

    if (characteristics->index_value->frequency[index_length - 1] < MAX_FREQUENCY) {

      // vue : ajustement intensite lissee

      plot.variable[index] = process;
      plot.viewpoint[index] = INTENSITY;

      smoothed_curves = new Curves(*(characteristics->index_value) , 's');

      title.str("");
      if (process > 0) {
        title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process << " - ";
      }
      title << SEQ_label[SEQL_SMOOTHED_OBSERVED_PROBABILITIES];
      plot[index].title = title.str();

      plot[index].xrange = Range(0 , index_length - 1);
      if (index_length - 1 < TIC_THRESHOLD) {
        plot[index].xtics = 1;
      }

      plot[index].yrange = Range(0. , 1.);

      plot[index].resize(index_value ? nb_value * 2 : nb_value);

      i = 0;
      for (j = 0;j < nb_value;j++) {
        legend.str("");
        legend << SEQ_label[SEQL_OBSERVED] << " "
               << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << j;
        plot[index][i].legend = legend.str();

        plot[index][i].style = "linespoints";

        smoothed_curves->plotable_write(j , plot[index][i]);
        i++;

        if (index_value) {
          legend.str("");
          legend << SEQ_label[SEQL_THEORETICAL] << " "
                 << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << j;
          plot[index][i].legend = legend.str();

          plot[index][i].style = "linespoints";

          index_value->plotable_write(j , plot[index][i]);
          i++;
        }
      }

      delete smoothed_curves;
      index++;
    }

    // vue : ajustement intensite

    plot.variable[index] = process;
    plot.viewpoint[index] = INTENSITY;

    if (process > 0) {
      title.str("");
      title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
      plot[index].title = title.str();
    }

    plot[index].xrange = Range(0 , index_length - 1);
    if (index_length - 1 < TIC_THRESHOLD) {
      plot[index].xtics = 1;
    }

    plot[index].yrange = Range(0. , 1.);

    plot[index].resize(index_value ? nb_value * 2 : nb_value);

    i = 0;
    for (j = 0;j < nb_value;j++) {
      legend.str("");
      legend << SEQ_label[SEQL_OBSERVED] << " "
             << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << j;
      plot[index][i].legend = legend.str();

      plot[index][i].style = "linespoints";

      characteristics->index_value->plotable_write(j , plot[index][i]);
      i++;

      if (index_value) {
        legend.str("");
        legend << SEQ_label[SEQL_THEORETICAL] << " "
               << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << j;
        plot[index][i].legend = legend.str();

        plot[index][i].style = "linespoints";

        index_value->plotable_write(j , plot[index][i]);
        i++;
      }
    }
    index++;

    // vue : loi empirique des longueurs des sequences

    plot.variable[index] = process;
    plot.viewpoint[index] = INTENSITY;

    plot[index].xrange = Range(0 , length_distribution->nb_value - 1);
    plot[index].yrange = Range(0 , ceil(length_distribution->max * YSCALE));

    if (length_distribution->nb_value - 1 < TIC_THRESHOLD) {
      plot[index].xtics = 1;
    }
    if (ceil(length_distribution->max * YSCALE) < TIC_THRESHOLD) {
      plot[index].ytics = 1;
    }

    plot[index].resize(1);

    legend.str("");
    legend << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[index][0].legend = legend.str();

    plot[index][0].style = "impulses";

    length_distribution->plotable_frequency_write(plot[index][0]);
    index++;
  }

  else {

    // vue : intensite theorique

    plot.variable[index] = process;
    plot.viewpoint[index] = INTENSITY;

    if (process > 0) {
      title.str("");
      title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
      plot[index].title = title.str();
    }

    plot[index].xrange = Range(0 , index_length - 1);
    if (index_length - 1 < TIC_THRESHOLD) {
      plot[index].xtics = 1;
    }

    plot[index].yrange = Range(0. , 1.);

    plot[index].resize(nb_value);

    for (i = 0;i < nb_value;i++) {
      legend.str("");
      legend << SEQ_label[SEQL_THEORETICAL] << " "
             << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i;
      plot[index][i].legend = legend.str();

      plot[index][i].style = "linespoints";

      index_value->plotable_write(i , plot[index][i]);
    }
    index++;
  }

  if ((first_occurrence) || (characteristics)) {
    for (i = 0;i < nb_value;i++) {
      if ((first_occurrence) && (first_occurrence[i])) {

        // vue : ajustement loi du temps avant la 1ere occurrence d'un etat/observation

        plot.variable[index] = process;
        plot.viewpoint[index] = FIRST_OCCURRENCE;

        if (process > 0) {
          title.str("");
          title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
          plot[index].title = title.str();
        }

        plot[index].xrange = Range(0 , MAX(first_occurrence[i]->nb_value , 2) - 1);
        if (MAX(first_occurrence[i]->nb_value , 2) - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }

        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->first_occurrence[i]->nb_element > 0)) {
          scale = characteristics->first_occurrence[i]->nb_element /
                  (1. - first_occurrence[i]->complement);
          plot[index].yrange = Range(0 , ceil(MAX(characteristics->first_occurrence[i]->max ,
                                                  first_occurrence[i]->max * scale) * YSCALE));

          plot[index].resize(2);

          legend.str("");
          legend << SEQ_label[SEQL_FIRST_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                 << " " << i << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
          plot[index][0].legend = legend.str();

          plot[index][0].style = "impulses";

          characteristics->first_occurrence[i]->plotable_frequency_write(plot[index][0]);
          j = 1;
        }

        else {
          scale = 1.;
          plot[index].yrange = Range(0. , MIN(first_occurrence[i]->max * YSCALE , 1.));

          plot[index].resize(1);
          j = 0;
        }

        legend.str("");
        legend << SEQ_label[SEQL_FIRST_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
               << " " << i << " " << STAT_label[STATL_DISTRIBUTION];
        plot[index][j].legend = legend.str();

        plot[index][j].style = "linespoints";

        first_occurrence[i]->plotable_mass_write(plot[index][j] , scale);
        index++;
      }

      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->first_occurrence[i]->nb_element > 0)) {

        // vue : loi empirique du temps avant la 1ere occurrence d'un etat/observation

        plot.variable[index] = process;
        plot.viewpoint[index] = FIRST_OCCURRENCE;

        if (process > 0) {
          title.str("");
          title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
          plot[index].title = title.str();
        }

        plot[index].xrange = Range(0 , MAX(characteristics->first_occurrence[i]->nb_value , 2) - 1);
        plot[index].yrange = Range(0 , ceil(characteristics->first_occurrence[i]->max * YSCALE));

        if (MAX(characteristics->first_occurrence[i]->nb_value , 2) - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }
        if (ceil(characteristics->first_occurrence[i]->max * YSCALE) < TIC_THRESHOLD) {
          plot[index].ytics = 1;
        }

        plot[index].resize(1);

        legend.str("");
        legend << SEQ_label[SEQL_FIRST_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
               << " " << i << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        plot[index][0].legend = legend.str();

        plot[index][0].style = "impulses";

        characteristics->first_occurrence[i]->plotable_frequency_write(plot[index][0]);
        index++;
      }
    }
  }

  if ((recurrence_time) || (characteristics)) {
    for (i = 0;i < nb_value;i++) {
      if ((recurrence_time) && (recurrence_time[i])) {

        // vue : ajustement loi empirique du temps de retour dans un etat/observation

        plot.variable[index] = process;
        plot.viewpoint[index] = RECURRENCE_TIME;

        if (process > 0) {
          title.str("");
          title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
          plot[index].title = title.str();
        }

        plot[index].xrange = Range(0 , recurrence_time[i]->nb_value - 1);
        if (recurrence_time[i]->nb_value - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }

        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->recurrence_time[i]->nb_element > 0)) {
          scale = characteristics->recurrence_time[i]->nb_element /
                  (1. - recurrence_time[i]->complement);
          plot[index].yrange = Range(0 , ceil(MAX(characteristics->recurrence_time[i]->max ,
                                                  recurrence_time[i]->max * scale) * YSCALE));

          plot[index].resize(2);

          legend.str("");
          legend << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
                 << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
          plot[index][0].legend = legend.str();

          plot[index][0].style = "impulses";

          characteristics->recurrence_time[i]->plotable_frequency_write(plot[index][0]);
          j = 1;
        }

        else {
          scale = 1.;
          plot[index].yrange = Range(0. , MIN(recurrence_time[i]->max * YSCALE , 1.));

          plot[index].resize(1);
          j = 0;
        }

        legend.str("");
        legend << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
               << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_DISTRIBUTION];
        plot[index][j].legend = legend.str();

        plot[index][j].style = "linespoints";

        recurrence_time[i]->plotable_mass_write(plot[index][j] , scale);
        index++;
      }

      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->recurrence_time[i]->nb_element > 0)) {

        // vue : loi empirique du temps de retour dans un etat/observation

        plot.variable[index] = process;
        plot.viewpoint[index] = RECURRENCE_TIME;

        if (process > 0) {
          title.str("");
          title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
          plot[index].title = title.str();
        }

        plot[index].xrange = Range(0 , characteristics->recurrence_time[i]->nb_value - 1);
        plot[index].yrange = Range(0 , ceil(characteristics->recurrence_time[i]->max * YSCALE));

        if (characteristics->recurrence_time[i]->nb_value - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }
        if (ceil(characteristics->recurrence_time[i]->max * YSCALE) < TIC_THRESHOLD) {
          plot[index].ytics = 1;
        }

        plot[index].resize(1);

        legend.str("");
        legend << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
               << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_DISTRIBUTION];
        plot[index][0].legend = legend.str();

        plot[index][0].style = "impulses";

        characteristics->recurrence_time[i]->plotable_frequency_write(plot[index][0]);
        index++;
      }
    }
  }

  if ((sojourn_time) || (characteristics)) {
    for (i = 0;i < nb_value;i++) {
      if ((sojourn_time) && (sojourn_time[i])) {

        // vue : ajustement loi du temps de sejour dans un etat/observation

        plot.variable[index] = process;
        plot.viewpoint[index] = SOJOURN_TIME;

        if (process > 0) {
          title.str("");
          title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
          plot[index].title = title.str();
        }

        plot[index].xrange = Range(0 , sojourn_time[i]->nb_value - 1);
        if (sojourn_time[i]->nb_value - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }

        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->sojourn_time[i]->nb_element > 0)) {
          scale = characteristics->sojourn_time[i]->nb_element /
                  (1. - sojourn_time[i]->complement);
          plot[index].yrange = Range(0 , ceil(MAX(characteristics->sojourn_time[i]->max ,
                                                  sojourn_time[i]->max * scale) * YSCALE));

          plot[index].resize(2);

          legend.str("");
          legend << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
                 << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
          plot[index][0].legend = legend.str();

          plot[index][0].style = "impulses";

          characteristics->sojourn_time[i]->plotable_frequency_write(plot[index][0]);
          j = 1;
        }

        else {
          scale = 1.;
          plot[index].yrange = Range(0. , MIN(sojourn_time[i]->max * YSCALE , 1.));

          plot[index].resize(1);
          j = 0;
        }

        legend.str("");
        legend << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
               << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION];
        sojourn_time[i]->plot_title_print(legend);
        plot[index][j].legend = legend.str();

        plot[index][j].style = "linespoints";

        sojourn_time[i]->plotable_mass_write(plot[index][j] , scale);
        index++;
      }

      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->sojourn_time[i]->nb_element > 0)) {

        // vue : loi empirique du temps de sejour dans un etat/observation

        plot.variable[index] = process;
        plot.viewpoint[index] = SOJOURN_TIME;

        if (process > 0) {
          title.str("");
          title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
          plot[index].title = title.str();
        }

        plot[index].xrange = Range(0 , characteristics->sojourn_time[i]->nb_value - 1);
        plot[index].yrange = Range(0 , ceil(characteristics->sojourn_time[i]->max * YSCALE));

        if (characteristics->sojourn_time[i]->nb_value - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }
        if (ceil(characteristics->sojourn_time[i]->max * YSCALE) < TIC_THRESHOLD) {
          plot[index].ytics = 1;
        }

        plot[index].resize(1);

        legend.str("");
        legend << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
               << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        plot[index][0].legend = legend.str();

        plot[index][0].style = "impulses";

        characteristics->sojourn_time[i]->plotable_frequency_write(plot[index][0]);
        index++;
      }

      if ((characteristics) && (i < characteristics->nb_value) &&
          (characteristics->initial_run) &&
          (characteristics->initial_run[i]->nb_element > 0)) {
        if ((forward) && (forward[i])) {

          // vue : ajustement loi du temps de sejour residuel

          plot.variable[index] = process;
          plot.viewpoint[index] = SOJOURN_TIME;

          if (process > 0) {
            title.str("");
            title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
            plot[index].title = title.str();
          }

          plot[index].xrange = Range(0 , forward[i]->nb_value - 1);
          plot[index].yrange = Range(0 , ceil(MAX(characteristics->initial_run[i]->max ,
                                                  forward[i]->max * characteristics->initial_run[i]->nb_element) * YSCALE));

          if (forward[i]->nb_value - 1 < TIC_THRESHOLD) {
            plot[index].xtics = 1;
          }

          plot[index].resize(2);

          legend.str("");
          legend << SEQ_label[SEQL_INITIAL_RUN] << " - "
                 << STAT_label[STATL_STATE] << " " << i << " "
                 << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
          plot[index][0].legend = legend.str();

          plot[index][0].style = "impulses";

          characteristics->initial_run[i]->plotable_frequency_write(plot[index][0]);

          legend.str("");
          legend << STAT_label[STATL_STATE] << " " << i << " " << SEQ_label[SEQL_FORWARD] << " "
                 << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION];
          plot[index][1].legend = legend.str();

          plot[index][1].style = "linespoints";

          forward[i]->plotable_mass_write(plot[index][1] , characteristics->initial_run[i]->nb_element);
          index++;
        }

        else {

          // vue : loi empirique du temps de sejour dans le 1ere etat rencontre

          plot.variable[index] = process;
          plot.viewpoint[index] = SOJOURN_TIME;

          if (process > 0) {
            title.str("");
            title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
            plot[index].title = title.str();
          }

          plot[index].xrange = Range(0 , characteristics->initial_run[i]->nb_value - 1);
          plot[index].yrange = Range(0 , ceil(characteristics->initial_run[i]->max * YSCALE));

          if (characteristics->initial_run[i]->nb_value - 1 < TIC_THRESHOLD) {
            plot[index].xtics = 1;
          }
          if (ceil(characteristics->initial_run[i]->max * YSCALE) < TIC_THRESHOLD) {
            plot[index].ytics = 1;
          }

          plot[index].resize(1);

          legend.str("");
          legend << SEQ_label[SEQL_INITIAL_RUN] << " - "
                 << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
                 << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
          plot[index][0].legend = legend.str();

          plot[index][0].style = "impulses";

          characteristics->initial_run[i]->plotable_frequency_write(plot[index][0]);
          index++;
        }
      }

      if ((forward) && (forward[i])) {

        // vue : ajustement loi du temps de sejour dans le dernier etat/observation rencontre

        plot.variable[index] = process;
        plot.viewpoint[index] = SOJOURN_TIME;

        if (process > 0) {
          title.str("");
          title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
          plot[index].title = title.str();
        }

        plot[index].xrange = Range(0 , forward[i]->nb_value - 1);
        if (forward[i]->nb_value - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }

        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->final_run[i]->nb_element > 0)) {
          scale = characteristics->final_run[i]->nb_element;
          plot[index].yrange = Range(0 , ceil(MAX(characteristics->final_run[i]->max ,
                                                  forward[i]->max * scale) * YSCALE));

          plot[index].resize(2);

          legend.str("");
          legend << SEQ_label[SEQL_FINAL_RUN] << " - "
                 << STAT_label[STATL_STATE] << " " << i << " "
                 << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
          plot[index][0].legend = legend.str();

          plot[index][0].style = "impulses";

          characteristics->final_run[i]->plotable_frequency_write(plot[index][0]);
          j = 1;

        }

        else {
          scale = 1.;
          plot[index].yrange = Range(0. , MIN(forward[i]->max * YSCALE , 1.));

          plot[index].resize(1);
          j = 0;
        }

        legend.str("");
        legend << STAT_label[STATL_STATE] << " " << i << " " << SEQ_label[SEQL_FORWARD] << " "
               << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION];
        plot[index][j].legend = legend.str();

        plot[index][j].style = "linespoints";

        forward[i]->plotable_mass_write(plot[index][j] , scale);
        index++;
      }

      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->final_run[i]->nb_element > 0)) {

        // vue : loi empirique du temps de sejour dans le dernier etat/observation rencontre

        plot.variable[index] = process;
        plot.viewpoint[index] = SOJOURN_TIME;

        if (process > 0) {
          title.str("");
          title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
          plot[index].title = title.str();
        }

        plot[index].xrange = Range(0 , characteristics->final_run[i]->nb_value - 1);
        plot[index].yrange = Range(0 , ceil(characteristics->final_run[i]->max * YSCALE));

        if (characteristics->final_run[i]->nb_value - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }
        if (ceil(characteristics->final_run[i]->max * YSCALE) < TIC_THRESHOLD) {
          plot[index].ytics = 1;
        }

        plot[index].resize(1);

        legend.str("");
        legend << SEQ_label[SEQL_FINAL_RUN] << " - "
               << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
               << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        plot[index][0].legend = legend.str();

        plot[index][0].style = "impulses";

        characteristics->final_run[i]->plotable_frequency_write(plot[index][0]);
        index++;
      }
    }
  }

  if ((nb_run) || (nb_occurrence) ||
      ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence))) {
    for (i = 0;i < nb_value;i++) {
      if (nb_run) {

        // vue : ajustement loi du nombre de series d'un etat/observation par sequence

        plot.variable[index] = process;
        plot.viewpoint[index] = COUNTING;

        if (process > 0) {
          title.str("");
          title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
          plot[index].title = title.str();
        }

        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->nb_run) &&
            (characteristics->nb_run[i]->nb_element > 0)) {
          dist_nb_value = nb_run[i]->plot_nb_value_computation(characteristics->nb_run[i]);
          scale = characteristics->nb_run[i]->nb_element;

          plot[index].yrange = Range(0 , ceil(MAX(characteristics->nb_run[i]->max ,
                                                  nb_run[i]->max * scale) * YSCALE));

          plot[index].resize(2);

          legend.str("");
          legend << SEQ_label[SEQL_NB_RUN_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                 << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
          plot[index][0].legend = legend.str();

          plot[index][0].style = "impulses";

          characteristics->nb_run[i]->plotable_frequency_write(plot[index][0]);
          j = 1;
        }

        else {
          dist_nb_value = nb_run[i]->plot_nb_value_computation();
          scale = 1.;

          plot[index].yrange = Range(0. , MIN(nb_run[i]->max * YSCALE , 1.));

          plot[index].resize(1);
          j = 0;
        }

        plot[index].xrange = Range(0 , dist_nb_value);
        if (dist_nb_value < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }

        legend.str("");
        if (length->variance == 0.) {
          legend << SEQ_label[SEQL_NB_RUN_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                 << " " << i << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
                 << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION];
        }
        else {
          legend << SEQ_label[SEQL_MIXTURE_OF] << SEQ_label[SEQL_NB_RUN_OF]
                 << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
                 << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION];
        }
        plot[index][j].legend = legend.str();

        plot[index][j].style = "linespoints";

        nb_run[i]->plotable_mass_write(plot[index][j] , scale);
        index++;
      }

      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->nb_run) &&
               (characteristics->nb_run[i]->nb_element > 0)) {

        // vue : loi empirique du nombre de series d'un etat/observation par sequence

        plot.variable[index] = process;
        plot.viewpoint[index] = COUNTING;

        if (process > 0) {
          title.str("");
          title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
          plot[index].title = title.str();
        }

        plot[index].xrange = Range(0 , MAX(characteristics->nb_run[i]->nb_value , 2) - 1);
        plot[index].yrange = Range(0 , ceil(characteristics->nb_run[i]->max * YSCALE));

        if (MAX(characteristics->nb_run[i]->nb_value , 2) - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }
        if (ceil(characteristics->nb_run[i]->max * YSCALE) < TIC_THRESHOLD) {
          plot[index].ytics = 1;
        }

        plot[index].resize(1);

        legend.str("");
        legend << SEQ_label[SEQL_NB_RUN_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
               << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        plot[index][0].legend = legend.str();

        plot[index][0].style = "impulses";

        characteristics->nb_run[i]->plotable_frequency_write(plot[index][0]);
        index++;
      }

      if (nb_occurrence) {

        // vue : ajustement loi du nombre d'occurrences d'un etat/observation par sequence

        plot.variable[index] = process;
        plot.viewpoint[index] = COUNTING;

        if (process > 0) {
          title.str("");
          title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
          plot[index].title = title.str();
        }

        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->nb_occurrence) &&
            (characteristics->nb_occurrence[i]->nb_element > 0)) {
          dist_nb_value = nb_occurrence[i]->plot_nb_value_computation(characteristics->nb_occurrence[i]);
          scale = characteristics->nb_occurrence[i]->nb_element;

          plot[index].yrange = Range(0 , ceil(MAX(characteristics->nb_occurrence[i]->max ,
                                                  nb_occurrence[i]->max * scale) * YSCALE));

          plot[index].resize(2);

          legend.str("");
          legend << SEQ_label[SEQL_NB_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                 << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
          plot[index][0].legend = legend.str();

          plot[index][0].style = "impulses";

          characteristics->nb_occurrence[i]->plotable_frequency_write(plot[index][0]);
          j = 1;
        }

        else {
          dist_nb_value = nb_occurrence[i]->plot_nb_value_computation();
          scale = 1.;

          plot[index].yrange = Range(0. , MIN(nb_occurrence[i]->max * YSCALE , 1.));

          plot[index].resize(1);
          j = 0;
        }

        plot[index].xrange = Range(0 , dist_nb_value);
        if (dist_nb_value < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }

        legend.str("");
        if (length->variance == 0.) {
          legend << SEQ_label[SEQL_NB_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                 << " " << i << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
                 << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION];
        }
        else {
          legend << SEQ_label[SEQL_MIXTURE_OF] << SEQ_label[SEQL_NB_OCCURRENCE_OF]
                 << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                 << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION];
        }
        plot[index][j].legend = legend.str();

        plot[index][j].style = "linespoints";

        nb_occurrence[i]->plotable_mass_write(plot[index][j] , scale);
        index++;
      }

      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->nb_occurrence) &&
               (characteristics->nb_occurrence[i]->nb_element > 0)) {

        // vue : loi empirique du nombre d'occurrences d'un etat/observation par sequence

        plot.variable[index] = process;
        plot.viewpoint[index] = COUNTING;

        if (process > 0) {
          title.str("");
          title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
          plot[index].title = title.str();
        }

        plot[index].xrange = Range(0 , MAX(characteristics->nb_occurrence[i]->nb_value , 2) - 1);
        plot[index].yrange = Range(0 , ceil(characteristics->nb_occurrence[i]->max * YSCALE));

        if (MAX(characteristics->nb_occurrence[i]->nb_value , 2) - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }
        if (ceil(characteristics->nb_occurrence[i]->max * YSCALE) < TIC_THRESHOLD) {
          plot[index].ytics = 1;
        }

        plot[index].resize(1);

        legend.str("");
        legend << SEQ_label[SEQL_NB_OCCURRENCE_OF] << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
               << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        plot[index][0].legend = legend.str();

        plot[index][0].style = "impulses";

        characteristics->nb_occurrence[i]->plotable_frequency_write(plot[index][0]);
        index++;
      }
    }

    if ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence)) {

      // vue : loi empirique des longueurs des sequences

      plot.variable[index] = process;
      plot.viewpoint[index] = COUNTING;

      plot[index].xrange = Range(0 , length_distribution->nb_value - 1);
      plot[index].yrange = Range(0 , ceil(length_distribution->max * YSCALE));

      if (length_distribution->nb_value - 1 < TIC_THRESHOLD) {
        plot[index].xtics = 1;
      }
      if (ceil(length_distribution->max * YSCALE) < TIC_THRESHOLD) {
        plot[index].ytics = 1;
      }

      plot[index].resize(1);

      legend.str("");
      legend << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      plot[index][0].legend = legend.str();

      plot[index][0].style = "impulses";

      length_distribution->plotable_frequency_write(plot[index][0]);
      index++;
    }
  }

  if (observation) {
    if (empirical_observation) {
      for (i = 0;i < nb_state;i++) {

        // vue : ajustement loi d'observation

        plot.variable[index] = process;
        plot.viewpoint[index] = OBSERVATION;

        title.str("");
        title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
        plot[index].title = title.str();

        plot[index].xrange = Range(0 , observation[i]->nb_value - 1);
        if (observation[i]->nb_value - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }

        if (empirical_observation[i]->nb_element > 0) {
          scale = empirical_observation[i]->nb_element;
          plot[index].yrange = Range(0 , ceil(MAX(empirical_observation[i]->max ,
                                                  observation[i]->max * scale) * YSCALE));

          plot[index].resize(2);

          legend.str("");
          legend << STAT_label[STATL_STATE] << " " << i << " "
                 << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
          plot[index][0].legend = legend.str();

          plot[index][0].style = "impulses";

          empirical_observation[i]->plotable_frequency_write(plot[index][0]);
          j = 1;
        }

        else {
          scale = 1;
          plot[index].yrange = Range(0 , MIN(observation[i]->max * YSCALE , 1.));

          plot[index].resize(1);
          j = 0;
        }

        legend.str("");
        legend << STAT_label[STATL_STATE] << " " << i << " "
               << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
        plot[index][j].legend = legend.str();

        plot[index][j].style = "linespoints";

        observation[i]->plotable_mass_write(plot[index][j] , scale);
        index++;
      }
    }

    else {

      // vue : lois d'observation

      plot.variable[index] = process;
      plot.viewpoint[index] = OBSERVATION;

      title.str("");
      title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
      plot[index].title = title.str();

      plot[index].xrange = Range(0 , nb_value - 1);
      if (nb_value - 1 < TIC_THRESHOLD) {
        plot[index].xtics = 1;
      }

      max = observation[0]->max;
      for (i = 1;i < nb_state;i++) {
        if (observation[i]->max > max) {
          max = observation[i]->max;
        }
      }
      plot[index].yrange = Range(0 , MIN(max * YSCALE , 1.));

      plot[index].resize(nb_state);

      for (i = 0;i < nb_state;i++) {
        legend.str("");
        legend << STAT_label[STATL_STATE] << " " << i << " "
               << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
        plot[index][i].legend = legend.str();

        plot[index][i].style = "linespoints";

        observation[i]->plotable_mass_write(plot[index][i]);
      }

      index++;
    }

    if (marginal_distribution) {
      if ((weight) && (mixture)) {

        // vue : ajustement melange de lois d'observation

        title.str("");
        title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process << " - "
              << STAT_label[STATL_THEORETICAL] << " " << STAT_label[STATL_WEIGHTS];
        plot[index].title = title.str();

        plot[index].xrange = Range(0 , nb_value - 1);
        if (nb_value - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }

        plot[index].yrange = Range(0 , ceil(MAX(marginal_distribution->max ,
                                                mixture->max * marginal_distribution->nb_element) * YSCALE));

        plot[index].resize(nb_state + 2);

        legend.str("");
        legend << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        plot[index][0].legend = legend.str();

        plot[index][0].style = "impulses";

        marginal_distribution->plotable_frequency_write(plot[index][0]);

        for (i = 0;i < nb_state;i++) {
          legend.str("");
          legend << STAT_label[STATL_STATE] << " " << i << " " << STAT_label[STATL_OBSERVATION] << " "
                 << STAT_label[STATL_DISTRIBUTION];
          observation[i]->plot_title_print(legend);
          plot[index][i + 1].legend = legend.str();

          plot[index][i + 1].style = "linespoints";

          observation[i]->plotable_mass_write(plot[index][i + 1] ,
                                              weight->mass[i] * marginal_distribution->nb_element);
        }

        plot[index][nb_state + 2].legend = STAT_label[STATL_MIXTURE];

        plot[index][nb_state + 2].style = "linespoints";

        mixture->plotable_mass_write(plot[index][nb_state + 2] ,
                                     marginal_distribution->nb_element);

        index++;
      }

      if ((restoration_weight) && (restoration_mixture)) {

        // vue : ajustement melange de lois d'observation (poids issus de la restauration)

        title.str("");
        if (process > 0) {
          title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
        }
        title << STAT_label[STATL_RESTORATION] << " " << STAT_label[STATL_WEIGHTS];
        plot[index].title = title.str();

        plot[index].xrange = Range(0 , nb_value - 1);
        if (nb_value - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }

        plot[index].yrange = Range(0 , ceil(MAX(marginal_distribution->max ,
                                                restoration_mixture->max * marginal_distribution->nb_element) * YSCALE));

        plot[index].resize(nb_state + 2);

        legend.str("");
        legend << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        plot[index][0].legend = legend.str();

        plot[index][0].style = "impulses";

        marginal_distribution->plotable_frequency_write(plot[index][0]);

        for (i = 0;i < nb_state;i++) {
          legend.str("");
          legend << STAT_label[STATL_STATE] << " " << i << " " << STAT_label[STATL_OBSERVATION] << " "
                 << STAT_label[STATL_DISTRIBUTION];
          observation[i]->plot_title_print(legend);
          plot[index][i + 1].legend = legend.str();

          plot[index][i + 1].style = "linespoints";

          observation[i]->plotable_mass_write(plot[index][i + 1] ,
                                              restoration_weight->mass[i] * marginal_distribution->nb_element);
        }

        plot[index][nb_state + 2].legend = STAT_label[STATL_MIXTURE];

        plot[index][nb_state + 2].style = "linespoints";

        restoration_mixture->plotable_mass_write(plot[index][nb_state + 2] ,
                                                 marginal_distribution->nb_element);

        index++;
      }
    }
  }
}


};  // namespace sequence_analysis
