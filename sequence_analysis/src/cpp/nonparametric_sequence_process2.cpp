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



#include <sstream>
#include <iomanip>
// #include <rw/vstream.h>
// #include <rw/rwfile.h>
#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequences.h"
#include "sequence_label.h"
#include "tool/config.h"

using namespace std;


extern int column_width(int value);
extern int column_width(int nb_value , const double *value , double scale = 1.);
extern char* label(const char *file_name);



/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Nonparametric_sequence_process.
 *
 *  arguments : stream, indice du processus d'observation,
 *              pointeurs sur les lois d'observation empiriques et
 *              sur les caracteristiques des sequences observees,
 *              flag niveau de detail, flag fichier ,
 *              pointeurs sur les lois de l'intervalle de temps residuel.
 *
 *--------------------------------------------------------------*/

ostream& Nonparametric_sequence_process::ascii_print(ostream &os , int process ,
                                                     Histogram **empirical_observation ,
                                                     const Sequence_characteristics *characteristics ,
                                                     bool exhaustive , bool file_flag , Forward **forward) const

{
  register int i , j;
  int buff , width[2];
  long old_adjust;
  double *pmass;


  old_adjust = os.setf(ios::left , ios::adjustfield);

  if (observation) {
    for (i = 0;i < nb_state;i++) {
      os << "\n" << STAT_word[STATW_STATE] << " " << i << " "
         << STAT_word[STATW_OBSERVATION_DISTRIBUTION] << endl;
      pmass = observation[i]->mass + observation[i]->offset;
      for (j = observation[i]->offset;j < observation[i]->nb_value;j++) {
        if (*pmass > 0.) {
          os << STAT_word[STATW_OUTPUT] << " " << j << " : " << *pmass << endl;
        }
        pmass++;
      }

      if ((empirical_observation) && (empirical_observation[i]->nb_element > 0) &&
          (exhaustive)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "   | " << STAT_label[STATL_STATE] << " " << i << " "
           << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_HISTOGRAM]
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
                               (characteristics ? characteristics->index_value : 0));
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
          os << SEQ_label[process == 0 ? SEQL_STATE_NO_OCCURRENCE : SEQL_OUTPUT_NO_OCCURRENCE]
             << " " << i << ": " << no_occurrence[i] << endl;
        }

        if (first_occurrence[i]) {
          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << SEQ_label[process == 0 ? SEQL_STATE_FIRST_OCCURRENCE : SEQL_OUTPUT_FIRST_OCCURRENCE]
             << " " << i << " " << STAT_label[STATL_DISTRIBUTION] << endl;
          first_occurrence[i]->ascii_characteristic_print(os , false , file_flag);
        }
      }

      if ((characteristics) && (i < characteristics->nb_value)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[process == 0 ? SEQL_STATE_FIRST_OCCURRENCE : SEQL_OUTPUT_FIRST_OCCURRENCE]
           << " " << i << " " << STAT_label[STATL_HISTOGRAM] << " - ";
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
          os << " | " << SEQ_label[process == 0 ? SEQL_STATE_FIRST_OCCURRENCE : SEQL_OUTPUT_FIRST_OCCURRENCE]
             << " " << i << " " << STAT_label[STATL_HISTOGRAM];
        }

        if ((first_occurrence) && (first_occurrence[i])) {
          os << " | " << SEQ_label[process == 0 ? SEQL_STATE_FIRST_OCCURRENCE : SEQL_OUTPUT_FIRST_OCCURRENCE]
             << " " << i << " " << STAT_label[STATL_DISTRIBUTION];
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->first_occurrence[i]->nb_element > 0)) {
            os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
               << STAT_label[STATL_FUNCTION];
          }
          os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION] << endl;

          first_occurrence[i]->ascii_print(os , file_flag , true , true ,
                                           (((characteristics) && (i < characteristics->nb_value) && (characteristics->first_occurrence[i]->nb_element > 0)) ? characteristics->first_occurrence[i] : 0));
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
          os << SEQ_label[process == 0 ? SEQL_LEAVING_STATE : SEQL_LEAVING_OUTPUT]
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
           << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
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
             << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_HISTOGRAM];
        }

        if ((recurrence_time) && (recurrence_time[i])) {
          os << " | " << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_DISTRIBUTION];
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->recurrence_time[i]->nb_element > 0)) {
            os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
               << STAT_label[STATL_FUNCTION];
          }
          os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION] << endl;

          recurrence_time[i]->ascii_print(os , file_flag , true , true ,
                                          (((characteristics) && (i < characteristics->nb_value) && (characteristics->recurrence_time[i]->nb_element > 0)) ?
                                           characteristics->recurrence_time[i] : 0));
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
          os << SEQ_label[process == 0 ? SEQL_STATE_ABSORPTION : SEQL_OUTPUT_ABSORPTION]
             << " " << i << ": " << absorption[i] << endl;
        }

        if (sojourn_time[i]) {
          if (sojourn_time[i]->ident != NONPARAMETRIC) {
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
          sojourn_time[i]->ascii_characteristic_print(os , false , file_flag);
        }
      }

      if ((characteristics) && (i < characteristics->nb_value)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
           << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
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
             << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_HISTOGRAM];
        }

        if ((sojourn_time) && (sojourn_time[i])) {
          os << " | " << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION];
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->sojourn_time[i]->nb_element > 0)) {
            os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
               << STAT_label[STATL_FUNCTION];
          }
          os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION] << endl;

          sojourn_time[i]->Distribution::ascii_print(os , file_flag , true , false ,
                                                     (((characteristics) && (i < characteristics->nb_value) && (characteristics->sojourn_time[i]->nb_element > 0)) ?
                                                      characteristics->sojourn_time[i] : 0));
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
           << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
        characteristics->initial_run[i]->ascii_characteristic_print(os , false , file_flag);

        if ((characteristics->initial_run[i]->nb_element > 0) && (exhaustive)) {
          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << "   | " << SEQ_label[SEQL_INITIAL_RUN] << " - "
             << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_HISTOGRAM];

          if ((forward) && (forward[i])) {
            os << " | " << STAT_label[STATL_STATE] << " " << i << " " << SEQ_label[SEQL_FORWARD] << " "
               << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION]
               << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
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
           << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
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
             << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_HISTOGRAM];
        }

        if ((forward) && (forward[i])) {
          os << " | " << STAT_label[STATL_STATE] << " " << i << " " << SEQ_label[SEQL_FORWARD] << " "
             << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION];
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->final_run[i]->nb_element > 0)) {
            os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
               << STAT_label[STATL_FUNCTION];
          }
          os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION] << endl;

          forward[i]->Distribution::ascii_print(os , file_flag , true , false ,
                                                (((characteristics) && (i < characteristics->nb_value) && (characteristics->final_run[i]->nb_element > 0)) ?
                                                 characteristics->final_run[i] : 0));
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
          os << SEQ_label[process == 0 ? SEQL_STATE_NB_RUN : SEQL_OUTPUT_NB_RUN]
             << " " << i << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
             << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
        }
        else {
          os << SEQ_label[SEQL_MIXTURE_OF]
             << SEQ_label[process == 0 ? SEQL_STATE_NB_RUN : SEQL_OUTPUT_NB_RUN] << " " << i << " "
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
        os << SEQ_label[process == 0 ? SEQL_STATE_NB_RUN : SEQL_OUTPUT_NB_RUN] << " " << i << " "
           << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
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
          os << " | " << SEQ_label[process == 0 ? SEQL_STATE_NB_RUN : SEQL_OUTPUT_NB_RUN]
             << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_HISTOGRAM];
        }

        if (nb_run) {
          if (length->variance == 0.) {
            os << " | " << SEQ_label[process == 0 ? SEQL_STATE_NB_RUN : SEQL_OUTPUT_NB_RUN]
               << " " << i << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
               << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION];
          }
          else {
            os << " | " << SEQ_label[SEQL_MIXTURE_OF]
               << SEQ_label[process == 0 ? SEQL_STATE_NB_RUN : SEQL_OUTPUT_NB_RUN] << " " << i << " "
               << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTIONS];
          }
          if ((characteristics) && (i < characteristics->nb_value) && (characteristics->nb_run) &&
              (characteristics->nb_run[i]->nb_element > 0)) {
            os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
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
                                 (((characteristics) && (i < characteristics->nb_value) && (characteristics->nb_run) && (characteristics->nb_run[i]->nb_element > 0)) ? characteristics->nb_run[i] : 0));
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
          os << SEQ_label[process == 0 ? SEQL_STATE_NB_OCCURRENCE : SEQL_OUTPUT_NB_OCCURRENCE]
             << " " << i << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
             << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
        }
        else {
          os << SEQ_label[SEQL_MIXTURE_OF]
             << SEQ_label[process == 0 ? SEQL_STATE_NB_OCCURRENCE : SEQL_OUTPUT_NB_OCCURRENCE]
             << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTIONS] << endl;
        }
        nb_occurrence[i]->ascii_characteristic_print(os , (length->variance > 0. ? false : true) , file_flag);
      }

      if ((characteristics) && (i < characteristics->nb_value) &&
          (characteristics->nb_occurrence)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[process == 0 ? SEQL_STATE_NB_OCCURRENCE : SEQL_OUTPUT_NB_OCCURRENCE]
           << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
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
          os << " | " << SEQ_label[process == 0 ? SEQL_STATE_NB_OCCURRENCE : SEQL_OUTPUT_NB_OCCURRENCE]
             << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_HISTOGRAM];
        }

        if (nb_occurrence) {
          if (length->variance == 0.) {
            os << " | " << SEQ_label[process == 0 ? SEQL_STATE_NB_OCCURRENCE : SEQL_OUTPUT_NB_OCCURRENCE]
               << " " << i << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
               << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION];
          }
          else {
            os << " | " << SEQ_label[SEQL_MIXTURE_OF]
               << SEQ_label[process == 0 ? SEQL_STATE_NB_OCCURRENCE : SEQL_OUTPUT_NB_OCCURRENCE]
               << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTIONS];
          }
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->nb_occurrence) &&
              (characteristics->nb_occurrence[i]->nb_element > 0)) {
            os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
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
                                        (((characteristics) && (i < characteristics->nb_value) && (characteristics->nb_occurrence) && (characteristics->nb_occurrence[i]->nb_element > 0)) ? characteristics->nb_occurrence[i] : 0));
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
 *  Ecriture d'un objet Nonparametric_sequence_process au format tableur.
 *
 *  arguments : stream, indice du processus d'observation,
 *              pointeurs sur les lois d'observation empiriques,
 *              sur les caracteristiques des sequences observees et
 *              sur les lois de l'intervalle de temps residuel.
 *
 *--------------------------------------------------------------*/

ostream& Nonparametric_sequence_process::spreadsheet_print(ostream &os , int process ,
                                                           Histogram **empirical_observation ,
                                                           const Sequence_characteristics *characteristics ,
                                                           Forward **forward) const

{
  register int i , j;
  double *pmass;
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
           << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_HISTOGRAM]
           << "\t" << STAT_label[STATL_STATE] << " " << i << " "
           << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION] << endl;

        observation[i]->spreadsheet_print(os , false , false , false , empirical_observation[i]);
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
      index_value->spreadsheet_print(os , (characteristics ? characteristics->index_value : 0));
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
          os << "\n" << SEQ_label[process == 0 ? SEQL_STATE_NO_OCCURRENCE : SEQL_OUTPUT_NO_OCCURRENCE]
             << " " << i << ": " << no_occurrence[i] << endl;
        }

        if (first_occurrence[i]) {
          os << "\n" << SEQ_label[process == 0 ? SEQL_STATE_FIRST_OCCURRENCE : SEQL_OUTPUT_FIRST_OCCURRENCE]
             << " " << i << " " << STAT_label[STATL_DISTRIBUTION] << endl;
          first_occurrence[i]->spreadsheet_characteristic_print(os);
        }
      }

      if ((characteristics) && (i < characteristics->nb_value)) {
        os << "\n" << SEQ_label[process == 0 ? SEQL_STATE_FIRST_OCCURRENCE : SEQL_OUTPUT_FIRST_OCCURRENCE]
           << " " << i << " " << STAT_label[STATL_HISTOGRAM] << "\t";
        characteristics->first_occurrence[i]->spreadsheet_characteristic_print(os);
      }

      if (((first_occurrence) && (first_occurrence[i])) ||
          ((characteristics) && (i < characteristics->nb_value) &&
           (characteristics->first_occurrence[i]->nb_element > 0))) {
        os << "\n";
        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->first_occurrence[i]->nb_element > 0)) {
          os << "\t" << SEQ_label[process == 0 ? SEQL_STATE_FIRST_OCCURRENCE : SEQL_OUTPUT_FIRST_OCCURRENCE]
             << " " << i << " " << STAT_label[STATL_HISTOGRAM];
        }

        if ((first_occurrence) && (first_occurrence[i])) {
          os << "\t" << SEQ_label[process == 0 ? SEQL_STATE_FIRST_OCCURRENCE : SEQL_OUTPUT_FIRST_OCCURRENCE]
             << " " << i << " " << STAT_label[STATL_DISTRIBUTION];
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->first_occurrence[i]->nb_element > 0)) {
            os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
               << STAT_label[STATL_FUNCTION];
          }
          os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION] << endl;

          first_occurrence[i]->spreadsheet_print(os , true , false , true ,
                                                 (((characteristics) && (i < characteristics->nb_value) && (characteristics->first_occurrence[i]->nb_element > 0)) ? characteristics->first_occurrence[i] : 0));
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
          os << "\n" << SEQ_label[process == 0 ? SEQL_LEAVING_STATE : SEQL_LEAVING_OUTPUT]
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
           << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
        characteristics->recurrence_time[i]->spreadsheet_characteristic_print(os);
      }

      if (((recurrence_time) && (recurrence_time[i])) ||
          ((characteristics) && (i < characteristics->nb_value) &&
           (characteristics->recurrence_time[i]->nb_element > 0))) {
        os << "\n";
        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->recurrence_time[i]->nb_element > 0)) {
          os << "\t" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_HISTOGRAM];
        }

        if ((recurrence_time) && (recurrence_time[i])) {
          os << "\t" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_DISTRIBUTION];
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->recurrence_time[i]->nb_element > 0)) {
            os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
               << STAT_label[STATL_FUNCTION];
          }
          os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION] << endl;

          recurrence_time[i]->spreadsheet_print(os , true , false , true ,
                                                (((characteristics) && (i < characteristics->nb_value) && (characteristics->recurrence_time[i]->nb_element > 0)) ?
                                                 characteristics->recurrence_time[i] : 0));
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
          os << "\n" << SEQ_label[process == 0 ? SEQL_STATE_ABSORPTION : SEQL_OUTPUT_ABSORPTION]
             << " " << i << ": " << absorption[i] << endl;
        }

        if (sojourn_time[i]) {
          if (sojourn_time[i]->ident != NONPARAMETRIC) {
            os << "\n" << STAT_word[STATW_STATE] << " " << i << "\t"
               << SEQ_word[SEQW_OCCUPANCY_DISTRIBUTION] << endl;
            sojourn_time[i]->spreadsheet_print(os);
          }
          else {
            os << "\n" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
               << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
          }
          sojourn_time[i]->spreadsheet_characteristic_print(os);
        }
      }

      if ((characteristics) && (i < characteristics->nb_value)) {
        os << "\n" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
           << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
        characteristics->sojourn_time[i]->spreadsheet_characteristic_print(os);
      }

      if (((sojourn_time) && (sojourn_time[i])) ||
          ((characteristics) && (i < characteristics->nb_value) &&
           (characteristics->sojourn_time[i]->nb_element > 0))) {
        os << "\n";
        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->sojourn_time[i]->nb_element > 0)) {
          os << "\t" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_HISTOGRAM];
        }

        if ((sojourn_time) && (sojourn_time[i])) {
          os << "\t" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION];
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->sojourn_time[i]->nb_element > 0)) {
            os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
               << STAT_label[STATL_FUNCTION];
          }
          os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION] << endl;

          sojourn_time[i]->Distribution::spreadsheet_print(os , true , false , false ,
                                                           (((characteristics) && (i < characteristics->nb_value) && (characteristics->sojourn_time[i]->nb_element > 0)) ?
                                                            characteristics->sojourn_time[i] : 0));
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
           << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
        characteristics->initial_run[i]->spreadsheet_characteristic_print(os);

        if (characteristics->initial_run[i]->nb_element > 0) {
          os << "\n\t" << SEQ_label[SEQL_INITIAL_RUN] << " - "
             << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_HISTOGRAM];

          if ((forward) && (forward[i])) {
            os << "\t" << STAT_label[STATL_STATE] << " " << i << " " << SEQ_label[SEQL_FORWARD] << " "
               << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION]
               << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
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
           << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
        characteristics->final_run[i]->spreadsheet_characteristic_print(os);
      }

      if (((forward) && (forward[i])) || ((characteristics) && (i < characteristics->nb_value) &&
           (characteristics->final_run[i]->nb_element > 0))) {
        os << "\n";
        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->final_run[i]->nb_element > 0)) {
          os << "\t" << SEQ_label[SEQL_FINAL_RUN] << " - "
             << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
             << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_HISTOGRAM];
        }

        if ((forward) && (forward[i])) {
          os << "\t" << STAT_label[STATL_STATE] << " " << i << " " << SEQ_label[SEQL_FORWARD] << " "
             << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_DISTRIBUTION];
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->final_run[i]->nb_element > 0)) {
            os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
               << STAT_label[STATL_FUNCTION];
          }
          os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION] << endl;

          forward[i]->Distribution::spreadsheet_print(os , true , false , false ,
                                                      (((characteristics) && (i < characteristics->nb_value) && (characteristics->final_run[i]->nb_element > 0)) ?
                                                       characteristics->final_run[i] : 0));

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
          os << "\n" << SEQ_label[process == 0 ? SEQL_STATE_NB_RUN : SEQL_OUTPUT_NB_RUN]
             << " " << i << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
             << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
        }
        else {
          os << "\n" << SEQ_label[SEQL_MIXTURE_OF]
             << SEQ_label[process == 0 ? SEQL_STATE_NB_RUN : SEQL_OUTPUT_NB_RUN] << " " << i << " "
             << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTIONS] << endl;
        }
        nb_run[i]->spreadsheet_characteristic_print(os , (length->variance > 0. ? false : true));
      }

      if ((characteristics) && (i < characteristics->nb_value) && (characteristics->nb_run)) {
        os << "\n" << SEQ_label[process == 0 ? SEQL_STATE_NB_RUN : SEQL_OUTPUT_NB_RUN]
           << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
        characteristics->nb_run[i]->spreadsheet_characteristic_print(os , (length->variance > 0. ? false : true));
      }

      os << "\n";
      if ((characteristics) && (i < characteristics->nb_value) && (characteristics->nb_run) &&
          (characteristics->nb_run[i]->nb_element > 0)) {
        os << "\t" << SEQ_label[process == 0 ? SEQL_STATE_NB_RUN : SEQL_OUTPUT_NB_RUN]
           << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_HISTOGRAM];
      }

      if (nb_run) {
        if (length->variance == 0.) {
          os << "\t" << SEQ_label[process == 0 ? SEQL_STATE_NB_RUN : SEQL_OUTPUT_NB_RUN]
             << " " << i << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
             << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION];
        }
        else {
          os << "\t" << SEQ_label[SEQL_MIXTURE_OF]
             << SEQ_label[process == 0 ? SEQL_STATE_NB_RUN : SEQL_OUTPUT_NB_RUN] << " " << i << " "
             << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTIONS];
        }
        if ((characteristics) && (i < characteristics->nb_value) && (characteristics->nb_run) &&
            (characteristics->nb_run[i]->nb_element > 0)) {
          os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
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
                                     (((characteristics) && (i < characteristics->nb_value) && (characteristics->nb_run) && (characteristics->nb_run[i]->nb_element > 0)) ? characteristics->nb_run[i] : 0));
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
          os << "\n" << SEQ_label[process == 0 ? SEQL_STATE_NB_OCCURRENCE : SEQL_OUTPUT_NB_OCCURRENCE]
             << " " << i << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
             << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
        }
        else {
          os << "\n" << SEQ_label[SEQL_MIXTURE_OF]
             << SEQ_label[process == 0 ? SEQL_STATE_NB_OCCURRENCE : SEQL_OUTPUT_NB_OCCURRENCE]
             << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTIONS] << endl;
        }
        nb_occurrence[i]->spreadsheet_characteristic_print(os , (length->variance > 0. ? false : true));
      }

      if ((characteristics) && (i < characteristics->nb_value) &&
          (characteristics->nb_occurrence)) {
        os << "\n" << SEQ_label[process == 0 ? SEQL_STATE_NB_OCCURRENCE : SEQL_OUTPUT_NB_OCCURRENCE]
           << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
        characteristics->nb_occurrence[i]->spreadsheet_characteristic_print(os , (length->variance > 0. ? false : true));
      }

      os << "\n";
      if ((characteristics) && (i < characteristics->nb_value) &&
          (characteristics->nb_occurrence) &&
          (characteristics->nb_occurrence[i]->nb_element > 0)) {
        os << "\t" << SEQ_label[process == 0 ? SEQL_STATE_NB_OCCURRENCE : SEQL_OUTPUT_NB_OCCURRENCE]
           << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_HISTOGRAM];
      }

      if (nb_occurrence) {
        if (length->variance == 0.) {
          os << "\t" << SEQ_label[process == 0 ? SEQL_STATE_NB_OCCURRENCE : SEQL_OUTPUT_NB_OCCURRENCE]
             << " " << i << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
             << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION];
        }
        else {
          os << "\t" << SEQ_label[SEQL_MIXTURE_OF]
             << SEQ_label[process == 0 ? SEQL_STATE_NB_OCCURRENCE : SEQL_OUTPUT_NB_OCCURRENCE]
             << " " << i << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTIONS];
        }
        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->nb_occurrence) &&
            (characteristics->nb_occurrence[i]->nb_element > 0)) {
          os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
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
                                            (((characteristics) && (i < characteristics->nb_value) && (characteristics->nb_occurrence) && (characteristics->nb_occurrence[i]->nb_element > 0)) ? characteristics->nb_occurrence[i] : 0));
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
 *  Sortie Gnuplot d'un objet Nonparametric_sequence_process.
 *
 *  arguments : prefixe des fichiers, titre des figures,
 *              indice du processus d'observation,
 *              pointeurs sur les lois d'observation empiriques,
 *              sur les caracteristiques des sequences observees,
 *              sur l'histogramme des longueurs des sequences et
 *              sur les lois de l'intervalle de temps residuel.
 *
 *--------------------------------------------------------------*/

bool Nonparametric_sequence_process::plot_print(const char *prefix , const char *title ,
                                                int process , Histogram **empirical_observation ,
                                                const Sequence_characteristics *characteristics ,
                                                const Histogram *hlength , Forward **forward) const

{
  bool status = false , start;
  register int i , j , k , m;
  int index_length , nb_histo , nb_dist , histo_index , dist_index , *dist_nb_value ,
      *index_dist;
  double *scale;
  Curves *smoothed_curves;
  const Distribution **pdist;
  const Histogram **phisto;
  ostringstream data_file_name[2];


  // ecriture des fichiers de donnees

  if ((index_value) || (characteristics)) {
    data_file_name[0] << prefix << process << 0 << ".dat";

    if (characteristics) {
      index_length = characteristics->index_value->plot_length_computation();
      if (characteristics->index_value->frequency[index_length - 1] < MAX_FREQUENCY) {
        smoothed_curves = new Curves(*(characteristics->index_value) , 's');
      }
      else {
        smoothed_curves = 0;
      }
    }

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
    data_file_name[1] << prefix << process << 1 << ".dat";

    pdist = new const Distribution*[6 * nb_value + nb_state];
    dist_nb_value = new int[6 * nb_value + nb_state];
    scale = new double[6 * nb_value + nb_state];
    phisto = new const Histogram*[1 + 7 * nb_value + nb_state];
    index_dist = new int[1 + 7 * nb_value + nb_state];

    nb_histo = 0;
    nb_dist = 0;

    if (hlength) {
      phisto[nb_histo] = hlength;
      index_dist[nb_histo++] = I_DEFAULT;
    }

    if ((first_occurrence) || (characteristics)) {
      for (i = 0;i < nb_value;i++) {
        if ((first_occurrence) && (first_occurrence[i])) {
          pdist[nb_dist] = first_occurrence[i];

          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->first_occurrence[i]->nb_element > 0)) {
            phisto[nb_histo] = characteristics->first_occurrence[i];
            index_dist[nb_histo] = nb_dist;
            dist_nb_value[nb_dist] = MIN(first_occurrence[i]->nb_value , phisto[nb_histo]->nb_value * 3);
            scale[nb_dist++] = phisto[nb_histo++]->nb_element /
                               first_occurrence[i]->cumul[first_occurrence[i]->nb_value - 1];
          }
          else {
            dist_nb_value[nb_dist] = first_occurrence[i]->nb_value;
            scale[nb_dist++] = 1.;
          }
        }

        else {
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->first_occurrence[i]->nb_element > 0)) {
            phisto[nb_histo] = characteristics->first_occurrence[i];
            index_dist[nb_histo++] = I_DEFAULT;
          }
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
            index_dist[nb_histo] = nb_dist;
            dist_nb_value[nb_dist] = MIN(recurrence_time[i]->nb_value , phisto[nb_histo]->nb_value * 3);
            scale[nb_dist++] = phisto[nb_histo++]->nb_element / (1. - leave[i]);
          }
          else {
            dist_nb_value[nb_dist] = recurrence_time[i]->nb_value;
            scale[nb_dist++] = 1.;
          }
        }

        else {
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->recurrence_time[i]->nb_element > 0)) {
            phisto[nb_histo] = characteristics->recurrence_time[i];
            index_dist[nb_histo++] = I_DEFAULT;
          }
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
            index_dist[nb_histo] = nb_dist;
            if (sojourn_time[i]->cumul[sojourn_time[i]->nb_value - 1] < CUMUL_THRESHOLD) {
              scale[nb_dist++] = phisto[nb_histo++]->nb_element / sojourn_time[i]->cumul[sojourn_time[i]->nb_value - 1];
            }
            else {
              scale[nb_dist++] = phisto[nb_histo++]->nb_element;
            }
          }
          else {
            scale[nb_dist++] = 1.;
          }
        }

        else {
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->sojourn_time[i]->nb_element > 0)) {
            phisto[nb_histo] = characteristics->sojourn_time[i];
            index_dist[nb_histo++] = I_DEFAULT;
          }
        }

        if ((characteristics) && (i < characteristics->nb_value) &&
            (characteristics->initial_run) &&
            (characteristics->initial_run[i]->nb_element > 0)) {
          if ((forward) && (forward[i])) {
            pdist[nb_dist] = forward[i];
            dist_nb_value[nb_dist] = forward[i]->nb_value;
            phisto[nb_histo] = characteristics->initial_run[i];
            index_dist[nb_histo] = nb_dist;
            scale[nb_dist++] = phisto[nb_histo++]->nb_element;
          }

          else {
            phisto[nb_histo] = characteristics->initial_run[i];
            index_dist[nb_histo++] = I_DEFAULT;
          }
        }

        if ((forward) && (forward[i])) {
          pdist[nb_dist] = forward[i];
          dist_nb_value[nb_dist] = forward[i]->nb_value;

          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->final_run[i]->nb_element > 0)) {
            phisto[nb_histo] = characteristics->final_run[i];
            index_dist[nb_histo] = nb_dist;
            scale[nb_dist++] = phisto[nb_histo++]->nb_element;
          }
          else {
            scale[nb_dist++] = 1.;
          }
        }

        else {
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->final_run[i]->nb_element > 0)) {
            phisto[nb_histo] = characteristics->final_run[i];
            index_dist[nb_histo++] = I_DEFAULT;
          }
        }
      }
    }

    if ((nb_run) || (nb_occurrence) ||
        ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence))) {
      for (i = 0;i < nb_value;i++) {
        if (nb_run) {
          pdist[nb_dist] = nb_run[i];

          if ((characteristics) && (i < characteristics->nb_value) && (characteristics->nb_run) &&
              (characteristics->nb_run[i]->nb_element > 0)) {
            phisto[nb_histo] = characteristics->nb_run[i];
            index_dist[nb_histo] = nb_dist;
            dist_nb_value[nb_dist] = nb_run[i]->plot_nb_value_computation(phisto[nb_histo]);
            scale[nb_dist++] = phisto[nb_histo++]->nb_element;
          }
          else {
            dist_nb_value[nb_dist] = nb_run[i]->plot_nb_value_computation();
            scale[nb_dist++] = 1.;
          }
        }

        else {
          if ((characteristics) && (i < characteristics->nb_value) && (characteristics->nb_run) &&
              (characteristics->nb_run[i]->nb_element > 0)) {
            phisto[nb_histo] = characteristics->nb_run[i];
            index_dist[nb_histo++] = I_DEFAULT;
          }
        }

        if (nb_occurrence) {
          pdist[nb_dist] = nb_occurrence[i];

          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->nb_occurrence) &&
              (characteristics->nb_occurrence[i]->nb_element > 0)) {
            phisto[nb_histo] = characteristics->nb_occurrence[i];
            index_dist[nb_histo] = nb_dist;
            dist_nb_value[nb_dist] = nb_occurrence[i]->plot_nb_value_computation(phisto[nb_histo]);
            scale[nb_dist++] = phisto[nb_histo++]->nb_element;
          }
          else {
            dist_nb_value[nb_dist] = nb_occurrence[i]->plot_nb_value_computation();
            scale[nb_dist++] = 1.;
          }
        }

        else {
          if ((characteristics) && (i < characteristics->nb_value) &&
              (characteristics->nb_occurrence) &&
              (characteristics->nb_occurrence[i]->nb_element > 0)) {
            phisto[nb_histo] = characteristics->nb_occurrence[i];
            index_dist[nb_histo++] = I_DEFAULT;
          }
        }
      }
    }

    if (observation) {
      for (i = 0;i < nb_state;i++) {
        pdist[nb_dist] = observation[i];
        dist_nb_value[nb_dist] = observation[i]->nb_value;

        if ((empirical_observation) && (empirical_observation[i]->nb_element > 0)) {
          phisto[nb_histo] = empirical_observation[i];
          index_dist[nb_histo] = nb_dist;
          scale[nb_dist++] = phisto[nb_histo++]->nb_element;
        }
        else {
          scale[nb_dist++] = 1.;
        }
      }
    }

    status = ::plot_print((data_file_name[1].str()).c_str() , nb_dist , pdist , scale ,
                          dist_nb_value , nb_histo , phisto , index_dist);

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

          out_file << "set title \"";
          if (title) {
            out_file << title;
            if (process > 0) {
              out_file << " - ";
            }
          }
          if (process > 0) {
            out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
          }
          out_file << "\"\n\n";

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

          if (hlength->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }
          if ((int)(hlength->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [0:" << hlength->nb_value - 1 << "] [0:"
                   << (int)(hlength->max * YSCALE) + 1 << "] \""
                   << label((data_file_name[1].str()).c_str()) << "\" using 1 title \""
                   << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM]
                   << "\" with impulses" << endl;

          if (hlength->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
          if ((int)(hlength->max * YSCALE) + 1 < TIC_THRESHOLD) {
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
                         << " title \"" << SEQ_label[process == 0 ? SEQL_STATE_FIRST_OCCURRENCE : SEQL_OUTPUT_FIRST_OCCURRENCE]
                         << " " << m << " " << STAT_label[STATL_HISTOGRAM] << "\" with impulses,\\" << endl;
                out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1
                         << " title \"" << SEQ_label[process == 0 ? SEQL_STATE_FIRST_OCCURRENCE : SEQL_OUTPUT_FIRST_OCCURRENCE]
                         << " " << m << " " << STAT_label[STATL_DISTRIBUTION] << "\" with linespoints" << endl;
                j++;
              }

              else {
                out_file << "plot [0:" << MAX(dist_nb_value[k] , 2) - 1 << "] [0:"
                         << MIN(pdist[k]->max * YSCALE , 1.) << "] \""
                         << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1
                         << " title \"" << SEQ_label[process == 0 ? SEQL_STATE_FIRST_OCCURRENCE : SEQL_OUTPUT_FIRST_OCCURRENCE]
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
                       << " title \"" << SEQ_label[process == 0 ? SEQL_STATE_FIRST_OCCURRENCE : SEQL_OUTPUT_FIRST_OCCURRENCE]
                       << " " << m << " " << STAT_label[STATL_HISTOGRAM] << "\" with impulses" << endl;

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
                         << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_HISTOGRAM]
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
                       << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_HISTOGRAM]
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
                         << " " << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_HISTOGRAM]
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
                       << " " << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_HISTOGRAM]
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
                         << " " << m << " " << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_HISTOGRAM]
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
                         << " " << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_HISTOGRAM]
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
                         << " " << m << " " << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_HISTOGRAM]
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
                       << " " << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_HISTOGRAM]
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
                         << " title \"" << SEQ_label[process == 0 ? SEQL_STATE_NB_RUN : SEQL_OUTPUT_NB_RUN]
                         << " " << m << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_HISTOGRAM]
                         << "\" with impulses,\\" << endl;
                out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1;
                if (length->variance == 0.) {
                  out_file << " title \"" << SEQ_label[process == 0 ? SEQL_STATE_NB_RUN : SEQL_OUTPUT_NB_RUN]
                           << " " << m << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
                           << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION];
                }
                else {
                  out_file << " title \"" << SEQ_label[SEQL_MIXTURE_OF]
                           << SEQ_label[process == 0 ? SEQL_STATE_NB_RUN : SEQL_OUTPUT_NB_RUN]
                           << " " << m << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION];
                }
                out_file << "\" with linespoints" << endl;
                j++;
              }

              else {
                out_file << "plot [0:" << dist_nb_value[k] - 1 << "] [0:"
                         << MIN(pdist[k]->max * YSCALE , 1.) << "] \""
                         << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1
                         << " title \"" << SEQ_label[process == 0 ? SEQL_STATE_NB_RUN : SEQL_OUTPUT_NB_RUN]
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
                       << " title \"" << SEQ_label[process == 0 ? SEQL_STATE_NB_RUN : SEQL_OUTPUT_NB_RUN]
                       << " " << m << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_HISTOGRAM]
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
                         << " title \"" << SEQ_label[process == 0 ? SEQL_STATE_NB_OCCURRENCE : SEQL_OUTPUT_NB_OCCURRENCE]
                         << " " << m << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_HISTOGRAM]
                         << "\" with impulses,\\" << endl;
                out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1;
                if (length->variance == 0.) {
                  out_file << " title \"" << SEQ_label[process == 0 ? SEQL_STATE_NB_OCCURRENCE : SEQL_OUTPUT_NB_OCCURRENCE]
                           << " " << m << " " << SEQ_label[SEQL_PER_LENGTH] << " " << length->offset << " "
                           << SEQ_label[SEQL_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION];
                }
                else {
                  out_file << " title \"" << SEQ_label[SEQL_MIXTURE_OF]
                           << SEQ_label[process == 0 ? SEQL_STATE_NB_OCCURRENCE : SEQL_OUTPUT_NB_OCCURRENCE]
                           << " " << m << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_DISTRIBUTION];
                }
                out_file << "\" with linespoints" << endl;
                j++;
              }

              else {
                out_file << "plot [0:" << dist_nb_value[k] - 1 << "] [0:"
                         << MIN(pdist[k]->max * YSCALE , 1.) << "] \""
                         << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo + k + 1
                         << " title \"" << SEQ_label[process == 0 ? SEQL_STATE_NB_OCCURRENCE : SEQL_OUTPUT_NB_OCCURRENCE]
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
                       << " title \"" << SEQ_label[process == 0 ? SEQL_STATE_NB_OCCURRENCE : SEQL_OUTPUT_NB_OCCURRENCE]
                       << " " << m << " " << SEQ_label[SEQL_PER_SEQUENCE] << " " << STAT_label[STATL_HISTOGRAM]
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

            if (hlength->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
            }
            if ((int)(hlength->max * YSCALE) + 1 < TIC_THRESHOLD) {
              out_file << "set ytics 0,1" << endl;
            }

            out_file << "plot [0:" << hlength->nb_value - 1 << "] [0:"
                     << (int)(hlength->max * YSCALE) + 1 << "] \""
                     << label((data_file_name[1].str()).c_str()) << "\" using 1 title \""
                     << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM]
                     << "\" with impulses" << endl;

            if (hlength->nb_value - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
            if ((int)(hlength->max * YSCALE) + 1 < TIC_THRESHOLD) {
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
                       << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_HISTOGRAM]
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
            out_file << endl;
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
      delete [] index_dist;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWspace Nonparametric_sequence_process::binaryStoreSize() const

{
  register int i;
  RWspace size;


  size = sizeof(nb_state) + sizeof(nb_value);

  size += sizeof(true);
  if (observation) {
    for (i = 0;i < nb_state;i++) {
      size += observation[i]->binaryStoreSize();
    }
  }

  size += sizeof(true);
  if (length) {
    size += length->binaryStoreSize();
  }

  size += sizeof(true);
  if (index_value) {
    size += index_value->binaryStoreSize();
  }

  size += sizeof(true);
  if (no_occurrence) {
    size += sizeof(*no_occurrence) * nb_value;
  }

  size += sizeof(true);
  if (first_occurrence) {
    for (i = 0;i < nb_value;i++) {
      size += first_occurrence[i]->binaryStoreSize();
    }
  }

  size += sizeof(true);
  if (leave) {
    size += sizeof(*leave) * nb_value;
  }

  size += sizeof(true);
  if (recurrence_time) {
    size += sizeof(true) * nb_value;
    for (i = 0;i < nb_value;i++) {
      if (recurrence_time[i]) {
        size += recurrence_time[i]->binaryStoreSize();
      }
    }
  }

  size += sizeof(true);
  if (absorption) {
    size += sizeof(*absorption) * nb_value;
  }

  size += sizeof(true);
  if (sojourn_time) {
    size += sizeof(true) * nb_value;
    for (i = 0;i < nb_value;i++) {
      if (sojourn_time[i]) {
        size += sojourn_time[i]->binaryStoreSize();
      }
    }
  }

  size += sizeof(true);
  if (nb_run) {
    for (i = 0;i < nb_value;i++) {
      size += nb_run[i]->binaryStoreSize();
    }
  }

  size += sizeof(true);
  if (nb_occurrence) {
    for (i = 0;i < nb_value;i++) {
      size += nb_occurrence[i]->binaryStoreSize();
    }
  }

  return size;
}


void Nonparametric_sequence_process::restoreGuts(RWvistream &is)

{
  bool status;
  register int i;


  remove();

  is >> nb_state >> nb_value;

  is >> status;
  if (status) {
    observation = new Distribution*[nb_state];
    for (i = 0;i < nb_state;i++) {
      observation[i] = new Distribution();
      observation[i]->restoreGuts(is);
    }
  }
  else {
    observation = 0;
  }

  is >> status;
  if (status) {
    length = new Distribution();
    length->restoreGuts(is);
  }
  else {
    length = 0;
  }

  is >> status;
  if (status) {
    index_value = new Curves();
    index_value->restoreGuts(is);
  }
  else {
    index_value = 0;
  }

  is >> status;
  if (status) {
    no_occurrence = new double[nb_value];
    for (i = 0;i < nb_value;i++) {
      is >> no_occurrence[i];
    }
  }
  else {
    no_occurrence = 0;
  }

  is >> status;
  if (status) {
    first_occurrence = new Distribution*[nb_value];
    for (i = 0;i < nb_value;i++) {
      first_occurrence[i] = new Distribution();
      first_occurrence[i]->restoreGuts(is);
    }
  }
  else {
    first_occurrence = 0;
  }

  is >> status;
  if (status) {
    leave = new double[nb_value];
    for (i = 0;i < nb_value;i++) {
      is >> leave[i];
    }
  }
  else {
    leave = 0;
  }

  is >> status;
  if (status) {
    recurrence_time = new Distribution*[nb_value];
    for (i = 0;i < nb_value;i++) {
      is >> status;
      if (status) {
        recurrence_time[i] = new Distribution();
        recurrence_time[i]->restoreGuts(is);
      }
      else {
        recurrence_time[i] = 0;
      }
    }
  }
  else {
    recurrence_time = 0;
  }

  is >> status;
  if (status) {
    absorption = new double[nb_value];
    for (i = 0;i < nb_value;i++) {
      is >> absorption[i];
    }
  }
  else {
    absorption = 0;
  }

  is >> status;
  if (status) {
    sojourn_time = new Parametric*[nb_value];
    for (i = 0;i < nb_value;i++) {
      is >> status;
      if (status) {
        sojourn_time[i] = new Parametric();
        sojourn_time[i]->restoreGuts(is);
      }
      else {
        sojourn_time[i] = 0;
      }
    }
  }
  else {
    sojourn_time = 0;
  }

  is >> status;
  if (status) {
    nb_run = new Distribution*[nb_value];
    for (i = 0;i < nb_value;i++) {
      nb_run[i] = new Distribution();
      nb_run[i]->restoreGuts(is);
    }
  }
  else {
    nb_run = 0;
  }

  is >> status;
  if (status) {
    nb_occurrence = new Distribution*[nb_value];
    for (i = 0;i < nb_value;i++) {
      nb_occurrence[i] = new Distribution();
      nb_occurrence[i]->restoreGuts(is);
    }
  }
  else {
    nb_occurrence = 0;
  }
}


void Nonparametric_sequence_process::restoreGuts(RWFile &file)

{
  bool status;
  register int i;


  remove();

  file.Read(nb_state);
  file.Read(nb_value);

  file.Read(status);
  if (status) {
    observation = new Distribution*[nb_state];
    for (i = 0;i < nb_state;i++) {
      observation[i] = new Distribution();
      observation[i]->restoreGuts(file);
    }
  }
  else {
    observation = 0;
  }

  file.Read(status);
  if (status) {
    length = new Distribution();
    length->restoreGuts(file);
  }
  else {
    length = 0;
  }

  file.Read(status);
  if (status) {
    index_value = new Curves();
    index_value->restoreGuts(file);
  }
  else {
    index_value = 0;
  }

  file.Read(status);
  if (status) {
    no_occurrence = new double[nb_value];
    file.Read(no_occurrence , nb_value);
  }
  else {
    no_occurrence = 0;
  }

  file.Read(status);
  if (status) {
    first_occurrence = new Distribution*[nb_value];
    for (i = 0;i < nb_value;i++) {
      first_occurrence[i] = new Distribution();
      first_occurrence[i]->restoreGuts(file);
    }
  }
  else {
    first_occurrence = 0;
  }

  file.Read(status);
  if (status) {
    leave = new double[nb_value];
    file.Read(leave , nb_value);
  }
  else {
    leave = 0;
  }

  file.Read(status);
  if (status) {
    recurrence_time = new Distribution*[nb_value];
    for (i = 0;i < nb_value;i++) {
      file.Read(status);
      if (status) {
        recurrence_time[i] = new Distribution();
        recurrence_time[i]->restoreGuts(file);
      }
      else {
        recurrence_time[i] = 0;
      }
    }
  }
  else {
    recurrence_time = 0;
  }

  file.Read(status);
  if (status) {
    absorption = new double[nb_value];
    file.Read(absorption , nb_value);
  }
  else {
    absorption = 0;
  }

  file.Read(status);
  if (status) {
    sojourn_time = new Parametric*[nb_value];
    for (i = 0;i < nb_value;i++) {
      file.Read(status);
      if (status) {
        sojourn_time[i] = new Parametric();
        sojourn_time[i]->restoreGuts(file);
      }
      else {
        sojourn_time[i] = 0;
      }
    }
  }
  else {
    sojourn_time = 0;
  }

  file.Read(status);
  if (status) {
    nb_run = new Distribution*[nb_value];
    for (i = 0;i < nb_value;i++) {
      nb_run[i] = new Distribution();
      nb_run[i]->restoreGuts(file);
    }
  }
  else {
    nb_run = 0;
  }

  file.Read(status);
  if (status) {
    nb_occurrence = new Distribution*[nb_value];
    for (i = 0;i < nb_value;i++) {
      nb_occurrence[i] = new Distribution();
      nb_occurrence[i]->restoreGuts(file);
    }
  }
  else {
    nb_occurrence = 0;
  }
}


void Nonparametric_sequence_process::saveGuts(RWvostream &os) const

{
  register int i;


  os << nb_state << nb_value;

  if (observation) {
    os << true;
    for (i = 0;i < nb_state;i++) {
      observation[i]->saveGuts(os);
    }
  }
  else {
    os << false;
  }

  if (length) {
    os << true;
    length->saveGuts(os);
  }
  else {
    os << false;
  }

  if (index_value) {
    os << true;
    index_value->saveGuts(os);
  }
  else {
    os << false;
  }

  if (no_occurrence) {
    os << true;
    for (i = 0;i < nb_value;i++) {
      os << no_occurrence[i];
    }
  }
  else {
    os << false;
  }

  if (first_occurrence) {
    os << true;
    for (i = 0;i < nb_value;i++) {
      first_occurrence[i]->saveGuts(os);
    }
  }
  else {
    os << false;
  }

  if (leave) {
    os << true;
    for (i = 0;i < nb_value;i++) {
      os << leave[i];
    }
  }
  else {
    os << false;
  }

  if (recurrence_time) {
    os << true;
    for (i = 0;i < nb_value;i++) {
      if (recurrence_time[i]) {
        os << true;
        recurrence_time[i]->saveGuts(os);
      }
      else {
        os << false;
      }
    }
  }
  else {
    os << false;
  }

  if (absorption) {
    os << true;
    for (i = 0;i < nb_value;i++) {
      os << absorption[i];
    }
  }
  else {
    os << false;
  }

  if (sojourn_time) {
    os << true;
    for (i = 0;i < nb_value;i++) {
      if (sojourn_time[i]) {
        os << true;
        sojourn_time[i]->saveGuts(os);
      }
      else {
        os << false;
      }
    }
  }
  else {
    os << false;
  }

  if (nb_run) {
    os << true;
    for (i = 0;i < nb_value;i++) {
      nb_run[i]->saveGuts(os);
    }
  }
  else {
    os << false;
  }

  if (nb_occurrence) {
    os << true;
    for (i = 0;i < nb_value;i++) {
      nb_occurrence[i]->saveGuts(os);
    }
  }
  else {
    os << false;
  }
}


void Nonparametric_sequence_process::saveGuts(RWFile &file) const

{
  register int i;


  file.Write(nb_state);
  file.Write(nb_value);

  if (observation) {
    file.Write(true);
    for (i = 0;i < nb_state;i++) {
      observation[i]->saveGuts(file);
    }
  }
  else {
    file.Write(false);
  }

  if (length) {
    file.Write(true);
    length->saveGuts(file);
  }
  else {
    file.Write(false);
  }

  if (index_value) {
    file.Write(true);
    index_value->saveGuts(file);
  }
  else {
    file.Write(false);
  }

  if (no_occurrence) {
    file.Write(true);
    file.Write(no_occurrence , nb_value);
  }
  else {
    file.Write(false);
  }

  if (first_occurrence) {
    file.Write(true);
    for (i = 0;i < nb_value;i++) {
      first_occurrence[i]->saveGuts(file);
    }
  }
  else {
    file.Write(false);
  }

  if (leave) {
    file.Write(true);
    file.Write(leave , nb_value);
  }
  else {
    file.Write(false);
  }

  if (recurrence_time) {
    file.Write(true);
    for (i = 0;i < nb_value;i++) {
      if (recurrence_time[i]) {
        file.Write(true);
        recurrence_time[i]->saveGuts(file);
      }
      else {
        file.Write(false);
      }
    }
  }
  else {
    file.Write(false);
  }

  if (absorption) {
    file.Write(true);
    file.Write(absorption , nb_value);
  }
  else {
    file.Write(false);
  }

  if (sojourn_time) {
    file.Write(true);
    for (i = 0;i < nb_value;i++) {
      if (sojourn_time[i]) {
        file.Write(true);
        sojourn_time[i]->saveGuts(file);
      }
      else {
        file.Write(false);
      }
    }
  }
  else {
    file.Write(false);
  }

  if (nb_run) {
    file.Write(true);
    for (i = 0;i < nb_value;i++) {
      nb_run[i]->saveGuts(file);
    }
  }
  else {
    file.Write(false);
  }

  if (nb_occurrence) {
    file.Write(true);
    for (i = 0;i < nb_value;i++) {
      nb_occurrence[i]->saveGuts(file);
    }
  }
  else {
    file.Write(false);
  }
} */
