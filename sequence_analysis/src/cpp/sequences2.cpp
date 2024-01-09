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
 *       $Id: sequences2.cpp 18074 2015-04-23 10:56:21Z guedon $
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



#include <limits.h>
#include <math.h>
#include <sstream>
#include <iomanip>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>

#include "tool/config.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

#include "sequences.h"
#include "sequence_label.h"

using namespace std;
using namespace boost::math;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Fusion d'objets Sequences.
 *
 *  arguments : reference sur un objet StatError, nombre d'objets Sequences,
 *              pointeurs sur les objets Sequences.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::merge(StatError &error , int nb_sample , const Sequences **iseq) const

{
  bool status = true;
  register int i , j , k , m , n , p , q;
  int inb_sequence , cumul_nb_sequence , *ilength , *iidentifier , **ivertex_identifier;
  const FrequencyDistribution **phisto;
  Sequences *seq;
  const Sequences **pseq;


  seq = NULL;
  error.init();

  for (i = 0;i < nb_sample;i++) {
    if (iseq[i]->index_parameter_type != index_parameter_type) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_INDEX_PARAMETER_TYPE];

      if (index_parameter_type == IMPLICIT_TYPE) {
        error.update((error_message.str()).c_str());
      }
      else {
        error.correction_update((error_message.str()).c_str() , SEQ_index_parameter_word[index_parameter_type]);
      }
    }
  }

  for (i = 0;i < nb_sample;i++) {
    if (iseq[i]->nb_variable != nb_variable) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << STAT_error[STATR_NB_VARIABLE];
      error.correction_update((error_message.str()).c_str() , nb_variable);
    }

    else {
      for (j = 0;j < nb_variable;j++) {
        if (iseq[i]->type[j] != type[j]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                        << STAT_label[STATL_VARIABLE] << " " << j + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          error.correction_update((error_message.str()).c_str() , STAT_variable_word[type[j]]);
        }
      }
    }
  }

  if (status) {
    nb_sample++;
    pseq = new const Sequences*[nb_sample];

    pseq[0] = this;
    for (i = 1;i < nb_sample;i++) {
      pseq[i] = iseq[i - 1];
    }

    // calcul du nombre de sequences

    inb_sequence = 0;
    for (i = 0;i < nb_sample;i++) {
      inb_sequence += pseq[i]->nb_sequence;
    }

    // comparaison des identificateurs des sequences

    iidentifier = new int[inb_sequence];

    cumul_nb_sequence = 0;
    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < pseq[j]->nb_sequence;k++) {
        iidentifier[i] = pseq[j]->identifier[k];

        for (m = 0;m < cumul_nb_sequence;m++) {
          if (iidentifier[i] == iidentifier[m]) {
            delete [] iidentifier;
            iidentifier = NULL;
            break;
          }
        }

        if (!iidentifier) {
          break;
        }
        i++;
      }

      if (!iidentifier) {
        break;
      }
      cumul_nb_sequence += pseq[j]->nb_sequence;
    }

    // copie des longueurs des sequences

    ilength = new int[inb_sequence];

    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < pseq[j]->nb_sequence;k++) {
        ilength[i++] = pseq[j]->length[k];
      }
    }

    // comparaison des identificateurs des vertex

    for (i = 0;i < nb_sample;i++) {
      if (!(pseq[i]->vertex_identifier)) {
        break;
      }
    }

    if (i == nb_sample) {
      ivertex_identifier = new int*[inb_sequence];

      cumul_nb_sequence = 0;
      i = 0;
      for (j = 0;j < nb_sample;j++) {
        for (k = 0;k < pseq[j]->nb_sequence;k++) {
          ivertex_identifier[i] = new int[pseq[j]->length[k]];
          for (m = 0;m < pseq[j]->length[k];m++) {
            ivertex_identifier[i][m] = pseq[j]->vertex_identifier[k][m];

            for (n = 0;n < cumul_nb_sequence;n++) {
              for (p = 0;p < ilength[n];p++) {
                if (ivertex_identifier[i][m] == ivertex_identifier[n][p]) {
                  for (q = 0;q <= i;q++) {
                    delete [] ivertex_identifier[q];
                  }
                  delete [] ivertex_identifier;
                  ivertex_identifier = NULL;
                  break;
                }
              }

              if (!ivertex_identifier) {
                break;
              }
            }

            if (!ivertex_identifier) {
              break;
            }
          }

          if (!ivertex_identifier) {
            break;
          }
          i++;
        }

        if (!ivertex_identifier) {
          break;
        }
        cumul_nb_sequence += pseq[j]->nb_sequence;
      }
    }

    else {
      ivertex_identifier = NULL;
    }

    seq = new Sequences(inb_sequence , iidentifier , ilength , ivertex_identifier ,
                        index_parameter_type , nb_variable , type);
    delete [] iidentifier;
    delete [] ilength;

    if (ivertex_identifier) {
      for (i = 0;i < inb_sequence;i++) {
        delete [] ivertex_identifier[i];
      }
      delete [] ivertex_identifier;
    }

    phisto = new const FrequencyDistribution*[nb_sample];

    // copy of index parameters

    if (seq->index_parameter) {
      i = 0;
      for (j = 0;j < nb_sample;j++) {
        for (k = 0;k < pseq[j]->nb_sequence;k++) {
          for (m = 0;m < (pseq[j]->index_parameter_type == POSITION ? pseq[j]->length[k] + 1 : pseq[j]->length[k]);m++) {
            seq->index_parameter[i][m] = pseq[j]->index_parameter[k][m];
          }
          i++;
        }
      }

      for (i = 0;i < nb_sample;i++) {
        phisto[i] = pseq[i]->index_parameter_distribution;
      }
      seq->index_parameter_distribution = new FrequencyDistribution(nb_sample , phisto);
    }

//    if ((seq->index_parameter_type == TIME) || ((seq->index_parameter_type == POSITION) &&
//         (seq->type[0] != NB_INTERNODE))) {
    if ((seq->index_parameter_type == TIME) || (seq->index_parameter_type == POSITION)) {
      for (i = 0;i < nb_sample;i++) {
        phisto[i] = pseq[i]->index_interval;
      }
      seq->index_interval = new FrequencyDistribution(nb_sample , phisto);
    }

    // copie des valeurs

    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < pseq[j]->nb_sequence;k++) {
        for (m = 0;m < pseq[j]->nb_variable;m++) {
          if ((pseq[j]->type[m] != REAL_VALUE) && (pseq[j]->type[m] != AUXILIARY)) {
            for (n = 0;n < pseq[j]->length[k];n++) {
              seq->int_sequence[i][m][n] = pseq[j]->int_sequence[k][m][n];
            }
          }

          else {
            for (n = 0;n < pseq[j]->length[k];n++) {
              seq->real_sequence[i][m][n] = pseq[j]->real_sequence[k][m][n];
            }
          }
        }
        i++;
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value[i] = pseq[0]->min_value[i];
      seq->max_value[i] = pseq[0]->max_value[i];
      for (j = 1;j < nb_sample;j++) {
        if (pseq[j]->min_value[i] < seq->min_value[i]) {
          seq->min_value[i] = pseq[j]->min_value[i];
        }
        if (pseq[j]->max_value[i] > seq->max_value[i]) {
          seq->max_value[i] = pseq[j]->max_value[i];
        }
      }

      if (seq->type[i] != AUXILIARY) {
        for (j = 0;j < nb_sample;j++) {
          if (pseq[j]->marginal_distribution[i]) {
            phisto[j] = pseq[j]->marginal_distribution[i];
          }
          else {
            break;
          }
        }

        if (j == nb_sample) {
          seq->marginal_distribution[i] = new FrequencyDistribution(nb_sample , phisto);
        }

        else {
          seq->build_marginal_histogram(i);
        }
      }
    }

    delete [] pseq;
    delete [] phisto;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Translation des valeurs d'une variable.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              parametre de translation.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::shift(StatError &error , int variable , int shift_param) const

{
  bool status = true;
  register int i , j;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (type[variable] == INT_VALUE) {
      if (shift_param + min_value[variable] < INT_MIN) {
        status = false;
        ostringstream correction_message;
        correction_message << STAT_error[STATR_GREATER_THAN] << " "
                           << INT_MIN - min_value[variable];
        error.correction_update(STAT_error[STATR_SHIFT_VALUE] , (correction_message.str()).c_str());
      }

      if (shift_param + max_value[variable] > INT_MAX) {
        status = false;
        ostringstream correction_message;
        correction_message << STAT_error[STATR_SMALLER_THAN] << " "
                           << INT_MAX - max_value[variable];
        error.correction_update(STAT_error[STATR_SHIFT_VALUE] , (correction_message.str()).c_str());
      }
    }

    else if (type[variable] != REAL_VALUE) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[REAL_VALUE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }
  }

  if (status) {
    seq = new Sequences(*this , variable , type[variable]);

    switch (seq->type[variable]) {

    // translation des valeurs entieres

    case INT_VALUE : {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->length[i];j++) {
          seq->int_sequence[i][variable][j] = int_sequence[i][variable][j] + shift_param;
        }
      }
      break;
    }

    // translation des valeurs reelles

    case REAL_VALUE : {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->length[i];j++) {
          seq->real_sequence[i][variable][j] = real_sequence[i][variable][j] + shift_param;
        }
      }
      break;
    }
    }

    seq->min_value[variable] = min_value[variable] + shift_param;
    seq->max_value[variable] = max_value[variable] + shift_param;

    if ((variable + 1 < seq->nb_variable) && (seq->type[variable + 1] == AUXILIARY)) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->length[i];j++) {
          seq->real_sequence[i][variable + 1][j] = real_sequence[i][variable + 1][j] + shift_param;
        }
      }

      seq->min_value[variable + 1] = min_value[variable + 1] + shift_param;
      seq->max_value[variable + 1] = max_value[variable + 1] + shift_param;
    }

    if ((seq->type[variable] == INT_VALUE) && (seq->min_value[variable] >= 0) &&
        (seq->max_value[variable] <= MARGINAL_DISTRIBUTION_MAX_VALUE)) {
      if (marginal_distribution[variable]) {
        seq->marginal_distribution[variable] = new FrequencyDistribution(*marginal_distribution[variable] , 's' , shift_param);
      }
      else {
        seq->build_marginal_frequency_distribution(variable);
      }
    }

    else {
      seq->build_marginal_histogram(variable);
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Translation des valeurs d'une variable reelle.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              parametre de translation.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::shift(StatError &error , int variable , double shift_param) const

{
  bool status = true;
  register int i , j;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (type[variable] != REAL_VALUE) {
      status = false;
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , STAT_variable_word[REAL_VALUE]);
    }
  }

  if (status) {
    seq = new Sequences(*this , variable , type[variable]);

    // translation des valeurs reelles

    for (i = 0;i < seq->nb_sequence;i++) {
      for (j = 0;j < seq->length[i];j++) {
        seq->real_sequence[i][variable][j] = real_sequence[i][variable][j] + shift_param;
      }
    }

    seq->min_value[variable] = min_value[variable] + shift_param;
    seq->max_value[variable] = max_value[variable] + shift_param;

    if ((variable + 1 < seq->nb_variable) && (seq->type[variable + 1] == AUXILIARY)) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->length[i];j++) {
          seq->real_sequence[i][variable + 1][j] = real_sequence[i][variable + 1][j] + shift_param;
        }
      }

      seq->min_value[variable + 1] = min_value[variable + 1] + shift_param;
      seq->max_value[variable + 1] = max_value[variable + 1] + shift_param;
    }

    seq->build_marginal_histogram(variable);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Seuillage des valeurs d'une variable.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              seuil, mode (ABOVE/BELOW).
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::thresholding(StatError &error , int variable , int threshold , int mode) const

{
  bool status = true;
  register int i , j;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != REAL_VALUE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[REAL_VALUE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    if ((type[variable] == INT_VALUE) && (variable + 1 < nb_variable) &&
        (type[variable + 1] == AUXILIARY)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      error.update((error_message.str()).c_str());
    }

    if (threshold <= min_value[variable]) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_error[STATR_GREATER_THAN] << " " << min_value[variable];
      error.correction_update(STAT_error[STATR_THRESHOLD_VALUE] , (correction_message.str()).c_str());
    }

    if (threshold >= max_value[variable]) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_error[STATR_SMALLER_THAN] << " " << max_value[variable];
      error.correction_update(STAT_error[STATR_THRESHOLD_VALUE] , (correction_message.str()).c_str());
    }
  }

  if (status) {
    seq = new Sequences(*this , variable , type[variable]);

    switch (seq->type[variable]) {

    // seuillage des valeurs entieres

    case INT_VALUE : {
      switch (mode) {

      case ABOVE : {
        for (i = 0;i < seq->nb_sequence;i++) {
          for (j = 0;j < seq->length[i];j++) {
            if (int_sequence[i][variable][j] > threshold) {
              seq->int_sequence[i][variable][j] = threshold;
            }
            else {
              seq->int_sequence[i][variable][j] = int_sequence[i][variable][j];
            }
          }
        }
        break;
      }

      case BELOW : {
        for (i = 0;i < seq->nb_sequence;i++) {
          for (j = 0;j < seq->length[i];j++) {
            if (int_sequence[i][variable][j] < threshold) {
              seq->int_sequence[i][variable][j] = threshold;
            }
            else {
              seq->int_sequence[i][variable][j] = int_sequence[i][variable][j];
            }
          }
        }
        break;
      }
      }

      break;
    }

    // seuillage des valeurs reelles

    case REAL_VALUE : {
      switch (mode) {

      case ABOVE : {
        for (i = 0;i < seq->nb_sequence;i++) {
          for (j = 0;j < seq->length[i];j++) {
            if (real_sequence[i][variable][j] > threshold) {
              seq->real_sequence[i][variable][j] = threshold;
            }
            else {
              seq->real_sequence[i][variable][j] = real_sequence[i][variable][j];
            }
          }
        }
        break;
      }

      case BELOW : {
        for (i = 0;i < seq->nb_sequence;i++) {
          for (j = 0;j < seq->length[i];j++) {
            if (real_sequence[i][variable][j] < threshold) {
              seq->real_sequence[i][variable][j] = threshold;
            }
            else {
              seq->real_sequence[i][variable][j] = real_sequence[i][variable][j];
            }
          }
        }
        break;
      }
      }

      break;
    }
    }

    switch (mode) {
    case ABOVE :
      seq->min_value[variable] = min_value[variable];
      seq->max_value[variable] = threshold;
      break;
    case BELOW :
      seq->min_value[variable] = threshold;
      seq->max_value[variable] = max_value[variable];
      break;
    }

    if ((seq->type[variable] == INT_VALUE) && (seq->min_value[variable] >= 0) &&
        (seq->max_value[variable] <= MARGINAL_DISTRIBUTION_MAX_VALUE)) {
      seq->build_marginal_frequency_distribution(variable);
    }
    else {
      seq->build_marginal_histogram(variable);
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Seuillage des valeurs d'une variable reelle.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              seuil, mode (ABOVE/BELOW).
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::thresholding(StatError &error , int variable , double threshold , int mode) const

{
  bool status = true;
  register int i , j;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (type[variable] != REAL_VALUE) {
      status = false;
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , STAT_variable_word[REAL_VALUE]);
    }

    if (threshold <= min_value[variable]) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_error[STATR_GREATER_THAN] << " " << min_value[variable];
      error.correction_update(STAT_error[STATR_THRESHOLD_VALUE] , (correction_message.str()).c_str());
    }

    if (threshold >= max_value[variable]) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_error[STATR_SMALLER_THAN] << " " << max_value[variable];
      error.correction_update(STAT_error[STATR_THRESHOLD_VALUE] , (correction_message.str()).c_str());
    }
  }

  if (status) {
    seq = new Sequences(*this , variable , type[variable]);

    // seuillage des valeurs reelles

    switch (mode) {

    case ABOVE : {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->length[i];j++) {
          if (real_sequence[i][variable][j] > threshold) {
            seq->real_sequence[i][variable][j] = threshold;
          }
          else {
            seq->real_sequence[i][variable][j] = real_sequence[i][variable][j];
          }
        }
      }

      seq->min_value[variable] = min_value[variable];
      seq->max_value[variable] = threshold;
      break;
    }

    case BELOW : {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->length[i];j++) {
          if (real_sequence[i][variable][j] < threshold) {
            seq->real_sequence[i][variable][j] = threshold;
          }
          else {
            seq->real_sequence[i][variable][j] = real_sequence[i][variable][j];
          }
        }
      }

      seq->min_value[variable] = threshold;
      seq->max_value[variable] = max_value[variable];
      break;
    }
    }

    seq->build_marginal_histogram(variable);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une variable.
 *
 *  arguments : reference sur un objet Sequences, indice de la variable,
 *              pas de regroupement, mode (FLOOR/ROUND/CEIL).
 *
 *--------------------------------------------------------------*/

void Sequences::cluster(const Sequences &seq , int variable , int step , int mode)

{
  register int i , j;


  switch (type[variable]) {

  // regroupement des valeurs entieres

  case INT_VALUE : {
    switch (mode) {

    case FLOOR : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          int_sequence[i][variable][j] = seq.int_sequence[i][variable][j] / step;
        }
      }

      min_value[variable] = (int)seq.min_value[variable] / step;
      max_value[variable] = (int)seq.max_value[variable] / step;
//      min_value[variable] = floor(seq.min_value[variable] / step);
//      max_value[variable] = floor(seq.max_value[variable] / step);
      break;
    }

    case ROUND : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          int_sequence[i][variable][j] = (seq.int_sequence[i][variable][j] + step / 2) / step;
//          int_sequence[i][variable][j] = (int)::round((double)seq.int_sequence[i][variable][j] / (double)step);
        }
      }

      min_value[variable] = ((int)seq.min_value[variable] + step / 2) / step;
      max_value[variable] = ((int)seq.max_value[variable] + step / 2) / step;
//      min_value[variable] = ::round(seq.min_value[variable] / step);
//      max_value[variable] = ::round(seq.max_value[variable] / step);
      break;
    }

    case CEIL : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          int_sequence[i][variable][j] = (seq.int_sequence[i][variable][j] + step - 1) / step;
//          int_sequence[i][variable][j] = (int)ceil((double)seq.int_sequence[i][variable][j] / (double)step);
        }
      }

      min_value[variable] = ((int)seq.min_value[variable] + step - 1) / step;
      max_value[variable] = ((int)seq.max_value[variable] + step - 1) / step;
//      min_value[variable] = ceil(seq.min_value[variable] / step);
//      max_value[variable] = ceil(seq.max_value[variable] / step);
      break;
    }
    }

    if (seq.marginal_distribution[variable]) {
      marginal_distribution[variable] = new FrequencyDistribution(*(seq.marginal_distribution[variable]) , 'c' , step , mode);
    }
    else {
      build_marginal_frequency_distribution(variable);
    }
    break;
  }

  // regroupement des valeurs reelles

  case REAL_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        real_sequence[i][variable][j] = seq.real_sequence[i][variable][j] / step;
      }
    }

    min_value[variable] = seq.min_value[variable] / step;
    max_value[variable] = seq.max_value[variable] / step;

    build_marginal_histogram(variable , seq.marginal_histogram[variable]->step / step);

    if ((variable + 1 < nb_variable) && (type[variable + 1] == AUXILIARY)) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          real_sequence[i][variable + 1][j] = seq.real_sequence[i][variable + 1][j] / step;
        }
      }

      min_value[variable + 1] = seq.min_value[variable + 1] / step;
      max_value[variable + 1] = seq.max_value[variable + 1] / step;
    }
    break;
  }
  }
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une variable.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              pas de regroupement, mode (FLOOR/ROUND/CEIL).
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::cluster(StatError &error , int variable , int step , int mode) const

{
  bool status = true;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != REAL_VALUE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[REAL_VALUE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    if ((type[variable] == INT_VALUE) && (variable + 1 < nb_variable) &&
        (type[variable + 1] == AUXILIARY)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      error.update((error_message.str()).c_str());
    }
  }

  if (step < 1) {
    status = false;
    error.update(STAT_error[STATR_CLUSTERING_STEP]);
  }

  if (status) {
    seq = new Sequences(*this , variable , type[variable]);
    seq->cluster(*this , variable , step , mode);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Transcodage des symboles d'une variable entiere.
 *
 *  arguments : reference sur un objet Sequences, indice de la variable,
 *              plus petit et plus grand symboles, table de transcodage des symboles,
 *              flag pour ajouter une variable.
 *
 *--------------------------------------------------------------*/

void Sequences::transcode(const Sequences &seq , int ivariable , int min_symbol ,
                          int max_symbol , int *symbol , bool add_flag)

{
  register int i , j , k;
  int variable , offset;


  // copy of index parameters

  if (seq.index_parameter_distribution) {
    index_parameter_distribution = new FrequencyDistribution(*(seq.index_parameter_distribution));
  }
  if (seq.index_interval) {
    index_interval = new FrequencyDistribution(*(seq.index_interval));
  }

  if (seq.index_parameter) {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < (index_parameter_type == POSITION ? length[i] + 1 : length[i]);j++) {
        index_parameter[i][j] = seq.index_parameter[i][j];
      }
    }
  }

  switch (add_flag) {
  case false :
    variable = ivariable;
    offset = 0;
    break;
  case true :
    variable = 0;
    offset = 1;
    break;
  }

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < nb_variable;j++) {
      if ((type[j] != REAL_VALUE) && (type[j] != AUXILIARY)) {

        // transcodage des symboles

        if (j == variable) {
          for (k = 0;k < length[i];k++) {
            int_sequence[i][j][k] = symbol[seq.int_sequence[i][ivariable][k] -
                                           (int)seq.min_value[variable]] + min_symbol;
          }
        }

        // copie des valeurs entieres

        else {
          for (k = 0;k < length[i];k++) {
            int_sequence[i][j][k] =  seq.int_sequence[i][j - offset][k];
          }
        }
      }

      // copie des valeurs reelles

      else {
        for (k = 0;k < length[i];k++) {
          real_sequence[i][j][k] = seq.real_sequence[i][j - offset][k];
        }
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    if (i == variable) {
      min_value[i] = min_symbol;
      max_value[i] = max_symbol;

      build_marginal_frequency_distribution(i);
    }

    else {
      min_value[i] = seq.min_value[i - offset];
      max_value[i] = seq.max_value[i - offset];

      if (seq.marginal_distribution[i - offset]) {
        marginal_distribution[i] = new FrequencyDistribution(*(seq.marginal_distribution[i - offset]));
      }
      if (seq.marginal_histogram[i - offset]) {
        marginal_histogram[i] = new Histogram(*(seq.marginal_histogram[i - offset]));
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Transcodage des symboles d'une variable entiere.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              table de transcodage des symboles.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::transcode(StatError &error , int variable , int *symbol) const

{
  bool status = true , *presence;
  register int i;
  int min_symbol , max_symbol;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != STATE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    if ((variable + 1 < nb_variable) && (type[variable + 1] == AUXILIARY)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      error.update((error_message.str()).c_str());
    }

    if (status) {
      min_symbol = symbol[0];
      max_symbol = symbol[0];

      for (i = 1;i <= (int)(max_value[variable] - min_value[variable]);i++) {
        if (symbol[i] < min_symbol) {
          min_symbol = symbol[i];
        }
        if (symbol[i] > max_symbol) {
          max_symbol = symbol[i];
        }
      }

      if (max_symbol - min_symbol == 0) {
        status = false;
        error.update(STAT_error[STATR_NB_SYMBOL]);
      }

      if (max_symbol - min_symbol > (int)(max_value[variable] - min_value[variable])) {
        status = false;
        error.update(STAT_error[STATR_NON_CONSECUTIVE_SYMBOLS]);
      }
    }

    if (status) {
      presence = new bool[max_symbol - min_symbol + 1];
      for (i = 0;i <= max_symbol - min_symbol;i++) {
        presence[i] = false;
      }

      for (i = 0;i <= (int)(max_value[variable] - min_value[variable]);i++) {
        presence[symbol[i] - min_symbol] = true;
      }

      for (i = 0;i <= max_symbol - min_symbol;i++) {
        if (!presence[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_error[STATR_MISSING_SYMBOL] << " " << i + min_symbol;
          error.update((error_message.str()).c_str());
        }
      }

      delete [] presence;
    }

    if (status) {
      for (i = 0;i <= (int)(max_value[variable] - min_value[variable]);i++) {
        symbol[i] -= min_symbol;
      }

      seq = new Sequences(nb_sequence , identifier , length , vertex_identifier ,
                          index_parameter_type , nb_variable , type);
      seq->transcode(*this , variable , min_symbol , max_symbol , symbol);
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une variable entiere.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              nombre de classes, bornes pour regrouper les valeurs.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::cluster(StatError &error , int variable ,
                              int nb_class , int *ilimit) const

{
  bool status = true;
  register int i , j , k;
  int *int_limit , *symbol , *itype;
  double *real_limit;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != STATE) &&
        (type[variable] != REAL_VALUE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    else if ((nb_class < 2) || (nb_class >= (int)(max_value[variable] - min_value[variable]) + 1)) {
      status = false;
      error.update(STAT_error[STATR_NB_CLASS]);
    }

    if ((variable + 1 < nb_variable) && (type[variable + 1] == AUXILIARY)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      error.update((error_message.str()).c_str());
    }
  }

  if (status) {
    if ((type[variable] == INT_VALUE) || (type[variable] == STATE)) {
      int_limit = new int[nb_class + 1];
      int_limit[0] = (int)min_value[variable];
      for (i = 1;i < nb_class;i++) {
        int_limit[i] = ilimit[i - 1];
      }
      int_limit[nb_class] = (int)max_value[variable] + 1;

      for (i = 0;i < nb_class;i++) {
        if (int_limit[i] >= int_limit[i + 1]) {
          status = false;
          error.update(STAT_error[STATR_CLUSTER_LIMIT]);
        }
      }

      if (status) {
        symbol = new int[(int)(max_value[variable] - min_value[variable]) + 1];

        i = 0;
        for (j = 0;j < nb_class;j++) {
          for (k = int_limit[j];k < int_limit[j + 1];k++) {
            symbol[i++] = j;
          }
        }

        seq = new Sequences(nb_sequence , identifier , length , vertex_identifier ,
                            index_parameter_type , nb_variable , type);
        seq->transcode(*this , variable , 0 , nb_class - 1 , symbol);

        delete [] symbol;
      }

      delete [] int_limit;
    }

    else {
      real_limit = new double[nb_class + 1];
      real_limit[0] = min_value[variable];
      for (i = 1;i < nb_class;i++) {
        real_limit[i] = ilimit[i - 1];
      }
      real_limit[nb_class] = max_value[variable] + 1;

      for (i = 0;i < nb_class;i++) {
        if (real_limit[i] >= real_limit[i + 1]) {
          status = false;
          error.update(STAT_error[STATR_CLUSTER_LIMIT]);
        }
      }

      if (status) {
        seq = new Sequences(*this , variable , INT_VALUE);
        seq->cluster(*this , variable , nb_class , real_limit);
      }

      delete [] real_limit;
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une variable reelle.
 *
 *  arguments : reference sur un objet Sequences, indice de la variable,
 *              nombre de classes, bornes pour regrouper les valeurs.
 *
 *--------------------------------------------------------------*/

void Sequences::cluster(const Sequences &seq , int variable , int nb_class , double *limit)

{
  register int i , j , k;


  // regroupement des valeurs reelles

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < length[i];j++) {
      for (k = 0;k < nb_class;k++) {
        if (seq.real_sequence[i][variable][j] < limit[k + 1]) {
          int_sequence[i][variable][j] = k;
          break;
        }
      }
    }
  }

  min_value_computation(variable);
  max_value_computation(variable);

  build_marginal_frequency_distribution(variable);
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une variable reelle.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              nombre de classes, bornes pour regrouper les valeurs.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::cluster(StatError &error , int variable ,
                              int nb_class , double *ilimit) const

{
  bool status = true;
  register int i;
  double *limit;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (type[variable] != REAL_VALUE) {
      status = false;
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , STAT_variable_word[REAL_VALUE]);
    }

    if ((variable + 1 < nb_variable) && (type[variable + 1] == AUXILIARY)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      error.update((error_message.str()).c_str());
    }
  }

  if (nb_class < 2) {
    status = false;
    error.update(STAT_error[STATR_NB_CLASS]);
  }

  if (status) {
    limit = new double[nb_class + 1];
    limit[0] = min_value[variable];
    for (i = 1;i < nb_class;i++) {
      limit[i] = ilimit[i - 1];
    }
    limit[nb_class] = max_value[variable] + DOUBLE_ERROR;

    for (i = 0;i < nb_class;i++) {
      if (limit[i] >= limit[i + 1]) {
        status = false;
        error.update(STAT_error[STATR_CLUSTER_LIMIT]);
      }
    }

    if (status) {
      seq = new Sequences(*this , variable , INT_VALUE);
      seq->cluster(*this , variable , nb_class , limit);
    }

    delete [] limit;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Changement d'unite d'une variable.
 *
 *  arguments : reference sur un objet StatError, variable, facteur d'echelle.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::scaling(StatError &error , int variable , int scaling_coeff) const

{
  bool status = true;
  register int i , j;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (type[variable] == INT_VALUE) {
      if ((min_value[variable] * scaling_coeff < INT_MIN) ||
          (max_value[variable] * scaling_coeff > INT_MAX)) {
        status = false;
        error.update(STAT_error[STATR_SCALING_COEFF]);
      }
    }

    else if (type[variable] != REAL_VALUE) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[REAL_VALUE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }
  }

  if (scaling_coeff <= 1) {
    status = false;
    error.update(STAT_error[STATR_SCALING_COEFF]);
  }

  if (status) {
    seq = new Sequences(*this , variable , type[variable]);

    switch (seq->type[variable]) {

    // mise a l'echelle des valeurs entieres

    case INT_VALUE : {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->length[i];j++) {
          seq->int_sequence[i][variable][j] = int_sequence[i][variable][j] * scaling_coeff;
        }
      }
      break;
    }

    // mise a l'echelle des valeurs reelles

    case REAL_VALUE : {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->length[i];j++) {
          seq->real_sequence[i][variable][j] = real_sequence[i][variable][j] * scaling_coeff;
        }
      }
      break;
    }
    }

    seq->min_value[variable] = min_value[variable] * scaling_coeff;
    seq->max_value[variable] = max_value[variable] * scaling_coeff;

    seq->build_marginal_frequency_distribution(variable);

    if ((variable + 1 < seq->nb_variable) && (seq->type[variable + 1] == AUXILIARY)) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->length[i];j++) {
          seq->real_sequence[i][variable + 1][j] = real_sequence[i][variable + 1][j] * scaling_coeff;
        }
      }

      seq->min_value[variable + 1] = min_value[variable + 1] * scaling_coeff;
      seq->max_value[variable + 1] = max_value[variable + 1] * scaling_coeff;
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Changement d'unite d'une variable.
 *
 *  arguments : reference sur un objet StatError, variable, facteur d'echelle.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::scaling(StatError &error , int variable , double scaling_coeff) const

{
  bool status = true;
  register int i , j;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != REAL_VALUE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[REAL_VALUE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }
  }

  if (scaling_coeff <= 0) {
    status = false;
    error.update(STAT_error[STATR_SCALING_COEFF]);
  }

  if (status) {
    seq = new Sequences(*this , variable , REAL_VALUE);

    switch (type[variable]) {

    // mise a l'echelle des valeurs entieres

    case INT_VALUE : {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->length[i];j++) {
          seq->real_sequence[i][variable][j] = int_sequence[i][variable][j] * scaling_coeff;
        }
      }
      break;
    }

    // mise a l'echelle des valeurs reelles

    case REAL_VALUE : {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->length[i];j++) {
          seq->real_sequence[i][variable][j] = real_sequence[i][variable][j] * scaling_coeff;
        }
      }
      break;
    }
    }

    seq->min_value[variable] = min_value[variable] * scaling_coeff;
    seq->max_value[variable] = max_value[variable] * scaling_coeff;

    seq->build_marginal_histogram(variable);

    if ((variable + 1 < seq->nb_variable) && (seq->type[variable + 1] == AUXILIARY)) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->length[i];j++) {
          seq->real_sequence[i][variable + 1][j] = real_sequence[i][variable + 1][j] * scaling_coeff;
        }
      }

      seq->min_value[variable + 1] = min_value[variable + 1] * scaling_coeff;
      seq->max_value[variable + 1] = max_value[variable + 1] * scaling_coeff;
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Arrondi des valeurs d'une variable reelle.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              mode (FLOOR/ROUND/CEIL).
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::round(StatError &error , int variable , int mode) const

{
  bool status = true;
  register int i , j , k;
  int *itype;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (variable != I_DEFAULT) {
    if ((variable < 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }

    else {
      variable--;

      if (type[variable] != REAL_VALUE) {
        status = false;
        error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , STAT_variable_word[REAL_VALUE]);
      }
    }
  }

  if (status) {
    for (i = 0;i < nb_variable;i++) {
      if (((variable == I_DEFAULT) && (type[i] == REAL_VALUE)) || (variable == i)) {
        if (min_value[i] < INT_MIN) {
          status = false;
          ostringstream error_message , correction_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_ROUNDED_VALUE];
          correction_message << STAT_error[STATR_GREATER_THAN] << " " << INT_MIN;
          error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
        }

        if (max_value[i] > INT_MAX) {
          status = false;
          ostringstream error_message , correction_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_ROUNDED_VALUE];
          correction_message << STAT_error[STATR_SMALLER_THAN] << " " << INT_MAX;
          error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
        }
      }
    }
  }

  if (status) {
    itype = new int[nb_variable];

    if (variable == I_DEFAULT) {
      for (i = 0;i < nb_variable;i++) {
        if (type[i] == REAL_VALUE) {
          itype[i] = INT_VALUE;
        }
        else {
          itype[i] = type[i];
        }
      }
    }

    else {
      for (i = 0;i < nb_variable;i++) {
        itype[i] = type[i];
      }
      itype[variable] = INT_VALUE;
    }

    seq = new Sequences(nb_sequence , identifier , length , vertex_identifier ,
                        index_parameter_type , nb_variable , itype);
    delete [] itype;

    // copy of index parameters

    if (index_parameter_distribution) {
      seq->index_parameter_distribution = new FrequencyDistribution(*index_parameter_distribution);
    }
    if (index_interval) {
      seq->index_interval = new FrequencyDistribution(*index_interval);
    }

    if (index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < (seq->index_parameter_type == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }
      }
    }

    for (i = 0;i < seq->nb_sequence;i++) {
      for (j = 0;j < seq->nb_variable;j++) {

        // copie des valeurs entieres

        if ((type[j] != REAL_VALUE) && (type[j] != AUXILIARY)) {
          for (k = 0;k < seq->length[i];k++) {
            seq->int_sequence[i][j][k] = int_sequence[i][j][k];
          }
        }

        else {

          // arrondi des valeurs reelles

          if (((variable == I_DEFAULT) && (type[j] == REAL_VALUE)) || (variable == j)) {
            switch (mode) {

            case FLOOR : {
              for (k = 0;k < seq->length[i];k++) {
                seq->int_sequence[i][j][k] = (int)floor(real_sequence[i][j][k]);
              }
              break;
            }

            case ROUND : {
              for (k = 0;k < seq->length[i];k++) {
                seq->int_sequence[i][j][k] = (int)::round(real_sequence[i][j][k]);
              }
              break;
            }

            case CEIL : {
              for (k = 0;k < seq->length[i];k++) {
                seq->int_sequence[i][j][k] = (int)ceil(real_sequence[i][j][k]);
              }
              break;
            }
            }
          }

          // copie des valeurs reelles

          else {
            for (k = 0;k < seq->length[i];k++) {
              seq->real_sequence[i][j][k] = real_sequence[i][j][k];
            }
          }
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      if (((variable == I_DEFAULT) && (type[i] == REAL_VALUE)) || (variable == i)) {
        switch (mode) {
        case FLOOR :
          seq->min_value[i] = floor(min_value[i]);
          seq->max_value[i] = floor(max_value[i]);
          break;
        case ROUND :
          seq->min_value[i] = ::round(min_value[i]);
          seq->max_value[i] = ::round(max_value[i]);
          break;
        case CEIL :
          seq->min_value[i] = ceil(min_value[i]);
          seq->max_value[i] = ceil(max_value[i]);
          break;
        }

        seq->build_marginal_frequency_distribution(i);
      }

      else {
        seq->min_value[i] = min_value[i];
        seq->max_value[i] = max_value[i];

        if (marginal_distribution[i]) {
          seq->marginal_distribution[i] = new FrequencyDistribution(*marginal_distribution[i]);
        }
        if (marginal_histogram[i]) {
          seq->marginal_histogram[i] = new Histogram(*marginal_histogram[i]);
        }
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Selection de sequences sur les valeurs prises par le parametre d'index.
 *
 *  arguments : reference sur un objet StatError, stream,
 *              bornes sur les parametres d'index,
 *              flag pour conserver ou rejeter les sequences selectionnees.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::index_parameter_select(StatError &error , ostream &os ,
                                             int min_index_parameter ,
                                             int max_index_parameter , bool keep) const

{
  bool status = true;
  register int i , j;
  int inb_sequence , *index , *iidentifier;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((index_parameter_type != TIME) && (index_parameter_type != POSITION)) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }

  else {
    if ((min_index_parameter < 0) || (min_index_parameter >= index_parameter_distribution->nb_value) ||
        (min_index_parameter > max_index_parameter)) {
      status = false;
      error.update(SEQ_error[SEQR_MIN_INDEX_PARAMETER]);
    }
    if ((max_index_parameter < index_parameter_distribution->offset) ||
        (max_index_parameter < min_index_parameter)) {
      status = false;
      error.update(SEQ_error[SEQR_MAX_INDEX_PARAMETER]);
    }
  }

  if (status) {

    // selection des sequences

    iidentifier = new int[nb_sequence];
    index = new int[nb_sequence];
    inb_sequence = 0;

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < (index_parameter_type == POSITION ? length[i] + 1 : length[i]);j++) {
        if ((index_parameter[i][j] >= min_index_parameter) &&
            (index_parameter[i][j] <= max_index_parameter)) {
          if (keep) {
            iidentifier[inb_sequence] = identifier[i];
            index[inb_sequence++] = i;
          }
          break;
        }
      }

      if ((!keep) && (j == (index_parameter_type == POSITION ? length[i] + 1 : length[i]))) {
        iidentifier[inb_sequence] = identifier[i];
        index[inb_sequence++] = i;
      }
    }

    if (inb_sequence == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    // copie des sequences

    if (status) {
      if (inb_sequence <= DISPLAY_NB_INDIVIDUAL) {
        os << "\n" << SEQ_label[inb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << ": ";
        for (i = 0;i < inb_sequence;i++) {
          os << iidentifier[i] << ", ";
        }
        os << endl;
      }

      seq = new Sequences(*this , inb_sequence , index);
    }

    delete [] iidentifier;
    delete [] index;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Selection de sequences sur les valeurs prises par une variable.
 *
 *  arguments : reference sur un objet StatError, stream, indice de la variable,
 *              bornes sur les valeurs, flag pour conserver ou rejeter
 *              les sequences selectionnees.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::value_select(StatError &error , ostream &os , int variable ,
                                   int imin_value , int imax_value , bool keep) const

{
  bool status = true;
  register int i , j;
  int inb_sequence , *index , *iidentifier;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != STATE) &&
        (type[variable] != REAL_VALUE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    else {
      if ((imin_value > max_value[variable]) || (imin_value > imax_value)) {
        status = false;
        error.update(STAT_error[STATR_MIN_VALUE]);
      }
      if ((imax_value < min_value[variable]) || (imax_value < imin_value)) {
        status = false;
        error.update(STAT_error[STATR_MAX_VALUE]);
      }
    }
  }

  if (status) {

    // selection des sequences

    iidentifier = new int[nb_sequence];
    index = new int[nb_sequence];
    inb_sequence = 0;

    if (type[variable] != REAL_VALUE) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          if ((int_sequence[i][variable][j] >= imin_value) &&
              (int_sequence[i][variable][j] <= imax_value)) {
            if (keep) {
              iidentifier[inb_sequence] = identifier[i];
              index[inb_sequence++] = i;
            }
            break;
          }
        }

        if ((!keep) && (j == length[i])) {
          iidentifier[inb_sequence] = identifier[i];
          index[inb_sequence++] = i;
        }
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          if ((real_sequence[i][variable][j] >= imin_value) &&
              (real_sequence[i][variable][j] <= imax_value)) {
            if (keep) {
              iidentifier[inb_sequence] = identifier[i];
              index[inb_sequence++] = i;
            }
            break;
          }
        }

        if ((!keep) && (j == length[i])) {
          iidentifier[inb_sequence] = identifier[i];
          index[inb_sequence++] = i;
        }
      }
    }

    if (inb_sequence == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    // copie des sequences

    if (status) {
      if (inb_sequence <= DISPLAY_NB_INDIVIDUAL) {
        os << "\n" << SEQ_label[inb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << ": ";
        for (i = 0;i < inb_sequence;i++) {
          os << iidentifier[i] << ", ";
        }
        os << endl;
      }

      seq = new Sequences(*this , inb_sequence , index);
    }

    delete [] iidentifier;
    delete [] index;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Selection de sequences sur les valeurs prises par une variable reelle.
 *
 *  arguments : reference sur un objet StatError, stream, indice de la variable,
 *              bornes sur les valeurs, flag pour conserver ou rejeter
 *              les sequences selectionnees.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::value_select(StatError &error , ostream &os , int variable ,
                                   double imin_value , double imax_value , bool keep) const

{
  bool status = true;
  register int i , j;
  int inb_sequence , *index , *iidentifier;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (type[variable] != REAL_VALUE) {
      status = false;
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , STAT_variable_word[REAL_VALUE]);
    }

    else {
      if ((imin_value > max_value[variable]) || (imin_value > imax_value)) {
        status = false;
        error.update(STAT_error[STATR_MIN_VALUE]);
      }
      if ((imax_value < min_value[variable]) || (imax_value < imin_value)) {
        status = false;
        error.update(STAT_error[STATR_MAX_VALUE]);
      }
    }
  }

  if (status) {

    // selection des sequences

    iidentifier = new int[nb_sequence];
    index = new int[nb_sequence];
    inb_sequence = 0;

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        if ((real_sequence[i][variable][j] >= imin_value) &&
            (real_sequence[i][variable][j] <= imax_value)) {
          if (keep) {
            iidentifier[inb_sequence] = identifier[i];
            index[inb_sequence++] = i;
          }
          break;
        }
      }

      if ((!keep) && (j == length[i])) {
        iidentifier[inb_sequence] = identifier[i];
        index[inb_sequence++] = i;
      }
    }

    if (inb_sequence == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    // copie des sequences

    if (status) {
      if (inb_sequence <= DISPLAY_NB_INDIVIDUAL) {
        os << "\n" << SEQ_label[inb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << ": ";
        for (i = 0;i < inb_sequence;i++) {
          os << iidentifier[i] << ", ";
        }
        os << endl;
      }

      seq = new Sequences(*this , inb_sequence , index);
    }

    delete [] iidentifier;
    delete [] index;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Selection de sequences par l'identificateur.
 *
 *  arguments : reference sur un objet StatError, nombre de sequences,
 *              identificateur des sequences, flag pour conserver ou rejeter
 *              les sequences selectionnees.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::select_individual(StatError &error , int inb_sequence ,
                                        int *iidentifier , bool keep) const

{
  bool status = true;
  int *index;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((inb_sequence < 1) || (inb_sequence > (keep ? nb_sequence : nb_sequence - 1))) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }

  else {
    status = selected_identifier_checking(error , nb_sequence , identifier , inb_sequence ,
                                          iidentifier , SEQ_label[SEQL_SEQUENCE]);
  }

  if (status) {
    index = identifier_select(nb_sequence , identifier , inb_sequence , iidentifier , keep);

    seq = new Sequences(*this , (keep ? inb_sequence : nb_sequence - inb_sequence) , index);

    delete [] index;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Suppression du parametre d'index.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::remove_index_parameter(StatError &error) const

{
  Sequences *seq;


  error.init();

  if (!index_parameter) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  else {
    seq = new Sequences(*this , 'm');
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Sequences avec transformation du parametre d'index implicite
 *  en parametre d'index explicite.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::explicit_index_parameter(StatError &error) const

{
  Sequences *seq;


  error.init();

  if (index_parameter) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  else {
    seq = new Sequences(*this , 'e');
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Selection de variables.
 *
 *  arguments : reference sur un objet Sequences, indices des variables.
 *
 *--------------------------------------------------------------*/

void Sequences::select_variable(const Sequences &seq , int *variable)

{
  register int i , j , k;


  // copy of index parameters

  if (seq.index_parameter_distribution) {
    index_parameter_distribution = new FrequencyDistribution(*(seq.index_parameter_distribution));
  }
  if (seq.index_interval) {
    index_interval = new FrequencyDistribution(*(seq.index_interval));
  }

  if (seq.index_parameter) {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < (index_parameter_type == POSITION ? length[i] + 1 : length[i]);j++) {
        index_parameter[i][j] = seq.index_parameter[i][j];
      }
    }
  }

  // copie des valeurs

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < nb_variable;j++) {
      if ((type[j] != REAL_VALUE) && (type[j] != AUXILIARY)) {
        for (k = 0;k < length[i];k++) {
          int_sequence[i][j][k] = seq.int_sequence[i][variable[j]][k];
        }
      }

      else {
        for (k = 0;k < length[i];k++) {
          real_sequence[i][j][k] = seq.real_sequence[i][variable[j]][k];
        }
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value[i] = seq.min_value[variable[i]];
    max_value[i] = seq.max_value[variable[i]];

    if (seq.marginal_distribution[variable[i]]) {
      marginal_distribution[i] = new FrequencyDistribution(*(seq.marginal_distribution[variable[i]]));
    }
    if (seq.marginal_histogram[variable[i]]) {
      marginal_histogram[i] = new Histogram(*(seq.marginal_histogram[variable[i]]));
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Selection de variables.
 *
 *  arguments : reference sur un objet StatError, nombre de variables,
 *              indices des variables, flag pour conserver ou rejeter
 *              les variables selectionnees.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::select_variable(StatError &error , int inb_variable ,
                                      int *ivariable , bool keep) const

{
  bool status = true , *selected_variable;
  register int i;
  int bnb_variable , *variable , *itype;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((inb_variable < 1) || (inb_variable > (keep ? nb_variable : nb_variable - 1))) {
    status = false;
    error.update(STAT_error[STATR_NB_SELECTED_VARIABLE]);
  }

  else {
    selected_variable = new bool[nb_variable + 1];
    for (i = 1;i <= nb_variable;i++) {
      selected_variable[i] = false;
    }

    for (i = 0;i < inb_variable;i++) {
      if ((ivariable[i] < 1) || (ivariable[i] > nb_variable)) {
        status = false;
        ostringstream error_message;
        error_message << ivariable[i] << ": " << STAT_error[STATR_VARIABLE_INDEX];
        error.update((error_message.str()).c_str());
      }

      else if (selected_variable[ivariable[i]]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << ivariable[i] << " "
                      << STAT_error[STATR_ALREADY_SELECTED];
        error.update((error_message.str()).c_str());
      }
      else {
        selected_variable[ivariable[i]] = true;
      }
    }

    delete [] selected_variable;
  }

  if (status) {
    variable = ::select_variable(nb_variable , inb_variable , ivariable , keep);

    bnb_variable = (keep ? inb_variable : nb_variable - inb_variable);

    for (i = 0;i < bnb_variable;i++) {
      if ((type[variable[i]] == AUXILIARY) &&
          ((i == 0) || (variable[i - 1] != variable[i] - 1))) {
        status = false;
        error.update(SEQ_error[SEQR_VARIABLE_INDICES]);
      }
    }

    if (status) {
      itype = new int[bnb_variable];
      for (i = 0;i < bnb_variable;i++) {
        itype[i] = type[variable[i]];
      }

      seq = new Sequences(nb_sequence , identifier , length , vertex_identifier ,
                          index_parameter_type , bnb_variable , itype);
      seq->select_variable(*this , variable);

      delete [] itype;
    }

    delete [] variable;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Concatenation des variables d'objets Sequences.
 *
 *  arguments : reference sur un objet StatError, nombre d'objets Sequences,
 *              pointeurs sur les objets Sequences, echantillon de reference pour
 *              les identificateurs.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::merge_variable(StatError &error , int nb_sample ,
                                     const Sequences **iseq , int ref_sample) const

{
  bool status = true;
  register int i , j , k , m;
  int inb_variable , *iidentifier , *itype , **ivertex_identifier;
  Sequences *seq;
  const Sequences **pseq;


  seq = NULL;
  error.init();

  for (i = 0;i < nb_sample;i++) {
    if ((iseq[i]->vertex_identifier) && (!vertex_identifier)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_SAMPLE_VERTEX_IDENTIFIER];
      error.update((error_message.str()).c_str());
    }

    if ((iseq[i]->index_parameter_type != IMPLICIT_TYPE) &&
        (iseq[i]->index_parameter_type != index_parameter_type)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_INDEX_PARAMETER_TYPE];

      if (index_parameter_type == IMPLICIT_TYPE) {
        error.update((error_message.str()).c_str());
      }
      else {
        error.correction_update((error_message.str()).c_str() , SEQ_index_parameter_word[index_parameter_type]);
      }
    }

    if (iseq[i]->nb_sequence != nb_sequence) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_NB_SEQUENCE];
      error.update((error_message.str()).c_str());
    }

    else {
      for (j = 0;j < nb_sequence;j++) {
        if (iseq[i]->length[j] != length[j]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                        << SEQ_label[SEQL_SEQUENCE] << " " << j + 1 << ": "
                        << SEQ_error[SEQR_SEQUENCE_LENGTH];
          error.update((error_message.str()).c_str());
        }

        else {
          if ((iseq[i]->vertex_identifier) && (vertex_identifier)) {
            for (k = 0;k < length[j];k++) {
              if (iseq[i]->vertex_identifier[j][k] != vertex_identifier[j][k]) {
                status = false;
                ostringstream error_message;
                error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                              << SEQ_label[SEQL_SEQUENCE] << " " << j + 1 << ": "
                              << SEQ_label[SEQL_VERTEX_IDENTIFIER] << " " << k << ": "
                              << SEQ_error[SEQR_VERTEX_IDENTIFIER];
                error.update((error_message.str()).c_str());
              }
            }
          }

          if ((iseq[i]->index_parameter_type != IMPLICIT_TYPE) &&
              (iseq[i]->index_parameter_type == index_parameter_type)) {
            for (k = 0;k < (index_parameter_type == POSITION ? length[j] + 1 : length[j]);k++) {
              if (iseq[i]->index_parameter[j][k] != index_parameter[j][k]) {
                status = false;
                ostringstream error_message;
                error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                              << SEQ_label[SEQL_SEQUENCE] << " " << j + 1 << ": "
                              << SEQ_label[SEQL_INDEX_PARAMETER] << " " << k << ": "
                              << SEQ_error[SEQR_INDEX_PARAMETER];
                error.update((error_message.str()).c_str());
              }
            }
          }
        }
      }
    }
  }

  if ((ref_sample != I_DEFAULT) && ((ref_sample < 1) || (ref_sample > nb_sample + 1))) {
    status = false;
    error.update(STAT_error[STATR_SAMPLE_INDEX]);
  }

  if (status) {
    nb_sample++;
    pseq = new const Sequences*[nb_sample];

    pseq[0] = this;
    inb_variable = nb_variable;
    for (i = 1;i < nb_sample;i++) {
      pseq[i] = iseq[i - 1];
      inb_variable += iseq[i - 1]->nb_variable;
    }

    // comparaison des identificateurs des sequences

    if (ref_sample == I_DEFAULT) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 1;j < nb_sample;j++) {
          if (pseq[j]->identifier[i] != pseq[0]->identifier[i]) {
            break;
          }
        }
        if (j < nb_sample) {
          break;
        }
      }

      if (i < nb_sequence) {
        iidentifier = NULL;
      }
      else {
        iidentifier = pseq[0]->identifier;
      }
    }

    else {
      ref_sample--;
      iidentifier = pseq[ref_sample]->identifier;
    }

    // comparaison des identificateurs des vertex

    if (ref_sample == I_DEFAULT) {
      for (i = 0;i < nb_sample;i++) {
        if (!(pseq[i]->vertex_identifier)) {
          break;
        }
      }

      if (i == nb_sample) {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 1;j < nb_sample;j++) {
            for (k = 0;k < pseq[j]->length[i];k++) {
              if (pseq[j]->vertex_identifier[i][k] != pseq[0]->vertex_identifier[i][k]) {
                break;
              }
            }

            if (k < pseq[j]->length[i]) {
              break;
            }
          }

          if (j < nb_sample) {
            break;
          }
        }

        if (i < nb_sequence) {
          ivertex_identifier = NULL;
        }
        else {
          ivertex_identifier = pseq[0]->vertex_identifier;
        }
      }

      else {
        ivertex_identifier = NULL;
      }
    }

    else {
      ivertex_identifier = pseq[ref_sample]->vertex_identifier;
    }

    itype = new int[inb_variable];
    inb_variable = 0;
    for (i = 0;i < nb_sample;i++) {
      for (j = 0;j < pseq[i]->nb_variable;j++) {
        itype[inb_variable] = pseq[i]->type[j];
        if ((inb_variable > 0) && (itype[inb_variable] == STATE)) {
          itype[inb_variable] = INT_VALUE;
        }
        inb_variable++;
      }
    }

    seq = new Sequences(nb_sequence , iidentifier , length , ivertex_identifier ,
                        index_parameter_type , inb_variable , itype);
    delete [] itype;

    // copy of index parameters

    if (index_parameter_distribution) {
      seq->index_parameter_distribution = new FrequencyDistribution(*index_parameter_distribution);
    }
    if (index_interval) {
      seq->index_interval = new FrequencyDistribution(*index_interval);
    }

    if (index_parameter) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < (index_parameter_type == POSITION ? length[i] + 1 : length[i]);j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }
      }
    }

    // copie des valeurs

    for (i = 0;i < nb_sequence;i++) {
      inb_variable = 0;
      for (j = 0;j < nb_sample;j++) {
        for (k = 0;k < pseq[j]->nb_variable;k++) {
          if ((seq->type[inb_variable] != REAL_VALUE) && (seq->type[inb_variable] != AUXILIARY)) {
            for (m = 0;m < length[i];m++) {
              seq->int_sequence[i][inb_variable][m] = pseq[j]->int_sequence[i][k][m];
            }
          }

          else {
            for (m = 0;m < length[i];m++) {
              seq->real_sequence[i][inb_variable][m] = pseq[j]->real_sequence[i][k][m];
            }
          }

          inb_variable++;
        }
      }
    }

    inb_variable = 0;
    for (i = 0;i < nb_sample;i++) {
      for (j = 0;j < pseq[i]->nb_variable;j++) {
        seq->min_value[inb_variable] = pseq[i]->min_value[j];
        seq->max_value[inb_variable] = pseq[i]->max_value[j];

        if (pseq[i]->marginal_distribution[j]) {
          seq->marginal_distribution[inb_variable] = new FrequencyDistribution(*(pseq[i]->marginal_distribution[j]));
        }
        if (pseq[i]->marginal_histogram[j]) {
          seq->marginal_histogram[inb_variable] = new Histogram(*(pseq[i]->marginal_histogram[j]));
        }

        inb_variable++;
      }
    }

    delete [] pseq;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Decalage d'une variable.
 *
 *  arguments : reference sur un objet StatError, indice de la variable, decalage.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::shift_variable(StatError &error , int variable , int lag) const

{
  bool status = true;
  register int i , j , k;
  int inb_sequence , *iidentifier , *ilength , *index , *pvertex_id , *cvertex_id ,
      *pindex_param , *cindex_param , *pisequence , *cisequence;
  double *prsequence , *crsequence;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (index_parameter_type == POSITION) {
    error.update(SEQ_error[SEQR_INDEX_PARAMETER]);
    status = false;
  }

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (type[variable] == AUXILIARY) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[REAL_VALUE] << " or "
                         << STAT_variable_word[STATE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }
  }

  if ((lag == 0) || (abs(lag) > length_distribution->mean)) {
    status = false;
    error.update(STAT_error[SEQR_VARIABLE_LAG]);
  }

  if (status) {

    // calcul de la longueur des sequences

    iidentifier = new int[nb_sequence];
    ilength = new int[nb_sequence];
    index = new int[nb_sequence];
    inb_sequence = 0;

    for (i = 0;i < nb_sequence;i++) {
      if (abs(lag) < length[i]) {
        iidentifier[inb_sequence] = identifier[i];
        ilength[inb_sequence] = length[i] - abs(lag);
        index[inb_sequence++] = i;
      }
    }

    seq = new Sequences(inb_sequence , iidentifier , ilength , vertex_identifier ,
                        index_parameter_type , nb_variable , type , false);

    // copy of vertex identifiers

    if (vertex_identifier) {
      for (i = 0;i < seq->nb_sequence;i++) {
        pvertex_id = seq->vertex_identifier[i];

        if (lag > 0) {
          cvertex_id = vertex_identifier[index[i]] + lag;
        }
        else {
          cvertex_id = vertex_identifier[index[i]];
        }

        for (j = 0;j < seq->length[i];j++) {
          *pvertex_id++ = *cvertex_id++;
        }
      }
    }

    // copy of index parameters

    if (index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        pindex_param = seq->index_parameter[i];

        if (lag > 0) {
          cindex_param = index_parameter[index[i]] + lag;
        }
        else {
          cindex_param = index_parameter[index[i]];
        }

        for (j = 0;j < seq->length[i];j++) {
          *pindex_param++ = *cindex_param++;
        }
      }

      seq->build_index_parameter_frequency_distribution();
      seq->index_interval_computation();
    }

    // copie des valeurs

    for (i = 0;i < seq->nb_sequence;i++) {
      for (j = 0;j < seq->nb_variable;j++) {
        if ((seq->type[j] != REAL_VALUE) && (seq->type[j] != AUXILIARY)) {
          pisequence = seq->int_sequence[i][j];

          if ((j != variable) && (lag > 0)) {
            cisequence = int_sequence[index[i]][j] + lag;
          }
          else if ((j == variable) && (lag < 0)) {
            cisequence = int_sequence[index[i]][j] - lag;
          }
          else {
            cisequence = int_sequence[index[i]][j];
          }

          for (k = 0;k < seq->length[i];k++) {
            *pisequence++ = *cisequence++;
          }
        }

        else {
          prsequence = seq->real_sequence[i][j];

          if ((j != variable) && (lag > 0))  {
            crsequence = real_sequence[index[i]][j] + lag;
          }
          else if ((j == variable) && (lag < 0)) {
            crsequence = real_sequence[index[i]][j] - lag;
          }
          else {
            crsequence = real_sequence[index[i]][j];
          }

          for (k = 0;k < seq->length[i];k++) {
            *prsequence++ = *crsequence++;
          }
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);

      if (seq->type[i] != AUXILIARY) {
        seq->build_marginal_frequency_distribution(i);
      }
    }

    delete [] iidentifier;
    delete [] ilength;
    delete [] index;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Inversion du sens de parcours des sequences.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::reverse(StatError &error) const

{
  Sequences *seq;


  error.init();

  if (index_parameter_type == TIME) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER]);
  }

  else {
    seq = new Sequences(*this , 'r');
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Selection des sequences sur un critere de longueur.
 *
 *  arguments : reference sur un objet StatError, stream, bornes sur la longueur,
 *              flag pour conserver ou rejeter les sequences selectionnees.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::length_select(StatError &error , ostream &os , int min_length ,
                                    int imax_length , bool keep) const

{
  bool status = true;
  register int i;
  int inb_sequence , *index , *iidentifier;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((min_length < 1) || (min_length >= length_distribution->nb_value) || (min_length > imax_length)) {
    status = false;
    error.update(SEQ_error[SEQR_MIN_SEQUENCE_LENGTH]);
  }
  if ((imax_length < length_distribution->offset) || (imax_length < min_length)) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_SEQUENCE_LENGTH]);
  }

  if (status) {

    // selection des sequences

    iidentifier = new int[nb_sequence];
    index = new int[nb_sequence];
    inb_sequence = 0;

    for (i = 0;i < nb_sequence;i++) {
      if ((length[i] >= min_length) && (length[i] <= imax_length)) {
        if (keep) {
          iidentifier[inb_sequence] = identifier[i];
          index[inb_sequence++] = i;
        }
      }

      else if (!keep) {
        iidentifier[inb_sequence] = identifier[i];
        index[inb_sequence++] = i;
      }
    }

    if (inb_sequence == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    // copie des sequences selectionnees

    if (status) {
      if (inb_sequence <= DISPLAY_NB_INDIVIDUAL) {
        os << "\n" << SEQ_label[inb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << ": ";
        for (i = 0;i < inb_sequence;i++) {
          os << iidentifier[i] << ", ";
        }
        os << endl;
      }

      seq = new Sequences(*this , inb_sequence , index);
    }

    delete [] iidentifier;
    delete [] index;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Suppression des premieres/dernieres series d'une valeur donne.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              valeur, position ('b' : begin, 'e' : end),
 *              longueur maximum de la serie supprimee.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::remove_run(StatError &error , int variable , int ivalue ,
                                 char position , int max_run_length) const

{
  bool status = true;
  register int i , j , k;
  int smax_length , inb_sequence , *iidentifier , *ilength , *index , *pvertex_id ,
      *cvertex_id , *pindex_param , *cindex_param , *pisequence , *cisequence;
  double *prsequence , *crsequence;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (index_parameter_type == POSITION) {
    error.update(SEQ_error[SEQR_INDEX_PARAMETER]);
    status = false;
  }

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != STATE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    else {
      if (!marginal_distribution[variable]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
        error.update((error_message.str()).c_str());
      }

      else if ((ivalue < marginal_distribution[variable]->offset) ||
               (ivalue >= marginal_distribution[variable]->nb_value) ||
               (marginal_distribution[variable]->frequency[ivalue] == 0)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << STAT_label[STATL_VALUE] << " " << ivalue << " "
                      << STAT_error[STATR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if ((max_run_length < I_DEFAULT) || (max_run_length == 0)) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_RUN_LENGTH]);
  }

  if (status) {

    // calcul de la longueur des sequences

    iidentifier = new int[nb_sequence];
    ilength = new int[nb_sequence];
    index = new int[nb_sequence];
    inb_sequence = 0;

    for (i = 0;i < nb_sequence;i++) {
      if (max_run_length == I_DEFAULT) {
        smax_length = length[i];
      }
      else {
        smax_length = MIN(max_run_length , length[i]);
      }

      switch (position) {

      case 'b' : {
        cisequence = int_sequence[i][variable];
        for (j = 0;j < smax_length;j++) {
          if (*cisequence++ != ivalue) {
            break;
          }
        }
        break;
      }

      case 'e' : {
        cisequence = int_sequence[i][variable] + length[i];
        for (j = 0;j < smax_length;j++) {
          if (*--cisequence != ivalue) {
            break;
          }
        }
        break;
      }
      }

      if (j < length[i]) {
        iidentifier[inb_sequence] = identifier[i];
        ilength[inb_sequence] = length[i] - j;
        index[inb_sequence++] = i;
      }
    }

    seq = new Sequences(inb_sequence , iidentifier , ilength , vertex_identifier ,
                        index_parameter_type , nb_variable , type , false);

    // copy of vertex identifiers

    if (vertex_identifier) {
      for (i = 0;i < seq->nb_sequence;i++) {
        pvertex_id = seq->vertex_identifier[i];

        switch (position) {
        case 'b' :
          cvertex_id = vertex_identifier[index[i]] + length[index[i]] - seq->length[i];
          break;
        case 'e' :
          cvertex_id = vertex_identifier[index[i]];
          break;
        }

        for (j = 0;j < seq->length[i];j++) {
          *pvertex_id++ = *cvertex_id++;
        }
      }
    }

    // copy of index parameters

    if (index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        pindex_param = seq->index_parameter[i];

        switch (position) {
        case 'b' :
          cindex_param = index_parameter[index[i]] + length[index[i]] - seq->length[i];
          break;
        case 'e' :
          cindex_param = index_parameter[index[i]];
          break;
        }

        for (j = 0;j < seq->length[i];j++) {
          *pindex_param++ = *cindex_param++;
        }
      }

      seq->build_index_parameter_frequency_distribution();
      seq->index_interval_computation();
    }

    // copie des valeurs

    for (i = 0;i < seq->nb_sequence;i++) {
      for (j = 0;j < seq->nb_variable;j++) {
        if ((seq->type[j] != REAL_VALUE) && (seq->type[j] != AUXILIARY)) {
          pisequence = seq->int_sequence[i][j];

          switch (position) {
          case 'b' :
            cisequence = int_sequence[index[i]][j] + length[index[i]] - seq->length[i];
            break;
          case 'e' :
            cisequence = int_sequence[index[i]][j];
            break;
          }

          for (k = 0;k < seq->length[i];k++) {
            *pisequence++ = *cisequence++;
          }
        }

        else {
          prsequence = seq->real_sequence[i][j];

          switch (position) {
          case 'b' :
            crsequence = real_sequence[index[i]][j] + length[index[i]] - seq->length[i];
            break;
          case 'e' :
            crsequence = real_sequence[index[i]][j];
            break;
          }

          for (k = 0;k < seq->length[i];k++) {
            *prsequence++ = *crsequence++;
          }
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);

      seq->build_marginal_frequency_distribution(i);
    }

    delete [] iidentifier;
    delete [] ilength;
    delete [] index;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de sous-sequences.
 *
 *  arguments : reference sur un objet StatError, parametres d'index minimum et
 *              maximum dans la sequence.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::index_parameter_extract(StatError &error , int min_index_parameter ,
                                              int max_index_parameter) const

{
  bool status = true;
  register int i , j , k;
  int inb_sequence , *iidentifier , *ilength , *index , *first_index , *pvertex_id ,
      *cvertex_id , *pindex_param , *cindex_param , *pisequence , *cisequence;
  double *prsequence , *crsequence;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((min_index_parameter < 0) || ((!index_parameter) && (min_index_parameter >= max_length)) ||
      ((max_index_parameter != I_DEFAULT) && (min_index_parameter > max_index_parameter))) {
    status = false;
    error.update(SEQ_error[SEQR_MIN_INDEX_PARAMETER]);
  }
  if ((max_index_parameter != I_DEFAULT) && ((max_index_parameter < 0) || ((!index_parameter) &&
        (max_index_parameter >= max_length)) || (max_index_parameter < min_index_parameter))) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_INDEX_PARAMETER]);
  }

  if (status) {

    // selection des sequences

    iidentifier = new int[nb_sequence];
    ilength = new int[nb_sequence];
    index = new int[nb_sequence];
    inb_sequence = 0;

    // parametre d'index explicite

    if (index_parameter) {
      first_index = new int[nb_sequence];

      if (max_index_parameter == I_DEFAULT) {
        for (i = 0;i < nb_sequence;i++) {
          if (index_parameter[i][length[i] - 1] >= min_index_parameter) {
            iidentifier[inb_sequence] = identifier[i];

            for (j = 0;j < length[i];j++) {
              if (index_parameter[i][j] >= min_index_parameter) {
                break;
              }
            }
            first_index[inb_sequence] = j;
            ilength[inb_sequence] = length[i] - j;

            index[inb_sequence++] = i;
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          if ((index_parameter[i][0] <= min_index_parameter) &&
              (index_parameter[i][length[i] - 1] >= max_index_parameter)) {
            iidentifier[inb_sequence] = identifier[i];

            for (j = 0;j < length[i];j++) {
              if (index_parameter[i][j] >= min_index_parameter) {
                break;
              }
            }
            first_index[inb_sequence] = j;
            ilength[inb_sequence] = length[i] - j;

            for (j = length[i] - 1;j >= 0;j--) {
              if (index_parameter[i][j] <= max_index_parameter) {
                break;
              }
            }
            ilength[inb_sequence] -= (length[i] - 1 - j);

            index[inb_sequence++] = i;
          }
        }
      }
    }

    // parametre d'index implicite

    else {
      if (max_index_parameter == I_DEFAULT) {
        for (i = 0;i < nb_sequence;i++) {
          if (length[i] > min_index_parameter) {
            iidentifier[inb_sequence] = identifier[i];
            ilength[inb_sequence] = length[i] - min_index_parameter;
            index[inb_sequence++] = i;
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          if (length[i] > max_index_parameter) {
            iidentifier[inb_sequence] = identifier[i];
            ilength[inb_sequence] = max_index_parameter - min_index_parameter + 1;
            index[inb_sequence++] = i;
          }
        }
      }
    }

    if (inb_sequence == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    // extraction des sous-sequences

    if (status) {
      seq = new Sequences(inb_sequence , iidentifier , ilength , vertex_identifier ,
                          index_parameter_type , nb_variable , type , false);

      // copy of vertex identifiers

      if (vertex_identifier) {
        for (i = 0;i < seq->nb_sequence;i++) {
          pvertex_id = seq->vertex_identifier[i];
          cvertex_id = vertex_identifier[index[i]] + first_index[i];
          for (k = 0;k < seq->length[i];k++) {
            *pvertex_id++ = *cvertex_id++;
          }
        }
      }

      // copy of index parameters

      if (index_parameter) {
        for (i = 0;i < seq->nb_sequence;i++) {
          pindex_param = seq->index_parameter[i];
          cindex_param = index_parameter[index[i]] + first_index[i];
          for (k = 0;k < seq->length[i];k++) {
            *pindex_param++ = *cindex_param++;
          }

          if (seq->index_parameter_type == POSITION) {
            if (max_index_parameter == I_DEFAULT) {
              *pindex_param = *cindex_param;
            }
            else {
              *pindex_param = max_index_parameter;
            }
          }
        }
      }

      // copie des valeurs

      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->nb_variable;j++) {
          if ((seq->type[j] != REAL_VALUE) && (seq->type[j] != AUXILIARY)) {
            pisequence = seq->int_sequence[i][j];
            if (seq->index_parameter) {
              cisequence = int_sequence[index[i]][j] + first_index[i];
            }
            else {
              cisequence = int_sequence[index[i]][j] + min_index_parameter;
            }

            for (k = 0;k < seq->length[i];k++) {
              *pisequence++ = *cisequence++;
            }
          }

          else {
            prsequence = seq->real_sequence[i][j];
            if (seq->index_parameter) {
              crsequence = real_sequence[index[i]][j] + first_index[i];
            }
            else {
              crsequence = real_sequence[index[i]][j] + min_index_parameter;
            }

            for (k = 0;k < seq->length[i];k++) {
              *prsequence++ = *crsequence++;
            }
          }
        }
      }

      if (index_parameter_distribution) {
        seq->build_index_parameter_frequency_distribution();
      }
      if (index_interval) {
        seq->index_interval_computation();
      }

      for (i = 0;i < seq->nb_variable;i++) {
        seq->min_value_computation(i);
        seq->max_value_computation(i);

        seq->build_marginal_frequency_distribution(i);
      }

      delete [] iidentifier;
      delete [] ilength;
      delete [] index;

      if (index_parameter) {
        delete [] first_index;
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Extraction par segmentation d'un objet Sequences.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              nombre de valeurs, valeurs, flag segments correspondant aux valeurs
 *              selectionnees ou non, segments concatenes par sequence ou non.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::segmentation_extract(StatError &error , int variable ,
                                           int nb_value , int *ivalue , bool keep ,
                                           bool concatenation) const

{
  bool status = true;
  register int i , j , k , m , n;
  int nb_present_value , nb_selected_value , nb_segment , inb_sequence , *pfrequency ,
      *selected_value , *itype , *pvertex_id , *cvertex_id , *pindex_param , *cindex_param ,
      *pisequence , *cisequence , *segment_length , *sequence_length , **segment_begin;
  double *prsequence , *crsequence;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (index_parameter_type == POSITION) {
    error.update(SEQ_error[SEQR_INDEX_PARAMETER]);
    status = false;
  }
  if (nb_variable == 1) {
    status = false;
    error.update(STAT_error[STATR_NB_VARIABLE]);
  }

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != STATE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    else {
      if (!marginal_distribution[variable]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
        error.update((error_message.str()).c_str());
      }

      else {
        pfrequency = marginal_distribution[variable]->frequency + marginal_distribution[variable]->offset;
        nb_present_value = 0;
        for (i = marginal_distribution[variable]->offset;i < marginal_distribution[variable]->nb_value;i++) {
          if (*pfrequency > 0) {
            nb_present_value++;
          }
          pfrequency++;
        }

        if ((nb_value < 1) || (nb_value > (keep ? nb_present_value : nb_present_value - 1))) {
          status = false;
          error.update(SEQ_error[SEQR_NB_SELECTED_VALUE]);
        }

        else {
          selected_value = new int[marginal_distribution[variable]->nb_value];
          for (i = marginal_distribution[variable]->offset;i < marginal_distribution[variable]->nb_value;i++) {
            selected_value[i] = false;
          }

          for (i = 0;i < nb_value;i++) {
            if ((ivalue[i] < marginal_distribution[variable]->offset) ||
                (ivalue[i] >= marginal_distribution[variable]->nb_value) ||
                (marginal_distribution[variable]->frequency[ivalue[i]] == 0)) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                            << STAT_label[STATL_VALUE] << " " << ivalue[i] << " "
                            << STAT_error[STATR_NOT_PRESENT];
              error.update((error_message.str()).c_str());
            }

            else if (selected_value[ivalue[i]]) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_VALUE] << " " << ivalue[i] << " "
                            << STAT_error[STATR_ALREADY_SELECTED];
              error.update((error_message.str()).c_str());
            }
            else {
              selected_value[ivalue[i]] = true;
            }
          }

          delete [] selected_value;
        }
      }
    }

    if ((variable + 1 < nb_variable) && ((type[variable + 1] != INT_VALUE) &&
         (type[variable + 1] != STATE) && (type[variable + 1] != REAL_VALUE))) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (status) {
    itype = new int[nb_variable - 1];
    for (i = 0;i < variable;i++) {
      itype[i] = type[i];
    }
    for (i = variable + 1;i < nb_variable;i++) {
      itype[i - 1] = type[i];
    }

    switch (keep) {

    case false : {
      nb_selected_value = nb_present_value - nb_value;
      selected_value = new int[nb_selected_value];

      pfrequency = marginal_distribution[variable]->frequency + marginal_distribution[variable]->offset;
      i = 0;

      for (j = marginal_distribution[variable]->offset;j < marginal_distribution[variable]->nb_value;j++) {
        if (*pfrequency > 0) {
          for (k = 0;k < nb_value;k++) {
            if (ivalue[k] == j) {
              break;
            }
          }

          if (k == nb_value) {
            selected_value[i++] = j;
          }
        }

        pfrequency++;
      }
      break;
    }

    case true : {
      nb_selected_value = nb_value;
      selected_value = ivalue;
      break;
    }
    }

    // initialisations

    segment_length = new int[cumul_length];

    segment_begin = new int*[cumul_length];
    for (i = 0;i < cumul_length;i++) {
      segment_begin[i] = 0;
    }

    // recherche des sequences a extraire

    i = -1;

    for (j = 0;j < nb_sequence;j++) {
      pisequence = int_sequence[j][variable];

      for (k = 0;k < nb_selected_value;k++) {
        if (*pisequence == selected_value[k]) {
          break;
        }
      }

      if (k < nb_selected_value) {
        i++;
        segment_begin[i] = new int[2];

        segment_begin[i][0] = j;
        segment_begin[i][1] = 0;
        segment_length[i] = 1;
      }

      for (k = 1;k < length[j];k++) {
        pisequence++;

        for (m = 0;m < nb_selected_value;m++) {
          if (*pisequence == selected_value[m]) {
            break;
          }
        }

        if (m  < nb_selected_value) {
          if (*pisequence != *(pisequence - 1)) {
            i++;
            segment_begin[i] = new int[2];

            segment_begin[i][0] = j;
            segment_begin[i][1] = k;
            segment_length[i] = 1;
          }
          else {
            segment_length[i]++;
          }
        }
      }
    }

    nb_segment = i + 1;

    switch (concatenation) {

    case false : {
      inb_sequence = nb_segment;
      sequence_length = segment_length;
      break;
    }

    case true : {
      sequence_length = new int[nb_sequence];

      i = 0;
      sequence_length[i] = segment_length[i];
      for (j = 1;j < nb_segment;j++) {
        if (segment_begin[j][0] != segment_begin[j - 1][0]) {
          i++;
          sequence_length[i] = 0;
        }
        sequence_length[i] += segment_length[j];
      }

      inb_sequence = i + 1;
      break;
    }
    }

    // creation de l'objet Sequences

    seq = new Sequences(inb_sequence , NULL , sequence_length , vertex_identifier ,
                        index_parameter_type , nb_variable - 1 , itype , false);

    // copy of vertex identifiers

    if (vertex_identifier) {
      i = -1;

      for (j = 0;j < nb_segment;j++) {
        switch (concatenation) {

        case false : {
          pvertex_id = seq->vertex_identifier[j];
          break;
        }

        case true : {
          if ((j == 0) || (segment_begin[j][0] != segment_begin[j - 1][0])) {
            pvertex_id = seq->vertex_identifier[++i];
          }
          break;
        }
        }

        cvertex_id = vertex_identifier[segment_begin[j][0]] + segment_begin[j][1];
        for (k = 0;k < segment_length[j];k++) {
          *pvertex_id++ = *cvertex_id++;
        }
      }
    }

    // copy of index parameters

    if (index_parameter) {
      i = -1;

      for (j = 0;j < nb_segment;j++) {
        switch (concatenation) {

        case false : {
          pindex_param = seq->index_parameter[j];
          break;
        }

        case true : {
          if ((j == 0) || (segment_begin[j][0] != segment_begin[j - 1][0])) {
            pindex_param = seq->index_parameter[++i];
          }
          break;
        }
        }

        cindex_param = index_parameter[segment_begin[j][0]] + segment_begin[j][1];
        for (k = 0;k < segment_length[j];k++) {
          *pindex_param++ = *cindex_param++;
        }
      }
    }

    if (index_parameter_distribution) {
      seq->build_index_parameter_frequency_distribution();
    }
    if (index_interval) {
      seq->index_interval_computation();
    }

    // copie des valeurs

    i = -1;
    for (j = 0;j < nb_segment;j++) {
      k = 0;
      for (m = 0;m < nb_variable;m++) {
        if (m != variable) {
          if ((type[m] != REAL_VALUE) && (type[m] != AUXILIARY)) {
            switch (concatenation) {

            case false : {
              pisequence = seq->int_sequence[j][k];
              break;
            }

            case true : {
              if ((j == 0) || (segment_begin[j][0] != segment_begin[j - 1][0])) {
                pisequence = seq->int_sequence[++i][k];
              }
              break;
            }
            }

            cisequence = int_sequence[segment_begin[j][0]][m] + segment_begin[j][1];
            for (n = 0;n < segment_length[j];n++) {
              *pisequence++ = *cisequence++;
            }
          }

          else {
            switch (concatenation) {

            case false : {
              prsequence = seq->real_sequence[j][k];
              break;
            }

            case true : {
              if ((j == 0) || (segment_begin[j][0] != segment_begin[j - 1][0])) {
                prsequence = seq->real_sequence[++i][k];
              }
              break;
            }
            }

            crsequence = real_sequence[segment_begin[j][0]][m] + segment_begin[j][1];
            for (n = 0;n < segment_length[j];n++) {
              *prsequence++ = *crsequence++;
            }
          }

          k++;
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);

      seq->build_marginal_frequency_distribution(i);
    }

    if (!keep) {
      delete [] selected_value;
    }

    delete [] segment_length;

    for (i = 0;i < nb_segment;i++) {
      delete [] segment_begin[i];
    }
    delete [] segment_begin;

    delete [] itype;

    if (concatenation) {
      delete [] sequence_length;
    }

#   ifdef MESSAGE
    if ((seq->index_parameter_type == TIME) && (seq->index_interval->variance > 0.)) {  // pour les suivis de croissance manguier
      double average_diff , individual_mean , global_mean , diff , individual_variance , global_variance;

      for (i = 0;i < seq->nb_variable;i++) {
        if (seq->type[i] == INT_VALUE) {
          average_diff = 0.;
          for (j = 0;j < seq->nb_sequence;j++) {
            average_diff += seq->int_sequence[j][i][seq->length[j] - 1] - seq->int_sequence[j][i][0];
          }
          average_diff /= seq->nb_sequence;

          global_mean = 0.;
          individual_variance = 0.;
          for (j = 0;j < seq->nb_sequence;j++) {
            individual_mean = 0.;
            for (k = 0;k < seq->length[j];k++) {
              individual_mean += seq->int_sequence[j][i][k];
            }
            global_mean += individual_mean;
            individual_mean /= seq->length[j];

            for (k = 0;k < seq->length[j];k++) {
              diff = seq->int_sequence[j][i][k] - individual_mean;
              individual_variance += diff * diff;
            }
          }
          global_mean /= seq->cumul_length;
          individual_variance /= (seq->cumul_length - 1);

          global_variance = 0.;
          for (j = 0;j < seq->nb_sequence;j++) {
            for (k = 0;k < seq->length[j];k++) {
              diff = seq->int_sequence[j][i][k] - global_mean;
              global_variance += diff * diff;
            }
          }
          global_variance /= (seq->cumul_length - 1);

          cout << "\n" << STAT_label[STATL_VARIABLE] << " " << i + 1 << " - "
               << "average difference: " <<  average_diff << " | "
               << "within-individual variance: " << individual_variance << " | "
               << "global variance: " << global_variance << " | "
               << individual_variance / global_variance;
        }
      }
      cout << endl;
    }
#   endif

  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Cumul des valeurs successives des sequences.
 *
 *  arguments : reference sur un objet StatError, indice de la variable.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::cumulate(StatError &error , int variable) const

{
  bool status = true;
  register int i , j , k , m;
  int inb_variable , *itype;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (variable != I_DEFAULT) {
    if ((variable < 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }
    else {
      variable--;
    }
  }

  for (i = 0;i < nb_variable;i++) {
    if (((variable == I_DEFAULT) || (variable == i)) &&
        ((type[i] != INT_VALUE) && (type[i] != STATE) && (type[i] != REAL_VALUE))) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (status) {
    if (variable == I_DEFAULT) {
      inb_variable = nb_variable;
      itype = type;
    }
    else {
      inb_variable = 1;
      itype = type + variable;
    }

    seq = new Sequences(nb_sequence , identifier , length , vertex_identifier ,
                        index_parameter_type , inb_variable , itype);

    // copy of index parameters

    if (index_parameter_distribution) {
      seq->index_parameter_distribution = new FrequencyDistribution(*index_parameter_distribution);
    }
    if (index_interval) {
      seq->index_interval = new FrequencyDistribution(*index_interval);
    }

    if (index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < (seq->index_parameter_type == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }
      }
    }

   // cumul des valeurs

    for (i = 0;i < nb_sequence;i++) {
      j = 0;
      for (k = 0;k < nb_variable;k++) {
        if ((variable == I_DEFAULT) || (variable == k)) {
          if (type[k] != REAL_VALUE) {
            seq->int_sequence[i][j][0] = int_sequence[i][k][0];
            for (m = 1;m < length[i];m++) {
              seq->int_sequence[i][j][m] = seq->int_sequence[i][j][m - 1] + int_sequence[i][k][m];
            }
          }

          else {
            seq->real_sequence[i][j][0] = real_sequence[i][k][0];
            for (m = 1;m < length[i];m++) {
              seq->real_sequence[i][j][m] = seq->real_sequence[i][j][m - 1] + real_sequence[i][k][m];
            }
          }

          j++;
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);

      seq->build_marginal_frequency_distribution(i);
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Differenciation au 1er ordre des sequences.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              premier element de la sequence garde ou pas.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::difference(StatError &error , int variable , bool first_element) const

{
  bool status = true;
  register int i , j , k , m;
  int inb_variable , *ilength , *itype , *pvertex_id , *cvertex_id ,
      *pindex_param , *cindex_param , *pisequence , *cisequence;
  double *prsequence , *crsequence;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (index_parameter_type == POSITION) {
    status = increasing_index_parameter_checking(error , true , SEQ_label[SEQL_SEQUENCE]);
  }

  if (variable != I_DEFAULT) {
    if ((variable < 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }
    else {
      variable--;
    }
  }

  for (i = 0;i < nb_variable;i++) {
    if (((variable == I_DEFAULT) || (variable == i)) &&
        ((type[i] != INT_VALUE) && (type[i] != STATE) && (type[i] != REAL_VALUE))) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if ((!first_element) && (length_distribution->offset < 2)) {
    status = false;
    ostringstream correction_message;
    correction_message << STAT_error[STATR_GREATER_THAN] << " " << 1;
    error.correction_update(SEQ_error[SEQR_MIN_SEQUENCE_LENGTH] , (correction_message.str()).c_str());
  }

  if (status) {
    if (first_element) {
      ilength = length;
    }
    else {
      ilength = new int[nb_sequence];
      for (i = 0;i < nb_sequence;i++) {
        ilength[i] = length[i] - 1;
      }
    }

    if (variable == I_DEFAULT) {
      inb_variable = nb_variable;
    }
    else {
      inb_variable = 1;
    }

    if ((index_parameter_type != IMPLICIT_TYPE) && ((index_interval->mean != 1.) ||
         (index_interval->variance > 0.))) {
      itype = new int[inb_variable];
      for (i = 0;i < inb_variable;i++) {
        itype[i] = REAL_VALUE;
      }
    }

    else {
      if (variable == I_DEFAULT) {
        itype = type;
      }
      else {
        itype = type + variable;
      }
    }

    seq = new Sequences(nb_sequence , identifier , ilength , vertex_identifier ,
                        index_parameter_type , inb_variable , itype , false);

    if (!first_element) {
      delete [] ilength;
    }
    if ((index_parameter_type != IMPLICIT_TYPE) && ((index_interval->mean != 1.) ||
         (index_interval->variance > 0.))) {
      delete [] itype;
    }

    // copy of vertex identifiers

    if (vertex_identifier) {
      for (i = 0;i < seq->nb_sequence;i++) {
        pvertex_id = seq->vertex_identifier[i];
        if (first_element) {
          cvertex_id = vertex_identifier[i];
        }
        else {
          cvertex_id = vertex_identifier[i] + 1;
        }

        for (j = 0;j < seq->length[i];j++) {
          *pvertex_id++ = *cvertex_id++;
        }
      }
    }

    // copy of index parameters

    if (index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        pindex_param = seq->index_parameter[i];
        if (first_element) {
          cindex_param = index_parameter[i];
        }
        else {
          cindex_param = index_parameter[i] + 1;
        }

        for (j = 0;j < (seq->index_parameter_type == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
          *pindex_param++ = *cindex_param++;
        }
      }
    }

    if (first_element) {
      if (index_parameter_distribution) {
        seq->index_parameter_distribution = new FrequencyDistribution(*index_parameter_distribution);
      }
      if (index_interval) {
        seq->index_interval = new FrequencyDistribution(*index_interval);
      }
    }

    else {
      seq->build_index_parameter_frequency_distribution();
      if (index_interval) {
        seq->index_interval_computation();
      }
    }

    // differenciation des sequences

    if ((index_parameter_type != IMPLICIT_TYPE) && ((index_interval->mean != 1.) ||
         (index_interval->variance > 0.))) {
      for (i = 0;i < nb_sequence;i++) {
        j = 0;
        for (k = 0;k < nb_variable;k++) {
          if ((variable == I_DEFAULT) || (variable == k)) {
            prsequence = seq->real_sequence[i][j];
            cindex_param = index_parameter[i];

            if (type[k] != REAL_VALUE) {
              cisequence = int_sequence[i][k];

              if (first_element) {
                if (*cindex_param > 0) {
                  *prsequence++ = (double)*cisequence / (double)*cindex_param;
                }
                else {
                  *prsequence++ = D_DEFAULT;
                }
              }

              for (m = 0;m < length[i] - 1;m++) {
                *prsequence++ = (double)(*(cisequence + 1) - *cisequence) /
                                (double)(*(cindex_param + 1) - *cindex_param);
                cindex_param++;
                cisequence++;
              }
            }

            else {
              crsequence = real_sequence[i][k];

              if (first_element) {
                if (*cindex_param > 0) {
                  *prsequence++ = *crsequence / *cindex_param;
                }
                else {
                  *prsequence++ = D_DEFAULT;
                }
              }

              for (m = 0;m < length[i] - 1;m++) {
                *prsequence++ = (*(crsequence + 1) - *crsequence) /
                                (*(cindex_param + 1) - *cindex_param);
                cindex_param++;
                crsequence++;
              }
            }

            j++;
          }
        }
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        j = 0;
        for (k = 0;k < nb_variable;k++) {
          if ((variable == I_DEFAULT) || (variable == k)) {
            if (type[k] != REAL_VALUE) {
              pisequence = seq->int_sequence[i][j];
              cisequence = int_sequence[i][k];

              if (first_element) {
                *pisequence++ = *cisequence;
              }
              for (m = 0;m < length[i] - 1;m++) {
                *pisequence++ = *(cisequence + 1) - *cisequence;
                cisequence++;
              }
            }

            else {
              prsequence = seq->real_sequence[i][j];
              crsequence = real_sequence[i][k];

              if (first_element) {
                *prsequence++ = *crsequence;
              }
              for (m = 0;m < length[i] - 1;m++) {
                *prsequence++ = *(crsequence + 1) - *crsequence;
                crsequence++;
              }
            }

            j++;
          }
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);

      seq->build_marginal_frequency_distribution(i);
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Computation of relative growth rates on the basis of cumulative dimensions.
 *
 *  arguments : reference on an object StatError, growth factor for
 *              computing the first relative growth rate.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::relative_growth_rate(StatError &error , double growth_factor) const

{
  bool status = true , begin;
  register int i , j , k;
  int *itype , *pindex_param , *cindex_param , *pisequence , *cisequence;
  double *prsequence , *crsequence;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (index_parameter_type == POSITION) {
    status = increasing_index_parameter_checking(error , true , SEQ_label[SEQL_SEQUENCE]);
  }

  for (i = 0;i < nb_variable;i++) {
    if ((type[i] != INT_VALUE) && (type[i] != STATE) && (type[i] != REAL_VALUE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }

/*    else {
      lstatus = increasing_sequence_checking(error , i , false , SEQ_label[SEQL_SEQUENCE] ,
                                             STAT_label[STATL_VALUE]);
      if (!lstatus) {
        status = false;
      }
    } */

    else if (min_value[i] < 0.) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_POSITIVE_MIN_VALUE];
      error.update((error_message.str()).c_str());
    }
  }

  if (status) {
    itype = new int[nb_variable];
    for (i = 0;i < nb_variable;i++) {
      itype[i] = REAL_VALUE;
    }

    seq = new Sequences(nb_sequence , identifier , length , vertex_identifier ,
                        index_parameter_type , nb_variable , itype);
    delete [] itype;

    // copy of index parameters

    if (index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < (seq->index_parameter_type == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }
      }
    }

    if (index_parameter_distribution) {
      seq->index_parameter_distribution = new FrequencyDistribution(*index_parameter_distribution);
    }
    if (index_interval) {
      seq->index_interval = new FrequencyDistribution(*index_interval);
    }

    // computaton of relative growth rate

    if ((index_parameter_type != IMPLICIT_TYPE) && ((index_interval->mean != 1.) ||
         (index_interval->variance > 0.))) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < nb_variable;j++) {
          prsequence = seq->real_sequence[i][j];
          cindex_param = index_parameter[i] + 1;
          begin = true;

          if (type[j] != REAL_VALUE) {
            cisequence = int_sequence[i][j] + 1;
            *prsequence++ = 0.;

            for (k = 1;k < length[i];k++) {
              if ((*cisequence > 0) && (*(cisequence - 1) > 0)) {
                *prsequence = (log(*cisequence) - log(*(cisequence - 1))) /
                              (*cindex_param - *(cindex_param - 1));

                if (begin) {
                  begin = false;
                  *(prsequence - 1) = *prsequence * growth_factor;
                }
                prsequence++;
              }

              else {
                *prsequence++ = 0.;
              }
              cindex_param++;
              cisequence++;
            }
          }

          else {
            crsequence = real_sequence[i][j] + 1;
            *prsequence++ = 0.;

            for (k = 1;k < length[i];k++) {
              if ((*crsequence > 0.) && (*(crsequence - 1) > 0.)) {
                *prsequence = (log(*crsequence) - log(*(crsequence - 1))) /
                              (*cindex_param - *(cindex_param - 1));

                if (begin) {
                  begin = false;
                  *(prsequence - 1) = *prsequence * growth_factor;
                }
                prsequence++;
              }

              else {
                *prsequence++ = 0.;
              }
              cindex_param++;
              crsequence++;
            }
          }
        }
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < nb_variable;j++) {
          prsequence = seq->real_sequence[i][j];
          begin = true;

          if (type[j] != REAL_VALUE) {
            cisequence = int_sequence[i][j] + 1;
            *prsequence++ = 0.;

            for (k = 1;k < length[i];k++) {
              if ((*cisequence > 0) && (*(cisequence - 1) > 0)) {
                *prsequence = log(*cisequence) - log(*(cisequence - 1));

                if (begin) {
                  begin = false;
                  *(prsequence - 1) = *prsequence * growth_factor;
                }
                prsequence++;
              }

              else {
                *prsequence++ = 0.;
              }
              cisequence++;
            }
          }

          else {
            crsequence = real_sequence[i][j] + 1;
            *prsequence++ = 0.;

            for (k = 1;k < length[i];k++) {
              if ((*crsequence > 0.) && (*(crsequence - 1) > 0.)) {
                *prsequence = log(*crsequence) - log(*(crsequence - 1));

                if (begin) {
                  begin = false;
                  *(prsequence - 1) = *prsequence * growth_factor;
                }
                prsequence++;
              }

              else {
                *prsequence++ = 0.;
              }
              crsequence++;
            }
          }
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);

      seq->build_marginal_frequency_distribution(i);
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Normalisation des sequences.
 *
 *  arguments : reference sur un objet StatError, indice de la variable.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::sequence_normalization(StatError &error , int variable) const

{
  bool status = true;
  register int i , j , k;
  int int_max , *itype;
  double real_max;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (variable != I_DEFAULT) {
    if ((variable < 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }
    else {
      variable--;
    }
  }

  for (i = 0;i < nb_variable;i++) {
    if ((variable == I_DEFAULT) || (variable == i)) {
      if ((type[i] != INT_VALUE) && (type[i] != STATE) && (type[i] != REAL_VALUE)) {
        status = false;
        ostringstream error_message , correction_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                      << STAT_error[STATR_VARIABLE_TYPE];
        correction_message << STAT_variable_word[INT_VALUE] << " or "
                           << STAT_variable_word[STATE] << " or "
                           << STAT_variable_word[REAL_VALUE];
        error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
      }

      if (min_value[i] < 0) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                      << STAT_error[STATR_POSITIVE_MIN_VALUE];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {
    itype = new int[nb_variable];

    if (variable == I_DEFAULT) {
      for (i = 0;i < nb_variable;i++) {
        itype[i] = REAL_VALUE;
      }
    }
    else {
      for (i = 0;i < nb_variable;i++) {
        itype[i] = type[i];
      }
      itype[variable] = REAL_VALUE;
    }

    seq = new Sequences(nb_sequence , identifier , length , vertex_identifier ,
                        index_parameter_type , nb_variable , itype);
    delete [] itype;

    // copy of index parameters

    if (index_parameter_distribution) {
      seq->index_parameter_distribution = new FrequencyDistribution(*index_parameter_distribution);
    }
    if (index_interval) {
      seq->index_interval = new FrequencyDistribution(*index_interval);
    }

    if (index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < (seq->index_parameter_type == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }
      }
    }

   // normalisation des sequences

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < nb_variable;j++) {
        if ((variable == I_DEFAULT) || (variable == j)) {
          if (type[j] != REAL_VALUE) {
            int_max = int_sequence[i][j][0];
            for (k = 1;k < length[i];k++) {
              if (int_sequence[i][j][k] >  int_max) {
                int_max = int_sequence[i][j][k];
              }
            }

            for (k = 0;k < length[i];k++) {
              seq->real_sequence[i][j][k] = (double)int_sequence[i][j][k] / (double)int_max;
            }
          }

          else {
            real_max = real_sequence[i][j][0];
            for (k = 1;k < length[i];k++) {
              if (real_sequence[i][j][k] >  real_max) {
                real_max = real_sequence[i][j][k];
              }
            }

            for (k = 0;k < length[i];k++) {
              seq->real_sequence[i][j][k] = real_sequence[i][j][k] / real_max;
            }
          }
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      if ((variable == I_DEFAULT) || (variable == i)) {
        seq->min_value_computation(i);
        seq->max_value[i] = 1.;

        seq->build_marginal_histogram(i);
      }

      else {
        seq->min_value[i] = min_value[i];
        seq->max_value[i] = max_value[i];

        if (marginal_distribution[i]) {
          seq->marginal_distribution[i] = new FrequencyDistribution(*marginal_distribution[i]);
        }
        if (marginal_histogram[i]) {
          seq->marginal_histogram[i] = new Histogram(*marginal_histogram[i]);
        }
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Filtrage de type moyenne mobile des sequences.
 *
 *  arguments : reference sur un objet StatError, demi-largeur du filtre,
 *              filtre, indice de la variable, debut/fin garde ou pas,
 *              tendance ou residus (par soustraction ou par division).
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::moving_average(StatError &error , int nb_point , double *filter ,
                                     int variable , bool begin_end , int output) const

{
  bool status = true;
  register int i , j , k , m , n;
  int inb_variable , *ilength , *itype , *pvertex_id , *cvertex_id ,
      *pindex_param , *cindex_param , *pisequence , *cisequence;
  double *prsequence , *crsequence , *ppoint;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((index_interval) && (index_interval->variance > 0.)) {
    status = false;
    error.update(SEQ_error[SEQR_UNEQUAL_INDEX_INTERVALS]);
  }

  if (variable != I_DEFAULT) {
    if ((variable < 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }
    else {
      variable--;
    }
  }

  for (i = 0;i < nb_variable;i++) {
    if (((variable == I_DEFAULT) || (variable == i)) &&
        ((type[i] != INT_VALUE) && (type[i] != STATE) && (type[i] != REAL_VALUE))) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if ((!begin_end) && (length_distribution->offset < 2 * nb_point + 1)) {
    status = false;
    ostringstream correction_message;
    correction_message << STAT_error[STATR_GREATER_THAN] << " " << 2 * nb_point;
    error.correction_update(SEQ_error[SEQR_MIN_SEQUENCE_LENGTH] , (correction_message.str()).c_str());
  }

  if (status) {
    if (begin_end) {
      ilength = length;
    }
    else {
      ilength = new int[nb_sequence];
      for (i = 0;i < nb_sequence;i++) {
        ilength[i] = length[i] - 2 * nb_point;
      }
    }

    if (variable == I_DEFAULT) {
      inb_variable = nb_variable;
    }
    else {
      inb_variable = 1;
    }

    if (output == SEQUENCE) {
      inb_variable *= 2;
    }

    itype = new int[inb_variable];

    if (output == SEQUENCE) {
      i = 0;
      for (j = 0;j < nb_variable;j++) {
        if ((variable == I_DEFAULT) || (variable == j)) {
          itype[i++] = type[j];
          itype[i++] = AUXILIARY;
        }
      }
    }

    else {
      for (i = 0;i < inb_variable;i++) {
        itype[i] = REAL_VALUE;
      }
    }

    seq = new Sequences(nb_sequence , identifier , ilength , vertex_identifier ,
                        index_parameter_type , inb_variable , itype , false);

    if (!begin_end) {
      delete [] ilength;
    }
    delete [] itype;

    // copy of vertex identifiers

    if (vertex_identifier) {
      for (i = 0;i < seq->nb_sequence;i++) {
        pvertex_id = seq->vertex_identifier[i];
        if (begin_end) {
          cvertex_id = vertex_identifier[i];
        }
        else {
          cvertex_id = vertex_identifier[i] + nb_point;
        }

        for (j = 0;j < seq->length[i];j++) {
          *pvertex_id++ = *cvertex_id++;
        }
      }
    }

    // copy of index parameters

    if (index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        pindex_param = seq->index_parameter[i];
        if (begin_end) {
          cindex_param = index_parameter[i];
        }
        else {
          cindex_param = index_parameter[i] + nb_point;
        }

        for (j = 0;j < (seq->index_parameter_type == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
          *pindex_param++ = *cindex_param++;
        }
      }
    }

    if (begin_end) {
      if (index_parameter_distribution) {
        seq->index_parameter_distribution = new FrequencyDistribution(*index_parameter_distribution);
      }
      if (index_interval) {
        seq->index_interval = new FrequencyDistribution(*index_interval);
      }
    }

    else {
      seq->build_index_parameter_frequency_distribution();
      if (index_interval) {
        seq->index_interval_computation();
      }
    }

    // filtrage de type moyenne mobile

    for (i = 0;i < nb_sequence;i++) {
      j = 0;
      for (k = 0;k < nb_variable;k++) {
        if ((variable == I_DEFAULT) || (variable == k)) {
          prsequence = seq->real_sequence[i][output == SEQUENCE ? j + 1 : j];

          if (type[k] != REAL_VALUE) {
            if (begin_end) {
              for (m = 0;m < MIN(nb_point , length[i]);m++) {
                cisequence = int_sequence[i][k];
                ppoint = filter;
                *prsequence = 0.;
                for (n = 0;n < 2 * nb_point + 1;n++) {
                  *prsequence += *cisequence * *ppoint++;
                  if ((m - nb_point + n >= 0) && (m - nb_point + n < length[i] - 1)) {
                    cisequence++;
                  }
                }
                prsequence++;
              }
            }

            for (m = nb_point;m < length[i] - nb_point;m++) {
              cisequence = int_sequence[i][k] + m - nb_point;
              ppoint = filter;
              *prsequence = 0.;
              for (n = 0;n < 2 * nb_point + 1;n++) {
                *prsequence += *cisequence++ * *ppoint++;
              }
              prsequence++;
            }

            if (begin_end) {
              for (m = MAX(length[i] - nb_point , nb_point);m < length[i];m++) {
                cisequence = int_sequence[i][k] + m - nb_point;
                ppoint = filter;
                *prsequence = 0.;
                for (n = 0;n < 2 * nb_point + 1;n++) {
                  *prsequence += *cisequence * *ppoint++;
                  if ((m - nb_point + n >= 0) && (m - nb_point + n < length[i] - 1)) {
                    cisequence++;
                  }
                }
                prsequence++;
              }
            }

            if (begin_end) {
              cisequence = int_sequence[i][k];
            }
            else {
              cisequence = int_sequence[i][k] + nb_point;
            }

            switch (output) {

            case SEQUENCE : {
              pisequence = seq->int_sequence[i][j];
              for (m = 0;m < seq->length[i];m++) {
                *pisequence++ = *cisequence++;
              }
              break;
            }

            case SUBTRACTION_RESIDUAL : {
              prsequence = seq->real_sequence[i][j];
              for (m = 0;m < seq->length[i];m++) {
                *prsequence = *cisequence++ - *prsequence;
                prsequence++;
              }
              break;
            }

            case DIVISION_RESIDUAL : {
              prsequence = seq->real_sequence[i][j];
              for (m = 0;m < seq->length[i];m++) {
                if (*prsequence != 0.) {
                  *prsequence = *cisequence / *prsequence;
                }
                prsequence++;
                cisequence++;
              }
              break;
            }
            }
          }

          else {
            if (begin_end) {
              for (m = 0;m < MIN(nb_point , length[i]);m++) {
                crsequence = real_sequence[i][k];
                ppoint = filter;
                *prsequence = 0.;
                for (n = 0;n < 2 * nb_point + 1;n++) {
                  *prsequence += *crsequence * *ppoint++;
                  if ((m - nb_point + n >= 0) && (m - nb_point + n < length[i] - 1)) {
                    crsequence++;
                  }
                }
                prsequence++;
              }
            }

            for (m = nb_point;m < length[i] - nb_point;m++) {
              crsequence = real_sequence[i][k] + m - nb_point;
              ppoint = filter;
              *prsequence = 0.;
              for (n = 0;n < 2 * nb_point + 1;n++) {
                *prsequence += *crsequence++ * *ppoint++;
              }
              prsequence++;
            }

            if (begin_end) {
              for (m = MAX(length[i] - nb_point , nb_point);m < length[i];m++) {
                crsequence = real_sequence[i][k] + m - nb_point;
                ppoint = filter;
                *prsequence = 0.;
                for (n = 0;n < 2 * nb_point + 1;n++) {
                  *prsequence += *crsequence * *ppoint++;
                  if ((m - nb_point + n >= 0) && (m - nb_point + n < length[i] - 1)) {
                    crsequence++;
                  }
                }
                prsequence++;
              }
            }

            prsequence = seq->real_sequence[i][j];
            if (begin_end) {
              crsequence = real_sequence[i][k];
            }
            else {
              crsequence = real_sequence[i][k] + nb_point;
            }

            switch (output) {

            case SEQUENCE : {
              for (m = 0;m < seq->length[i];m++) {
                *prsequence++ = *crsequence++;
              }
              break;
            }

            case SUBTRACTION_RESIDUAL : {
              for (m = 0;m < seq->length[i];m++) {
                *prsequence = *crsequence++ - *prsequence;
                prsequence++;
              }
              break;
            }

            case DIVISION_RESIDUAL : {
              for (m = 0;m < seq->length[i];m++) {
                if (*prsequence != 0.) {
                  *prsequence = *crsequence / *prsequence;
                }
                prsequence++;
                crsequence++;
              }
              break;
            }
            }
          }

          if (output == SEQUENCE) {
            j++;
          }
          j++;
        }
      }
    }

    if (output == SEQUENCE) {
      i = 0;

      if (begin_end) {
        for (j = 0;j < nb_variable;j++) {
          if ((variable == I_DEFAULT) || (variable == j)) {
            seq->min_value[i] = min_value[j];
            seq->max_value[i] = max_value[j];

            if (marginal_distribution[j]) {
              seq->marginal_distribution[i] = new FrequencyDistribution(*marginal_distribution[j]);
            }
            if (marginal_histogram[j]) {
              seq->marginal_histogram[i] = new Histogram(*marginal_histogram[j]);
            }
            i++;

            seq->min_value_computation(i);
            seq->max_value_computation(i);
            i++;
          }
        }
      }

      else {
        for (j = 0;j < nb_variable;j++) {
          if ((variable == I_DEFAULT) || (variable == j)) {
            seq->min_value_computation(i);
            seq->max_value_computation(i);

            seq->build_marginal_frequency_distribution(i);
            i++;

            seq->min_value_computation(i);
            seq->max_value_computation(i);
            i++;
          }
        }
      }
    }

    else {
      for (i = 0;i < seq->nb_variable;i++) {
        seq->min_value_computation(i);
        seq->max_value_computation(i);

        seq->build_marginal_histogram(i);
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Filtrage de type moyenne mobile des sequences.
 *
 *  arguments : reference sur un objet StatError, loi symmetrique,
 *              indice de la variable, debut/fin supprime ou pas,
 *              tendance ou residus (par soustraction ou par division).
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::moving_average(StatError &error , const Distribution &dist ,
                                     int variable , bool begin_end , int output) const

{
  bool status = true;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((dist.offset != 0) || ((dist.nb_value - dist.offset) % 2 == 0)) {
    status = false;
    error.correction_update(STAT_error[STATR_NB_VALUE] , STAT_error[STATR_ODD]);
  }
  if (fabs(dist.skewness_computation()) > SKEWNESS_ROUNDNESS) {
    status = false;
    error.update(STAT_error[STATR_NON_SYMMETRICAL_DISTRIBUTION]);
  }
  if (dist.complement > 0.) {
    status = false;
    error.update(STAT_error[STATR_UNPROPER_DISTRIBUTION]);
  }

  if (status) {
    seq = moving_average(error , dist.nb_value / 2 , dist.mass , variable ,
                         begin_end , output);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des sequences de moyennes et d'ecart-types.
 *
 *  arguments : reference sur un objet StatError, path,
 *              tailles d'echantillons, flag calcul des ecart-types,
 *              sortie (sequences, residus ou residus standardisees).
 *
 *--------------------------------------------------------------*/

bool Sequences::pointwise_average_ascii_print(StatError &error , const char *path ,
                                              int *size , bool standard_deviation ,
                                              int output) const

{
  bool status;
  register int i , j , k , m;
  int buff , inb_sequence , *width;
  long old_adjust;
  double standard_normal_value , half_confidence_interval , *t_value;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    old_adjust = out_file.flags(ios::adjustfield);

    if (standard_deviation) {

#     ifdef MESSAGE
      normal normal_dist;
      standard_normal_value = quantile(complement(normal_dist , 0.025));
      cout << "\nTEST: " << standard_normal_value;
#     endif

      t_value = new double[length[0]];
      for (i = 0;i < length[0];i++) {
        if (size[i] > 1) {
//          t_value[i] = t_value_computation(false , size[i] - 1 , 0.05);
          students_t students_dist(size[i] - 1);
          t_value[i] = quantile(complement(students_dist , 0.025));
        }
      }

#     ifdef MESSAGE
      cout << " | " << t_value[0] << " " << size[0] << endl;
#     endif

    }

    else {
      t_value = NULL;
    }

    // calcul des largeurs des colonnes

    inb_sequence = nb_sequence;
    if (standard_deviation) {
      inb_sequence--;
    }

    width = new int[2 * nb_variable + 2];

    for (i = 0;i < nb_variable;i++) {
      width[i] = column_width(length[0] , real_sequence[0][i]);
      if (standard_deviation) {
        buff = column_width(length[nb_sequence - 1] , real_sequence[nb_sequence - 1][i]);
        if (buff > width[i]) {
          width[i] = buff;
        }
      }
      width[i] += ASCII_SPACE;
    }

    if (index_parameter) {
      width[nb_variable] = column_width(index_parameter_distribution->nb_value - 1);
    }
    else {
      width[nb_variable] = column_width(max_length);
    }

    width[nb_variable + 1] = column_width(nb_sequence) + ASCII_SPACE;

    if (nb_sequence < POINTWISE_AVERAGE_NB_SEQUENCE) {
      for (i = 0;i < nb_variable;i++) {
        width[nb_variable + i + 2] = 0;
        for (j = 1;j < inb_sequence;j++) {
          buff = column_width(length[j] , real_sequence[j][i]);
          if (buff > width[nb_variable + i + 2]) {
            width[nb_variable + i + 2] = buff;
          }
        }
        width[nb_variable + i + 2] += ASCII_SPACE;
      }
    }

    out_file.setf(ios::right , ios::adjustfield);

    switch (output) {
    case SUBTRACTION_RESIDUAL :
      out_file << STAT_label[STATL_RESIDUAL] << "\n" << endl;
      break;
    case STANDARDIZED_RESIDUAL :
      out_file << STAT_label[STATL_STANDARDIZED_RESIDUAL] << "\n" << endl;
      break;
    }

    for (i = 0;i < nb_variable;i++) {
      out_file << STAT_label[STATL_VARIABLE] << " " << i + 1 << endl;

      if (index_parameter_type == TIME) {
        out_file << "\n" << SEQ_label[SEQL_TIME];
      }
      else {
        out_file << "\n" << SEQ_label[SEQL_INDEX];
      }

      out_file << "   " << STAT_label[STATL_MEAN];
      if (standard_deviation) {
        out_file << "   " << STAT_label[STATL_MEAN_CONFIDENCE_INTERVAL]
                 << "   " << STAT_label[STATL_STANDARD_DEVIATION];
      }
      out_file << "   " << STAT_label[STATL_SAMPLE_SIZE];

      if (nb_sequence < POINTWISE_AVERAGE_NB_SEQUENCE) {
        out_file << "   ";
        for (j = 1;j < inb_sequence;j++) {
          out_file << "   " << SEQ_label[SEQL_SEQUENCE] << " " << identifier[j];
        }
      }
      out_file << endl;

      for (j = 0;j < length[0];j++) {
        out_file << setw(width[nb_variable]) << (index_parameter ? index_parameter[0][j] : j) << "  "
                 << setw(width[i]) << real_sequence[0][i][j];

        if (standard_deviation) {
          if (size[j] > 1) {
//            half_confidence_interval = standard_normal_value * real_sequence[nb_sequence - 1][i][j] / sqrt((double)size[j]);
            half_confidence_interval = t_value[j] * real_sequence[nb_sequence - 1][i][j] / sqrt((double)size[j]);
            out_file << setw(width[i]) << real_sequence[0][i][j] - half_confidence_interval
                     << setw(width[i]) << real_sequence[0][i][j] + half_confidence_interval;
          }

          else {
            out_file << setw(width[i]) << " "
                     << setw(width[i]) << " ";
          }

          out_file << setw(width[i]) << real_sequence[nb_sequence - 1][i][j];
        }

        out_file << setw(width[nb_variable + 1]) << size[j];
 
        if (inb_sequence - 1 < POINTWISE_AVERAGE_NB_SEQUENCE) {
          out_file << "   ";

          if (index_parameter) {
            for (k = 1;k < inb_sequence;k++) {
              for (m = 0;m < length[k];m++) {
                if (index_parameter[k][m] == index_parameter[0][j]) {
                  out_file << setw(width[nb_variable + i + 2]) << real_sequence[k][i][m];
                  break;
                }
              }

              if (m == length[k]) {
                out_file << setw(width[nb_variable + i + 2]) << " ";
              }
            }
          }

          else {
            for (k = 1;k < inb_sequence;k++) {
              if (j < length[k]) {
                out_file << setw(width[nb_variable + i + 2]) << real_sequence[k][i][j];
              }
              else {
                out_file << setw(width[nb_variable + i + 2]) << " ";
              }
            }
          }
        }

        out_file << endl;
      }

      out_file << endl;
    }

    delete [] width;
    delete [] t_value;

    out_file.setf((FMTFLAGS)old_adjust , ios::adjustfield);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des sequences de moyennes et d'ecart-types au format tableur.
 *
 *  arguments : reference sur un objet StatError, path,
 *              tailles d'echantillons, flag calcul des ecart-types,
 *              sortie (sequences, residus ou residus standardisees).
 *
 *--------------------------------------------------------------*/

bool Sequences::pointwise_average_spreadsheet_print(StatError &error , const char *path ,
                                                    int *size , bool standard_deviation ,
                                                    int output) const

{
  bool status;
  register int i , j , k , m;
  int inb_sequence;
  double standard_normal_value , half_confidence_interval , *t_value;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    if (standard_deviation) {
//      normal normal_dist;
//      standard_normal_value = quantile(complement(normal_dist , 0.025));

      t_value = new double[length[0]];
      for (i = 0;i < length[0];i++) {
        if (size[i] > 1) {
//          t_value[i] = t_value_computation(false , size[i] - 1 , 0.05);
          students_t students_dist(size[i] - 1);
          t_value[i] = quantile(complement(students_dist , 0.025));
        }
      }
    }

    else {
      t_value = NULL;
    }

    switch (output) {
    case SUBTRACTION_RESIDUAL :
      out_file << STAT_label[STATL_RESIDUAL] << "\n" << endl;
      break;
    case STANDARDIZED_RESIDUAL :
      out_file << STAT_label[STATL_STANDARDIZED_RESIDUAL] << "\n" << endl;
      break;
    }

    for (i = 0;i < nb_variable;i++) {
      out_file << STAT_label[STATL_VARIABLE] << "\t" << i + 1 << endl;

      if (index_parameter_type == TIME) {
        out_file << "\n" << SEQ_label[SEQL_TIME];
      }
      else {
        out_file << "\n" << SEQ_label[SEQL_INDEX];
      }

      out_file << "\t" << STAT_label[STATL_MEAN];
      if (standard_deviation) {
        out_file << "\t" << STAT_label[STATL_MEAN_CONFIDENCE_INTERVAL]
                 << "\t\t" << STAT_label[STATL_STANDARD_DEVIATION];
      }
      out_file << "\t" << STAT_label[STATL_SAMPLE_SIZE];

      if (nb_sequence < POINTWISE_AVERAGE_NB_SEQUENCE) {
        inb_sequence = nb_sequence;
        if (standard_deviation) {
          inb_sequence--;
        }

        out_file << "\t";
        for (j = 1;j < inb_sequence;j++) {
          out_file << "\t" << SEQ_label[SEQL_SEQUENCE] << " " << identifier[j];
        }
      }
      out_file << endl;

      for (j = 0;j < length[0];j++) {
        out_file << (index_parameter ? index_parameter[0][j] : j)
                 << "\t" << real_sequence[0][i][j];

        if (standard_deviation) {
          if (size[j] > 1) {
//            half_confidence_interval = standard_normal_value * real_sequence[nb_sequence - 1][i][j] / sqrt((double)size[j]);
            half_confidence_interval = t_value[j] * real_sequence[nb_sequence - 1][i][j] / sqrt((double)size[j]);
            out_file << "\t" << real_sequence[0][i][j] - half_confidence_interval
                     << "\t" << real_sequence[0][i][j] + half_confidence_interval;
          }

          else {
            out_file << "\t\t";
          }

          out_file << "\t" << real_sequence[nb_sequence - 1][i][j];
        }

        out_file << "\t" << size[j];
 
        if (inb_sequence - 1 < POINTWISE_AVERAGE_NB_SEQUENCE) {
          out_file << "\t";

          if (index_parameter) {
            for (k = 1;k < inb_sequence;k++) {
              out_file << "\t";
              for (m = 0;m < length[k];m++) {
                if (index_parameter[k][m] == index_parameter[0][j]) {
                  out_file << real_sequence[k][i][m];
                  break;
                }
              }
            }
          }

          else {
            for (k = 1;k < inb_sequence;k++) {
              out_file << "\t";
              if (j < length[k]) {
                out_file << real_sequence[k][i][j];
              }
            }
          }
        }

        out_file << endl;
      }

      out_file << endl;
    }

    delete [] t_value;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la sequence des moyennes (et des ecart-types).
 *
 *  arguments : reference sur un objet StatError, flag donnees circulaires,
 *              flag calcul des ecart-types,
 *              sortie (sequences, residus ou residus standardisees),
 *              path, format ('a' : ASCII / 's' : Spreadsheet).
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::pointwise_average(StatError &error , bool circular ,
                                        bool standard_deviation , int output ,
                                        const char *path , char format) const

{
  bool status = true;
  register int i , j , k;
  int inb_sequence , min_identifier , max_identifier , unit , *iidentifier ,
      *ilength , *itype , *pindex_param , *size;
  double diff , *prsequence , *pmean , *pstandard_deviation , *pmean_direction1 ,
         *pmean_direction2 , ***mean_direction;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (nb_sequence == 1) {
    status = false;
    error.update(SEQ_error[SEQR_SINGLE_SEQUENCE]);
  }
  if ((index_parameter_type != IMPLICIT_TYPE) && (index_parameter_type != TIME)) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }

  for (i = 0;i < nb_variable;i++) {
    if ((type[i] != INT_VALUE) && (type[i] != STATE) && (type[i] != REAL_VALUE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (status) {
    if ((output == STANDARDIZED_RESIDUAL) && (!standard_deviation)) {
      standard_deviation = true;
    }

    inb_sequence = nb_sequence + (standard_deviation ? 2 : 1);

    iidentifier = new int[inb_sequence];

    min_identifier = identifier[0];
    for (i = 0;i < nb_sequence;i++) {
      if (identifier[i] < min_identifier) {
        min_identifier = identifier[i];
      }

      iidentifier[i + 1] = identifier[i];
    }

    iidentifier[0] = min_identifier - 1;

    if (standard_deviation) {
      max_identifier = identifier[nb_sequence - 1];
      for (i = 0;i < nb_sequence - 1;i++) {
        if (identifier[i] > max_identifier) {
          max_identifier = identifier[i];
        }
      }

      iidentifier[inb_sequence - 1] = max_identifier + 1;
    }

    ilength = new int[inb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      ilength[i + 1] = length[i];
    }

    if (index_parameter) {
      ilength[0] = 0;
      for (i = index_parameter_distribution->offset;i < index_parameter_distribution->nb_value;i++) {
        if (index_parameter_distribution->frequency[i] > 0) {
          ilength[0]++;
        }
      }
    }

    else {
      ilength[0] = max_length;
    }

    if (standard_deviation) {
      ilength[inb_sequence - 1] = ilength[0];
    }

    itype = new int[nb_variable];
    for (i = 0;i < nb_variable;i++) {
      itype[i] = REAL_VALUE;
    }

    seq = new Sequences(inb_sequence , iidentifier , ilength , NULL ,
                        index_parameter_type , nb_variable , itype);

    delete [] iidentifier;
    delete [] ilength;
    delete [] itype;

    if (index_parameter) {
      pindex_param = seq->index_parameter[0];
      for (i = index_parameter_distribution->offset;i < index_parameter_distribution->nb_value;i++) {
        if (index_parameter_distribution->frequency[i] > 0) {
          *pindex_param++ = i;
        }
      }

      if (standard_deviation) {
        for (i = 0;i < seq->length[seq->nb_sequence - 1];i++) {
          seq->index_parameter[seq->nb_sequence - 1][i] = seq->index_parameter[0][i];
        }
      }

      // copy of index parameters

      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          seq->index_parameter[i + 1][j] = index_parameter[i][j];
        }
      }

      seq->build_index_parameter_frequency_distribution();
      seq->index_interval_computation();
    }

    // calcul des tailles d'echantillons

    size = new int[seq->length[0]];

    if (index_parameter) {
      pindex_param = seq->index_parameter[0];
      i = 0;
      for (j = index_parameter_distribution->offset;j < index_parameter_distribution->nb_value;j++) {
        if (index_parameter_distribution->frequency[j] > 0) {
          size[i++] = index_parameter_distribution->frequency[j];
        }
      }
    }

    else {
      size[0] = nb_sequence;
      for (i = 1;i < max_length;i++) {
        size[i] = size[i - 1] - length_distribution->frequency[i];
      }
    }

    switch (circular) {

    case false : {

      // calcul des moyennes

      for (i = 0;i < seq->nb_variable;i++) {
        for (j = 0;j < seq->length[0];j++) {
          seq->real_sequence[0][i][j] = 0.;
        }
      }

      if (index_parameter) {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < nb_variable;j++) {
            pindex_param = seq->index_parameter[0];
            prsequence = seq->real_sequence[0][j];

            if (type[j] != REAL_VALUE) {
              for (k = 0;k < length[i];k++) {
                while (*pindex_param < index_parameter[i][k]) {
                  pindex_param++;
                  prsequence++;
                }
                pindex_param++;
                *prsequence++ += int_sequence[i][j][k];
              }
            }

            else {
              for (k = 0;k < length[i];k++) {
                while (*pindex_param < index_parameter[i][k]) {
                  pindex_param++;
                  prsequence++;
                }
                pindex_param++;
                *prsequence++ += real_sequence[i][j][k];
              }
            }
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < nb_variable;j++) {
            if (type[j] != REAL_VALUE) {
              for (k = 0;k < length[i];k++) {
                seq->real_sequence[0][j][k] += int_sequence[i][j][k];
              }
            }

            else {
              for (k = 0;k < length[i];k++) {
                seq->real_sequence[0][j][k] += real_sequence[i][j][k];
              }
            }
          }
        }
      }

      for (i = 0;i < seq->nb_variable;i++) {
        for (j = 0;j < seq->length[0];j++) {
          seq->real_sequence[0][i][j] /= size[j];
        }
      }

      // calcul des ecart-types

      if (standard_deviation) {
        for (i = 0;i < seq->nb_variable;i++) {
          for (j = 0;j < seq->length[seq->nb_sequence - 1];j++) {
            seq->real_sequence[seq->nb_sequence - 1][i][j] = 0.;
          }
        }

        if (index_parameter) {
          for (i = 0;i < nb_sequence;i++) {
            for (j = 0;j < nb_variable;j++) {
              pindex_param = seq->index_parameter[seq->nb_sequence - 1];
              prsequence = seq->real_sequence[seq->nb_sequence - 1][j];
              pmean = seq->real_sequence[0][j];

              if (type[j] != REAL_VALUE) {
                for (k = 0;k < length[i];k++) {
                  while (*pindex_param < index_parameter[i][k]) {
                    pindex_param++;
                    prsequence++;
                    pmean++;
                  }
                  pindex_param++;
                  diff = int_sequence[i][j][k] - *pmean++;
                  *prsequence++ += diff * diff;
                }
              }

              else {
                for (k = 0;k < length[i];k++) {
                  while (*pindex_param < index_parameter[i][k]) {
                    pindex_param++;
                    prsequence++;
                    pmean++;
                  }
                  pindex_param++;
                  diff = real_sequence[i][j][k] - *pmean++;
                  *prsequence++ += diff * diff;
                }
              }
            }
          }
        }

        else {
          for (i = 0;i < nb_sequence;i++) {
            for (j = 0;j < nb_variable;j++) {
              if (type[j] != REAL_VALUE) {
                for (k = 0;k < length[i];k++) {
                  diff = int_sequence[i][j][k] - seq->real_sequence[0][j][k];
                  seq->real_sequence[seq->nb_sequence - 1][j][k] += diff * diff;
                }
              }

              else {
                for (k = 0;k < length[i];k++) {
                  diff = real_sequence[i][j][k] - seq->real_sequence[0][j][k];
                  seq->real_sequence[seq->nb_sequence - 1][j][k] += diff * diff;
                }
              }
            }
          }
        }

        for (i = 0;i < seq->nb_variable;i++) {
          for (j = 0;j < seq->length[seq->nb_sequence - 1];j++) {
            if (size[j] > 1) {
              seq->real_sequence[seq->nb_sequence - 1][i][j] = sqrt(seq->real_sequence[seq->nb_sequence - 1][i][j] /
                                                                    (size[j] - 1));
            }
          }
        }
      }
      break;
    }

    case true : {

      // choix de l'unite

      unit = RADIAN;
      for (i = 0;i < nb_variable;i++) {
        if (max_value[i] - min_value[i] > 2 * M_PI) {
          unit = DEGREE;
          break;
        }
      }

      // calcul des directions moyennes

      mean_direction = new double**[seq->nb_variable];
      for (i = 0;i < seq->nb_variable;i++) {
        mean_direction[i] = new double*[3];
        for (j = 0;j < 3;j++) {
          mean_direction[i][j] = new double[seq->length[0]];
        }
      }

      for (i = 0;i < seq->nb_variable;i++) {
        for (j = 0;j < seq->length[0];j++) {
          mean_direction[i][0][j] = 0.;
          mean_direction[i][1][j] = 0.;
        }
      }

      if (index_parameter) {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < nb_variable;j++) {
            pindex_param = seq->index_parameter[0];
            pmean_direction1 = mean_direction[j][0];
            pmean_direction2 = mean_direction[j][1];

            if (type[j] != REAL_VALUE) {
              for (k = 0;k < length[i];k++) {
                while (*pindex_param < index_parameter[i][k]) {
                  pindex_param++;
                  pmean_direction1++;
                  pmean_direction2++;
                }
                pindex_param++;
                *pmean_direction1++ += cos(int_sequence[i][j][k] * M_PI / 180);
                *pmean_direction2++ += sin(int_sequence[i][j][k] * M_PI / 180);
              }
            }

            else {
              switch (unit) {

              case DEGREE : {
                for (k = 0;k < length[i];k++) {
                  while (*pindex_param < index_parameter[i][k]) {
                    pindex_param++;
                    pmean_direction1++;
                    pmean_direction2++;
                  }
                  pindex_param++;
                  *pmean_direction1++ += cos(real_sequence[i][j][k] * M_PI / 180);
                  *pmean_direction2++ += sin(real_sequence[i][j][k] * M_PI / 180);
                }
                break;
              }

              case RADIAN : {
                for (k = 0;k < length[i];k++) {
                  while (*pindex_param < index_parameter[i][k]) {
                    pindex_param++;
                    pmean_direction1++;
                    pmean_direction2++;
                  }
                  pindex_param++;
                  *pmean_direction1++ += cos(real_sequence[i][j][k]);
                  *pmean_direction2++ += sin(real_sequence[i][j][k]);
                }
                break;
              }
              }
            }
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < nb_variable;j++) {
            if (type[j] != REAL_VALUE) {
              for (k = 0;k < length[i];k++) {
                mean_direction[j][0][k] += cos(int_sequence[i][j][k] * M_PI / 180);
                mean_direction[j][1][k] += sin(int_sequence[i][j][k] * M_PI / 180);
              }
            }

            else {
              switch (unit) {

              case DEGREE : {
                for (k = 0;k < length[i];k++) {
                  mean_direction[j][0][k] += cos(real_sequence[i][j][k] * M_PI / 180);
                  mean_direction[j][1][k] += sin(real_sequence[i][j][k] * M_PI / 180);
                }
                break;
              }

              case RADIAN : {
                for (k = 0;k < length[i];k++) {
                  mean_direction[j][0][k] += cos(real_sequence[i][j][k]);
                  mean_direction[j][1][k] += sin(real_sequence[i][j][k]);
                }
                break;
              }
              }
            }
          }
        }
      }

      for (i = 0;i < seq->nb_variable;i++) {
        for (j = 0;j < seq->length[0];j++) {
          mean_direction[i][0][j] /= size[j];
          mean_direction[i][1][j] /= size[j];

          mean_direction[i][2][j] = sqrt(mean_direction[i][0][j] * mean_direction[i][0][j] +
                                         mean_direction[i][1][j] * mean_direction[i][1][j]);

          if (mean_direction[i][2][j] > 0.) {
            seq->real_sequence[0][i][j] = atan(mean_direction[i][1][j] / mean_direction[i][0][j]);

            if (mean_direction[i][0][j] < 0.) {
              seq->real_sequence[0][i][j] += M_PI;
            }
            if (unit == DEGREE) {
              seq->real_sequence[0][i][j] *= (180 / M_PI);
            }
          }

          else {
            seq->real_sequence[0][i][j] = D_DEFAULT;
          }
        }
      }

      // calcul des ecart-types

      if (standard_deviation) {
        for (i = 0;i < seq->nb_variable;i++) {
          for (j = 0;j < seq->length[0];j++) {
            if (mean_direction[i][2][j] > 0.) {
              seq->real_sequence[seq->nb_sequence - 1][i][j] = sqrt(-2 * log(mean_direction[i][2][j]));
              if (unit == DEGREE) {
                seq->real_sequence[seq->nb_sequence - 1][i][j] *= (180 / M_PI);
              }
            }

            else {
              seq->real_sequence[seq->nb_sequence - 1][i][j] = D_DEFAULT;
            }
          }
        }
      }

      for (i = 0;i < seq->nb_variable;i++) {
        for (j = 0;j < 3;j++) {
          delete [] mean_direction[i][j];
        }
        delete [] mean_direction[i];
      }
      delete [] mean_direction;
      break;
    }
    }

    switch (output) {

    // copie des sequences

    case SEQUENCE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < nb_variable;j++) {
          if (type[j] != REAL_VALUE) {
            for (k = 0;k < length[i];k++) {
              seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k];
            }
          }

          else {
            for (k = 0;k < length[i];k++) {
              seq->real_sequence[i + 1][j][k] = real_sequence[i][j][k];
            }
          }
        }
      }
      break;
    }

    // calcul des residus

    case SUBTRACTION_RESIDUAL : {
      switch (circular) {

      case false : {
        if (index_parameter) {
          for (i = 0;i < nb_sequence;i++) {
            for (j = 0;j < nb_variable;j++) {
              pindex_param = seq->index_parameter[0];
              pmean = seq->real_sequence[0][j];

              if (type[j] != REAL_VALUE) {
                for (k = 0;k < length[i];k++) {
                  while (*pindex_param < index_parameter[i][k]) {
                    pindex_param++;
                    pmean++;
                  }
                  pindex_param++;
                  seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k] - *pmean++;
                }
              }

              else {
                for (k = 0;k < length[i];k++) {
                  while (*pindex_param < index_parameter[i][k]) {
                    pindex_param++;
                    pmean++;
                  }
                  pindex_param++;
                  seq->real_sequence[i + 1][j][k] = real_sequence[i][j][k] - *pmean++;
                }
              }
            }
          }
        }

        else {
          for (i = 0;i < nb_sequence;i++) {
            for (j = 0;j < nb_variable;j++) {
              if (type[j] != REAL_VALUE) {
                for (k = 0;k < length[i];k++) {
                  seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k] - seq->real_sequence[0][j][k];
                }
              }

              else {
                for (k = 0;k < length[i];k++) {
                  seq->real_sequence[i + 1][j][k] = real_sequence[i][j][k] - seq->real_sequence[0][j][k];
                }
              }
            }
          }
        }
        break;
      }

      case true : {
        if (index_parameter) {
          for (i = 0;i < nb_sequence;i++) {
            for (j = 0;j < nb_variable;j++) {
              pindex_param = seq->index_parameter[0];
              pmean = seq->real_sequence[0][j];

              if (type[j] != REAL_VALUE) {
                for (k = 0;k < length[i];k++) {
                  while (*pindex_param < index_parameter[i][k]) {
                    pindex_param++;
                    pmean++;
                  }

                  pindex_param++;
                  if (fabs(int_sequence[i][j][k] - *pmean) <= 180) {
                    seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k] - *pmean++;
                  }
                  else if (int_sequence[i][j][k] - *pmean > 180) {
                    seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k] - *pmean++ - 360;
                  }
                  else {
                    seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k] - *pmean++ + 360;
                  }
                }
              }

              else {
                switch (unit) {

                case DEGREE : {
                  for (k = 0;k < length[i];k++) {
                    while (*pindex_param < index_parameter[i][k]) {
                      pindex_param++;
                      pmean++;
                    }

                    pindex_param++;
                    if (fabs(int_sequence[i][j][k] - *pmean) <= 180) {
                      seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k] - *pmean++;
                    }
                    else if (int_sequence[i][j][k] - *pmean > 180) {
                      seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k] - *pmean++ - 360;
                    }
                    else {
                      seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k] - *pmean++ + 360;
                    }
                  }
                  break;
                }

                case RADIAN : {
                  for (k = 0;k < length[i];k++) {
                    while (*pindex_param < index_parameter[i][k]) {
                      pindex_param++;
                      pmean++;
                    }

                    pindex_param++;
                    if (fabs(real_sequence[i][j][k] - *pmean) <= M_PI) {
                      seq->real_sequence[i + 1][j][k] = real_sequence[i][j][k] - *pmean++;
                    }
                    else if (real_sequence[i][j][k] - *pmean > M_PI) {
                      seq->real_sequence[i + 1][j][k] = real_sequence[i][j][k] - *pmean++ - 2 * M_PI;
                    }
                    else {
                      seq->real_sequence[i + 1][j][k] = real_sequence[i][j][k] - *pmean++ + 2 * M_PI;
                    }
                  }
                  break;
                }
                }
              }
            }
          }
        }

        else {
          for (i = 0;i < nb_sequence;i++) {
            for (j = 0;j < nb_variable;j++) {
              if (type[j] != REAL_VALUE) {
                for (k = 0;k < length[i];k++) {
                  if (fabs(int_sequence[i][j][k] - seq->real_sequence[0][j][k]) <= 180) {
                    seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k] - seq->real_sequence[0][j][k];
                  }
                  else if (int_sequence[i][j][k] - seq->real_sequence[0][j][k] > 180) {
                    seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k] - seq->real_sequence[0][j][k] - 360;
                  }
                  else {
                    seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k] - seq->real_sequence[0][j][k] + 360;
                  }
                }
              }

              else {
                switch (unit) {

                case DEGREE : {
                  for (k = 0;k < length[i];k++) {
                    if (fabs(int_sequence[i][j][k] - seq->real_sequence[0][j][k]) <= 180) {
                      seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k] - seq->real_sequence[0][j][k];
                    }
                    else if (int_sequence[i][j][k] - seq->real_sequence[0][j][k] > 180) {
                      seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k] - seq->real_sequence[0][j][k] - 360;
                    }
                    else {
                      seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k] - seq->real_sequence[0][j][k] + 360;
                    }
                  }
                  break;
                }

                case RADIAN : {
                  for (k = 0;k < length[i];k++) {
                    if (fabs(real_sequence[i][j][k] - seq->real_sequence[0][j][k]) <= M_PI) {
                      seq->real_sequence[i + 1][j][k] = real_sequence[i][j][k] - seq->real_sequence[0][j][k];
                    }
                    else if (int_sequence[i][j][k] - seq->real_sequence[0][j][k] > M_PI) {
                      seq->real_sequence[i + 1][j][k] = real_sequence[i][j][k] - seq->real_sequence[0][j][k] - 2 * M_PI;
                    }
                    else {
                      seq->real_sequence[i + 1][j][k] = real_sequence[i][j][k] - seq->real_sequence[0][j][k] + 2 * M_PI;
                    }
                  }
                  break;
                }
                }
              }
            }
          }
        }
        break;
      }
      }
      break;
    }

    // calcul des residus standardises

    case STANDARDIZED_RESIDUAL : {
      if (index_parameter) {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < nb_variable;j++) {
            pindex_param = seq->index_parameter[0];
            pmean = seq->real_sequence[0][j];
            pstandard_deviation = seq->real_sequence[seq->nb_sequence - 1][j];

            if (type[j] != REAL_VALUE) {
              for (k = 0;k < length[i];k++) {
                while (*pindex_param < index_parameter[i][k]) {
                  pindex_param++;
                  pmean++;
                  pstandard_deviation++;
                }
                pindex_param++;
                if (*pstandard_deviation > 0.) {
                  seq->real_sequence[i + 1][j][k] = (int_sequence[i][j][k] - *pmean) /
                                                    *pstandard_deviation;
                }
                else {
                  seq->real_sequence[i + 1][j][k] = 0.;
                }
                pmean++;
                pstandard_deviation++;
              }
            }

            else {
              for (k = 0;k < length[i];k++) {
                while (*pindex_param < index_parameter[i][k]) {
                  pindex_param++;
                  pmean++;
                  pstandard_deviation++;
                }
                pindex_param++;
                if (*pstandard_deviation > 0.) {
                  seq->real_sequence[i + 1][j][k] = (real_sequence[i][j][k] - *pmean) /
                                                    *pstandard_deviation;
                }
                else {
                  seq->real_sequence[i + 1][j][k] = 0.;
                }
                pmean++;
                pstandard_deviation++;
              }
            }
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < nb_variable;j++) {
            if (type[j] != REAL_VALUE) {
              for (k = 0;k < length[i];k++) {
                if (seq->real_sequence[seq->nb_sequence - 1][j][k] > 0.) {
                   seq->real_sequence[i + 1][j][k] = (int_sequence[i][j][k] - seq->real_sequence[0][j][k]) /
                                                     seq->real_sequence[seq->nb_sequence - 1][j][k];
                }
                else {
                   seq->real_sequence[i + 1][j][k] = 0.;
                }
              }
            }

            else {
              for (k = 0;k < length[i];k++) {
                if (seq->real_sequence[seq->nb_sequence - 1][j][k] > 0.) {
                   seq->real_sequence[i + 1][j][k] = (real_sequence[i][j][k] - seq->real_sequence[0][j][k]) /
                                                     seq->real_sequence[seq->nb_sequence - 1][j][k];
                }
                else {
                   seq->real_sequence[i + 1][j][k] = 0.;
                }
              }
            }
          }
        }
      }
      break;
    }
    }

    if ((output == SEQUENCE) && (!standard_deviation)) {
      for (i = 0;i < seq->nb_variable;i++) {
        seq->min_value[i] = min_value[i];
        seq->max_value[i] = max_value[i];
      }
    }
    else {
      for (i = 0;i < seq->nb_variable;i++) {
        seq->min_value_computation(i);
        seq->max_value_computation(i);
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      seq->build_marginal_histogram(i);
    }

    // ecriture des sequences de moyennes et d'ecart-types

    if (path) {
      switch (format) {
      case 'a' :
        status = seq->pointwise_average_ascii_print(error , path , size ,
                                                    standard_deviation , output);
        break;
      case 's' :
        status = seq->pointwise_average_spreadsheet_print(error , path , size ,
                                                          standard_deviation , output);
        break;
      }

      if (!status) {

#       ifdef MESSAGE
        cout << error;
#       endif

      }
    }

    delete [] size;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des sequences des temps de retour pour une valeur prise
 *  par une variable entiere.
 *
 *  arguments : reference sur un objet StatError, indice de la variable, valeur.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::recurrence_time_sequences(StatError &error , int variable , int value) const

{
  bool status = true;
  register int i , j;
  int inb_sequence , ilength , previous_index , *psequence;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != STATE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    else {
      if (!marginal_distribution[variable]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
        error.update((error_message.str()).c_str());
      }

      else if ((value < marginal_distribution[variable]->offset) ||
               (value >= marginal_distribution[variable]->nb_value) ||
               (marginal_distribution[variable]->frequency[value] == 0)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << STAT_label[STATL_VALUE] << " " << value << " "
                      << STAT_error[STATR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {
    seq = new Sequences(nb_sequence , 1);

    // calcul des sequences des temps de retour

    inb_sequence = 0;

    for (i = 0;i < nb_sequence;i++) {
      previous_index = 0;
      ilength = 0;
      for (j = 0;j < length[i];j++) {
        if (int_sequence[i][variable][j] == value) {
          if (ilength == 0) {
            seq->int_sequence[inb_sequence][0] = new int[length[i]];
            psequence = seq->int_sequence[inb_sequence][0];
          }

          *psequence++ = j - previous_index;
          previous_index = j;
          ilength++;
        }
      }

      if (ilength > 0) {
        seq->length[inb_sequence] = ilength;
        seq->identifier[inb_sequence++] = identifier[i];
      }
    }

    seq->nb_sequence = inb_sequence;

    seq->max_length_computation();
    seq->cumul_length_computation();
    seq->build_length_frequency_distribution();

    seq->min_value_computation(0);
    seq->max_value_computation(0);

    seq->build_marginal_frequency_distribution(0);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des sequences des temps de sejour par une variable entiere.
 *
 *  arguments : reference sur un objet StatError, indice de la variable.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::sojourn_time_sequences(StatError &error , int variable) const

{
  bool status = true;
  register int i , j;
  int ilength , begin_run , *pstate , *psequence , itype[2];
//  int run_length;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != STATE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    else if (!marginal_distribution[variable]) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
      error.update((error_message.str()).c_str());
    }
  }

  if (status) {
    itype[0] = type[variable];
    itype[1] = INT_VALUE;

    seq = new Sequences(nb_sequence , identifier , length , NULL ,
                        IMPLICIT_TYPE , 2 , itype);

    // calcul des sequences de temps de sejour

    if ((index_parameter_type == TIME) && (index_interval->variance > 0.)) {  // pour les suivis de croissance manguier
      for (i = 0;i < nb_sequence;i++) {                                       // et les pousses d'arabidopsis
        pstate = seq->int_sequence[i][0];
        psequence = seq->int_sequence[i][1];
        begin_run = index_parameter_distribution->offset;
//        begin_run = 0;
        ilength = 0;

        for (j = 0;j < length[i] - 1;j++) {
          if (int_sequence[i][variable][j + 1] != int_sequence[i][variable][j]) {
            *pstate++ = int_sequence[i][variable][j];
            *psequence++ = index_parameter[i][j + 1] - begin_run;
            begin_run = index_parameter[i][j + 1];
            ilength++;
          }
        }

        *pstate = int_sequence[i][variable][length[i] - 1];
        *psequence = index_parameter[i][j] + 1 - begin_run;

        seq->length[i] = ilength + 1;
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        pstate = seq->int_sequence[i][0];
        psequence = seq->int_sequence[i][1];
//        run_length = 1;
        begin_run = 0;
        ilength = 0;

        for (j = 0;j < length[i] - 1;j++) {
          if (int_sequence[i][variable][j + 1] != int_sequence[i][variable][j]) {
            *pstate++ = int_sequence[i][variable][j];
//            *psequence++ = run_length;
//            run_length = 0;
            *psequence++ = j + 1 - begin_run;
            begin_run = j + 1;
            ilength++;
          }

//          run_length++;
        }

        *pstate = int_sequence[i][variable][length[i] - 1];
//        *psequence = run_length;
        *psequence = length[i] - begin_run;

        seq->length[i] = ilength + 1;
      }
    }

    seq->max_length_computation();
    seq->cumul_length_computation();
    delete seq->length_distribution;
    seq->build_length_frequency_distribution();

    for (i = 0;i < 2;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);

      seq->build_marginal_frequency_distribution(i);
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Discretisation des positions.
 *
 *  arguments : reference sur un objet StatError,
 *              pas de discretisation.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::transform_position(StatError &error , int step) const

{
  bool status = true;
  register int i , j , k , m;
  int inter_position , nb_unit , *ilength , **pisequence;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (index_parameter_type != POSITION) {
    status = false;
    error.correction_update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE] , SEQ_index_parameter_word[POSITION]);
  }

  else if ((step < 1) || ((index_interval) && (step > index_interval->mean))) {
    status = false;
    error.update(SEQ_error[SEQR_POSITION_STEP]);
  }

  for (i = 0;i < nb_variable;i++) {
    if ((type[i] != INT_VALUE) && (type[i] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (status) {
    ilength = new int[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      ilength[i] = index_parameter[i][length[i]] / step + 1 + length[i];
    }

    seq = new Sequences(nb_sequence , identifier , ilength , nb_variable);
    delete [] ilength;

    // extraction des sequences

    pisequence = new int*[nb_variable];

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < nb_variable;j++) {
        pisequence[j] = seq->int_sequence[i][j];
      }
      seq->length[i] = 0;

      for (j = 0;j <= length[i];j++) {
        if (j == 0) {
          inter_position = index_parameter[i][j];
        }
        else {
          inter_position = index_parameter[i][j] - index_parameter[i][j - 1];
        }

        nb_unit = (inter_position % step == 0 ? inter_position / step : inter_position / step + 1);
        if ((nb_unit > 0) && (j < length[i])) {
          nb_unit--;
        }

        for (k = 0;k < nb_variable;k++) {
          for (m = 0;m < nb_unit;m++) {
            *pisequence[k]++ = (int)min_value[k] - 1;
          }
          if (j < length[i]) {
            *pisequence[k]++ = int_sequence[i][k][j];
          }
        }

        if (j < length[i]) {
          seq->length[i] += nb_unit + 1;
        }
        else {
          seq->length[i] += nb_unit;
        }
      }
    }

    seq->max_length_computation();
    seq->cumul_length_computation();
    delete seq->length_distribution;
    seq->build_length_frequency_distribution();

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value[i] = min_value[i] - 1;
      seq->max_value[i] = max_value[i];

      seq->build_marginal_frequency_distribution(i);
    }

    delete [] pisequence;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Croisement des sequences.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::cross(StatError &error) const

{
  bool status = true;
  register int i , j , k , m;
  int sense = 0 , *ilength;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (type[0] != INT_VALUE) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_VARIABLE] << " " << 1 << ": "
                  << STAT_error[STATR_VARIABLE_TYPE];
    error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
  }

  for (i = 1;i < nb_sequence;i++) {
    if (length[i] > length[i - 1]) {
      if (sense == 0) {
        sense++;
      }
      else if (sense == -1) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEQUENCE] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_LENGTH];
        error.update((error_message.str()).c_str());
      }
    }

    else if (length[i] < length[i - 1]) {
      if (sense == 0) {
        sense--;
      }
      else if (sense == 1) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEQUENCE] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_LENGTH];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {
    ilength = new int[max_length];
    for (i = 0;i < max_length;i++) {
      ilength[i] = nb_sequence;
    }

    seq = new Sequences(max_length , NULL , ilength , NULL , IMPLICIT_TYPE ,
                        nb_variable , type);
    delete [] ilength;

    // constitution des sequences croisees

    for (i = 0;i < seq->nb_sequence;i++) {
      j = 0;
      while (length[j] <= i) {
        j++;
      }

      k = 0;
      do {
        for (m = 0;m < seq->nb_variable;m++) {
          if (seq->type[m] != REAL_VALUE) {
            seq->int_sequence[i][m][k] = int_sequence[j][m][i];
          }
          else {
            seq->real_sequence[i][m][k] = real_sequence[j][m][i];
          }
        }
        j++;
        k++;
      }
      while ((j < nb_sequence) && (length[j] > i));

      seq->length[i] = k;
    }

    seq->max_length = nb_sequence;
    seq->cumul_length = cumul_length;
    delete seq->length_distribution;
    seq->build_length_frequency_distribution();

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value[i] = min_value[i];
      seq->max_value[i] = max_value[i];

      if (marginal_distribution[i]) {
        seq->marginal_distribution[i] = new FrequencyDistribution(*marginal_distribution[i]);
      }
      if (marginal_histogram[i]) {
        seq->marginal_histogram[i] = new Histogram(*marginal_histogram[i]);
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la longueur maximum des sequences.
 *
 *--------------------------------------------------------------*/

void Sequences::max_length_computation()

{
  register int i;


  max_length = length[0];
  for (i = 1;i < nb_sequence;i++) {
    if (length[i] > max_length) {
      max_length = length[i];
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la longueur cumulee des sequences.
 *
 *--------------------------------------------------------------*/

void Sequences::cumul_length_computation()

{
  register int i;


  cumul_length = 0;
  for (i = 0;i < nb_sequence;i++) {
    cumul_length += length[i];
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction de la loi empirique des longueurs des sequences.
 *
 *--------------------------------------------------------------*/

void Sequences::build_length_frequency_distribution()

{
  register int i;


  length_distribution = new FrequencyDistribution(max_length + 1);

  length_distribution->nb_element = nb_sequence;
  for (i = 0;i < nb_sequence;i++) {
    (length_distribution->frequency[length[i]])++;
  }

  length_distribution->nb_value_computation();
  length_distribution->offset_computation();
  length_distribution->max_computation();
  length_distribution->mean_computation();
  length_distribution->variance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Calcul des temps a partir des intervalles de temps /
 *  calcul des positions a partir des intervalles inter-positions.
 *
 *--------------------------------------------------------------*/

void Sequences::index_parameter_computation()

{
  if ((index_parameter_type == TIME_INTERVAL) || (index_parameter_type == POSITION_INTERVAL)) {
    register int i , j;


    switch (index_parameter_type) {
    case TIME_INTERVAL :
      index_parameter_type = TIME;
      break;
    case POSITION_INTERVAL :
      index_parameter_type = POSITION;
      break;
    }

    for (i = 0;i < nb_sequence;i++) {
      for (j = 1;j < (index_parameter_type == POSITION ? length[i] + 1 : length[i]);j++) {
        index_parameter[i][j] += index_parameter[i][j - 1];
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur minimum prise par le parametre d'index.
 *
 *--------------------------------------------------------------*/

int Sequences::min_index_parameter_computation() const

{
  register int i;
  int min_index_parameter = I_DEFAULT;


  if (index_parameter) {
    min_index_parameter = index_parameter[0][0];
    for (i = 1;i < nb_sequence;i++) {
      if (index_parameter[i][0] < min_index_parameter) {
        min_index_parameter = index_parameter[i][0];
      }
    }
  }

  return min_index_parameter;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur maximum prise par le parametre d'index.
 *
 *  argument : flag derniere position.
 *
 *--------------------------------------------------------------*/

int Sequences::max_index_parameter_computation(bool last_position) const

{
  register int i;
  int max_index_parameter = I_DEFAULT;


  if (index_parameter) {
    if ((index_parameter_type == TIME) || (last_position)) {
      max_index_parameter = index_parameter[0][length[0] - 1];
      for (i = 1;i < nb_sequence;i++) {
        if (index_parameter[i][length[i] - 1] > max_index_parameter) {
          max_index_parameter = index_parameter[i][length[i] - 1];
        }
      }
    }

    else {
      max_index_parameter = index_parameter[0][length[0]];
      for (i = 1;i < nb_sequence;i++) {
        if (index_parameter[i][length[i]] > max_index_parameter) {
          max_index_parameter = index_parameter[i][length[i]];
        }
      }
    }
  }

  return max_index_parameter;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la loi empirique des parametres d'index.
 *
 *--------------------------------------------------------------*/

void Sequences::build_index_parameter_frequency_distribution()

{
  if (index_parameter) {
    register int i , j;


    index_parameter_distribution = new FrequencyDistribution(max_index_parameter_computation() + 1);

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < (index_parameter_type == POSITION ? length[i] + 1 : length[i]);j++) {
        (index_parameter_distribution->frequency[index_parameter[i][j]])++;
      }
    }

    index_parameter_distribution->offset_computation();
    index_parameter_distribution->nb_element = cumul_length;
    if (index_parameter_type == POSITION) {
      index_parameter_distribution->nb_element += nb_sequence;
    }
    index_parameter_distribution->max_computation();
    index_parameter_distribution->mean_computation();
    index_parameter_distribution->variance_computation();
  }
}


/*--------------------------------------------------------------*
 *
 *  Extraction de la loi empirique des intervalles entre index successifs.
 *
 *--------------------------------------------------------------*/

void Sequences::index_interval_computation()

{
//  if ((index_parameter_type == TIME) || ((index_parameter_type == POSITION) &&
//       (type[0] != NB_INTERNODE))) {
  if ((index_parameter_type == TIME) || (index_parameter_type == POSITION)) {
    register int i , j;


    index_interval = new FrequencyDistribution(max_index_parameter_computation(true) + 1);

    // constitution de la loi empirique des intervalles entre index successifs

    for (i = 0;i < nb_sequence;i++) {
      for (j = 1;j < length[i];j++) {
        (index_interval->frequency[index_parameter[i][j] - index_parameter[i][j - 1]])++;
      }
    }

    index_interval->nb_value_computation();
    index_interval->offset_computation();
    index_interval->nb_element = cumul_length - nb_sequence;
    index_interval->max_computation();
    index_interval->mean_computation();
    index_interval->variance_computation();
  }
}


/*--------------------------------------------------------------*
 *
 *  Extraction de la loi empirique des intervalles entre index successifs
 *  pour une valeur d'une variable entiere.
 *
 *  arguments : reference sur un objet StatError, indice de la variable, valeur.
 *
 *--------------------------------------------------------------*/

FrequencyDistribution* Sequences::value_index_interval_computation(StatError &error , int variable ,
                                                                   int value) const

{
  bool status = true;
  register int i , j;
  int previous_index_param , *pindex_param , *pisequence;
  FrequencyDistribution *value_index_interval;


  value_index_interval = NULL;
  error.init();

  if ((index_parameter_type != TIME) && (index_parameter_type != POSITION)) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (!marginal_distribution[variable]) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
      error.update((error_message.str()).c_str());
    }

    else if ((value < marginal_distribution[variable]->offset) ||
             (value >= marginal_distribution[variable]->nb_value) ||
             (marginal_distribution[variable]->frequency[value] <= 1)) {
      status = false;
      error.update(SEQ_error[SEQR_VALUE]);
    }
  }

  if (status) {
    value_index_interval = new FrequencyDistribution(max_index_parameter_computation(true) + 1);

    for (i = 0;i < nb_sequence;i++) {
      pindex_param = index_parameter[i];
      pisequence = int_sequence[i][variable];
      previous_index_param = I_DEFAULT;

      for (j = 0;j < length[i];j++) {
        if (*pisequence == value) {
          if (previous_index_param != I_DEFAULT) {
            (value_index_interval->frequency[*pindex_param - previous_index_param])++;
          }
          previous_index_param = *pindex_param;
        }

        pindex_param++;
        pisequence++;
      }
    }

    // extraction des caracteristiques de la loi empirique

    value_index_interval->nb_value_computation();
    value_index_interval->offset_computation();
    value_index_interval->nb_element_computation();

    if (value_index_interval->nb_element == 0) {
      delete value_index_interval;
      value_index_interval = NULL;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    else {
      value_index_interval->max_computation();
      value_index_interval->mean_computation();
      value_index_interval->variance_computation();
    }
  }

  return value_index_interval;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur minimum prise par une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Sequences::min_value_computation(int variable)

{
  register int i , j;


  if ((type[variable] != REAL_VALUE) && (type[variable] != AUXILIARY)) {
    min_value[variable] = int_sequence[0][variable][0];

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        if (int_sequence[i][variable][j] < min_value[variable]) {
          min_value[variable] = int_sequence[i][variable][j];
        }
      }
    }
  }

  else {
    min_value[variable] = real_sequence[0][variable][0];

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        if (real_sequence[i][variable][j] < min_value[variable]) {
          min_value[variable] = real_sequence[i][variable][j];
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur maximum prise par une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Sequences::max_value_computation(int variable)

{
  register int i , j;


  if ((type[variable] != REAL_VALUE) && (type[variable] != AUXILIARY)) {
    max_value[variable] = int_sequence[0][variable][0];

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        if (int_sequence[i][variable][j] > max_value[variable]) {
          max_value[variable] = int_sequence[i][variable][j];
        }
      }
    }
  }

  else {
    max_value[variable] = real_sequence[0][variable][0];

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        if (real_sequence[i][variable][j] > max_value[variable]) {
          max_value[variable] = real_sequence[i][variable][j];
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la loi marginale empirique pour une variable entiere positive.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Sequences::marginal_frequency_distribution_computation(int variable)

{
  register int i , j;


  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    marginal_distribution[variable]->frequency[i] = 0;
  }

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < length[i];j++) {
      (marginal_distribution[variable]->frequency[int_sequence[i][variable][j]])++;
    }
  }

  marginal_distribution[variable]->offset = (int)min_value[variable];
  marginal_distribution[variable]->nb_element_computation();
//  marginal_distribution[variable]->nb_element = cumul_length;
  marginal_distribution[variable]->max_computation();
  marginal_distribution[variable]->mean_computation();
  marginal_distribution[variable]->variance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Construction de la loi marginale empirique pour une variable entiere positive.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Sequences::build_marginal_frequency_distribution(int variable)

{
  if (type[variable] != AUXILIARY) {
    if ((type[variable] != REAL_VALUE) && (min_value[variable] >= 0) &&
        (max_value[variable] <= MARGINAL_DISTRIBUTION_MAX_VALUE)) {
      marginal_distribution[variable] = new FrequencyDistribution((int)max_value[variable] + 1);
      marginal_frequency_distribution_computation(variable);
    }

    else {
      build_marginal_histogram(variable);
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction de l'histogramme marginal pour une variable.
 *
 *  arguments : indice de la variable, pas de regroupement, valeur minimum.
 *
 *--------------------------------------------------------------*/

void Sequences::build_marginal_histogram(int variable , double step , double imin_value)

{
  if ((!marginal_histogram[variable]) || (step != marginal_histogram[variable]->step) ||
      (imin_value != D_INF)) {
    register int i , j;
    int *pisequence;
    double *prsequence;


    // construction de l'histogramme

    if (step == D_DEFAULT) {
      step = MAX(::round((max_value[variable] - min_value[variable]) * HISTOGRAM_FREQUENCY / cumul_length) , 1);

#     ifdef MESSAGE
      cout << "\n" << STAT_label[STATL_VARIABLE] << " " << variable + 1 << " - "
           << STAT_label[STATL_STEP] << ": " << step << endl;
//           << " (" << min_value[variable] << ", " << max_value[variable] << ")"
#     endif

    }

    if (imin_value == D_INF) {
      imin_value = floor(min_value[variable] / step) * step;
    }

    if (marginal_histogram[variable]) {
      marginal_histogram[variable]->nb_category = (int)floor((max_value[variable] - imin_value) / step) + 1;

      delete [] marginal_histogram[variable]->frequency;
      marginal_histogram[variable]->frequency = new int[marginal_histogram[variable]->nb_category];
    }

    else {
      marginal_histogram[variable] = new Histogram((int)floor((max_value[variable] - imin_value) / step) + 1 , false);

      marginal_histogram[variable]->nb_element = cumul_length;
      marginal_histogram[variable]->type = type[variable];
    }

    marginal_histogram[variable]->step = step;
    marginal_histogram[variable]->min_value = imin_value;
    marginal_histogram[variable]->max_value = ceil(max_value[variable] / step) * step;

    // calcul des frequences

    for (i = 0;i < marginal_histogram[variable]->nb_category;i++) {
      marginal_histogram[variable]->frequency[i] = 0;
    }

    if ((type[variable] != REAL_VALUE) && (type[variable] != AUXILIARY)) {
      for (i = 0;i < nb_sequence;i++) {
        pisequence = int_sequence[i][variable];
        for (j = 0;j < length[i];j++) {
//          (marginal_histogram[variable]->frequency[(int)((*pisequence++ - imin_value) / step)])++;
          (marginal_histogram[variable]->frequency[(int)floor((*pisequence++ - imin_value) / step)])++;
        }
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        prsequence = real_sequence[i][variable];
        for (j = 0;j < length[i];j++) {
//          (marginal_histogram[variable]->frequency[(int)((*prsequence++ - imin_value) / step)])++;
          (marginal_histogram[variable]->frequency[(int)floor((*prsequence++ - imin_value) / step)])++;
        }
      }
    }

    marginal_histogram[variable]->max_computation();
  }
}


/*--------------------------------------------------------------*
 *
 *  Changement du pas de regroupement de l'histogramme marginal.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              pas de regroupement, valeur minimum.
 *
 *--------------------------------------------------------------*/

bool Sequences::select_step(StatError &error , int variable ,
                            double step , double imin_value)

{
  bool status = true;


  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (!marginal_histogram[variable]) {
      status = false;
      error.update(STAT_error[STATR_MARGINAL_HISTOGRAM]);
    }
    if ((step <= 0.) || ((type[variable] != REAL_VALUE) && ((int)step != step))) {
      status = false;
      error.update(STAT_error[STATR_HISTOGRAM_STEP]);
    }
    if ((imin_value != D_INF) && ((imin_value <= min_value[variable] - step) ||
         (imin_value > min_value[variable]) || ((type[variable] != REAL_VALUE) &&
          ((int)imin_value != imin_value)))) {
      status = false;
      error.update(STAT_error[STATR_HISTOGRAM_MIN_VALUE]);
    }
  }

  if (status) {
    build_marginal_histogram(variable , step , imin_value);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la moyenne d'une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

double Sequences::mean_computation(int variable) const

{
  register int i , j;
  double mean;


  if (marginal_distribution[variable]) {
    mean = marginal_distribution[variable]->mean;
  }

  else {
    mean = 0.;

    if ((type[variable] != REAL_VALUE) && (type[variable] != AUXILIARY)) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          mean += int_sequence[i][variable][j];
        }
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          mean += real_sequence[i][variable][j];
        }
      }
    }

    mean /= cumul_length;
  }

  return mean;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la variance d'une variable.
 *
 *  arguments : indice de la variable, moyenne.
 *
 *--------------------------------------------------------------*/

double Sequences::variance_computation(int variable , double mean) const

{
  register int i , j;
  double variance , diff;


  if (marginal_distribution[variable]) {
    variance = marginal_distribution[variable]->variance;
  }

  else {
    variance = 0.;

    if (cumul_length > 1) {
      if ((type[variable] != REAL_VALUE) && (type[variable] != AUXILIARY)) {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            diff = int_sequence[i][variable][j] - mean;
            variance += diff * diff;
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            diff = real_sequence[i][variable][j] - mean;
            variance += diff * diff;
          }
        }
      }

      variance /= (cumul_length - 1);
    }
  }

  return variance;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'ecart absolu moyen d'une variable.
 *
 *  arguments : indice de la variable, moyenne.
 *
 *--------------------------------------------------------------*/

double Sequences::mean_absolute_deviation_computation(int variable , double mean) const

{
  register int i , j;
  double mean_absolute_deviation ;


  if (marginal_distribution[variable]) {
    mean_absolute_deviation = marginal_distribution[variable]->mean_absolute_deviation_computation();
  }

  else {
    mean_absolute_deviation = 0.;

    if ((type[variable] != REAL_VALUE) && (type[variable] != AUXILIARY)) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          mean_absolute_deviation += fabs(int_sequence[i][variable][j] - mean);
        }
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          mean_absolute_deviation += fabs(real_sequence[i][variable][j] - mean);
        }
      }
    }

    mean_absolute_deviation /= cumul_length;
  }

  return mean_absolute_deviation;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la difference absolue moyenne d'une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

double Sequences::mean_absolute_difference_computation(int variable) const

{
  register int i , j , k , m;
  double mean_absolute_difference;


  mean_absolute_difference = 0.;

  if (cumul_length > 1) {
    if ((type[variable] != REAL_VALUE) && (type[variable] != AUXILIARY)) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          for (k = j + 1;k < length[i];k++) {
            mean_absolute_difference += abs(int_sequence[i][variable][j] -
                                            int_sequence[i][variable][k]);
          }
        }

        for (j = i + 1;j < nb_sequence;j++) {
          for (k = 0;k < length[i];k++) {
            for (m = 0;m < length[j];m++) {
              mean_absolute_difference += abs(int_sequence[i][variable][k] -
                                              int_sequence[j][variable][m]);
            }
          }
        }
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          for (k = j + 1;k < length[i];k++) {
            mean_absolute_difference += fabs(real_sequence[i][variable][j] -
                                             real_sequence[i][variable][k]);
          }
        }

        for (j = i + 1;j < nb_sequence;j++) {
          for (k = 0;k < length[i];k++) {
            for (m = 0;m < length[j];m++) {
              mean_absolute_difference += fabs(real_sequence[i][variable][k] -
                                               real_sequence[j][variable][m]);
            }
          }
        }
      }
    }

    mean_absolute_difference = 2 * mean_absolute_difference /
                               (cumul_length * (double)(cumul_length - 1));
  }

  return mean_absolute_difference;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du coefficient d'asymetrie d'une variable.
 *
 *  arguments : indice de la variable, moyenne, variance.
 *
 *--------------------------------------------------------------*/

double Sequences::skewness_computation(int variable , double mean , double variance) const

{
  register int i , j;
  double skewness , diff;


  if (marginal_distribution[variable]) {
    skewness = marginal_distribution[variable]->skewness_computation();
  }

  else {
    skewness = 0.;

    if ((cumul_length > 2) && (variance > 0.)) {
      if ((type[variable] != REAL_VALUE) && (type[variable] != AUXILIARY)) {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            diff = int_sequence[i][variable][j] - mean;
            skewness += diff * diff * diff;
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            diff = real_sequence[i][variable][j] - mean;
            skewness += diff * diff * diff;
          }
        }
      }

      skewness = skewness * cumul_length / ((cumul_length - 1) *
                  (double)(cumul_length - 2) * pow(variance , 1.5));
    }
  }

  return skewness;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'exces d'applatissement d'une variable :
 *  exces d'applatissement = coefficient d'applatissement - 3.
 *
 *  arguments : indice de la variable, moyenne, variance.
 *
 *--------------------------------------------------------------*/

double Sequences::kurtosis_computation(int variable , double mean , double variance) const

{
  register int i , j;
  double kurtosis , diff;


  if (marginal_distribution[variable]) {
    kurtosis = marginal_distribution[variable]->kurtosis_computation();
  }

  else {
    if (variance == 0.) {
      kurtosis = -2.;
    }

    else {
      kurtosis = 0.;

      if ((type[variable] != REAL_VALUE) && (type[variable] != AUXILIARY)) {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            diff = int_sequence[i][variable][j] - mean;
            kurtosis += diff * diff * diff * diff;
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            diff = real_sequence[i][variable][j] - mean;
            kurtosis += diff * diff * diff * diff;
          }
        }
      }

      kurtosis = kurtosis / ((cumul_length - 1) * variance * variance) - 3.;
    }
  }

  return kurtosis;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la direction moyenne d'une variable circulaire.
 *
 *  arguments : indice de la variable, unite (DEGREE/RADIAN).
 *
 *--------------------------------------------------------------*/

double* Sequences::mean_direction_computation(int variable , int unit) const

{
  register int i , j;
  double *mean_direction;


  mean_direction = new double[4];
//  mean_direction = new double[2];

  mean_direction[0] = 0.;
  mean_direction[1] = 0.;

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        mean_direction[0] += cos(int_sequence[i][variable][j] * M_PI / 180);
        mean_direction[1] += sin(int_sequence[i][variable][j] * M_PI / 180);
      }
    }
    break;
  }

  case REAL_VALUE : {
    switch (unit) {

    case DEGREE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          mean_direction[0] += cos(real_sequence[i][variable][j] * M_PI / 180);
          mean_direction[1] += sin(real_sequence[i][variable][j] * M_PI / 180);
        }
      }
      break;
    }

    case RADIAN : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          mean_direction[0] += cos(real_sequence[i][variable][j]);
          mean_direction[1] += sin(real_sequence[i][variable][j]);
        }
      }
      break;
    }
    }
    break;
  }
  }

  mean_direction[0] /= cumul_length;
  mean_direction[1] /= cumul_length;

  mean_direction[2] = sqrt(mean_direction[0] * mean_direction[0] +
                           mean_direction[1] * mean_direction[1]);

  if (mean_direction[2] > 0.) {
    mean_direction[3] = atan(mean_direction[1] / mean_direction[0]);

    if (mean_direction[0] < 0.) {
      mean_direction[3] += M_PI;
    }
    if (unit == DEGREE) {
      mean_direction[3] *= (180 / M_PI);
    }
  }

  else {
    mean_direction[3] = D_DEFAULT;
  }

  return mean_direction;
}


};  // namespace sequence_analysis
