/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       StructureAnalysis: Exploring and Analyzing Plant Architecture
 *
 *       Copyright 1995-2018 CIRAD AGAP
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
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



#include <limits.h>
#include <math.h>

#include <string>
#include <vector>
#include <sstream>
#include <iomanip>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>

#include "stat_tool/stat_label.h"

#include "stat_tool/quantile_computation.hpp"

#include "sequences.h"
#include "sequence_label.h"

using namespace std;
using namespace boost::math;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Merging of Sequences objects.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] nb_sample number of Sequences objects,
 *  \param[in] iseq      pointer on the Sequences objects.
 *
 *  \return              Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::merge(StatError &error , int nb_sample , const Sequences **iseq) const

{
  bool status = true;
  int i , j , k , m , n , p , q;
  int inb_sequence , cumul_nb_sequence , *ilength , *iidentifier , **ivertex_identifier;
  const FrequencyDistribution **phisto;
  Sequences *seq;
  const Sequences **pseq;


  seq = NULL;
  error.init();

  for (i = 0;i < nb_sample;i++) {
    if (iseq[i]->index_param_type != index_param_type) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_INDEX_PARAMETER_TYPE];

      if (index_param_type == IMPLICIT_TYPE) {
        error.update((error_message.str()).c_str());
      }
      else {
        error.correction_update((error_message.str()).c_str() , SEQ_index_parameter_word[index_param_type]);
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

    // computation of the number of sequences

    inb_sequence = 0;
    for (i = 0;i < nb_sample;i++) {
      inb_sequence += pseq[i]->nb_sequence;
    }

    // comparison of the sequence identifiers

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

    // copy of sequence lengths

    ilength = new int[inb_sequence];

    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < pseq[j]->nb_sequence;k++) {
        ilength[i++] = pseq[j]->length[k];
      }
    }

    // comparison of vertex identifiers

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
                        index_param_type , nb_variable , type);
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
          for (m = 0;m < (pseq[j]->index_param_type == POSITION ? pseq[j]->length[k] + 1 : pseq[j]->length[k]);m++) {
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

//    if ((seq->index_param_type == TIME) || ((seq->index_param_type == POSITION) &&
//         (seq->type[0] != NB_INTERNODE))) {
    if ((seq->index_param_type == TIME) || (seq->index_param_type == POSITION)) {
      for (i = 0;i < nb_sample;i++) {
        phisto[i] = pseq[i]->index_interval;
      }
      seq->index_interval = new FrequencyDistribution(nb_sample , phisto);
    }

    // copy of values

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


/*--------------------------------------------------------------*/
/**
 *  \brief Merging of Sequences objects.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] nb_sample number of Sequences objects,
 *  \param[in] iseq      pointer on the Sequences objects.
 *
 *  \return              Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::merge(StatError &error , int nb_sample , const vector<Sequences> iseq) const

{
  int i;
  Sequences *seq;
  const Sequences **pseq;


  pseq = new const Sequences*[nb_sample];
  for (i = 0;i < nb_sample;i++) {
    pseq[i] = new Sequences(iseq[i]);
  }

  seq = merge(error , nb_sample , pseq);

  for (i = 0;i < nb_sample;i++) {
    delete pseq[i];
  }
  delete [] pseq;

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Shifting of values of a variable.
 *
 *  \param[in] error       reference on a StatError object,
 *  \param[in] variable    variable index,
 *  \param[in] shift_param integer shifting parameter.
 *
 *  \return                Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::shift(StatError &error , int variable , int shift_param) const

{
  bool status = true;
  int i , j;
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

    // shifting of integer values

    case INT_VALUE : {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->length[i];j++) {
          seq->int_sequence[i][variable][j] = int_sequence[i][variable][j] + shift_param;
        }
      }
      break;
    }

    // shifting of real values

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
        seq->marginal_distribution[variable] = new FrequencyDistribution(*marginal_distribution[variable] ,
                                                                         SHIFT , shift_param);
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


/*--------------------------------------------------------------*/
/**
 *  \brief Shifting of values of a real-valued variable.
 *
 *  \param[in] error       reference on a StatError object,
 *  \param[in] variable    variable index,
 *  \param[in] shift_param real shifting parameter.
 *
 *  \return                Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::shift(StatError &error , int variable , double shift_param) const

{
  bool status = true;
  int i , j;
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

    // shifting of real values

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


/*--------------------------------------------------------------*/
/**
 *  \brief Thresholding of values of a variable.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] variable  variable index,
 *  \param[in] threshold integer threshold,
 *  \param[in] mode      mode (ABOVE/BELOW).
 *
 *  \return              Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::thresholding(StatError &error , int variable , int threshold ,
                                   threshold_direction mode) const

{
  bool status = true;
  int i , j;
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

    // thresholding of integer values

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

    // thresholding of real values

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


/*--------------------------------------------------------------*/
/**
 *  \brief Thresholding of values of a real-valued variable.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] variable  variable index,
 *  \param[in] threshold real threshold,
 *  \param[in] mode      mode (ABOVE/BELOW).
 *
 *  \return              Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::thresholding(StatError &error , int variable , double threshold ,
                                   threshold_direction mode) const

{
  bool status = true;
  int i , j;
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

    // thresholding of real values

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


/*--------------------------------------------------------------*/
/**
 *  \brief Clustering of values of a variable.
 *
 *  \param[in] seq      reference on a Sequences object,
 *  \param[in] variable variable index,
 *  \param[in] step     clustering step,
 *  \param[in] mode     mode (FLOOR/ROUND/CEIL).
 */
/*--------------------------------------------------------------*/

void Sequences::cluster(const Sequences &seq , int variable , int step , rounding mode)

{
  int i , j;


  switch (type[variable]) {

  // clustering of integer values

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
      marginal_distribution[variable] = new FrequencyDistribution(*(seq.marginal_distribution[variable]) ,
                                                                  CLUSTER , step , mode);
    }
    else {
      build_marginal_frequency_distribution(variable);
    }
    break;
  }

  // clustering of real values

  case REAL_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        real_sequence[i][variable][j] = seq.real_sequence[i][variable][j] / step;
      }
    }

    min_value[variable] = seq.min_value[variable] / step;
    max_value[variable] = seq.max_value[variable] / step;

    build_marginal_histogram(variable , seq.marginal_histogram[variable]->bin_width / step);

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


/*--------------------------------------------------------------*/
/**
 *  \brief Clustering of values of a variable.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] step     clustering step,
 *  \param[in] mode     mode (FLOOR/ROUND/CEIL).
 *
 *  \return             Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::cluster(StatError &error , int variable , int step , rounding mode) const

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


/*--------------------------------------------------------------*/
/**
 *  \brief Transcoding of categories of an integer-valued variable.
 *
 *  \param[in] seq          reference on a Sequences object,
 *  \param[in] ivariable    variable index,
 *  \param[in] min_category lowest category,
 *  \param[in] max_category highest category,
 *  \param[in] category     transcoding table,
 *  \param[in] add_variable flag for adding a variable.
 */
/*--------------------------------------------------------------*/

void Sequences::transcode(const Sequences &seq , int ivariable , int min_category ,
                          int max_category , int *category , bool add_variable)

{
  int i , j , k;
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
      for (j = 0;j < (index_param_type == POSITION ? length[i] + 1 : length[i]);j++) {
        index_parameter[i][j] = seq.index_parameter[i][j];
      }
    }
  }

  if (add_variable) {
    variable = 0;
    offset = 1;
  }
  else {
    variable = ivariable;
    offset = 0;
  }

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < nb_variable;j++) {
      if ((type[j] != REAL_VALUE) && (type[j] != AUXILIARY)) {

        // transcoding of categories

        if (j == variable) {
          for (k = 0;k < length[i];k++) {
            int_sequence[i][j][k] = category[seq.int_sequence[i][ivariable][k] -
                                             (int)seq.min_value[variable]] + min_category;
          }
        }

        // copy of integer values

        else {
          for (k = 0;k < length[i];k++) {
            int_sequence[i][j][k] = seq.int_sequence[i][j - offset][k];
          }
        }
      }

      // copy of real values

      else {
        for (k = 0;k < length[i];k++) {
          real_sequence[i][j][k] = seq.real_sequence[i][j - offset][k];
        }
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    if (i == variable) {
      min_value[i] = min_category;
      max_value[i] = max_category;

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


/*--------------------------------------------------------------*/
/**
 *  \brief Transcoding of categories of an integer-valued variable.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] category transcoding table.
 *
 *  \return             Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::transcode(StatError &error , int variable , int *category) const

{
  bool status = true , *presence;
  int i;
  int min_category , max_category;
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
      min_category = category[0];
      max_category = category[0];

      for (i = 1;i <= (int)(max_value[variable] - min_value[variable]);i++) {
        if (category[i] < min_category) {
          min_category = category[i];
        }
        if (category[i] > max_category) {
          max_category = category[i];
        }
      }

      if (max_category - min_category == 0) {
        status = false;
        error.update(STAT_error[STATR_NB_CATEGORY]);
      }

      if (max_category - min_category > (int)(max_value[variable] - min_value[variable])) {
        status = false;
        error.update(STAT_error[STATR_NON_CONSECUTIVE_CATEGORIES]);
      }
    }

    if (status) {
      presence = new bool[max_category - min_category + 1];
      for (i = 0;i <= max_category - min_category;i++) {
        presence[i] = false;
      }

      for (i = 0;i <= (int)(max_value[variable] - min_value[variable]);i++) {
        presence[category[i] - min_category] = true;
      }

      for (i = 0;i <= max_category - min_category;i++) {
        if (!presence[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_error[STATR_MISSING_CATEGORY] << " " << i + min_category;
          error.update((error_message.str()).c_str());
        }
      }

      delete [] presence;
    }

    if (status) {
      for (i = 0;i <= (int)(max_value[variable] - min_value[variable]);i++) {
        category[i] -= min_category;
      }

      seq = new Sequences(nb_sequence , identifier , length , vertex_identifier ,
                          index_param_type , nb_variable , type);
      seq->transcode(*this , variable , min_category , max_category , category);
    }
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Transcoding of categories of an integer-valued variable.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] category transcoding table.
 *
 *  \return             Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::transcode(StatError &error , int variable , vector<int> category) const

{
  return transcode(error , variable , category.data());
}


/*--------------------------------------------------------------*/
/**
 *  \brief Partitioning of values of a variable.
 *
 *  \param[in] error     reference on a StatError object, 
 *  \param[in] variable  variable index,
 *  \param[in] nb_class  number of classes,
 *  \param[in] ilimit    integer limits between classes (beginning of classes).
 *
 *  \return              Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::cluster(StatError &error , int variable ,
                              int nb_class , int *ilimit) const

{
  bool status = true;
  int i , j , k;
  int *int_limit , *category;
  variable_nature *itype;
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
        category = new int[(int)(max_value[variable] - min_value[variable]) + 1];

        i = 0;
        for (j = 0;j < nb_class;j++) {
          for (k = int_limit[j];k < int_limit[j + 1];k++) {
            category[i++] = j;
          }
        }

        seq = new Sequences(nb_sequence , identifier , length , vertex_identifier ,
                            index_param_type , nb_variable , type);
        seq->transcode(*this , variable , 0 , nb_class - 1 , category);

        delete [] category;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Partitioning of values of a variable.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] nb_class number of classes,
 *  \param[in] ilimit   integer limits between classes (beginning of classes).
 *
 *  \return             Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::cluster(StatError &error , int variable ,
                              int nb_class , vector<int> ilimit) const

{
  return cluster(error , variable , nb_class , ilimit.data());
}


/*--------------------------------------------------------------*/
/**
 *  \brief Partitioning of values of a real-valued variable.
 *
 *  \param[in] seq      reference on a Sequences object,
 *  \param[in] variable variable index,
 *  \param[in] nb_class number of classes,
 *  \param[in] limit    real limits between classes (beginning of classes).
 */
/*--------------------------------------------------------------*/

void Sequences::cluster(const Sequences &seq , int variable , int nb_class , double *limit)

{
  int i , j , k;


  // grouping of real values

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


/*--------------------------------------------------------------*/
/**
 *  \brief Partitioning of values of a real-valued variable.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] nb_class number of classes,
 *  \param[in] ilimit   real limits between classes (beginning of classes).
 *
 *  \return             Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::cluster(StatError &error , int variable ,
                              int nb_class , double *ilimit) const

{
  bool status = true;
  int i;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Partitioning of values of a real-valued variable.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] nb_class number of classes,
 *  \param[in] ilimit   real limits between classes (beginning of classes).
 *
 *  \return             Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::cluster(StatError &error , int variable ,
                              int nb_class , vector<double> ilimit) const

{
  return cluster(error , variable , nb_class , ilimit.data());
}


/*--------------------------------------------------------------*/
/**
 *  \brief Scaling of a variable.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] variable      variable index,
 *  \param[in] scaling_coeff integer scaling factor.
 *
 *  \return                  Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::scaling(StatError &error , int variable , int scaling_coeff) const

{
  bool status = true;
  int i , j;
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

    // scaling of integer values

    case INT_VALUE : {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->length[i];j++) {
          seq->int_sequence[i][variable][j] = int_sequence[i][variable][j] * scaling_coeff;
        }
      }
      break;
    }

    // scaling of real values

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


/*--------------------------------------------------------------*/
/**
 *  \brief Scaling of a variable.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] variable      variable index,
 *  \param[in] scaling_coeff real scaling factor.
 *
 *  \return                  Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::scaling(StatError &error , int variable , double scaling_coeff) const

{
  bool status = true;
  int i , j;
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

    // scaling of integer values

    case INT_VALUE : {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->length[i];j++) {
          seq->real_sequence[i][variable][j] = int_sequence[i][variable][j] * scaling_coeff;
        }
      }
      break;
    }

    // scaling of real values

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


/*--------------------------------------------------------------*/
/**
 *  \brief Rounding of values of a real-valued variable.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] mode     mode (FLOOR/ROUND/CEIL).
 *
 *  \return             Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::round(StatError &error , int variable , rounding mode) const

{
  bool status = true;
  int i , j , k;
  variable_nature *itype;
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
    itype = new variable_nature[nb_variable];

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
                        index_param_type , nb_variable , itype);
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
        for (j = 0;j < (seq->index_param_type == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }
      }
    }

    for (i = 0;i < seq->nb_sequence;i++) {
      for (j = 0;j < seq->nb_variable;j++) {

        // copy of integer values

        if ((type[j] != REAL_VALUE) && (type[j] != AUXILIARY)) {
          for (k = 0;k < seq->length[i];k++) {
            seq->int_sequence[i][j][k] = int_sequence[i][j][k];
          }
        }

        else {

          // rounding of real values

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

          // copy of real values

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


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of sequences taking values in a given range for the index parameter.
 *
 *  \param[in] error               reference on a StatError object,
 *  \param[in] display             flag for displaying the selected individuals,
 *  \param[in] min_index_parameter lowest index parameter,
 *  \param[in] max_index_parameter highest index parameter,
 *  \param[in] keep                flag for keeping or rejecting the selected sequences.
 *
 *  \return                        Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::index_parameter_select(StatError &error , bool display ,
                                             int min_index_parameter ,
                                             int max_index_parameter , bool keep) const

{
  bool status = true;
  int i , j;
  int inb_sequence , *index , *iidentifier;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((index_param_type != TIME) && (index_param_type != POSITION)) {
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

    // selection of sequences

    iidentifier = new int[nb_sequence];
    index = new int[nb_sequence];
    inb_sequence = 0;

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < (index_param_type == POSITION ? length[i] + 1 : length[i]);j++) {
        if ((index_parameter[i][j] >= min_index_parameter) &&
            (index_parameter[i][j] <= max_index_parameter)) {
          if (keep) {
            iidentifier[inb_sequence] = identifier[i];
            index[inb_sequence++] = i;
          }
          break;
        }
      }

      if ((!keep) && (j == (index_param_type == POSITION ? length[i] + 1 : length[i]))) {
        iidentifier[inb_sequence] = identifier[i];
        index[inb_sequence++] = i;
      }
    }

    if (inb_sequence == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    // copy of sequences

    if (status) {
      if ((display) && (inb_sequence <= DISPLAY_NB_INDIVIDUAL)) {
        cout << "\n" << SEQ_label[inb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << ": ";
        for (i = 0;i < inb_sequence;i++) {
          cout << iidentifier[i] << ", ";
        }
        cout << endl;
      }

      seq = new Sequences(*this , inb_sequence , index);
    }

    delete [] iidentifier;
    delete [] index;
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of sequences taking values in a given range for a variable.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] display    flag for displaying the selected individuals,
 *  \param[in] variable   variable index,
 *  \param[in] imin_value lowest integer value,
 *  \param[in] imax_value highest integer value,
 *  \param[in] keep       flag for keeping or rejecting the selected sequences.
 *
 *  \return               Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::value_select(StatError &error , bool display , int variable ,
                                   int imin_value , int imax_value , bool keep) const

{
  bool status = true;
  int i , j;
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

    // selection of sequences

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

    // copy of sequences

    if (status) {
      if ((display) && (inb_sequence <= DISPLAY_NB_INDIVIDUAL)) {
        cout << "\n" << SEQ_label[inb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << ": ";
        for (i = 0;i < inb_sequence;i++) {
          cout << iidentifier[i] << ", ";
        }
        cout << endl;
      }

      seq = new Sequences(*this , inb_sequence , index);
    }

    delete [] iidentifier;
    delete [] index;
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of sequences taking values in a given range for a real-valued variable.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] display    flag for displaying the selected individuals,
 *  \param[in] variable   variable index,
 *  \param[in] imin_value lowest real value,
 *  \param[in] imax_value highest real value,
 *  \param[in] keep       flag for keeping or rejecting the selected sequences.
 *
 *  \return               Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::value_select(StatError &error , bool display , int variable ,
                                   double imin_value , double imax_value , bool keep) const

{
  bool status = true;
  int i , j;
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

    // selection of sequences

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

    // copy of sequences

    if (status) {
      if ((display) && (inb_sequence <= DISPLAY_NB_INDIVIDUAL)) {
        cout << "\n" << SEQ_label[inb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << ": ";
        for (i = 0;i < inb_sequence;i++) {
          cout << iidentifier[i] << ", ";
        }
        cout << endl;
      }

      seq = new Sequences(*this , inb_sequence , index);
    }

    delete [] iidentifier;
    delete [] index;
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of sequences by their identifiers.
 *
 *  \param[in] error        reference on a StatError object,
 *  \param[in] inb_sequence number of sequences,
 *  \param[in] iidentifier  sequence identifiers,
 *  \param[in] keep         flag for keeping or rejecting the selected individuals.
 *
 *  \return                 Sequences object.
 */
/*--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of sequences by their identifiers.
 *
 *  \param[in] error        reference on a StatError object,
 *  \param[in] inb_sequence number of sequences,
 *  \param[in] iidentifier  sequence identifiers,
 *  \param[in] keep         flag for keeping or rejecting the selected sequences.
 *
 *  \return                 Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::select_individual(StatError &error , int inb_sequence ,
                                        vector<int> iidentifier , bool keep) const

{
  return select_individual(error , inb_sequence , iidentifier.data() , keep);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a Sequences object transforming the implicit index parameters in
 *         explicit index parameters.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::explicit_index_parameter(StatError &error) const

{
  Sequences *seq;


  error.init();

  if (index_parameter) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  else {
    seq = new Sequences(*this , EXPLICIT_INDEX_PARAMETER);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Removing of the index parameters.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::remove_index_parameter(StatError &error) const

{
  Sequences *seq;


  error.init();

  if (!index_parameter) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  else {
    seq = new Sequences(*this , REMOVE_INDEX_PARAMETER);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of variables.
 *
 *  \param[in] seq      reference on a Sequences object,
 *  \param[in] variable variable indices.
 */
/*--------------------------------------------------------------*/

void Sequences::select_variable(const Sequences &seq , int *variable)

{
  int i , j , k;


  // copy of index parameters

  if (seq.index_parameter_distribution) {
    index_parameter_distribution = new FrequencyDistribution(*(seq.index_parameter_distribution));
  }
  if (seq.index_interval) {
    index_interval = new FrequencyDistribution(*(seq.index_interval));
  }

  if (seq.index_parameter) {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < (index_param_type == POSITION ? length[i] + 1 : length[i]);j++) {
        index_parameter[i][j] = seq.index_parameter[i][j];
      }
    }
  }

  // copy of values

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


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of variables.
 *
 *  \param[in] error        reference on a StatError object,
 *  \param[in] inb_variable number of variables,
 *  \param[in] ivariable    variable indices,
 *  \param[in] keep         flag for keeping or rejecting the selected variables.
 *
 *  \return                 Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::select_variable(StatError &error , int inb_variable ,
                                      int *ivariable , bool keep) const

{
  bool status = true , *selected_variable;
  int i;
  int bnb_variable , *variable;
  variable_nature *itype;
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
      itype = new variable_nature[bnb_variable];
      for (i = 0;i < bnb_variable;i++) {
        itype[i] = type[variable[i]];
      }

      seq = new Sequences(nb_sequence , identifier , length , vertex_identifier ,
                          index_param_type , bnb_variable , itype);
      seq->select_variable(*this , variable);

      delete [] itype;
    }

    delete [] variable;
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of variables.
 *
 *  \param[in] error        reference on a StatError object,
 *  \param[in] inb_variable number of variables,
 *  \param[in] ivariable    variable indices,
 *  \param[in] keep         flag for keeping or rejecting the selected variables.
 *
 *  \return                 Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::select_variable(StatError &error , int inb_variable ,
                                      vector<int> ivariable , bool keep) const

{
  return select_variable(error , inb_variable , ivariable.data() , keep);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Summation of variables.
 *
 *  \param[in] error              reference on a StatError object,
 *  \param[in] nb_summed_variable number of variables to be summed,
 *  \param[in] variable           variable indices.
 *
 *  \return                       Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::sum_variable(StatError &error , int nb_summed_variable , int *ivariable) const

{
  bool status = true , *selected_variable;
  int i , j , k , m , n;
  int inb_variable , *copied_variable , *summed_variable;
  variable_nature *itype;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (nb_variable == 1) {
    status = false;
    error.update(STAT_error[STATR_NB_VARIABLE]);
  }

  if ((nb_summed_variable < 2) || (nb_summed_variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_NB_SUMMED_VARIABLE]);
  }

  else {
    selected_variable = new bool[nb_variable + 1];
    for (i = 1;i <= nb_variable;i++) {
      selected_variable[i] = false;
    }

    for (i = 0;i < nb_summed_variable;i++) {
      if ((ivariable[i] < 1) || (ivariable[i] > nb_variable)) {
        status = false;
        ostringstream error_message;
        error_message << ivariable[i] << ": " << STAT_error[STATR_VARIABLE_INDEX];
        error.update((error_message.str()).c_str());
      }

      else {
        if (selected_variable[ivariable[i]]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << ivariable[i] << " "
                        << STAT_error[STATR_ALREADY_SELECTED];
          error.update((error_message.str()).c_str());
        }

        else {
          selected_variable[ivariable[i]] = true;

          if ((type[ivariable[i] - 1] != INT_VALUE) && (type[ivariable[i] - 1] != REAL_VALUE)) {
            status = false;
            ostringstream error_message , correction_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << ivariable[i] << ": "
                          << STAT_error[STATR_VARIABLE_TYPE];
            correction_message << STAT_variable_word[INT_VALUE] << " or "
                               << STAT_variable_word[REAL_VALUE];
            error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
          }

          else if ((ivariable[i] < nb_variable) && (type[ivariable[i]] == AUXILIARY)) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << ivariable[i] + 1 << ": "
                          << STAT_error[STATR_VARIABLE_TYPE];
            error.update((error_message.str()).c_str());
          }
        }
      }
    }

    delete [] selected_variable;
  }

  if (status) {
    inb_variable = nb_variable - nb_summed_variable + 1;
    for (i = 0;i < nb_summed_variable;i++) {
      ivariable[i]--;
    }

    summed_variable = new int[nb_summed_variable];
    copied_variable = new int[nb_variable - nb_summed_variable];
    i = 0;
    j = 0;
    for (k = 0;k < nb_variable;k++) {
      for (m = 0;m < nb_summed_variable;m++) {
        if (ivariable[m] == k) {
          summed_variable[i++] = k;
          break;
        }
      }
      if (m == nb_summed_variable) {
        copied_variable[j++] = k;
      }
    }

    itype = new variable_nature[inb_variable];
    i = 0;
    for (j = 0;j < inb_variable;j++) {
      if (j == summed_variable[0]) {
        itype[j] = type[j];
        for (k = 1;k < nb_summed_variable;k++) {
          if (type[summed_variable[k]] == REAL_VALUE) {
           itype[j] = REAL_VALUE;
          }
        }
      }
      else {
        itype[j] = type[copied_variable[i++]];
      }
    }

    seq = new Sequences(nb_sequence , identifier , length , vertex_identifier ,
                        index_param_type , inb_variable , itype);
    delete [] itype;

    // copy of index parameters

    if (index_parameter_distribution) {
      seq->index_parameter_distribution = new FrequencyDistribution(*(index_parameter_distribution));
    }
    if (index_interval) {
      seq->index_interval = new FrequencyDistribution(*(index_interval));
    }

    if (index_parameter) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < (index_param_type == POSITION ? length[i] + 1 : length[i]);j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }
      }
    }

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        k = 0;
        for (m = 0;m < inb_variable;m++) {

          // summation of values

          if (m == summed_variable[0]) {
            switch (seq->type[m]) {

            case INT_VALUE : {
              seq->int_sequence[i][m][j] = 0.;

              for (n = 0;n < nb_summed_variable;n++) {
                switch (type[summed_variable[n]]) {
                case INT_VALUE :
                  seq->int_sequence[i][m][j] += int_sequence[i][summed_variable[n]][j];
                  break;
                case REAL_VALUE :
                  seq->int_sequence[i][m][j] += real_sequence[i][summed_variable[n]][j];
                  break;
                }
              }
              break;
            }

            case REAL_VALUE : {
              seq->real_sequence[i][m][j] = 0.;

              for (n = 0;n < nb_summed_variable;n++) {
                switch (type[summed_variable[n]]) {
                case INT_VALUE :
                  seq->real_sequence[i][m][j] += int_sequence[i][summed_variable[n]][j];
                  break;
                case REAL_VALUE :
                  seq->real_sequence[i][m][j] += real_sequence[i][summed_variable[n]][j];
                  break;
                }
              }
              break;
            }
            }
          }

          // copy of values

          else {
            switch (seq->type[m]) {
            case INT_VALUE :
              seq->int_sequence[i][m][j] = int_sequence[i][copied_variable[k++]][j];
              break;
            case REAL_VALUE :
              seq->real_sequence[i][m][j] = real_sequence[i][copied_variable[k++]][j];
              break;
            }
          }
        }
      }
    }

    i = 0;
    for (j = 0;j < inb_variable;j++) {
      if (j == summed_variable[0]) {
        seq->min_value_computation(j);
        seq->max_value_computation(j);

        seq->build_marginal_frequency_distribution(j);
      }

      else {
        seq->min_value[j] = min_value[copied_variable[i]];
        seq->max_value[j] = max_value[copied_variable[i]];

        if (marginal_distribution[copied_variable[i]]) {
          seq->marginal_distribution[j] = new FrequencyDistribution(*(marginal_distribution[copied_variable[i]]));
        }
        if (marginal_histogram[copied_variable[i]]) {
          seq->marginal_histogram[j] = new Histogram(*(marginal_histogram[copied_variable[i]]));
        }
        i++;
      }
    }

    delete [] summed_variable;
    delete [] copied_variable;
  }

  return seq;
}



/*--------------------------------------------------------------*/
/**
 *  \brief Summation of variables.
 *
 *  \param[in] error              reference on a StatError object,
 *  \param[in] nb_summed_variable number of variables to be summed,
 *  \param[in] variable           variable indices.
 *
 *  \return                       Vectors object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::sum_variable(StatError &error , int nb_summed_variable , vector<int> ivariable) const

{
  return sum_variable(error , nb_summed_variable , ivariable.data());
}


/*--------------------------------------------------------------*/
/**
 *  \brief Merging of variables of Sequences objects.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] nb_sample  number of Sequences objects,
 *  \param[in] iseq       pointer on the Sequences objects,
 *  \param[in] ref_sample reference Sequences object for the identifiers.
 *
 *  \return               Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::merge_variable(StatError &error , int nb_sample ,
                                     const Sequences **iseq , int ref_sample) const

{
  bool status = true;
  int i , j , k , m;
  int inb_variable , *iidentifier , **ivertex_identifier;
  variable_nature *itype;
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

    if ((iseq[i]->index_param_type != IMPLICIT_TYPE) &&
        (iseq[i]->index_param_type != index_param_type)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_INDEX_PARAMETER_TYPE];

      if (index_param_type == IMPLICIT_TYPE) {
        error.update((error_message.str()).c_str());
      }
      else {
        error.correction_update((error_message.str()).c_str() , SEQ_index_parameter_word[index_param_type]);
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

          if ((iseq[i]->index_param_type != IMPLICIT_TYPE) &&
              (iseq[i]->index_param_type == index_param_type)) {
            for (k = 0;k < (index_param_type == POSITION ? length[j] + 1 : length[j]);k++) {
              if (iseq[i]->index_parameter[j][k] != index_parameter[j][k]) {
                status = false;
                ostringstream error_message;
                error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                              << SEQ_label[SEQL_SEQUENCE] << " " << j + 1 << ": "
                              << SEQ_label[SEQL_INDEX] << " " << k << ": "
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

    // comparison of sequence identifiers

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

    // comparison of vertex identifiers

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

    itype = new variable_nature[inb_variable];
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
                        index_param_type , inb_variable , itype);
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
        for (j = 0;j < (index_param_type == POSITION ? length[i] + 1 : length[i]);j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }
      }
    }

    // copy of values

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


/*--------------------------------------------------------------*/
/**
 *  \brief Merging of variables of Sequences objects.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] nb_sample  number of Sequences objects,
 *  \param[in] iseq       pointer on the Sequences objects,
 *  \param[in] ref_sample reference Sequences object for the identifiers.
 *
 *  \return               Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::merge_variable(StatError &error , int nb_sample ,
                                     const vector<Sequences> iseq , int ref_sample) const

{
  int i;
  Sequences *seq;
  const Sequences **pseq;


  pseq = new const Sequences*[nb_sample];
  for (i = 0;i < nb_sample;i++) {
    pseq[i] = new Sequences(iseq[i]);
  }

  seq = merge_variable(error , nb_sample , pseq , ref_sample);

  for (i = 0;i < nb_sample;i++) {
    delete pseq[i];
  }
  delete [] pseq;

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Differences between data and residuals in order to build auxiliary variables.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] residual reference on a Sequences object.
 *
 *  \return             Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::difference_variable(StatError &error , const Sequences &residual) const

{
  bool status = true;
  int i , j , k , m;
  int offset , inb_variable;
  variable_nature *itype;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (((residual.vertex_identifier) && (!vertex_identifier)) ||
      ((!(residual.vertex_identifier)) && (vertex_identifier))) {
    status = false;
    error.update(SEQ_error[SEQR_SAMPLE_VERTEX_IDENTIFIER]);
  }

  if (residual.index_param_type != index_param_type) {
    status = false;
    ostringstream error_message;
    error_message << SEQ_error[SEQR_INDEX_PARAMETER_TYPE];

    if (index_param_type == IMPLICIT_TYPE) {
      error.update((error_message.str()).c_str());
    }
    else {
      error.correction_update((error_message.str()).c_str() , SEQ_index_parameter_word[index_param_type]);
    }
  }

  if (residual.nb_sequence != nb_sequence) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }

  else {
    for (i = 0;i < nb_sequence;i++) {
      if (residual.identifier[i] != identifier[i]) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEQUENCE] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_SEQUENCE_IDENTIFIER];
        error.update((error_message.str()).c_str());
      }

      if (residual.length[i] != length[i]) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEQUENCE] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_SEQUENCE_LENGTH];
        error.update((error_message.str()).c_str());
      }

      else {
        if ((residual.vertex_identifier) && (vertex_identifier)) {
          for (j = 0;j < length[i];j++) {
            if (residual.vertex_identifier[i][j] != vertex_identifier[i][j]) {
              status = false;
              ostringstream error_message;
              error_message << SEQ_label[SEQL_SEQUENCE] << " " << i + 1 << ": "
                            << SEQ_label[SEQL_VERTEX_IDENTIFIER] << " " << j << ": "
                            << SEQ_error[SEQR_VERTEX_IDENTIFIER];
              error.update((error_message.str()).c_str());
            }
          }
        }

        if ((residual.index_param_type != IMPLICIT_TYPE) &&
            (residual.index_param_type == index_param_type)) {
          for (j = 0;j < (index_param_type == POSITION ? length[i] + 1 : length[i]);j++) {
            if (residual.index_parameter[i][j] != index_parameter[i][j]) {
              status = false;
              ostringstream error_message;
              error_message << SEQ_label[SEQL_SEQUENCE] << " " << i + 1 << ": "
                            << SEQ_label[SEQL_INDEX] << " " << j << ": "
                            << SEQ_error[SEQR_INDEX_PARAMETER];
              error.update((error_message.str()).c_str());
            }
          }
        }
      }
    }
  }

  if (residual.nb_variable != nb_variable) {
    status = false;
    error.correction_update(STAT_error[STATR_NB_VARIABLE] , nb_variable);
  }

  else {
    if (type[0] == STATE) {
      offset = 1;

      if (nb_variable == 1) {
        status = false;
        error.update(STAT_error[STATR_NB_VARIABLE]);
      }

      if (residual.type[0] != STATE) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " 1: "
                      << STAT_error[STATR_VARIABLE_TYPE];
        error.correction_update((error_message.str()).c_str() , STAT_variable_word[STATE]);
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            if (residual.int_sequence[i][0][j] != int_sequence[i][0][j]) {
              status = false;
              ostringstream error_message , correction_message;
              error_message << SEQ_label[SEQL_SEQUENCE] << " " << i + 1 << ": "
                            << SEQ_label[SEQL_INDEX_PARAMETER] << " " << j << ": "
                            << SEQ_error[SEQR_STATE];
              correction_message << STAT_label[STATL_STATE] << " " << int_sequence[i][0][j];
              error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
            }
          }
        }
      }
    }

    else {
      offset = 0;
    }

    for (i = offset;i < nb_variable;i++) {
      if (residual.type[i] != REAL_VALUE) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                      << STAT_error[STATR_VARIABLE_TYPE];
        error.correction_update((error_message.str()).c_str() , STAT_variable_word[REAL_VALUE]);
      }
    }
  }

  if (status) {
    inb_variable = offset + (nb_variable - offset) * 2;

    itype = new variable_nature[inb_variable];

    if (type[0] == STATE) {
      itype[0] = type[0];
    }
    i = offset;
    for (j = offset;j < nb_variable;j++) {
      itype[i++] = type[j];
      itype[i++] = AUXILIARY;
    }

    seq = new Sequences(nb_sequence , identifier , length , vertex_identifier ,
                        index_param_type , inb_variable , itype);
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
        for (j = 0;j < (index_param_type == POSITION ? length[i] + 1 : length[i]);j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }
      }
    }

    // copy of variables and building of auxiliary variables

    for (i = 0;i < nb_sequence;i++) {
      if (type[0] == STATE) {
        for (j = 0;j < length[i];j++) {
          seq->int_sequence[i][0][j] = int_sequence[i][0][j];
        }
      }

      j = offset;
      for (k = offset;k < nb_variable;k++) {
        switch (type[k]) {

        case INT_VALUE : {
          for (m = 0;m < length[i];m++) {
            seq->int_sequence[i][j][m] = int_sequence[i][k][m];
          }
          j++;
          for (m = 0;m < length[i];m++) {
            seq->real_sequence[i][j][m] = int_sequence[i][k][m] - residual.real_sequence[i][k][m];
          }
          j++;
          break;
        }

        case REAL_VALUE : {
          for (m = 0;m < length[i];m++) {
            seq->real_sequence[i][j][m] = real_sequence[i][k][m];
          }
          j++;
          for (m = 0;m < length[i];m++) {
            seq->real_sequence[i][j][m] = real_sequence[i][k][m] - residual.real_sequence[i][k][m];
          }
          j++;
          break;
        }
        }
      }
    }

    if (type[0] == STATE) {
      seq->min_value[0] = min_value[0];
      seq->max_value[0] = max_value[0];
      seq->marginal_distribution[0] = new FrequencyDistribution(*marginal_distribution[0]);
    }

    i = offset;
    for (j = offset;j < nb_variable;j++) {
      seq->min_value[i] = min_value[j];
      seq->max_value[i] = max_value[j];

      if (marginal_distribution[j]) {
        seq->marginal_distribution[i] = new FrequencyDistribution(*marginal_distribution[j]);
      }
      else {
        seq->marginal_distribution[i] = NULL;
      }

      if (marginal_histogram[j]) {
        seq->marginal_histogram[i] = new Histogram(*marginal_histogram[j]);
      }
      else {
        seq->marginal_histogram[i] = NULL;
      }
      i++;

      seq->min_value_computation(i);
      seq->max_value_computation(i);
      i++;
    }
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Variable shift.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] lag      lag.
 *
 *  \return             Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::shift_variable(StatError &error , int variable , int lag) const

{
  bool status = true;
  int i , j , k;
  int inb_sequence , *iidentifier , *ilength , *index , *pvertex_id , *cvertex_id ,
      *pindex_param , *cindex_param , *pisequence , *cisequence;
  double *prsequence , *crsequence;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (index_param_type == POSITION) {
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

    // computation of the sequence lengths

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
                        index_param_type , nb_variable , type , false);

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

    // copy of values

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

          if ((j != variable) && (lag > 0)) {
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


/*--------------------------------------------------------------*/
/**
 *  \brief Reversing of the direction of sequences.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::reverse(StatError &error) const

{
  Sequences *seq;


  error.init();

  if (index_param_type == TIME) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER]);
  }

  else {
    seq = new Sequences(*this , REVERSE);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of sequences on a sequence length criterion.
 *
 *  \param[in] error       reference on a StatError object,
 *  \param[in] display     flag for displaying the selected individuals,
 *  \param[in] min_length  lowest sequence length,
 *  \param[in] imax_length highest sequence length,
 *  \param[in] keep        flag for keeping or rejecting the selected sequences.
 *
 *  \return                Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::length_select(StatError &error , bool display , int min_length ,
                                    int imax_length , bool keep) const

{
  bool status = true;
  int i;
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

    // selection of sequences

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

    // copy of selected sequences

    if (status) {
      if ((display) && (inb_sequence <= DISPLAY_NB_INDIVIDUAL)) {
        cout << "\n" << SEQ_label[inb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << ": ";
        for (i = 0;i < inb_sequence;i++) {
          cout << iidentifier[i] << ", ";
        }
        cout << endl;
      }

      seq = new Sequences(*this , inb_sequence , index);
    }

    delete [] iidentifier;
    delete [] index;
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Removing of the first/last run for a given value.
 *
 *  \param[in] error          reference on a StatError object,
 *  \param[in] variable       variable index,
 *  \param[in] ivalue         value,
 *  \param[in] position       position (BEGIN_RUN/END_RUN),
 *  \param[in] max_run_length maximum length of the removed runs.
 *
 *  \return                   Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::remove_run(StatError &error , int variable , int ivalue ,
                                 run_position position , int max_run_length) const

{
  bool status = true;
  int i , j , k;
  int smax_length , inb_sequence , *iidentifier , *ilength , *index , *pvertex_id ,
      *cvertex_id , *pindex_param , *cindex_param , *pisequence , *cisequence;
  double *prsequence , *crsequence;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (index_param_type == POSITION) {
    error.update(SEQ_error[SEQR_INDEX_PARAMETER]);
    status = false;
  }

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
      if (!marginal_distribution[variable]) {
        if ((ivalue < min_value[variable]) || (ivalue > max_value[variable])) {
          status = false;
          error.update(STAT_error[STATR_VALUE]);
        }
/*        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
        error.update((error_message.str()).c_str()); */
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

    // computation of the sequence lengths

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

      case BEGIN_RUN : {
        if (type[variable] != REAL_VALUE) {
          cisequence = int_sequence[i][variable];
          for (j = 0;j < smax_length;j++) {
            if (*cisequence++ != ivalue) {
              break;
            }
          }
        }

        else {
          crsequence = real_sequence[i][variable];
          for (j = 0;j < smax_length;j++) {
            if (*crsequence++ != (double)ivalue) {
              break;
            }
          }
        }
        break;
      }

      case END_RUN : {
        if (type[variable] != REAL_VALUE) {
          cisequence = int_sequence[i][variable] + length[i];
          for (j = 0;j < smax_length;j++) {
            if (*--cisequence != ivalue) {
              break;
            }
          }
        }

        else {
          crsequence = real_sequence[i][variable] + length[i];
          for (j = 0;j < smax_length;j++) {
            if (*--crsequence != (double)ivalue) {
              break;
            }
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
                        index_param_type , nb_variable , type , false);

    // copy of vertex identifiers

    if (vertex_identifier) {
      for (i = 0;i < seq->nb_sequence;i++) {
        pvertex_id = seq->vertex_identifier[i];

        switch (position) {
        case BEGIN_RUN :
          cvertex_id = vertex_identifier[index[i]] + length[index[i]] - seq->length[i];
          break;
        case END_RUN :
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
        case BEGIN_RUN :
          cindex_param = index_parameter[index[i]] + length[index[i]] - seq->length[i];
          break;
        case END_RUN :
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

    // copy of values

    for (i = 0;i < seq->nb_sequence;i++) {
      for (j = 0;j < seq->nb_variable;j++) {
        if ((seq->type[j] != REAL_VALUE) && (seq->type[j] != AUXILIARY)) {
          pisequence = seq->int_sequence[i][j];

          switch (position) {
          case BEGIN_RUN :
            cisequence = int_sequence[index[i]][j] + length[index[i]] - seq->length[i];
            break;
          case END_RUN :
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
          case BEGIN_RUN :
            crsequence = real_sequence[index[i]][j] + length[index[i]] - seq->length[i];
            break;
          case END_RUN :
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


/*--------------------------------------------------------------*/
/**
 *  \brief Truncation of sequences.
 *
 *  \param[in] error               reference on a StatError object,
 *  \param[in] max_index_parameter highest index parameter.
 *
 *  \return                        Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::truncate(StatError &error , int max_index_parameter) const

{
  int i , j , k;
  int *ilength;
  Sequences *seq;


  error.init();

  if ((max_index_parameter < (index_parameter ? index_parameter_distribution->offset : 1)) ||
      ((!index_parameter) && (max_index_parameter >= max_length))) {
    seq = NULL;
    error.update(SEQ_error[SEQR_MAX_INDEX_PARAMETER]);
  }

  else {
    ilength = new int[nb_sequence];

    // explicit index parameters

    if (index_parameter) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = length[i] - 1;j >= 0;j--) {
          if (index_parameter[i][j] <= max_index_parameter) {
            break;
          }
        }
        ilength[i] = MAX(j + 1 , 1);

#       ifdef MESSAGE
        if (ilength[i] == 0) {
          cout << "\n" << identifier[i] << " " << j << " " << length[i] << endl;
        }
#       endif

      }
    }

    // implicit index parameters

    else {
      for (i = 0;i < nb_sequence;i++) {
        ilength[i] = MIN(max_index_parameter + 1 , length[i]);
      }
    }

    // extraction of truncated sequences

    seq = new Sequences(nb_sequence , identifier , ilength , vertex_identifier ,
                        index_param_type , nb_variable , type);
    delete [] ilength;

    // copy of index parameters

    if (index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->length[i];j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }

        if (seq->index_param_type == POSITION) {
          seq->index_parameter[i][seq->length[i]] = (index_interval->variance == 0. ? index_parameter[i][seq->length[i]] : index_parameter[i][seq->length[i] - 1]);
        }
      }

      seq->build_index_parameter_frequency_distribution();
      seq->index_interval_computation();
    }

    // copy of values

    for (i = 0;i < seq->nb_sequence;i++) {
      for (j = 0;j < seq->nb_variable;j++) {
        if ((seq->type[j] != REAL_VALUE) && (seq->type[j] != AUXILIARY)) {
          for (k = 0;k < seq->length[i];k++) {
            seq->int_sequence[i][j][k] = int_sequence[i][j][k];
          }
        }

        else {
          for (k = 0;k < seq->length[i];k++) {
            seq->real_sequence[i][j][k] = real_sequence[i][j][k];
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
  }
 
  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of sub-sequences.
 *
 *  \param[in] error               reference on a StatError object,
 *  \param[in] min_index_parameter lowest index parameter,
 *  \param[in] max_index_parameter highest index parameter.
 *
 *  \return                        Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::index_parameter_extract(StatError &error , int min_index_parameter ,
                                              int max_index_parameter) const

{
  bool status = true;
  int i , j , k;
  int inb_sequence , *iidentifier , *ilength , *index , *first_index , *pvertex_id ,
      *cvertex_id , *pindex_param , *cindex_param , *pisequence , *cisequence;
  double *prsequence , *crsequence;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((min_index_parameter < (index_parameter ? index_parameter_distribution->offset : 1)) ||
      ((!index_parameter) && (min_index_parameter >= max_length)) ||
      ((max_index_parameter != I_DEFAULT) && (min_index_parameter > max_index_parameter))) {
    status = false;
    error.update(SEQ_error[SEQR_MIN_INDEX_PARAMETER]);
  }
  if ((max_index_parameter != I_DEFAULT) && ((max_index_parameter < (index_parameter ? index_parameter_distribution->offset : 1)) ||
       ((!index_parameter) && (max_index_parameter >= max_length)) || (max_index_parameter < min_index_parameter))) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_INDEX_PARAMETER]);
  }

  if (status) {

    // selection of sequences

    iidentifier = new int[nb_sequence];
    ilength = new int[nb_sequence];
    index = new int[nb_sequence];
    inb_sequence = 0;

    // explicit index parameters

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

    // implicit index parameters

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

    // extraction of sub-sequences

    if (status) {
      seq = new Sequences(inb_sequence , iidentifier , ilength , vertex_identifier ,
                          index_param_type , nb_variable , type , false);

      // copy of vertex identifiers

      if (vertex_identifier) {
        for (i = 0;i < seq->nb_sequence;i++) {
          pvertex_id = seq->vertex_identifier[i];
          cvertex_id = vertex_identifier[index[i]] + first_index[i];
          for (j = 0;j < seq->length[i];j++) {
            *pvertex_id++ = *cvertex_id++;
          }
        }
      }

      // copy of index parameters

      if (index_parameter) {
        for (i = 0;i < seq->nb_sequence;i++) {
          pindex_param = seq->index_parameter[i];
          cindex_param = index_parameter[index[i]] + first_index[i];
          for (j = 0;j < seq->length[i];j++) {
            *pindex_param++ = *cindex_param++;
          }

          if (seq->index_param_type == POSITION) {
            if (max_index_parameter == I_DEFAULT) {
              *pindex_param = *cindex_param;
            }
            else {
              *pindex_param = (index_interval->variance == 0. ? max_index_parameter + index_interval->mean : max_index_parameter);
            }
          }
        }
      }

      // copy of values

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


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction by segmentation of a Sequences object.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] variable      variable index,
 *  \param[in] nb_value      number of values,
 *  \param[in] ivalue        values,
 *  \param[in] keep          flag for keeping or rejecting the selected segments,
 *  \param[in] concatenation segments merged by sequence or not.
 *
 *  \return                  Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::segmentation_extract(StatError &error , int variable ,
                                           int nb_value , int *ivalue , bool keep ,
                                           bool concatenation) const

{
  bool status = true;
  int i , j , k , m , n;
  int nb_present_value , nb_selected_value , nb_segment , inb_sequence ,
      *selected_value , *pvertex_id , *cvertex_id , *pindex_param , *cindex_param ,
      *pisequence , *cisequence , *segment_length , *sequence_length , **segment_begin;
  variable_nature *itype;
  double *prsequence , *crsequence;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (index_param_type == POSITION) {
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
        nb_present_value = 0;
        for (i = marginal_distribution[variable]->offset;i < marginal_distribution[variable]->nb_value;i++) {
          if (marginal_distribution[variable]->frequency[i] > 0) {
            nb_present_value++;
          }
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
    itype = new variable_nature[nb_variable - 1];
    for (i = 0;i < variable;i++) {
      itype[i] = type[i];
    }
    for (i = variable + 1;i < nb_variable;i++) {
      itype[i - 1] = type[i];
    }

    if (keep) {
      nb_selected_value = nb_value;
      selected_value = ivalue;
    }

    else {
      nb_selected_value = nb_present_value - nb_value;
      selected_value = new int[nb_selected_value];
      i = 0;

      for (j = marginal_distribution[variable]->offset;j < marginal_distribution[variable]->nb_value;j++) {
        if (marginal_distribution[variable]->frequency[j] > 0) {
          for (k = 0;k < nb_value;k++) {
            if (ivalue[k] == j) {
              break;
            }
          }

          if (k == nb_value) {
            selected_value[i++] = j;
          }
        }
      }
    }

    // initializations

    segment_length = new int[cumul_length];

    segment_begin = new int*[cumul_length];
    for (i = 0;i < cumul_length;i++) {
      segment_begin[i] = 0;
    }

    // search for the sub-sequences to be extracted

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

    if (concatenation) {
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
    }

    else {
      inb_sequence = nb_segment;
      sequence_length = segment_length;
    }

    // construction of the Sequences object

    seq = new Sequences(inb_sequence , NULL , sequence_length , vertex_identifier ,
                        index_param_type , nb_variable - 1 , itype , false);

    // copy of vertex identifiers

    if (vertex_identifier) {
      i = -1;

      for (j = 0;j < nb_segment;j++) {
        if (concatenation) {
          if ((j == 0) || (segment_begin[j][0] != segment_begin[j - 1][0])) {
            pvertex_id = seq->vertex_identifier[++i];
          }
        }
        else {
          pvertex_id = seq->vertex_identifier[j];
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
        if (concatenation) {
          if ((j == 0) || (segment_begin[j][0] != segment_begin[j - 1][0])) {
            pindex_param = seq->index_parameter[++i];
          }
        }
        else {
          pindex_param = seq->index_parameter[j];
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

    // copy of values

    i = -1;
    for (j = 0;j < nb_segment;j++) {
      k = 0;
      for (m = 0;m < nb_variable;m++) {
        if (m != variable) {
          if ((type[m] != REAL_VALUE) && (type[m] != AUXILIARY)) {
            if (concatenation) {
              if ((j == 0) || (segment_begin[j][0] != segment_begin[j - 1][0])) {
                pisequence = seq->int_sequence[++i][k];
              }
            }
            else {
              pisequence = seq->int_sequence[j][k];
            }

            cisequence = int_sequence[segment_begin[j][0]][m] + segment_begin[j][1];
            for (n = 0;n < segment_length[j];n++) {
              *pisequence++ = *cisequence++;
            }
          }

          else {
            if (concatenation) {
              if ((j == 0) || (segment_begin[j][0] != segment_begin[j - 1][0])) {
                prsequence = seq->real_sequence[++i][k];
              }
            }
            else {
              prsequence = seq->real_sequence[j][k];
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
    if ((seq->index_param_type == TIME) && (seq->index_interval->variance > 0.)) {  // for the mango growth follow-ups
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
               << "average difference: " << average_diff << " | "
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


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction by segmentation of a Sequences object.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] variable      variable index,
 *  \param[in] nb_value      number of values,
 *  \param[in] ivalue        values,
 *  \param[in] keep          flag for keeping or rejecting the selected segments,
 *  \param[in] concatenation segments merged by sequence or not.
 *
 *  \return                  Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::segmentation_extract(StatError &error , int variable ,
                                           int nb_value , vector<int> ivalue , bool keep ,
                                           bool concatenation) const

{
  return segmentation_extract(error , variable , nb_value , ivalue.data() , keep , concatenation);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Summation of successive values along sequences.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index.
 *
 *  \return             Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::cumulate(StatError &error , int variable) const

{
  bool status = true;
  int i , j , k , m;
  int offset , inb_variable;
  variable_nature *itype;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (type[0] == STATE) {
    offset = 1;
    if (nb_variable == 1) {
      status = false;
      error.update(STAT_error[STATR_NB_VARIABLE]);
    }
  }
  else {
    offset = 0;
  }

  if (variable != I_DEFAULT) {
    if ((variable < offset + 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }
    else {
      variable--;
    }
  }

  for (i = offset;i < nb_variable;i++) {
    if (((variable == I_DEFAULT) || (variable == i)) && (type[i] != INT_VALUE) && (type[i] != REAL_VALUE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
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
      inb_variable = offset + 1;
      itype = new variable_nature[inb_variable];
      if (type[0] == STATE) {
        itype[0] = type[0];
      }
      itype[offset] = type[variable];
    }

    seq = new Sequences(nb_sequence , identifier , length , vertex_identifier ,
                        index_param_type , inb_variable , itype);

    if (variable != I_DEFAULT) {
      delete [] itype;
    }

    // copy of index parameters

    if (index_parameter_distribution) {
      seq->index_parameter_distribution = new FrequencyDistribution(*index_parameter_distribution);
    }
    if (index_interval) {
      seq->index_interval = new FrequencyDistribution(*index_interval);
    }

    if (index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < (seq->index_param_type == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }
      }
    }

    // copy of the state variable

    if (type[0] == STATE) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->length[i];j++) {
          seq->int_sequence[i][0][j] = int_sequence[i][0][j];
        }
      }

      seq->min_value[0] = min_value[0];
      seq->max_value[0] = max_value[0];
      seq->marginal_distribution[0] = new FrequencyDistribution(*marginal_distribution[0]);
    }

    // summation of values

    for (i = 0;i < nb_sequence;i++) {
      j = offset;
      for (k = offset;k < nb_variable;k++) {
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

    for (i = offset;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);

      seq->build_marginal_frequency_distribution(i);
    }
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief First-order differencing of sequences.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] variable      variable index,
 *  \param[in] first_element 1st element of the sequence kept or not.
 *
 *  \return                  Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::difference(StatError &error , int variable , bool first_element) const

{
  bool status = true;
  int i , j , k , m;
  int offset , inb_variable , *ilength , *pvertex_id , *cvertex_id , *pindex_param , *cindex_param ,
      *pisequence , *cisequence;
  variable_nature *itype;
  double *prsequence , *crsequence;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (index_param_type == POSITION) {
    status = increasing_index_parameter_checking(error , true , SEQ_label[SEQL_SEQUENCE]);
  }

  if (type[0] == STATE) {
    offset = 1;
    if (nb_variable == 1) {
      status = false;
      error.update(STAT_error[STATR_NB_VARIABLE]);
    }
  }
  else {
    offset = 0;
  }

  if (variable != I_DEFAULT) {
    if ((variable < offset + 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }
    else {
      variable--;
    }
  }

  for (i = offset;i < nb_variable;i++) {
    if (((variable == I_DEFAULT) || (variable == i)) && (type[i] != INT_VALUE) && (type[i] != REAL_VALUE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
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
      inb_variable = offset + 1;
    }

    if ((index_param_type != IMPLICIT_TYPE) && ((index_interval->mean != 1.) ||
         (index_interval->variance > 0.))) {
      itype = new variable_nature[inb_variable];
      if (type[0] == STATE) {
        itype[0] = type[0];
      }
      for (i = offset;i < inb_variable;i++) {
        itype[i] = REAL_VALUE;
      }
    }

    else {
      if (variable == I_DEFAULT) {
        itype = type;
      }
      else {
        itype = new variable_nature[inb_variable];
        if (type[0] == STATE) {
          itype[0] = type[0];
        }
        itype[offset] = type[variable];
      }
    }

    seq = new Sequences(nb_sequence , identifier , ilength , vertex_identifier ,
                        index_param_type , inb_variable , itype , false);

    if (!first_element) {
      delete [] ilength;
    }
    if (((index_param_type != IMPLICIT_TYPE) && ((index_interval->mean != 1.) ||
          (index_interval->variance > 0.))) || (variable != I_DEFAULT)) {
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

        for (j = 0;j < (seq->index_param_type == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
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

    // copy of the state variable

    if (type[0] == STATE) {
      for (i = 0;i < seq->nb_sequence;i++) {
        pisequence = seq->int_sequence[i][0];
        if (first_element) {
          cisequence = int_sequence[i][0];
        }
        else {
          cisequence = int_sequence[i][0] + 1;
        }
        for (j = 0;j < seq->length[i];j++) {
          *pisequence++ = *cisequence++;
        }
      }

      if (first_element) {
        seq->min_value[0] = min_value[0];
        seq->max_value[0] = max_value[0];

        seq->marginal_distribution[0] = new FrequencyDistribution(*marginal_distribution[0]);
      }

      else {
        seq->min_value_computation(0);
        seq->max_value_computation(0);

        seq->build_marginal_frequency_distribution(0);
      }
    }

    // differencing of sequences

    if ((index_param_type != IMPLICIT_TYPE) && ((index_interval->mean != 1.) ||
         (index_interval->variance > 0.))) {
      for (i = 0;i < nb_sequence;i++) {
        j = offset;
        for (k = offset;k < nb_variable;k++) {
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
        j = offset;
        for (k = offset;k < nb_variable;k++) {
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

    for (i = offset;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);

      seq->build_marginal_frequency_distribution(i);
    }
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Log-transform of values.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] base     base of the logarithm.
 *
 *  \return             Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::log_transform(StatError &error , int variable , log_base base) const

{
  bool status = true;
  int i , j , k , m;
  int offset , inb_variable;
  variable_nature *itype;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (type[0] == STATE) {
    offset = 1;
    if (nb_variable == 1) {
      status = false;
      error.update(STAT_error[STATR_NB_VARIABLE]);
    }
  }
  else {
    offset = 0;
  }

  if (variable != I_DEFAULT) {
    if ((variable < offset + 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }
    else {
      variable--;
    }
  }

  for (i = offset;i < nb_variable;i++) {
    if ((variable == I_DEFAULT) || (variable == i)) {
      if ((type[i] != INT_VALUE) && (type[i] != REAL_VALUE)) {
        status = false;
        ostringstream error_message , correction_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                      << STAT_error[STATR_VARIABLE_TYPE];
        correction_message << STAT_variable_word[INT_VALUE] << " or "
                           << STAT_variable_word[REAL_VALUE];
        error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
      }

      else if (min_value[i] <= 0.) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                      << STAT_error[STATR_POSITIVE_MIN_VALUE];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {
    if (variable == I_DEFAULT) {
      inb_variable = nb_variable;
    }
    else {
      inb_variable = offset + 1;
    }

    itype = new variable_nature[inb_variable];
    if (type[0] == STATE) {
      itype[0] = type[0];
    }
    for (i = offset;i < inb_variable;i++) {
      itype[i] = REAL_VALUE;
    }

    seq = new Sequences(nb_sequence , identifier , length , vertex_identifier ,
                        index_param_type , inb_variable , itype);
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
        for (j = 0;j < (seq->index_param_type == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }
      }
    }

    // copy of the state variable

    if (type[0] == STATE) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->length[i];j++) {
          seq->int_sequence[i][0][j] = int_sequence[i][0][j];
        }
      }

      seq->min_value[0] = min_value[0];
      seq->max_value[0] = max_value[0];

      seq->marginal_distribution[0] = new FrequencyDistribution(*marginal_distribution[0]);
    }

    // log-transform of values

    switch (base) {

    case NATURAL : {
      for (i = 0;i < nb_sequence;i++) {
        j = offset;
        for (k = offset;k < nb_variable;k++) {
          if ((variable == I_DEFAULT) || (variable == k)) {
            if (type[k] != REAL_VALUE) {
              for (m = 0;m < length[i];m++) {
                seq->real_sequence[i][j][m] = log(int_sequence[i][k][m]);
              }
            }

            else {
              for (m = 0;m < length[i];m++) {
                seq->real_sequence[i][j][m] = log(real_sequence[i][k][m]);
              }
            }

            j++;
          }
        }
      }
      break;
    }

    case TWO : {
      for (i = 0;i < nb_sequence;i++) {
        j = offset;
        for (k = offset;k < nb_variable;k++) {
          if ((variable == I_DEFAULT) || (variable == k)) {
            if (type[k] != REAL_VALUE) {
              for (m = 0;m < length[i];m++) {
                seq->real_sequence[i][j][m] = log2(int_sequence[i][k][m]);
              }
            }

            else {
              for (m = 0;m < length[i];m++) {
                seq->real_sequence[i][j][m] = log2(real_sequence[i][k][m]);
              }
            }

            j++;
          }
        }
      }
      break;
    }

    case TEN : {
      for (i = 0;i < nb_sequence;i++) {
        j = offset;
        for (k = offset;k < nb_variable;k++) {
          if ((variable == I_DEFAULT) || (variable == k)) {
            if (type[k] != REAL_VALUE) {
              for (m = 0;m < length[i];m++) {
                seq->real_sequence[i][j][m] = log10(int_sequence[i][k][m]);
              }
            }

            else {
              for (m = 0;m < length[i];m++) {
                seq->real_sequence[i][j][m] = log10(real_sequence[i][k][m]);
              }
            }

            j++;
          }
        }
      }
      break;
    }
    }

    for (i = offset;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);

      seq->build_marginal_histogram(i);
    }
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of relative growth rates on the basis of cumulative dimensions.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] growth_factor growth factor for computing the first relative growth rate.
 *
 *  \return                  Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::relative_growth_rate(StatError &error , double growth_factor) const

{
  bool status = true , begin;
  int i , j , k;
  int offset;
  variable_nature *itype;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (index_param_type == POSITION) {
    status = increasing_index_parameter_checking(error , true , SEQ_label[SEQL_SEQUENCE]);
  }

  if (type[0] == STATE) {
    offset = 1;
    if (nb_variable == 1) {
      status = false;
      error.update(STAT_error[STATR_NB_VARIABLE]);
    }
  }
  else {
    offset = 0;
  }

  for (i = offset;i < nb_variable;i++) {
    if ((type[i] != INT_VALUE) && (type[i] != REAL_VALUE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
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
    itype = new variable_nature[nb_variable];
    if (type[0] == STATE) {
      itype[0] = type[0];
    }
    for (i = offset;i < nb_variable;i++) {
      itype[i] = REAL_VALUE;
    }

    seq = new Sequences(nb_sequence , identifier , length , vertex_identifier ,
                        index_param_type , nb_variable , itype);
    delete [] itype;

    // copy of index parameters

    if (index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < (seq->index_param_type == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
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

    // copy of the state variable

    if (type[0] == STATE) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->length[i];j++) {
          seq->int_sequence[i][0][j] = int_sequence[i][0][j];
        }
      }

      seq->min_value[0] = min_value[0];
      seq->max_value[0] = max_value[0];

      seq->marginal_distribution[0] = new FrequencyDistribution(*marginal_distribution[0]);
    }

    // computaton of relative growth rates

    if ((index_param_type != IMPLICIT_TYPE) && ((index_interval->mean != 1.) ||
         (index_interval->variance > 0.))) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = offset;j < nb_variable;j++) {
          seq->real_sequence[i][j][0] = 0.;
          begin = true;

          if (type[j] != REAL_VALUE) {
            for (k = 1;k < length[i];k++) {
              if ((int_sequence[i][j][k] > 0) && (int_sequence[i][j][k - 1] > 0)) {
                seq->real_sequence[i][j][k] = (log(int_sequence[i][j][k]) - log(int_sequence[i][j][k - 1])) /
                                              (index_parameter[i][k] - index_parameter[i][k - 1]);

                if (begin) {
                  begin = false;
                  seq->real_sequence[i][j][k - 1] = seq->real_sequence[i][j][k] * growth_factor;
                }
              }

              else {
                seq->real_sequence[i][j][k] = 0.;
              }
            }
          }

          else {
            for (k = 1;k < length[i];k++) {
              if ((real_sequence[i][j][k] > 0.) && (real_sequence[i][j][k - 1] > 0.)) {
                seq->real_sequence[i][j][k] = (log(real_sequence[i][j][k]) - log(real_sequence[i][j][k - 1])) /
                                              (index_parameter[i][k] - index_parameter[i][k - 1]);

                if (begin) {
                  begin = false;
                  seq->real_sequence[i][j][k - 1] = seq->real_sequence[i][j][k] * growth_factor;
                }
              }

              else {
                seq->real_sequence[i][j][k] = 0.;
              }
            }
          }
        }
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        for (j = offset;j < nb_variable;j++) {
          seq->real_sequence[i][j][0] = 0.;
          begin = true;

          if (type[j] != REAL_VALUE) {
            for (k = 1;k < length[i];k++) {
              if ((int_sequence[i][j][k] > 0) && (int_sequence[i][j][k - 1] > 0)) {
                seq->real_sequence[i][j][k] = log(int_sequence[i][j][k]) - log(int_sequence[i][j][k - 1]);

                if (begin) {
                  begin = false;
                  seq->real_sequence[i][j][k - 1] = seq->real_sequence[i][j][k] * growth_factor;
                }
              }

              else {
                seq->real_sequence[i][j][k] = 0.;
              }
            }
          }

          else {
            for (k = 1;k < length[i];k++) {
              if ((real_sequence[i][j][k] > 0.) && (real_sequence[i][j][k - 1] > 0.)) {
                seq->real_sequence[i][j][k] = log(real_sequence[i][j][k]) - log(real_sequence[i][j][k - 1]);

                if (begin) {
                  begin = false;
                  seq->real_sequence[i][j][k - 1] = seq->real_sequence[i][j][k] * growth_factor;
                }
              }

              else {
                seq->real_sequence[i][j][k] = 0.;
              }
            }
          }
        }
      }
    }

    for (i = offset;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);

      seq->build_marginal_histogram(i);
    }
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Normalization of sequences (for mangoo GU growth profiles).
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index.
 *
 *  \return             Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::sequence_normalization(StatError &error , int variable) const

{
  bool status = true;
  int i , j , k;
  int offset , int_max;
  variable_nature *itype;
  double real_max;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (type[0] == STATE) {
    offset = 1;
    if (nb_variable == 1) {
      status = false;
      error.update(STAT_error[STATR_NB_VARIABLE]);
    }
  }
  else {
    offset = 0;
  }

  if (variable != I_DEFAULT) {
    if ((variable < offset + 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }
    else {
      variable--;
    }
  }

  for (i = offset;i < nb_variable;i++) {
    if ((variable == I_DEFAULT) || (variable == i)) {
      if ((type[i] != INT_VALUE) && (type[i] != REAL_VALUE)) {
        status = false;
        ostringstream error_message , correction_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                      << STAT_error[STATR_VARIABLE_TYPE];
        correction_message << STAT_variable_word[INT_VALUE] << " or "
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
    itype = new variable_nature[nb_variable];

    if (type[0] == STATE) {
      itype[0] = type[0];
    }

    if (variable == I_DEFAULT) {
      for (i = offset;i < nb_variable;i++) {
        itype[i] = REAL_VALUE;
      }
    }
    else {
      for (i = offset;i < nb_variable;i++) {
        itype[i] = type[i];
      }
      itype[variable] = REAL_VALUE;
    }

    seq = new Sequences(nb_sequence , identifier , length , vertex_identifier ,
                        index_param_type , nb_variable , itype);
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
        for (j = 0;j < (seq->index_param_type == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }
      }
    }

    // copy of the state variable

    if (type[0] == STATE) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->length[i];j++) {
          seq->int_sequence[i][0][j] = int_sequence[i][0][j];
        }
      }

      seq->min_value[0] = min_value[0];
      seq->max_value[0] = max_value[0];

      seq->marginal_distribution[0] = new FrequencyDistribution(*marginal_distribution[0]);
    }

   // normalization of sequences

    for (i = 0;i < nb_sequence;i++) {
      for (j = offset;j < nb_variable;j++) {
        if ((variable == I_DEFAULT) || (variable == j)) {
          if (type[j] != REAL_VALUE) {
            int_max = int_sequence[i][j][0];
            for (k = 1;k < length[i];k++) {
              if (int_sequence[i][j][k] > int_max) {
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
              if (real_sequence[i][j][k] > real_max) {
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

    for (i = offset;i < seq->nb_variable;i++) {
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


/*--------------------------------------------------------------*/
/**
 *  \brief Filtering of sequences using a symmetric smoothing filter.
 *
 *  \param[in] error        reference on a StatError object,
 *  \param[in] nb_point     filter half width,
 *  \param[in] filter       filter,
 *  \param[in] variable     variable index,
 *  \param[in] begin_end    begin and end kept or not,
 *  \param[in] segmentation smoothing by segment or not using the state variable,
 *  \param[in] output       trend, substraction residuals or division residuals.
 *
 *  \return                 Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::moving_average(StatError &error , int nb_point , double *filter ,
                                     int variable , bool begin_end , bool segmentation ,
                                     sequence_type output) const

{
  bool status = true;
  int i , j , k , m , n , p;
  int offset , inb_variable , nb_segment , *ilength , *pvertex_id , *cvertex_id , *pindex_param , *cindex_param ,
      *pisequence , *cisequence , *change_point;
  variable_nature *itype;
  double *prsequence , *crsequence , *pfilter;
  Sequences *seq;


  seq = NULL;
  error.init();

  if ((index_interval) && (index_interval->variance > 0.)) {
    status = false;
    error.update(SEQ_error[SEQR_UNEQUAL_INDEX_INTERVALS]);
  }

  if (type[0] == STATE) {
    offset = 1;
    if (nb_variable == 1) {
      status = false;
      error.update(STAT_error[STATR_NB_VARIABLE]);
    }
  }
  else {
    offset = 0;
  }

  if (variable != I_DEFAULT) {
    if ((variable < offset + 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }
    else {
      variable--;
    }
  }

  for (i = offset;i < nb_variable;i++) {
    if (((variable == I_DEFAULT) || (variable == i)) && (type[i] != INT_VALUE) && (type[i] != REAL_VALUE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
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
      inb_variable = nb_variable - offset;
    }
    else {
      inb_variable = 1;
    }

    if (output == SEQUENCE) {
      inb_variable *= 2;
    }
    inb_variable += offset;

    itype = new variable_nature[inb_variable];

    if (type[0] == STATE) {
      itype[0] = type[0];
    }

    if (output == SEQUENCE) {
      i = offset;
      for (j = offset;j < nb_variable;j++) {
        if ((variable == I_DEFAULT) || (variable == j)) {
          itype[i++] = type[j];
          itype[i++] = AUXILIARY;
        }
      }
    }

    else {
      for (i = offset;i < inb_variable;i++) {
        itype[i] = REAL_VALUE;
      }
    }

    seq = new Sequences(nb_sequence , identifier , ilength , vertex_identifier ,
                        index_param_type , inb_variable , itype , false);

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

        for (j = 0;j < (seq->index_param_type == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
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

    // copy of the state variable

    if (type[0] == STATE) {
      for (i = 0;i < seq->nb_sequence;i++) {
        pisequence = seq->int_sequence[i][0];
        if (begin_end) {
          cisequence = int_sequence[i][0];
        }
        else {
          cisequence = int_sequence[i][0] + nb_point;
        }

        for (j = 0;j < seq->length[i];j++) {
          *pisequence++ = *cisequence++;
        }
      }

      if (begin_end) {
        seq->min_value[0] = min_value[0];
        seq->max_value[0] = max_value[0];
        seq->marginal_distribution[0] = new FrequencyDistribution(*marginal_distribution[0]);
      }

      else {
        seq->min_value_computation(0);
        seq->max_value_computation(0);
        seq->build_marginal_frequency_distribution(0);
      }
    }

    // filtering using a symmetric smoothing filter

    if ((type[0] == STATE) && (begin_end) && (segmentation)) {
      change_point = new int[max_length];

      for (i = 0;i < nb_sequence;i++) {
        change_point[0] = 0;
        j = 1;
        for (k = 1;k < length[i];k++) {
          if (int_sequence[i][0][k] != int_sequence[i][0][k - 1]) {
            change_point[j++] = k;
          }
        }
        change_point[j] = length[i];
        nb_segment = j;

        j = 1;
        for (k = 1;k < nb_variable;k++) {
          if ((variable == I_DEFAULT) || (variable == k)) {
            prsequence = seq->real_sequence[i][output == SEQUENCE ? j + 1 : j];

            switch (type[k]) {

            case INT_VALUE : {
              for (m = 0;m < nb_segment;m++) {
                for (n = change_point[m];n < MIN(change_point[m] + nb_point , change_point[m + 1]);n++) {
                  cisequence = int_sequence[i][k] + change_point[m];
                  pfilter = filter;
                  *prsequence = 0.;
                  for (p = 0;p < 2 * nb_point + 1;p++) {
                    *prsequence += *cisequence * *pfilter++;
                    if ((n - nb_point + p >= change_point[m]) && (n - nb_point + p < change_point[m + 1] - 1)) {
                      cisequence++;
                    }
                  }
                  prsequence++;
                }

                for (n = change_point[m] + nb_point;n < change_point[m + 1] - nb_point;n++) {
                  cisequence = int_sequence[i][k] + n - nb_point;
                  pfilter = filter;
                  *prsequence = 0.;
                  for (p = 0;p < 2 * nb_point + 1;p++) {
                    *prsequence += *cisequence++ * *pfilter++;
                  }
                  prsequence++;
                }

                for (n = MAX(change_point[m + 1] - nb_point , change_point[m] + nb_point);n < change_point[m + 1];n++) {
                  cisequence = int_sequence[i][k] + n - nb_point;
                  pfilter = filter;
                  *prsequence = 0.;
                  for (p = 0;p < 2 * nb_point + 1;p++) {
                    *prsequence += *cisequence * *pfilter++;
                    if (n - nb_point + p < change_point[m + 1] - 1) {
                      cisequence++;
                    }
                  }
                  prsequence++;
                }
              }

              cisequence = int_sequence[i][k];

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
              break;
            }

            case REAL_VALUE : {
              for (m = 0;m < nb_segment;m++) {
                for (n = change_point[m];n < MIN(change_point[m] + nb_point , change_point[m + 1]);n++) {
                  crsequence = real_sequence[i][k] + change_point[m];
                  pfilter = filter;
                  *prsequence = 0.;
                  for (p = 0;p < 2 * nb_point + 1;p++) {
                    *prsequence += *crsequence * *pfilter++;
                    if ((n - nb_point + p >= change_point[m]) && (n - nb_point + p < change_point[m + 1] - 1)) {
                      crsequence++;
                    }
                  }
                  prsequence++;
                }

                for (n = change_point[m] + nb_point;n < change_point[m + 1] - nb_point;n++) {
                  crsequence = real_sequence[i][k] + n - nb_point;
                  pfilter = filter;
                  *prsequence = 0.;
                  for (p = 0;p < 2 * nb_point + 1;p++) {
                    *prsequence += *crsequence++ * *pfilter++;
                  }
                  prsequence++;
                }

                for (n = MAX(change_point[m + 1] - nb_point , change_point[m] + nb_point);n < change_point[m + 1];n++) {
                  crsequence = real_sequence[i][k] + n - nb_point;
                  pfilter = filter;
                  *prsequence = 0.;
                  for (p = 0;p < 2 * nb_point + 1;p++) {
                    *prsequence += *crsequence * *pfilter++;
                    if (n - nb_point + p < change_point[m + 1] - 1) {
                      crsequence++;
                    }
                  }
                  prsequence++;
                }
              }

              prsequence = seq->real_sequence[i][j];
              crsequence = real_sequence[i][k];

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
              break;
            }
            }

            if (output == SEQUENCE) {
              j++;
            }
            j++;
          }
        }
      }

      delete [] change_point;
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        j = offset;
        for (k = offset;k < nb_variable;k++) {
          if ((variable == I_DEFAULT) || (variable == k)) {
            prsequence = seq->real_sequence[i][output == SEQUENCE ? j + 1 : j];

            switch (type[k]) {

            case INT_VALUE : {
              if (begin_end) {
                for (m = 0;m < MIN(nb_point , length[i]);m++) {
                  cisequence = int_sequence[i][k];
                  pfilter = filter;
                  *prsequence = 0.;
                  for (n = 0;n < 2 * nb_point + 1;n++) {
                    *prsequence += *cisequence * *pfilter++;
                    if ((m - nb_point + n >= 0) && (m - nb_point + n < length[i] - 1)) {
                      cisequence++;
                    }
                  }
                  prsequence++;
                }
              }

              for (m = nb_point;m < length[i] - nb_point;m++) {
                cisequence = int_sequence[i][k] + m - nb_point;
                pfilter = filter;
                *prsequence = 0.;
                for (n = 0;n < 2 * nb_point + 1;n++) {
                  *prsequence += *cisequence++ * *pfilter++;
                }
                prsequence++;
              }

              if (begin_end) {
                for (m = MAX(length[i] - nb_point , nb_point);m < length[i];m++) {
                  cisequence = int_sequence[i][k] + m - nb_point;
                  pfilter = filter;
                  *prsequence = 0.;
                  for (n = 0;n < 2 * nb_point + 1;n++) {
                    *prsequence += *cisequence * *pfilter++;
                    if (m - nb_point + n < length[i] - 1) {
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
              break;
            }

            case REAL_VALUE : {
              if (begin_end) {
                for (m = 0;m < MIN(nb_point , length[i]);m++) {
                  crsequence = real_sequence[i][k];
                  pfilter = filter;
                  *prsequence = 0.;
                  for (n = 0;n < 2 * nb_point + 1;n++) {
                  *prsequence += *crsequence * *pfilter++;
                    if ((m - nb_point + n >= 0) && (m - nb_point + n < length[i] - 1)) {
                      crsequence++;
                    }
                  }
                  prsequence++;
                }
              }

              for (m = nb_point;m < length[i] - nb_point;m++) {
                crsequence = real_sequence[i][k] + m - nb_point;
                pfilter = filter;
                *prsequence = 0.;
                for (n = 0;n < 2 * nb_point + 1;n++) {
                  *prsequence += *crsequence++ * *pfilter++;
                }
                prsequence++;
              }

              if (begin_end) {
                for (m = MAX(length[i] - nb_point , nb_point);m < length[i];m++) {
                  crsequence = real_sequence[i][k] + m - nb_point;
                  pfilter = filter;
                  *prsequence = 0.;
                  for (n = 0;n < 2 * nb_point + 1;n++) {
                    *prsequence += *crsequence * *pfilter++;
                    if (m - nb_point + n < length[i] - 1) {
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
              break;
            }
            }

            if (output == SEQUENCE) {
              j++;
            }
            j++;
          }
        }
      }
    }

    if (output == SEQUENCE) {
      i = offset;

      if (begin_end) {
        for (j = offset;j < nb_variable;j++) {
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
        for (j = offset;j < nb_variable;j++) {
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
      for (i = offset;i < seq->nb_variable;i++) {
        seq->min_value_computation(i);
        seq->max_value_computation(i);

        seq->build_marginal_histogram(i);
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Filtering of sequences using a symmetric smoothing filter.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] nb_point  filter half width,
 *  \param[in] filter    filter,
 *  \param[in] variable  variable index,
 *  \param[in] begin_end begin and end kept or not,
 *  \param[in] state     smoothing by segment or not using the state variable,
 *  \param[in] output    trend, substraction residuals or division residuals.
 *
 *  \return              Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::moving_average(StatError &error , int nb_point , vector<double> filter ,
                                     int variable , bool begin_end , bool segmentation ,
                                     sequence_type output) const

{
  return moving_average(error , nb_point , filter.data() , variable , begin_end , segmentation , output);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Filtering of sequences using a symmetric smoothing filter.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] dist      symmetric discrete distribution,
 *  \param[in] variable  variable index,
 *  \param[in] begin_end begin and end kept or not,
 *  \param[in] state     smoothing by segment or not using the state variable,
 *  \param[in] output    trend, substraction residuals or division residuals.
 *
 *  \return              Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::moving_average(StatError &error , const Distribution &dist ,
                                     int variable , bool begin_end , bool segmentation ,
                                     sequence_type output) const

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
                         begin_end , segmentation , output);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of pointwise mean, median or mean direction of sequences and
 *         associated dispersion measures.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] frequency  frequencies for successive index parameter values,
 *  \param[in] dispersion flag computation of dispersion measures,
 *  \param[in] output     output (sequences, residuals or standardized residuals).
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool Sequences::pointwise_average_ascii_print(StatError &error , const string path ,
                                              int *frequency , bool dispersion ,
                                              sequence_type output) const

{
  bool status;
  int i , j , k , m;
  int buff , inb_sequence , *width;
  double standard_normal_value , half_confidence_interval , *t_value;
  ios_base::fmtflags format_flags;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    format_flags = out_file.setf(ios::right , ios::adjustfield);

    if (dispersion) {

#     ifdef MESSAGE
      normal normal_dist;
      standard_normal_value = quantile(complement(normal_dist , 0.025));
      cout << "\nTEST: " << standard_normal_value;
#     endif

      t_value = new double[length[0]];
      for (i = 0;i < length[0];i++) {
        if (frequency[i] > 1) {
//          t_value[i] = t_value_computation(false , frequency[i] - 1 , 0.05);
          students_t students_dist(frequency[i] - 1);
          t_value[i] = quantile(complement(students_dist , 0.025));
        }
      }

#     ifdef MESSAGE
      cout << " | " << t_value[0] << " " << frequency[0] << endl;
#     endif

    }

    else {
      t_value = NULL;
    }

    // computation of the column widths

    inb_sequence = nb_sequence;
    if (dispersion) {
      inb_sequence--;
    }

    width = new int[2 * nb_variable + 2];

    for (i = 0;i < nb_variable;i++) {
      width[i] = column_width(length[0] , real_sequence[0][i]);
      if (dispersion) {
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

      if (index_param_type == TIME) {
        out_file << "\n" << SEQ_label[SEQL_TIME];
      }
      else {
        out_file << "\n" << SEQ_label[SEQL_INDEX];
      }

      out_file << "   " << STAT_label[STATL_MEAN];
      if (dispersion) {
        out_file << "   " << STAT_label[STATL_MEAN_CONFIDENCE_INTERVAL]
                 << "   " << STAT_label[STATL_STANDARD_DEVIATION];
      }
      out_file << "   " << STAT_label[STATL_FREQUENCY];

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

        if (dispersion) {
          if (frequency[j] > 1) {
//            half_confidence_interval = standard_normal_value * real_sequence[nb_sequence - 1][i][j] / sqrt((double)frequency[j]);
            half_confidence_interval = t_value[j] * real_sequence[nb_sequence - 1][i][j] / sqrt((double)frequency[j]);
            out_file << setw(width[i]) << real_sequence[0][i][j] - half_confidence_interval
                     << setw(width[i]) << real_sequence[0][i][j] + half_confidence_interval;
          }

          else {
            out_file << setw(width[i]) << " "
                     << setw(width[i]) << " ";
          }

          out_file << setw(width[i]) << real_sequence[nb_sequence - 1][i][j];
        }

        out_file << setw(width[nb_variable + 1]) << frequency[j];
 
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

    out_file.setf(format_flags , ios::adjustfield);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the pointwise mean, median or mean direction of sequences and
 *         associated dispersion measures at the spreadsheet format.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] frequency  frequencies for successive index parameter values,
 *  \param[in] dispersion flag computation of dispersion measures,
 *  \param[in] output     output (sequences, residuals or standardized residuals).
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool Sequences::pointwise_average_spreadsheet_print(StatError &error , const string path ,
                                                    int *frequency , bool dispersion ,
                                                    sequence_type output) const

{
  bool status;
  int i , j , k , m;
  int inb_sequence;
  double standard_normal_value , half_confidence_interval , *t_value;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    if (dispersion) {
//      normal normal_dist;
//      standard_normal_value = quantile(complement(normal_dist , 0.025));

      t_value = new double[length[0]];
      for (i = 0;i < length[0];i++) {
        if (frequency[i] > 1) {
//          t_value[i] = t_value_computation(false , frequency[i] - 1 , 0.05);
          students_t students_dist(frequency[i] - 1);
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

      if (index_param_type == TIME) {
        out_file << "\n" << SEQ_label[SEQL_TIME];
      }
      else {
        out_file << "\n" << SEQ_label[SEQL_INDEX];
      }

      out_file << "\t" << STAT_label[STATL_MEAN];
      if (dispersion) {
        out_file << "\t" << STAT_label[STATL_MEAN_CONFIDENCE_INTERVAL]
                 << "\t\t" << STAT_label[STATL_STANDARD_DEVIATION];
      }
      out_file << "\t" << STAT_label[STATL_FREQUENCY];

      if (nb_sequence < POINTWISE_AVERAGE_NB_SEQUENCE) {
        inb_sequence = nb_sequence;
        if (dispersion) {
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

        if (dispersion) {
          if (frequency[j] > 1) {
//            half_confidence_interval = standard_normal_value * real_sequence[nb_sequence - 1][i][j] / sqrt((double)frequency[j]);
            half_confidence_interval = t_value[j] * real_sequence[nb_sequence - 1][i][j] / sqrt((double)frequency[j]);
            out_file << "\t" << real_sequence[0][i][j] - half_confidence_interval
                     << "\t" << real_sequence[0][i][j] + half_confidence_interval;
          }

          else {
            out_file << "\t\t";
          }

          out_file << "\t" << real_sequence[nb_sequence - 1][i][j];
        }

        out_file << "\t" << frequency[j];
 
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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the pointwise mean, median or mean direction of sequences and
 *         associated dispersion measures.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] circular   flag circular variables,
 *  \param[in] robust     flag computation of robust location and dispersion measures,
 *  \param[in] dispersion flag computation of dispersion measures,
 *  \param[in] output     output (sequences, residuals or standardized residuals),
 *  \param[in] path       file path,
 *  \param[in] format     format (ASCII/SPREADSHEET).
 *
 *  \return               Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::pointwise_average(StatError &error , bool circular , bool robust ,
                                        bool dispersion , sequence_type output ,
                                        const string path , output_format format) const

{
  bool status = true;
  int i , j , k;
  int inb_sequence , min_identifier , max_identifier , *iidentifier , *ilength ,
      *pindex_param , *frequency , *index , *int_sample , *pisample;
  variable_nature *itype;
  angle_unit unit;
  double diff , *prsequence , *plocation , *pdispersion , *real_sample , *prsample ,
         *pmean_direction1 , *pmean_direction2 , ***mean_direction;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (nb_sequence == 1) {
    status = false;
    error.update(SEQ_error[SEQR_SINGLE_SEQUENCE]);
  }
  if ((index_param_type != IMPLICIT_TYPE) && (index_param_type != TIME)) {
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
    if ((output == STANDARDIZED_RESIDUAL) && (!dispersion)) {
      dispersion = true;
    }

    inb_sequence = nb_sequence + (dispersion ? 2 : 1);

    iidentifier = new int[inb_sequence];

    min_identifier = identifier[0];
    for (i = 0;i < nb_sequence;i++) {
      if (identifier[i] < min_identifier) {
        min_identifier = identifier[i];
      }

      iidentifier[i + 1] = identifier[i];
    }

    iidentifier[0] = min_identifier - 1;

    if (dispersion) {
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

    if (dispersion) {
      ilength[inb_sequence - 1] = ilength[0];
    }

    itype = new variable_nature[nb_variable];
    for (i = 0;i < nb_variable;i++) {
      itype[i] = REAL_VALUE;
    }

    seq = new Sequences(inb_sequence , iidentifier , ilength , NULL ,
                        index_param_type , nb_variable , itype);

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

      if (dispersion) {
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

    // computation of frequencies for each index parameter value

    frequency = new int[seq->length[0]];

    if (index_parameter) {
      pindex_param = seq->index_parameter[0];
      i = 0;
      for (j = index_parameter_distribution->offset;j < index_parameter_distribution->nb_value;j++) {
        if (index_parameter_distribution->frequency[j] > 0) {
          frequency[i++] = index_parameter_distribution->frequency[j];
        }
      }
    }

    else {
      frequency[0] = nb_sequence;
      for (i = 1;i < max_length;i++) {
        frequency[i] = frequency[i - 1] - length_distribution->frequency[i];
      }
    }

    if (robust) {
      if (index_parameter) {
        index = new int[nb_sequence];
      }
      else {
        index = NULL;
      }
      int_sample = new int[nb_sequence];
      real_sample = new double[nb_sequence];
    }

    if (circular) {

      // choice of the angle unit

      unit = RADIAN;
      for (i = 0;i < nb_variable;i++) {
        if (max_value[i] - min_value[i] > 2 * M_PI) {
          unit = DEGREE;
          break;
        }
      }

      // computation of mean directions

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
          mean_direction[i][0][j] /= frequency[j];
          mean_direction[i][1][j] /= frequency[j];

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

      // computation of circular standard deviations

      if (dispersion) {
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
    }

    else {
      if (robust) {

        // computation of medians

        if (index_parameter) {
          for (i = 0;i < nb_variable;i++) {
            for (j = 0;j < nb_sequence;j++) {
              index[j] = 0;
            }

            if (type[i] != REAL_VALUE) {
              for (j = 0;j < seq->length[0];j++) {
                pisample = int_sample;
                for (k = 0;k < nb_sequence;k++) {
                  while ((index[k] < length[k]) && (index_parameter[k][index[k]] < seq->index_parameter[0][j])) {
                    index[k]++;
                  }
                  if ((index[k] < length[k]) && (index_parameter[k][index[k]] == seq->index_parameter[0][j])) {
                    *pisample++ = int_sequence[k][i][index[k]];
                  }
                }

                seq->real_sequence[0][i][j] = quantile_computation(frequency[j] , int_sample , 0.5);
              }
            }

            else {
              for (j = 0;j < seq->length[0];j++) {
                prsample = real_sample;
                for (k = 0;k < nb_sequence;k++) {
                  while ((index[k] < length[k]) && (index_parameter[k][index[k]] < seq->index_parameter[0][j])) {
                    index[k]++;
                  }
                  if ((index[k] < length[k]) && (index_parameter[k][index[k]] == seq->index_parameter[0][j])) {
                    *prsample++ = real_sequence[k][i][index[k]];
                  }
                }

                seq->real_sequence[0][i][j] = quantile_computation(frequency[j] , real_sample , 0.5);
              }
            }
          }
        }

        else {
          for (i = 0;i < nb_variable;i++) {
            if (type[i] != REAL_VALUE) {
              for (j = 0;j < max_length;j++) {
                pisample = int_sample;
                for (k = 0;k < nb_sequence;k++) {
                  if (j < length[k]) {
                    *pisample++ = int_sequence[k][i][j];
                  }
                }

                seq->real_sequence[0][i][j] = quantile_computation(frequency[j] , int_sample , 0.5);
              }
            }

            else {
              for (j = 0;j < max_length;j++) {
                prsample = real_sample;
                for (k = 0;k < nb_sequence;k++) {
                  if (j < length[k]) {
                    *prsample++ = real_sequence[k][i][j];
                  }
                }

                seq->real_sequence[0][i][j] = quantile_computation(frequency[j] , real_sample , 0.5);
              }
            }
          }
        }

        // computation of mean absolute deviations from the median

        if (dispersion) {
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
                plocation = seq->real_sequence[0][j];

                if (type[j] != REAL_VALUE) {
                  for (k = 0;k < length[i];k++) {
                    while (*pindex_param < index_parameter[i][k]) {
                      pindex_param++;
                      prsequence++;
                      plocation++;
                    }
                    pindex_param++;
                    *prsequence++ += fabs(int_sequence[i][j][k] - *plocation++);
                  }
                }

                else {
                  for (k = 0;k < length[i];k++) {
                    while (*pindex_param < index_parameter[i][k]) {
                      pindex_param++;
                      prsequence++;
                      plocation++;
                    }
                    pindex_param++;
                    *prsequence++ += fabs(real_sequence[i][j][k] - *plocation++);
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
                    seq->real_sequence[seq->nb_sequence - 1][j][k] += fabs(int_sequence[i][j][k] - seq->real_sequence[0][j][k]);
                  }
                }

                else {
                  for (k = 0;k < length[i];k++) {
                    seq->real_sequence[seq->nb_sequence - 1][j][k] += fabs(real_sequence[i][j][k] - seq->real_sequence[0][j][k]);
                  }
                }
              }
            }
          }

          for (i = 0;i < seq->nb_variable;i++) {
            for (j = 0;j < seq->length[seq->nb_sequence - 1];j++) {
              if (frequency[j] > 1) {
                seq->real_sequence[seq->nb_sequence - 1][i][j] = seq->real_sequence[seq->nb_sequence - 1][i][j] /
                                                                 (frequency[j] - 1);
              }
            }
          }
        }
      }

      else {

        // computation of means

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
            seq->real_sequence[0][i][j] /= frequency[j];
          }
        }

        // computation of standard deviations

        if (dispersion) {
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
                plocation = seq->real_sequence[0][j];

                if (type[j] != REAL_VALUE) {
                  for (k = 0;k < length[i];k++) {
                    while (*pindex_param < index_parameter[i][k]) {
                      pindex_param++;
                      prsequence++;
                      plocation++;
                    }
                    pindex_param++;
                    diff = int_sequence[i][j][k] - *plocation++;
                    *prsequence++ += diff * diff;
                  }
                }

                else {
                  for (k = 0;k < length[i];k++) {
                    while (*pindex_param < index_parameter[i][k]) {
                      pindex_param++;
                      prsequence++;
                      plocation++;
                    }
                    pindex_param++;
                    diff = real_sequence[i][j][k] - *plocation++;
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
              if (frequency[j] > 1) {
                seq->real_sequence[seq->nb_sequence - 1][i][j] = sqrt(seq->real_sequence[seq->nb_sequence - 1][i][j] /
                                                                      (frequency[j] - 1));
              }
            }
          }
        }
      }
    }

    switch (output) {

    // copy of sequences

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

    // computation of residuals

    case SUBTRACTION_RESIDUAL : {
      if (circular) {
        if (index_parameter) {
          for (i = 0;i < nb_sequence;i++) {
            for (j = 0;j < nb_variable;j++) {
              pindex_param = seq->index_parameter[0];
              plocation = seq->real_sequence[0][j];

              if (type[j] != REAL_VALUE) {
                for (k = 0;k < length[i];k++) {
                  while (*pindex_param < index_parameter[i][k]) {
                    pindex_param++;
                    plocation++;
                  }

                  pindex_param++;
                  if (fabs(int_sequence[i][j][k] - *plocation) <= 180) {
                    seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k] - *plocation++;
                  }
                  else if (int_sequence[i][j][k] - *plocation > 180) {
                    seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k] - *plocation++ - 360;
                  }
                  else {
                    seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k] - *plocation++ + 360;
                  }
                }
              }

              else {
                switch (unit) {

                case DEGREE : {
                  for (k = 0;k < length[i];k++) {
                    while (*pindex_param < index_parameter[i][k]) {
                      pindex_param++;
                      plocation++;
                    }

                    pindex_param++;
                    if (fabs(int_sequence[i][j][k] - *plocation) <= 180) {
                      seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k] - *plocation++;
                    }
                    else if (int_sequence[i][j][k] - *plocation > 180) {
                      seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k] - *plocation++ - 360;
                    }
                    else {
                      seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k] - *plocation++ + 360;
                    }
                  }
                  break;
                }

                case RADIAN : {
                  for (k = 0;k < length[i];k++) {
                    while (*pindex_param < index_parameter[i][k]) {
                      pindex_param++;
                      plocation++;
                    }

                    pindex_param++;
                    if (fabs(real_sequence[i][j][k] - *plocation) <= M_PI) {
                      seq->real_sequence[i + 1][j][k] = real_sequence[i][j][k] - *plocation++;
                    }
                    else if (real_sequence[i][j][k] - *plocation > M_PI) {
                      seq->real_sequence[i + 1][j][k] = real_sequence[i][j][k] - *plocation++ - 2 * M_PI;
                    }
                    else {
                      seq->real_sequence[i + 1][j][k] = real_sequence[i][j][k] - *plocation++ + 2 * M_PI;
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
      }

      else {
        if (index_parameter) {
          for (i = 0;i < nb_sequence;i++) {
            for (j = 0;j < nb_variable;j++) {
              pindex_param = seq->index_parameter[0];
              plocation = seq->real_sequence[0][j];

              if (type[j] != REAL_VALUE) {
                for (k = 0;k < length[i];k++) {
                  while (*pindex_param < index_parameter[i][k]) {
                    pindex_param++;
                    plocation++;
                  }
                  pindex_param++;
                  seq->real_sequence[i + 1][j][k] = int_sequence[i][j][k] - *plocation++;
                }
              }

              else {
                for (k = 0;k < length[i];k++) {
                  while (*pindex_param < index_parameter[i][k]) {
                    pindex_param++;
                    plocation++;
                  }
                  pindex_param++;
                  seq->real_sequence[i + 1][j][k] = real_sequence[i][j][k] - *plocation++;
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
      }
      break;
    }

    // computation of standardized residuals

    case STANDARDIZED_RESIDUAL : {
      if (index_parameter) {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < nb_variable;j++) {
            pindex_param = seq->index_parameter[0];
            plocation = seq->real_sequence[0][j];
            pdispersion = seq->real_sequence[seq->nb_sequence - 1][j];

            if (type[j] != REAL_VALUE) {
              for (k = 0;k < length[i];k++) {
                while (*pindex_param < index_parameter[i][k]) {
                  pindex_param++;
                  plocation++;
                  pdispersion++;
                }
                pindex_param++;
                if (*pdispersion > 0.) {
                  seq->real_sequence[i + 1][j][k] = (int_sequence[i][j][k] - *plocation) /
                                                    *pdispersion;
                }
                else {
                  seq->real_sequence[i + 1][j][k] = 0.;
                }
                plocation++;
                pdispersion++;
              }
            }

            else {
              for (k = 0;k < length[i];k++) {
                while (*pindex_param < index_parameter[i][k]) {
                  pindex_param++;
                  plocation++;
                  pdispersion++;
                }
                pindex_param++;
                if (*pdispersion > 0.) {
                  seq->real_sequence[i + 1][j][k] = (real_sequence[i][j][k] - *plocation) /
                                                    *pdispersion;
                }
                else {
                  seq->real_sequence[i + 1][j][k] = 0.;
                }
                plocation++;
                pdispersion++;
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

    if ((output == SEQUENCE) && (!dispersion)) {
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

    // writing of pointwise mean, median or mean direction of sequences and associated dispersion measures

    if (!path.empty()) {
      switch (format) {
      case ASCII :
        status = seq->pointwise_average_ascii_print(error , path , frequency ,
                                                    dispersion , output);
        break;
      case SPREADSHEET :
        status = seq->pointwise_average_spreadsheet_print(error , path , frequency ,
                                                          dispersion , output);
        break;
      }

      if (!status) {

#       ifdef MESSAGE
        cout << error;
#       endif

      }
    }

    delete [] frequency;

    if (robust) {
      delete [] index;
      delete [] int_sample;
      delete [] real_sample;
    }
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the recurrence time sequences for a value taken by
 *         an integer-valued variable.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] value    value.
 *
 *  \return             Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::recurrence_time_sequences(StatError &error , int variable , int value) const

{
  bool status = true;
  int i , j;
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

    // computation of the recurrence time sequences

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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of sojourn time sequences for an integer-valued variable.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index.
 *
 *  \return             Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::sojourn_time_sequences(StatError &error , int variable) const

{
  bool status = true;
  int i , j;
  int ilength , begin_run , *pstate , *psequence;
  variable_nature itype[2];
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

    // computation of sojourn time sequences

    if ((index_param_type == TIME) && (index_interval->variance > 0.)) {  // for the mango growth follow-ups and
      for (i = 0;i < nb_sequence;i++) {                                   // the Arabidopsis rosettes
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


/*--------------------------------------------------------------*/
/**
 *  \brief Discretization of positions.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] step  discretization step.
 *
 *  \return          Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::transform_position(StatError &error , int step) const

{
  bool status = true;
  int i , j , k , m;
  int inter_position , nb_unit , *ilength , **pisequence;
  Sequences *seq;


  seq = NULL;
  error.init();

  if (index_param_type != POSITION) {
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

    // extraction of sequences

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


/*--------------------------------------------------------------*/
/**
 *  \brief Crossing of sequences.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::cross(StatError &error) const

{
  bool status = true;
  int i , j , k , m;
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

    // construction of crossed sequences

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


};  // namespace sequence_analysis
