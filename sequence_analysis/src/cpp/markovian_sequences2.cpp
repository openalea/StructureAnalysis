/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2017 CIRAD/INRA/Inria Virtual Plants
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



#include <string>
// #include <vector>
#include <sstream>
#include <iomanip>

#include <boost/math/distributions/normal.hpp>

#include "stat_tool/stat_label.h"

#include "stat_tool/quantile_computation.hpp"

#include "sequences.h"
#include "sequence_label.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a MarkovianSequences object.
 *
 *  \param[in,out] os           stream,
 *  \param[in]     exhaustive   flag detail level,
 *  \param[in]     comment_flag flag file.
 */
/*--------------------------------------------------------------*/

ostream& MarkovianSequences::ascii_write(ostream &os , bool exhaustive , bool comment_flag) const

{
  int i , j , k;
  int *int_value , *pint_value;
  double mean , variance , median , lower_quartile , upper_quartile , *real_value , *preal_value;


  if (index_param_type == TIME) {
    os << SEQ_word[SEQW_INDEX_PARAMETER] << " : "
       << SEQ_index_parameter_word[index_param_type] << "   ";
    if (comment_flag) {
      os << "# ";
    }
    os << "(" << SEQ_label[SEQL_MIN_INDEX_PARAMETER] << ": " << index_parameter_distribution->offset << ", "
       << SEQ_label[SEQL_MAX_INDEX_PARAMETER] << ": " << index_parameter_distribution->nb_value - 1 << ")" << endl;

    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    index_parameter_distribution->ascii_characteristic_print(os , false , comment_flag);

    if (exhaustive) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << "   | " << SEQ_label[SEQL_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      index_parameter_distribution->ascii_print(os , comment_flag);
    }

    if (index_interval) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_TIME_INTERVAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
      index_interval->ascii_characteristic_print(os , false , comment_flag);

      if (exhaustive) {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << "   | " << SEQ_label[SEQL_TIME_INTERVAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
        index_interval->ascii_print(os , comment_flag);
      }
    }

    os << "\n";
  }

  os << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;

  if ((self_transition) && (exhaustive)) {
    for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
      if (self_transition[i]) {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << "   | " << STAT_label[STATL_STATE] << " " << i << " - "
           << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION] << endl;

        self_transition[i]->ascii_print(os , comment_flag);
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    os << "\n" << STAT_word[STATW_VARIABLE] << " " << i + 1 << " : "
       << STAT_variable_word[type[i]];

    if (type[i] != AUXILIARY) {
      os  << "   ";
      if (comment_flag) {
        os << "# ";
      }

      if (type[i] == STATE) {
        os << "(" << marginal_distribution[i]->nb_value << " "
           << STAT_label[marginal_distribution[i]->nb_value == 1 ? STATL_STATE : STATL_STATES] << ")" << endl;
      }
      else {
        os << "(" << STAT_label[STATL_MIN_VALUE] << ": " << min_value[i] << ", "
           << STAT_label[STATL_MAX_VALUE] << ": " << max_value[i] << ")" << endl;
      }

      if (marginal_distribution[i]) {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << STAT_label[type[i] == STATE ? STATL_STATE : STATL_MARGINAL] << " "
           << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
        marginal_distribution[i]->ascii_characteristic_print(os , false , comment_flag);

        if ((marginal_distribution[i]->nb_value <= ASCII_NB_VALUE) || (exhaustive)) {
          os << "\n";
          if (comment_flag) {
            os << "# ";
          }
          os << "   | " << STAT_label[STATL_FREQUENCY] << endl;
          marginal_distribution[i]->ascii_print(os , comment_flag);
        }

#       ifdef DEBUG
        ContinuousParametric *dist;
;
        marginal_histogram[i] = new Histogram(*marginal_distribution[i]);

        dist = new ContinuousParametric(VON_MISES , 137.5 , 180 * 180 / (50. * 50. * M_PI * M_PI) , DEGREE);
        os << "\n";
        dist->ascii_parameter_print(os);
        dist->ascii_print(os , comment_flag , true , NULL , marginal_distribution[i]);
        os << "\n";
        dist->ascii_print(os , comment_flag , true , marginal_histogram[i]);
        delete dist;

        dist = new ContinuousParametric(GAUSSIAN , 137.5 , 50.);
        os << "\n";
        dist->ascii_parameter_print(os);
        dist->ascii_print(os , comment_flag , true , NULL , marginal_distribution[i]);
//        os << "\n";
//        dist->ascii_print(os , comment_flag , true , marginal_histogram[i]);
        delete dist;

        delete marginal_histogram[i];
        marginal_histogram[i] = NULL;
#       endif

      }

      else {
        mean = mean_computation(i);
        variance = variance_computation(i , mean);

        if (variance > 0.) {
          switch (type[i]) {

          case INT_VALUE : {
            int_value = new int[cumul_length];
            pint_value = int_value;
            for (j = 0;j < nb_sequence;j++) {
              for (k = 0;k < length[j];k++) {
                *pint_value++ = int_sequence[j][i][k];
              }
            }

            lower_quartile = quantile_computation(cumul_length , int_value , 0.25);
            median = quantile_computation(cumul_length , int_value , 0.5);
            upper_quartile = quantile_computation(cumul_length , int_value , 0.75);

            delete [] int_value;
            break;
          }

          case REAL_VALUE : {
            real_value = new double[cumul_length];
            preal_value = real_value;
            for (j = 0;j < nb_sequence;j++) {
              for (k = 0;k < length[j];k++) {
                *preal_value++ = real_sequence[j][i][k];
              }
            }

            lower_quartile = quantile_computation(cumul_length , real_value , 0.25);
            median = quantile_computation(cumul_length , real_value , 0.5);
            upper_quartile = quantile_computation(cumul_length , real_value , 0.75);

            delete [] real_value;
            break;
          }
          }
        }

        else {
          median = mean;
        }

        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_SAMPLE_SIZE] << ": " << cumul_length << endl;

        if (comment_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_MEAN] << ": " << mean << "   "
           << STAT_label[STATL_MEDIAN] << ": " << median << endl;

        if (comment_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
           << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance);
        if (variance > 0.) {
          os << "   " << STAT_label[STATL_LOWER_QUARTILE] << ": " << lower_quartile
             << "   " << STAT_label[STATL_UPPER_QUARTILE] << ": " << upper_quartile;
        }
        os << endl;

        if ((variance > 0.) && (exhaustive)) {
          if (comment_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << skewness_computation(i , mean , variance) << "   "
             << STAT_label[STATL_KURTOSIS_COEFF] << ": " << kurtosis_computation(i , mean , variance) << endl;
        }

        if ((marginal_histogram[i]) && (exhaustive)) {
          os << "\n";
          if (comment_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM] << endl;

          os << "\n";
          if (comment_flag) {
            os << "# ";
          }
          os << " " << STAT_label[STATL_VALUE] << "  | " << STAT_label[STATL_FREQUENCY] << endl;
          marginal_histogram[i]->ascii_print(os , comment_flag);
        }
      }

      if (characteristics[i]) {
        characteristics[i]->ascii_print(os , type[i] , *length_distribution , exhaustive , comment_flag);
      }
    }

    else {

#     ifdef MESSAGE
      mean = mean_computation(i);
      variance = variance_computation(i , mean);

      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_MEAN] << ": " << mean << "   "
         << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
         << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance) << endl;
#     endif

//      os << endl;
    }
  }

  os << "\n";
  if (comment_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
  length_distribution->ascii_characteristic_print(os , false , comment_flag);

  if (exhaustive) {
    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << "   | " << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    length_distribution->ascii_print(os , comment_flag);
  }

  os << "\n";
  if (comment_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_CUMUL_LENGTH] << ": " << cumul_length << endl;

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a MarkovianSequences object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& MarkovianSequences::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , exhaustive , false);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a MarkovianSequences object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool MarkovianSequences::ascii_write(StatError &error , const string path ,
                                     bool exhaustive) const

{
  bool status;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    ascii_write(out_file , exhaustive , false);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a MarkovianSequences object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     format     format (LINE/COLUMN/VECTOR/POSTERIOR_PROBABILITY),
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& MarkovianSequences::ascii_data_write(ostream &os , output_sequence_format format ,
                                              bool exhaustive) const

{
  ascii_write(os , exhaustive , false);
  ascii_print(os , format , false);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a MarkovianSequences object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] format     format (LINE/COLUMN/VECTOR/POSTERIOR_PROBABILITY),
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool MarkovianSequences::ascii_data_write(StatError &error , const string path ,
                                          output_sequence_format format , bool exhaustive) const

{
  bool status = false;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    if (format != 'a') {
      ascii_write(out_file , exhaustive , true);
    }
    ascii_print(out_file , format , true);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a MarkovianSequences object in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool MarkovianSequences::spreadsheet_write(StatError &error , const string path) const

{
  bool status;
  int i , j , k;
  int *int_value , *pint_value;
  double mean , variance , median , lower_quartile , upper_quartile , *real_value , *preal_value;
  Curves *smoothed_curves;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    if (index_param_type == TIME) {
      out_file << SEQ_word[SEQW_INDEX_PARAMETER] << "\t"
               << SEQ_index_parameter_word[index_param_type] << "\t\t"
               << SEQ_label[SEQL_MIN_INDEX_PARAMETER] << "\t" << index_parameter_distribution->offset << "\t\t"
               << SEQ_label[SEQL_MAX_INDEX_PARAMETER] << "\t" << index_parameter_distribution->nb_value - 1 << endl;

      out_file << "\n" << SEQ_label[SEQL_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      index_parameter_distribution->spreadsheet_characteristic_print(out_file);

      out_file << "\n\t" << SEQ_label[SEQL_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      index_parameter_distribution->spreadsheet_print(out_file);
      out_file << endl;
    }

    out_file << nb_variable << "\t" << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;

    if (self_transition) {
      for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
        if (self_transition[i]) {
          out_file << "\n\t" << STAT_label[STATL_STATE] << " " << i << " - "
                   << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION] << endl;
          self_transition[i]->spreadsheet_print(out_file);

          smoothed_curves = new Curves(*(self_transition[i]) , SMOOTHING);

          out_file << "\n\t" << STAT_label[STATL_STATE] << " " << i << " - "
                   << SEQ_label[SEQL_SMOOTHED] << " " << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION] << endl;
          smoothed_curves->spreadsheet_print(out_file);

          delete smoothed_curves;
        }
      }
    }

    for (i = 0;i < nb_variable;i++) {
      out_file << "\n" << STAT_word[STATW_VARIABLE] << "\t" << i + 1 << "\t"
               << STAT_variable_word[type[i]];

      if (type[i] != AUXILIARY) {
        if (type[i] == STATE) {
          out_file << "\t\t" << marginal_distribution[i]->nb_value << "\t"
                   << STAT_label[marginal_distribution[i]->nb_value == 1 ? STATL_STATE : STATL_STATES] << endl;
        }
        else {
          out_file << "\t\t" << STAT_label[STATL_MIN_VALUE] << "\t" << min_value[i]
                   << "\t\t" << STAT_label[STATL_MAX_VALUE] << "\t" << max_value[i] << endl;
        }

        if (marginal_distribution[i]) {
          out_file << "\n" << STAT_label[type[i] == STATE ? STATL_STATE : STATL_MARGINAL] << " "
                   << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
          marginal_distribution[i]->spreadsheet_characteristic_print(out_file);

          out_file << "\n\t" << STAT_label[STATL_FREQUENCY] << endl;
          marginal_distribution[i]->spreadsheet_print(out_file);
        }

        else {
          mean = mean_computation(i);
          variance = variance_computation(i , mean);

          if (variance > 0.) {
            switch (type[i]) {

            case INT_VALUE : {
              int_value = new int[cumul_length];
              pint_value = int_value;
              for (j = 0;j < nb_sequence;j++) {
                for (k = 0;k < length[j];k++) {
                  *pint_value++ = int_sequence[j][i][k];
                }
              }

              lower_quartile = quantile_computation(cumul_length , int_value , 0.25);
              median = quantile_computation(cumul_length , int_value , 0.5);
              upper_quartile = quantile_computation(cumul_length , int_value , 0.75);

              delete [] int_value;
              break;
            }

            case REAL_VALUE : {
              real_value = new double[cumul_length];
              preal_value = real_value;
              for (j = 0;j < nb_sequence;j++) {
                for (k = 0;k < length[j];k++) {
                  *preal_value++ = real_sequence[j][i][k];
                }
              }

              lower_quartile = quantile_computation(cumul_length , real_value , 0.25);
              median = quantile_computation(cumul_length , real_value , 0.5);
              upper_quartile = quantile_computation(cumul_length , real_value , 0.75);

              delete [] real_value;
              break;
            }
            }
          }

          else {
            median = mean;
          }

          out_file << "\n" << STAT_label[STATL_SAMPLE_SIZE] << "\t" << cumul_length << endl;

          out_file << STAT_label[STATL_MEAN] << "\t" << mean << "\t\t"
                   << STAT_label[STATL_MEDIAN] << "\t" << median << endl;

          out_file << STAT_label[STATL_VARIANCE] << "\t" << variance << "\t\t"
                   << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(variance);
          if (variance > 0.) {
            out_file << "\t\t" << STAT_label[STATL_LOWER_QUARTILE] << "\t" << lower_quartile
                     << "\t\t" << STAT_label[STATL_UPPER_QUARTILE] << "\t" << upper_quartile;
          }
          out_file << endl;

          if (variance > 0.) {
            out_file << STAT_label[STATL_SKEWNESS_COEFF] << "\t" << skewness_computation(i , mean , variance) << "\t\t"
                     << STAT_label[STATL_KURTOSIS_COEFF] << "\t" << kurtosis_computation(i , mean , variance) << endl;
          }

          if (marginal_histogram[i]) {
            out_file << "\n" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM] << endl;
            out_file << "\n" << STAT_label[STATL_VALUE] << "\t" << STAT_label[STATL_FREQUENCY] << endl;
            marginal_histogram[i]->spreadsheet_print(out_file);
          }
        }

        if (characteristics[i]) {
          characteristics[i]->spreadsheet_print(out_file , type[i] , *length_distribution);
        }
      }

      else {
        out_file << endl;
      }
    }

    out_file << "\n" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    length_distribution->spreadsheet_characteristic_print(out_file);

    out_file << "\n\t" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    length_distribution->spreadsheet_print(out_file);

    out_file << "\n" << SEQ_label[SEQL_CUMUL_LENGTH] << "\t" << cumul_length << endl;
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a MarkovianSequences object using Gnuplot for a variable
 *         in the case of the absence of the characteristic distributions.
 *
 *  \param[in] prefix      file prefix,
 *  \param[in] title       figure title,
 *  \param[in] variable    variable index,
 *  \param[in] nb_variable number of variables.
 *
 *  \return                error status.
 */
/*--------------------------------------------------------------*/

bool MarkovianSequences::plot_print(const char *prefix , const char *title , int variable ,
                                    int nb_variable) const

{
  bool status;
  int i;
  int nb_histo;
  const FrequencyDistribution *phisto[1];
  ostringstream data_file_name[2];


  // writing of the data files

  data_file_name[0] << prefix << variable + 1 << 0 << ".dat";

  nb_histo = 0;
  if (index_parameter_distribution) {
    phisto[nb_histo++] = index_parameter_distribution;
  }

  status = length_distribution->plot_print((data_file_name[0].str()).c_str() , nb_histo , phisto);

  if (status) {
    if (marginal_distribution[variable]) {
      data_file_name[1] << prefix << variable + 1 << 1 << ".dat";
      marginal_distribution[variable]->plot_print((data_file_name[1].str()).c_str());
    }
    else if (marginal_histogram[variable]) {
      data_file_name[1] << prefix << variable + 1 << 1 << ".dat";
      marginal_histogram[variable]->plot_print((data_file_name[1].str()).c_str());
    }

    // writing of the script files

    for (i = 0;i < 2;i++) {
      ostringstream file_name[2];

      switch (i) {

      case 0 : {
        if (nb_variable == 1) {
          file_name[0] << prefix << ".plot";
        }
        else {
          file_name[0] << prefix << variable + 1 << ".plot";
        }
        break;
      }

      case 1 : {
        if (nb_variable == 1) {
          file_name[0] << prefix << ".print";
        }
        else {
          file_name[0] << prefix << variable + 1 << ".print";
        }
        break;
      }
      }

      ofstream out_file((file_name[0].str()).c_str());

      if (i == 1) {
        out_file << "set terminal postscript" << endl;

        if (nb_variable == 1) {
          file_name[1] << label(prefix) << ".ps";
        }
        else {
          file_name[1] << label(prefix) << variable + 1 << ".ps";
        }
        out_file << "set output \"" << file_name[1].str() << "\"\n\n";
      }

      out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
               << "set title";
      if (title) {
        out_file << " \"" << title << "\"";
      }
      out_file << "\n\n";

      if (marginal_distribution[variable]) {
        if (marginal_distribution[variable]->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }
        if ((int)(marginal_distribution[variable]->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [0:" << MAX(marginal_distribution[variable]->nb_value - 1 , 1) << "] [0:"
                 << (int)(marginal_distribution[variable]->max * YSCALE) + 1 << "] \""
                 << label((data_file_name[1].str()).c_str()) << "\" using 1 title \""
                 << STAT_label[STATL_VARIABLE] << " " << variable + 1 << " - "
                 << STAT_label[type[variable] == STATE ? STATL_STATE : STATL_MARGINAL] << " "
                 << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses" << endl;

        if (marginal_distribution[variable]->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
        if ((int)(marginal_distribution[variable]->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics autofreq" << endl;
        }

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;
      }

      else if (marginal_histogram[variable]) {
        if ((int)(marginal_histogram[variable]->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [" << marginal_histogram[variable]->min_value - marginal_histogram[variable]->bin_width << ":"
                 << marginal_histogram[variable]->max_value + marginal_histogram[variable]->bin_width << "] [0:"
                 << (int)(marginal_histogram[variable]->max * YSCALE) + 1 << "] \""
                 << label((data_file_name[1].str()).c_str()) << "\" using 1:2 title \""
                 << STAT_label[STATL_VARIABLE] << " " << variable + 1 << " "
                 << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM]
                 << "\" with histeps" << endl;

        if ((int)(marginal_histogram[variable]->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics autofreq" << endl;
        }

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;
      }

      if (length_distribution->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics 0,1" << endl;
      }
      if ((int)(length_distribution->max * YSCALE) + 1 < TIC_THRESHOLD) {
        out_file << "set ytics 0,1" << endl;
      }

      out_file << "plot [0:" << length_distribution->nb_value - 1 << "] [0:"
               << (int)(length_distribution->max * YSCALE) + 1 << "] \""
               << label((data_file_name[0].str()).c_str()) << "\" using 1 title \""
               << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
               << "\" with impulses" << endl;

      if (length_distribution->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
      }
      if ((int)(length_distribution->max * YSCALE) + 1 < TIC_THRESHOLD) {
        out_file << "set ytics autofreq" << endl;
      }

      if (index_parameter_distribution) {
        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        if (index_parameter_distribution->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }
        if ((int)(index_parameter_distribution->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [" << index_parameter_distribution->offset << ":"
                 << index_parameter_distribution->nb_value - 1 << "] [0:"
                 << (int)(index_parameter_distribution->max * YSCALE) + 1 << "] \""
                 << label((data_file_name[0].str()).c_str()) << "\" using 2 title \""
                 << SEQ_label[SEQL_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                 << "\" with impulses" << endl;

        if (index_parameter_distribution->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
        if ((int)(index_parameter_distribution->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics autofreq" << endl;
        }
      }

      if (i == 1) {
        out_file << "\nset terminal x11" << endl;
      }

      out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a MarkovianSequences object using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool MarkovianSequences::plot_write(StatError &error , const char *prefix ,
                                    const char *title) const

{
  bool status , start;
  int i , j;
  int max_frequency[NB_OUTPUT];
  ostringstream data_file_name[NB_OUTPUT];


  error.init();

  if (characteristics[0]) {
    status = characteristics[0]->plot_print(prefix , title , 0 , nb_variable , type[0] , *length_distribution);
  }
  else {
    status = plot_print(prefix , title , 0 , nb_variable);
  }

  if (status) {
    for (i = 1;i < nb_variable;i++) {
      if (characteristics[i]) {
        characteristics[i]->plot_print(prefix , title , i , nb_variable , type[i] , *length_distribution);
      }
      else {
        plot_print(prefix , title , i , nb_variable);
      }
    }

    if (self_transition) {
      for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
        if (self_transition[i]) {
          max_frequency[i] = self_transition[i]->max_frequency_computation();

          data_file_name[i] << prefix << i << ".dat";
          self_transition[i]->plot_print_standard_residual((data_file_name[i].str()).c_str());
        }
      }

      // writing of the script files

      for (i = 0;i < 2;i++) {
        ostringstream file_name[2];

        switch (i) {
        case 0 :
          file_name[0] << prefix << 1 << 0 << ".plot";
          break;
        case 1 :
          file_name[0] << prefix << 1 << 0 << ".print";
          break;
        }

        ofstream out_file((file_name[0].str()).c_str());

        if (i == 1) {
          out_file << "set terminal postscript" << endl;
          file_name[1] << label(prefix) << 1 << 0 << ".ps";
          out_file << "set output \"" << file_name[1].str() << "\"\n\n";
        }

        out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                 << "set title";
        if ((title) || (nb_variable > 1)) {
          out_file << " \"";
          if (title) {
            out_file << title;
            if (nb_variable > 1) {
              out_file << " - ";
            }
          }
          if (nb_variable > 1) {
            out_file << STAT_label[STATL_VARIABLE] << " " << 1;
          }
          out_file << "\"";
        }
        out_file << "\n\n";

        start = true;
        for (j = 0;j < marginal_distribution[0]->nb_value;j++) {
          if (self_transition[j]) {
            if (!start) {
              if (i == 0) {
                out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
              }
              out_file << endl;
            }
            else {
              start = false;
            }

            if (self_transition[j]->length - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }

            out_file << "plot [0:" << self_transition[j]->length - 1 << "] [0:1] \""
                     << label((data_file_name[j].str()).c_str()) << "\" using 1:2 title \""
                     << STAT_label[STATL_STATE] << " " << j << " - " << SEQ_label[SEQL_OBSERVED] << " "
                     << SEQ_label[SEQL_SELF_TRANSITION] << "\" with points" << endl;

            if (self_transition[j]->length - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set xlabel \"" << SEQ_label[SEQL_INDEX] << "\"" << endl;
            out_file << "set ylabel \"" << STAT_label[STATL_FREQUENCY] << "\"" << endl;

            if (self_transition[j]->length - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }
            if ((int)(max_frequency[j] * YSCALE) + 1 < TIC_THRESHOLD) {
              out_file << "set ytics 0,1" << endl;
            }

            out_file << "plot [0:" << self_transition[j]->length - 1
                     << "] [0:" << (int)(max_frequency[j] * YSCALE) + 1 << "] \""
                     << label((data_file_name[j].str()).c_str())
                     << "\" using 1:3 title \"" << SEQ_label[SEQL_TRANSITION_COUNTS]
                     << "\" with impulses" << endl;

            if (self_transition[j]->length - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
            if ((int)(max_frequency[j] * YSCALE) + 1 < TIC_THRESHOLD) {
              out_file << "set ytics autofreq" << endl;
            }

            out_file << "set xlabel" << endl;
            out_file << "set ylabel" << endl;
          }
        }

        if (i == 1) {
          out_file << "\nset terminal x11" << endl;
        }

        out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a MarkovianSequences object for a variable
 *         in the case of the absence of the characteristic distributions.
 *
 *  \param[in] plot     reference on a MultiPlotSet object,
 *  \param[in] index    MultiPlot index,
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::plotable_write(MultiPlotSet &plot , int &index , int variable) const

{
  ostringstream legend;


  /*  nb_plot_set = 1;
  if ((marginal_distribution[variable]) || (marginal_histogram[variable])) {
    nb_plot_set++;
  }
  if (index_parameter_distribution) {
    nb_plot_set++;
  } */

  plot.variable_nb_viewpoint[variable] = 1;

  if (marginal_distribution[variable]) {

    // marginal frequency distribution

    plot.variable[index] = variable;

    plot[index].xrange = Range(0 , MAX(marginal_distribution[variable]->nb_value - 1 , 1));
    plot[index].yrange = Range(0 , ceil(marginal_distribution[variable]->max * YSCALE));

    if (marginal_distribution[variable]->nb_value - 1 < TIC_THRESHOLD) {
      plot[index].xtics = 1;
    }
    if (ceil(marginal_distribution[variable]->max * YSCALE) < TIC_THRESHOLD) {
      plot[index].ytics = 1;
    }

    plot[index].resize(1);

    legend.str("");
    legend << STAT_label[STATL_VARIABLE] << " " << variable + 1 << " - "
           << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[index][0].legend = legend.str();

    plot[index][0].style = "impulses";

    marginal_distribution[variable]->plotable_frequency_write(plot[index][0]);
  }

  else if (marginal_histogram[variable]) {

    // marginal histogram

    plot[index].xrange = Range(marginal_histogram[variable]->min_value - marginal_histogram[variable]->bin_width ,
                               marginal_histogram[variable]->max_value + marginal_histogram[variable]->bin_width);
    plot[index].yrange = Range(0 , ceil(marginal_histogram[variable]->max * YSCALE));

    if (ceil(marginal_histogram[variable]->max * YSCALE) < TIC_THRESHOLD) {
      plot[index].ytics = 1;
    }

    plot[index].resize(1);

    legend.str("");
    legend << STAT_label[STATL_VARIABLE] << " " << variable + 1 << " "
           << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM];
    plot[index][0].legend = legend.str();

    plot[index][0].style = "histeps";

    marginal_histogram[variable]->plotable_write(plot[index][0]);
  }

  index++;

  // sequence length frequency distribution

  plot.variable[index] = variable;

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

  if (index_parameter_distribution) {

    // index parameter frequency distribution

    plot.variable[index] = variable;

    plot[index].xrange = Range(index_parameter_distribution->offset , index_parameter_distribution->nb_value - 1);
    plot[index].yrange = Range(0 , ceil(index_parameter_distribution->max * YSCALE));

    if (index_parameter_distribution->nb_value - 1 < TIC_THRESHOLD) {
      plot[index].xtics = 1;
    }
    if (ceil(index_parameter_distribution->max * YSCALE) < TIC_THRESHOLD) {
      plot[index].ytics = 1;
    }

    plot[index].resize(1);

    legend.str("");
    legend << SEQ_label[SEQL_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[index][0].legend = legend.str();

    plot[index][0].style = "impulses";

    index_parameter_distribution->plotable_frequency_write(plot[index][0]);
    index++;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a MarkovianSequences object.
 *
 *  \return MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* MarkovianSequences::get_plotable() const

{
  int i , j;
  int nb_plot_set , index_length , index , max_frequency;
  ostringstream title , legend;
  MultiPlotSet *plot_set;


  // computation of the number of plots

  nb_plot_set = 0;

  for (i = 0;i < nb_variable;i++) {
    if (characteristics[i]) {
      index_length = characteristics[i]->index_value->plot_length_computation();

      nb_plot_set += 2;
      if (characteristics[i]->index_value->frequency[index_length - 1] < MAX_FREQUENCY) {
        nb_plot_set++;
      }

      nb_plot_set++;
      for (j = 0;j < characteristics[i]->nb_value;j++) {
        if (characteristics[i]->first_occurrence[j]->nb_element > 0) {
          nb_plot_set++;
        }
      }

      nb_plot_set++;
      for (j = 0;j < characteristics[i]->nb_value;j++) {
        if (characteristics[i]->recurrence_time[j]->nb_element > 0) {
          nb_plot_set++;
        }
      }

      nb_plot_set++;
      for (j = 0;j < characteristics[i]->nb_value;j++) {
        if (characteristics[i]->sojourn_time[j]->nb_element > 0) {
          nb_plot_set++;
        }
        if ((characteristics[i]->initial_run) &&
            (characteristics[i]->initial_run[j]->nb_element > 0)) {
          nb_plot_set++;
        }
        if (characteristics[i]->final_run[j]->nb_element > 0) {
          nb_plot_set++;
        }
      }

      if ((characteristics[i]->nb_run) && (characteristics[i]->nb_occurrence)) {
        nb_plot_set += 3;
        for (j = 0;j < characteristics[i]->nb_value;j++) {
          if ((characteristics[i]->nb_run[j]->nb_element > 0) &&
              (characteristics[i]->nb_occurrence[j]->nb_element > 0)) {
            nb_plot_set += 2;
          }
        }
      }
    }

    else if (type[i] != AUXILIARY) {
      nb_plot_set++;
      if ((marginal_distribution[i]) || (marginal_histogram[i])) {
        nb_plot_set++;
      }
      if (index_parameter_distribution) {
        nb_plot_set++;
      }
    }
  }

  if (self_transition) {
    for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
      if (self_transition[i]) {
        nb_plot_set += 2;
      }
    }
  }

  plot_set = new MultiPlotSet(nb_plot_set , nb_variable);

  MultiPlotSet &plot = *plot_set;

  plot.border = "15 lw 0";

  for (i = 0;i < nb_variable;i++) {
    plot.variable_nb_viewpoint[i] = 0;
  }

  index = 0;
  if (self_transition) {
    plot.variable_nb_viewpoint[0]++;

    for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
      if (self_transition[i]) {

        // self-transition probability as a function of the index parameter

        if (nb_variable > 1) {
          title.str("");
          title << STAT_label[STATL_VARIABLE] << " " << 1;
          plot[index].title = title.str();
        }

        plot[index].xrange = Range(0 , self_transition[i]->length - 1);
        plot[index].yrange = Range(0. , 1.);

        if (self_transition[i]->length - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }

        plot[index].resize(1);

        legend.str("");
        legend << STAT_label[STATL_STATE] << " " << i << " - "
               << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION];
        plot[index][0].legend = legend.str();

        plot[index][0].style = "linespoint";

        self_transition[i]->plotable_write(0 , plot[index][0]);
        index++;

        // frequency distributions of indexed transition counts

        if (nb_variable > 1) {
          title.str("");
          title << STAT_label[STATL_VARIABLE] << " " << 1;
          plot[index].title = title.str();
        }

        plot[index].xrange = Range(0 , self_transition[i]->length - 1);
        max_frequency = self_transition[i]->max_frequency_computation();
        plot[index].yrange = Range(0 , ceil(max_frequency * YSCALE));

        if (self_transition[i]->length - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }
        if (ceil(max_frequency * YSCALE) < TIC_THRESHOLD) {
          plot[index].ytics = 1;
        }

        plot[index].xlabel = SEQ_label[SEQL_INDEX];
        plot[index].ylabel = STAT_label[STATL_FREQUENCY];

        plot[index].resize(1);

        legend.str("");
        legend << STAT_label[STATL_STATE] << " " << i << " - "
               << SEQ_label[SEQL_TRANSITION_COUNTS];
        plot[index][0].legend = legend.str();

        plot[index][0].style = "impulses";

        self_transition[i]->plotable_frequency_write(plot[index][0]);
        index++;
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    if (characteristics[i]) {
      characteristics[i]->plotable_write(plot , index , i , type[i] , *length_distribution);
    }
    else if (type[i] != AUXILIARY) {
      plotable_write(plot , index , i);
    }
  }

  return plot_set;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a fitted observation linear trend model at the spreadsheet format.
 *
 *  \param[in,out] os       stream,
 *  \param[in]     variable variable index,
 *  \param[in]     process  pointer on a continuous observation process.
 */
/*--------------------------------------------------------------*/

ostream& MarkovianSequences::linear_model_spreadsheet_print(ostream &os , int variable ,
                                                            ContinuousParametricProcess *process) const

{
  bool *used_sequence;
  int i , j , k , m , n , r;
  int frequency , *index;
  double buff;


  if (type[variable] == INT_VALUE) {
    used_sequence = new bool[nb_sequence];
  }
  if (index_param_type == TIME) {
    index = new int[nb_sequence];
  }

  os << "\n" << SEQ_label[SEQL_INDEX] << "\t" << STAT_label[STATL_OBSERVATION];
  for (i = 0;i < process->nb_state;i++) {
    os << "\t" << STAT_label[STATL_STATE] << " " << i << " " << STAT_label[STATL_LINEAR_MODEL];
  }
  if (type[variable] == INT_VALUE) {
    os << "\t" << STAT_label[STATL_FREQUENCY];
  }
  os << endl;

  switch (type[0]) {

  case STATE : {
    for (i = 0;i < process->nb_state;i++) {
      switch (index_param_type) {

      case IMPLICIT_TYPE : {
        switch (type[variable]) {

        case INT_VALUE : {
          for (j = 0;j < max_length;j++) {
            for (k = 0;k < nb_sequence;k++) {
              used_sequence[k] = false;
            }
            for (k = 0;k < nb_sequence;k++) {
              if ((j < length[k]) && (int_sequence[k][0][j] == i) && (!used_sequence[k])) {
                used_sequence[k] = true;
                os << j << "\t" << int_sequence[k][variable][j];
                for (m = 0;m <= i;m++) {
                  os << "\t";
                }
                os << process->observation[i]->intercept + process->observation[i]->slope * j;

                frequency = 1;
                for (m = k + 1;m < nb_sequence;m++) {
                  if ((j < length[m]) && (int_sequence[m][0][j] == i) &&
                      (int_sequence[m][variable][j] == int_sequence[k][variable][j])) {
                    used_sequence[m] = true;
                    frequency++;
                  }
                }
                for (m = i;m < process->nb_state;m++) {
                  os << "\t";
                }
                os << frequency << endl;
              }
            }
          }
          break;
        }

        case REAL_VALUE : {
          for (j = 0;j < max_length;j++) {
            for (k = 0;k < nb_sequence;k++) {
              if ((j < length[k]) && (int_sequence[k][0][j] == i)) {
                os << j << "\t" << real_sequence[k][variable][j];
                for (m = 0;m <= i;m++) {
                  os << "\t";
                }
                os << process->observation[i]->intercept + process->observation[i]->slope * j << endl;
              }
            }
          }
          break;
        }
        }
        break;
      }

      case TIME : {
        for (j = 0;j < nb_sequence;j++) {
          index[j] = 0;
        }

        switch (type[variable]) {

        case INT_VALUE : {
          for (j = index_parameter_distribution->offset;j < index_parameter_distribution->nb_value;j++) {
            if (index_parameter_distribution->frequency[j] > 0) {
              for (k = 0;k < nb_sequence;k++) {
                used_sequence[k] = false;
              }
              for (k = 0;k < nb_sequence;k++) {
                m = index[k];
                while ((m < length[k] - 1) && (index_parameter[k][m] < j)) {
                  m++;
                }

                if ((index_parameter[k][m] == j) && (int_sequence[k][0][m] == i) && (!used_sequence[k])) {
                  index[k] = m;
                  used_sequence[k] = true;
                  os << j << "\t" << int_sequence[k][variable][m];
                  for (n = 0;n <= i;n++) {
                    os << "\t";
                  }
                  os << process->observation[i]->intercept + process->observation[i]->slope * j;

                  frequency = 1;
                  for (n = k + 1;n < nb_sequence;n++) {
                    r = index[n];
                    while ((r < length[n] - 1) && (index_parameter[n][r] < j)) {
                      r++;
                    }

                    if ((index_parameter[n][r] == j) && (int_sequence[n][0][r] == i) &&
                        (int_sequence[n][variable][r] == int_sequence[k][variable][m])) {
                      index[n] = r;
                      used_sequence[n] = true;
                      frequency++;
                    }
                  }

                  for (n = i;n < process->nb_state;n++) {
                    os << "\t";
                  }
                  os << frequency << endl;
                }
              }
            }
          }
          break;
        }

        case REAL_VALUE : {
          for (j = index_parameter_distribution->offset;j < index_parameter_distribution->nb_value;j++) {
            if (index_parameter_distribution->frequency[j] > 0) {
              for (k = 0;k < nb_sequence;k++) {
                m = index[k];
                while ((m < length[k] - 1) && (index_parameter[k][m] < j)) {
                  m++;
                }

                if ((index_parameter[k][m] == j) && (int_sequence[k][0][m] == i)) {
                  index[k] = m;
                  os << j << "\t" << real_sequence[k][variable][m];
                  for (n = 0;n <= i;n++) {
                    os << "\t";
                  }
                  os << process->observation[i]->intercept + process->observation[i]->slope * j << endl;
                }
              }
            }
          }
          break;
        }
        }
        break;
      }
      }
    }
    break;
  }

  default : {
    switch (index_param_type) {

    case IMPLICIT_TYPE : {
      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < max_length;i++) {
          for (j = 0;j < nb_sequence;j++) {
            used_sequence[j] = false;
          }
          for (j = 0;j < nb_sequence;j++) {
            if ((i < length[j]) && (!used_sequence[j])) {
              used_sequence[j] = true;
              os << i << "\t" << int_sequence[j][variable][i];
              for (k = 0;k < process->nb_state;k++) {
                os << "\t";
                buff = process->observation[k]->intercept + process->observation[k]->slope * i;
                if ((buff >= min_value[variable]) && (buff <= max_value[variable])) {
                  os << buff;
                }
              }

              frequency = 1;
              for (k = j + 1;k < nb_sequence;k++) {
                if ((i < length[k]) && (int_sequence[k][variable][i] == int_sequence[j][variable][i])) {
                  used_sequence[k] = true;
                  frequency++;
                }
              }
              os << "\t" << frequency << endl;
            }
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < max_length;i++) {
          for (j = 0;j < nb_sequence;j++) {
            if (i < length[j]) {
              os << i << "\t" << real_sequence[j][variable][i];
              for (k = 0;k < process->nb_state;k++) {
                os << "\t";
                buff = process->observation[k]->intercept + process->observation[k]->slope * i;
                if ((buff >= min_value[variable]) && ( buff <= max_value[variable])) {
                  os << buff;
                }
              }
              os << endl;
            }
          }
        }
        break;
      }
      }
      break;
    }

    case TIME : {
      for (i = 0;i < nb_sequence;i++) {
        index[i] = 0;
      }

      switch (type[variable]) {

      case INT_VALUE : {
        for (i = index_parameter_distribution->offset;i < index_parameter_distribution->nb_value;i++) {
          if (index_parameter_distribution->frequency[i] > 0) {
            for (j = 0;j < nb_sequence;j++) {
              used_sequence[j] = false;
            }
            for (j = 0;j < nb_sequence;j++) {
              k = index[j];
              while ((k < length[j] - 1) && (index_parameter[j][k] < i)) {
                k++;
              }

              if ((index_parameter[j][k] == i) && (!used_sequence[j])) {
                index[j] = k;
                used_sequence[j] = true;
                os << i << "\t" << int_sequence[j][variable][k];
                for (m = 0;m < process->nb_state;m++) {
                  os << "\t";
                  buff = process->observation[m]->intercept + process->observation[m]->slope * i;
                  if ((buff >= min_value[variable]) && (buff <= max_value[variable])) {
                    os << buff;
                  }
                }

                frequency = 1;
                for (m = j + 1;m < nb_sequence;m++) {
                  n = index[m];
                  while ((n < length[m] - 1) && (index_parameter[m][n] < i)) {
                    n++;
                  }

                  if ((index_parameter[m][n] == i) && (int_sequence[m][variable][n] == int_sequence[j][variable][k])) {
                    index[m] = n;
                    used_sequence[m] = true;
                    frequency++;
                  }
                }
                os << "\t" << frequency << endl;
              }
            }
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = index_parameter_distribution->offset;i < index_parameter_distribution->nb_value;i++) {
          if (index_parameter_distribution->frequency[i] > 0) {
            for (j = 0;j < nb_sequence;j++) {
              k = index[j];
              while ((k < length[j] - 1) && (index_parameter[j][k] < i)) {
                k++;
              }

              if (index_parameter[j][k] == i) {
                index[j] = k;
                os << i << "\t" << real_sequence[j][variable][k];
                for (m = 0;m < process->nb_state;m++) {
                  os << "\t";
                  buff = process->observation[m]->intercept + process->observation[m]->slope * i;
                  if ((buff >= min_value[variable]) && (buff <= max_value[variable])) {
                    os << buff;
                  }
                }
                os << endl;
              }
            }
          }
        }
        break;
      }
      }
      break;
    }
    }
    break;
  }
  }

  if (type[variable] == INT_VALUE) {
    delete [] used_sequence;
  }
  if (index_param_type == TIME) {
    delete [] index;
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a fitted observation linear trend model using Gnuplot.
 *
 *  \param[in] prefix   file prefix,
 *  \param[in] title    figure title,
 *  \param[in] variable variable index,
 *  \param[in] process  pointer on a continuous observation process.
 *
 *  \return             error status.
 */
/*--------------------------------------------------------------*/

bool MarkovianSequences::linear_model_plot_print(const char *prefix , const char *title , int variable ,
                                                 ContinuousParametricProcess *process) const

{
  bool status = false;
  int i , j;
  int process_index , *state_min_index_parameter , *state_max_index_parameter , *pstate;
  double buff , *state_min_value , *state_max_value;
  ostringstream data_file_name[NB_STATE * 2 + 1];
  ofstream *out_data_file[NB_STATE + 1];


  // writing of data files

  state_min_index_parameter = new int[process->nb_state + 1];
  state_max_index_parameter = new int[process->nb_state + 1];

  switch (type[0]) {

  case STATE : {
    process_index = variable;

    switch (index_param_type) {

    case IMPLICIT_TYPE : {
      for (i = 0;i < process->nb_state;i++) {
        state_min_index_parameter[i] = max_length - 1;
        state_max_index_parameter[i] = 0;
      }

      for (i = 0;i < nb_sequence;i++) {
        pstate = int_sequence[i][0];

        for (j = 0;j < length[i];j++) {
          if (j < state_min_index_parameter[*pstate]) {
            state_min_index_parameter[*pstate] = j;
          }
          if (j > state_max_index_parameter[*pstate]) {
            state_max_index_parameter[*pstate] = j;
          }
          pstate++;
        }
      }

      state_min_index_parameter[process->nb_state] = 0;
      state_max_index_parameter[process->nb_state] = max_length - 1;
      break;
    }

    case TIME : {
      for (i = 0;i < process->nb_state;i++) {
        state_min_index_parameter[i] = index_parameter_distribution->nb_value - 1;
        state_max_index_parameter[i] = index_parameter_distribution->offset;
      }

      for (i = 0;i < nb_sequence;i++) {
        pstate = int_sequence[i][0];

        for (j = 0;j < length[i];j++) {
          if (index_parameter[i][j] < state_min_index_parameter[*pstate]) {
            state_min_index_parameter[*pstate] = index_parameter[i][j];
          }
          if (index_parameter[i][j] > state_max_index_parameter[*pstate]) {
            state_max_index_parameter[*pstate] = index_parameter[i][j];
          }
          pstate++;
        }
      }

      state_min_index_parameter[process->nb_state] = index_parameter_distribution->offset;
      state_max_index_parameter[process->nb_state] = index_parameter_distribution->nb_value - 1;
      break;
    }
    }

    for (i = 0;i < process->nb_state;i++) {
      if (state_max_index_parameter[i] == state_min_index_parameter[i]) {
        state_max_index_parameter[i]++;
      }

      if (marginal_distribution[0]->frequency[i] == 0) {
        switch (index_param_type) {

        case IMPLICIT_TYPE : {
          state_min_index_parameter[i] = 0;
          state_max_index_parameter[i] = max_length - 1;
          break;
        }

        case TIME : {
          state_min_index_parameter[i] = index_parameter_distribution->offset;
          state_max_index_parameter[i] = index_parameter_distribution->nb_value - 1;
          break;
        }
        }

        buff = (min_value[variable] - process->observation[i]->intercept) / process->observation[i]->slope;
        if ((process->observation[i]->slope > 0.) && (buff > state_min_index_parameter[i])) {
          state_min_index_parameter[i] = ceil(buff);
        }
        if ((process->observation[i]->slope < 0.) && (buff < state_max_index_parameter[i])) {
          state_max_index_parameter[i] = floor(buff);
        }

        buff = (max_value[variable] - process->observation[i]->intercept) / process->observation[i]->slope;
        if ((process->observation[i]->slope < 0.) && (buff > state_min_index_parameter[i])) {
          state_min_index_parameter[i] = ceil(buff);
        }
        if ((process->observation[i]->slope > 0.) && (buff < state_max_index_parameter[i])) {
          state_max_index_parameter[i] = floor(buff);
        }
      }
    }

    state_min_value = new double[process->nb_state];
    state_max_value = new double[process->nb_state];
 
    for (i = 0;i < process->nb_state;i++) {
      state_min_value[i] = max_value[variable];
      state_max_value[i] = min_value[variable];
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        pstate = int_sequence[i][0];

        for (j = 0;j < length[i];j++) {
          if (int_sequence[i][variable][j] < state_min_value[*pstate]) {
            state_min_value[*pstate] = int_sequence[i][variable][j];
          }
          if (int_sequence[i][variable][j] > state_max_value[*pstate]) {
            state_max_value[*pstate] = int_sequence[i][variable][j];
          }
          pstate++;
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        pstate = int_sequence[i][0];

        for (j = 0;j < length[i];j++) {
          if (real_sequence[i][variable][j] < state_min_value[*pstate]) {
            state_min_value[*pstate] = real_sequence[i][variable][j];
          }
          if (real_sequence[i][variable][j] > state_max_value[*pstate]) {
            state_max_value[*pstate] = real_sequence[i][variable][j];
          }
          pstate++;
        }
      }
      break;
    }
    }
    break;
  }

  default : {
    process_index = variable + 1;

    switch (index_param_type) {

    case IMPLICIT_TYPE : {
      for (i = 0;i <= process->nb_state;i++) {
        state_min_index_parameter[i] = 0;
        state_max_index_parameter[i] = max_length - 1;
      }
      break;
    }

    case TIME : {
      for (i = 0;i <= process->nb_state;i++) {
        state_min_index_parameter[i] = index_parameter_distribution->offset;
        state_max_index_parameter[i] = index_parameter_distribution->nb_value - 1;
      }
      break;
    }
    }

    for (i = 0;i < process->nb_state;i++) {
      buff = (min_value[variable] - process->observation[i]->intercept) / process->observation[i]->slope;
      if ((process->observation[i]->slope > 0.) && (buff > state_min_index_parameter[i])) {
        state_min_index_parameter[i] = ceil(buff);
      }
      if ((process->observation[i]->slope < 0.) && (buff < state_max_index_parameter[i])) {
        state_max_index_parameter[i] = floor(buff);
      }

      buff = (max_value[variable] - process->observation[i]->intercept) / process->observation[i]->slope;
      if ((process->observation[i]->slope < 0.) && (buff > state_min_index_parameter[i])) {
        state_min_index_parameter[i] = ceil(buff);
      }
      if ((process->observation[i]->slope > 0.) && (buff < state_max_index_parameter[i])) {
        state_max_index_parameter[i] = floor(buff);
      }

//      cout << "\n" << STAT_label[STATL_STATE] << " " << i << ": " << state_min_index_parameter[i]
//           << ", " << state_max_index_parameter[i] << endl;
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    data_file_name[i * 2] << prefix << process_index << i * 2 << ".dat";
    out_data_file[0]= new ofstream((data_file_name[i * 2].str()).c_str());

    if (out_data_file[0]) {
      status = true;

      *out_data_file[0] << state_min_index_parameter[i] << " "
                        << process->observation[i]->intercept + process->observation[i]->slope * state_min_index_parameter[i] << endl;
      *out_data_file[0] << state_max_index_parameter[i] << " "
                        << process->observation[i]->intercept + process->observation[i]->slope * state_max_index_parameter[i] << endl;

      out_data_file[0]->close();
      delete out_data_file[0];
    }
  }

  if (type[0] == STATE) {
    for (i = 0;i < process->nb_state;i++) {
      if (marginal_distribution[0]->frequency[i] > 0) {
        data_file_name[i * 2 + 1] << prefix << process_index << i * 2 + 1 << ".dat";
        out_data_file[i] = new ofstream ((data_file_name[i * 2 + 1].str()).c_str());
      }
    }
    switch (index_param_type) {

    case IMPLICIT_TYPE : {
      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          pstate = int_sequence[i][0];
          for (j = 0;j < length[i];j++) {
            *out_data_file[*pstate++] << j << " " << int_sequence[i][variable][j] << endl;
           }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          pstate = int_sequence[i][0];
          for (j = 0;j < length[i];j++) {
            *out_data_file[*pstate++] << j << " " << real_sequence[i][variable][j] << endl;
          }
        }
        break;
      }
      }
      break;
    }

    case TIME : {
      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          pstate = int_sequence[i][0];
          for (j = 0;j < length[i];j++) {
            *out_data_file[*pstate++] << index_parameter[i][j] << " " << int_sequence[i][variable][j] << endl;
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          pstate = int_sequence[i][0];
          for (j = 0;j < length[i];j++) {
            *out_data_file[*pstate++] << index_parameter[i][j] << " " << real_sequence[i][variable][j] << endl;
          }
        }
        break;
      }
      }
      break;
    }
    }

    for (i = 0;i < process->nb_state;i++) {
      if (marginal_distribution[0]->frequency[i] > 0) {
        out_data_file[i]->close();
        delete out_data_file[i];
      }
    }
  }

  data_file_name[process->nb_state * 2] << prefix << process_index << process->nb_state * 2 << ".dat";
  out_data_file[process->nb_state] = new ofstream ((data_file_name[process->nb_state * 2].str()).c_str());

  if (out_data_file[process->nb_state]) {
    switch (index_param_type) {

    case IMPLICIT_TYPE : {
      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            *out_data_file[process->nb_state] << j << " " << int_sequence[i][variable][j] << endl;
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            *out_data_file[process->nb_state] << j << " " << real_sequence[i][variable][j] << endl;
          }
        }
        break;
      }
      }
      break;
    }

    case TIME : {
      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            *out_data_file[process->nb_state] << index_parameter[i][j] << " " << int_sequence[i][variable][j] << endl;
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            *out_data_file[process->nb_state] << index_parameter[i][j] << " " << real_sequence[i][variable][j] << endl;
          }
        }
        break;
      }
      }
      break;
    }
    }

    out_data_file[process->nb_state]->close();
    delete out_data_file[process->nb_state];
  }

  if (status) {

    // writing of the script files

    for (i = 0;i < 2;i++) {
      ostringstream file_name[2];

      switch (i) {
      case 0 :
        file_name[0] << prefix << process_index << ".plot";
        break;
      case 1 :
        file_name[0] << prefix << process_index << ".print";
        break;
      }

      ofstream out_file((file_name[0].str()).c_str());

      if (i == 1) {
        out_file << "set terminal postscript" << endl;
        file_name[1] << label(prefix) << process << ".ps";
        out_file << "set output \"" << file_name[1].str() << "\"\n\n";
      }

      out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
               << "set title \"";
      if (title) {
        out_file << title << " - ";
      }
      out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process_index << "\"\n\n";

      out_file << "set xlabel \"" << SEQ_label[SEQL_INDEX] << "\"" << endl;
      out_file << "set ylabel \"" << STAT_label[STATL_OBSERVATION] << "\"" << endl;

      if (type[0] == STATE) {
        for (j = 0;j < process->nb_state;j++) {
          if (marginal_distribution[0]->frequency[j] > 0) {
            if (state_max_index_parameter[j] - state_min_index_parameter[j] < TIC_THRESHOLD) {
              out_file << "set xtics " << state_min_index_parameter[j] << ",1" << endl;
            }
            if (state_max_value[j] - state_min_value[j] < TIC_THRESHOLD) {
              out_file << "set ytics " << MIN(state_min_value[j] , 0) << ",1" << endl;
            }

            out_file << "plot [" << state_min_index_parameter[j] << ":" << state_max_index_parameter[j] << "] [";
/*            if ((state_min_value[j] >= 0.) && (state_max_value[j] - state_min_value[j] > state_min_value[j] * PLOT_RANGE_RATIO)) {
              out_file << 0;
            }
            else { */
              out_file << state_min_value[j];
//            }
            out_file << ":" << MAX(state_max_value[j] , state_min_value[j] + 1) << "] \""
                     << label((data_file_name[2 * j + 1].str()).c_str()) << "\" using 1:2 notitle with points,\\" << endl;
            out_file << "\"" << label((data_file_name[2 * j].str()).c_str()) << "\" using 1:2 title \""
                     << STAT_label[STATL_STATE] << " " << j << " " << STAT_label[STATL_OBSERVATION] << " "
                     << STAT_label[STATL_MODEL] << "\" with lines" << endl;

            if (state_max_index_parameter[j] - state_min_index_parameter[j] < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
            if (state_max_value[j] - state_min_value[j] < TIC_THRESHOLD) {
              out_file << "set ytics autofreq" << endl;
            }

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;
          }
        }
      }

      if (state_max_index_parameter[process->nb_state] - state_min_index_parameter[process->nb_state] < TIC_THRESHOLD) {
        out_file << "set xtics " << state_min_index_parameter[process->nb_state] << ",1" << endl;
      }
      if (max_value[variable] - min_value[variable] < TIC_THRESHOLD) {
        out_file << "set ytics " << MIN(min_value[variable] , 0) << ",1" << endl;
      }

      out_file << "plot [" << state_min_index_parameter[process->nb_state] << ":" << state_max_index_parameter[process->nb_state] << "] [";
/*      if ((min_value[variable] >= 0.) && (max_value[variable] - min_value[variable] > min_value[variable] * PLOT_RANGE_RATIO)) {
        out_file << 0;
      }
      else { */
        out_file << min_value[variable];
//      }
      out_file << ":" << MAX(max_value[variable] , min_value[variable] + 1) << "] \""
               << label((data_file_name[2 * process->nb_state].str()).c_str()) << "\" using 1:2 notitle with points,\\" << endl;
      for (j = 0;j < process->nb_state;j++) {
        out_file << "\"" << label((data_file_name[2 * j].str()).c_str()) << "\" using 1:2 title \""
                 << STAT_label[STATL_STATE] << " " << j << " " << STAT_label[STATL_OBSERVATION] << " "
                 << STAT_label[STATL_MODEL] << "\" with lines";
        if (j < process->nb_state - 1) {
          out_file << ",\\";
        }
        out_file << endl;
      }

      if (state_max_index_parameter[process->nb_state] - state_min_index_parameter[process->nb_state] < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
      }
      if (max_value[variable] - min_value[variable] < TIC_THRESHOLD) {
        out_file << "set ytics autofreq" << endl;
      }

      out_file << "set xlabel" << endl;
      out_file << "set ylabel" << endl;

      if (i == 1) {
        out_file << "\nset terminal x11" << endl;
      }

      out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
    }
  }

  delete [] state_min_index_parameter;
  delete [] state_max_index_parameter;
  if (type[0] == STATE) {
    delete [] state_min_value;
    delete [] state_max_value;
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a fitted observation linear trend model.
 *
 *  \param[in] plot     file prefix,
 *  \param[in] index    MultiPlot index,
 *  \param[in] variable variable index,
 *  \param[in] process  pointer on a continuous observation process.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::linear_model_plotable_write(MultiPlotSet &plot , int &index , int variable ,
                                                     ContinuousParametricProcess *process) const

{
  bool status = false;
  int i , j;
  int process_index , plot_offset , *state_min_index_parameter , *state_max_index_parameter , *pstate;
  double buff , *state_min_value , *state_max_value;
  ostringstream title , legend;


  // computation of bounds

  state_min_index_parameter = new int[process->nb_state + 1];
  state_max_index_parameter = new int[process->nb_state + 1];

  switch (type[0]) {

  case STATE : {
    process_index = variable;
    plot_offset = process->nb_state;

    switch (index_param_type) {

    case IMPLICIT_TYPE : {
      for (i = 0;i < process->nb_state;i++) {
        state_min_index_parameter[i] = max_length - 1;
        state_max_index_parameter[i] = 0;
      }

      for (i = 0;i < nb_sequence;i++) {
        pstate = int_sequence[i][0];

        for (j = 0;j < length[i];j++) {
          if (j < state_min_index_parameter[*pstate]) {
            state_min_index_parameter[*pstate] = j;
          }
          if (j > state_max_index_parameter[*pstate]) {
            state_max_index_parameter[*pstate] = j;
          }
          pstate++;
        }
      }

      state_min_index_parameter[process->nb_state] = 0;
      state_max_index_parameter[process->nb_state] = max_length - 1;
      break;
    }

    case TIME : {
      for (i = 0;i < process->nb_state;i++) {
        state_min_index_parameter[i] = index_parameter_distribution->nb_value - 1;
        state_max_index_parameter[i] = index_parameter_distribution->offset;
      }

      for (i = 0;i < nb_sequence;i++) {
        pstate = int_sequence[i][0];

        for (j = 0;j < length[i];j++) {
          if (index_parameter[i][j] < state_min_index_parameter[*pstate]) {
            state_min_index_parameter[*pstate] = index_parameter[i][j];
          }
          if (index_parameter[i][j] > state_max_index_parameter[*pstate]) {
            state_max_index_parameter[*pstate] = index_parameter[i][j];
          }
          pstate++;
        }
      }

      state_min_index_parameter[process->nb_state] = index_parameter_distribution->offset;
      state_max_index_parameter[process->nb_state] = index_parameter_distribution->nb_value - 1;
      break;
    }
    }

    for (i = 0;i < process->nb_state;i++) {
      if (state_max_index_parameter[i] == state_min_index_parameter[i]) {
        state_max_index_parameter[i]++;
      }

      if (marginal_distribution[0]->frequency[i] == 0) {
        switch (index_param_type) {

        case IMPLICIT_TYPE : {
          state_min_index_parameter[i] = 0;
          state_max_index_parameter[i] = max_length - 1;
          break;
        }

        case TIME : {
          state_min_index_parameter[i] = index_parameter_distribution->offset;
          state_max_index_parameter[i] = index_parameter_distribution->nb_value - 1;
          break;
        }
        }

        buff = (min_value[variable] - process->observation[i]->intercept) / process->observation[i]->slope;
        if ((process->observation[i]->slope > 0.) && (buff > state_min_index_parameter[i])) {
          state_min_index_parameter[i] = ceil(buff);
        }
        if ((process->observation[i]->slope < 0.) && (buff < state_max_index_parameter[i])) {
          state_max_index_parameter[i] = floor(buff);
        }

        buff = (max_value[variable] - process->observation[i]->intercept) / process->observation[i]->slope;
        if ((process->observation[i]->slope < 0.) && (buff > state_min_index_parameter[i])) {
          state_min_index_parameter[i] = ceil(buff);
        }
        if ((process->observation[i]->slope > 0.) && (buff < state_max_index_parameter[i])) {
          state_max_index_parameter[i] = floor(buff);
        }
      }
    }

    state_min_value = new double[process->nb_state];
    state_max_value = new double[process->nb_state];

    for (i = 0;i < process->nb_state;i++) {
      state_min_value[i] = max_value[variable];
      state_max_value[i] = min_value[variable];
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        pstate = int_sequence[i][0];

        for (j = 0;j < length[i];j++) {
          if (int_sequence[i][variable][j] < state_min_value[*pstate]) {
            state_min_value[*pstate] = int_sequence[i][variable][j];
          }
          if (int_sequence[i][variable][j] > state_max_value[*pstate]) {
            state_max_value[*pstate] = int_sequence[i][variable][j];
          }
          pstate++;
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        pstate = int_sequence[i][0];

        for (j = 0;j < length[i];j++) {
          if (real_sequence[i][variable][j] < state_min_value[*pstate]) {
            state_min_value[*pstate] = real_sequence[i][variable][j];
          }
          if (real_sequence[i][variable][j] > state_max_value[*pstate]) {
            state_max_value[*pstate] = real_sequence[i][variable][j];
          }
          pstate++;
        }
      }
      break;
    }
    }
    break;
  }

  default : {
    process_index = variable + 1;
    plot_offset = 0;

    switch (index_param_type) {

    case IMPLICIT_TYPE : {
      for (i = 0;i <= process->nb_state;i++) {
        state_min_index_parameter[i] = 0;
        state_max_index_parameter[i] = max_length - 1;
      }
      break;
    }

    case TIME : {
      for (i = 0;i <= process->nb_state;i++) {
        state_min_index_parameter[i] = index_parameter_distribution->offset;
        state_max_index_parameter[i] = index_parameter_distribution->nb_value - 1;
      }
      break;
    }
    }

    for (i = 0;i < process->nb_state;i++) {
      buff = (min_value[variable] - process->observation[i]->intercept) / process->observation[i]->slope;
      if ((process->observation[i]->slope > 0.) && (buff > state_min_index_parameter[i])) {
        state_min_index_parameter[i] = ceil(buff);
      }
      if ((process->observation[i]->slope < 0.) && (buff < state_max_index_parameter[i])) {
        state_max_index_parameter[i] = floor(buff);
      }

      buff = (max_value[variable] - process->observation[i]->intercept) / process->observation[i]->slope;
      if ((process->observation[i]->slope < 0.) && (buff > state_min_index_parameter[i])) {
        state_min_index_parameter[i] = ceil(buff);
      }
      if ((process->observation[i]->slope > 0.) && (buff < state_max_index_parameter[i])) {
        state_max_index_parameter[i] = floor(buff);
      }
    }
    break;
  }
  }

  plot.variable_nb_viewpoint[variable] = 1;

  title.str("");
  title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process_index;

  // linear function and observations for each state

  if (type[0] == STATE) {
    for (i = 0;i < process->nb_state;i++) {
      plot[index + i].title = title.str();

      plot[index + i].xrange = Range(state_min_index_parameter[i] , state_max_index_parameter[i]);
      plot[index + i].yrange = Range(state_min_value[i] , MAX(state_max_value[i] , state_min_value[i] + 1));

      if (state_max_index_parameter[i] - state_min_index_parameter[i] < TIC_THRESHOLD) {
        plot[index + i].xtics = 1;
      }
      if (state_max_value[i] - state_min_value[i] < TIC_THRESHOLD) {
        plot[index + i].ytics = 1;
      }

      plot[index + i].xlabel = SEQ_label[SEQL_INDEX];
      plot[index + i].ylabel = STAT_label[STATL_OBSERVATION];

      plot[index + i].resize(2);

/*      legend.str("");
      legend << STAT_label[STATL_STATE] << " " << i << " " << STAT_label[STATL_OBSERVATION];
      plot[index + i][0].legend = legend.str(); */

      legend.str("");
      legend << STAT_label[STATL_STATE] << " " << i << " " << STAT_label[STATL_OBSERVATION] << " "
             << STAT_label[STATL_MODEL];
      plot[index + i][1].legend = legend.str();

      plot[index + i][0].style = "points";
      plot[index + i][1].style = "lines";
    }

    switch (index_param_type) {

    case IMPLICIT_TYPE : {
      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          pstate = int_sequence[i][0];
          for (j = 0;j < length[i];j++) {
            plot[index + *pstate++][0].add_point(j , int_sequence[i][variable][j]);
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          pstate = int_sequence[i][0];
          for (j = 0;j < length[i];j++) {
            plot[index + *pstate++][0].add_point(j , real_sequence[i][variable][j]);
          }
        }
        break;
      }
      }
      break;
    }

    case TIME : {
      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          pstate = int_sequence[i][0];
          for (j = 0;j < length[i];j++) {
            plot[index + *pstate++][0].add_point(index_parameter[i][j] , int_sequence[i][variable][j]);
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          pstate = int_sequence[i][0];
          for (j = 0;j < length[i];j++) {
            plot[index + *pstate++][0].add_point(index_parameter[i][j] , real_sequence[i][variable][j]);
          }
        }
        break;
      }
      }
      break;
    }
    }

    for (i = 0;i < process->nb_state;i++) {
      plot[index + i][1].add_point(state_min_index_parameter[i] , process->observation[i]->intercept +
                                   process->observation[i]->slope * state_min_index_parameter[i]);
      plot[index + i][1].add_point(state_max_index_parameter[i] , process->observation[i]->intercept +
                                   process->observation[i]->slope * state_max_index_parameter[i]);
    }
  }

  // linear functions and pooled observations

  plot[index + plot_offset].title = title.str();

  plot[index + plot_offset].xrange = Range(state_min_index_parameter[process->nb_state] ,
                                           state_max_index_parameter[process->nb_state]);
  plot[index + plot_offset].yrange = Range(min_value[variable] , MAX(max_value[variable] , min_value[variable]));

  if (state_max_index_parameter[process->nb_state] - state_min_index_parameter[process->nb_state] < TIC_THRESHOLD) {
    plot[index + plot_offset].xtics = 1;
  }
  if (max_value[variable] - min_value[variable] < TIC_THRESHOLD) {
    plot[index + plot_offset].ytics = 1;
  }

  plot[index + plot_offset].xlabel = SEQ_label[SEQL_INDEX];
  plot[index + plot_offset].ylabel = STAT_label[STATL_OBSERVATION];

  plot[index + plot_offset].resize(process->nb_state + 1);

/*  legend.str("");
  legend << STAT_label[STATL_OBSERVATION];
  plot[index + plot_offset][0].legend = legend.str(); */

  plot[index + plot_offset][0].style = "points";

  for (i = 0;i < process->nb_state;i++) {
    legend.str("");
    legend << STAT_label[STATL_STATE] << " " << i << " " << STAT_label[STATL_OBSERVATION] << " "
           << STAT_label[STATL_MODEL];
    plot[index + plot_offset][i + 1].legend = legend.str();

    plot[index + plot_offset][i + 1].style = "lines";
  }

  switch (index_param_type) {

  case IMPLICIT_TYPE : {
    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          plot[index + plot_offset][0].add_point(j , int_sequence[i][variable][j]);
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          plot[index + plot_offset][0].add_point(j , real_sequence[i][variable][j]);
        }
      }
      break;
    }
    }
    break;
  }

  case TIME : {
    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          plot[index + plot_offset][0].add_point(index_parameter[i][j] , int_sequence[i][variable][j]);
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          plot[index + plot_offset][0].add_point(index_parameter[i][j] , real_sequence[i][variable][j]);
        }
      }
      break;
    }
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    plot[index + plot_offset][i + 1].add_point(state_min_index_parameter[i] , process->observation[i]->intercept +
                                               process->observation[i]->slope * state_min_index_parameter[i]);
    plot[index + plot_offset][i + 1].add_point(state_max_index_parameter[i] , process->observation[i]->intercept +
                                               process->observation[i]->slope * state_max_index_parameter[i]);
  }

  switch (type[0]) {
  case STATE :
    index += process->nb_state + 1;
    break;
  default :
    index++;
    break;
  }

  delete [] state_min_index_parameter;
  delete [] state_max_index_parameter;
  if (type[0] == STATE) {
    delete [] state_min_value;
    delete [] state_max_value;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a fitted autocorrelation function for an observation autoregressive model.
 *
 *  \param[in,out] os        stream,
 *  \param[in]     variable  variable index,
 *  \param[in]     process   pointer on a continuous observation process,
 *  \param[in]     file_flag flag file.
 */
/*--------------------------------------------------------------*/

ostream& MarkovianSequences::autoregressive_model_ascii_print(ostream &os , int variable ,
                                                              ContinuousParametricProcess *process ,
                                                              bool file_flag) const

{
  int i , j;
  int max_lag , width[5];
  double standard_normal_value , *confidence_limit;
  Correlation *correl;
  normal dist;
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::right , ios::adjustfield);

  max_lag = max_length * (1. - AUTOCORRELATION_FREQUENCY_RATIO);
  correl = new Correlation(2 , max_lag + 1 , true , PEARSON);
  confidence_limit = new double[max_lag + 1];
  standard_normal_value = quantile(complement(dist , 0.025));

  os << "\n";
  for (i = 0;i < process->nb_state;i++) {
    autocorrelation_computation(*correl , i , variable);

    if (correl->length > 0) {
      correl->point[1][0] = 1.;
      for (j = 1;j < correl->length;j++) {
        correl->point[1][j] = correl->point[1][j - 1] * process->observation[i]->autoregressive_coeff;
      }
      for (j = 0;j < correl->length;j++) {
        confidence_limit[j] = standard_normal_value / sqrt((double)(correl->frequency[j]));
      }

      // computation of the column widths

      width[0] = column_width(correl->length - 1);
      width[1] = column_width(correl->length , correl->point[0]) + ASCII_SPACE;
      width[2] = column_width(correl->length , correl->point[1]) + ASCII_SPACE;
      width[3] = column_width(correl->length , confidence_limit) + ASCII_SPACE;
      width[4] = column_width(correl->frequency[0]) + ASCII_SPACE;

      if (file_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_LAG] << " | " << STAT_label[STATL_STATE] << " " << i << " "
         << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_AUTO] << SEQ_label[SEQL_CORRELATION_FUNCTION]
         << " | " << STAT_label[STATL_STATE] << " " << i << " "
         << SEQ_label[SEQL_THEORETICAL] << " " << SEQ_label[SEQL_AUTO] << SEQ_label[SEQL_CORRELATION_FUNCTION]
         << " | " << SEQ_label[SEQL_RANDOMNESS_95_CONFIDENCE_LIMIT]
         << " | " << STAT_label[STATL_FREQUENCY] << endl;

      for (j = 0;j < correl->length;j++) {
        if (file_flag) {
          os << "# ";
        }
        os << setw(width[0]) << j;
        os << setw(width[1]) << correl->point[0][j];
        os << setw(width[2]) << correl->point[1][j];
        os << setw(width[3]) << confidence_limit[j];
        os << setw(width[4]) << correl->frequency[j] << endl;
      }
      os << endl;
    }

    correl->length = max_lag + 1;
  }

  delete correl;
  delete [] confidence_limit;

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a fitted autocorrelation function for
 *         an observation autoregressive model at the spreadsheet format.
 *
 *  \param[in,out] os       stream,
 *  \param[in]     variable variable index,
 *  \param[in]     process  pointer on a continuous observation process.
 */
/*--------------------------------------------------------------*/

ostream& MarkovianSequences::autoregressive_model_spreadsheet_print(ostream &os , int variable ,
                                                                    ContinuousParametricProcess *process) const

{
  int i , j;
  int max_lag;
  double standard_normal_value , confidence_limit;
  Correlation *correl;
  normal dist;


  max_lag = max_length * (1. - AUTOCORRELATION_FREQUENCY_RATIO);
  correl = new Correlation(2 , max_lag + 1 , true , PEARSON);
  standard_normal_value = quantile(complement(dist , 0.025));

  os << "\n";
  for (i = 0;i < process->nb_state;i++) {
    autocorrelation_computation(*correl , i , variable);

    if (correl->length > 0) {
      correl->point[1][0] = 1.;
      for (j = 1;j < correl->length;j++) {
        correl->point[1][j] = correl->point[1][j - 1] * process->observation[i]->autoregressive_coeff;
      }

      os << SEQ_label[SEQL_LAG] << "\t" << STAT_label[STATL_STATE] << " " << i << " "
         << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_AUTO] << SEQ_label[SEQL_CORRELATION_FUNCTION]
         << "\t" << STAT_label[STATL_STATE] << " " << i << " "
         << SEQ_label[SEQL_THEORETICAL] << " " << SEQ_label[SEQL_AUTO] << SEQ_label[SEQL_CORRELATION_FUNCTION]
         << "\t" << SEQ_label[SEQL_RANDOMNESS_95_CONFIDENCE_LIMIT]
         << "\t" << SEQ_label[SEQL_RANDOMNESS_95_CONFIDENCE_LIMIT]
         << "\t" << STAT_label[STATL_FREQUENCY] << endl;

      for (j = 0;j < correl->length;j++) {
        confidence_limit = standard_normal_value / sqrt((double)(correl->frequency[j]));

        os << j << "\t" << correl->point[0][j] << "\t" << correl->point[1][j] << "\t"
           << confidence_limit << "\t" << -confidence_limit << "\t" << correl->frequency[j] << endl;
      }
      os << endl;
    }

    correl->length = max_lag + 1;
  }

  delete correl;

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a fitted autocorrelation function for
 *         an observation autoregressive model using Gnuplot.
 *
 *  \param[in] prefix   file prefix,
 *  \param[in] title    figure title,
 *  \param[in] variable variable index,
 *  \param[in] process  pointer on a continuous observation process.
 *
 *  \return             error status.
 */
/*--------------------------------------------------------------*/

bool MarkovianSequences::autoregressive_model_plot_print(const char *prefix , const char *title , int variable ,
                                                         ContinuousParametricProcess *process) const

{
  bool status = false , start;
  int i , j;
  int max_lag;
  double standard_normal_value , confidence_limit;
  Correlation **correl;
  normal dist;
  ostringstream data_file_name[NB_STATE];
  ofstream *out_data_file[NB_STATE];


  max_lag = max_length * (1. - AUTOCORRELATION_FREQUENCY_RATIO);
  standard_normal_value = quantile(complement(dist , 0.025));
  correl = new Correlation*[process->nb_state];

  for (i = 0;i < process->nb_state;i++) {
    correl[i] = new Correlation(2 , max_lag + 1 , true , PEARSON);

    autocorrelation_computation(*correl[i] , i , variable);

    if (correl[i]->length > 0) {

      // writing of data files

      correl[i]->point[1][0] = 1.;
      for (j = 1;j < correl[i]->length;j++) {
        correl[i]->point[1][j] = correl[i]->point[1][j - 1] * process->observation[i]->autoregressive_coeff;
      }

      data_file_name[i] << prefix << variable << i << ".dat";
      out_data_file[i] = new ofstream((data_file_name[i].str()).c_str());

      if (out_data_file[i]) {
        status = true;

        for (j = 0;j < correl[i]->length;j++) {
          confidence_limit = standard_normal_value / sqrt((double)(correl[i]->frequency[j]));

          *out_data_file[i] << j << " " << correl[i]->point[0][j] << " " << correl[i]->point[1][j] << " "
                            << confidence_limit << " " << -confidence_limit << " " << correl[i]->frequency[j] << endl;
        }
      }
    }
  }

  if (status) {

    // writing of the script files

    for (i = 0;i < 2;i++) {
      ostringstream file_name[2];

      switch (i) {
      case 0 :
        file_name[0] << prefix << variable << ".plot";
        break;
      case 1 :
        file_name[0] << prefix << variable << ".print";
        break;
      }

      ofstream out_file((file_name[0].str()).c_str());

      if (i == 1) {
        out_file << "set terminal postscript" << endl;
        file_name[1] << label(prefix) << process << ".ps";
        out_file << "set output \"" << file_name[1].str() << "\"\n\n";
      }

      out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
               << "set title \"";
      if (title) {
        out_file << title << " - ";
      }
      out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << variable << "\"\n\n";

      start = true;
      for (j = 0;j < process->nb_state;j++) {
        if (correl[j]->length > 0) {
          if (!start) {
            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;
          }
          else {
            start = false;
          }

          out_file << "set xlabel \"" << SEQ_label[SEQL_LAG] << "\"" << endl;

          if (correl[j]->length - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }

          out_file << "plot [" << 0 << ":" << correl[j]->length - 1 << "] [-1:1] "
                   << "\"" << label((data_file_name[j].str()).c_str()) << "\" using 1:2 title \""
                   << STAT_label[STATL_STATE] << " " << j << " " << SEQ_label[SEQL_OBSERVED] << " "
                   << SEQ_label[SEQL_AUTO] << SEQ_label[SEQL_CORRELATION_FUNCTION]
                   << "\" with linespoints,\\" << endl;
          out_file << "\"" << label((data_file_name[j].str()).c_str()) << "\" using 1:3 title \""
                   << STAT_label[STATL_STATE] << " " << j << " " << SEQ_label[SEQL_THEORETICAL] << " "
                   << SEQ_label[SEQL_AUTO] << SEQ_label[SEQL_CORRELATION_FUNCTION]
                   << "\" with linespoints,\\" << endl;
          out_file << "\"" << label((data_file_name[j].str()).c_str()) << "\" using 1:4 notitle with lines,\\" << endl;
          out_file << "\"" << label((data_file_name[j].str()).c_str()) << "\" using 1:5 notitle with lines" << endl;

          out_file << "set xlabel" << endl;

          if (correl[j]->length - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
        }
      }

      out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
    }
  }

  for (i = 0;i < process->nb_state;i++) {
    delete correl[i];
  }
  delete [] correl;

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a fitted autocorrelation function for an observation autoregressive model.
 *
 *  \param[in] plot     file prefix,
 *  \param[in] index    MultiPlot index,
 *  \param[in] variable variable index,
 *  \param[in] process  pointer on a continuous observation process.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::autoregressive_model_plotable_write(MultiPlotSet &plot , int &index , int variable ,
                                                             ContinuousParametricProcess *process) const

{
  bool status = false;
  int i , j;
  int max_lag;
  double standard_normal_value , confidence_limit;
  Correlation *correl;
  normal dist;
  ostringstream title , legend;


  max_lag = max_length * (1. - AUTOCORRELATION_FREQUENCY_RATIO);
  correl = new Correlation(2 , max_lag + 1 , true , PEARSON);
  standard_normal_value = quantile(complement(dist , 0.025));

  plot.variable_nb_viewpoint[variable] = 1;

  title.str("");
  title << STAT_label[STATL_OUTPUT_PROCESS] << " " << variable;

  for (i = 0;i < process->nb_state;i++) {
    autocorrelation_computation(*correl , i , variable);

    if (correl->length > 0) {
      correl->point[1][0] = 1.;
      for (j = 1;j < correl->length;j++) {
        correl->point[1][j] = correl->point[1][j - 1] * process->observation[i]->autoregressive_coeff;
      }

      plot[index].title = title.str();

      plot[index].xrange = Range(0 , correl->length - 1);
      plot[index].yrange = Range(-1 , 1);

      if (correl->length - 1 < TIC_THRESHOLD) {
        plot[index].xtics = 1;
      }

      plot[index].xlabel = SEQ_label[SEQL_LAG];
      plot[index].ylabel = STAT_label[STATL_CORRELATION_COEFF];

      plot[index].resize(4);

      legend.str("");
      legend << STAT_label[STATL_STATE] << " " << i << SEQ_label[SEQL_OBSERVED] << " "
             << SEQ_label[SEQL_AUTO] << SEQ_label[SEQL_CORRELATION_FUNCTION];
      plot[index][0].legend = legend.str();

      plot[index][0].style = "linespoints";

      legend.str("");
      legend << STAT_label[STATL_STATE] << " " << i << SEQ_label[SEQL_THEORETICAL] << " "
             << SEQ_label[SEQL_AUTO] << SEQ_label[SEQL_CORRELATION_FUNCTION];
      plot[index][1].legend = legend.str();

      plot[index][1].style = "linespoints";

      plot[index][2].style = "lines";
      plot[index][3].style = "lines";

      for (j = 0;j < correl->length;j++) {
        confidence_limit = standard_normal_value / sqrt((double)(correl->frequency[j]));

        plot[index][0].add_point(j , correl->point[0][j]);
        plot[index][1].add_point(j , correl->point[1][j]);
        plot[index][2].add_point(j , confidence_limit);
        plot[index][3].add_point(j , -confidence_limit);
      }
 
      index++;
    }

    correl->length = max_lag + 1;
  }

  delete correl;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a MarkovianSequences object in a file at the MTG format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path,
 *  \param[in] itype variable types (NOMINAL/NUMERIC).
 */
/*--------------------------------------------------------------*/

bool MarkovianSequences::mtg_write(StatError &error , const string path , variable_type *itype) const

{
  bool status;
  int i , j , k , m;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    // writing of the header

    out_file << "CODE:\tFORM-A" << endl;

    out_file << "\nCLASSES:\nSYMBOL\tSCALE\tDECOMPOSITION\tINDEXATION\tDEFINITION" << endl;
    out_file << "$\t0\tFREE\tFREE\tIMPLICIT" << endl;
    out_file << "U\t1\tFREE\tFREE\tIMPLICIT" << endl;
    out_file << "E\t2\tFREE\tFREE\tIMPLICIT" << endl;

    for (i = 0;i < nb_variable;i++) {
      switch (itype[i]) {

      case NOMINAL : {
        for (j = 1;j < marginal_distribution[i]->nb_value;j++) {
          out_file << (char)('F' + j) << "\t2\tFREE\tFREE\tIMPLICIT" << endl;
        }
        break;
      }

      case NUMERIC : {
        out_file << "F\t2\tFREE\tFREE\tIMPLICIT" << endl;
        break;
      }
      }
    }

    out_file << "\nDESCRIPTION:\nLEFT\tRIGHT\tRELTYPE\tMAX" << endl;
    out_file << "E\tE\t<\t1" << endl;

    for (i = 0;i < nb_variable;i++) {
      switch (itype[i]) {

      case NOMINAL : {
        out_file << "E\t";
        for (j = 1;j < marginal_distribution[i]->nb_value;j++) {
          out_file << (char)('F' + j);
          if (j < marginal_distribution[i]->nb_value - 1) {
            out_file << ",";
          }
        }
        out_file << "\t+\t1" << endl;
        break;
      }

      case NUMERIC : {
        out_file << "E\tF\t+\t?" << endl;
        break;
      }
      }
    }

    out_file << "\nFEATURES:\nNAME\tTYPE" << endl;

    // writing of the topological code

    out_file << "\nMTG:\nENTITY-CODE\n" << endl;

    for (i = 0;i < nb_sequence;i++) {
      out_file << "/U" << i + 1 << endl;

      for (j = 0;j < length[i];j++) {
        if (j == 0) {
          out_file << "\t/";
        }
        else {
          out_file << "\t^<";
        }
        out_file << 'E' << j + 1 << endl;

        for (k = 0;k < nb_variable;k++) {
          switch (itype[k]) {

          case NOMINAL : {
            if (int_sequence[i][k][j] > 0) {
              out_file <<"\t\t+" << (char)('F' + int_sequence[i][k][j]) << 1 << endl;
            }
            break;
          }

          case NUMERIC : {
            for (m = 0;m < int_sequence[i][k][j];m++) {
              out_file <<"\t\t+F" << m + 1 << endl;
            }
            break;
          }
          }
        }
      }
    }
  }

  return status;
}


};  // namespace sequence_analysis
