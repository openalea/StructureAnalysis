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



#include <limits.h>
#include <math.h>

#include <string>
#include <vector>
#include <sstream>
#include <iomanip>

#include <boost/tokenizer.hpp>

#include "distribution.h"
#include "stat_label.h"

using namespace std;
using namespace boost;


namespace stat_tool {


extern bool cumul_matching_plot_print(const char *path , int nb_cumul , int *offset ,
                                      int *nb_value , double **cumul);



/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the results of a comparison of frequency distributions.
 *
 *  \param[in,out] os            stream,
 *  \param[in]     nb_histo      number of frequency distributions,
 *  \param[in]     ihisto        pointer on the frequency distributions,
 *  \param[in]     type          variable type (NOMINAL/ORDINAL/NUMERIC),
 *  \param[in]     dissimilarity dissimilarities.
 */
/*--------------------------------------------------------------*/

ostream& FrequencyDistribution::dissimilarity_ascii_write(ostream &os , int nb_histo ,
                                                          const FrequencyDistribution **ihisto ,
                                                          variable_type type , double **dissimilarity) const

{
  int i , j;
  int max_nb_value , buff , width[3];
  double information , **cumul;
  Test *test;
  const FrequencyDistribution **histo;
  ios_base::fmtflags format_flags;


  nb_histo++;
  histo = new const FrequencyDistribution*[nb_histo];

  histo[0] = this;
  for (i = 1;i < nb_histo;i++) {
    histo[i] = ihisto[i - 1];
  }

  format_flags = os.setf(ios::right , ios::adjustfield);

  // writing of the frequency distribution characteristics

  for (i = 0;i < nb_histo;i++) {
    os << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i + 1  << " - ";

    if (type == NOMINAL) {
      os << STAT_label[STATL_SAMPLE_SIZE] << ": " << histo[i]->nb_element << endl;
    }

    else {
      histo[i]->ascii_characteristic_print(os , true);

      os << STAT_label[STATL_MEAN_ABSOLUTE_DEVIATION] << ": " << histo[i]->mean_absolute_deviation_computation(histo[i]->mean);
      if (histo[i]->mean > 0.) {
        os << "   " << STAT_label[STATL_CONCENTRATION_COEFF] << ": " << histo[i]->concentration_computation();
      }
      os << endl;
    }

    information = histo[i]->information_computation();

    os << STAT_label[STATL_INFORMATION] << ": " << information << " ("
       << information / histo[i]->nb_element << ")\n" << endl;
  }

  cumul = new double*[nb_histo];
  for (i = 0;i < nb_histo;i++) {
    cumul[i] = histo[i]->cumul_computation();
  }

  // computation of the column widths

  max_nb_value = histo[0]->nb_value;
  for (i = 1;i < nb_histo;i++) {
    if (histo[i]->nb_value > max_nb_value) {
      max_nb_value = histo[i]->nb_value;
    }
  }

  width[0] = column_width(max_nb_value - 1);

  width[1] = 0;
  for (i = 0;i < nb_histo;i++) {
    buff = column_width(histo[i]->max);
    if (buff > width[1]) {
      width[1] = buff;
    }
  }
  width[1] += ASCII_SPACE;

  width[2] = 0;
  for (i = 0;i < nb_histo;i++) {
    buff = column_width(histo[i]->nb_value - histo[i]->offset , cumul[i] + histo[i]->offset);
    if (buff > width[2]) {
      width[2] = buff;
    }
  }
  width[2] += ASCII_SPACE;

  // writing of the frequency distributions and the cumulative distribution functions

  os << "  ";
  for (i = 0;i < nb_histo;i++) {
    os << " | " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i + 1;
  }
  for (i = 0;i < nb_histo;i++) {
    os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
       << i + 1 << " " << STAT_label[STATL_FUNCTION];
  }
  os << endl;

  for (i = 0;i < max_nb_value;i++) {
    os << setw(width[0]) << i;
    for (j = 0;j < nb_histo;j++) {
      if (i < histo[j]->nb_value) {
        os << setw(width[1]) << histo[j]->frequency[i];
      }
      else {
        os << setw(width[1]) << " ";
      }
    }
    for (j = 0;j < nb_histo;j++) {
      if (i < histo[j]->nb_value) {
        os << setw(width[2]) << cumul[j][i];
      }
      else {
        os << setw(width[2]) << " ";
      }
    }
    os << endl;
  }

  for (i = 0;i < nb_histo;i++) {
    delete [] cumul[i];
  }
  delete [] cumul;

  // computation of the column widths

  width[0] = column_width(nb_histo);

  width[1] = 0;
  for (i = 0;i < nb_histo;i++) {
    buff = column_width(nb_histo , dissimilarity[i]);
    if (buff > width[1]) {
      width[1] = buff;
    }
  }
  width[1] += ASCII_SPACE;

  // writing of the pairwise dissimilarity matrix between frequency distributions

  os << "\n" << STAT_label[STATL_DISSIMILARITIES] << " between "
     << STAT_label[STATL_FREQUENCY_DISTRIBUTIONS] << endl;

  os << "\n           ";
  for (i = 0;i < nb_histo;i++) {
    os << " | " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i + 1;
  }
  os << endl;

  for (i = 0;i < nb_histo;i++) {
    os << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " ";
    os << setw(width[0]) << i + 1 << " ";
    for (j = 0;j < nb_histo;j++) {
      os << setw(width[1]) << dissimilarity[i][j];
    }
    os << endl;
  }

  // analysis of variance

  switch (type) {

  case ORDINAL : {
    test = histo[0]->kruskal_wallis_test(nb_histo - 1 , histo + 1);

    os << "\n" << STAT_label[STATL_KRUSKAL_WALLIS_TEST] << endl;
    test->ascii_print(os);
    break;
  }

  case NUMERIC : {
    int df[3];
    double diff , square_sum[3] , mean_square[3];
    FrequencyDistribution *merged_histo;


    merged_histo = new FrequencyDistribution(nb_histo , histo);

    square_sum[0] = 0.;
    square_sum[1] = 0.;

    for (i = 0;i < nb_histo;i++) {
      diff = histo[i]->mean - merged_histo->mean;
      square_sum[0] += diff * diff * histo[i]->nb_element;
      square_sum[1] += histo[i]->variance * (histo[i]->nb_element - 1);
    }

    square_sum[2] = merged_histo->variance * (merged_histo->nb_element - 1);

    df[0] = nb_histo - 1;
    df[1] = merged_histo->nb_element - nb_histo;
    df[2] = merged_histo->nb_element - 1;

    for (i = 0;i < 3;i++) {
      mean_square[i] = square_sum[i] / df[i];
    }

    delete merged_histo;

#   ifdef DEBUG
    os << "\ntest: " << square_sum[0] + square_sum[1] << " | " << square_sum[2] << endl;
#   endif

    width[0] = column_width(merged_histo->nb_element - 1) + ASCII_SPACE;
    width[1] = column_width(3 , square_sum) + ASCII_SPACE;
    width[2] = column_width(3 , mean_square) + ASCII_SPACE;

    os << "\n" << STAT_label[STATL_VARIANCE_ANALYSIS] << endl;
    os << "\n" << STAT_label[STATL_VARIATION_SOURCE] << " | " << STAT_label[STATL_FREEDOM_DEGREES]
       << " | " << STAT_label[STATL_SQUARE_SUM] << " | " << STAT_label[STATL_MEAN_SQUARE] << endl;
    for (i = 0;i < 3;i++) {
      switch (i) {
      case 0 :
        os << STAT_label[STATL_BETWEEN_SAMPLES];
        break;
      case 1 :
        os << STAT_label[STATL_WITHIN_SAMPLES] << " ";
        break;
      case 2 :
        os << STAT_label[STATL_TOTAL] << "          ";
        break;
      }

      os << setw(width[0]) << df[i];
      os << setw(width[1]) << square_sum[i];
      os << setw(width[2]) << mean_square[i] << endl;
    }
    os << endl;

    if ((df[0] > 0) && (df[1] > 0)) {
      test = new Test(FISHER , true , df[0] , df[1] , mean_square[0] / mean_square[1]);
      test->F_critical_probability_computation();

      test->ascii_print(os , false , (df[0] == 1 ? false : true));

      delete test;
    }
    break;
  }
  }

  os.setf(format_flags , ios::adjustfield);

  delete [] histo;

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the results of a comparison of frequency distributions in a file.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] path          file path,
 *  \param[in] nb_histo      number of frequency distributions,
 *  \param[in] ihisto        pointer on the frequency distributions,
 *  \param[in] type          variable type (NOMINAL/ORDINAL/NUMERIC),
 *  \param[in] dissimilarity dissimilarities.
 *
 *  \return                  error status.
 */
/*--------------------------------------------------------------*/

bool FrequencyDistribution::dissimilarity_ascii_write(StatError &error , const string path ,
                                                      int nb_histo , const FrequencyDistribution **ihisto ,
                                                      variable_type type , double **dissimilarity) const

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
    dissimilarity_ascii_write(out_file , nb_histo , ihisto , type , dissimilarity);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the results of a comparison of frequency distributions in a file
 *         at the spreadsheet format.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] path          file path,
 *  \param[in] nb_histo      number of frequency distributions,
 *  \param[in] ihisto        pointer on the frequency distributions,
 *  \param[in] type          variable type (NOMINAL/ORDINAL/NUMERIC),
 *  \param[in] dissimilarity dissimilarities.
 *
 *  \return                  error status.
 */
/*--------------------------------------------------------------*/

bool FrequencyDistribution::dissimilarity_spreadsheet_write(StatError &error , const string path ,
                                                            int nb_histo , const FrequencyDistribution **ihisto ,
                                                            variable_type type , double **dissimilarity) const

{
  bool status;
  int i , j;
  int max_nb_value;
  double information , **cumul , **concentration;
  Test *test;
  const FrequencyDistribution **histo;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    nb_histo++;
    histo = new const FrequencyDistribution*[nb_histo];

    histo[0] = this;
    for (i = 1;i < nb_histo;i++) {
      histo[i] = ihisto[i - 1];
    }

    // writing of the frequency distribution characteristics

    for (i = 0;i < nb_histo;i++) {
      out_file << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i + 1  << "\t";

      if (type == NOMINAL) {
        out_file << STAT_label[STATL_SAMPLE_SIZE] << "\t" << histo[i]->nb_element << endl;
      }

      else {
        histo[i]->spreadsheet_characteristic_print(out_file , true);

        out_file << STAT_label[STATL_MEAN_ABSOLUTE_DEVIATION] << "\t" << histo[i]->mean_absolute_deviation_computation(histo[i]->mean);
        if (histo[i]->mean > 0.) {
          out_file << "\t\t" << STAT_label[STATL_CONCENTRATION_COEFF] << "\t" << histo[i]->concentration_computation();
        }
        out_file << endl;
      }

      information = histo[i]->information_computation();

      out_file << STAT_label[STATL_INFORMATION] << "\t" << information << "\t"
               << information / histo[i]->nb_element << "\n" << endl;
    }

    // writing of the frequency distributions, cumulative distribution functions and concentration curves

    cumul = new double*[nb_histo];
    for (i = 0;i < nb_histo;i++) {
      cumul[i] = histo[i]->cumul_computation();
    }

    if (type != NOMINAL) {
      concentration = new double*[nb_histo];
      for (i = 0;i < nb_histo;i++) {
        concentration[i] = histo[i]->concentration_function_computation();
      }
    }

    max_nb_value = histo[0]->nb_value;
    for (i = 1;i < nb_histo;i++) {
      if (histo[i]->nb_value > max_nb_value) {
        max_nb_value = histo[i]->nb_value;
      }
    }

    for (i = 0;i < nb_histo;i++) {
      out_file << "\t" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i + 1;
    }
    for (i = 0;i < nb_histo;i++) {
      out_file << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
               << i + 1 << " " << STAT_label[STATL_FUNCTION];
    }

    if (type != NOMINAL) {
      for (i = 0;i < nb_histo;i++) {
        if (histo[i]->variance > 0.) {
          out_file << "\t" << STAT_label[STATL_CONCENTRATION] << " " << STAT_label[STATL_FUNCTION] << " "
                   << i + 1;
        }
      }
    }
    out_file << endl;

    for (i = 0;i < max_nb_value;i++) {
      out_file << i;
      for (j = 0;j < nb_histo;j++) {
        out_file << "\t";
        if (i < histo[j]->nb_value) {
          out_file << histo[j]->frequency[i];
        }
      }
      for (j = 0;j < nb_histo;j++) {
        out_file << "\t";
        if (i < histo[j]->nb_value) {
          out_file << cumul[j][i];
        }
      }

      if (type != NOMINAL) {
        for (j = 0;j < nb_histo;j++) {
          if (histo[j]->variance > 0.) {
            out_file << "\t";
            if (i < histo[j]->nb_value) {
              out_file << concentration[j][i];
            }
          }
        }
      }
      out_file << endl;
    }

    for (i = 0;i < nb_histo;i++) {
      delete [] cumul[i];
    }
    delete [] cumul;

    if (type != NOMINAL) {
      for (i = 0;i < nb_histo;i++) {
        delete [] concentration[i];
      }
      delete [] concentration;
    }

    // writing of the pairwise dissimilarity matrix between frequency distributions

    out_file << "\n" << STAT_label[STATL_DISSIMILARITIES] << " between "
             << STAT_label[STATL_FREQUENCY_DISTRIBUTIONS] << "\n\n";

    for (i = 0;i < nb_histo;i++) {
      out_file << "\t" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i + 1;
    }
    out_file << endl;

    for (i = 0;i < nb_histo;i++) {
      out_file << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i + 1;
      for (j = 0;j < nb_histo;j++) {
        out_file << "\t" << dissimilarity[i][j];
      }
      out_file << endl;
    }

    // analysis of variance

    switch (type) {

    case ORDINAL : {
      test = histo[0]->kruskal_wallis_test(nb_histo - 1 , histo + 1);

      out_file << "\n" << STAT_label[STATL_KRUSKAL_WALLIS_TEST] << endl;
      test->spreadsheet_print(out_file);
      break;
    }

    case NUMERIC : {
      int df[3];
      double diff , square_sum[3] , mean_square[3];
      FrequencyDistribution *merged_histo;


      merged_histo = new FrequencyDistribution(nb_histo , histo);

      square_sum[0] = 0.;
      square_sum[1] = 0.;

      for (i = 0;i < nb_histo;i++) {
        diff = histo[i]->mean - merged_histo->mean;
        square_sum[0] += diff * diff * histo[i]->nb_element;
        square_sum[1] += histo[i]->variance * (histo[i]->nb_element - 1);
      }

      square_sum[2] = merged_histo->variance * (merged_histo->nb_element - 1);

      df[0] = nb_histo - 1;
      df[1] = merged_histo->nb_element - nb_histo;
      df[2] = merged_histo->nb_element - 1;

      for (i = 0;i < 3;i++) {
        mean_square[i] = square_sum[i] / df[i];
      }

      delete merged_histo;

      out_file << "\n" << STAT_label[STATL_VARIANCE_ANALYSIS] << endl;
      out_file << "\n" << STAT_label[STATL_VARIATION_SOURCE] << "\t" << STAT_label[STATL_FREEDOM_DEGREES]
               << "\t" << STAT_label[STATL_SQUARE_SUM] << "\t" << STAT_label[STATL_MEAN_SQUARE] << endl;
      for (i = 0;i < 3;i++) {
        switch (i) {
        case 0 :
          out_file << STAT_label[STATL_BETWEEN_SAMPLES];
          break;
        case 1 :
          out_file << STAT_label[STATL_WITHIN_SAMPLES];
          break;
        case 2 :
          out_file << STAT_label[STATL_TOTAL];
          break;
        }

        out_file << "\t" << df[i] << "\t" << square_sum[i] << "\t" << mean_square[i] << endl;
      }
      out_file << endl;

      if ((df[0] > 0) && (df[1] > 0)) {
        test = new Test(FISHER , true , df[0] , df[1] , mean_square[0] / mean_square[1]);
        test->F_critical_probability_computation();

        test->spreadsheet_print(out_file , (df[0] == 1 ? false : true));

        delete test;
      }
      break;
    }
    }

    delete [] histo;
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Kruskal-Wallis test (analysis of variance on ranks).
 *
 *  \param[in] nb_histo number of frequency distributions,
 *  \param[in] ihisto   pointer on the frequency distributions.
 *
 *  \return             Test object.
 */
/*--------------------------------------------------------------*/

Test* FrequencyDistribution::kruskal_wallis_test(int nb_histo , const FrequencyDistribution **ihisto) const

{
  int i , j;
  int *pfrequency;
  double correction , value , sum , *rank;
  const FrequencyDistribution **histo;
  FrequencyDistribution *merged_histo;
  Test *test;


  nb_histo++;
  histo = new const FrequencyDistribution*[nb_histo];

  histo[0] = this;
  for (i = 1;i < nb_histo;i++) {
    histo[i] = ihisto[i - 1];
  }

  merged_histo = new FrequencyDistribution(nb_histo , histo);

  // computation of the correction term for ties

  pfrequency = merged_histo->frequency + merged_histo->offset;
  correction = 0.;
  for (i = merged_histo->offset;i < merged_histo->nb_value;i++) {
    if (*pfrequency > 1) {
      correction += *pfrequency * ((double)*pfrequency * (double)*pfrequency - 1);
    }
    pfrequency++;
  }

  // rank computation

  rank = merged_histo->rank_computation();

  // computation of the Kruskal-Wallis statistic

  value = 0.;
  for (i = 0;i < nb_histo;i++) {
    sum = 0.;
    for (j = histo[i]->offset;j < histo[i]->nb_value;j++) {
      sum += histo[i]->frequency[j] * rank[j];
    }
    value += sum * sum / histo[i]->nb_element;
  }

  value = (12 * value / (merged_histo->nb_element * ((double)merged_histo->nb_element + 1)) -
           3 * (merged_histo->nb_element + 1)) / (1. - correction / (merged_histo->nb_element *
           ((double)merged_histo->nb_element * (double)merged_histo->nb_element - 1)));

  test = new Test(CHI2 , true , nb_histo - 1 , I_DEFAULT , value);

  test->chi2_critical_probability_computation();

  delete [] histo;
  delete merged_histo;
  delete [] rank;

  return test;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Comparison of frequency distributions.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] display  flag for displaying comparison outputs,
 *  \param[in] nb_histo number of frequency distributions,
 *  \param[in] ihisto   pointer on the frequency distributions,
 *  \param[in] type     variable type (NOMINAL/ORDINAL/NUMERIC),
 *  \param[in] path     file path,
 *  \param[in] format   file format (ASCII/SPREADSHEET).
 *
 *  \return             error status.
 */
/*--------------------------------------------------------------*/

bool FrequencyDistribution::comparison(StatError &error , bool display , int nb_histo ,
                                       const FrequencyDistribution **ihisto , variable_type type ,
                                       const string path , output_format format) const

{
  bool status = true;
  int i , j , k , m;
  int max_nb_value , distance = 1;
  double **dissimilarity;
  Distribution **dist;
  const FrequencyDistribution **histo;


  nb_histo++;
  histo = new const FrequencyDistribution*[nb_histo];

  histo[0] = this;
  max_nb_value = nb_value;
  for (i = 1;i < nb_histo;i++) {
    histo[i] = ihisto[i - 1];
    if (histo[i]->nb_value > max_nb_value) {
      max_nb_value = histo[i]->nb_value;
    }
  }

  dissimilarity = new double*[nb_histo];
  dist = new Distribution*[nb_histo];

  for (i = 0;i < nb_histo;i++) {
    dissimilarity[i] = new double[nb_histo];

    dist[i] = new Distribution(max_nb_value);
    histo[i]->distribution_estimation(dist[i]);
    for (j = histo[i]->nb_value;j < max_nb_value;j++) {
      dist[i]->mass[j] = 0.;
    }
  }

# ifdef DEBUG
  if (type == ORDINAL) {
    int buff , width[2];
    double *rank , *rank_mean;
    FrequencyDistribution *merged_histo;


    merged_histo = new FrequencyDistribution(nb_histo , histo);

    // rank computation

    rank = merged_histo->rank_computation();

    // computation of the mean rank for each frequency distribution

    rank_mean = new double[nb_histo];

    for (i = 0;i < nb_histo;i++) {
      rank_mean[i] = 0.;
      for (j = histo[i]->offset;j < histo[i]->nb_value;j++) {
        rank_mean[i] += histo[i]->frequency[j] * rank[j];
      }
      rank_mean[i] /= histo[i]->nb_element;
    }

    for (i = 0;i < nb_histo;i++) {
      dissimilarity[i][i] = 0.;
      for (j = i + 1;j < nb_histo;j++) {
        dissimilarity[i][j] = rank_mean[j] - rank_mean[i];
        dissimilarity[j][i] = -dissimilarity[i][j];
      }
    }

    // computation of the column widths

    width[0] = column_width(nb_histo);

    width[1] = 0;
    for (i = 0;i < nb_histo;i++) {
      buff = column_width(nb_histo , dissimilarity[i]);
      if (buff > width[1]) {
        width[1] = buff;
      }
    }
    width[1] += ASCII_SPACE;

    // writing of the pairwise dissimilarity matrix between frequency distributions

    cout << "\n" << STAT_label[STATL_DISSIMILARITIES] << " between "
         << STAT_label[STATL_FREQUENCY_DISTRIBUTIONS] << endl;

    cout << "\n           ";
    for (i = 0;i < nb_histo;i++) {
      cout << " | " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i + 1;
    }
    cout << endl;

    for (i = 0;i < nb_histo;i++) {
      cout << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " ";
      cout << setw(width[0]) << i + 1 << " ";
      for (j = 0;j < nb_histo;j++) {
        cout << setw(width[1]) << dissimilarity[i][j];
      }
      cout << endl;
    }
    cout << endl;

    delete merged_histo;
    delete [] rank;
    delete [] rank_mean;
  }
# endif

  for (i = 0;i < nb_histo;i++) {
    dissimilarity[i][i] = 0.;

    for (j = i + 1;j < nb_histo;j++) {
      dissimilarity[i][j] = 0.;

      for (k = MIN(histo[i]->offset , histo[j]->offset);k < MAX(histo[i]->nb_value , histo[j]->nb_value);k++) {
        for (m = k + 1;m < MAX(histo[i]->nb_value , histo[j]->nb_value);m++) {
          if (type == NOMINAL) {
            dissimilarity[i][j] += fabs(dist[i]->mass[k] * dist[j]->mass[m] -
                                        dist[i]->mass[m] * dist[j]->mass[k]);
          }

          else {
            if (type == NUMERIC) {
              distance = m - k;
            }
            dissimilarity[i][j] += (dist[i]->mass[k] * dist[j]->mass[m] -
                                    dist[i]->mass[m] * dist[j]->mass[k]) * distance;
          }
        }
      }

      if (type == NOMINAL) {
        dissimilarity[j][i] = dissimilarity[i][j];
      }
      else {
        dissimilarity[j][i] = -dissimilarity[i][j];
      }
    }
  }

  if (display) {
    dissimilarity_ascii_write(cout , nb_histo - 1 , ihisto , type , dissimilarity);
  }

  if (!path.empty()) {
    switch (format) {
    case ASCII :
      status = dissimilarity_ascii_write(error , path , nb_histo - 1 , ihisto ,
                                         type , dissimilarity);
      break;
    case SPREADSHEET :
      status = dissimilarity_spreadsheet_write(error , path , nb_histo - 1 , ihisto ,
                                               type , dissimilarity);
      break;
    }
  }

  delete [] histo;

  for (i = 0;i < nb_histo;i++) {
    delete [] dissimilarity[i];
    delete dist[i];
  }
  delete [] dissimilarity;
  delete [] dist;

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Comparison of frequency distributions.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] display  flag for displaying comparison outputs,
 *  \param[in] nb_histo number of frequency distributions,
 *  \param[in] ihisto   pointer on the frequency distributions,
 *  \param[in] type     variable type (NOMINAL/ORDINAL/NUMERIC),
 *  \param[in] path     file path,
 *  \param[in] format   file format (ASCII/SPREADSHEET).
 *
 *  \return             error status.
 */
/*--------------------------------------------------------------*/

bool FrequencyDistribution::comparison(StatError &error , bool display , int nb_histo ,
                                       const vector<FrequencyDistribution> ihisto , variable_type type ,
                                       const string path , output_format format) const

{
  bool status;
  int i;
  const FrequencyDistribution **histo;


  histo = new const FrequencyDistribution*[nb_histo];
  for (i = 0;i < nb_histo;i++) {
    histo[i] = new FrequencyDistribution(ihisto[i]);
  }

  status = comparison(error , display , nb_histo , histo , type , path , format);

  for (i = 0;i < nb_histo;i++) {
    delete histo[i];
  }
  delete [] histo;

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief F test of variance comparison.
 *
 *  \param[in] display flag for displaying the test results,
 *  \param[in] histo   reference on a frequency distribution.
 */
/*--------------------------------------------------------------*/

void FrequencyDistribution::F_comparison(bool display , const FrequencyDistribution &histo) const

{
  if ((nb_element > 1) && (histo.nb_element > 1)) {
    Test *test;


    if (variance > histo.variance) {
      test = new Test(FISHER , true , nb_element - 1 , histo.nb_element - 1 ,
                      variance / histo.variance);
    }
    else {
      test = new Test(FISHER , true , histo.nb_element - 1 , nb_element - 1 ,
                      histo.variance / variance);
    }

    test->F_critical_probability_computation();

    if (display) {
      cout << *test;
    }

    delete test;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Student's t test of mean comparison.
 *
 *  \param[in] display flag for displaying the test results,
 *  \param[in] histo   reference on a frequency distribution.
 */
/*--------------------------------------------------------------*/

void FrequencyDistribution::t_comparison(bool display , const FrequencyDistribution &histo) const

{
  int df;
  double value , buff1 , buff2;
  Test *test;


  buff1 = variance / nb_element;
  buff2 = histo.variance / histo.nb_element;

  value = fabs(mean - histo.mean) / sqrt(buff1 + buff2);
  df = (int)round((buff1 + buff2) * (buff1 + buff2) /
                  (buff1 * buff1 / (nb_element - 1) + buff2 * buff2 / (histo.nb_element - 1)));

  if (df > 0) {
    test = new Test(STUDENT , false , df , I_DEFAULT , value);

    test->t_critical_probability_computation();

    if (display) {
      cout << *test;
    }

    delete test;
  }

# ifdef DEBUG
  value = fabs(mean - histo.mean) /
          sqrt(((nb_element - 1) * variance + (histo.nb_element - 1) * histo.variance) /
               (nb_element + histo.nb_element - 2) *
               (1. / (double)nb_element + 1. / (double)histo.nb_element));

  if (nb_element + histo.nb_element > 2) {
    test = new Test(STUDENT , false , nb_element + histo.nb_element - 2 , I_DEFAULT , value);

    test->t_critical_probability_computation();

    cout << "\n" << *test;

    delete test;
  }
# endif

}


/*--------------------------------------------------------------*/
/**
 *  \brief Wilcoxon-Mann-Whitney test of distribution comparison.
 *
 *  \param[in] error   reference on a StatError object,
 *  \param[in] display flag for displaying the test results,
 *  \param[in] ihisto  reference on a frequency distribution.
 *
 *  \return            error status.
 */
/*--------------------------------------------------------------*/

bool FrequencyDistribution::wilcoxon_mann_whitney_comparison(StatError &error , bool display ,
                                                             const FrequencyDistribution &ihisto) const

{
  bool status;
  int i;
  int nb_equal , min , max , *pfrequency;
  double correction , value , nb_sup , mean , variance , *rank;
  const FrequencyDistribution **histo;
  FrequencyDistribution *merged_histo;
  Test *test;


  if (nb_element * (double)ihisto.nb_element > INT_MAX) {
    status = false;
    error.update(STAT_error[STATR_SAMPLE_SIZES]);
  }

  else {
    status = true;

    histo = new const FrequencyDistribution*[2];

    histo[0] = this;
    histo[1] = &ihisto;

    merged_histo = new FrequencyDistribution(2 , histo);

    // computation of the correction term for ties

    pfrequency = merged_histo->frequency + merged_histo->offset;
    correction = 0.;
    for (i = merged_histo->offset;i < merged_histo->nb_value;i++) {
      if (*pfrequency > 1) {
        correction += *pfrequency * ((double)*pfrequency * (double)*pfrequency - 1);
      }
      pfrequency++;
    }

    // rank computation

    rank = merged_histo->rank_computation();

    // computation of the Wilcoxon-Mann-Whitney statistic

    value = 0.;
    for (i = histo[0]->offset;i < histo[0]->nb_value;i++) {
      value += histo[0]->frequency[i] * rank[i];
    }

    nb_sup = value - histo[0]->nb_element * ((double)histo[0]->nb_element + 1) / 2.;

    mean = histo[0]->nb_element * ((double)merged_histo->nb_element + 1) / 2.;

    // continuity correction

    if (value < mean) {
      value += 0.5;
    }
    else {
      value -= 0.5;
    }

    variance = histo[0]->nb_element * (double)histo[1]->nb_element *
               ((merged_histo->nb_element * ((double)merged_histo->nb_element *
                  (double)merged_histo->nb_element - 1)) - correction) /
               (merged_histo->nb_element * ((double)merged_histo->nb_element - 1) * 12.);

    value = fabs(value - mean) / sqrt(variance);

    test = new Test(STANDARD_NORMAL , false , I_DEFAULT , I_DEFAULT , value);

    test->standard_normal_critical_probability_computation();

    min = MAX(histo[0]->offset , histo[1]->offset);
    max = MIN(histo[0]->nb_value , histo[1]->nb_value) - 1;
    nb_equal = 0;

    if (max >= min) {
      for (i = min;i <= max;i++) {
        nb_equal += histo[0]->frequency[i] * histo[1]->frequency[i];
      }
    }

    if (display) {
      cout << STAT_label[STATL_TWO_SIDED] << " " << STAT_label[STATL_WILCOXON_MANN_WHITNEY_TEST];
      cout << *test;

      cout << STAT_label[STATL_MANN_WHITNEY_INFERIOR_PROBABILITY] << " = "
           << (histo[0]->nb_element * (double)histo[1]->nb_element - nb_equal / 2. - nb_sup) /
              (histo[0]->nb_element * (double)histo[1]->nb_element) << "   "
           << STAT_label[STATL_MANN_WHITNEY_EQUAL_PROBABILITY] << " = "
           << nb_equal / (histo[0]->nb_element * (double)histo[1]->nb_element) << "   "
           << STAT_label[STATL_MANN_WHITNEY_SUPERIOR_PROBABILITY] << " = "
           << (nb_sup - nb_equal / 2.) / (histo[0]->nb_element * (double)histo[1]->nb_element) << endl;
    }

    delete [] histo;
    delete merged_histo;
    delete [] rank;

    delete test;
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a family of frequency distributions using Gnuplot:
 *         - frequency distributions and cumulative frequencies,
 *         - probability mass functions and cumulative distribution functions,
 *         - matching of cumulative distribution functions,
 *         - concentration curves.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] prefix   file prefix,
 *  \param[in] nb_histo number of frequency distributions,
 *  \param[in] ihisto   pointer on the frequency distributions,
 *  \param[in] title    figure title.
 *
 *  \return             error status.
 */
/*--------------------------------------------------------------*/

bool FrequencyDistribution::plot_write(StatError &error , const char *prefix , int nb_histo ,
                                       const FrequencyDistribution **ihisto , const char *title) const

{
  bool status = true;


  error.init();

  if (nb_histo > PLOT_NB_HISTOGRAM) {
    status = false;
    error.update(STAT_error[STATR_PLOT_NB_HISTOGRAM]);
  }

  else {
    bool cumul_concentration_flag;
    int i , j , k;
    int max_nb_value , max_frequency , max_range , reference_matching ,
        reference_concentration , *poffset , *pnb_value;
    double max_mass , shift , **cumul , **concentration;
    const FrequencyDistribution **histo , *phisto[2] , **merged_histo;
    ostringstream *data_file_name;


    nb_histo++;
    histo = new const FrequencyDistribution*[nb_histo];

    histo[0] = this;
    for (i = 1;i < nb_histo;i++) {
      histo[i] = ihisto[i - 1];
    }

    // writing of the data files

    data_file_name = new ostringstream[nb_histo + 2];

    if (nb_histo > 1) {
      data_file_name[0] << prefix << 0 << ".dat";

      merged_histo = new const FrequencyDistribution*[nb_histo];
      merged_histo[nb_histo - 1] = new FrequencyDistribution(*histo[nb_histo - 1]);

      for (i = nb_histo - 2;i >= 0;i--) {
        phisto[0] = merged_histo[i + 1];
        phisto[1] = histo[i];
        merged_histo[i] = new FrequencyDistribution(2 , phisto);
      }

      status = merged_histo[0]->plot_print((data_file_name[0].str()).c_str() ,
                                            nb_histo - 1 , merged_histo + 1);
    }

    if (status) {
      poffset = new int[nb_histo];
      pnb_value = new int[nb_histo];

      cumul = new double*[nb_histo];
      concentration = new double*[nb_histo];
      for (i = 0;i < nb_histo;i++) {
        cumul[i] = NULL;
        concentration[i] = NULL;
      }

      max_nb_value = 0;
      max_frequency = 0;
      max_mass = 0.;
      cumul_concentration_flag = false;
      max_range = 0;
      shift = 0.;

      for (i = 0;i < nb_histo;i++) {
        data_file_name[i + 1] << prefix << i + 1 << ".dat";

        poffset[i] = histo[i]->offset;
        pnb_value[i] = histo[i]->nb_value;

        // computation of the cumulative distribution functions and concentration curves

        cumul[i] = histo[i]->cumul_computation();
        concentration[i] = histo[i]->concentration_function_computation();

        // computation of the maximum number of values, the largest support,
        // the maximum frequency and the maximum probability

        if (histo[i]->nb_value > max_nb_value) {
          max_nb_value = histo[i]->nb_value;
        }
        if (histo[i]->max > max_frequency) {
          max_frequency = histo[i]->max;
        }
        if ((double)histo[i]->max / (double)histo[i]->nb_element > max_mass) {
          max_mass = (double)histo[i]->max / (double)histo[i]->nb_element;
        }

        if (((histo[i]->variance == D_DEFAULT) && (histo[i]->nb_value > histo[i]->offset + 1)) ||
            (histo[i]->variance > 0.)) {
          cumul_concentration_flag = true;
          if (histo[i]->nb_value - histo[i]->offset > max_range) {
            max_range = histo[i]->nb_value - histo[i]->offset;
            reference_matching = i;
          }
          reference_concentration = i;
        }

        status = histo[i]->plot_print((data_file_name[i + 1].str()).c_str() , cumul[i] ,
                                      concentration[i] , shift);

        if (!status) {
          break;
        }

        if (nb_histo > 1) {
          if (PLOT_SHIFT * (nb_histo - 1) < PLOT_MAX_SHIFT) {
            shift += PLOT_SHIFT;
          }
          else {
            shift += PLOT_MAX_SHIFT / (nb_histo - 1);
          }
        }
      }

      if ((cumul_concentration_flag) && (nb_histo > 1)) {
        data_file_name[nb_histo + 1] << prefix << nb_histo + 1 << ".dat";
        cumul_matching_plot_print((data_file_name[nb_histo + 1].str()).c_str() , nb_histo ,
                                  poffset , pnb_value , cumul);
      }

      delete [] poffset;
      delete [] pnb_value;

      for (i = 0;i < nb_histo;i++) {
        delete [] cumul[i];
        delete [] concentration[i];
      }
      delete [] cumul;
      delete [] concentration;
    }

    if (status) {

      // writing of the script files

      for (i = 0;i < 2;i++) {
        ostringstream file_name[2];

        switch (i) {
        case 0 :
          file_name[0] << prefix << ".plot";
          break;
        case 1 :
          file_name[0] << prefix << ".print";
          break;
        }

        ofstream out_file((file_name[0].str()).c_str());

        if (i == 1) {
          out_file << "set terminal postscript" << endl;
          file_name[1] << label(prefix) << ".ps";
          out_file << "set output \"" << file_name[1].str() << "\"\n\n";
        }

        out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                 << "set title";
        if (title) {
          out_file << " \"" << title << "\"";
        }
        out_file << "\n\n";

        // shifted frequency distributions

        if (MAX(max_nb_value , 2) < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }
        if ((int)(max_frequency * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [0:" << MAX(max_nb_value , 2) << "] [0:"
                 << (int)(max_frequency * YSCALE) + 1 << "] ";
        for (j = 0;j < nb_histo;j++) {
          out_file << "\"" << label((data_file_name[j + 1].str()).c_str()) << "\" using 2:3 title \""
                   << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
          if (nb_histo > 1) {
            out_file << " " << j + 1;
          }
          out_file << "\" with impulses";
          if (j < nb_histo - 1) {
            out_file << ",\\";
          }
          out_file << endl;
        }

        if (MAX(max_nb_value , 2) < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
        if ((int)(max_frequency * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics autofreq" << endl;
        }

        if (nb_histo > 1) {

          // cumulative frequencies

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;

          if (MAX(merged_histo[0]->nb_value , 2) - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }
          if ((int)(merged_histo[0]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [0:" << MAX(merged_histo[0]->nb_value , 2) - 1 << "] [0:"
                   << (int)(merged_histo[0]->max * YSCALE) + 1 << "] ";
          for (j = 0;j < nb_histo;j++) {
            out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                     << j + 1 << " title \"" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << j + 1
                     << "\" with impulses";
            if (j < nb_histo - 1) {
              out_file << ",\\";
            }
            out_file << endl;
          }

          if (MAX(merged_histo[0]->nb_value , 2) - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
          if ((int)(merged_histo[0]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics autofreq" << endl;
          }

          // probability mass functions

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;

          if (MAX(max_nb_value , 2) - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }

          out_file << "plot [0:" << MAX(max_nb_value , 2) - 1 << "] [0:"
                   << MIN(max_mass * YSCALE , 1.) << "] ";
          for (j = 0;j < nb_histo;j++) {
            out_file << "\"" << label((data_file_name[j + 1].str()).c_str()) << "\" using 1:4 title \""
                     << STAT_label[STATL_DISTRIBUTION] << " " << j + 1 << "\" with linespoints";
            if (j < nb_histo - 1) {
              out_file << ",\\";
            }
            out_file << endl;
          }

          if (MAX(max_nb_value , 2) - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
        }

        if (cumul_concentration_flag) {

          // cumulative distribution functions

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;

          if (max_nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }

          out_file << "plot [0:" << max_nb_value - 1 << "] [0:1] ";
          j = 0;
          for (k = 0;k < nb_histo;k++) {
            if (((histo[k]->variance == D_DEFAULT) && (histo[k]->nb_value > histo[k]->offset + 1)) ||
                (histo[k]->variance > 0.)) {
              if (j > 0) {
                out_file << ",\\\n";
              }
              j++;
              out_file << "\"" << label((data_file_name[k + 1].str()).c_str()) << "\" using 1:5 title \""
                       << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION];
              if (nb_histo > 1) {
                out_file << " " << k + 1;
              }
              out_file << " " << STAT_label[STATL_FUNCTION] << "\" with linespoints";
            }
          }
          out_file << endl;

          if (max_nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }

          // matching of cumulative distribution functions taking as reference
          // the distribution with the largest support

          if (nb_histo > 1) {
            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set title";
            out_file << " \"";
            if (title) {
              out_file << title << " - ";
            }
            out_file << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION]
                     << " " << STAT_label[STATL_FUNCTION] << " " << STAT_label[STATL_MATCHING] << "\"\n\n";

            out_file << "set grid\n" << "set xtics 0,0.1\n" << "set ytics 0,0.1" << endl;

            out_file << "plot [0:1] [0:1] ";
            j = 0;
            for (k = 0;k < nb_histo;k++) {
              if (((histo[k]->variance == D_DEFAULT) && (histo[k]->nb_value > histo[k]->offset + 1)) ||
                  (histo[k]->variance > 0.)) {
                if (j > 0) {
                  out_file << ",\\\n";
                }
                j++;
                out_file << "\"" << label((data_file_name[nb_histo + 1].str()).c_str()) << "\" using "
                         << reference_matching + 1 << ":" << k + 1 << " title \""
                         << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
                         << k + 1 << " " << STAT_label[STATL_FUNCTION] << "\" with linespoints";
              }
            }
            out_file << endl;

            out_file << "unset grid\n" << "set xtics autofreq\n" << "set ytics autofreq" << endl;
          }

          // concentration curves

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;

          out_file << "set title";
          if (title) {
            out_file << " \"" << title << "\"";
          }
          out_file << "\n\n";

          out_file << "set grid\n" << "set xtics 0,0.1\n" << "set ytics 0,0.1" << endl;

          out_file << "plot [0:1] [0:1] ";
          for (j = 0;j < nb_histo;j++) {
            if (((histo[j]->variance == D_DEFAULT) && (histo[j]->nb_value > histo[j]->offset + 1)) ||
                (histo[j]->variance > 0.)) {
              out_file << "\"" << label((data_file_name[j + 1].str()).c_str()) << "\" using 5:6 title \""
                       << STAT_label[STATL_CONCENTRATION] << " " << STAT_label[STATL_CURVE];
              if (nb_histo > 1) {
                out_file << " " << j + 1;
              }
              out_file << "\" with linespoints,\\" << endl;
            }
          }
          out_file << "\"" << label((data_file_name[reference_concentration + 1].str()).c_str())
                   << "\" using 5:5 notitle with lines" << endl;

          out_file << "unset grid\n" << "set xtics autofreq\n" << "set ytics autofreq" << endl;
        }

        if (i == 1) {
          out_file << "\nset terminal x11" << endl;
        }

        out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
    }

    if (nb_histo > 1) {
      for (i = 0;i < nb_histo;i++) {
        delete merged_histo[i];
      }
      delete [] merged_histo;
    }

    delete [] histo;
    delete [] data_file_name;

    if (!status) {
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a family of frequency distributions:
 *         - frequency distributions and cumulative frequencies,
 *         - probability mass functions and cumulative distribution functions,
 *         - matching of cumulative distribution functions,
 *         - concentration curves.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] nb_histo number of frequency distributions,
 *  \param[in] ihisto   pointer on the frequency distributions.
 *
 *  \return             plots.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* FrequencyDistribution::get_plotable_frequency_distributions(StatError &error , int nb_histo ,
                                                                          const FrequencyDistribution **ihisto) const

{
  MultiPlotSet *plot_set;


  error.init();

  if (nb_histo > PLOT_NB_HISTOGRAM) {
    plot_set = NULL;
    error.update(STAT_error[STATR_PLOT_NB_HISTOGRAM]);
  }

  else {
    int i , j , k;
    int max_nb_value , max_frequency , cumul_concentration_nb_histo , max_range ,
        reference_matching , nb_plot_set;
    double max_mass , shift , **cumul;
    const FrequencyDistribution **histo , *phisto[2] , **merged_histo;
    ostringstream legend , title;


    nb_histo++;
    histo = new const FrequencyDistribution*[nb_histo];

    histo[0] = this;
    for (i = 1;i < nb_histo;i++) {
      histo[i] = ihisto[i - 1];
    }

    cumul = new double*[nb_histo];

    max_nb_value = 0;
    max_frequency = 0;
    max_mass = 0.;
    cumul_concentration_nb_histo = 0;
    max_range = 0;

    for (i = 0;i < nb_histo;i++) {

      // computation of the cumulative distribution functions

      cumul[i] = histo[i]->cumul_computation();

      // computation of the maximum number of values, the largest support,
      // the maximum frequency and the maximum probability

      if (histo[i]->nb_value > max_nb_value) {
        max_nb_value = histo[i]->nb_value;
      }
      if (histo[i]->max > max_frequency) {
        max_frequency = histo[i]->max;
      }
      if ((double)histo[i]->max / (double)histo[i]->nb_element > max_mass) {
        max_mass = (double)histo[i]->max / (double)histo[i]->nb_element;
      }

      if (((histo[i]->variance == D_DEFAULT) && (histo[i]->nb_value > histo[i]->offset + 1)) ||
          (histo[i]->variance > 0.)) {
        cumul_concentration_nb_histo++;
        if (histo[i]->nb_value - histo[i]->offset > max_range) {
          max_range = histo[i]->nb_value - histo[i]->offset;
          reference_matching = i;
        }
      }
    }

    nb_plot_set = 3;
    if (nb_histo == 1) {
      nb_plot_set -= 2;
    }
    if (cumul_concentration_nb_histo > 0) {
      nb_plot_set += 3;
      if (nb_histo == 1) {
        nb_plot_set -= 1;
      }
    }

    // number of views

    plot_set = new MultiPlotSet(nb_plot_set);
    MultiPlotSet &plot = *plot_set;

    title.str("");
    title << STAT_label[STATL_FREQUENCY];
    if (nb_histo == 1) {
      title << " " << STAT_label[STATL_DISTRIBUTION];
    }
    else {
      title << " " << STAT_label[STATL_DISTRIBUTIONS];
    }
    plot.title = title.str();

    plot.border = "15 lw 0";

    // shifted frequency distributions

    plot[0].xrange = Range(0 , MAX(max_nb_value , 2));
    plot[0].yrange = Range(0 , ceil(max_frequency * YSCALE));

    if (MAX(max_nb_value , 2) < TIC_THRESHOLD) {
      plot[0].xtics = 1;
    }
    if (ceil(max_frequency * YSCALE) < TIC_THRESHOLD) {
      plot[0].ytics = 1;
    }

    plot[0].resize(nb_histo);

    shift = 0.;

    for (i = 0;i < nb_histo;i++) {
      legend.str("");
      legend << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i + 1;
      plot[0][i].legend = legend.str();

      plot[0][i].style = "impulses";

      for (j = histo[i]->offset;j < histo[i]->nb_value;j++) {
        if (histo[i]->frequency[j] > 0) {
          plot[0][i].add_point(j + shift , histo[i]->frequency[j]);
        }
      }

      if (PLOT_SHIFT * (nb_histo - 1) < PLOT_MAX_SHIFT) {
        shift += PLOT_SHIFT;
      }
      else {
        shift += PLOT_MAX_SHIFT / (nb_histo - 1);
      }
    }

    if (nb_histo > 1) {

      // cumulative frequencies

      merged_histo = new const FrequencyDistribution*[nb_histo];
      merged_histo[nb_histo - 1] = new FrequencyDistribution(*histo[nb_histo - 1]);

      for (i = nb_histo - 2;i >= 0;i--) {
        phisto[0] = merged_histo[i + 1];
        phisto[1] = histo[i];
        merged_histo[i] = new FrequencyDistribution(2 , phisto);
      }

      plot[1].xrange = Range(0 , MAX(merged_histo[0]->nb_value , 2) - 1);
      plot[1].yrange = Range(0 , ceil(merged_histo[0]->max * YSCALE));

      if (MAX(merged_histo[0]->nb_value , 2) - 1 < TIC_THRESHOLD) {
        plot[1].xtics = 1;
      }
      if (ceil(merged_histo[0]->max * YSCALE) < TIC_THRESHOLD) {
        plot[1].ytics = 1;
      }

      plot[1].resize(nb_histo);

      for (i = 0;i < nb_histo;i++) {
        legend.str("");
        legend << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i + 1;
        plot[1][i].legend = legend.str();

        plot[1][i].style = "impulses";

        merged_histo[i]->plotable_frequency_write(plot[1][i]);
      }

      for (i = 0;i < nb_histo;i++) {
        delete merged_histo[i];
      }
      delete [] merged_histo;

      // probability mass functions

      plot[2].xrange = Range(0 , MAX(max_nb_value , 2) - 1);
      plot[2].yrange = Range(0 , MIN(max_mass * YSCALE , 1.));

      if (MAX(max_nb_value , 2) - 1 < TIC_THRESHOLD) {
        plot[2].xtics = 1;
      }

      plot[2].resize(nb_histo);

      for (i = 0;i < nb_histo;i++) {
        legend.str("");
        legend << STAT_label[STATL_DISTRIBUTION] << " " << i + 1;
        plot[2][i].legend = legend.str();

        plot[2][i].style = "linespoints";

        histo[i]->plotable_mass_write(plot[2][i]);
      }

      i = 3;
    }

    else {
      i = 1;
    }

    if (cumul_concentration_nb_histo > 0) {

      // cumulative distribution functions

      plot[i].xrange = Range(0 , max_nb_value - 1);
      plot[i].yrange = Range(0. , 1.);

      if (max_nb_value - 1 < TIC_THRESHOLD) {
        plot[i].xtics = 1;
      }

      plot[i].resize(cumul_concentration_nb_histo);

      j = 0;
      for (k = 0;k < nb_histo;k++) {
        if (((histo[k]->variance == D_DEFAULT) && (histo[k]->nb_value > histo[k]->offset + 1)) ||
            (histo[k]->variance > 0.)) {
          legend.str("");
          legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION];
          if (nb_histo > 1) {
            legend << " " << k + 1;
          }
          legend << " " << STAT_label[STATL_FUNCTION];
          plot[i][j].legend = legend.str();

          plot[i][j].style = "linespoints";

          histo[k]->plotable_cumul_write(plot[i][j] , cumul[k]);
          j++;
        }
      }

      i++;

      // matching of cumulative distribution functions taking as reference
      // the distribution with the largest support

      if (nb_histo > 1) {
        title.str("");
        title << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION]
              << " " << STAT_label[STATL_FUNCTION] << " " << STAT_label[STATL_MATCHING];
        plot[i].title = title.str();

        plot[i].xrange = Range(0. , 1.);
        plot[i].yrange = Range(0. , 1.);

        plot[i].grid = true;

        plot[i].xtics = 0.1;
        plot[i].ytics = 0.1;

        plot[i].resize(cumul_concentration_nb_histo);

        j = 0;
        for (k = 0;k < nb_histo;k++) {
          if (((histo[k]->variance == D_DEFAULT) && (histo[k]->nb_value > histo[k]->offset + 1)) ||
              (histo[k]->variance > 0.)) {
            legend.str("");
            legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION]
                   << " " << k + 1 << " " << STAT_label[STATL_FUNCTION];
            plot[i][j].legend = legend.str();

            plot[i][j].style = "linespoints";

            histo[k]->plotable_cumul_matching_write(plot[i][j] , histo[reference_matching]->offset ,
                                                    histo[reference_matching]->nb_value ,
                                                    cumul[reference_matching] , cumul[k]);
            j++;
          }
        }

        i++;
      }

      // concentration curves

      plot[i].xrange = Range(0. , 1.);
      plot[i].yrange = Range(0. , 1.);

      plot[i].grid = true;

      plot[i].xtics = 0.1;
      plot[i].ytics = 0.1;

      plot[i].resize(cumul_concentration_nb_histo + 1);

      j = 0;
      for (k = 0;k < nb_histo;k++) {
        if (((histo[k]->variance == D_DEFAULT) && (histo[k]->nb_value > histo[k]->offset + 1)) ||
            (histo[k]->variance > 0.)) {
          legend.str("");
          legend << STAT_label[STATL_CONCENTRATION] << " " << STAT_label[STATL_CURVE];
          if (nb_histo > 1) {
            legend << " " << k + 1;
          }
          plot[i][j].legend = legend.str();

          plot[i][j].style = "linespoints";

          histo[k]->plotable_concentration_write(plot[i][j] , cumul[k]);
          j++;
        }
      }

      plot[i][j].style = "lines";

      plot[i][j].add_point(0. , 0.);
      plot[i][j].add_point(1. , 1.);
    }

    for (i = 0;i < nb_histo;i++) {
      delete [] cumul[i];
    }
    delete [] cumul;

    delete [] histo;
  }

  return plot_set;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a frequency distribution.
 *
 *  \return MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* FrequencyDistribution::get_plotable() const

{
  StatError error;

  return get_plotable_frequency_distributions(error , 0 , NULL);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a DiscreteDistributionData object from
 *         a FrequencyDistribution and a Distribution objects.
 *
 *  \param[in] histo reference on a FrequencyDistribution object,
 *  \param[in] dist  pointer on a Distribution object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData::DiscreteDistributionData(const FrequencyDistribution &histo ,
                                                   const Distribution *dist)
:FrequencyDistribution(histo)

{
  if (dist) {
    distribution = new DiscreteParametricModel(*dist);
  }
  else {
    distribution = NULL;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a DiscreteDistributionData object from
 *         a FrequencyDistribution and a DiscreteParametric objects.
 *
 *  \param[in] histo reference on a FrequencyDistribution object,
 *  \param[in] dist  pointer on a DiscreteParametric object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData::DiscreteDistributionData(const FrequencyDistribution &histo ,
                                                   const DiscreteParametric *dist)
:FrequencyDistribution(histo)

{
  if (dist) {
    distribution = new DiscreteParametricModel(*dist);
  }
  else {
    distribution = NULL;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor by copy of the DiscreteDistributionData class.
 *
 *  \param[in] histo      reference on a DiscreteDistributionData object,
 *  \param[in] model_flag flag copy of the included DiscreteParametricModel object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData::DiscreteDistributionData(const DiscreteDistributionData &histo ,
                                                   bool model_flag)
:FrequencyDistribution(histo)

{
  if ((model_flag) && (histo.distribution)) {
    distribution = new DiscreteParametricModel(*(histo.distribution) , false);
  }
  else {
    distribution = NULL;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the DiscreteDistributionData class.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData::~DiscreteDistributionData()

{
  delete distribution;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the DiscreteDistributionData class.
 *
 *  \param[in] histo reference on a DiscreteDistributionData object.
 *
 *  \return          DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData& DiscreteDistributionData::operator=(const DiscreteDistributionData &histo)

{
  if (&histo != this) {
    delete distribution;
    delete [] frequency;

    FrequencyDistribution::copy(histo);
    if (histo.distribution) {
      distribution = new DiscreteParametricModel(*(histo.distribution) , false);
    }
    else {
      distribution = NULL;
    }
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of the DiscreteParametricModel object included in
 *         a DiscreteDistributionData object.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          DiscreteParametricModel object.
 */
/*--------------------------------------------------------------*/

DiscreteParametricModel* DiscreteDistributionData::extract_model(StatError &error) const

{
  DiscreteParametricModel *dist;


  error.init();

  if (!distribution) {
    dist = NULL;
    error.update(STAT_error[STATR_NON_EXISTING_DISTRIBUTION]);
  }

  else {
    dist = new DiscreteParametricModel(*distribution);
    dist->frequency_distribution = new DiscreteDistributionData(*this , false);
  }

  return dist;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a FrequencyDistribution object from a file.
 *         Format: n rows with ordered value and associated frequency.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData* DiscreteDistributionData::ascii_read(StatError &error , const string path)

{
  string buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  bool status , lstatus;
  int i;
  int line , nb_element , value , index , max_index;
  DiscreteDistributionData *histo;
  ifstream in_file(path.c_str());


  histo = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {

    // 1st pass: search for the number of possible values and format checking

    status = true;
    line = 0;
    max_index = -1;
    nb_element = 0;

    while (getline(in_file , buffer)) {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.find('#');
      if (position != string::npos) {
        buffer.erase(position);
      }
      i = 0;

      tokenizer tok_buffer(buffer , separator);

      for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
        if (i <= 1) {
          lstatus = true;

/*          try {
            value = stoi(*token);   in C++ 11
          }
          catch(invalid_argument &arg) {
            lstatus = false;
          } */
          value = atoi(token->c_str());

          if ((lstatus) && (value < 0)) {
            lstatus = false;
          }

          if (!lstatus) {
            status = false;
            error.update(STAT_parsing[STATP_DATA_TYPE] , line , i + 1);
          }

          else {
            switch (i) {

            // test ordered values

            case 0 : {
              if (value <= max_index) {
                status = false;
                error.update(STAT_parsing[STATP_VALUE_ORDER] , line , i + 1);
              }
              else {
                max_index = value;
              }

//              if (value >= SAMPLE_NB_VALUE) {
              if (value >= SAMPLE_NB_VALUE * 5) {
                status = false;
                error.update(STAT_parsing[STATP_MAX_VALUE] , line , i + 1);
              }
              break;
            }

            case 1 : {
              if (value > 0) {
                nb_element += value;
              }
              break;
            }
            }
          }
        }

        i++;
      }

      // test 2 items per line

      if ((i > 0) && (i != 2)) {
        status = false;
        error.correction_update(STAT_parsing[STATP_NB_TOKEN] , 2 , line);
      }
    }

    if (nb_element == 0) {
      status = false;
      error.update(STAT_parsing[STATP_EMPTY_SAMPLE]);
    }

    // 2nd pass: file reading

    if (status) {
//      in_file.close();
//      in_file.open(path.c_str() , ios::in);

      in_file.clear();
      in_file.seekg(0 , ios::beg);

      histo = new DiscreteDistributionData(max_index + 1);

      while (getline(in_file , buffer)) {
        position = buffer.find('#');
        if (position != string::npos) {
          buffer.erase(position);
        }
        i = 0;

        tokenizer tok_buffer(buffer , separator);

        for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
          switch (i) {
          case 0 :
//            index = stoi(*token);   in C++ 11
            index = atoi(token->c_str());
            break;
          case 1 :
//            histo->frequency[index] = stoi(*token);   in C++ 11
            histo->frequency[index] = atoi(token->c_str());
            break;
          }

          i++;
        }

        if ((i > 0) && (index == max_index)) {
          break;
        }
      }

      histo->nb_value_computation();
      histo->offset_computation();
      histo->nb_element = nb_element;
      histo->max_computation();
      histo->mean_computation();
      histo->variance_computation();
    }
  }

  return histo;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing on a single line of a DiscreteDistributionData object.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& DiscreteDistributionData::line_write(ostream &os) const

{
  os << STAT_label[STATL_SAMPLE_SIZE] << ": " << nb_element;

  if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
    os << "   " << STAT_label[STATL_MEAN] << ": " << mean
       << "   " << STAT_label[STATL_VARIANCE] << ": " << variance
       << "   " << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance);
   }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a DiscreteDistributionData object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level,
 *  \param[in]     file_flag  flag file.
 */
/*--------------------------------------------------------------*/

ostream& DiscreteDistributionData::ascii_write(ostream &os , bool exhaustive ,
                                               bool file_flag) const

{
  if (distribution) {
    distribution->ascii_write(os , this , exhaustive , file_flag);
  }
  else {
    FrequencyDistribution::ascii_write(os , exhaustive , file_flag);
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a DiscreteDistributionData object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& DiscreteDistributionData::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , exhaustive , false);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a DiscreteDistributionData object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool DiscreteDistributionData::ascii_write(StatError &error , const string path ,
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
    ascii_write(out_file , exhaustive , true);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a DiscreteDistributionData object in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool DiscreteDistributionData::spreadsheet_write(StatError &error , const string path) const

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

    if (distribution) {
      distribution->spreadsheet_write(out_file , this);
    }

    else {
      double information = information_computation();


      out_file << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      spreadsheet_characteristic_print(out_file , true);

      out_file << STAT_label[STATL_MEAN_ABSOLUTE_DEVIATION] << "\t" << mean_absolute_deviation_computation(mean);
      if (mean > 0.) {
        out_file << "\t\t" << STAT_label[STATL_CONCENTRATION_COEFF] << "\t" << concentration_computation();
      }
      out_file << endl;

      out_file << STAT_label[STATL_INFORMATION] << "\t" << information << "\t"
               << information / nb_element << endl;

      out_file << "\n\t" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t" << STAT_label[STATL_CUMULATIVE]
               << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION];
      if (mean > 0.) {
        out_file << "\t" << STAT_label[STATL_CONCENTRATION] << " " << STAT_label[STATL_FUNCTION];
      }
      out_file << endl;
      spreadsheet_print(out_file , true , (mean == 0. ? false : true));
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a DiscreteDistributionData object using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool DiscreteDistributionData::plot_write(StatError &error , const char *prefix ,
                                          const char *title) const

{
  bool status;


  if (distribution) {
    status = distribution->plot_write(error , prefix , title , this);
  }

  else {
    status = FrequencyDistribution::plot_write(error , prefix , 0 , NULL , title);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a DiscreteDistributionData object.
 *
 *  \return MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* DiscreteDistributionData::get_plotable() const

{
  MultiPlotSet *plot_set;


  if (distribution) {
    plot_set = distribution->get_plotable(this);
  }
  else {
    StatError error;
    plot_set = FrequencyDistribution::get_plotable_frequency_distributions(error , 0 , NULL);
  }

  return plot_set;
}


};  // namespace stat_tool
