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
 *       $Id: frequency_distribution2.cpp 18002 2015-04-23 06:57:18Z guedon $
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

#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"
#include "tool/config.h"

#include "stat_tools.h"
#include "distribution.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {


extern bool cumul_matching_plot_print(const char *path , int nb_cumul , int *offset ,
                                      int *nb_value , double **cumul);



/*--------------------------------------------------------------*
 *
 *  Ecriture des resultats d'une comparaison de lois empiriques.
 *
 *  arguments : stream, nombre de lois empiriques, pointeurs sur les lois empiriques,
 *              type de variable (SYMBOLIC/ORDINAL/NUMERIC), dissimilarites.
 *
 *--------------------------------------------------------------*/

ostream& FrequencyDistribution::dissimilarity_ascii_write(ostream &os , int nb_histo ,
                                                          const FrequencyDistribution **ihisto ,
                                                          int type , double **dissimilarity) const

{
  register int i , j;
  int max_nb_value , buff , width[3];
  long old_adjust;
  double information , **cumul;
  Test *test;
  const FrequencyDistribution **histo;


  nb_histo++;
  histo = new const FrequencyDistribution*[nb_histo];

  histo[0] = this;
  for (i = 1;i < nb_histo;i++) {
    histo[i] = ihisto[i - 1];
  }

  old_adjust = os.setf(ios::right , ios::adjustfield);

  // ecriture des caracteristiques des lois empiriques

  for (i = 0;i < nb_histo;i++) {
    os << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i + 1  << " - ";

    if (type == SYMBOLIC) {
      os << STAT_label[STATL_SAMPLE_SIZE] << ": " << histo[i]->nb_element << endl;
    }

    else {
      histo[i]->ascii_characteristic_print(os , true);

      os << STAT_label[STATL_MEAN_ABSOLUTE_DEVIATION] << ": " << histo[i]->mean_absolute_deviation_computation();
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

  // calcul des largeurs des colonnes

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

  // ecriture des lois empiriques et des fonctions de repartition

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

  // calcul des largeurs des colonnes

  width[0] = column_width(nb_histo);

  width[1] = 0;
  for (i = 0;i < nb_histo;i++) {
    buff = column_width(nb_histo , dissimilarity[i]);
    if (buff > width[1]) {
      width[1] = buff;
    }
  }
  width[1] += ASCII_SPACE;

  // ecriture de la matrice des dissimilarites entre lois empiriques

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

  // analyse de variance

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

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  delete [] histo;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des resultats d'une comparaison de lois empiriques dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path, nombre de lois empiriques,
 *              pointeurs sur les lois empiriques, type de variable (SYMBOLIC/ORDINAL/NUMERIC),
 *              dissimilarites.
 *
 *--------------------------------------------------------------*/

bool FrequencyDistribution::dissimilarity_ascii_write(StatError &error , const char *path ,
                                                      int nb_histo , const FrequencyDistribution **ihisto ,
                                                      int type , double **dissimilarity) const

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
    dissimilarity_ascii_write(out_file , nb_histo , ihisto , type , dissimilarity);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des resultats d'une comparaison de lois empiriques dans un fichier
 *  au format tableur.
 *
 *  arguments : reference sur un objet StatError, path, nombre de lois empiriques,
 *              pointeurs sur les lois empiriques, type de variable (SYMBOLIC/ORDINAL/NUMERIC),
 *              dissimilarites.
 *
 *--------------------------------------------------------------*/

bool FrequencyDistribution::dissimilarity_spreadsheet_write(StatError &error , const char *path ,
                                                            int nb_histo , const FrequencyDistribution **ihisto ,
                                                            int type , double **dissimilarity) const

{
  bool status;
  register int i , j;
  int max_nb_value;
  double information , **cumul , **concentration;
  Test *test;
  const FrequencyDistribution **histo;
  ofstream out_file(path);


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

    // ecriture des caracteristiques des lois empiriques

    for (i = 0;i < nb_histo;i++) {
      out_file << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i + 1  << "\t";

      if (type == SYMBOLIC) {
        out_file << STAT_label[STATL_SAMPLE_SIZE] << "\t" << histo[i]->nb_element << endl;
      }

      else {
        histo[i]->spreadsheet_characteristic_print(out_file , true);

        out_file << STAT_label[STATL_MEAN_ABSOLUTE_DEVIATION] << "\t" << histo[i]->mean_absolute_deviation_computation();
        if (histo[i]->mean > 0.) {
          out_file << "\t\t" << STAT_label[STATL_CONCENTRATION_COEFF] << "\t" << histo[i]->concentration_computation();
        }
        out_file << endl;
      }

      information = histo[i]->information_computation();

      out_file << STAT_label[STATL_INFORMATION] << "\t" << information << "\t"
               << information / histo[i]->nb_element << "\n" << endl;
    }

    // ecriture des lois empiriques, des fonctions de repartition et des courbes de concentration

    cumul = new double*[nb_histo];
    for (i = 0;i < nb_histo;i++) {
      cumul[i] = histo[i]->cumul_computation();
    }

    if (type != SYMBOLIC) {
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

    if (type != SYMBOLIC) {
      for (i = 0;i < nb_histo;i++) {
        if (histo[i]->mean > 0.) {
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

      if (type != SYMBOLIC) {
        for (j = 0;j < nb_histo;j++) {
          if (histo[j]->mean > 0.) {
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

    if (type != SYMBOLIC) {
      for (i = 0;i < nb_histo;i++) {
        delete [] concentration[i];
      }
      delete [] concentration;
    }

    // ecriture de la matrice des dissimilarites entre lois empiriques

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

    // analyse de variance

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


/*--------------------------------------------------------------*
 *
 *  Test de Kruskal-Wallis (analyse de variance dans le cas ordinal).
 *
 *  arguments : nombre de lois empiriques, pointeurs sur les lois empiriques.
 *
 *--------------------------------------------------------------*/

Test* FrequencyDistribution::kruskal_wallis_test(int nb_histo , const FrequencyDistribution **ihisto) const

{
  register int i , j;
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

  // calcul du terme de correction pour les ex-aequo

  pfrequency = merged_histo->frequency + merged_histo->offset;
  correction = 0.;
  for (i = merged_histo->offset;i < merged_histo->nb_value;i++) {
    if (*pfrequency > 1) {
      correction += *pfrequency * ((double)*pfrequency * (double)*pfrequency - 1);
    }
    pfrequency++;
  }

  // calcul des rangs

  rank = merged_histo->rank_computation();

  // calcul de la statistique de Kruskal-Wallis

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


/*--------------------------------------------------------------*
 *
 *  Comparaison de lois empiriques.
 *
 *  arguments : reference sur un objet StatError, stream, nombre de lois empiriques,
 *              pointeurs sur les lois empiriques, type de variable (SYMBOLIC/ORDINAL/NUMERIC),
 *              path, format de fichier ('a' : ASCII, 's' : Spreadsheet).
 *
 *--------------------------------------------------------------*/

bool FrequencyDistribution::comparison(StatError &error , ostream &os , int nb_histo ,
                                       const FrequencyDistribution **ihisto , int type ,
                                       const char *path , char format) const

{
  bool status = true;
  register int i , j , k , m;
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

    // calcul des rangs

    rank = merged_histo->rank_computation();

    // calcul des rangs moyens pour chaque loi empirique

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

    // calcul des largeurs des colonnes

    width[0] = column_width(nb_histo);

    width[1] = 0;
    for (i = 0;i < nb_histo;i++) {
      buff = column_width(nb_histo , dissimilarity[i]);
      if (buff > width[1]) {
        width[1] = buff;
      }
    }
    width[1] += ASCII_SPACE;

    // ecriture de la matrice des dissimilarites entre lois empiriques

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
    os << endl;

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
          if (type == SYMBOLIC) {
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

      if (type == SYMBOLIC) {
        dissimilarity[j][i] = dissimilarity[i][j];
      }
      else {
        dissimilarity[j][i] = -dissimilarity[i][j];
      }
    }
  }

# ifdef MESSAGE
  dissimilarity_ascii_write(os , nb_histo - 1 , ihisto , type , dissimilarity);
# endif

  if (path) {
    switch (format) {
    case 'a' :
      status = dissimilarity_ascii_write(error , path , nb_histo - 1 , ihisto ,
                                         type , dissimilarity);
      break;
    case 's' :
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


/*--------------------------------------------------------------*
 *
 *  Test F de comparaison des variances de 2 lois empiriques.
 *
 *  arguments : stream, reference sur une loi empirique.
 *
 *--------------------------------------------------------------*/

void FrequencyDistribution::F_comparison(ostream &os , const FrequencyDistribution &histo) const

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

#   ifdef MESSAGE
    os << *test;
#   endif

    delete test;
  }
}


/*--------------------------------------------------------------*
 *
 *  Test t de comparaison des moyennes de 2 lois empiriques.
 *
 *  arguments : stream, reference sur une loi empirique.
 *
 *--------------------------------------------------------------*/

void FrequencyDistribution::t_comparison(ostream &os , const FrequencyDistribution &histo) const

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

#   ifdef MESSAGE
    os << *test;
#   endif

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

    os << "\n" << *test;

    delete test;
  }
# endif

}


/*--------------------------------------------------------------*
 *
 *  Test de comparaison de Wilcoxon-Mann-Whitney.
 *
 *  arguments : reference sur un objet StatError, stream,
 *              reference sur une loi empirique.
 *
 *--------------------------------------------------------------*/

bool FrequencyDistribution::wilcoxon_mann_whitney_comparison(StatError &error , ostream &os ,
                                                             const FrequencyDistribution &ihisto) const

{
  bool status;
  register int i;
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

    // calcul du terme de correction pour les ex-aequo

    pfrequency = merged_histo->frequency + merged_histo->offset;
    correction = 0.;
    for (i = merged_histo->offset;i < merged_histo->nb_value;i++) {
      if (*pfrequency > 1) {
        correction += *pfrequency * ((double)*pfrequency * (double)*pfrequency - 1);
      }
      pfrequency++;
    }

    // calcul des rangs

    rank = merged_histo->rank_computation();

    // calcul de la statistique de Wilcoxon

    value = 0.;
    for (i = histo[0]->offset;i < histo[0]->nb_value;i++) {
      value += histo[0]->frequency[i] * rank[i];
    }

    nb_sup = value - histo[0]->nb_element * ((double)histo[0]->nb_element + 1) / 2.;

    mean = histo[0]->nb_element * ((double)merged_histo->nb_element + 1) / 2.;

    // correction de continuite

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

#   ifdef MESSAGE
    os << STAT_label[STATL_TWO_SIDED] << " " << STAT_label[STATL_WILCOXON_MANN_WHITNEY_TEST];
    os << *test;

    os << STAT_label[STATL_MANN_WHITNEY_INFERIOR_PROBABILITY] << " = "
       << (histo[0]->nb_element * (double)histo[1]->nb_element - nb_equal / 2. - nb_sup) /
          (histo[0]->nb_element * (double)histo[1]->nb_element) << "   "
       << STAT_label[STATL_MANN_WHITNEY_EQUAL_PROBABILITY] << " = "
       << nb_equal / (histo[0]->nb_element * (double)histo[1]->nb_element) << "   "
       << STAT_label[STATL_MANN_WHITNEY_SUPERIOR_PROBABILITY] << " = "
       << (nb_sup - nb_equal / 2.) / (histo[0]->nb_element * (double)histo[1]->nb_element) << endl;
#   endif

    delete [] histo;
    delete merged_histo;
    delete [] rank;

    delete test;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'une famille de lois empiriques :
 *  - histogrammes seuls et histogrammes etages,
 *  - lois et fonctions de repartition,
 *  - mise en correspondance des fonctions de repartition,
 *  - courbes de concentration.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              nombre de lois empiriques, pointeurs sur les lois empiriques, titre des figures.
 *
 *--------------------------------------------------------------*/

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
    register int i , j , k;
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

    // ecriture des fichiers de donnees

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

        // calcul des fonctions de repartition et de concentration

        cumul[i] = histo[i]->cumul_computation();
        concentration[i] = histo[i]->concentration_function_computation();

        // calcul du nombre de valeurs maximum, de l'etendue maximum,
        // de la frequence maximum et de la probabilite maximum

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

      // ecriture du fichier de commandes et du fichier d'impression

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

        // histogrammes decales

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

          // histogrammes etages

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

          // lois

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

          // fonctions de repartition

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

          // mise en correspondance des fonctions de repartition en prenant
          // comme reference la loi empirique dont l'etendue est la plus grande

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

          // courbes de concentration

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


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'une famille de lois empiriques :
 *  - histogrammes seuls et histogrammes etages,
 *  - lois et fonctions de repartition,
 *  - mise en correspondance des fonctions de repartition,
 *  - courbes de concentration.
 *
 *  arguments : reference sur un objet StatError, nombre de lois empiriques,
 *              pointeurs sur les lois empiriques, titre des figures.
 *
 *--------------------------------------------------------------*/

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
    register int i , j , k;
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

      // calcul des fonctions de repartition

      cumul[i] = histo[i]->cumul_computation();

      // calcul du nombre de valeurs maximum, de l'etendue maximum,
      // de la frequence maximum et de la probabilite maximum

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

    // nombre de fenetres

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

    // 1ere vue : histogrammes decales

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

      // 2eme vue : histogrammes etages

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

      // 3eme vue : lois

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

      // 4eme vue : fonctions de repartition

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

      // 5eme vue : mise en correspondance des fonctions de repartition en prenant
      // comme reference la loi empirique dont l'etendue est la plus grande

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

      // 6eme vue : courbes de concentration

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


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'une loi empirique.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* FrequencyDistribution::get_plotable() const

{
  StatError error;

  return get_plotable_frequency_distributions(error , 0 , NULL);
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet DiscreteDistributionData a partir
 *  d'un objet FrequencyDistribution et d'un objet Distribution.
 *
 *  arguments : reference sur un objet FrequencyDistribution,
 *              pointeur sur un objet Distribution.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet DiscreteDistributionData a partir
 *  d'un objet FrequencyDistribution et d'un objet DiscreteParametric.
 *
 *  arguments : reference sur un objet FrequencyDistribution,
 *              pointeur sur un objet DiscreteParametric.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe DiscreteDistributionData.
 *
 *  arguments : reference sur un objet DiscreteDistributionData,
 *              flag copie de l'objet DiscreteParametricModel.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe DiscreteDistributionData.
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData::~DiscreteDistributionData()

{
  delete distribution;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe DiscreteDistributionData.
 *
 *  argument : reference sur un objet DiscreteDistributionData.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Extraction de la partie "modele" d'un objet DiscreteDistributionData.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet FrequencyDistribution a partir d'un fichier.
 *  Format : n lignes de la forme (valeur) (frequence).
 *  Les valeurs sont ordonnees.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData* frequency_distribution_ascii_read(StatError &error , const char *path)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status , lstatus;
  register int i;
  int line , nb_element;
  long value , index , max_index;
  DiscreteDistributionData *histo;
  ifstream in_file(path);


  histo = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {

    // 1ere passe : analyse des lignes et recherche du nombre de valeurs

    status = true;
    line = 0;
    max_index = -1;
    nb_element = 0;

    while (buffer.readLine(in_file , false)) {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.first('#');
      if (position != RW_NPOS) {
        buffer.remove(position);
      }
      i = 0;

      RWCTokenizer next(buffer);

      while (!((token = next()).isNull())) {
          if (i <= 1) {
          lstatus = locale.stringToNum(token , &value);
          if ((lstatus) && (value < 0)) {
            lstatus = false;
          }

          if (!lstatus) {
            status = false;
            error.update(STAT_parsing[STATP_DATA_TYPE] , line , i + 1);
          }

          else {
            switch (i) {

            // test valeurs ordonnees

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

      // test deux valeurs par ligne

      if ((i > 0) && (i != 2)) {
        status = false;
        error.correction_update(STAT_parsing[STATP_NB_TOKEN] , 2 , line);
      }
    }

    if (nb_element == 0) {
      status = false;
      error.update(STAT_parsing[STATP_EMPTY_SAMPLE]);
    }

    // 2eme passe : lecture du fichier

    if (status) {
//      in_file.close();
//      in_file.open(path , ios::in);

      in_file.clear();
      in_file.seekg(0 , ios::beg);

      histo = new DiscreteDistributionData(max_index + 1);

      while (buffer.readLine(in_file , false)) {
        position = buffer.first('#');
        if (position != RW_NPOS) {
          buffer.remove(position);
        }
        i = 0;

        RWCTokenizer next(buffer);

        while (!((token = next()).isNull())) {
          switch (i) {
          case 0 :
            locale.stringToNum(token , &index);
            break;
          case 1 :
            locale.stringToNum(token , &value);
            histo->frequency[index] = value;
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


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet DiscreteDistributionData.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet DiscreteDistributionData.
 *
 *  arguments : stream, flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet DiscreteDistributionData.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& DiscreteDistributionData::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet DiscreteDistributionData dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool DiscreteDistributionData::ascii_write(StatError &error , const char *path ,
                                           bool exhaustive) const

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
    ascii_write(out_file , exhaustive , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet DiscreteDistributionData dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool DiscreteDistributionData::spreadsheet_write(StatError &error , const char *path) const

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

    if (distribution) {
      distribution->spreadsheet_write(out_file , this);
    }

    else {
      double information = information_computation();


      out_file << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      spreadsheet_characteristic_print(out_file , true);

      out_file << STAT_label[STATL_MEAN_ABSOLUTE_DEVIATION] << "\t" << mean_absolute_deviation_computation();
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


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet DiscreteDistributionData.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet DiscreteDistributionData.
 *
 *--------------------------------------------------------------*/

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
