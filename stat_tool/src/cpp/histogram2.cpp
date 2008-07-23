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



#include <limits.h>
#include <math.h>
#include <sstream>
#include <iomanip>
#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"
#include "tool/config.h"

// #include <rw/vstream.h>
// #include <rw/rwfile.h>
#include "stat_tools.h"
#include "distribution.h"
#include "stat_label.h"

using namespace std;


extern int column_width(int value);
extern int column_width(int nb_value , const double *value , double scale = 1.);
extern bool cumul_matching_plot_print(const char *path , int nb_cumul , int *offset ,
                                      int *nb_value , double **cumul);
extern char* label(const char *file_name);



/*--------------------------------------------------------------*
 *
 *  Ecriture des resultats d'une comparaison d'histogrammes.
 *
 *  arguments : stream, nombre d'histogrammes, pointeurs sur les histogrammes,
 *              type de variable (SYMBOLIC/ORDINAL/NUMERIC), dissimilarites.
 *
 *--------------------------------------------------------------*/

ostream& dissimilarity_ascii_write(ostream &os , int nb_histo , const Histogram **histo ,
                                   int type , double **dissimilarity)

{
  register int i , j;
  int nb_value , buff , width[3];
  long old_adjust;
  double information , **cumul;
  Test *test;


  old_adjust = os.setf(ios::right , ios::adjustfield);

  // ecriture des caracteristiques des histogrammes

  for (i = 0;i < nb_histo;i++) {
    os << STAT_label[STATL_HISTOGRAM] << " " << i + 1  << " - ";

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

  nb_value = histo[0]->nb_value;
  for (i = 1;i < nb_histo;i++) {
    if (histo[i]->nb_value > nb_value) {
      nb_value = histo[i]->nb_value;
    }
  }

  width[0] = column_width(nb_value - 1);

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

  // ecriture des histogrammes et des fonctions de repartition

  os << "  ";
  for (i = 0;i < nb_histo;i++) {
    os << " | " << STAT_label[STATL_HISTOGRAM] << " " << i + 1;
  }
  for (i = 0;i < nb_histo;i++) {
    os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
       << i + 1 << " " << STAT_label[STATL_FUNCTION];
  }
  os << endl;

  for (i = 0;i < nb_value;i++) {
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

  // ecriture de la matrice des dissimilarites entre lois deduites des histogrammes

  os << "\n" << STAT_label[STATL_DISSIMILARITIES] << " between " << STAT_label[STATL_HISTOGRAMS] << endl;

  os << "\n           ";
  for (i = 0;i < nb_histo;i++) {
    os << " | " << STAT_label[STATL_HISTOGRAM] << " " << i + 1;
  }
  os << endl;

  for (i = 0;i < nb_histo;i++) {
    os << STAT_label[STATL_HISTOGRAM] << " ";
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
    Histogram *merged_histo;


    merged_histo = new Histogram(nb_histo , histo);

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

    test = new Test(FISHER , true , df[0] , df[1] , mean_square[0] / mean_square[1]);
    test->F_critical_probability_computation();

    test->ascii_print(os , false , (df[0] == 1 ? false : true));

    delete test;
    break;
  }
  }

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des resultats d'une comparaison d'histogrammes dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path, nombre d'histogrammes,
 *              pointeurs sur les histogrammes, type de variable (SYMBOLIC/ORDINAL/NUMERIC),
 *              dissimilarites.
 *
 *--------------------------------------------------------------*/

bool dissimilarity_ascii_write(Format_error &error , const char *path , int nb_histo ,
                               const Histogram **histo , int type , double **dissimilarity)

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
    dissimilarity_ascii_write(out_file , nb_histo , histo , type , dissimilarity);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des resultats d'une comparaison d'histogrammes dans un fichier
 *  au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path, nombre d'histogrammes,
 *              pointeurs sur les histogrammes, type de variable (SYMBOLIC/ORDINAL/NUMERIC),
 *              dissimilarites.
 *
 *--------------------------------------------------------------*/

bool dissimilarity_spreadsheet_write(Format_error &error , const char *path , int nb_histo ,
                                     const Histogram **histo , int type , double **dissimilarity)

{
  bool status;
  register int i , j;
  int nb_value;
  double information , **cumul , **concentration;
  Test *test;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    // ecriture des caracteristiques des histogrammes

    for (i = 0;i < nb_histo;i++) {
      out_file << STAT_label[STATL_HISTOGRAM] << " " << i + 1  << "\t";

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

    // ecriture des histogrammes, des fonctions de repartition et des courbes de concentration

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

    nb_value = histo[0]->nb_value;
    for (i = 1;i < nb_histo;i++) {
      if (histo[i]->nb_value > nb_value) {
        nb_value = histo[i]->nb_value;
      }
    }

    for (i = 0;i < nb_histo;i++) {
      out_file << "\t" << STAT_label[STATL_HISTOGRAM] << " " << i + 1;
    }
    for (i = 0;i < nb_histo;i++) {
      out_file << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
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

    for (i = 0;i < nb_value;i++) {
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

    // ecriture de la matrice des dissimilarites entre lois deduites des histogrammes

    out_file << "\n" << STAT_label[STATL_DISSIMILARITIES] << " between " << STAT_label[STATL_HISTOGRAMS] << "\n\n";

    for (i = 0;i < nb_histo;i++) {
      out_file << "\t" << STAT_label[STATL_HISTOGRAM] << " " << i + 1;
    }
    out_file << endl;

    for (i = 0;i < nb_histo;i++) {
      out_file << STAT_label[STATL_HISTOGRAM] << " " << i + 1;
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
      Histogram *merged_histo;


      merged_histo = new Histogram(nb_histo , histo);

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

      test = new Test(FISHER , true , df[0] , df[1] , mean_square[0] / mean_square[1]);
      test->F_critical_probability_computation();

      test->spreadsheet_print(out_file , (df[0] == 1 ? false : true));

      delete test;
      break;
    }
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Test de Kruskal-Wallis (analyse de variance dans le cas ordinal).
 *
 *  arguments : nombre d'histogrammes, pointeurs sur les histogrammes.
 *
 *--------------------------------------------------------------*/

Test* Histogram::kruskal_wallis_test(int nb_histo , const Histogram **ihisto) const

{
  register int i , j;
  int *pfrequency;
  double correction , value , sum , *rank;
  const Histogram **histo;
  Histogram *merged_histo;
  Test *test;


  nb_histo++;
  histo = new const Histogram*[nb_histo];

  histo[0] = this;
  for (i = 1;i < nb_histo;i++) {
    histo[i] = ihisto[i - 1];
  }

  merged_histo = new Histogram(nb_histo , histo);

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
    pfrequency = histo[i]->frequency + histo[i]->offset;
    sum = 0.;
    for (j = histo[i]->offset;j < histo[i]->nb_value;j++) {
      sum += *pfrequency++ * rank[j];
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
 *  Comparaison d'histogrammes.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre d'histogrammes,
 *              pointeurs sur les histogrammes, type de variable (SYMBOLIC/ORDINAL/NUMERIC),
 *              path, format de fichier ('a' : ASCII, 's' : Spreadsheet).
 *
 *--------------------------------------------------------------*/

bool Histogram::comparison(Format_error &error , ostream &os , int nb_histo ,
                           const Histogram **ihisto , int type , const char *path ,
                           char format) const

{
  bool status = true;
  register int i , j , k , m;
  int max_nb_value , distance = 1;
  double **dissimilarity;
  Distribution **dist;
  const Histogram **histo;


  nb_histo++;
  histo = new const Histogram*[nb_histo];

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
    int buff , *pfrequency , width[2];
    double *rank , *rank_mean;
    Histogram *merged_histo;


    merged_histo = new Histogram(nb_histo , histo);

    // calcul des rangs

    rank = merged_histo->rank_computation();

    // calcul des rangs moyens pour chaque histogramme

    rank_mean = new double[nb_histo];

    for (i = 0;i < nb_histo;i++) {
      pfrequency = histo[i]->frequency + histo[i]->offset;
      rank_mean[i] = 0.;
      for (j = histo[i]->offset;j < histo[i]->nb_value;j++) {
        rank_mean[i] += *pfrequency++ * rank[j];
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

    // ecriture de la matrice des dissimilarites entre lois deduites des histogrammes

    os << "\n" << STAT_label[STATL_DISSIMILARITIES] << " between " << STAT_label[STATL_HISTOGRAMS] << endl;

    os << "\n           ";
    for (i = 0;i < nb_histo;i++) {
      os << " | " << STAT_label[STATL_HISTOGRAM] << " " << i + 1;
    }
    os << endl;

    for (i = 0;i < nb_histo;i++) {
      os << STAT_label[STATL_HISTOGRAM] << " ";
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
  dissimilarity_ascii_write(os , nb_histo , histo , type , dissimilarity);
# endif

  if (path) {
    switch (format) {
    case 'a' :
      status = dissimilarity_ascii_write(error , path , nb_histo , histo , type , dissimilarity);
      break;
    case 's' :
      status = dissimilarity_spreadsheet_write(error , path , nb_histo , histo , type , dissimilarity);
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
 *  Test F de comparaison des variances de 2 histogrammes.
 *
 *  arguments : stream, reference sur un histogramme.
 *
 *--------------------------------------------------------------*/

void Histogram::F_comparison(ostream &os , const Histogram &histo) const

{
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

# ifdef MESSAGE
  os << *test;
# endif

  delete test;
}


/*--------------------------------------------------------------*
 *
 *  Test t de comparaison des moyennes de 2 histogrammes.
 *
 *  arguments : stream, reference sur un histogramme.
 *
 *--------------------------------------------------------------*/

void Histogram::t_comparison(ostream &os , const Histogram &histo) const

{
  int df;
  double value , buff1 , buff2;
  Test *test;


  buff1 = variance / nb_element;
  buff2 = histo.variance / histo.nb_element;

  value = fabs(mean - histo.mean) / sqrt(buff1 + buff2);
  df = (int)round((buff1 + buff2) * (buff1 + buff2) /
                  (buff1 * buff1 / (nb_element - 1) + buff2 * buff2 / (histo.nb_element - 1)));

  test = new Test(STUDENT , false , df , I_DEFAULT , value);

  test->t_critical_probability_computation();

# ifdef MESSAGE
  os << *test;
# endif

  delete test;

# ifdef DEBUG
  value = fabs(mean - histo.mean) /
          sqrt(((nb_element - 1) * variance + (histo.nb_element - 1) * histo.variance) /
               (nb_element + histo.nb_element - 2) *
               (1. / (double)nb_element + 1. / (double)histo.nb_element));

  test = new Test(STUDENT , false , nb_element + histo.nb_element - 2 , I_DEFAULT , value);

  test->t_critical_probability_computation();

  os << "\n" << *test;

  delete test;
# endif

}


/*--------------------------------------------------------------*
 *
 *  Test de comparaison de Wilcoxon-Mann-Whitney.
 *
 *  arguments : reference sur un objet Format_error, stream,
 *              reference sur un histogramme.
 *
 *--------------------------------------------------------------*/

bool Histogram::wilcoxon_mann_whitney_comparison(Format_error &error , ostream &os ,
                                                 const Histogram &ihisto) const

{
  bool status;
  register int i;
  int nb_equal , min , max , *pfrequency0 , *pfrequency1;
  double correction , value , nb_sup , mean , variance , *rank;
  const Histogram **histo;
  Histogram *merged_histo;
  Test *test;


  if (nb_element * (double)ihisto.nb_element > INT_MAX) {
    status = false;
    error.update(STAT_error[STATR_SAMPLE_SIZES]);
  }

  else {
    status = true;

    histo = new const Histogram*[2];

    histo[0] = this;
    histo[1] = &ihisto;

    merged_histo = new Histogram(2 , histo);

    // calcul du terme de correction pour les ex-aequo

    pfrequency0 = merged_histo->frequency + merged_histo->offset;
    correction = 0.;
    for (i = merged_histo->offset;i < merged_histo->nb_value;i++) {
      if (*pfrequency0 > 1) {
        correction += *pfrequency0 * ((double)*pfrequency0 * (double)*pfrequency0 - 1);
      }
      pfrequency0++;
    }

    // calcul des rangs

    rank = merged_histo->rank_computation();

    // calcul de la statistique de Wilcoxon

    pfrequency0 = histo[0]->frequency + histo[0]->offset;
    value = 0.;
    for (i = histo[0]->offset;i < histo[0]->nb_value;i++) {
      value += *pfrequency0++ * rank[i];
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

    variance = histo[0]->nb_element * (double)histo[1]->nb_element * ((merged_histo->nb_element *
                 ((double)merged_histo->nb_element * (double)merged_histo->nb_element - 1)) - correction) /
               (merged_histo->nb_element * ((double)merged_histo->nb_element - 1) * 12.);

    value = fabs(value - mean) / sqrt(variance);

    test = new Test(STANDARD_NORMAL , false , I_DEFAULT , I_DEFAULT , value);

    test->standard_normal_critical_probability_computation();

    min = MAX(histo[0]->offset , histo[1]->offset);
    max = MIN(histo[0]->nb_value , histo[1]->nb_value) - 1;
    nb_equal = 0;

    if (max >= min) {
      pfrequency0 = histo[0]->frequency + min;
      pfrequency1 = histo[1]->frequency + min;

      for (i = min;i <= max;i++) {
        nb_equal += *pfrequency0++ * *pfrequency1++;
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
 *  Sortie Gnuplot d'une famille d'histogrammes :
 *  - histogrammes seuls et histogrammes etages,
 *  - lois et fonctions de repartition deduites des histogrammes,
 *  - mise en correspondance des fonctions de repartition,
 *  - courbes de concentration.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              nombre d'histogrammes, pointeurs sur les histogrammes, titre des figures.
 *
 *--------------------------------------------------------------*/

bool plot_write(Format_error &error , const char *prefix , int nb_histo ,
                const Histogram **histo , const char *title)

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
        reference_concentration , *offset , *nb_value;
    double max_mass , shift , **cumul , **concentration;
    const Histogram *phisto[2];
    const Histogram **merged_histo;
    ostringstream *data_file_name;


    // ecriture des fichiers de donnees

    data_file_name = new ostringstream[nb_histo + 2];

    if (nb_histo > 1) {
      data_file_name[0] << prefix << 0 << ".dat";

      merged_histo = new const Histogram*[nb_histo];
      merged_histo[nb_histo - 1] = new Histogram(*histo[nb_histo - 1]);

      for (i = nb_histo - 2;i >= 0;i--) {
        phisto[0] = merged_histo[i + 1];
        phisto[1] = histo[i];
        merged_histo[i] = new Histogram(2 , phisto);
      }

      status = merged_histo[0]->plot_print((data_file_name[0].str()).c_str() ,
                                            nb_histo - 1 , merged_histo + 1);
    }

    if (status) {
      offset = new int[nb_histo];
      nb_value = new int[nb_histo];

      cumul = new double*[nb_histo];
      concentration = new double*[nb_histo];
      for (i = 0;i < nb_histo;i++) {
        cumul[i] = 0;
        concentration[i] = 0;
      }

      max_nb_value = 0;
      max_frequency = 0;
      max_mass = 0.;
      cumul_concentration_flag = false;
      max_range = 0;
      shift = 0.;

      for (i = 0;i < nb_histo;i++) {
        data_file_name[i + 1] << prefix << i + 1 << ".dat";

        offset[i] = histo[i]->offset;
        nb_value[i] = histo[i]->nb_value;

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

        if (histo[i]->variance > 0.) {
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
                                  offset , nb_value , cumul);
      }

      delete [] offset;
      delete [] nb_value;

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
                   << STAT_label[STATL_HISTOGRAM];
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
                     << j + 1 << " title \"" << STAT_label[STATL_HISTOGRAM] << " " << j + 1
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

          // lois deduite des histogrammes

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

          // fonctions de repartition deduites des histogrammes

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;

          if (max_nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }

          out_file << "plot [0:" << max_nb_value - 1 << "] [0:1] ";
          k = 0;
          for (j = 0;j < nb_histo;j++) {
            if (histo[j]->variance > 0.) {
              if (k > 0) {
                out_file << ",\\\n";
              }
              k++;
              out_file << "\"" << label((data_file_name[j + 1].str()).c_str()) << "\" using 1:5 title \""
                       << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM];
              if (nb_histo > 1) {
                out_file << " " << j + 1;
              }
              out_file << " " << STAT_label[STATL_FUNCTION] << "\" with linespoints";
            }
          }
          out_file << endl;

          if (max_nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }

          // mise en correspondance des fonctions de repartition en prenant
          // comme reference l'histogramme dont l'etendue est la plus grande

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
            out_file << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM]
                     << " " << STAT_label[STATL_FUNCTION] << " " << STAT_label[STATL_MATCHING] << "\"\n\n";

            out_file << "set grid\n" << "set xtics 0,0.1\n" << "set ytics 0,0.1" << endl;

            out_file << "plot [0:1] [0:1] ";
            k = 0;
            for (j = 0;j < nb_histo;j++) {
              if (histo[j]->variance > 0.) {
                if (k > 0) {
                  out_file << ",\\\n";
                }
                k++;
                out_file << "\"" << label((data_file_name[nb_histo + 1].str()).c_str()) << "\" using "
                         << reference_matching + 1 << ":" << j + 1 << " title \""
                         << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
                         << j + 1 << " " << STAT_label[STATL_FUNCTION] << "\" with linespoints";
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
            if (histo[j]->variance > 0.) {
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

    delete [] data_file_name;

    if (!status) {
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Distribution_data a partir
 *  d'un objet Histogram et d'un objet Distribution.
 *
 *  arguments : reference sur un objet Histogram,
 *              pointeur sur un objet Distribution.
 *
 *--------------------------------------------------------------*/

Distribution_data::Distribution_data(const Histogram &histo , const Distribution *dist)
:Histogram(histo)

{
  if (dist) {
    distribution = new Parametric_model(*dist);
  }
  else {
    distribution = 0;
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Distribution_data a partir
 *  d'un objet Histogram et d'un objet Parametric.
 *
 *  arguments : reference sur un objet Histogram,
 *              pointeur sur un objet Parametric.
 *
 *--------------------------------------------------------------*/

Distribution_data::Distribution_data(const Histogram &histo , const Parametric *dist)
:Histogram(histo)

{
  if (dist) {
    distribution = new Parametric_model(*dist);
  }
  else {
    distribution = 0;
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe Distribution_data.
 *
 *  arguments : reference sur un objet Distribution_data,
 *              flag copie de l'objet Parametric_model.
 *
 *--------------------------------------------------------------*/

Distribution_data::Distribution_data(const Distribution_data &histo , bool model_flag)
:Histogram(histo)

{
  if ((model_flag) && (histo.distribution)) {
    distribution = new Parametric_model(*(histo.distribution) , false);
  }
  else {
    distribution = 0;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Distribution_data.
 *
 *--------------------------------------------------------------*/

Distribution_data::~Distribution_data()

{
  delete distribution;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Distribution_data.
 *
 *  argument : reference sur un objet Distribution_data.
 *
 *--------------------------------------------------------------*/

Distribution_data& Distribution_data::operator=(const Distribution_data &histo)

{
  if (&histo != this) {
    delete distribution;
    delete [] frequency;

    Histogram::copy(histo);
    if (histo.distribution) {
      distribution = new Parametric_model(*(histo.distribution) , false);
    }
    else {
      distribution = 0;
    }
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de la partie "modele" d'un objet Distribution_data.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Parametric_model* Distribution_data::extract_model(Format_error &error) const

{
  Parametric_model *dist;


  error.init();

  if (!distribution) {
    dist = 0;
    error.update(STAT_error[STATR_NON_EXISTING_DISTRIBUTION]);
  }

  else {
    dist = new Parametric_model(*distribution);
    dist->histogram = new Distribution_data(*this , false);
  }

  return dist;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Histogram a partir d'un fichier.
 *  Format : n lignes de la forme (valeur) (frequence).
 *  Les valeurs sont ordonnees.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

Distribution_data* histogram_ascii_read(Format_error &error , const char *path)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status , lstatus;
  register int i;
  int line , nb_element;
  long value , index , max_index;
  Distribution_data *histo;
  ifstream in_file(path);


  histo = 0;
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

      histo = new Distribution_data(max_index + 1);

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
 *  Ecriture sur une ligne d'un objet Distribution_data.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Distribution_data::line_write(ostream &os) const

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
 *  Ecriture d'un objet Distribution_data.
 *
 *  arguments : stream, flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& Distribution_data::ascii_write(ostream &os , bool exhaustive ,
                                        bool file_flag) const

{
  if (distribution) {
    distribution->ascii_write(os , this , exhaustive , file_flag);
  }
  else {
    Histogram::ascii_write(os , exhaustive , file_flag);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Distribution_data.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Distribution_data::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Distribution_data dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Distribution_data::ascii_write(Format_error &error , const char *path ,
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
 *  Ecriture d'un objet Distribution_data dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Distribution_data::spreadsheet_write(Format_error &error , const char *path) const

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


      out_file << STAT_label[STATL_HISTOGRAM] << "\t";
      spreadsheet_characteristic_print(out_file , true);

      out_file << STAT_label[STATL_MEAN_ABSOLUTE_DEVIATION] << "\t" << mean_absolute_deviation_computation();
      if (mean > 0.) {
        out_file << "\t\t" << STAT_label[STATL_CONCENTRATION_COEFF] << "\t" << concentration_computation();
      }
      out_file << endl;

      out_file << STAT_label[STATL_INFORMATION] << "\t" << information << "\t"
               << information / nb_element << endl;

      out_file << "\n\t" << STAT_label[STATL_HISTOGRAM] << "\t" << STAT_label[STATL_CUMULATIVE]
               << " " << STAT_label[STATL_HISTOGRAM] << " " << STAT_label[STATL_FUNCTION];
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
 *  Sortie Gnuplot d'un objet Distribution_data.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Distribution_data::plot_write(Format_error &error , const char *prefix ,
                                   const char *title) const

{
  bool status;


  if (distribution) {
    status = distribution->plot_write(error , prefix , title , this);
  }

  else {
    const Histogram *histo[1];

    histo[0] = this;
    status = ::plot_write(error , prefix , 1 , histo , title);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWDEFINE_COLLECTABLE(Distribution_data , STATI_DISTRIBUTION_DATA);


RWspace Distribution_data::binaryStoreSize() const

{
  RWspace size = Histogram::binaryStoreSize();
  if (distribution) {
    size += distribution->recursiveStoreSize();
  }

  return size;
}


void Distribution_data::restoreGuts(RWvistream &is)

{
  delete distribution;

  Histogram::restoreGuts(is);

  is >> distribution;
  if (distribution == RWnilCollectable) {
    distribution = 0;
  }
}


void Distribution_data::restoreGuts(RWFile &file)

{
  delete distribution;

  Histogram::restoreGuts(file);

  file >> distribution;
  if (distribution == RWnilCollectable) {
    distribution = 0;
  }
}


void Distribution_data::saveGuts(RWvostream &os) const

{
  Histogram::saveGuts(os);

  if (distribution) {
    os << distribution;
  }
  else {
    os << RWnilCollectable;
  }
}


void Distribution_data::saveGuts(RWFile &file) const

{
  Histogram::saveGuts(file);

  if (distribution) {
    file << distribution;
  }
  else {
    file << RWnilCollectable;
  }
} */
