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
 *       $Id: correlation.cpp 18049 2015-04-23 09:42:18Z guedon $
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



#include <math.h>
#include <sstream>
#include <iomanip>

#include <boost/math/distributions/normal.hpp>

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

using namespace std;
using namespace boost::math;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Correlation.
 *
 *--------------------------------------------------------------*/

Correlation::Correlation()

{
  type = PEARSON;

  variable_type = NULL;

  variable1 = NULL;
  variable2 = NULL;

  white_noise = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Correlation.
 *
 *  arguments : type, decalage maximum, variables.
 *
 *--------------------------------------------------------------*/

Correlation::Correlation(int itype , int max_lag , int ivariable1 , int ivariable2)
:Curves(1 , max_lag + 1 , true , false)

{
  type = itype;

  variable_type = new int[1];
  variable_type[0] = OBSERVED_VALUE;

  variable1 = new int[1];
  variable2 = new int[1];
  variable1[0] = ivariable1;
  variable2[0] = ivariable2;

  white_noise = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Correlation.
 *
 *  arguments : nombre de fonctions de correlation, nombre de points,
 *              flag sur les frequences, type.
 *
 *--------------------------------------------------------------*/

Correlation::Correlation(int inb_curve , int ilength , bool frequency_flag , int itype)
:Curves(inb_curve , ilength , frequency_flag , false)

{
  type = itype;

  variable_type = new int[nb_curve];

  variable1 = new int[nb_curve];
  variable2 = new int[nb_curve];

  white_noise = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Correlation.
 *
 *  argument : reference sur un objet Correlation.
 *
 *--------------------------------------------------------------*/

void Correlation::copy(const Correlation &correl)

{
  register int i;


  type = correl.type;

  variable_type = new int[nb_curve];
  for (i = 0;i < nb_curve;i++) {
    variable_type[i] = correl.variable_type[i];
  }

  variable1 = new int[nb_curve];
  variable2 = new int[nb_curve];
  for (i = 0;i < nb_curve;i++) {
    variable1[i] = correl.variable1[i];
    variable2[i] = correl.variable2[i];
  }

  if (correl.white_noise) {
    white_noise = new double[length];
    for (i = 0;i < length;i++) {
      white_noise[i] = correl.white_noise[i];
    }
  }
  else {
    white_noise = NULL;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Correlation.
 *
 *--------------------------------------------------------------*/

void Correlation::remove()

{
  delete [] variable_type;

  delete [] variable1;
  delete [] variable2;

  delete [] white_noise;
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Correlation.
 *
 *--------------------------------------------------------------*/

Correlation::~Correlation()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Correlation.
 *
 *  argument : reference sur un objet Correlation.
 *
 *--------------------------------------------------------------*/

Correlation& Correlation::operator=(const Correlation &correl)

{
  if (&correl != this) {
    remove();
    Curves::remove();

    Curves::copy(correl);
    copy(correl);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Fusion d'objets Correlation.
 *
 *  arguments : reference sur un objet StatError, nombre d'objets Correlation,
 *              pointeurs sur les objets Correlation.
 *
 *--------------------------------------------------------------*/

Correlation* Correlation::merge(StatError &error , int nb_correl ,
                                const Correlation **icorrel) const

{
  bool status = true;
  register int i , j , k , m;
  int inb_curve , *pfrequency;
  Correlation *correl;
  const Correlation **pcorrel;


  correl = NULL;
  error.init();

  pfrequency = frequency;

  for (i = 0;i < nb_correl;i++) {
    if ((icorrel[i]->type != type) || (icorrel[i]->offset != offset)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << SEQ_label[SEQL_CORRELATION_FUNCTION] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_CORRELATION_COEFF_TYPE];

      switch (type) {
      case SPEARMAN :
        correction_message << SEQ_label[SEQL_SPEARMAN] << " ";
        break;
      case KENDALL :
        correction_message << SEQ_label[SEQL_KENDALL] << " ";
        break;
      }

      if (offset == 1) {
        correction_message << SEQ_label[SEQL_PARTIAL] << " ";
      }

      if ((type == SPEARMAN) || (type == KENDALL)) {
        correction_message << SEQ_label[SEQL_RANK] << " ";
      }

      correction_message << SEQ_label[SEQL_CORRELATION_FUNCTION];

      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }

    if (icorrel[i]->length != length) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_CORRELATION_FUNCTION] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_MAX_LAG];
      error.correction_update((error_message.str()).c_str() , length - 1);
    }

    else if (icorrel[i]->frequency) {
      if (!pfrequency) {
        pfrequency = icorrel[i]->frequency;
      }

      else {
        for (j = offset;j < length;j++) {
          if (icorrel[i]->frequency[j] != pfrequency[j]) {
            status = false;
            ostringstream error_message;
            error_message << SEQ_label[SEQL_CORRELATION_FUNCTION] << " " << i + 2 << ": "
                          << SEQ_label[SEQL_LAG] << " " << j << ": " << SEQ_error[SEQR_FREQUENCY];
            error.update((error_message.str()).c_str());
          }
        }
      }
    }
  }

  if (status) {
    nb_correl++;
    pcorrel = new const Correlation*[nb_correl];

    pcorrel[0] = this;
    for (i = 1;i < nb_correl;i++) {
      pcorrel[i] = icorrel[i - 1];
    }

    inb_curve = 0;
    for (i = 0;i < nb_correl;i++) {
      inb_curve += pcorrel[i]->nb_curve;
    }

    correl = new Correlation(inb_curve , length , (pfrequency ? true : false) , type);

    correl->offset = offset;

    i = 0;
    for (j = 0;j < nb_correl;j++) {
      for (k = 0;k < pcorrel[j]->nb_curve;k++) {
        correl->variable_type[i] = pcorrel[j]->variable_type[k];
        correl->variable1[i] = pcorrel[j]->variable1[k];
        correl->variable2[i++] = pcorrel[j]->variable2[k];
      }
    }

    if (pfrequency) {
      for (i = 0;i < length;i++) {
        correl->frequency[i] = pfrequency[i];
      }
    }

    i = 0;
    for (j = 0;j < nb_correl;j++) {
      for (k = 0;k < pcorrel[j]->nb_curve;k++) {
        for (m = 0;m < length;m++) {
          correl->point[i][m] = pcorrel[j]->point[k][m];
        }
        i++;
      }
    }

    delete [] pcorrel;
  }

  return correl;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet Correlation.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Correlation::line_write(ostream &os) const

{
  register int i;
  int autocorrelation , cross_correlation;


  switch (type) {
  case SPEARMAN :
    os << SEQ_label[SEQL_SPEARMAN] << " ";
    break;
  case KENDALL :
    os << SEQ_label[SEQL_KENDALL] << " ";
    break;
  }

  if (offset == 1) {
    os << SEQ_label[SEQL_PARTIAL] << " ";
  }

  if ((type == SPEARMAN) || (type == KENDALL)) {
    os << SEQ_label[SEQL_RANK] << " ";
  }

  autocorrelation = true;
  cross_correlation = true;
  for (i = 0;i < nb_curve;i++) {
    if (variable_type[i] == OBSERVED_VALUE) {
      if (variable1[i] != variable2[i]) {
        autocorrelation = false;
      }
      else if (variable1[i] == variable2[i]) {
        cross_correlation = false;
      }
    }

    else {
      cross_correlation = false;
    }
  }

  if (autocorrelation) {
    os << SEQ_label[SEQL_AUTO];
  }
  else if (cross_correlation) {
    os << SEQ_label[SEQL_CROSS];
  }

  os << SEQ_label[SEQL_CORRELATION_FUNCTION] << "   " << SEQ_label[SEQL_MAX_LAG] << ": " << length - 1;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Correlation.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Correlation::ascii_write(ostream &os , bool exhaustive) const

{
  bool autocorrelation , cross_correlation;
  register int i , j;
  int *width;
  double standard_normal_value , *confidence_limit;
  long old_adjust;


  old_adjust = os.setf(ios::right , ios::adjustfield);

  // calcul des limites de confiance

  if (frequency) {
    normal dist;
    standard_normal_value = quantile(complement(dist , 0.025));

    confidence_limit = new double[length];

    for (i = 0;i < length;i++) {
      switch (type) {
      case PEARSON :
        confidence_limit[i] = standard_normal_value / sqrt((double)frequency[i]);
        break;
      case SPEARMAN :
        confidence_limit[i] = standard_normal_value / sqrt((double)frequency[i]);
        break;
      case KENDALL :
        confidence_limit[i] = standard_normal_value * sqrt((2 * (2 * (double)frequency[i] + 5)) /
                               (9 * (double)frequency[i] * (double)(frequency[i] - 1)));
        break;
      }
    }
  }

  os << "\n";
  switch (type) {
  case SPEARMAN :
    os << SEQ_label[SEQL_SPEARMAN] << " ";
    break;
  case KENDALL :
    os << SEQ_label[SEQL_KENDALL] << " ";
    break;
  }

  if (offset == 1) {
    os << SEQ_label[SEQL_PARTIAL] << " ";
  }

  if ((type == SPEARMAN) || (type == KENDALL)) {
    os << SEQ_label[SEQL_RANK] << " ";
  }

  autocorrelation = true;
  cross_correlation = true;
  for (i = 0;i < nb_curve;i++) {
    if (variable_type[i] == OBSERVED_VALUE) {
      if (variable1[i] != variable2[i]) {
        autocorrelation = false;
      }
      else if (variable1[i] == variable2[i]) {
        cross_correlation = false;
      }
    }

    else {
      cross_correlation = false;
    }
  }

  if (autocorrelation) {
    os << SEQ_label[SEQL_AUTO];
  }
  else if (cross_correlation) {
    os << SEQ_label[SEQL_CROSS];
  }
  os << SEQ_label[SEQL_CORRELATION_FUNCTION];

  // calcul des largeurs des colonnes

  width = new int[nb_curve + 4];

  width[0] = column_width(length - 1);
  for (i = 0;i < nb_curve;i++) {
    width[i + 1] = column_width(length - offset , point[i] + offset) + ASCII_SPACE;
  }

  i = nb_curve + 1;
  if (white_noise) {
    width[i++] = column_width(length , white_noise) + ASCII_SPACE;
  }
  if (frequency) {
    width[i++] = column_width(length - offset , confidence_limit + offset) + ASCII_SPACE;
    width[i++] = column_width(frequency[offset]) + ASCII_SPACE;
  }

  os << "\n  ";
  for (i = 0;i < nb_curve;i++) {
    if (variable_type[i] == OBSERVED_VALUE) {
      os << " | " << STAT_label[STATL_VARIABLE] << " " << variable1[i];
      if (variable1[i] != variable2[i]) {
        os << " " << STAT_label[STATL_VARIABLE] << " " << variable2[i];
      }
    }

    else {
      switch (variable_type[i]) {
      case OBSERVED_STATE :
        os << " | " << SEQ_label[SEQL_OBSERVED] << " " << STAT_label[STATL_STATE];
        break;
      case THEORETICAL_STATE :
        os << " | " << SEQ_label[SEQL_THEORETICAL] << " " << STAT_label[STATL_STATE];
        break;
      case OBSERVED_OUTPUT :
        os << " | " << SEQ_label[SEQL_OBSERVED] << " " << STAT_label[STATL_OUTPUT];
        break;
      case THEORETICAL_OUTPUT :
        os << " | " << SEQ_label[SEQL_THEORETICAL] << " " << STAT_label[STATL_OUTPUT];
        break;
      }

      os << " " << variable1[i];
    }
  }

  if (white_noise) {
    os << " | " << SEQ_label[SEQL_WHITE_NOISE];
  }
  if (frequency) {
    os << " | " << SEQ_label[SEQL_RANDOMNESS_95_CONFIDENCE_LIMIT]
       << " | " << STAT_label[STATL_FREQUENCY];
  }
  os << endl;

  for (i = offset;i < length;i++) {
    os << setw(width[0]) << i;
    for (j = 0;j < nb_curve;j++) {
      os << setw(width[j + 1]) << point[j][i];
    }

    j = nb_curve + 1;
    if (white_noise) {
      os << setw(width[j++]) << white_noise[i];
    }
    if (frequency) {
      os << setw(width[j++]) << confidence_limit[i];
      os << setw(width[j++]) << frequency[i];
    }
    os << endl;
  }

  if (frequency) {
    delete [] confidence_limit;
  }
  delete [] width;

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Correlation dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Correlation::ascii_write(StatError &error , const char *path , bool exhaustive) const

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
    ascii_write(out_file , exhaustive);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Correlation dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool Correlation::spreadsheet_write(StatError &error , const char *path) const

{
  bool status , autocorrelation , cross_correlation;
  register int i , j;
  double standard_normal_value , confidence_limit;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    switch (type) {
    case SPEARMAN :
      out_file << SEQ_label[SEQL_SPEARMAN] << " ";
      break;
    case KENDALL :
      out_file << SEQ_label[SEQL_KENDALL] << " ";
      break;
    }

    if (offset == 1) {
      out_file << SEQ_label[SEQL_PARTIAL] << " ";
    }

    if ((type == SPEARMAN) || (type == KENDALL)) {
      out_file << SEQ_label[SEQL_RANK] << " ";
    }

    autocorrelation = true;
    cross_correlation = true;
    for (i = 0;i < nb_curve;i++) {
      if (variable_type[i] == OBSERVED_VALUE) {
        if (variable1[i] != variable2[i]) {
          autocorrelation = false;
        }
        else if (variable1[i] == variable2[i]) {
          cross_correlation = false;
        }
      }

      else {
        cross_correlation = false;
      }
    }

    if (autocorrelation) {
      out_file << SEQ_label[SEQL_AUTO];
    }
    else if (cross_correlation) {
      out_file << SEQ_label[SEQL_CROSS];
    }
    out_file << SEQ_label[SEQL_CORRELATION_FUNCTION];

    out_file << "\n";
    for (i = 0;i < nb_curve;i++) {
      if (variable_type[i] == OBSERVED_VALUE) {
        out_file << "\t" << STAT_label[STATL_VARIABLE] << " " << variable1[i];
        if (variable1[i] != variable2[i]) {
          out_file << " " << STAT_label[STATL_VARIABLE] << " " << variable2[i];
        }
      }

      else {
        switch (variable_type[i]) {
        case OBSERVED_STATE :
          out_file << "\t" << SEQ_label[SEQL_OBSERVED] << " " << STAT_label[STATL_STATE];
          break;
        case THEORETICAL_STATE :
          out_file << "\t" << SEQ_label[SEQL_THEORETICAL] << " " << STAT_label[STATL_STATE];
          break;
        case OBSERVED_OUTPUT :
          out_file << "\t" << SEQ_label[SEQL_OBSERVED] << " " << STAT_label[STATL_OUTPUT];
          break;
        case THEORETICAL_OUTPUT :
          out_file << "\t" << SEQ_label[SEQL_THEORETICAL] << " " << STAT_label[STATL_OUTPUT];
          break;
        }

        out_file << " " << variable1[i];
      }
    }

    if (white_noise) {
      out_file << "\t" << SEQ_label[SEQL_WHITE_NOISE];
    }
    if (frequency) {
      out_file << "\t" << SEQ_label[SEQL_RANDOMNESS_95_CONFIDENCE_LIMIT]
               << "\t" << SEQ_label[SEQL_RANDOMNESS_95_CONFIDENCE_LIMIT]
               << "\t" << STAT_label[STATL_FREQUENCY];
    }
    out_file << endl;

    if (frequency) {
      normal dist;
      standard_normal_value = quantile(complement(dist , 0.025));
    }

    for (i = offset;i < length;i++) {
      out_file << i;
      for (j = 0;j < nb_curve;j++) {
        out_file << "\t" << point[j][i];
      }

      if (white_noise) {
        out_file << "\t" << white_noise[i];
      }

      if (frequency) {
        switch (type) {
        case PEARSON :
          confidence_limit = standard_normal_value / sqrt((double)frequency[i]);
          break;
        case SPEARMAN :
          confidence_limit = standard_normal_value / sqrt((double)frequency[i]);
          break;
        case KENDALL :
          confidence_limit = standard_normal_value * sqrt((2 * (2 * (double)frequency[i] + 5)) /
                              (9 * (double)frequency[i] * (double)(frequency[i] - 1)));
          break;
        }

        out_file << "\t" << confidence_limit << "\t" << -confidence_limit
                 << "\t" << frequency[i];
      }

      out_file << endl;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Correlation et des limites de confiance au format Gnuplot.
 *
 *  arguments : path, limites de confiance.
 *
 *--------------------------------------------------------------*/

bool Correlation::plot_print(const char *path , double *confidence_limit) const

{
  bool status = false;
  register int i , j;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    for (i = 0;i < length;i++) {
      for (j = 0;j < nb_curve;j++) {
        out_file << point[j][i] << " ";
      }
      if (white_noise) {
        out_file << white_noise[i] << " ";
      }
      if (frequency) {
        out_file << confidence_limit[i] << " " << -confidence_limit[i] << " "
                 << frequency[i];
      }
      out_file << endl;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Correlation.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Correlation::plot_write(StatError &error , const char *prefix ,
                             const char *title) const

{
  bool status , autocorrelation , cross_correlation;
  register int i , j;
  double standard_normal_value , *confidence_limit = NULL;
  ostringstream data_file_name;


  error.init();

  // calcul des limites de confiance

  if (frequency) {
    normal dist;
    standard_normal_value = quantile(complement(dist , 0.025));

    confidence_limit = new double[length];

    for (i = 0;i < length;i++) {
      switch (type) {
      case PEARSON :
        confidence_limit[i] = standard_normal_value / sqrt((double)frequency[i]);
        break;
      case SPEARMAN :
        confidence_limit[i] = standard_normal_value / sqrt((double)frequency[i]);
        break;
      case KENDALL :
        confidence_limit[i] = standard_normal_value * sqrt((2 * (2 * (double)frequency[i] + 5)) /
                               (9 * (double)frequency[i] * (double)(frequency[i] - 1)));
        break;
      }
    }
  }

  // ecriture du fichier de donnees

  data_file_name << prefix << ".dat";
  status = plot_print((data_file_name.str()).c_str() , confidence_limit);

  if (frequency) {
    delete [] confidence_limit;
  }

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  else {

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
               << "set title" << " \"";
      if (title) {
        out_file << title << " - ";
      }

      switch (type) {
      case SPEARMAN :
        out_file << SEQ_label[SEQL_SPEARMAN] << " ";
        break;
      case KENDALL :
        out_file << SEQ_label[SEQL_KENDALL] << " ";
        break;
      }

      if (offset == 1) {
        out_file << SEQ_label[SEQL_PARTIAL] << " ";
      }

      if ((type == SPEARMAN) || (type == KENDALL)) {
        out_file << SEQ_label[SEQL_RANK] << " ";
      }

      autocorrelation = true;
      cross_correlation = true;
      for (j = 0;j < nb_curve;j++) {
        if (variable_type[j] == OBSERVED_VALUE) {
           if (variable1[j] != variable2[j]) {
            autocorrelation = false;
          }
          else if (variable1[j] == variable2[j]) {
            cross_correlation = false;
          }
        }

        else {
          cross_correlation = false;
        }
      }

      if (autocorrelation) {
        out_file << SEQ_label[SEQL_AUTO];
      }
      else if (cross_correlation) {
        out_file << SEQ_label[SEQL_CROSS];
      }
      out_file << SEQ_label[SEQL_CORRELATION_FUNCTION] << "\"\n\n";

      out_file << "set xlabel \"" << SEQ_label[SEQL_LAG] << "\"" << endl;

      if (length - 1 < TIC_THRESHOLD) {
        out_file << "set xtics 0,1" << endl;
      }

      out_file << "plot [" << offset << ":" << length - 1 << "] [-1:1] ";
      for (j = 0;j < nb_curve;j++) {
        out_file << "\"" << label((data_file_name.str()).c_str()) << "\" using " << j + 1 << " title \"";

        if (variable_type[j] == OBSERVED_VALUE) {
          out_file << STAT_label[STATL_VARIABLE] << " " << variable1[j];
          if (variable1[j] != variable2[j]) {
            out_file << " " << STAT_label[STATL_VARIABLE] << " " << variable2[j];
          }
        }

        else {
          switch (variable_type[j]) {
          case OBSERVED_STATE :
            out_file << SEQ_label[SEQL_OBSERVED] << " " << STAT_label[STATL_STATE];
            break;
          case THEORETICAL_STATE :
            out_file << SEQ_label[SEQL_THEORETICAL] << " " << STAT_label[STATL_STATE];
            break;
          case OBSERVED_OUTPUT :
            out_file << SEQ_label[SEQL_OBSERVED] << " " << STAT_label[STATL_OUTPUT];
            break;
          case THEORETICAL_OUTPUT :
            out_file << SEQ_label[SEQL_THEORETICAL] << " " << STAT_label[STATL_OUTPUT];
            break;
          }

          out_file << " " << variable1[j];
        }

        out_file << "\" with linespoints";
        if (j < nb_curve - 1) {
          out_file << ",\\" << endl;
        }
      }

      j = nb_curve + 1;
      if (white_noise) {
        out_file << ",\\\n\"" << label((data_file_name.str()).c_str()) << "\" using " << j++ << " title \""
                 << SEQ_label[SEQL_WHITE_NOISE] << "\" with linespoints";
      }
      if (frequency) {
        out_file << ",\\\n\"" << label((data_file_name.str()).c_str()) << "\" using " << j++
                 << " notitle with lines";
        out_file << ",\\\n\"" << label((data_file_name.str()).c_str()) << "\" using " << j++
                 << " notitle with lines";
      }
      out_file << endl;

      out_file << "set xlabel" << endl;

      if (length - 1 < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
      }

      if (frequency) {
        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        out_file << "set xlabel \"" << SEQ_label[SEQL_LAG] << "\"" << endl;
        out_file << "set ylabel \"" << SEQ_label[SEQL_PAIR_FREQUENCY] << "\"" << endl;

        if (length - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }
        if (frequency[0] < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [0:" << length - 1 << "] [0:" << frequency[0] << "] \""
                 << label((data_file_name.str()).c_str()) << "\" using " << j++
                 << " notitle with impulses" << endl;

        out_file << "set xlabel" << endl;
        out_file << "set ylabel" << endl;

        if (length - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
        if (frequency[0] < TIC_THRESHOLD) {
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


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet Correlation.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Correlation::get_plotable() const

{
  bool autocorrelation , cross_correlation;
  register int i , j;
  int nb_plot;
  double standard_normal_value , *confidence_limit = NULL;
  ostringstream title , legend;
  MultiPlotSet *plot_set;


  plot_set = new MultiPlotSet(frequency ? 2 : 1);
  MultiPlotSet &plot = *plot_set;

  title.str("");

  switch (type) {
  case SPEARMAN :
    title << SEQ_label[SEQL_SPEARMAN] << " ";
    break;
  case KENDALL :
    title << SEQ_label[SEQL_KENDALL] << " ";
    break;
  }

  if (offset == 1) {
    title << SEQ_label[SEQL_PARTIAL] << " ";
  }

  if ((type == SPEARMAN) || (type == KENDALL)) {
    title << SEQ_label[SEQL_RANK] << " ";
  }

  autocorrelation = true;
  cross_correlation = true;
  for (j = 0;j < nb_curve;j++) {
    if (variable_type[j] == OBSERVED_VALUE) {
      if (variable1[j] != variable2[j]) {
        autocorrelation = false;
      }
      else if (variable1[j] == variable2[j]) {
        cross_correlation = false;
      }
    }

    else {
      cross_correlation = false;
    }
  }

  if (autocorrelation) {
    title << SEQ_label[SEQL_AUTO];
  }
  else if (cross_correlation) {
    title << SEQ_label[SEQL_CROSS];
  }
  title << SEQ_label[SEQL_CORRELATION_FUNCTION];

  plot.title = title.str();

  plot.border = "15 lw 0";

  // 1ere vue : fonctions de correlation

  plot[0].xrange = Range(0 , length - 1);
  plot[0].yrange = Range(-1., 1.);

  plot[0].xlabel = SEQ_label[SEQL_LAG];

  if (length - 1 < TIC_THRESHOLD) {
    plot[0].xtics = 1;
  }

  nb_plot = nb_curve;
  if (white_noise) {
    nb_plot++;
  }

  // calcul des limites de confiance

  if (frequency) {
    normal dist;
    standard_normal_value = quantile(complement(dist , 0.025));

    confidence_limit = new double[length];

    for (i = 0;i < length;i++) {
      switch (type) {
      case PEARSON :
        confidence_limit[i] = standard_normal_value / sqrt((double)frequency[i]);
        break;
      case SPEARMAN :
        confidence_limit[i] = standard_normal_value / sqrt((double)frequency[i]);
        break;
      case KENDALL :
        confidence_limit[i] = standard_normal_value * sqrt((2 * (2 * (double)frequency[i] + 5)) /
                               (9 * (double)frequency[i] * (double)(frequency[i] - 1)));
        break;
      }
    }

    nb_plot += 2;
  }

  plot[0].resize(nb_plot);

  for (i = 0;i < nb_curve;i++) {
    legend.str("");

    if (variable_type[i] == OBSERVED_VALUE) {
      legend << STAT_label[STATL_VARIABLE] << " " << variable1[i];
      if (variable1[i] != variable2[i]) {
        legend << " " << STAT_label[STATL_VARIABLE] << " " << variable2[i];
      }
    }

    else {
      switch (variable_type[i]) {
      case OBSERVED_STATE :
        legend << SEQ_label[SEQL_OBSERVED] << " " << STAT_label[STATL_STATE];
        break;
      case THEORETICAL_STATE :
        legend << SEQ_label[SEQL_THEORETICAL] << " " << STAT_label[STATL_STATE];
        break;
      case OBSERVED_OUTPUT :
        legend << SEQ_label[SEQL_OBSERVED] << " " << STAT_label[STATL_OUTPUT];
        break;
      case THEORETICAL_OUTPUT :
        legend << SEQ_label[SEQL_THEORETICAL] << " " << STAT_label[STATL_OUTPUT];
        break;
      }

      legend << " " << variable1[i];
    }

    plot[0][i].legend = legend.str();

    plot[0][i].style = "linespoints";

    plotable_write(i , plot[0][i]);
  }

  i = nb_curve;
  if (white_noise) {
    plot[0][i].legend = SEQ_label[SEQL_WHITE_NOISE];

    plot[0][i].style = "linespoints";

    for (j = 0;j < length;j++) {
      plot[0][i].add_point(j , white_noise[j]);
    }
    i++;
  }

  if (frequency) {
    plot[0][i].style = "lines";

    for (j = 0;j < length;j++) {
      plot[0][i].add_point(j , confidence_limit[j]);
    }
    i++;

    plot[0][i].style = "lines";

    for (j = 0;j < length;j++) {
      plot[0][i].add_point(j , -confidence_limit[j]);
    }

    // 2eme vue : frequences

    plot[1].xrange = Range(0 , length - 1);
    plot[1].yrange = Range(0 , frequency[0]);

    plot[1].xlabel = SEQ_label[SEQL_LAG];
    plot[1].ylabel = SEQ_label[SEQL_PAIR_FREQUENCY];

    if (length - 1 < TIC_THRESHOLD) {
      plot[1].xtics = 1;
    }
    if (frequency[0] < TIC_THRESHOLD) {
      plot[1].ytics = 1;
    }

    plot[1].resize(1);

    plot[1][0].style = "impulses";

    plotable_frequency_write(plot[1][0]);

    delete [] confidence_limit;
  }

  return plot_set;
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'une fonction de correlation a partir d'un objet Sequences.
 *
 *  arguments : reference sur un objet Correlation, indices des variables,
 *              normalisation (APPROXIMATED/EXACT), calcul des moyennes par individu ou non.
 *
 *--------------------------------------------------------------*/

void Sequences::correlation_computation(Correlation &correl , int variable1 , int variable2 ,
                                        int normalization , bool individual_mean) const

{
  if (correl.type == PEARSON) {
    register int i , j , k;
    int max_lag = correl.length - 1 , *pisequence1 , *pisequence2 , *pfrequency;
    double *prsequence1 , *prsequence2;
    double variance1 , variance2 , diff , norm , *mean1 , *mean2 , *ppoint;


    // calcul des moyennes et des variances

/*    mean1 = mean_computation(variable1);

    if (variable1 == variable2) {
      mean2 = mean1;
      norm = variance_computation(variable1 , mean1);
    }
    else {
      mean2 = mean_computation(variable2);
      norm = sqrt(variance_computation(variable1 , mean1) *
                  variance_computation(variable2 , mean2));
    }

    norm *= (cumul_length - 1); */

    mean1 = new double[nb_sequence];
    mean2 = new double[nb_sequence];

    if (individual_mean) {
      variance1 = 0.;

      if (type[variable1] != REAL_VALUE) {
        for (i = 0;i < nb_sequence;i++) {
          pisequence1 = int_sequence[i][variable1];
          mean1[i] = 0.;
          for (j = 0;j < length[i];j++) {
            mean1[i] += *pisequence1++;
          }
          mean1[i] /= length[i];

          pisequence1 = int_sequence[i][variable1];
          for (j = 0;j < length[i];j++) {
            diff = *pisequence1++ - mean1[i];
            variance1 += diff * diff;
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          prsequence1 = real_sequence[i][variable1];
          mean1[i] = 0.;
          for (j = 0;j < length[i];j++) {
            mean1[i] += *prsequence1++;
          }
          mean1[i] /= length[i];

          prsequence1 = real_sequence[i][variable1];
          for (j = 0;j < length[i];j++) {
            diff = *prsequence1++ - mean1[i];
            variance1 += diff * diff;
          }
        }
      }

      if (variable1 == variable2) {
        for (i = 0;i < nb_sequence;i++) {
          mean2[i] = mean1[i];
        }
        norm = variance1;
      }

      else {
        variance2 = 0.;

        if (type[variable2] != REAL_VALUE) {
          for (i = 0;i < nb_sequence;i++) {
            pisequence2 = int_sequence[i][variable2];
            mean2[i] = 0.;
            for (j = 0;j < length[i];j++) {
              mean2[i] += *pisequence2++;
            }
            mean2[i] /= length[i];

            pisequence2 = int_sequence[i][variable2];
            for (j = 0;j < length[i];j++) {
              diff = *pisequence2++ - mean2[i];
              variance2 += diff * diff;
            }
          }
        }

        else {
          for (i = 0;i < nb_sequence;i++) {
            prsequence2 = real_sequence[i][variable2];
            mean2[i] = 0.;
            for (j = 0;j < length[i];j++) {
              mean2[i] += *prsequence2++;
            }
            mean2[i] /= length[i];

            prsequence2 = real_sequence[i][variable2];
            for (j = 0;j < length[i];j++) {
              diff = *prsequence2++ - mean2[i];
              variance2 += diff * diff;
            }
          }
        }

        norm = sqrt(variance1 * variance2);
      }
    }

    else {
      mean1[0] = mean_computation(variable1);
      for (i = 1;i < nb_sequence;i++) {
        mean1[i] = mean1[0];
      }

      if (variable1 == variable2) {
        for (i = 0;i < nb_sequence;i++) {
          mean2[i] = mean1[i];
        }
        norm = variance_computation(variable1 , mean1[0]);
      }

      else {
        mean2[0] = mean_computation(variable2);
        for (i = 1;i < nb_sequence;i++) {
          mean2[i] = mean2[0];
        }
        norm = sqrt(variance_computation(variable1 , mean1[0]) *
                    variance_computation(variable2 , mean2[0]));
      }

      norm *= (cumul_length - 1);
    }

    // calcul des coefficients de correlation

    ppoint = correl.point[0];
    pfrequency = correl.frequency;

    for (i = 0;i <= max_lag;i++) {
      *ppoint = 0.;
      *pfrequency = 0;

      for (j = 0;j < nb_sequence;j++) {
        if (length[j] > i) {
          if (type[variable1] != REAL_VALUE) {
            pisequence1 = int_sequence[j][variable1];
          }
          else {
            prsequence1 = real_sequence[j][variable1];
          }

          if (type[variable2] != REAL_VALUE) {
            pisequence2 = int_sequence[j][variable2] + i;
          }
          else {
            prsequence2 = real_sequence[j][variable2] + i;
          }

          if ((type[variable1] != REAL_VALUE) && (type[variable2] != REAL_VALUE)) {
            for (k = 0;k < length[j] - i;k++) {
              *ppoint += (*pisequence1++ - mean1[j]) * (*pisequence2++ - mean2[j]);
            }
          }
          else if ((type[variable1] != REAL_VALUE) && (type[variable2] == REAL_VALUE)) {
            for (k = 0;k < length[j] - i;k++) {
              *ppoint += (*pisequence1++ - mean1[j]) * (*prsequence2++ - mean2[j]);
            }
          }
          else if ((type[variable1] == REAL_VALUE) && (type[variable2] != REAL_VALUE)) {
            for (k = 0;k < length[j] - i;k++) {
              *ppoint += (*prsequence1++ - mean1[j]) * (*pisequence2++ - mean2[j]);
            }
          }
//          else if ((type[variable1] == REAL_VALUE) && (type[variable2] == REAL_VALUE)) {
          else {
            for (k = 0;k < length[j] - i;k++) {
              *ppoint += (*prsequence1++ - mean1[j]) * (*prsequence2++ - mean2[j]);
            }
          }
     
          *pfrequency += length[j] - i;
        }
      }

      switch (normalization) {
      case APPROXIMATED :
        *ppoint++ /= norm;
        break;
      case EXACT :
        *ppoint++ *= cumul_length / (*pfrequency * norm);
        break;
      }

//      if (*pfrequency++ <= CORRELATION_MIN_FREQUENCY) {
      if (*pfrequency++ <= cumul_length * FREQUENCY_RATIO) {
        correl.length = i + 1;
        break;
      }
    }

    if (normalization == APPROXIMATED) {
      pfrequency = correl.frequency;
      for (i = 0;i < correl.length;i++) {
        *pfrequency++ = cumul_length;
      }
    }

    delete [] mean1;
    delete [] mean2;
  }

  else if (correl.type == SPEARMAN) {
    register int i , j , k;
    int max_lag = correl.length - 1 , *pisequence1 , *pisequence2 , *pfrequency;
    double main_term , correction , norm , rank_mean , *ppoint , *rank[2];


    // calcul du terme principal et des termes de correction pour les ex-aequo

    main_term = cumul_length * ((double)cumul_length * (double)cumul_length - 1);

    pfrequency = marginal_distribution[variable1]->frequency + marginal_distribution[variable1]->offset;
    correction = 0.;
    for (i = marginal_distribution[variable1]->offset;i < marginal_distribution[variable1]->nb_value;i++) {
      if (*pfrequency > 1) {
        correction += *pfrequency * ((double)*pfrequency * (double)*pfrequency - 1);
      }
      pfrequency++;
    }

    if (variable1 == variable2) {
      norm = (main_term - correction) / 12.;
    }

    else {
      norm = sqrt(main_term - correction) / 12.;

      pfrequency = marginal_distribution[variable2]->frequency + marginal_distribution[variable2]->offset;
      correction = 0.;
      for (i = marginal_distribution[variable2]->offset;i < marginal_distribution[variable2]->nb_value;i++) {
        if (*pfrequency > 1) {
          correction += *pfrequency * ((double)*pfrequency * (double)*pfrequency - 1);
        }
        pfrequency++;
      }

      norm *= sqrt(main_term - correction);
    }

    // calcul des rangs

    rank_mean = (double)(cumul_length + 1) / 2.;

    rank[0] = marginal_distribution[variable1]->rank_computation();

    if (variable1 == variable2) {
      rank[1] = new double[marginal_distribution[variable1]->nb_value];

      for (i = marginal_distribution[variable1]->offset;i < marginal_distribution[variable1]->nb_value;i++) {
        rank[1][i] = rank[0][i];
      }
    }

    else {
      rank[1] = marginal_distribution[variable2]->rank_computation();
    }

    // calcul des coefficients de correlation

    ppoint = correl.point[0];
    pfrequency = correl.frequency;

    for (i = 0;i <= max_lag;i++) {
      *ppoint = 0.;
      *pfrequency = 0;

      // calcul des differences de rangs centrees

      for (j = 0;j < nb_sequence;j++) {
        if (length[j] > i) {
          pisequence1 = int_sequence[j][variable1];
          pisequence2 = int_sequence[j][variable2] + i;
          for (k = 0;k < length[j] - i;k++) {
            *ppoint += (rank[0][*pisequence1++] - rank_mean) * (rank[1][*pisequence2++] - rank_mean);
          }
          *pfrequency += length[j] - i;
        }
      }

      switch (normalization) {
      case APPROXIMATED :
        *ppoint++ /= norm;
        break;
      case EXACT :
        *ppoint++ *= cumul_length / (*pfrequency * norm);
        break;
      }

//      if (*pfrequency++ <= CORRELATION_MIN_FREQUENCY) {
      if (*pfrequency++ <= cumul_length * FREQUENCY_RATIO) {
        correl.length = i + 1;
        break;
      }
    }

    if (normalization == APPROXIMATED) {
      pfrequency = correl.frequency;
      for (i = 0;i < correl.length;i++) {
        *pfrequency++ = cumul_length;
      }
    }

    for (i = 0;i < 2;i++) {
      delete [] rank[i];
    }
  }

  else {
    register int i , j , k;
    int max_lag = correl.length - 1 , nb_vector , *pisequence1 , *pisequence2 ,
        *pfrequency , **int_vector;
    double *ppoint;
    Vectors *vec;


    int_vector = new int*[cumul_length];
    for (i = 0;i < cumul_length;i++) {
      int_vector[i] = new int[2];
    }

    ppoint = correl.point[0];
    pfrequency = correl.frequency;

    for (i = 0;i <= max_lag;i++) {

      // constitution de l'echantillon de vecteurs

      nb_vector = 0;
      for (j = 0;j < nb_sequence;j++) {
        if (length[j] > i) {
          pisequence1 = int_sequence[j][variable1];
          pisequence2 = int_sequence[j][variable2] + i;
          for (k = 0;k < length[j] - i;k++) {
            int_vector[nb_vector][0] = *pisequence1++;
            int_vector[nb_vector][1] = *pisequence2++;
            nb_vector++;
          }
        }
      }

      *pfrequency = nb_vector;

      vec = new Vectors(nb_vector , NULL , 2 , int_vector);

      // calcul du coefficient de correlation

      switch (correl.type) {
      case SPEARMAN2 :
        *ppoint++ = vec->spearman_rank_single_correlation_computation();
        break;
      case KENDALL :
        *ppoint++ = vec->kendall_rank_single_correlation_computation();
        break;
      }

      delete vec;

//      if (*pfrequency++ <= CORRELATION_MIN_FREQUENCY) {
      if (*pfrequency++ <= cumul_length * FREQUENCY_RATIO) {
        correl.length = i + 1;
        break;
      }
    }

    for (i = 0;i < cumul_length;i++) {
      delete [] int_vector[i];
    }
    delete [] int_vector;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'une fonction de correlation a partir d'un objet Sequences.
 *
 *  arguments : reference sur un objet StatError, indices des variables,
 *              type de coefficient (PEARSON/SPEARMAN/KENDALL), decalage maximum,
 *              normalisation (APPROXIMATED/EXACT), calcul des moyennes par individu ou non.
 *
 *--------------------------------------------------------------*/

Correlation* Sequences::correlation_computation(StatError &error , int variable1 ,
                                                int variable2 , int itype , int max_lag ,
                                                int normalization , bool individual_mean) const

{
  bool status = true;
  Correlation *correl;


  correl = NULL;
  error.init();

  if ((variable1 < 1) || (variable1 > nb_variable)) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_VARIABLE] << " " << variable1 << ": "
                  << STAT_error[STATR_VARIABLE_INDEX];
    error.update((error_message.str()).c_str());
  }

  else {
    variable1--;

    if ((itype == PEARSON) && (type[variable1] != INT_VALUE) && (type[variable1] != STATE) &&
        (type[variable1] != REAL_VALUE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable1 + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }

//    if (((itype == SPEARMAN) || (itype == KENDALL)) && (!marginal_distribution[variable1])) {
    if (((itype == SPEARMAN) || (itype == SPEARMAN2) || (itype == KENDALL)) &&
        (!marginal_distribution[variable1])) {
      status = false;
      ostringstream error_message;
      error_message << STAT_error[STATR_RANK_CORRELATION_COMPUTATION] << ": "
                    << STAT_label[STATL_VARIABLE] << " " << variable1 + 1 << " "
                    << STAT_error[STATR_SHIFTED_SCALED];
      error.update((error_message.str()).c_str());
    }
  }

  if ((variable2 < 1) || (variable2 > nb_variable)) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_VARIABLE] << " " << variable2 << ": "
                  << STAT_error[STATR_VARIABLE_INDEX];
    error.update((error_message.str()).c_str());
  }

  else {
    variable2--;

    if ((itype == PEARSON) && (type[variable2] != INT_VALUE) && (type[variable2] != STATE) &&
        (type[variable2] != REAL_VALUE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable2 + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }

//    if (((itype == SPEARMAN) || (itype == KENDALL)) && (!marginal_distribution[variable2])) {
    if (((itype == SPEARMAN) || (itype == SPEARMAN2) || (itype == KENDALL)) &&
        (!marginal_distribution[variable2])) {
      status = false;
      ostringstream error_message;
      error_message << STAT_error[STATR_RANK_CORRELATION_COMPUTATION] << ": "
                    << STAT_label[STATL_VARIABLE] << " " << variable2 + 1 << " "
                    << STAT_error[STATR_SHIFTED_SCALED];
      error.update((error_message.str()).c_str());
    }
  }

  if ((max_lag < I_DEFAULT) || (max_lag >= max_length)) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_LAG]);
  }

  if (status) {
    if (max_lag == I_DEFAULT) {
      max_lag = max_length - 1;
    }

    // creation de la fonction de correlation

//    correl = new Correlation(itype , max_lag , variable1 + 1 , variable2 + 1);
    correl = new Correlation((itype == SPEARMAN2 ? SPEARMAN : itype) , max_lag , variable1 + 1 , variable2 + 1);
    correlation_computation(*correl , variable1 , variable2 , normalization , individual_mean);
  }

  return correl;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la fonction de correlation theorique d'un bruit blanc
 *  pour une filtre donne.
 *
 *  arguments : reference sur un objet StatError, taille du filtre, filtre,
 *              flag calcul du filtre correspondant aux residus.
 *
 *--------------------------------------------------------------*/

bool Correlation::white_noise_correlation(StatError &error , int nb_point , double *filter ,
                                          int residual)

{
  bool status = true;
  register int i , j;
  double variance;


  error.init();

  if ((type != PEARSON) || (offset != 0)) {
    status = false;
    ostringstream correction_message;
    correction_message << SEQ_label[SEQL_PEARSON] << " "
                       << SEQ_label[SEQL_CORRELATION_FUNCTION];
    error.correction_update(SEQ_error[SEQR_CORRELATION_COEFF_TYPE] , (correction_message.str()).c_str());
  }

  for (i = 1;i < nb_curve;i++) {
    if ((point[i][0] < point[0][0] - DOUBLE_ERROR) || (point[i][0] > point[0][0] + DOUBLE_ERROR)) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_CORRELATION_FUNCTION] << " " << i + 1  << " "
                    << SEQ_error[SEQR_INCOMPATIBLE_CORRELATION_FUNCTION];
      error.update((error_message.str()).c_str());
    }
  }

  if (status) {
    if (residual) {
      for (i = 0;i < nb_point;i++) {
        filter[i] = -filter[i];
      }
      filter[nb_point / 2]++;
    }

    delete [] white_noise;
    white_noise = new double[length];

    variance = 0.;
    for (i = 0;i < nb_point;i++) {
      variance += filter[i] * filter[i];
    }

    white_noise[0] = point[0][0];
    for (i = 1;i < MIN(nb_point , length);i++) {
      white_noise[i] = 0.;
      for (j = 0;j < nb_point - i;j++) {
        white_noise[i] += filter[i + j] * filter[j];
      }
      white_noise[i] = white_noise[i] * point[0][0] / variance;
    }

    for (i = nb_point;i < length;i++) {
      white_noise[i] = 0.;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la fonction de correlation theorique d'un bruit blanc
 *  pour une filtre donne.
 *
 *  arguments : reference sur un objet StatError, loi symmetrique.
 *
 *--------------------------------------------------------------*/

bool Correlation::white_noise_correlation(StatError &error , const Distribution &dist)

{
  bool status = true;


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
    status = white_noise_correlation(error , dist.nb_value , dist.mass);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la fonction de correlation theorique d'un bruit blanc
 *  pour une differenciation.
 *
 *  arguments : reference sur un objet StatError, ordre de la differenciation.
 *
 *--------------------------------------------------------------*/

bool Correlation::white_noise_correlation(StatError &error , int order)

{
  bool status = true;
  register int i;
  double *filter;


  error.init();

  if ((order < 1) || (order > MAX_DIFFERENCING_ORDER)) {
    status = false;
    error.update(SEQ_error[SEQR_DIFFERENCING_ORDER]);
  }

  if (status) {
    filter = new double[order + 1];

    filter[0] = 1.;
    for (i = 1;i <= order;i++) {
      filter[i] = -filter[i - 1] * (order - i + 1) / i;
    }

#   ifdef DEBUG
    cout << "\nfiltre : ";
    for (i = 0;i <= order;i++) {
      cout << filter[i] << " ";
    }
    cout << endl;
#   endif

    status = white_noise_correlation(error , order + 1 , filter , false);

    delete [] filter;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'une fonction d'autocorrelation partielle a partir d'un objet Sequences.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              type de coefficient (PEARSON/KENDALL), decalage maximum.
 *
 *--------------------------------------------------------------*/

Correlation* Sequences::partial_autocorrelation_computation(StatError &error , int variable ,
                                                            int itype , int max_lag) const

{
  bool status = true;
  register int i , j;
  int *pfrequency , *cfrequency;
  double sum , denom , *ppoint , *cpoint1 , *cpoint2 , *aux_correl , *paux_correl ,
         *aux , *paux1 , *paux2;
  Correlation *correl , *partial_correl;


  partial_correl = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((itype == PEARSON) && (type[variable] != INT_VALUE) && (type[variable] != STATE) &&
        (type[variable] != REAL_VALUE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }

    if ((itype == KENDALL) && (!marginal_distribution[variable])) {
      status = false;
      ostringstream error_message;
      error_message << STAT_error[STATR_RANK_CORRELATION_COMPUTATION] << ": "
                    << STAT_label[STATL_VARIABLE] << " " << variable + 1 << " "
                    << STAT_error[STATR_SHIFTED_SCALED];
      error.update((error_message.str()).c_str());
    }
  }

  if ((max_lag != I_DEFAULT) && ((max_lag < 1) || (max_lag >= max_length))) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_LAG]);
  }

  if (status) {
    if (max_lag == I_DEFAULT) {
      max_lag = max_length - 1;
    }

    // construction de la fonction de d'autocorrelation

    correl = correlation_computation(error , variable + 1 , variable + 1 , itype ,
                                     max_lag , APPROXIMATED);
    max_lag = correl->length - 1;

    // creation de la fonction d'autocorrelation partielle

    partial_correl = new Correlation(itype , max_lag , variable + 1 , variable + 1);
    partial_correl->offset = 1;

    // calcul des coefficients de correlation partiels

    paux_correl = new double[max_lag + 1];
    aux_correl = new double[max_lag];

    ppoint = partial_correl->point[0];
    cpoint2 = correl->point[0] + 1;
    pfrequency = partial_correl->frequency;
    cfrequency = correl->frequency;

    *ppoint++ = 0.;
    *pfrequency++ = *cfrequency++;

    denom = 1.;
    for (i = 1;i <= max_lag;i++) {
      *pfrequency++ = *cfrequency++;

      cpoint1 = correl->point[0] + i;
//      cpoint2 = correl->point[0] + 1;
      paux1 = paux_correl + 1;
      sum = 0.;
//      double sum2 = 0.;
      for (j = 1;j < i;j++) {
//        sum2 += *paux1 * *cpoint2++;
        sum += *paux1++ * *--cpoint1;
      }

//      *ppoint = (*cpoint2 - sum) / (1. - sum2);

      *ppoint = (*cpoint2++ - sum) / denom;
      denom *= (1. - *ppoint * *ppoint);

      aux = aux_correl + 1;
      paux1 = paux_correl + 1;
      paux2 = paux_correl + i;
      for (j = 1;j < i;j++) {
        *aux++ = *paux1++ - *ppoint * *--paux2;
      }

      paux1 = paux_correl + 1;
      aux = aux_correl + 1;
      for (j = 1;j < i;j++) {
        *paux1++ = *aux++;
      }

      *paux1 = *ppoint++;
    }

    delete correl;
    delete [] paux_correl;
    delete [] aux_correl;
  }

  return partial_correl;
}


};  // namespace sequence_analysis
