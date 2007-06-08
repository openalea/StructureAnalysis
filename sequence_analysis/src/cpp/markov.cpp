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
#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"
// #include <rw/vstream.h>
// #include <rw/rwfile.h>
#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"
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
extern char* label(const char *file_name);



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Function.
 *
 *--------------------------------------------------------------*/

Function::Function()

{
  residual = 0;
  frequency = 0;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Function.
 *
 *  arguments : identificateur, longueur des courbes, parametres.
 *
 *--------------------------------------------------------------*/

Function::Function(int iident , int length , double *iparameter)
:Regression_kernel(iident , 0 , length - 1)

{
  register int i;


  for (i = 0;i < nb_parameter;i++) {
    parameter[i] = iparameter[i];
  }

  residual = 0;
  frequency = 0;

  computation();
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Function.
 *
 *  arguments : identificateur, longueur des courbes.
 *
 *--------------------------------------------------------------*/

Function::Function(int iident , int length)
:Regression_kernel(iident , 0 , length - 1)

{
  register int i;


  residual = new double[max_value + 1];
  for (i = 0;i <= max_value;i++) {
    residual[i] = -D_INF;
  }

  frequency = new int[max_value + 1];
  for (i = 0;i <= max_value;i++) {
    frequency[i] = 0;
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Function.
 *
 *  argument : reference sur un objet Function.
 *
 *--------------------------------------------------------------*/

void Function::copy(const Function &function)

{
  if ((function.residual) && (function.frequency)) {
    register int i;


    residual = new double[max_value + 1];
    for (i = 0;i <= max_value;i++) {
      residual[i] = function.residual[i];
    }

    frequency = new int[max_value + 1];
    for (i = 0;i <= max_value;i++) {
      frequency[i] = function.frequency[i];
    }
  }

  else {
    residual = 0;
    frequency = 0;
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe Function.
 *
 *  argument : reference sur un objet Function.
 *
 *--------------------------------------------------------------*/

Function::Function(const Function &function)

{
  Regression_kernel::copy(function);
  copy(function);
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Function.
 *
 *--------------------------------------------------------------*/

void Function::remove()

{
  delete [] residual;
  delete [] frequency;
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Function.
 *
 *--------------------------------------------------------------*/

Function::~Function()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Function.
 *
 *  argument : reference sur un objet Function.
 *
 *--------------------------------------------------------------*/

Function& Function::operator=(const Function &function)

{
  if (&function != this) {
    remove();
    Regression_kernel::remove();

    Regression_kernel::copy(function);
    copy(function);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Analyse du format d'un objet Function.
 *
 *  arguments : reference sur un objet Format_error, stream,
 *              reference sur l'indice de la ligne lue, longueur,
 *              bornes sur les valeurs de la fonction.
 *
 *--------------------------------------------------------------*/

Function* function_parsing(Format_error &error , ifstream &in_file , int &line ,
                           int length , double min , double max)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status = true , lstatus;
  register int i , j;
  int ident = I_DEFAULT , nb_parameter = 0;
  long index;
  double parameter[3];
  Function *function;


  function = 0;

  while (buffer.readLine(in_file , false)) {
    line++;

#   ifdef DEBUG
    cout << line << "  " << buffer << endl;
#   endif

    position = buffer.first('#');
    if (position != RW_NPOS) {
      buffer.remove(position);
    }
    i = 0;

    RWCTokenizer next(buffer);

    while (!((token = next()).isNull())) {
      if (i <= 1) {
        switch (i) {

        // test mot cle MONOMOLECULAR / LOGISTIC

        case 0 : {
          for (j = 1;j < 3;j++) {
            if (token == STAT_function_word[j]) {
              ident = j;
              break;
            }
          }

          if (j == 3) {
            status = false;
            error.update(STAT_parsing[STATP_KEY_WORD] , line , i + 1);
          }
          else {
            nb_parameter = 3;
            parameter[0] = D_DEFAULT;
          }
          break;
        }

        // test mot cle FUNCTION

        case 1 : {
          if (token != STAT_word[STATW_FUNCTION]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_FUNCTION] , line , i + 1);
          }
          break;
        }
        }
      }

      else {
        switch ((i - 2) % 4) {

        // test mot cle PARAMETER

        case 0 : {
          if (token != STAT_word[STATW_PARAMETER]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_PARAMETER] , line , i + 1);
          }
          break;
        }

        // test indice du parametre

        case 1 : {
          lstatus = locale.stringToNum(token , &index);
          if ((lstatus) && (index != (i - 2) / 4 + 1)) {
            lstatus = false;
          }

          if (!lstatus) {
            status = false;
            error.correction_update(STAT_parsing[STATP_PARAMETER_INDEX] , (i - 2) / 4 + 1 , line , i + 1);
          }
          break;
        }

        // test separateur

        case 2 : {
          if (token != ":") {
            status = false;
            error.update(STAT_parsing[STATP_SEPARATOR] , line , i + 1);
          }
          break;
        }

        // test valeur du parametre

        case 3 : {
          if ((i - 2) / 4 < nb_parameter) {
            lstatus = locale.stringToNum(token , parameter + (i - 2) / 4);

            if (lstatus) {
              switch (ident) {

              case STAT_MONOMOLECULAR : {
                switch ((i - 2) / 4) {

                case 0 : {
                  if ((parameter[0] < min) || (parameter[0] > max)) {
                    lstatus = false;
                  }
                  break;
                }

                case 1 : {
                  if ((parameter[0] != D_DEFAULT) &&
                      ((parameter[0] + parameter[1] < min) || (parameter[0] + parameter[1] > max))) {
                    lstatus = false;
                  }
                  break;
                }

                case 2 : {
                  if (parameter[2] <= 0.) {
                    lstatus = false;
                  }
                  break;
                }
                }

                break;
              }

              case STAT_LOGISTIC : {
                switch ((i - 2) / 4) {

                case 0 : {
                  if ((parameter[0] < min) || (parameter[0] > max)) {
                    lstatus = false;
                  }
                  break;
                }

                case 1 : {
                  if ((parameter[0] != D_DEFAULT) &&
                      ((parameter[0] / (1. + parameter[1]) < min) || (parameter[0] / (1. + parameter[1]) > max))) {
                    lstatus = false;
                  }
                  break;
                }

                case 2 : {
                  if (parameter[2] <= 0.) {
                    lstatus = false;
                  }
                  break;
                }
                }

                break;
              }
              }
            }

            if (!lstatus) {
              status = false;
              error.update(STAT_parsing[STATP_PARAMETER_VALUE] , line , i + 1);
            }
          }
          break;
        }
        }
      }

      i++;
    }

    if (i > 0) {
      if (i != 14) {
        status = false;
        error.update(STAT_parsing[STATP_FORMAT] , line);
      }

      break;
    }
  }

  if (ident == I_DEFAULT) {
    status = false;
    error.update(STAT_parsing[STATP_FORMAT] , line);
  }

  if (status) {
    function = new Function(ident , length , parameter);
  }

  return function;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une fonction et d'une courbe au format ASCII.
 *
 *  arguments : stream, flag niveau de detail, flag fichier,
 *              pointeur sur un objet Curves.
 *
 *--------------------------------------------------------------*/

ostream& Function::ascii_print(ostream &os , bool exhaustive , bool file_flag ,
                               const Curves *curves) const

{
  register int i;
  int *pfrequency , width[6];
  double self_transition_mean , residual_mean , residual_standard_deviation ,
         *standard_residual , *presidual , *pstandard_residual , *ppoint , *cpoint ,
         square_sum[3];
  long old_adjust;


  old_adjust = os.setf(ios::right , ios::adjustfield);

  ascii_parameter_print(os);
  os << endl;

  if (curves) {
    self_transition_mean = curves->mean_computation(0);
    square_sum[0] = regression_square_sum_computation(self_transition_mean);
    square_sum[1] = residual_square_sum_computation();
    square_sum[2] = curves->total_square_sum_computation(0 , self_transition_mean);

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_REGRESSION_VARIATION_TOTAL_VARIATION] << ": "
       << 1. - square_sum[1] / square_sum[2] << endl;

    if (file_flag) {
      os << "# ";
    }
    os << regression_df << " " << STAT_label[STATL_REGRESSION] << " " << STAT_label[STATL_FREEDOM_DEGREES] << "   "
       << residual_df << " " << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_FREEDOM_DEGREES] << endl;

#   ifdef DEBUG
    os << "\n" << STAT_label[STATL_REGRESSION] << " " << STAT_label[STATL_SQUARE_SUM] << ": " << square_sum[0]
       << "   " << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_SQUARE_SUM] << ": " << square_sum[1]
       << "   " << STAT_label[STATL_TOTAL] << " " << STAT_label[STATL_SQUARE_SUM] << ": " << square_sum[2] << endl;
#   endif

    // ecriture moyenne et ecart-type des residus

    residual_mean = residual_mean_computation();
    residual_standard_deviation = sqrt(residual_variance_computation(residual_mean));

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_MEAN] << ": " << residual_mean << "   "
       << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_STANDARD_DEVIATION] << ": "
       << residual_standard_deviation << endl;
  }

  if (exhaustive) {
    if (curves) {

      // calcul des residus reduits

      standard_residual = new double[max_value + 1];

      pfrequency = frequency;
      pstandard_residual = standard_residual;
      presidual = residual;

      for (i = 0;i <= max_value;i++) {
        if (*pfrequency++ > 0) {
          *pstandard_residual = *presidual / residual_standard_deviation;
        }
        pstandard_residual++;
        presidual++;
      }
    }

    // calcul des largeurs des colonnes

    width[0] = column_width(max_value);
    width[2] = column_width(max_value + 1 , point) + ASCII_SPACE;
    if (curves) {
      width[1] = column_width(curves->length , curves->point[0]) + ASCII_SPACE;
      width[3] = column_width(max_value + 1 , residual) + ASCII_SPACE;
      width[4] = column_width(max_value + 1 , standard_residual) + ASCII_SPACE;
      width[5] = column_width(curves->max_frequency_computation()) + ASCII_SPACE;
    }

    // ecriture des reponses observees et theoriques, des residus, des residus reduits et
    // des frequences

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "  ";
    if (curves) {
      os << " | " << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION];
    }
    os << " | " << SEQ_label[SEQL_THEORETICAL] << " " << SEQ_label[SEQL_SELF_TRANSITION];
    if (curves) {
      os << " | " << STAT_label[STATL_RESIDUAL] << " | " << STAT_label[STATL_STANDARDIZED_RESIDUAL]
         << " | " << STAT_label[STATL_FREQUENCY];
    }
    os << endl;

    ppoint = point;
    if (curves) {
      cpoint = curves->point[0];
      presidual = residual;
      pstandard_residual = standard_residual;
      pfrequency = frequency;
    }

    for (i = 0;i <= max_value;i++) {
      if (file_flag) {
        os << "# ";
      }
      os << setw(width[0]) << i;

      if (curves) {
        if (*pfrequency > 0) {
          os << setw(width[1]) << *cpoint;
        }
        else {
          os << setw(width[1]) << " ";
        }
        cpoint++;
      }

      os << setw(width[2]) << *ppoint++;

      if (curves) {
        if (*pfrequency > 0) {
          os << setw(width[3]) << *presidual;
          os << setw(width[4]) << *pstandard_residual;
        }
        else {
          os << setw(width[3]) << " ";
          os << setw(width[4]) << " ";
        }
        presidual++;
        pstandard_residual++;

        os << setw(width[5]) << *pfrequency++;
      }

      os << endl;
    }

    if (curves) {
      delete [] standard_residual;
    }
  }

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une fonction et d'une courbe au format tableur.
 *
 *  arguments : stream, pointeur sur un objet Curves.
 *
 *--------------------------------------------------------------*/

ostream& Function::spreadsheet_print(ostream &os , const Curves *curves) const

{
  register int i;
  int *pfrequency;
  double self_transition_mean , residual_mean , residual_standard_deviation , *presidual ,
         *ppoint , *cpoint , square_sum[3];


  os << STAT_function_word[ident] << " " << STAT_word[STATW_FUNCTION];
  for (i = 0;i < nb_parameter;i++) {
    os << "\t\t" << STAT_word[STATW_PARAMETER] << " " << i + 1 << "\t" << parameter[i];
  }
  os << endl;

  if (curves) {
    self_transition_mean = curves->mean_computation(0);
    square_sum[0] = regression_square_sum_computation(self_transition_mean);
    square_sum[1] = residual_square_sum_computation();
    square_sum[2] = curves->total_square_sum_computation(0 , self_transition_mean);

    os << "\n" << STAT_label[STATL_REGRESSION_VARIATION_TOTAL_VARIATION] << "\t"
       << 1. - square_sum[1] / square_sum[2] << endl;

    os << regression_df << "\t" << STAT_label[STATL_REGRESSION] << " " << STAT_label[STATL_FREEDOM_DEGREES] << "\t\t"
       << residual_df << "\t" << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_FREEDOM_DEGREES] << endl;

#   ifdef DEBUG
    os << "\n" << STAT_label[STATL_REGRESSION] << " " << STAT_label[STATL_SQUARE_SUM] << "\t" << square_sum[0]
       << "\t\t" << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_SQUARE_SUM] << "\t" << square_sum[1]
       << "\t\t" << STAT_label[STATL_TOTAL] << " " << STAT_label[STATL_SQUARE_SUM] << "\t" << square_sum[2] << endl;
#   endif

    // ecriture moyenne et ecart-type des residus

    residual_mean = residual_mean_computation();
    residual_standard_deviation = sqrt(residual_variance_computation(residual_mean));

    os << "\n" << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_MEAN] << "\t" << residual_mean
       << "\t\t" << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_STANDARD_DEVIATION] << "\t"
       << residual_standard_deviation << endl;
  }

  // ecriture des reponses observees et theoriques, des residus des residus reduits et
  // des frequences

  os << "\n";
  if (curves) {
    os << "\t" << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION];
  }
  os << "\t" << SEQ_label[SEQL_THEORETICAL] << " " << SEQ_label[SEQL_SELF_TRANSITION];
  if (curves) {
    os << "\t" << STAT_label[STATL_RESIDUAL] << "\t" << STAT_label[STATL_STANDARDIZED_RESIDUAL]
       << "\t" << SEQ_label[SEQL_ASYMPTOTE] << "\t" << SEQ_label[SEQL_ASYMPTOTE]
       << "\t" << STAT_label[STATL_FREQUENCY];
  }
  os << endl;

  ppoint = point;
  if (curves) {
    cpoint = curves->point[0];
    presidual = residual;
    pfrequency = frequency;
  }

  for (i = 0;i <= max_value;i++) {
    os << i;

    if (curves) {
      os << "\t";
      if (*pfrequency > 0) {
        os << *cpoint;
      }
      cpoint++;
    }

    os << "\t" << *ppoint;

    if (curves) {
      if (*pfrequency > 0) {
        os << "\t" << *presidual << "\t" << *presidual / residual_standard_deviation;
      }
      else {
        os << "\t\t";
      }
      presidual++;

      os << "\t" << (1. - *ppoint) / residual_standard_deviation
         << "\t" << -*ppoint / residual_standard_deviation << "\t" << *pfrequency++;
    }
    ppoint++;

    os << endl;
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une fonction et des bornes sur les residus standardises
 *  au format Gnuplot.
 *
 *  arguments : path, ecart-type des residus.
 *
 *--------------------------------------------------------------*/

bool Function::plot_print(const char *path , double residual_standard_deviation) const

{
  bool status = false;
  register int i;
  double *ppoint;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    ppoint = point;
    for (i = min_value;i <= max_value;i++) {
      out_file << *ppoint;
      if (residual_standard_deviation != D_DEFAULT) {
        out_file << " " << (1. - *ppoint) / residual_standard_deviation
                 << " " << -*ppoint / residual_standard_deviation;
      }
      out_file << endl;
      ppoint++;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une courbe et des residus reduits correspondants
 *  au format Gnuplot.
 *
 *  arguments : path, pointeur sur les residus reduits.
 *
 *--------------------------------------------------------------*/

bool Curves::plot_print_0(const char *path , double *standard_residual) const

{
  bool status = false;
  register int i;
  int *pfrequency;
  double *ppoint , *pstandard_residual;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    // ecriture des reponses observees et des residus reduits

    pfrequency = frequency;
    ppoint = point[0];
    if (standard_residual) {
      pstandard_residual = standard_residual;
    }

    for (i = 0;i < length;i++) {
      if (*pfrequency > 0) {
        out_file << i << " " << *ppoint;
        if (standard_residual) {
          out_file << " " << *pstandard_residual;
        }
        out_file << " " << *pfrequency << endl;
      }
      pfrequency++;
      ppoint++;
      if (standard_residual) {
        pstandard_residual++;
      }
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWspace Function::binaryStoreSize() const

{
  RWspace size;


  size = Regression_kernel::binaryStoreSize();

  size += sizeof(true);
  if ((residual) && (frequency)) {
    size += sizeof(residual) * (max_value + 1) + sizeof(frequency) * (max_value + 1);
  }

  return size;
}


void Function::restoreGuts(RWvistream &is)

{
  bool status;
  register int i;


  remove();

  Regression_kernel::restoreGuts(is);

  is >> status;

  if (status) {
    residual = new double[max_value + 1];
    for (i = 0;i <= max_value;i++) {
      is >> residual[i];
    }

    frequency = new int[max_value + 1];
    for (i = 0;i <= max_value;i++) {
      is >> frequency[i];
    }
  }

  else {
    residual = 0;
    frequency = 0;
  }
}


void Function::restoreGuts(RWFile &file)

{
  bool status;


  remove();

  Regression_kernel::restoreGuts(file);

  file.Read(status);

  if (status) {
    residual = new double[max_value + 1];
    file.Read(residual , max_value + 1);

    frequency = new int[max_value + 1];
    file.Read(frequency , max_value + 1);
  }

  else {
    residual = 0;
    frequency = 0;
  }
}


void Function::saveGuts(RWvostream &os) const

{
  register int i;


  Regression_kernel::saveGuts(os);

  if ((residual) && (frequency)) {
    os << true;
    for (i = 0;i <= max_value;i++) {
      os << residual[i];
    }
    for (i = 0;i <= max_value;i++) {
      os << frequency[i];
    }
  }

  else {
    os << false;
  }
}


void Function::saveGuts(RWFile &file) const

{
  Regression_kernel::saveGuts(file);

  if ((residual) && (frequency)) {
    file.Write(true);
    file.Write(residual , max_value + 1);
    file.Write(frequency , max_value + 1);
  }

  else {
    file.Write(false);
  }
} */


/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Markov.
 *
 *--------------------------------------------------------------*/

Markov::Markov()

{
  nb_iterator = 0;
  markov_data = 0;

  order = 1;
  self_row = 0;

  homogeneous = true;

  homogeneity = 0;
  self_transition = 0;

  nb_output_process = 0;
  process = 0;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Markov (non-homogene).
 *
 *  arguments : nombre d'etats, identificateurs evolution
 *              des probabilites de transition.
 *
 *--------------------------------------------------------------*/

Markov::Markov(int inb_state , int *ident)
:Chain('o' , inb_state)

{
  register int i;


  nb_iterator = 0;
  markov_data = 0;

  order = 1;

  self_row = new int[nb_state];
  for (i = 0;i < nb_state;i++) {
    self_row[i] = i;
  }

  homogeneous = true;

  homogeneity = new bool[nb_state];
  self_transition = new Function*[nb_state];

  for (i = 0;i < nb_state;i++) {
    if (ident[i] == I_DEFAULT) {
      homogeneity[i] = true;
    }
    else {
      homogeneous = false;
      homogeneity[i] = false;
    }
    self_transition[i] = 0;
  }

  nb_output_process = 0;
  process = new Nonparametric_sequence_process*[1];
  process[0] = new Nonparametric_sequence_process(nb_state , nb_state); // , false);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Markov.
 *
 *  arguments : nombre d'etats, ordre, nombre de processus d'observation,
 *              nombre de valeurs observees par processus.
 *
 *--------------------------------------------------------------*/

Markov::Markov(int inb_state , int iorder , int inb_output_process , int *nb_value)
:Chain('o' , inb_state , (int)pow((double)inb_state , iorder) , true)

{
  register int i;


  nb_iterator = 0;
  markov_data = 0;

  order = iorder;

  self_row = new int[nb_state];
  self_row_computation();

  homogeneous = true;

  homogeneity = new bool[nb_state];
  self_transition = new Function*[nb_state];

  for (i = 0;i < nb_state;i++) {
    homogeneity[i] = true;
    self_transition[i] = 0;
  }

  nb_output_process = inb_output_process;
  process = new Nonparametric_sequence_process*[nb_output_process + 1];
  process[0] = new Nonparametric_sequence_process(nb_state , nb_state , false);
  for (i = 1;i <= nb_output_process;i++) {
    process[i] = new Nonparametric_sequence_process(nb_state , *nb_value++ , true);
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Markov (non-homogene).
 *
 *  arguments : pointeur sur un objet Chain, ordre, pointeurs sur des objets Function et
 *              sur un objet Nonparametric_process, longueur des sequences.
 *
 *--------------------------------------------------------------*/

Markov::Markov(const Chain *pchain , int iorder , const Function **pself_transition ,
               const Nonparametric_process *pobservation , int length)
:Chain(*pchain)

{
  register int i;


  nb_iterator = 0;
  markov_data = 0;

  order = iorder;

  self_row = new int[nb_state];
  self_row_computation();

  homogeneous = true;

  homogeneity = new bool[nb_state];
  self_transition = new Function*[nb_state];

  for (i = 0;i < nb_state;i++) {
    if (pself_transition[i]) {
      homogeneous = false;
      homogeneity[i] = false;
      self_transition[i] = new Function(*pself_transition[i]);
    }

    else {
      homogeneity[i] = true;
      self_transition[i] = 0;
    }
  }

  nb_output_process = (pobservation ? 1 : 0);
  process = new Nonparametric_sequence_process*[nb_output_process + 1];
  process[0] = new Nonparametric_sequence_process(nb_state , nb_state , false);
  if (pobservation) {
    process[1] = new Nonparametric_sequence_process(*pobservation);
  }

  characteristic_computation(length , true);
}


/*--------------------------------------------------------------*
 *
 *  Construction d'une chaine de Markov d'ordre superieur a 1 a partir
 *  d'une chaine de Markov.
 *
 *  arguments : ordre, reference sur un objet Markov.
 *
 *--------------------------------------------------------------*/

Markov::Markov(int iorder , const Markov &markov)
:Chain('o' , markov.nb_state , (int)pow((double)markov.nb_state , iorder) , true)

{
  register int i , j;


  accessibility = new bool*[nb_state];
  for (i = 0;i < nb_state;i++) {
    accessibility[i] = new bool[nb_state];
    for (j = 0;j < nb_state;j++) {
      accessibility[i][j] = markov.accessibility[i][j];
    }
  }

  nb_component = markov.nb_component;
  component_nb_state = new int[nb_component];
  component = new int*[nb_component];

  for (i = 0;i < nb_component;i++) {
    component_nb_state[i] = markov.component_nb_state[i];
    component[i] = new int[component_nb_state[i]];
    for (j = 0;j < component_nb_state[i];j++) {
      component[i][j] = markov.component[i][j];
    }
  }

  state_type = new char[nb_state];
  for (i = 0;i < nb_state;i++) {
    state_type[i] = markov.state_type[i];
  }

  nb_iterator = 0;
  markov_data = 0;

  order = iorder;

  self_row = new int[nb_state];
  self_row_computation();

  homogeneous = true;

  homogeneity = new bool[nb_state];
  self_transition = new Function*[nb_state];

  for (i = 0;i < nb_state;i++) {
    homogeneity[i] = true;
    self_transition[i] = 0;
  }

  nb_output_process = markov.nb_output_process;

  process = new Nonparametric_sequence_process*[nb_output_process + 1];
  for (i = 0;i <= nb_output_process;i++) {
    process[i] = new Nonparametric_sequence_process(*(markov.process[i]) , 'c' , false);
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Markov.
 *
 *  arguments : reference sur un objet Markov,
 *              flag copie de l'objet Markov_data,
 *              flag copie des lois caracteristiques.
 *
 *--------------------------------------------------------------*/

void Markov::copy(const Markov &markov , bool data_flag , bool characteristic_flag)

{
  register int i;


  nb_iterator = 0;

  if ((data_flag) && (markov.markov_data)) {
    markov_data = new Markov_data(*(markov.markov_data) , false);
  }
  else {
    markov_data = 0;
  }

  order = markov.order;

  self_row = new int[nb_state];
  for (i = 0;i < nb_state;i++) {
    self_row[i] = markov.self_row[i];
  }

  homogeneous = markov.homogeneous;

  homogeneity = new bool[nb_state];
  self_transition = new Function*[nb_state];

  for (i = 0;i < nb_state;i++) {
    homogeneity[i] = markov.homogeneity[i];
    if (homogeneity[i]) {
      self_transition[i] = 0;
    }
    else {
      self_transition[i] = new Function(*(markov.self_transition[i]));
    }
  }

  nb_output_process = markov.nb_output_process;

  process = new Nonparametric_sequence_process*[nb_output_process + 1];
  for (i = 0;i <= nb_output_process;i++) {
    process[i] = new Nonparametric_sequence_process(*(markov.process[i]) , 'c' , characteristic_flag);
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Markov.
 *
 *--------------------------------------------------------------*/

void Markov::remove()

{
  register int i;


  delete markov_data;

  delete [] self_row;

  if (self_transition) {
    for (i = 0;i < nb_state;i++) {
      if (!homogeneity[i]) {
        delete self_transition[i];
      }
    }
    delete [] self_transition;
  }

  delete [] homogeneity;

  if (process) {
    for (i = 0;i <= nb_output_process;i++) {
      delete process[i];
    }
    delete [] process;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Markov.
 *
 *--------------------------------------------------------------*/

Markov::~Markov()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Markov.
 *
 *  argument : reference sur un objet Markov.
 *
 *--------------------------------------------------------------*/

Markov& Markov::operator=(const Markov &markov)

{
  if (&markov != this) {
    remove();
    Chain::remove();

    Chain::copy(markov);
    copy(markov);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'une loi.
 *
 *  arguments : reference sur un objet Format_error, type de loi,
 *              variable, valeur (etat dans le cas des lois d'observation).
 *
 *--------------------------------------------------------------*/

Parametric_model* Markov::extract(Format_error &error , int type ,
                                  int variable , int value) const

{
  bool status = true;
  int hvariable;
  Distribution *pdist;
  Parametric *pparam;
  Parametric_model *dist;
  Histogram *phisto;


  dist = 0;
  error.init();

  pdist = 0;
  pparam = 0;

  if (type == OBSERVATION) {
    if ((variable < 1) || (variable > nb_output_process)) {
      status = false;
      error.update(SEQ_error[SEQR_OUTPUT_PROCESS_INDEX]);
    }

    else {
      if ((value < 0) || (value >= nb_state)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_STATE] << " " << value << " "
                      << SEQ_error[SEQR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      else {
        pdist = process[variable]->observation[value];
      }
    }
  }

  else {
    if ((variable < 0) || (variable > nb_output_process)) {
      status = false;
      error.update(SEQ_error[SEQR_OUTPUT_PROCESS_INDEX]);
    }

    else {
      if ((value < 0) || (value >= process[variable]->nb_value)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[variable == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                      << value << " " << SEQ_error[SEQR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      else {
        switch (type) {
        case FIRST_OCCURRENCE :
          pdist = process[variable]->first_occurrence[value];
          break;
        case RECURRENCE_TIME :
          if (process[variable]->recurrence_time) {
            pdist = process[variable]->recurrence_time[value];
          }
          break;
        case SOJOURN_TIME :
          if (process[variable]->sojourn_time) {
            pparam = process[variable]->sojourn_time[value];
          }
          break;
        case NB_RUN :
          pdist = process[variable]->nb_run[value];
          break;
        case NB_OCCURRENCE :
          pdist = process[variable]->nb_occurrence[value];
          break;
        }

        if ((!pdist) && (!pparam)) {
          status = false;
          error.update(SEQ_error[SEQR_NON_EXISTING_CHARACTERISTIC_DISTRIBUTION]);
        }
      }
    }
  }

  if (status) {
    phisto = 0;

    if (markov_data) {
      switch (markov_data->type[0]) {
      case INT_VALUE :
        hvariable = variable - 1;
        break;
      case STATE :
        hvariable = variable;
        break;
      }

      if (hvariable >= 0) {
        switch (type) {

        case OBSERVATION : {
          if ((markov_data->observation) && (markov_data->observation[hvariable])) {
            phisto = markov_data->observation[hvariable][value];
          }
          break;
        }

        case FIRST_OCCURRENCE : {
          phisto = markov_data->characteristics[hvariable]->first_occurrence[value];
          break;
        }

        case RECURRENCE_TIME : {
          if (markov_data->characteristics[hvariable]->recurrence_time[value]->nb_element > 0) {
            phisto = markov_data->characteristics[hvariable]->recurrence_time[value];
          }
          break;
        }

        case SOJOURN_TIME : {
          if (markov_data->characteristics[hvariable]->sojourn_time[value]->nb_element > 0) {
            phisto = markov_data->characteristics[hvariable]->sojourn_time[value];
          }
          break;
        }

        case NB_RUN : {
          phisto = markov_data->characteristics[hvariable]->nb_run[value];
          break;
        }

        case NB_OCCURRENCE : {
          phisto = markov_data->characteristics[hvariable]->nb_occurrence[value];
          break;
        }
        }
      }
    }

    if (pdist) {
      dist = new Parametric_model(*pdist , phisto);
    }
    else if (pparam) {
      dist = new Parametric_model(*pparam , phisto);
    }
  }

  return dist;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de la partie "donnees" d'un objet Markov.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Markov_data* Markov::extract_data(Format_error &error) const

{
  bool status = true;
  Markov_data *seq;


  seq = 0;
  error.init();

  if (!markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }
  else if (nb_output_process + 1 != markov_data->nb_variable) {
    status = false;
    error.update(SEQ_error[SEQR_STATE_SEQUENCES]);
  }

  if (status) {
    seq = new Markov_data(*markov_data);
    seq->markov = new Markov(*this , false);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Application d'un seuil sur les parametres d'une chaine de Markov.
 *
 *  argument : probabilite minimum.
 *
 *--------------------------------------------------------------*/

Markov* Markov::thresholding(double min_probability) const

{
  register int i;
  Markov *markov;


  markov = new Markov(*this , false , false);

  markov->Chain::thresholding(min_probability);
  markov->component_computation();

  for (i = 1;i <= markov->nb_output_process;i++) {
    markov->process[i]->thresholding(min_probability);
  }

  return markov;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un "Mixture Transition Distribution" modele a partir
 *  d'une chaine de Markov d'ordre 1 et de poids.
 *
 *  arguments : reference sur un objet Format_error, stream, ordre, poids,
 *              longueur des sequences.
 *
 *--------------------------------------------------------------*/

Markov* Markov::MTD_model_building(Format_error &error , ostream &os , int iorder ,
                                   double *weight , int length) const

{
  bool status = true;
  register int i , j , k;
  int state_index[ORDER];
  double negative_weight_sum , positive_weight_sum , min_transition , max_transition ,
         *pinitial , *cinitial , *ptransition , *pweight;
  Markov *markov;


  markov = 0;
  error.init();

  if (order != 1) {
    status = false;
    error.update(SEQ_error[SEQR_INPUT_MARKOV_ORDER]);
  }
  if (!homogeneous) {
    status = false;
    error.update(SEQ_error[SEQR_NON_HOMOGENEOUS_INPUT_MARKOV]);
  }

  if ((iorder < 2) || (iorder > ORDER)) {
    status = false;
    error.update(SEQ_error[SEQR_ORDER]);
  }

  else {
    negative_weight_sum = 0.;
    positive_weight_sum = 0.;
    for (i = 0;i < iorder;i++) {
      if (weight[i] < 0.) {
        negative_weight_sum += weight[i];
      }
      else {
        positive_weight_sum += weight[i];
      }
    }

    if ((negative_weight_sum + positive_weight_sum < 1. - DOUBLE_ERROR) ||
        (negative_weight_sum + positive_weight_sum > 1. + DOUBLE_ERROR)) {
      status = false;
      error.update(SEQ_error[SEQR_WEIGHT_VALUES]);
    }

    // verification des conditions de Raftery et Tavare (cas de poids negatifs)

    else if (negative_weight_sum < 0.) {

#     ifdef MESSAGE
      os << "\nConditions of Raftery and Tavare on the original transition probabilities and the weights" << endl;
#     endif

      for (i = 0;i < nb_state;i++) {
        min_transition = 1.;
        max_transition = 0.;
        for (j = 0;j < nb_state;j++) {
          if (transition[j][i] < min_transition) {
            min_transition = transition[j][i];
          }
          if (transition[j][i] > max_transition) {
            max_transition = transition[j][i];
          }
        }

#       ifdef MESSAGE
        os << STAT_label[STATL_STATE] << " " << i << ": "
           << positive_weight_sum * min_transition + negative_weight_sum * max_transition << endl;
#       endif

        if (positive_weight_sum * min_transition + negative_weight_sum * max_transition < 0.) {
          status = false;
          error.update(SEQ_error[SEQR_WEIGHTS_INCOMPATIBLE_WITH_TRANSITIONS]);
        }
      }
    }
  }

  if (status) {
    markov = new Markov(iorder , *this);

    pinitial = markov->initial;
    cinitial = initial;
    for (i = 0;i < markov->nb_state;i++) {
      *pinitial++ = *cinitial++;
    }

    // calcul des probabilites de transition au sens du modele MTD

    for (i = 0;i < markov->order;i++) {
      state_index[i] = 0;
    }

    for (i = 0;i < markov->nb_row;i++) {
      ptransition = markov->transition[i];
      for (j = 0;j < markov->nb_state;j++) {
        *ptransition = 0.;
        pweight = weight;
        for (k = 0;k < markov->order;k++) {
          *ptransition += *pweight++ * transition[state_index[markov->order - 1 - k]][j];
        }
        ptransition++;
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

    markov->characteristic_computation(length , true);
  }

  return markov;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Markov a partir d'un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              longueur des sequences.
 *
 *--------------------------------------------------------------*/

Markov* markov_ascii_read(Format_error &error , const char *path , int length)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status , lstatus , homogeneous;
  register int i , j;
  int line , order , homogeneity , nb_state;
  double index;
  const Chain *chain;
  const Function **self_transition;
  const Nonparametric_process *observation;
  Markov *markov;
  ifstream in_file(path);


  markov = 0;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;

    if (length < 2) {
      status = false;
      error.update(SEQ_error[SEQR_SHORT_SEQUENCE_LENGTH]);
    }
    if (length > MAX_LENGTH) {
      status = false;
      error.update(SEQ_error[SEQR_LONG_SEQUENCE_LENGTH]);
    }

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

        // test mot cle MARKOV_CHAIN / NON-HOMOGENEOUS_MARKOV_CHAIN

        if (i == 0) {
          if (token == SEQ_word[SEQW_MARKOV_CHAIN]) {
            homogeneous = true;
          }
          else {
            if (token == SEQ_word[SEQW_NON_HOMOGENEOUS_MARKOV_CHAIN]) {
              homogeneous = false;
              order = 1;
            }
            else {
              status = false;
              error.update(STAT_parsing[STATP_KEY_WORD] , line);
            }
          }
        }

        i++;
      }

      if (i > 0) {
        if (i != 1) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }
        break;
      }
    }

    // analyse du format et lecture de la chaine de Markov

    chain = chain_parsing(error , in_file , line , 'o' , homogeneous , order);

    if (chain) {
      nb_state = chain->nb_state;
      self_transition = new const Function*[nb_state];
      for (i = 0;i < nb_state;i++) {
        self_transition[i] = 0;
      }

      if (!homogeneous) {

        // analyse format des fonctions d'evolution
        // des probabilites de rester dans un etat

        for (i = 0;i < nb_state;i++) {
          homogeneity = I_DEFAULT;

          while (buffer.readLine(in_file , false)) {
            line++;

#           ifdef DEBUG
            cout << line << "  " << buffer << endl;
#           endif

            position = buffer.first('#');
            if (position != RW_NPOS) {
              buffer.remove(position);
            }
            j = 0;

            RWCTokenizer next(buffer);

            while (!((token = next()).isNull())) {
              switch (j) {

              // test mot cle STATE

              case 0 : {
                if (token != STAT_word[STATW_STATE]) {
                  status = false;
                  error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_STATE] , line , j + 1);
                }
                break;
              }

              // test indice etat

              case 1 : {
                lstatus = locale.stringToNum(token , &index);
                if ((lstatus) && (index != i)) {
                  lstatus = false;
                }

                if (!lstatus) {
                  status = false;
                  error.correction_update(STAT_parsing[STATP_STATE_INDEX] , i , line , j + 1);
                }
                break;
              }

              // test mot cle HOMOGENEOUS / NON-HOMOGENEOUS

              case 2 : {
                if (token == SEQ_word[SEQW_HOMOGENEOUS]) {
                  homogeneity = true;
                }
                else {
                  if (token == SEQ_word[SEQW_NON_HOMOGENEOUS]) {
                    homogeneity = false;
                  }
                  else {
                    status = false;
                    error.update(STAT_parsing[STATP_KEY_WORD] , line , j + 1);
                  }
                }
                break;
              }
              }

              j++;
            }

            if (j > 0) {
              if (j != 3) {
                status = false;
                error.update(STAT_parsing[STATP_FORMAT] , line);
              }

              if (!homogeneity) {
                self_transition[i] = function_parsing(error , in_file , line , length);
                if (!self_transition[i]) {
                  status = false;
                }
              }

              break;
            }
          }

          if (homogeneity == I_DEFAULT) {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line);
          }
        }
      }

      // analyse du format et lecture des lois d'observation

      observation = 0;

      if (homogeneous) {
        while (buffer.readLine(in_file , false)) {
          line++;

#         ifdef DEBUG
          cout << line << "  " << buffer << endl;
#         endif

          position = buffer.first('#');
          if (position != RW_NPOS) {
            buffer.remove(position);
          }
          i = 0;

          RWCTokenizer next(buffer);

          while (!((token = next()).isNull())) {

            // test mot cle OUTPUT_PROCESS

            if (i == 0) {
              if (token != STAT_word[STATW_OUTPUT_PROCESS]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_OUTPUT_PROCESS] , line);
              }
            }

            i++;
          }

          if (i > 0) {
            if (i != 1) {
              status = false;
              error.update(STAT_parsing[STATP_FORMAT] , line);
            }

            // analyse du format et lecture des lois d'observation

            observation = observation_parsing(error , in_file , line , nb_state , false);
            if (!observation) {
              status = false;
            }

            break;
          }
        }
      }

      while (buffer.readLine(in_file , false)) {
        line++;

#       ifdef DEBUG
        cout << line << "  " << buffer << endl;
#       endif

        position = buffer.first('#');
        if (position != RW_NPOS) {
          buffer.remove(position);
        }
        if (!(buffer.isNull())) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }
      }

      if (status) {
        markov = new Markov(chain , order , self_transition , observation , length);
      }

      delete chain;

      for (i = 0;i < nb_state;i++) {
        delete self_transition[i];
      }
      delete [] self_transition;

      delete observation;
    }
  }

  return markov;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet Markov.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Markov::line_write(ostream &os) const

{
  if (!homogeneous) {
    os << SEQ_label[SEQL_NON_HOMOGENEOUS] << " - ";
  }

  os << nb_state << " " << STAT_word[STATW_STATES];

  if (homogeneous) {
    os << "   " << STAT_word[STATW_ORDER] << " " << order;
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Markov et de la structure de donnees associee.
 *
 *  arguments : stream, pointeur sur un objet Markov_data,
 *              flag niveau de detail, flag fichier, flag Markov cache.
 *
 *--------------------------------------------------------------*/

ostream& Markov::ascii_write(ostream &os , const Markov_data *seq , bool exhaustive ,
                             bool file_flag , bool hidden) const

{
  register int i;
  int variable;
  Histogram **observation = 0;
  Sequence_characteristics *characteristics = 0;


  switch (hidden) {

  case false : {
    switch (homogeneous) {
    case true :
      os << SEQ_word[SEQW_MARKOV_CHAIN] << endl;
      break;
    case false :
      os << SEQ_word[SEQW_NON_HOMOGENEOUS_MARKOV_CHAIN] << endl;
      break;
    }
    break;
  }

  case true : {
    os << SEQ_word[SEQW_HIDDEN_MARKOV_CHAIN] << endl;
    break;
  }
  }

  // ecriture des parametres de la chaine de Markov

  ascii_print(os , file_flag , homogeneous , order);

  // ecriture des parametres des fonctions d'evolution
  // des probabilites de rester dans un etat

  if (!homogeneous) {
    for (i = 0;i < nb_state;i++) {
      os << "\n" << STAT_word[STATW_STATE] << " " << i << " ";

      switch (homogeneity[i]) {
      case true :
        os << SEQ_word[SEQW_HOMOGENEOUS] << endl;
        break;
      case false :
        os << SEQ_word[SEQW_NON_HOMOGENEOUS] << endl;
        self_transition[i]->ascii_print(os , exhaustive , file_flag ,
                                        (seq ? seq->self_transition[i] : 0));
        break;
      }
    }
  }

  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
  }
  process[0]->ascii_print(os , 0 , 0 , characteristics , exhaustive , file_flag);

  // ecriture des lois associees a chaque processus d'observation

  if (hidden) {
    os << "\n" << nb_output_process << " "
       << STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES] << endl;
  }

  for (i = 1;i <= nb_output_process;i++) {
    os << "\n" << STAT_word[STATW_OUTPUT_PROCESS];
    if (hidden) {
      os << " " << i;
    }
    os << endl;

    if (seq) {
      switch (seq->type[0]) {
      case INT_VALUE :
        variable = i - 1;
        break;
      case STATE :
        variable = i;
        break;
      }

      if (seq->observation) {
        observation = seq->observation[variable];
      }
      characteristics = seq->characteristics[variable];
    }

    process[i]->ascii_print(os , i , observation , characteristics ,
                            exhaustive , file_flag);
  }

  if (seq) {
    int nb_parameter = nb_parameter_computation(hidden ? MIN_PROBABILITY : 0.);
    double information , likelihood;


    // ecriture de l'histogramme des longueurs des sequences

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
    seq->hlength->ascii_characteristic_print(os , false , file_flag);

    if (exhaustive) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << "   | " << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << endl;
      seq->hlength->ascii_print(os , file_flag);
    }

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_CUMUL_LENGTH] << ": " << seq->cumul_length << endl;

    // ecriture de la quantite d'information des sequences dans le cas iid

    information = seq->iid_information_computation();

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_IID_INFORMATION] << ": " << information << " ("
       << information / seq->cumul_length << ")" << endl;

    // ecriture des vraisemblances des sequences

    switch (hidden) {

    case false : {
      likelihood = seq->likelihood;

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   ("
         << STAT_label[STATL_NORMALIZED] << ": " << likelihood / seq->cumul_length << ")" << endl;
      break;
    }

    case true : {
      if (seq->likelihood != D_INF) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << seq->likelihood << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << seq->likelihood / seq->cumul_length << ")" << endl;
      }

      likelihood = seq->hidden_likelihood;

      if (likelihood != D_INF) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << likelihood << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << likelihood / seq->cumul_length << ")" << endl;
      }
      break;
    }
    }

    if ((likelihood != D_INF) && (nb_component == 1)) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << "): "
         << 2 * (likelihood - nb_parameter) << endl;

      if (nb_parameter < seq->cumul_length - 1) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
           << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << "): "
           << 2 * (likelihood - (double)(nb_parameter * seq->cumul_length) /
             (double)(seq->cumul_length - nb_parameter - 1)) << endl;
      }

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << "): "
         << 2 * likelihood - nb_parameter * log((double)seq->cumul_length) << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Markov.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Markov::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , markov_data , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Markov dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Markov::ascii_write(Format_error &error , const char *path ,
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
    ascii_write(out_file , markov_data , exhaustive , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Markov et de la structure de donnees associee
 *  dans un fichier au format tableur.
 *
 *  arguments : stream, pointeur sur un objet Markov_data, flag Markov cache.
 *
 *--------------------------------------------------------------*/

ostream& Markov::spreadsheet_write(ostream &os , const Markov_data *seq ,
                                   bool hidden) const

{
  register int i;
  int variable;
  Histogram **observation = 0;
  Sequence_characteristics *characteristics = 0;


  switch (hidden) {

  case false : {
    switch (homogeneous) {
    case true :
      os << SEQ_word[SEQW_MARKOV_CHAIN] << endl;
      break;
    case false :
      os << SEQ_word[SEQW_NON_HOMOGENEOUS_MARKOV_CHAIN] << endl;
      break;
    }
    break;
  }

  case true : {
    os << SEQ_word[SEQW_HIDDEN_MARKOV_CHAIN] << endl;
    break;
  }
  }

  // ecriture des parametres de la chaine de Markov

  spreadsheet_print(os , homogeneous , order);

  // ecriture des parametres des fonctions d'evolution
  // des probabilites de rester dans un etat

  if (!homogeneous) {
    for (i = 0;i < nb_state;i++) {
      os << "\n" << STAT_word[STATW_STATE] << "\t" << i << "\t";

      switch (homogeneity[i]) {
      case true :
        os << SEQ_word[SEQW_HOMOGENEOUS] << endl;
        break;
      case false :
        os << SEQ_word[SEQW_NON_HOMOGENEOUS] << endl;
        self_transition[i]->spreadsheet_print(os , (seq ? seq->self_transition[i] : 0));
        break;
      }
    }
  }

  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
  }
  process[0]->spreadsheet_print(os , 0 , 0 , characteristics);

  // ecriture des lois associees a chaque processus d'observation

  if (hidden) {
    os << "\n" << nb_output_process << "\t"
       << STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES] << endl;
  }

  for (i = 1;i <= nb_output_process;i++) {
    os << "\n" << STAT_word[STATW_OUTPUT_PROCESS];
    if (hidden) {
      os << "\t" << i;
    }
    os << endl;

    if (seq) {
      switch (seq->type[0]) {
      case INT_VALUE :
        variable = i - 1;
        break;
      case STATE :
        variable = i;
        break;
      }

      if (seq->observation) {
        observation = seq->observation[variable];
      }
      characteristics = seq->characteristics[variable];
    }

    process[i]->spreadsheet_print(os , i , observation , characteristics);
  }

  if (seq) {
    int nb_parameter = nb_parameter_computation(hidden ? MIN_PROBABILITY : 0.);
    double information , likelihood;


    // ecriture de l'histogramme des longueurs des sequences

    os << "\n" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
    seq->hlength->spreadsheet_characteristic_print(os);

    os << "\n\t" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << endl;
    seq->hlength->spreadsheet_print(os);

    os << "\n" << SEQ_label[SEQL_CUMUL_LENGTH] << "\t" << seq->cumul_length << endl;

    // ecriture de la quantite d'information des sequences dans le cas iid

    information = seq->iid_information_computation();

    os << "\n" << SEQ_label[SEQL_IID_INFORMATION] << "\t" << information << "\t"
       << information / seq->cumul_length << endl;

    // ecriture des vraisemblances des sequences

    switch (hidden) {

    case false : {
      likelihood = seq->likelihood;

      os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << likelihood << "\t"
         << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / seq->cumul_length << endl;
      break;
    }

    case true : {
      if (seq->likelihood != D_INF) {
        os << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << "\t" << seq->likelihood << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << seq->likelihood / seq->cumul_length << endl;
      }

      likelihood = seq->hidden_likelihood;

      if (likelihood != D_INF) {
        os << "\n" << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << "\t" << likelihood << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / seq->cumul_length << endl;
      }
      break;
    }
    }

    if ((likelihood != D_INF) && (nb_component == 1)) {
      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << ")\t"
         << 2 * (likelihood - nb_parameter) << endl;

      if (nb_parameter < seq->cumul_length - 1) {
        os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
           << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << ")\t"
           << 2 * (likelihood - (double)(nb_parameter * seq->cumul_length) /
              (double)(seq->cumul_length - nb_parameter - 1)) << endl;
      }

      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << ")\t"
         << 2 * likelihood - nb_parameter * log((double)seq->cumul_length) << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Markov dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Markov::spreadsheet_write(Format_error &error , const char *path) const

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
    spreadsheet_write(out_file , markov_data);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Markov et de la structure de donnees associee.
 *
 *  arguments : prefixe des fichiers, titre des figures,
 *              pointeur sur les sequences observees.
 *
 *--------------------------------------------------------------*/

bool Markov::plot_write(const char *prefix , const char *title ,
                        const Markov_data *seq) const

{
  bool status;
  register int i , j;
  int variable , start , *pfrequency , max_frequency[NB_STATE];
  double residual_mean , residual_standard_deviation , *standard_residual , *presidual ,
         *pstandard_residual , min_standard_residual[NB_STATE] , max_standard_residual[NB_STATE];
  Histogram *hlength = 0 , **observation = 0;
  Sequence_characteristics *characteristics = 0;
  ostringstream data_file_name[NB_STATE * 2];


  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
    hlength = seq->hlength;
  }

  status = process[0]->plot_print(prefix , title , 0 , 0 , characteristics , hlength);

  if (status) {
    if (!homogeneous) {

      // ecriture des fichiers de donnees

      for (i = 0;i < nb_state;i++) {
        if (!homogeneity[i]) {
          if (seq) {
            max_frequency[i] = seq->self_transition[i]->max_frequency_computation();

            // calcul des residus reduits

            residual_mean = self_transition[i]->residual_mean_computation();
            residual_standard_deviation = sqrt(self_transition[i]->residual_variance_computation(residual_mean));

            standard_residual = new double[self_transition[i]->max_value + 1];

            pfrequency = self_transition[i]->frequency;
            pstandard_residual = standard_residual;
            presidual = self_transition[i]->residual;
            min_standard_residual[i] = 0.;
            max_standard_residual[i] = 0.;

            for (j = 0;j <= self_transition[i]->max_value;j++) {
              if (*pfrequency++ > 0) {
                *pstandard_residual = *presidual / residual_standard_deviation;
                if (*pstandard_residual < min_standard_residual[i]) {
                  min_standard_residual[i] = *pstandard_residual;
                }
                if (*pstandard_residual > max_standard_residual[i]) {
                  max_standard_residual[i] = *pstandard_residual;
                }
                presidual++;
                pstandard_residual++;
              }
            }
          }

          data_file_name[i * 2] << prefix << i * 2 << ".dat";
          self_transition[i]->plot_print((data_file_name[i * 2].str()).c_str() ,
                                         (seq ? residual_standard_deviation : D_DEFAULT));

          if (seq) {
            data_file_name[i * 2 + 1] << prefix << i * 2 + 1 << ".dat";
            seq->self_transition[i]->plot_print_0((data_file_name[i * 2 + 1].str()).c_str() ,
                                                  standard_residual);
            delete [] standard_residual;
          }
        }
      }

      // ecriture du fichier de commandes et du fichier d'impression

      for (i = 0;i < 2;i++) {
        ostringstream file_name[2];

        switch (i) {
        case 0 :
          file_name[0] << prefix << 0 << 0 << ".plot";
          break;
        case 1 :
          file_name[0] << prefix << 0 << 0 << ".print";
          break;
        }

        ofstream out_file((file_name[0].str()).c_str());

        if (i == 1) {
          out_file << "set terminal postscript" << endl;
          file_name[1] << label(prefix) << 0 << 0 << ".ps";
          out_file << "set output \"" << file_name[1].str() << "\"\n\n";
        }

        out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                 << "set title";
        if (title) {
          out_file << " \"" << title << "\"";
        }
        out_file << "\n\n";

        start = true;
        for (j = 0;j < nb_state;j++) {
          if (!homogeneity[j]) {
            if (!start) {
              if (i == 0) {
                out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
              }
              out_file << endl;
            }
            else {
              start = false;
            }

            if (self_transition[j]->max_value < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }

            out_file << "plot [0:" << self_transition[j]->max_value << "] [0:1] ";
            if (seq) {
              out_file << "\""<< label((data_file_name[j * 2 + 1].str()).c_str())
                       << "\" using 1:2 title \"" << STAT_label[STATL_STATE] << " " << j << " - "
                       << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION]
                       << "\" with points,\\" << endl;
            }
            out_file << "\"" << label((data_file_name[j * 2].str()).c_str())
                     << "\" using 1 title \"" << STAT_label[STATL_STATE] << " " << j << " - "
                     << SEQ_label[SEQL_THEORETICAL] << " " << SEQ_label[SEQL_SELF_TRANSITION]
                     << "\" with linespoints" << endl;

            if (self_transition[j]->max_value < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }

            if (seq) {
              if (i == 0) {
                out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
              }
              out_file << endl;

              out_file << "set xlabel \"" << SEQ_label[SEQL_INDEX] << "\"" << endl;
              out_file << "set ylabel \"" << STAT_label[STATL_STANDARDIZED_RESIDUAL] << "\"" << endl;

              if (seq->self_transition[j]->length - 1 < TIC_THRESHOLD) {
                out_file << "set xtics 0,1" << endl;
              }

              out_file << "plot [0:" << seq->self_transition[j]->length - 1 << "] ["
                       << min_standard_residual[j] << ":" << max_standard_residual[j] << "] \""
                       << label((data_file_name[j * 2 + 1].str()).c_str())
                       << "\" using 1:3 notitle with points";
              if (((1. - self_transition[j]->point[0]) / residual_standard_deviation <= max_standard_residual[j]) ||
                  ((1. - self_transition[j]->point[seq->self_transition[j]->length - 1]) / residual_standard_deviation <=
                   max_standard_residual[j])) {
                out_file << ",\\\n\"" << label((data_file_name[j * 2].str()).c_str())
                         << "\" using 2 title \"" << SEQ_label[SEQL_ASYMPTOTE] << "\" with lines";
              }
              if ((-self_transition[j]->point[0] / residual_standard_deviation >= min_standard_residual[j]) ||
                  (-self_transition[j]->point[seq->self_transition[j]->length - 1] / residual_standard_deviation >=
                   min_standard_residual[j])) {
                out_file << ",\\\n\"" << label((data_file_name[j * 2].str()).c_str())
                         << "\" using 3 title \"" << SEQ_label[SEQL_ASYMPTOTE] << "\" with lines";
              }
              out_file << endl;

              out_file << "set xlabel" << endl;
              out_file << "set ylabel" << endl;

              if (seq->self_transition[j]->length - 1 < TIC_THRESHOLD) {
                out_file << "set xtics autofreq" << endl;
              }

              if (i == 0) {
                out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
              }
              out_file << endl;

              out_file << "set xlabel \"" << SEQ_label[SEQL_INDEX] << "\"" << endl;
              out_file << "set ylabel \"" << STAT_label[STATL_FREQUENCY] << "\"" << endl;

              if (seq->self_transition[j]->length - 1 < TIC_THRESHOLD) {
                out_file << "set xtics 0,1" << endl;
              }
              if (max_frequency[j] < TIC_THRESHOLD) {
                out_file << "set ytics 0,1" << endl;
              }

              out_file << "plot [0:" << seq->self_transition[j]->length - 1
                       << "] [0:" << max_frequency[j] << "] \""
                       << label((data_file_name[j * 2 + 1].str()).c_str())
                       << "\" using 1:4 notitle with impulses" << endl;

              out_file << "set xlabel" << endl;
              out_file << "set ylabel" << endl;

              if (seq->self_transition[j]->length - 1 < TIC_THRESHOLD) {
                out_file << "set xtics autofreq" << endl;
              }
              if (max_frequency[j] < TIC_THRESHOLD) {
                out_file << "set ytics autofreq" << endl;
              }
            }
          }
        }

        if (i == 1) {
          out_file << "\nset terminal x11" << endl;
        }

        out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
    }

    if (seq) {
      hlength = seq->hlength;
    }

    for (i = 1;i <= nb_output_process;i++) {
      if (seq) {
        switch (seq->type[0]) {
        case INT_VALUE :
          variable = i - 1;
          break;
        case STATE :
          variable = i;
          break;
        }

        if (seq->observation) {
          observation = seq->observation[variable];
        }
        characteristics = seq->characteristics[variable];
      }

      process[i]->plot_print(prefix , title , i , observation ,
                             characteristics , hlength);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Markov.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Markov::plot_write(Format_error &error , const char *prefix ,
                        const char *title) const

{
  bool status = plot_write(prefix , title , markov_data);

  error.init();

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWDEFINE_COLLECTABLE(Markov , STATI_MARKOV);


RWspace Markov::binaryStoreSize() const

{
  register int i;
  RWspace size;


  size = Chain::binaryStoreSize() + sizeof(order) + sizeof(*self_row) * nb_state +
         sizeof(homogeneous) + sizeof(*homogeneity) * nb_state;

  for (i = 0;i < nb_state;i++) {
    if (!homogeneity[i]) {
      size += self_transition[i]->binaryStoreSize();
    }
  }

  size += sizeof(nb_output_process);
  for (i = 0;i <= nb_output_process;i++) {
    size += process[i]->binaryStoreSize();
  }

  if (markov_data) {
    size += markov_data->recursiveStoreSize();
  }

  return size;
}


void Markov::restoreGuts(RWvistream &is)

{
  register int i;


  remove();

  Chain::restoreGuts(is);

  is >> order;

  self_row = new int[nb_state];
  for (i = 0;i < nb_state;i++) {
    is >> self_row[i];
  }

  is >> homogeneous;

  homogeneity = new bool[nb_state];
  for (i = 0;i < nb_state;i++) {
    is >> homogeneity[i];
  }

  self_transition = new Function*[nb_state];
  for (i = 0;i < nb_state;i++) {
    if (homogeneity[i]) {
      self_transition[i] = 0;
    }
    else {
      self_transition[i] = new Function();
      self_transition[i]->restoreGuts(is);
    }
  }

  is >> nb_output_process;

  process = new Nonparametric_sequence_process*[nb_output_process + 1];
  for (i = 0;i <= nb_output_process;i++) {
    process[i] = new Nonparametric_sequence_process();
    process[i]->restoreGuts(is);
  }

  is >> markov_data;
  if (markov_data == RWnilCollectable) {
    markov_data = 0;
  }
}


void Markov::restoreGuts(RWFile &file)

{
  register int i;


  remove();

  Chain::restoreGuts(file);

  file.Read(order);

  self_row = new int[nb_state];
  file.Read(self_row , nb_state);

  file.Read(homogeneous);

  homogeneity = new bool[nb_state];
  file.Read(homogeneity , nb_state);

  self_transition = new Function*[nb_state];
  for (i = 0;i < nb_state;i++) {
    if (homogeneity[i]) {
      self_transition[i] = 0;
    }
    else {
      self_transition[i] = new Function();
      self_transition[i]->restoreGuts(file);
    }
  }

  file.Read(nb_output_process);

  process = new Nonparametric_sequence_process*[nb_output_process + 1];
  for (i = 0;i <= nb_output_process;i++) {
    process[i] = new Nonparametric_sequence_process();
    process[i]->restoreGuts(file);
  }

  file >> markov_data;
  if (markov_data == RWnilCollectable) {
    markov_data = 0;
  }
}


void Markov::saveGuts(RWvostream &os) const

{
  register int i;


  Chain::saveGuts(os);

  os << order;
  for (i = 0;i < nb_state;i++) {
    os << self_row[i];
  }

  os << homogeneous;
  for (i = 0;i < nb_state;i++) {
    os << homogeneity[i];
  }

  for (i = 0;i < nb_state;i++) {
    if (!homogeneity[i]) {
      self_transition[i]->saveGuts(os);
    }
  }

  os << nb_output_process;
  for (i = 0;i <= nb_output_process;i++) {
    process[i]->saveGuts(os);
  }

  if (markov_data) {
    os << markov_data;
  }
  else {
    os << RWnilCollectable;
  }
}


void Markov::saveGuts(RWFile &file) const

{
  register int i;


  Chain::saveGuts(file);

  file.Write(order);
  file.Write(self_row , nb_state);

  file.Write(homogeneous);
  file.Write(homogeneity , nb_state);

  for (i = 0;i < nb_state;i++) {
    if (!homogeneity[i]) {
      self_transition[i]->saveGuts(file);
    }
  }

  file.Write(nb_output_process);
  for (i = 0;i <= nb_output_process;i++) {
    process[i]->saveGuts(file);
  }

  if (markov_data) {
    file << markov_data;
  }
  else {
    file << RWnilCollectable;
  }
} */


/*--------------------------------------------------------------*
 *
 *  Calcul des indices des lignes de la matrice des probabilites
 *  de transition correspondant aux probabilites de rester dans un etat.
 *
 *--------------------------------------------------------------*/

void Markov::self_row_computation()

{
  register int i , j;
  int state_index[ORDER] , *pself_row;


  for (i = 0;i < order;i++) {
    state_index[i] = 0;
  }

  pself_row = self_row;
  for (i = 0;i < nb_row;i++) {
    for (j = 1;j < order;j++) {
      if (state_index[j] != state_index[j - 1]) {
        break;
      }
    }

    if (j == order) {
      *pself_row++ = i;
    }

    for (j = 0;j < order;j++) {
      if (state_index[j] < nb_state - 1) {
        state_index[j]++;
        break;
      }
      else {
        state_index[j] = 0;
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Extraction des classes d'une chaine de Markov a partir
 *  de l'accessibilite des etats.
 *
 *--------------------------------------------------------------*/

void Markov::component_computation()

{
  bool **logic_transition;
  register int i , j , k;
  int power = nb_row / nb_state;
  double sum , **ptransition;


  // calcul de la matrice des transitions possibles entre etats

  logic_transition = new bool*[nb_state];

  for (i = 0;i < nb_state;i++) {
    logic_transition[i] = new bool[nb_state];

    for (j = 0;j < nb_state;j++) {
      if (j == i) {
        logic_transition[i][j] = false;
      }

      else {
        ptransition = transition + i * power;
        sum = 0.;
        for (k = 0;k < power;k++) {
          sum += *(*ptransition + j);
          ptransition++;
        }
        logic_transition[i][j] = (sum == 0. ? false : true);
      }
    }
  }

  Chain::component_computation(logic_transition);

  for (i = 0;i < nb_state;i++) {
    delete [] logic_transition[i];
  }
  delete [] logic_transition;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de parametres independants d'un objet Markov.
 *
 *  argument : probabilite minimum.
 *
 *--------------------------------------------------------------*/

int Markov::nb_parameter_computation(double min_probability) const

{
  register int i;
  int nb_parameter = Chain::nb_parameter_computation(min_probability);


  for (i = 1;i <= nb_output_process;i++) {
    nb_parameter += process[i]->nb_parameter_computation(min_probability);
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Markov_data.
 *
 *--------------------------------------------------------------*/

Markov_data::Markov_data()

{
  markov = 0;
  chain_data = 0;
  likelihood = D_INF;
  hidden_likelihood = D_INF;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Markov_data.
 *
 *  arguments : nombre de variables, histogramme des longueurs des sequences,
 *              flag initialisation.
 *
 *--------------------------------------------------------------*/

Markov_data::Markov_data(int inb_variable , const Histogram &ihlength , bool init_flag)
:Markovian_sequences(inb_variable , ihlength , init_flag)

{
  markov = 0;
  chain_data = 0;
  likelihood = D_INF;
  hidden_likelihood = D_INF;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Markov_data.
 *
 *  arguments : nombre de sequences, longueurs des sequences, sequences,
 *              pointeur sur un objet Markov.
 *
 *--------------------------------------------------------------*/

Markov_data::Markov_data(int inb_sequence , int *ilength , int ***isequence ,
                         const Markov &imarkov)
:Markovian_sequences(imarkov.nb_output_process + 1 , inb_sequence , ilength , isequence)

{
  register int i;


  type[0] = STATE;
  markov = new Markov(imarkov , false);

  for (i = 0;i < markov->nb_state;i++) {
    if (!(markov->homogeneity[i])) {
      if (markov->self_transition[i]->max_value < max_length - 1) {
        delete [] markov->self_transition[i]->point;
        markov->self_transition[i]->max_value = max_length - 1;
        markov->self_transition[i]->point = new double[markov->self_transition[i]->max_value + 1];
        markov->self_transition[i]->computation();
      }
    }
  }

  // extraction des caracteristiques des sequences simulees

  if (!(markov->homogeneous)) {
    self_transition_computation(markov->homogeneity);
  }

  build_transition_count(markov->order);
  build_observation_histogram();

  markov->characteristic_computation(*this , true);

  // calcul de la vraisemblance

  likelihood = markov->likelihood_computation(*this);
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Markov_data a partir d'un objet Markovian_sequences.
 *
 *  argument : reference sur un objet Markovian_sequences.
 *
 *--------------------------------------------------------------*/

Markov_data::Markov_data(const Markovian_sequences &seq)
:Markovian_sequences(seq)

{
  markov = 0;
  chain_data = 0;
  likelihood = D_INF;
  hidden_likelihood = D_INF;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Markov_data a partir d'un objet Markovian_sequences
 *  avec ajout d'une variable.
 *
 *  arguments : reference sur un objet Markovian_sequences, indice de la variable.
 *
 *--------------------------------------------------------------*/

Markov_data::Markov_data(const Markovian_sequences &seq , int variable)
:Markovian_sequences(seq , 'a' , variable)

{
  markov = 0;
  chain_data = 0;
  likelihood = D_INF;
  hidden_likelihood = D_INF;
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Markov_data.
 *
 *  arguments : reference sur un objet Markov_data,
 *              flag copie de l'objet Markov.
 *
 *--------------------------------------------------------------*/

void Markov_data::copy(const Markov_data &seq , bool model_flag)

{
  if ((model_flag) && (seq.markov)) {
    markov = new Markov(*(seq.markov) , false);
  }
  else {
    markov = 0;
  }

  if (seq.chain_data) {
    chain_data = new Chain_data(*(seq.chain_data));
  }
  else {
    chain_data = 0;
  }

  likelihood = seq.likelihood;
  hidden_likelihood = seq.hidden_likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Markov_data.
 *
 *--------------------------------------------------------------*/

Markov_data::~Markov_data()

{
  delete markov;
  delete chain_data;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Markov_data.
 *
 *  argument : reference sur un objet Markov_data.
 *
 *--------------------------------------------------------------*/

Markov_data& Markov_data::operator=(const Markov_data &seq)

{
  if (&seq != this) {
    delete markov;
    delete chain_data;

    remove();
    Sequences::remove();

    Sequences::copy(seq);
    Markovian_sequences::copy(seq);
    copy(seq);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'un histogramme.
 *
 *  arguments : reference sur un objet Format_error, type d'histogramme,
 *              variable, valeur (etat dans le cas des histogrammes d'observation).
 *
 *--------------------------------------------------------------*/

Distribution_data* Markov_data::extract(Format_error &error , int type ,
                                        int variable , int value) const

{
  bool status = true;
  Distribution *pdist;
  Parametric *pparam;
  Histogram *phisto;
  Distribution_data *histo;


  histo = 0;
  error.init();

  if (type == OBSERVATION) {
    if ((variable < 2) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }

    else {
      variable--;

      if ((value < 0) || (value >= marginal[0]->nb_value)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_STATE] << " " << value << " "
                      << SEQ_error[SEQR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      else {
        phisto = observation[variable][value];

        if (phisto->nb_element == 0) {
          status = false;
          error.update(STAT_error[STATR_EMPTY_HISTOGRAM]);
        }
      }
    }
  }

  else {
    if ((variable < 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }

    else {
      variable--;

      if ((value < 0) || (value >= marginal[variable]->nb_value)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[variable == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                      << value << " " << SEQ_error[SEQR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      else {
        switch (type) {
        case FIRST_OCCURRENCE :
          phisto = characteristics[variable]->first_occurrence[value];
          break;
        case RECURRENCE_TIME :
          phisto = characteristics[variable]->recurrence_time[value];
          break;
        case SOJOURN_TIME :
          phisto = characteristics[variable]->sojourn_time[value];
          break;
        case FINAL_RUN :
          phisto = characteristics[variable]->final_run[value];
          break;
        case NB_RUN :
          phisto = characteristics[variable]->nb_run[value];
          break;
        case NB_OCCURRENCE :
          phisto = characteristics[variable]->nb_occurrence[value];
          break;
        }

        if (phisto->nb_element == 0) {
          status = false;
          error.update(STAT_error[STATR_EMPTY_HISTOGRAM]);
        }
      }
    }
  }

  if (status) {
    pdist = 0;
    pparam = 0;

    switch (type) {
    case OBSERVATION :
      pdist = markov->process[variable]->observation[value];
      break;
    case FIRST_OCCURRENCE :
      pdist = markov->process[variable]->first_occurrence[value];
      break;
    case RECURRENCE_TIME :
      if (markov->process[variable]->recurrence_time) {
        pdist = markov->process[variable]->recurrence_time[value];
      }
      break;
    case SOJOURN_TIME :
      if (markov->process[variable]->sojourn_time) {
        pparam = markov->process[variable]->sojourn_time[value];
      }
      break;
    case NB_RUN :
      pdist = markov->process[variable]->nb_run[value];
      break;
    case NB_OCCURRENCE :
      pdist = markov->process[variable]->nb_occurrence[value];
      break;
    }

    if (pdist) {
      histo = new Distribution_data(*phisto , pdist);
    }
    else {
      histo = new Distribution_data(*phisto , pparam);
    }
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Markov_data.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Markov_data::ascii_write(ostream &os , bool exhaustive) const

{
  if (markov) {
    markov->ascii_write(os , this , exhaustive , false ,
                        ::test_hidden(markov->nb_output_process , markov->process));
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Markov_data dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Markov_data::ascii_write(Format_error &error , const char *path ,
                              bool exhaustive) const

{
  bool status = false;


  if (markov) {
    ofstream out_file(path);

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      markov->ascii_write(out_file , this , exhaustive , true ,
                          ::test_hidden(markov->nb_output_process , markov->process));
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Markov_data dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Markov_data::spreadsheet_write(Format_error &error , const char *path) const

{
  bool status = false;


  if (markov) {
    ofstream out_file(path);

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      markov->spreadsheet_write(out_file , this ,
                                ::test_hidden(markov->nb_output_process , markov->process));
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Markov_data.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Markov_data::plot_write(Format_error &error , const char *prefix ,
                             const char *title) const

{
  bool status = false;


  if (markov) {
    status = markov->plot_write(prefix , title , this);

    error.init();

    if (!status) {
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWDEFINE_COLLECTABLE(Markov_data , STATI_MARKOV_DATA);


RWspace Markov_data::binaryStoreSize() const

{
  RWspace size = Markovian_sequences::binaryStoreSize() + sizeof(likelihood) + sizeof(hidden_likelihood);

  size += sizeof(true);
  if (chain_data) {
    size += chain_data->binaryStoreSize();
  }

  if (markov) {
    size += markov->recursiveStoreSize();
  }

  return size;
}


void Markov_data::restoreGuts(RWvistream &is)

{
  bool status;


  delete markov;
  delete chain_data;

  Markovian_sequences::restoreGuts(is);

  is >> status;
  if (status) {
    chain_data = new Chain_data();
    chain_data->restoreGuts(is);
  }
  else {
    chain_data = 0;
  }

  is >> likelihood >> hidden_likelihood;

  is >> markov;
  if (markov == RWnilCollectable) {
    markov = 0;
  }
}


void Markov_data::restoreGuts(RWFile &file)

{
  bool status;


  delete markov;
  delete chain_data;

  Markovian_sequences::restoreGuts(file);

  file.Read(status);
  if (status) {
    chain_data = new Chain_data();
    chain_data->restoreGuts(file);
  }
  else {
    chain_data = 0;
  }

  file.Read(likelihood);
  file.Read(hidden_likelihood);

  file >> markov;
  if (markov == RWnilCollectable) {
    markov = 0;
  }
}


void Markov_data::saveGuts(RWvostream &os) const

{
  Markovian_sequences::saveGuts(os);

  if (chain_data) {
    os << true;
    chain_data->saveGuts(os);
  }
  else {
    os << false;
  }

  os << likelihood << hidden_likelihood;

  if (markov) {
    os << markov;
  }
  else {
    os << RWnilCollectable;
  }
}


void Markov_data::saveGuts(RWFile &file) const

{
  Markovian_sequences::saveGuts(file);

  if (chain_data) {
    file.Write(true);
    chain_data->saveGuts(file);
  }
  else {
    file.Write(false);
  }

  file.Write(likelihood);
  file.Write(hidden_likelihood);

  if (markov) {
    file << markov;
  }
  else {
    file << RWnilCollectable;
  }
} */
