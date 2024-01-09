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
 *       $Id: nonhomogeneous_markov.cpp 3257 2007-06-06 12:56:12Z dufourko $
 *
 *       Forum for V-Plants developers: amldevlp@cirad.fr
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
#include "tool/config.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/regression.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

#include "sequences.h"
#include "nonhomogeneous_markov.h"
#include "sequence_label.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Ecriture d'une courbe et des residus standardises correspondants au format Gnuplot.
 *
 *  arguments : path, pointeur sur les residus standardises.
 *
 *--------------------------------------------------------------*/

bool Curves::plot_print_standard_residual(const char *path , double *standard_residual) const

{
  bool status = false;
  register int i;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    // ecriture des reponses observees et des residus standardises

    for (i = 0;i < length;i++) {
      if (frequency[i] > 0) {
        out_file << i << " " << point[0][i];
        if (standard_residual) {
          out_file << " " << standard_residual[i];
        }
        out_file << " " << frequency[i] << endl;
      }
    }
  }

  return status;
}


};  // namespace stat_tool



using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Function.
 *
 *--------------------------------------------------------------*/

Function::Function()

{
  residual = NULL;
  frequency = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Function.
 *
 *  arguments : identificateur, longueur des courbes, parametres.
 *
 *--------------------------------------------------------------*/

Function::Function(int iident , int length , double *iparameter)
:RegressionKernel(iident , 0 , length - 1)

{
  register int i;


  for (i = 0;i < nb_parameter;i++) {
    parameter[i] = iparameter[i];
  }

  residual = NULL;
  frequency = NULL;

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
:RegressionKernel(iident , 0 , length - 1)

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
    residual = NULL;
    frequency = NULL;
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
  RegressionKernel::copy(function);
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
    RegressionKernel::remove();

    RegressionKernel::copy(function);
    copy(function);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Analyse du format d'un objet Function.
 *
 *  arguments : reference sur un objet StatError, stream,
 *              reference sur l'indice de la ligne lue, longueur,
 *              bornes sur les valeurs de la fonction.
 *
 *--------------------------------------------------------------*/

Function* function_parsing(StatError &error , ifstream &in_file , int &line ,
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


  function = NULL;

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
         *standard_residual , square_sum[3];
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

      // calcul des residus standardises

      standard_residual = new double[max_value + 1];

      for (i = 0;i <= max_value;i++) {
        if (frequency[i] > 0) {
          standard_residual[i] = residual[i] / residual_standard_deviation;
        }
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

    // ecriture des reponses observees et theoriques, des residus,
    // des residus standardises et des frequences

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

    for (i = 0;i <= max_value;i++) {
      if (file_flag) {
        os << "# ";
      }
      os << setw(width[0]) << i;

      if (curves) {
        if (frequency[i] > 0) {
          os << setw(width[1]) << curves->point[0][i];
        }
        else {
          os << setw(width[1]) << " ";
        }
      }

      os << setw(width[2]) << point[i];

      if (curves) {
        if (frequency[i] > 0) {
          os << setw(width[3]) << residual[i];
          os << setw(width[4]) << standard_residual[i];
        }
        else {
          os << setw(width[3]) << " ";
          os << setw(width[4]) << " ";
        }

        os << setw(width[5]) << frequency[i];
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
  double self_transition_mean , residual_mean , residual_standard_deviation , square_sum[3];


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

  // ecriture des reponses observees et theoriques, des residus,
  // des residus standardises et des frequences

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

  for (i = 0;i <= max_value;i++) {
    os << i;

    if (curves) {
      os << "\t";
      if (frequency[i] > 0) {
        os << curves->point[0][i];
      }
    }

    os << "\t" << point[i];

    if (curves) {
      if (frequency[i] > 0) {
        os << "\t" << residual[i] << "\t" << residual[i] / residual_standard_deviation;
      }
      else {
        os << "\t\t";
      }

      os << "\t" << (1. - point[i]) / residual_standard_deviation
         << "\t" << -point[i] / residual_standard_deviation << "\t" << frequency[i];
    }

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
  ofstream out_file(path);


  if (out_file) {
    status = true;

    for (i = 0;i <= max_value - min_value;i++) {
      out_file << point[i];
      if (residual_standard_deviation != D_DEFAULT) {
        out_file << " " << (1. - point[i]) / residual_standard_deviation
                 << " " << -point[i] / residual_standard_deviation;
      }
      out_file << endl;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe NonhomogeneousMarkov.
 *
 *--------------------------------------------------------------*/

NonhomogeneousMarkov::NonhomogeneousMarkov()

{
  markov_data = NULL;

  homogeneity = NULL;
  self_transition = NULL;

  process = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe NonhomogeneousMarkov.
 *
 *  arguments : nombre d'etats, identificateurs evolution
 *              des probabilites de transition.
 *
 *--------------------------------------------------------------*/

NonhomogeneousMarkov::NonhomogeneousMarkov(int inb_state , int *ident)
:Chain('o' , inb_state)

{
  register int i;


  markov_data = NULL;

  homogeneity = new bool[nb_state];
  self_transition = new Function*[nb_state];

  for (i = 0;i < nb_state;i++) {
    if (ident[i] == I_DEFAULT) {
      homogeneity[i] = true;
    }
    else {
      homogeneity[i] = false;
    }
    self_transition[i] = NULL;
  }

  process = new CategoricalSequenceProcess(nb_state , nb_state);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe NonhomogeneousMarkov.
 *
 *  arguments : pointeur sur un objet Chain, pointeurs sur des objets Function,
 *              longueur des sequences.
 *
 *--------------------------------------------------------------*/

NonhomogeneousMarkov::NonhomogeneousMarkov(const Chain *pchain , const Function **pself_transition ,
                                           int length)
:Chain(*pchain)

{
  register int i;


  markov_data = NULL;

  homogeneity = new bool[nb_state];
  self_transition = new Function*[nb_state];

  for (i = 0;i < nb_state;i++) {
    if (pself_transition[i]) {
      homogeneity[i] = false;
      self_transition[i] = new Function(*pself_transition[i]);
    }

    else {
      homogeneity[i] = true;
      self_transition[i] = NULL;
    }
  }

  process = new CategoricalSequenceProcess(nb_state , nb_state);

  characteristic_computation(length , true);
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet NonhomogeneousMarkov.
 *
 *  arguments : reference sur un objet NonhomogeneousMarkov,
 *              flag copie de l'objet NonhomogeneousMarkovData,
 *              flag copie des lois caracteristiques.
 *
 *--------------------------------------------------------------*/

void NonhomogeneousMarkov::copy(const NonhomogeneousMarkov &markov , bool data_flag ,
                                bool characteristic_flag)

{
  register int i;


  if ((data_flag) && (markov.markov_data)) {
    markov_data = new NonhomogeneousMarkovData(*(markov.markov_data) , false);
  }
  else {
    markov_data = NULL;
  }

  homogeneity = new bool[nb_state];
  self_transition = new Function*[nb_state];

  for (i = 0;i < nb_state;i++) {
    homogeneity[i] = markov.homogeneity[i];
    if (homogeneity[i]) {
      self_transition[i] = NULL;
    }
    else {
      self_transition[i] = new Function(*(markov.self_transition[i]));
    }
  }

  process = new CategoricalSequenceProcess(*(markov.process) , 'c' , characteristic_flag);
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet NonhomogeneousMarkov.
 *
 *--------------------------------------------------------------*/

void NonhomogeneousMarkov::remove()

{
  register int i;


  delete markov_data;

  if (self_transition) {
    for (i = 0;i < nb_state;i++) {
      if (!homogeneity[i]) {
        delete self_transition[i];
      }
    }
    delete [] self_transition;
  }

  delete [] homogeneity;

  delete process;
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe NonhomogeneousMarkov.
 *
 *--------------------------------------------------------------*/

NonhomogeneousMarkov::~NonhomogeneousMarkov()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe NonhomogeneousMarkov.
 *
 *  argument : reference sur un objet NonhomogeneousMarkov.
 *
 *--------------------------------------------------------------*/

NonhomogeneousMarkov& NonhomogeneousMarkov::operator=(const NonhomogeneousMarkov &markov)

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
 *  arguments : reference sur un objet StatError, type de loi, etat.
 *
 *--------------------------------------------------------------*/

DiscreteParametricModel* NonhomogeneousMarkov::extract(StatError &error , int type ,
                                                       int state) const

{
  bool status = true;
  Distribution *pdist;
  DiscreteParametric *pparam;
  DiscreteParametricModel *dist;
  FrequencyDistribution *phisto;


  dist = NULL;
  error.init();

  pdist = NULL;
  pparam = NULL;

  if ((state < 0) || (state >= process->nb_value)) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_STATE] << " " << state << " "
                  << STAT_error[STATR_NOT_PRESENT];
    error.update((error_message.str()).c_str());
  }

  else {
    switch (type) {
    case FIRST_OCCURRENCE :
      pdist = process->first_occurrence[state];
      break;
    case RECURRENCE_TIME :
      if (process->recurrence_time) {
        pdist = process->recurrence_time[state];
      }
      break;
    case SOJOURN_TIME :
      if (process->sojourn_time) {
        pparam = process->sojourn_time[state];
      }
      break;
    case NB_RUN :
      pdist = process->nb_run[state];
      break;
    case NB_OCCURRENCE :
      pdist = process->nb_occurrence[state];
      break;
    }

    if ((!pdist) && (!pparam)) {
      status = false;
      error.update(SEQ_error[SEQR_NON_EXISTING_CHARACTERISTIC_DISTRIBUTION]);
    }
  }

  if (status) {
    phisto = NULL;

    if (markov_data) {
      switch (type) {

      case FIRST_OCCURRENCE : {
        phisto = markov_data->characteristics[0]->first_occurrence[state];
        break;
      }

      case RECURRENCE_TIME : {
        if (markov_data->characteristics[0]->recurrence_time[state]->nb_element > 0) {
          phisto = markov_data->characteristics[0]->recurrence_time[state];
        }
        break;
      }

      case SOJOURN_TIME : {
        if (markov_data->characteristics[0]->sojourn_time[state]->nb_element > 0) {
          phisto = markov_data->characteristics[0]->sojourn_time[state];
        }
        break;
      }

      case NB_RUN : {
        phisto = markov_data->characteristics[0]->nb_run[state];
        break;
      }

      case NB_OCCURRENCE : {
        phisto = markov_data->characteristics[0]->nb_occurrence[state];
        break;
      }
      }
    }

    if (pdist) {
      dist = new DiscreteParametricModel(*pdist , phisto);
    }
    else if (pparam) {
      dist = new DiscreteParametricModel(*pparam , phisto);
    }
  }

  return dist;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet NonhomogeneousMarkov a partir d'un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              longueur des sequences.
 *
 *--------------------------------------------------------------*/

NonhomogeneousMarkov* nonhomogeneous_markov_ascii_read(StatError &error , const char *path ,
                                                       int length)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status , lstatus;
  register int i , j;
  int line , homogeneity , nb_state;
  double index;
  const Chain *chain;
  const Function **self_transition;
  NonhomogeneousMarkov *markov;
  ifstream in_file(path);


  markov = NULL;
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

        // test mot cle NONHOMOGENEOUS_MARKOV_CHAIN

        if (i == 0) {
          if (token != SEQ_word[SEQW_NONHOMOGENEOUS_MARKOV_CHAIN]) {
            status = false;
            error.update(STAT_parsing[STATP_KEY_WORD] , line);
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

    chain = chain_parsing(error , in_file , line , 'o');

    if (chain) {
      nb_state = chain->nb_state;
      self_transition = new const Function*[nb_state];
      for (i = 0;i < nb_state;i++) {
        self_transition[i] = NULL;
      }

      // analyse format des fonctions d'evolution
      // des probabilites de rester dans un etat

      for (i = 0;i < nb_state;i++) {
        homogeneity = I_DEFAULT;

        while (buffer.readLine(in_file , false)) {
          line++;

#         ifdef DEBUG
          cout << line << "  " << buffer << endl;
#         endif

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

            // test mot cle HOMOGENEOUS / NONHOMOGENEOUS

            case 2 : {
              if (token == SEQ_word[SEQW_HOMOGENEOUS]) {
                homogeneity = true;
              }
              else {
                if (token == SEQ_word[SEQW_NONHOMOGENEOUS]) {
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
        markov = new NonhomogeneousMarkov(chain , self_transition , length);
      }

      delete chain;

      for (i = 0;i < nb_state;i++) {
        delete self_transition[i];
      }
      delete [] self_transition;
    }
  }

  return markov;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet NonhomogeneousMarkov.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& NonhomogeneousMarkov::line_write(ostream &os) const

{
  os << nb_state << " " << STAT_word[STATW_STATES];

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet NonhomogeneousMarkov et de la structure de donnees associee.
 *
 *  arguments : stream, pointeur sur un objet NonhomogeneousMarkovData,
 *              flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& NonhomogeneousMarkov::ascii_write(ostream &os , const NonhomogeneousMarkovData *seq ,
                                           bool exhaustive , bool file_flag) const

{
  register int i;


  os << SEQ_word[SEQW_NONHOMOGENEOUS_MARKOV_CHAIN] << endl;
 
  // ecriture des parametres de la chaine de Markov

  ascii_print(os , file_flag);

  // ecriture des parametres des fonctions d'evolution
  // des probabilites de rester dans un etat

  for (i = 0;i < nb_state;i++) {
    os << "\n" << STAT_word[STATW_STATE] << " " << i << " ";

    switch (homogeneity[i]) {
    case true :
      os << SEQ_word[SEQW_HOMOGENEOUS] << endl;
      break;
    case false :
      os << SEQ_word[SEQW_NONHOMOGENEOUS] << endl;
      self_transition[i]->ascii_print(os , exhaustive , file_flag ,
                                      (seq ? seq->self_transition[i] : NULL));
      break;
    }
  }

  process->ascii_print(os , 0 , NULL , NULL , (seq ? seq->characteristics[0] : NULL) ,
                       exhaustive , file_flag);

  if (seq) {
    int nb_parameter = nb_parameter_computation();
    double information , likelihood;


    // ecriture de la loi empirique des longueurs des sequences

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    seq->length_distribution->ascii_characteristic_print(os , false , file_flag);

    if (exhaustive) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << "   | " << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      seq->length_distribution->ascii_print(os , file_flag);
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

    likelihood = seq->likelihood;

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   ("
       << STAT_label[STATL_NORMALIZED] << ": " << likelihood / seq->cumul_length << ")" << endl;

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

# ifdef DEBUG
  MultiPlotSet *plot_set;

  plot_set = get_plotable(seq);
  delete plot_set;
# endif

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet NonhomogeneousMarkov.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& NonhomogeneousMarkov::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , markov_data , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet NonhomogeneousMarkov dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool NonhomogeneousMarkov::ascii_write(StatError &error , const char *path ,
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
 *  Ecriture d'un objet NonhomogeneousMarkov et de la structure de donnees associee
 *  dans un fichier au format tableur.
 *
 *  arguments : stream, pointeur sur un objet NonhomogeneousMarkovData.
 *
 *--------------------------------------------------------------*/

ostream& NonhomogeneousMarkov::spreadsheet_write(ostream &os , const NonhomogeneousMarkovData *seq) const

{
  register int i;


  os << SEQ_word[SEQW_NONHOMOGENEOUS_MARKOV_CHAIN] << endl;

  // ecriture des parametres de la chaine de Markov

  spreadsheet_print(os);

  // ecriture des parametres des fonctions d'evolution
  // des probabilites de rester dans un etat

  for (i = 0;i < nb_state;i++) {
    os << "\n" << STAT_word[STATW_STATE] << "\t" << i << "\t";

    switch (homogeneity[i]) {
    case true :
      os << SEQ_word[SEQW_HOMOGENEOUS] << endl;
      break;
    case false :
      os << SEQ_word[SEQW_NONHOMOGENEOUS] << endl;
      self_transition[i]->spreadsheet_print(os , (seq ? seq->self_transition[i] : NULL));
      break;
    }
  }

  process->spreadsheet_print(os , 0 , NULL , NULL , (seq ? seq->characteristics[0] : NULL));

  if (seq) {
    int nb_parameter = nb_parameter_computation();
    double information , likelihood;


    // ecriture de la loi empirique des longueurs des sequences

    os << "\n" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    seq->length_distribution->spreadsheet_characteristic_print(os);

    os << "\n\t" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    seq->length_distribution->spreadsheet_print(os);

    os << "\n" << SEQ_label[SEQL_CUMUL_LENGTH] << "\t" << seq->cumul_length << endl;

    // ecriture de la quantite d'information des sequences dans le cas iid

    information = seq->iid_information_computation();

    os << "\n" << SEQ_label[SEQL_IID_INFORMATION] << "\t" << information << "\t"
       << information / seq->cumul_length << endl;

    // ecriture des vraisemblances des sequences

    likelihood = seq->likelihood;

    os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << likelihood << "\t"
       << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / seq->cumul_length << endl;

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
 *  Ecriture d'un objet NonhomogeneousMarkov dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool NonhomogeneousMarkov::spreadsheet_write(StatError &error , const char *path) const

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
 *  Sortie Gnuplot d'un objet NonhomogeneousMarkov et de la structure de donnees associee.
 *
 *  arguments : prefixe des fichiers, titre des figures,
 *              pointeur sur les sequences observees.
 *
 *--------------------------------------------------------------*/

bool NonhomogeneousMarkov::plot_write(const char *prefix , const char *title ,
                                      const NonhomogeneousMarkovData *seq) const

{
  bool status;
  register int i , j;
  int variable , start , *pfrequency , max_frequency[NB_STATE];
  double residual_mean , residual_standard_deviation , *standard_residual , *presidual ,
         min_standard_residual[NB_STATE] , max_standard_residual[NB_STATE];
  ostringstream data_file_name[NB_STATE * 2];


  if (seq) {
    status = process->plot_print(prefix , title , 0 , NULL , NULL ,
                                 seq->characteristics[0] , seq->length_distribution);
  }
  else {
    status = process->plot_print(prefix , title , 0);
  }

  if (status) {

    // ecriture des fichiers de donnees

    for (i = 0;i < nb_state;i++) {
      if (!homogeneity[i]) {
        if (seq) {
          max_frequency[i] = seq->self_transition[i]->max_frequency_computation();

          // calcul des residus standardises

          residual_mean = self_transition[i]->residual_mean_computation();
          residual_standard_deviation = sqrt(self_transition[i]->residual_variance_computation(residual_mean));

          standard_residual = new double[self_transition[i]->max_value + 1];

          pfrequency = self_transition[i]->frequency;
          presidual = self_transition[i]->residual;
          min_standard_residual[i] = 0.;
          max_standard_residual[i] = 0.;

          for (j = 0;j <= self_transition[i]->max_value;j++) {
            if (*pfrequency++ > 0) {
              standard_residual[j] = *presidual / residual_standard_deviation;
              if (standard_residual[j] < min_standard_residual[i]) {
                min_standard_residual[i] = standard_residual[j];
              }
              if (standard_residual[j] > max_standard_residual[i]) {
                max_standard_residual[i] = standard_residual[j];
              }
              presidual++;
            }
          }
        }

        data_file_name[i * 2] << prefix << i * 2 << ".dat";
        self_transition[i]->plot_print((data_file_name[i * 2].str()).c_str() ,
                                       (seq ? residual_standard_deviation : D_DEFAULT));

        if (seq) {
          data_file_name[i * 2 + 1] << prefix << i * 2 + 1 << ".dat";
          seq->self_transition[i]->plot_print_standard_residual((data_file_name[i * 2 + 1].str()).c_str() ,
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
            if ((int)(max_frequency[j] * YSCALE) + 1 < TIC_THRESHOLD) {
              out_file << "set ytics 0,1" << endl;
            }

            out_file << "plot [0:" << seq->self_transition[j]->length - 1
                     << "] [0:" << (int)(max_frequency[j] * YSCALE) + 1 << "] \""
                     << label((data_file_name[j * 2 + 1].str()).c_str())
                     << "\" using 1:4 title \"" << STAT_label[STATL_STATE] << " "
                     << j << " - "<< SEQ_label[SEQL_TRANSITION_COUNTS]
                     << "\" with impulses" << endl;

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

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet NonhomogeneousMarkov.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool NonhomogeneousMarkov::plot_write(StatError &error , const char *prefix ,
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
 *  Sortie graphique d'un objet NonhomogeneousMarkov et
 *  de la structure de donnees associee.
 *
 *  argument : pointeur sur les sequences observees.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* NonhomogeneousMarkov::get_plotable(const NonhomogeneousMarkovData *seq) const

{
  register int i , j , k;
  int nb_plot_set , index_length , index , nb_plot , max_frequency , *pfrequency;
  double residual_mean , residual_standard_deviation , min_standard_residual ,
         max_standard_residual , *standard_residual , *presidual;
  ostringstream legend;
  FrequencyDistribution *length_distribution;
  SequenceCharacteristics *characteristics;
  MultiPlotSet *plot_set;


  if (seq) {
    characteristics = seq->characteristics[0];
    length_distribution = seq->length_distribution;
  }
  else {
    characteristics = NULL;
    length_distribution = NULL;
  }

  // calcul du nombre de vues

  nb_plot_set = 0;

  if ((process->index_value) || (characteristics)) {
    nb_plot_set++;

    if (characteristics) {
      index_length = characteristics->index_value->plot_length_computation();

      if (characteristics->index_value->frequency[index_length - 1] < MAX_FREQUENCY) {
        nb_plot_set++;
      }
      nb_plot_set++;
    }
  }

  if ((process->first_occurrence) || (characteristics)) {
    for (i = 0;i < nb_state;i++) {
      if ((process->first_occurrence) && (process->first_occurrence[i])) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->first_occurrence[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }
  }

  if ((process->recurrence_time) || (characteristics)) {
    for (i = 0;i < nb_state;i++) {
      if ((process->recurrence_time) && (process->recurrence_time[i])) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->recurrence_time[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }
  }

  if ((process->sojourn_time) || (characteristics)) {
    for (i = 0;i < nb_state;i++) {
      if ((process->sojourn_time) && (process->sojourn_time[i])) {
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

      if ((characteristics) && (i < characteristics->nb_value) &&
          (characteristics->final_run[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }
  }

  if ((process->nb_run) || (process->nb_occurrence) ||
      ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence))) {
    for (i = 0;i < nb_state;i++) {
      if (process->nb_run) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->nb_run) && (characteristics->nb_run[i]->nb_element > 0)) {
        nb_plot_set++;
      }

      if (process->nb_occurrence) {
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

  for (i = 0;i < nb_state;i++) {
    if (!homogeneity[i]) {
      nb_plot_set++;
      if (seq) {
        nb_plot_set += 2;
      }
    }
  }

  plot_set = new MultiPlotSet(nb_plot_set , 1);

  MultiPlotSet &plot = *plot_set;

  plot.border = "15 lw 0";

  plot.variable_nb_viewpoint[0] = 1;

  index = 0;
  for (i = 0;i < nb_state;i++) {
    if (!homogeneity[i]) {

      // vue : probabilites de rester dans l'etat i indexees

      plot.variable[index] = 0;
      plot.viewpoint[index] = SELF_TRANSITION;

      plot[index].xrange = Range(0 , self_transition[i]->max_value);
      plot[index].yrange = Range(0. , 1.);

      if (self_transition[i]->max_value < TIC_THRESHOLD) {
        plot[index].xtics = 1;
      }

      if (seq) {
        plot[index].resize(2);

        legend.str("");
        legend << STAT_label[STATL_STATE] << " " << i << " - "
               << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION];
        plot[index][0].legend = legend.str();

        plot[index][0].style = "points";

        seq->self_transition[i]->plotable_write(0 , plot[index][0]);
        j = 1;
      }

      else {
        plot[index].resize(1);
        j = 0;
      }

      legend.str("");
      legend << STAT_label[STATL_STATE] << " " << i << " - "
             << SEQ_label[SEQL_THEORETICAL] << " " << SEQ_label[SEQL_SELF_TRANSITION];
      plot[index][j].legend = legend.str();

      plot[index][j].style = "linespoint";

      self_transition[i]->plotable_write(plot[index][j]);
      index++;

      if (seq) {

        // calcul des residus standardises

        residual_mean = self_transition[i]->residual_mean_computation();
        residual_standard_deviation = sqrt(self_transition[i]->residual_variance_computation(residual_mean));

        standard_residual = new double[self_transition[i]->max_value + 1];

        pfrequency = self_transition[i]->frequency;
        presidual = self_transition[i]->residual;
        min_standard_residual = 0.;
        max_standard_residual = 0.;

        for (j = 0;j <= self_transition[i]->max_value;j++) {
          if (*pfrequency++ > 0) {
            standard_residual[j] = *presidual / residual_standard_deviation;
            if (standard_residual[j] < min_standard_residual) {
              min_standard_residual = standard_residual[j];
            }
            if (standard_residual[j] > max_standard_residual) {
              max_standard_residual = standard_residual[j];
            }
            presidual++;
          }
        }

        // vue : residus standardises

        plot.variable[index] = 0;
        plot.viewpoint[index] = SELF_TRANSITION;

        plot[index].xrange = Range(0 , seq->self_transition[i]->length - 1);
        if (seq->self_transition[i]->length - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }

        plot[index].yrange = Range(min_standard_residual , max_standard_residual);

        plot[index].xlabel = SEQ_label[SEQL_INDEX];
        plot[index].ylabel = STAT_label[STATL_STANDARDIZED_RESIDUAL];

        nb_plot = 1;
        if (((1. - self_transition[i]->point[0]) / residual_standard_deviation <= max_standard_residual) ||
            ((1. - self_transition[i]->point[seq->self_transition[i]->length - 1]) / residual_standard_deviation <=
             max_standard_residual)) {
          nb_plot++;
        }
        if ((-self_transition[i]->point[0] / residual_standard_deviation >= min_standard_residual) ||
            (-self_transition[i]->point[seq->self_transition[i]->length - 1] / residual_standard_deviation >=
             min_standard_residual)) {
          nb_plot++;
        }
        plot[index].resize(nb_plot);

        plot[index][0].style = "points";

        pfrequency = self_transition[i]->frequency;
        for (j = 0;j <= self_transition[i]->max_value;j++) {
          if (*pfrequency++ > 0) {
            plot[index][0].add_point(j , standard_residual[j]);
          }
        }

        j = 1;
        if (((1. - self_transition[i]->point[0]) / residual_standard_deviation <= max_standard_residual) ||
            ((1. - self_transition[i]->point[seq->self_transition[i]->length - 1]) / residual_standard_deviation <=
             max_standard_residual)) {
          plot[index][j].legend = SEQ_label[SEQL_ASYMPTOTE];

          plot[index][j].style = "lines";

          for (k = 0;k <= self_transition[i]->max_value;k++) {
            plot[index][j].add_point(k , (1. - self_transition[i]->point[k]) / residual_standard_deviation);
          }
          j++;
        }

        if ((-self_transition[i]->point[0] / residual_standard_deviation >= min_standard_residual) ||
            (-self_transition[i]->point[seq->self_transition[i]->length - 1] / residual_standard_deviation >=
             min_standard_residual)) {
          plot[index][j].legend = SEQ_label[SEQL_ASYMPTOTE];

          plot[index][j].style = "lines";

          for (k = 0;k <= self_transition[i]->max_value;k++) {
            plot[index][j].add_point(k , -self_transition[i]->point[k] / residual_standard_deviation);
          }
        }
        index++;

        // vue : loi empirique des comptages de transition indexes

        plot.variable[index] = 0;
        plot.viewpoint[index] = SELF_TRANSITION;

        plot[index].xrange = Range(0 , seq->self_transition[i]->length - 1);
        max_frequency = seq->self_transition[i]->max_frequency_computation();
        plot[index].yrange = Range(0 , ceil(max_frequency * YSCALE));

        if (seq->self_transition[i]->length - 1 < TIC_THRESHOLD) {
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

        seq->self_transition[i]->plotable_frequency_write(plot[index][0]);
        index++;

        delete [] standard_residual;
      }
    }
  }

  process->plotable_write(*plot_set , index , 0 , NULL , NULL ,
                          characteristics , length_distribution);

  return plot_set;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet NonhomogeneousMarkov.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* NonhomogeneousMarkov::get_plotable() const

{
  return get_plotable(markov_data);
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de parametres independants d'un objet NonhomogeneousMarkov.
 *
 *--------------------------------------------------------------*/

int NonhomogeneousMarkov::nb_parameter_computation() const

{
  register int i;
  int nb_parameter = Chain::nb_parameter_computation();


  for (i = 0;i < nb_state;i++) {
    if (!homogeneity[i]) {
      nb_parameter += self_transition[i]->nb_parameter - 1;
    }
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe NonhomogeneousMarkovData.
 *
 *--------------------------------------------------------------*/

NonhomogeneousMarkovData::NonhomogeneousMarkovData()

{
  markov = NULL;
  chain_data = NULL;
  likelihood = D_INF;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe NonhomogeneousMarkovData.
 *
 *  argument : loi empirique des longueurs des sequences.
 *
 *--------------------------------------------------------------*/

NonhomogeneousMarkovData::NonhomogeneousMarkovData(const FrequencyDistribution &ilength_distribution)
:MarkovianSequences(ilength_distribution , 1 , NULL , false)

{
  markov = NULL;
  chain_data = NULL;
  likelihood = D_INF;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet NonhomogeneousMarkovData a partir
 *  d'un objet MarkovianSequences.
 *
 *  argument : reference sur un objet MarkovianSequences.
 *
 *--------------------------------------------------------------*/

NonhomogeneousMarkovData::NonhomogeneousMarkovData(const MarkovianSequences &seq)
:MarkovianSequences(seq)

{
  markov = NULL;
  chain_data = NULL;
  likelihood = D_INF;
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet NonhomogeneousMarkovData.
 *
 *  arguments : reference sur un objet NonhomogeneousMarkovData,
 *              flag copie de l'objet NonhomogeneousMarkov.
 *
 *--------------------------------------------------------------*/

void NonhomogeneousMarkovData::copy(const NonhomogeneousMarkovData &seq , bool model_flag)

{
  if ((model_flag) && (seq.markov)) {
    markov = new NonhomogeneousMarkov(*(seq.markov) , false);
  }
  else {
    markov = NULL;
  }

  if (seq.chain_data) {
    chain_data = new ChainData(*(seq.chain_data));
  }
  else {
    chain_data = NULL;
  }

  likelihood = seq.likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe NonhomogeneousMarkovData.
 *
 *--------------------------------------------------------------*/

NonhomogeneousMarkovData::~NonhomogeneousMarkovData()

{
  delete markov;
  delete chain_data;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe NonhomogeneousMarkovData.
 *
 *  argument : reference sur un objet NonhomogeneousMarkovData.
 *
 *--------------------------------------------------------------*/

NonhomogeneousMarkovData& NonhomogeneousMarkovData::operator=(const NonhomogeneousMarkovData &seq)

{
  if (&seq != this) {
    delete markov;
    delete chain_data;

    remove();
    Sequences::remove();

    Sequences::copy(seq);
    MarkovianSequences::copy(seq);
    copy(seq);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'une loi empirique.
 *
 *  arguments : reference sur un objet StatError, type de loi empirique, etat.
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData* NonhomogeneousMarkovData::extract(StatError &error , int type ,
                                                            int state) const

{
  bool status = true;
  Distribution *pdist;
  DiscreteParametric *pparam;
  FrequencyDistribution *phisto;
  DiscreteDistributionData *histo;


  histo = NULL;
  error.init();

  if ((state < 0) || (state >= marginal_distribution[0]->nb_value)) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_STATE] << " " << state << " "
                  << STAT_error[STATR_NOT_PRESENT];
    error.update((error_message.str()).c_str());
  }

  else {
    switch (type) {
    case FIRST_OCCURRENCE :
      phisto = characteristics[0]->first_occurrence[state];
      break;
    case RECURRENCE_TIME :
      phisto = characteristics[0]->recurrence_time[state];
      break;
    case SOJOURN_TIME :
      phisto = characteristics[0]->sojourn_time[state];
      break;
    case FINAL_RUN :
      phisto = characteristics[0]->final_run[state];
      break;
    case NB_RUN :
      phisto = characteristics[0]->nb_run[state];
      break;
    case NB_OCCURRENCE :
      phisto = characteristics[0]->nb_occurrence[state];
      break;
    }

    if (phisto->nb_element == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }
  }

  if (status) {
    pdist = NULL;
    pparam = NULL;

    switch (type) {
    case FIRST_OCCURRENCE :
      pdist = markov->process->first_occurrence[state];
      break;
    case RECURRENCE_TIME :
      if (markov->process->recurrence_time) {
        pdist = markov->process->recurrence_time[state];
      }
      break;
    case SOJOURN_TIME :
      if (markov->process->sojourn_time) {
        pparam = markov->process->sojourn_time[state];
      }
      break;
    case NB_RUN :
      pdist = markov->process->nb_run[state];
      break;
    case NB_OCCURRENCE :
      pdist = markov->process->nb_occurrence[state];
      break;
    }

    if (pdist) {
      histo = new DiscreteDistributionData(*phisto , pdist);
    }
    else {
      histo = new DiscreteDistributionData(*phisto , pparam);
    }
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Suppression du parametre d'index.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

NonhomogeneousMarkovData* NonhomogeneousMarkovData::remove_index_parameter(StatError &error) const

{
  NonhomogeneousMarkovData *seq;


  error.init();

  if (!index_parameter) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  else {
    seq = new NonhomogeneousMarkovData(*this , true , 'm');
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet NonhomogeneousMarkovData avec transformation du parametre d'index implicite
 *  en parametre d'index explicite.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

NonhomogeneousMarkovData* NonhomogeneousMarkovData::explicit_index_parameter(StatError &error) const

{
  NonhomogeneousMarkovData *seq;


  error.init();

  if (index_parameter) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  else {
    seq = new NonhomogeneousMarkovData(*this , true , 'e');
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet NonhomogeneousMarkovData.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& NonhomogeneousMarkovData::ascii_write(ostream &os , bool exhaustive) const

{
  if (markov) {
    markov->ascii_write(os , this , exhaustive , false);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet NonhomogeneousMarkovData dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool NonhomogeneousMarkovData::ascii_write(StatError &error , const char *path ,
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
      markov->ascii_write(out_file , this , exhaustive , true);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet NonhomogeneousMarkovData dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool NonhomogeneousMarkovData::spreadsheet_write(StatError &error , const char *path) const

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
      markov->spreadsheet_write(out_file , this);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet NonhomogeneousMarkovData.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool NonhomogeneousMarkovData::plot_write(StatError &error , const char *prefix ,
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
 *  Sortie graphique d'un objet NonhomogeneousMarkovData.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* NonhomogeneousMarkovData::get_plotable() const

{
  MultiPlotSet *plot_set;


  if (markov) {
    plot_set = markov->get_plotable(this);
  }
  else {
    plot_set = NULL;
  }

  return plot_set;
}


};  // namespace sequence_analysis
