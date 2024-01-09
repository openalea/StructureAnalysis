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
 *       $Id: discrete_parametric.cpp 17993 2015-04-23 06:46:53Z guedon $
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

#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"

#include "stat_tools.h"
#include "distribution.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {


extern bool cumul_matching_plot_print(const char *path , int nb_cumul , int *offset ,
                                      int *nb_value , double **cumul);



/*--------------------------------------------------------------*
 *
 *  Initialisation des parametres d'une loi discrete.
 *
 *  arguments : bornes inferieure et superieure, parametre, probabilite.
 *
 *--------------------------------------------------------------*/

void DiscreteParametric::init(int iinf_bound , int isup_bound ,
                              double iparameter , double iprobability)

{
  inf_bound = iinf_bound;
  sup_bound = isup_bound;
  parameter = iparameter;
  probability = iprobability;
}


/*--------------------------------------------------------------*
 *
 *  Initialisation de l'identificateur et des parametres d'une loi discrete.
 *
 *  arguments : identificateur, bornes inferieure et superieure,
 *              parametre, probabilite.
 *
 *--------------------------------------------------------------*/

void DiscreteParametric::init(int iident , int iinf_bound , int isup_bound ,
                              double iparameter , double iprobability)

{
  ident = iident;

  inf_bound = iinf_bound;
  sup_bound = isup_bound;
  parameter = iparameter;
  probability = iprobability;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe DiscreteParametric.
 *
 *  arguments : nombre de valeurs, identificateur, bornes inferieure et
 *              superieure, parametre, probabilite.
 *
 *--------------------------------------------------------------*/

DiscreteParametric::DiscreteParametric(int inb_value , int iident ,
                                       int iinf_bound , int isup_bound ,
                                       double iparameter , double iprobability)
:Distribution(inb_value)

{
  ident = iident;

  inf_bound = iinf_bound;
  sup_bound = isup_bound;
  parameter = iparameter;
  probability = iprobability;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe DiscreteParametric.
 *
 *  arguments : identificateur, bornes inferieure et superieure,
 *              parametre, probabilite, seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

DiscreteParametric::DiscreteParametric(int iident , int iinf_bound ,
                                       int isup_bound , double iparameter ,
                                       double iprobability , double cumul_threshold)

{
  ident = iident;

  inf_bound = iinf_bound;
  sup_bound = isup_bound;
  parameter = iparameter;
  probability = iprobability;

  nb_value = 0;

  if ((ident == BINOMIAL) || (ident == UNIFORM)) {
    nb_value = sup_bound + 1;
  }
  else {
    if ((ident == POISSON) || (ident == NEGATIVE_BINOMIAL)) {
      nb_value = (int)round(inf_bound + (parametric_mean_computation() - inf_bound +
                                         sqrt(parametric_variance_computation())) * 20.);
      if (nb_value == inf_bound) {
        nb_value++;
      }
    }
  }

  Distribution::init(nb_value);

  computation(1 , cumul_threshold);
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet DiscreteParametric a partir d'un objet Distribution.
 *
 *  arguments : reference sur un objet Distribution, nombre de valeurs allouees.
 *
 *--------------------------------------------------------------*/

DiscreteParametric::DiscreteParametric(const Distribution &dist , int ialloc_nb_value)
:Distribution(dist , 'c' , ialloc_nb_value)

{
  ident = CATEGORICAL;

  inf_bound = I_DEFAULT;
  sup_bound = I_DEFAULT;
  parameter = D_DEFAULT;
  probability = D_DEFAULT;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet DiscreteParametric a partir d'un objet Distribution
 *  avec changement d'echelle.
 *
 *  arguments : reference sur un objet Distribution, facteur d'echelle.
 *
 *--------------------------------------------------------------*/

DiscreteParametric::DiscreteParametric(const Distribution &dist , double scaling_coeff)
:Distribution(dist , scaling_coeff)

{
  ident = CATEGORICAL;

  inf_bound = I_DEFAULT;
  sup_bound = I_DEFAULT;
  parameter = D_DEFAULT;
  probability = D_DEFAULT;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'une loi inter-evenement a partir d'une loi
 *  inter-evenement initiale par dilatation/retraction de l'echelle des temps.
 *
 *  arguments : reference sur une loi inter-evenement, facteur d'echelle.
 *
 *--------------------------------------------------------------*/

DiscreteParametric::DiscreteParametric(const DiscreteParametric &dist , double scaling_coeff)
:Distribution((int)floor(dist.nb_value * scaling_coeff) + 1)

{
  double scaled_mean , scaled_variance , ratio , shifted_mean;


  // scaled_mean = scaling_coeff * dist.mean;
  // scaled_variance = scaling_coeff * scaling_coeff * dist.variance;
  scaled_mean = scaling_coeff * dist.parametric_mean_computation();
  scaled_variance = scaling_coeff * scaling_coeff * dist.parametric_variance_computation();

  inf_bound = (int)floor(dist.inf_bound * scaling_coeff);
  sup_bound = I_DEFAULT;
  parameter = D_DEFAULT;
  probability = D_DEFAULT;

  shifted_mean = scaled_mean - inf_bound;
  ratio = scaled_variance / shifted_mean;

  // cas binomiale

  if (ratio < 1. - POISSON_RANGE) {
    ident = BINOMIAL;
    sup_bound = (int)ceil(inf_bound + shifted_mean * shifted_mean /
                (shifted_mean - scaled_variance));
    if (sup_bound <= inf_bound) {
      sup_bound = inf_bound + 1;
    }
    probability = shifted_mean / (sup_bound - inf_bound);
  }

  // cas binomiale negative

  else if (ratio > 1. + POISSON_RANGE) {
    ident = NEGATIVE_BINOMIAL;
    parameter = shifted_mean * shifted_mean / (scaled_variance - shifted_mean);
    probability = shifted_mean / scaled_variance;
  }

  // cas Poisson

  else {
    ident = POISSON;
    parameter = shifted_mean;
  }

  computation();
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet DiscreteParametric a partir
 *  d'un objet FrequencyDistribution.
 *
 *  argument : reference sur un objet FrequencyDistribution.
 *
 *--------------------------------------------------------------*/

DiscreteParametric::DiscreteParametric(const FrequencyDistribution &histo)
:Distribution(histo)

{
  ident = CATEGORICAL;

  inf_bound = I_DEFAULT;
  sup_bound = I_DEFAULT;
  parameter = D_DEFAULT;
  probability = D_DEFAULT;
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet DiscreteParametric.
 *
 *  argument : reference sur un objet DiscreteParametric.
 *
 *--------------------------------------------------------------*/

void DiscreteParametric::copy(const DiscreteParametric &dist)

{
  ident = dist.ident;

  inf_bound = dist.inf_bound;
  sup_bound = dist.sup_bound;
  parameter = dist.parameter;
  probability = dist.probability;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe DiscreteParametric.
 *
 *  arguments : reference sur un objet DiscreteParametric, type de transformation
 *              ('c' : copie , 'n' : copie avec renormalisation),
 *              nombre de valeurs allouees.
 *
 *--------------------------------------------------------------*/

DiscreteParametric::DiscreteParametric(const DiscreteParametric &dist ,
                                       char transform , int ialloc_nb_value)
:Distribution(dist , transform , ialloc_nb_value)

{
  copy(dist);
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe DiscreteParametric.
 *
 *  argument : reference sur un objet DiscreteParametric.
 *
 *--------------------------------------------------------------*/

DiscreteParametric& DiscreteParametric::operator=(const DiscreteParametric &dist)

{
  if (&dist != this) {
    delete [] mass;
    delete [] cumul;

    Distribution::copy(dist);
    copy(dist);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Analyse du format d'un objet DiscreteParametric.
 *
 *  arguments : reference sur un objet StatError, stream,
 *              reference sur l'indice de la ligne lue, identificateur
 *              de la derniere loi dans la liste, seuil sur la fonction
 *              de repartition, borne inferieure minimum.
 *
 *--------------------------------------------------------------*/

DiscreteParametric* discrete_parametric_parsing(StatError &error , ifstream &in_file ,
                                                int &line , int last_ident ,
                                                double cumul_threshold , int min_inf_bound)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status = true , lstatus;
  register int i , j;
  int ident = I_DEFAULT;
  long inf_bound , sup_bound = I_DEFAULT;
  double parameter = D_DEFAULT , probability = D_DEFAULT;
  DiscreteParametric *dist;


  dist = NULL;

  while (buffer.readLine(in_file , false)) {
    line++;

#   ifdef DEBUG
    cout << line << " " << buffer << endl;
#   endif

    position = buffer.first('#');
    if (position != RW_NPOS) {
      buffer.remove(position);
    }
    i = 0;

    RWCTokenizer next(buffer);

    while (!((token = next()).isNull())) {

      // test nom de la loi

      if (i == 0) {
        for (j = BINOMIAL;j <= last_ident;j++) {
          if (token == STAT_discrete_distribution_word[j]) {
            ident = j;
            break;
          }
        }

        if (j == last_ident + 1) {
          status = false;
          error.update(STAT_parsing[STATP_DISTRIBUTION_NAME] , line , i + 1);
        }
      }

      // test nom du parametre

      else {
        switch ((i - 1) % 3) {

        case 0 : {
          switch ((i - 1) / 3) {

          // 1er parametre : borne inferieure

          case 0 : {
            if (token != STAT_word[STATW_INF_BOUND]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_INF_BOUND] , line , i + 1);
            }
            break;
          }

          // 2eme parametre : borne superieure (binomiale , uniforme)
          // ou parametre (Poisson , binomiale negative)

          case 1 : {
            if (((ident == BINOMIAL) || (ident == UNIFORM)) &&
                (token != STAT_word[STATW_SUP_BOUND])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_SUP_BOUND] , line , i + 1);
            }

            if (((ident == POISSON) || (ident == NEGATIVE_BINOMIAL)) &&
                (token != STAT_word[STATW_PARAMETER])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_PARAMETER] , line , i + 1);
            }
            break;
          }

          // 3eme parametre : probabilite (binomiale , binomiale negative)

          case 2 : {
            if (((ident == BINOMIAL) || (ident == NEGATIVE_BINOMIAL)) &&
                (token != STAT_word[STATW_PROBABILITY])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_PROBABILITY] , line , i + 1);
            }
            break;
          }
          }

          break;
        }

        // test separateur

        case 1 : {
          if (token != ":") {
            status = false;
            error.update(STAT_parsing[STATP_SEPARATOR] , line , i + 1);
          }
          break;
        }

        // test valeur du parametre

        case 2 : {
          switch ((i - 1) / 3) {

          // 1er parametre : borne inferieure

          case 0 : {
            lstatus = locale.stringToNum(token , &inf_bound);
            if ((lstatus) && ((inf_bound < min_inf_bound) || (inf_bound > MAX_INF_BOUND))) {
              lstatus = false;
            }
            break;
          }

          // 2eme parametre : borne superieure (binomiale , uniforme)
          // ou parametre (Poisson , binomiale negative)

          case 1 : {
            if ((ident == BINOMIAL) || (ident == UNIFORM)) {
              lstatus = locale.stringToNum(token , &sup_bound);

              if (lstatus) {
                switch (ident) {

                case BINOMIAL : {
                  if ((sup_bound <= inf_bound) || (sup_bound - inf_bound > MAX_DIFF_BOUND)) {
                    lstatus = false;
                  }
                  break;
                }

                case UNIFORM : {
                  if ((sup_bound < inf_bound) || (sup_bound - inf_bound > MAX_DIFF_BOUND)) {
                    lstatus = false;
                  }
                  break;
                }
                }
              }
            }

            if ((ident == POISSON) || (ident == NEGATIVE_BINOMIAL)) {
              lstatus = locale.stringToNum(token , &parameter);
              if ((lstatus) && ((parameter <= 0.) ||
                   ((ident == POISSON) && (parameter > MAX_MEAN)))) {
                lstatus = false;
              }
            }
            break;
          }

          // 3eme parametre : probabilite (binomiale , binomiale negative)

          case 2 : {
            if ((ident == BINOMIAL) || (ident == NEGATIVE_BINOMIAL)) {
              lstatus = locale.stringToNum(token , &probability);

              if (lstatus) {
                switch (ident) {

                case BINOMIAL : {
                  if ((probability < 0.) || (probability > 1.)) {
                    lstatus = false;
                  }
                  break;
                }

                case NEGATIVE_BINOMIAL : {
                  if ((probability <= 0.) || (probability >= 1.)) {
                    lstatus = false;
                  }
                  else if (parameter * (1. - probability) / probability > MAX_MEAN) {
                    lstatus = false;
                  }
                  break;
                }
                }
              }
            }
            break;
          }
          }

          if (!lstatus) {
            status = false;
            error.update(STAT_parsing[STATP_PARAMETER_VALUE] , line , i + 1);
          }
          break;
        }
        }
      }

      i++;
    }

    if (i > 0) {
      if ((((ident == BINOMIAL) || (ident == NEGATIVE_BINOMIAL)) && (i != 10)) ||
          (((ident == POISSON) || (ident == UNIFORM)) && (i != 7))) {
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
    dist = new DiscreteParametric(ident , inf_bound , sup_bound ,
                                  parameter , probability , cumul_threshold);
  }

  return dist;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des parametres d'une loi discrete.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& DiscreteParametric::ascii_print(ostream &os) const

{
  os << STAT_discrete_distribution_word[ident] << "   ";

  if (inf_bound != I_DEFAULT) {
    os << STAT_word[STATW_INF_BOUND] << " : " << inf_bound << "   ";
  }
  if (sup_bound != I_DEFAULT) {
    os << STAT_word[STATW_SUP_BOUND] << " : " << sup_bound << "   ";
  }
  if (parameter != D_DEFAULT) {
    os << STAT_word[STATW_PARAMETER] << " : " << parameter << "   ";
  }
  if (probability != D_DEFAULT) {
    os << STAT_word[STATW_PROBABILITY] << " : " << probability;
  }
  os << endl;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des parametres d'une loi discrete au format tableur.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& DiscreteParametric::spreadsheet_print(ostream &os) const

{
  os << STAT_discrete_distribution_word[ident] << "\t";

  if (inf_bound != I_DEFAULT) {
    os << STAT_word[STATW_INF_BOUND] << "\t" << inf_bound << "\t";
  }
  if (sup_bound != I_DEFAULT) {
    os << STAT_word[STATW_SUP_BOUND] << "\t" << sup_bound << "\t";
  }
  if (parameter != D_DEFAULT) {
    os << STAT_word[STATW_PARAMETER] << "\t" << parameter << "\t";
  }
  if (probability != D_DEFAULT) {
    os << STAT_word[STATW_PROBABILITY] << "\t" << probability;
  }
  os << endl;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des caracteristiques d'une loi.
 *
 *  arguments : stream, flag ecriture des parametres de forme, flag commentaire.
 *
 *--------------------------------------------------------------*/

ostream& DiscreteParametric::ascii_parametric_characteristic_print(ostream &os , bool shape , bool comment_flag) const

{
  if (ident == CATEGORICAL) {
    ascii_characteristic_print(os , shape , comment_flag);
  }

  else {
    double variance = parametric_variance_computation();


    if (comment_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_MEAN] << ": " << parametric_mean_computation() << "   "
       << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
       << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance) << endl;

    if ((shape) && (variance > 0.)) {
      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << parametric_skewness_computation() << "   "
         << STAT_label[STATL_KURTOSIS_COEFF] << ": " << parametric_kurtosis_computation() << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des caracteristiques d'une loi au format tableur.
 *
 *  arguments : stream, flag ecriture des parametres de forme.
 *
 *--------------------------------------------------------------*/

ostream& DiscreteParametric::spreadsheet_parametric_characteristic_print(ostream &os , bool shape) const

{
  if (ident == CATEGORICAL) {
    spreadsheet_characteristic_print(os , shape);
  }

  else {
    double variance = parametric_variance_computation();


    os << STAT_label[STATL_MEAN] << "\t" << parametric_mean_computation() << "\t"
       << STAT_label[STATL_VARIANCE] << "\t" << variance << "\t"
       << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(variance) << endl;

    if ((shape) && (variance > 0.)) {
      os << STAT_label[STATL_SKEWNESS_COEFF] << "\t" << parametric_skewness_computation() << "\t"
         << STAT_label[STATL_KURTOSIS_COEFF] << "\t" << parametric_kurtosis_computation() << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des parametres d'une loi discrete au format Gnuplot.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& DiscreteParametric::plot_title_print(ostream &os) const

{
  if (ident != CATEGORICAL) {
    os << " " << STAT_discrete_distribution_letter[ident] << "(";

    os << inf_bound;
    if (sup_bound != I_DEFAULT) {
      os << ", " << sup_bound;
    }
    if (parameter != D_DEFAULT) {
      os << ", " << parameter;
    }
    if (probability != D_DEFAULT) {
      os << ", " << probability;
    }
    os << ")";
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Visualisation d'une loi discrete parametrique.
 *
 *  arguments : stream, reference sur un objet DiscreteParametric.
 *
 *--------------------------------------------------------------*/

ostream& operator<<(ostream &os , const DiscreteParametric &dist)

{
  os.precision(5);

  os << endl;
  dist.ascii_print(os);
  dist.Distribution::print(os);

  os.precision(6);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de parametres d'une loi discrete.
 *
 *--------------------------------------------------------------*/

int DiscreteParametric::nb_parameter_computation()

{
  int bnb_parameter;


  switch (ident) {
  case BINOMIAL :
    bnb_parameter = 3;
    break;
  case POISSON :
    bnb_parameter = 2;
    break;
  case NEGATIVE_BINOMIAL :
    bnb_parameter = 3;
    break;
  case UNIFORM :
    bnb_parameter = 2;
    break;
  default :
    bnb_parameter = 0;
    break;
  }

  return bnb_parameter;
}


/*--------------------------------------------------------------*
 *
 *  Mise a jour du nombre de parametres d'une loi discrete.
 *
 *--------------------------------------------------------------*/

void DiscreteParametric::nb_parameter_update()

{
  nb_parameter = nb_parameter_computation();
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la moyenne d'une loi discrete parametrique.
 *
 *--------------------------------------------------------------*/

double DiscreteParametric::parametric_mean_computation() const

{
  double parametric_mean;


  switch (ident) {
  case BINOMIAL :
    parametric_mean = inf_bound + (sup_bound - inf_bound) * probability;
    break;
  case POISSON :
    parametric_mean = inf_bound + parameter;
    break;
  case NEGATIVE_BINOMIAL :
    parametric_mean = inf_bound + parameter * (1. - probability) / probability;
    break;
  case UNIFORM :
    parametric_mean = (double)(inf_bound + sup_bound) / 2.;
    break;
  default :
    parametric_mean = mean;
    break;
  }

  return parametric_mean;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la variance d'une loi discrete parametrique.
 *
 *--------------------------------------------------------------*/

double DiscreteParametric::parametric_variance_computation() const

{
  double parametric_variance;


  switch (ident) {
  case BINOMIAL :
    parametric_variance = (sup_bound - inf_bound) * probability * (1. - probability);
    break;
  case POISSON :
    parametric_variance = parameter;
    break;
  case NEGATIVE_BINOMIAL :
    parametric_variance = parameter * (1. - probability) / (probability * probability);
    break;
  case UNIFORM :
    parametric_variance = (double)((sup_bound - inf_bound + 2) *
                          (sup_bound - inf_bound)) / 12.;
    break;
  default :
    parametric_variance = variance;
    break;
  }

  return parametric_variance;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du coefficient d'asymetrie d'une loi discrete parametrique.
 *
 *--------------------------------------------------------------*/

double DiscreteParametric::parametric_skewness_computation() const

{
  double parametric_skewness;


  switch (ident) {

  case BINOMIAL : {
    if ((probability > 0.) && (probability < 1.)) {
      parametric_skewness = (1. - 2. * probability) /
                            sqrt((sup_bound - inf_bound) * probability * (1. - probability));
    }
    else {
      parametric_skewness = 0.;
    }
    break;
  }

  case POISSON : {
    parametric_skewness = 1. / sqrt(parameter);
    break;
  }

  case NEGATIVE_BINOMIAL : {
    parametric_skewness = (2. - probability) / sqrt(parameter * (1. - probability));
    break;
  }

  case UNIFORM : {
    parametric_skewness = 0.;
    break;
  }

  default : {
    parametric_skewness = D_INF;
    break;
  }
  }

  return parametric_skewness;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'exces d'applatissement d'une loi discrete parametrique:
 *  exces d'applatissement = coefficient d'applatissement - 3..
 *
 *--------------------------------------------------------------*/

double DiscreteParametric::parametric_kurtosis_computation() const

{
  double parametric_kurtosis;


  switch (ident) {

  case BINOMIAL : {
    if ((probability > 0.) && (probability < 1.)) {
      parametric_kurtosis = (1. - 6. * probability * (1. - probability)) /
                            ((sup_bound - inf_bound) * probability * (1. - probability));
    }
    else {
      parametric_kurtosis = -D_INF;
    }
    break;
  }

  case POISSON : {
    parametric_kurtosis = 1. / parameter;
    break;
  }

  case NEGATIVE_BINOMIAL : {
    parametric_kurtosis = (probability * probability + 6. * (1. - probability)) /
                          (parameter * (1. - probability));
    break;
  }

  case UNIFORM : {
    parametric_kurtosis = -0.5;
    break;
  }

  default : {
    parametric_kurtosis = -D_INF;
    break;
  }
  }

  return parametric_kurtosis;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la distance entre deux lois discretes (sup de la difference
 *  absolue des fonctions de repartition dans le cas de fonctions de repartition
 *  ne se croisant pas; sinon somme des sup avant et apres croisement
 *  de la difference des fonctions de repartition).
 *
 *--------------------------------------------------------------*/

double DiscreteParametric::sup_norm_distance_computation(const DiscreteParametric &dist) const

{
  bool crossing;
  register int i;
  int inf , sup;
  double buff , distance , max_absolute_diff;


  inf = MAX(offset , dist.offset);
  sup = MIN(nb_value , dist.nb_value) - 1;
  if (sup < inf) {
    distance = 1.;
  }

  else {
    distance = 0.;
    max_absolute_diff = fabs(cumul[inf] - dist.cumul[inf]);
    crossing = false;

    for (i = inf + 1;i <= sup;i++) {
      buff = fabs(cumul[i] - dist.cumul[i]);
      if (buff > max_absolute_diff) {
        max_absolute_diff = buff;
      }

      if ((!crossing) && (((cumul[i] > dist.cumul[i]) && (cumul[i - 1] <= dist.cumul[i - 1])) ||
          ((cumul[i] <= dist.cumul[i]) && (cumul[i - 1] > dist.cumul[i - 1])))) {
        crossing = true;
        distance = max_absolute_diff;
        max_absolute_diff = 0.;
      }
    }
    distance += max_absolute_diff;

#   ifdef DEBUG
    double overlap;

    overlap = 0.;
    for (i = inf;i <= sup;i++) {
      overlap += MIN(mass[i] , dist.mass[i]);
    }

    cout << "\nSup norm distance: " << distance << " " << 1. - overlap << endl;
#   endif
  }

  return distance;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet DiscreteParametricModel a partir
 *  d'un objet FrequencyDistribution.
 *
 *  argument : reference sur un objet FrequencyDistribution.
 *
 *--------------------------------------------------------------*/

DiscreteParametricModel::DiscreteParametricModel(const FrequencyDistribution &histo)
:DiscreteParametric(histo)

{
  frequency_distribution = new DiscreteDistributionData(histo);
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet DiscreteParametricModel a partir
 *  d'un objet Distribution et d'un objet FrequencyDistribution.
 *
 *  arguments : reference sur un objet Distribution,
 *              pointeur sur un objet FrequencyDistribution.
 *
 *--------------------------------------------------------------*/

DiscreteParametricModel::DiscreteParametricModel(const Distribution &dist ,
                                                 const FrequencyDistribution *histo)
:DiscreteParametric(dist)

{
  if ((histo) && (histo->nb_element > 0)) {
    frequency_distribution = new DiscreteDistributionData(*histo);
  }
  else {
    frequency_distribution = NULL;
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet DiscreteParametricModel a partir
 *  d'un objet DiscreteParametric et d'un objet FrequencyDistribution.
 *
 *  arguments : reference sur un objet DiscreteParametric,
 *              pointeur sur un objet FrequencyDistribution.
 *
 *--------------------------------------------------------------*/

DiscreteParametricModel::DiscreteParametricModel(const DiscreteParametric &dist ,
                                                 const FrequencyDistribution *histo)
:DiscreteParametric(dist)

{
  if ((histo) && (histo->nb_element > 0)) {
    frequency_distribution = new DiscreteDistributionData(*histo);
  }
  else {
    frequency_distribution = NULL;
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe DiscreteParametricModel.
 *
 *  arguments : reference sur un objet DiscreteParametricModel,
 *              flag copie de l'objet DiscreteDistributionData.
 *
 *--------------------------------------------------------------*/

DiscreteParametricModel::DiscreteParametricModel(const DiscreteParametricModel &dist ,
                                                 bool data_flag)
:DiscreteParametric(dist)

{
  if ((data_flag) && (dist.frequency_distribution)) {
    frequency_distribution = new DiscreteDistributionData(*(dist.frequency_distribution) , false);
  }
  else {
    frequency_distribution = NULL;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe DiscreteParametricModel.
 *
 *--------------------------------------------------------------*/

DiscreteParametricModel::~DiscreteParametricModel()

{
  if (frequency_distribution) {
    delete frequency_distribution; 
  }
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe DiscreteParametricModel.
 *
 *  argument : reference sur un objet DiscreteParametricModel.
 *
 *--------------------------------------------------------------*/

DiscreteParametricModel& DiscreteParametricModel::operator=(const DiscreteParametricModel &dist)

{
  if (&dist != this) {
    delete frequency_distribution;

    delete [] mass;
    delete [] cumul;

    Distribution::copy(dist);
    DiscreteParametric::copy(dist);

    if (dist.frequency_distribution) {
      frequency_distribution = new DiscreteDistributionData(*(dist.frequency_distribution) , false);
    }
    else {
      frequency_distribution = NULL;
    }
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de la partie "donnees" d'un objet DiscreteParametricModel.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData* DiscreteParametricModel::extract_data(StatError &error) const

{
  DiscreteDistributionData *histo;


  error.init();

  if (!frequency_distribution) {
    histo = NULL;
    error.update(STAT_error[STATR_NON_EXISTING_FREQUENCY_DISTRIBUTION]);
  }

  else {
    histo = new DiscreteDistributionData(*frequency_distribution);
    histo->distribution = new DiscreteParametricModel(*this , false);
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet DiscreteParametricModel a partir d'un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

DiscreteParametricModel* discrete_parametric_ascii_read(StatError &error , const char *path ,
                                                        double cumul_threshold)

{
  RWCString buffer;
  size_t position;
  bool status;
  int line;
  DiscreteParametric *pdist;
  DiscreteParametricModel *dist;
  ifstream in_file(path);


  dist = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;

    pdist = discrete_parametric_parsing(error , in_file , line ,
                                        UNIFORM , cumul_threshold);

    if (!pdist) {
      status = false;
    }

    while (buffer.readLine(in_file , false)) {
      line++;

#     ifdef DEBUG
      cout << line << " " << buffer << endl;
#     endif

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
      dist = new DiscreteParametricModel(*pdist);
    }

    delete pdist;
  }

  return dist;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet DiscreteParametricModel.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& DiscreteParametricModel::line_write(ostream &os) const

{
  os << STAT_discrete_distribution_word[ident];

  if (ident == CATEGORICAL) {
    if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
      os << "   " << STAT_label[STATL_MEAN] << ": " << mean
         << "   " << STAT_label[STATL_VARIANCE] << ": " << variance;
    }
  }

  else {
    os << "   " << STAT_label[STATL_MEAN] << ": " << parametric_mean_computation()
       << "   " << STAT_label[STATL_VARIANCE] << ": " << parametric_variance_computation();
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une loi discrete parametrique et d'une loi discrete empirique.
 *
 *  arguments : stream, pointeur sur une loi discrete empirique,
 *              flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& DiscreteParametricModel::ascii_write(ostream &os , const DiscreteDistributionData *histo ,
                                              bool exhaustive , bool file_flag) const

{
  if (ident == CATEGORICAL) {
    file_flag = false;
  }

  ascii_print(os);
  if (complement > 0.) {
    os << STAT_label[STATL_UNPROPER] << " " << STAT_label[STATL_DISTRIBUTION] << " ("
       << STAT_label[STATL_COMPLEMENTARY_PROBABILITY] << ": " << complement << ")" << endl;
  }

  ascii_parametric_characteristic_print(os , true , file_flag);

  if (file_flag) {
    os << "# ";
  }
  os << STAT_label[STATL_MEAN_ABSOLUTE_DEVIATION] << ": " << mean_absolute_deviation_computation();
  if (mean > 0.) {
    os << "   " << STAT_label[STATL_CONCENTRATION_COEFF] << ": " << concentration_computation();
  }
  os << endl;

  if (histo) {
    double likelihood , information;
    Test test(CHI2);


    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    histo->ascii_characteristic_print(os , true , file_flag);

    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_MEAN_ABSOLUTE_DEVIATION] << ": " << histo->mean_absolute_deviation_computation();
    if (histo->mean > 0.) {
      os << "   " << STAT_label[STATL_CONCENTRATION_COEFF] << ": " << histo->concentration_computation();
    }
    os << endl;

    likelihood = likelihood_computation(*histo);
    information = histo->information_computation();

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   ("
       << STAT_label[STATL_NORMALIZED] << ": " << likelihood / histo->nb_element << ")" << endl;

    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_MAX_LIKELIHOOD] << ": " << information << "   ("
       << STAT_label[STATL_INFORMATION] << ": " << information / histo->nb_element << ")" << endl;

    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_DEVIANCE] << ": " << 2 * (information - likelihood) << endl;

    chi2_fit(*histo , test);
    os << "\n";
    test.ascii_print(os , file_flag);
  }

  else {
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_INFORMATION] << ": " << information_computation() << endl;
  }

  if (exhaustive) {
    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "  ";
    if (histo) {
      os << " | " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    }
    os << " | " << STAT_label[STATL_DISTRIBUTION];
    if (histo) {
      os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
         << STAT_label[STATL_FUNCTION];
    }
    os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
       << STAT_label[STATL_FUNCTION] << endl;

    Distribution::ascii_print(os , file_flag , true , false , histo);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet DiscreteParametricModel.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& DiscreteParametricModel::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , frequency_distribution , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet DiscreteParametricModel dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool DiscreteParametricModel::ascii_write(StatError &error , const char *path ,
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
    ascii_write(out_file , frequency_distribution , exhaustive , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une loi discrete parametrique et d'une loi discrete empirique
 *  dans un fichier au format tableur.
 *
 *  arguments : stream, pointeur sur une loi discrete empirique.
 *
 *--------------------------------------------------------------*/

ostream& DiscreteParametricModel::spreadsheet_write(ostream &os ,
                                                    const DiscreteDistributionData *histo) const

{


  spreadsheet_print(os);
  if (complement > 0.) {
    os << STAT_label[STATL_UNPROPER] << " " << STAT_label[STATL_DISTRIBUTION] << "\t"
       << STAT_label[STATL_COMPLEMENTARY_PROBABILITY] << "\t" << complement << endl;
  }

  spreadsheet_parametric_characteristic_print(os , true);

  os << STAT_label[STATL_MEAN_ABSOLUTE_DEVIATION] << "\t" << mean_absolute_deviation_computation();
  if (mean > 0.) {
    os << "\t\t" << STAT_label[STATL_CONCENTRATION_COEFF] << "\t" << concentration_computation();
  }
  os << endl;

  if (histo) {
    double likelihood , information;
    Test test(CHI2);


    os << "\n" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    histo->spreadsheet_characteristic_print(os , true);

    os << STAT_label[STATL_MEAN_ABSOLUTE_DEVIATION] << "\t" << histo->mean_absolute_deviation_computation();
    if (histo->mean > 0.) {
      os << "\t\t" << STAT_label[STATL_CONCENTRATION_COEFF] << "\t" << histo->concentration_computation();
    }
    os << endl;

    likelihood = likelihood_computation(*histo);
    information = histo->information_computation();

    os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << likelihood << "\t"
       << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / histo->nb_element << endl;
    os << STAT_label[STATL_MAX_LIKELIHOOD] << "\t" << information << "\t"
       << STAT_label[STATL_INFORMATION] << "\t" << information / histo->nb_element << endl;
    os << STAT_label[STATL_DEVIANCE] << "\t" << 2 * (information - likelihood) << endl;

    chi2_fit(*histo , test);
    os << "\n";
    test.spreadsheet_print(os);
  }

  else {
    os << STAT_label[STATL_INFORMATION] << "\t" << information_computation() << endl;
  }

  os << "\n";
  if (histo) {
    os << "\t" << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
  }
  os << "\t" << STAT_label[STATL_DISTRIBUTION];
  if (histo) {
    os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
       << STAT_label[STATL_FUNCTION];
  }
  os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
     << STAT_label[STATL_FUNCTION];
  if (mean > 0.) {
    if ((histo) && (histo->mean > 0.)) {
      os << "\t" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_CONCENTRATION] << " "
         << STAT_label[STATL_FUNCTION] << "\t" << STAT_label[STATL_DISTRIBUTION] << " ";
    }
    else {
      os << "\t";
    }
    os << STAT_label[STATL_CONCENTRATION] << " " << STAT_label[STATL_FUNCTION];
  }
  os << endl;

  Distribution::spreadsheet_print(os , true , true , false , histo);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet DiscreteParametricModel dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool DiscreteParametricModel::spreadsheet_write(StatError &error , const char *path) const

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
    spreadsheet_write(out_file , frequency_distribution);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'une loi discrete parametrique et d'une loi discrete empirique.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures, pointeur sur une loi discrete empirique.
 *
 *--------------------------------------------------------------*/

bool DiscreteParametricModel::plot_write(StatError &error , const char *prefix , const char *title ,
                                         const DiscreteDistributionData *histo) const

{
  bool status;


  if (histo) {
    register int i;
    int plot_nb_value , *poffset , *pnb_value;
    double scale , **pcumul , **pconcentration;
    ostringstream *data_file_name;


    error.init();

    // ecriture des fichiers de donnees

    data_file_name = new ostringstream[3];

    poffset = new int[2];
    pnb_value = new int[2];
    pcumul = new double*[2];

    pconcentration = new double*[2];
    pconcentration[1] = NULL;

    data_file_name[0] << prefix << 0 << ".dat";

    poffset[0] = histo->offset;
    pnb_value[0] = histo->nb_value;

    // calcul des fonctions de repartition et de concentration de la loi discrete empirique

    scale = histo->nb_element / (1. - complement);
    pcumul[0] = histo->cumul_computation(scale);
    pconcentration[0] = histo->concentration_function_computation(scale);

    status = histo->plot_print((data_file_name[0].str()).c_str() , pcumul[0] , pconcentration[0]);

    if (status) {
      data_file_name[1] << prefix << 1 << ".dat";

      poffset[1] = offset;
      pnb_value[1] = nb_value;
      plot_nb_value = plot_nb_value_computation(histo);
      pcumul[1] = cumul;

      // calcul de la fonction de concentration

      pconcentration[1] = concentration_function_computation();
      plot_print((data_file_name[1].str()).c_str() , pconcentration[1] , scale);

      if ((variance > 0.) && (histo->variance > 0.)) {
        data_file_name[2] << prefix << 2 << ".dat";
        cumul_matching_plot_print((data_file_name[2].str()).c_str() , 2 , poffset , pnb_value , pcumul);
      }
    }

    delete [] poffset;
    delete [] pnb_value;

    delete [] pcumul[0];
    delete [] pcumul;

    delete [] pconcentration[0];
    delete [] pconcentration[1];
    delete [] pconcentration;

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

        // ajustement

        if (plot_nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }

        out_file << "plot [0:" << MAX(plot_nb_value - 1 , 1) << "] [0:"
                 << (int)(MAX(histo->max , max * scale) * YSCALE) + 1
                 << "] \"" << label((data_file_name[0].str()).c_str()) << "\" using 1:3 title \""
                 << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses,\\" << endl;
        out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using 1:2 title \""
                 << STAT_label[STATL_DISTRIBUTION];
        plot_title_print(out_file);
        out_file << "\" with linespoints" << endl;

        if (plot_nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }

        if (variance > 0.) {

          // fonctions de repartition

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;

          if (plot_nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }

          out_file << "plot [0:" << plot_nb_value - 1 << "] [0:" << 1. - complement << "] ";
          if (histo->variance > 0.) {
            out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using 1:5 title \""
                     << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
                     << STAT_label[STATL_FUNCTION] << "\" with linespoints,\\" << endl;
          }
          out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using 1:3 title \""
                   << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
                   << STAT_label[STATL_FUNCTION] << "\" with linespoints" << endl;

          if (plot_nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }

          if (histo->variance > 0.) {

            // mise en correspondance des fonctions de repartition en prenant
            // comme reference celle de la loi

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

            out_file << "plot [0:" << 1. - complement << "] [0:" << 1. - complement << "] \""
                     << label((data_file_name[2].str()).c_str())
                     << "\" using 2:1 notitle with linespoints,\\" << endl;
            out_file << "\"" << label((data_file_name[2].str()).c_str())
                     << "\" using 2:2 notitle with lines" << endl;

/*            out_file << "plot [0:" << 1. - complement << "] [0:" << 1. - complement << "] \""
                     << label((data_file_name[2].str()).c_str()) << "\" using 2:1 title \""
                     << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
                     << STAT_label[STATL_FUNCTION] << "\" with linespoints,\\" << endl;
            out_file << "\"" << label((data_file_name[2].str()).c_str()) << "\" using 2:2 title \""
                     << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
                     << STAT_label[STATL_FUNCTION] << "\" with linespoints" << endl; */

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

          out_file << "plot [0:" << 1. - complement << "] [0:" << 1. - complement << "] ";
          if (histo->variance > 0.) {
            out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using 5:6 title \""
                     << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_CONCENTRATION] << " "
                     << STAT_label[STATL_CURVE] << "\" with linespoints,\\" << endl;
          }
          out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using 3:4 title \""
                   << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_CONCENTRATION] << " "
                   << STAT_label[STATL_CURVE] << "\" with linespoints,\\" << endl;
          out_file << "\"" << label((data_file_name[1].str()).c_str())
                   << "\" using 3:3 notitle with lines" << endl;

          out_file << "unset grid\n" << "set xtics autofreq\n" << "set ytics autofreq" << endl;
        }

        if (i == 1) {
          out_file << "\nset terminal x11" << endl;
        }

        out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
    }

    delete [] data_file_name;

    if (!status) {
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }
  }

  else {
    status = Distribution::plot_write(error , prefix , 0 , NULL , title);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet DiscreteParametricModel.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool DiscreteParametricModel::plot_write(StatError &error , const char *prefix ,
                                         const char *title) const

{
  return plot_write(error , prefix , title , frequency_distribution);
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'une loi discrete parametrique et d'une loi discrete empirique.
 *
 *  argument : pointeur sur une loi discrete empirique.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* DiscreteParametricModel::get_plotable(const DiscreteDistributionData *histo) const

{
  MultiPlotSet *plot_set;


  if (histo) {
    register int i , j;
    int nb_plot_set , plot_nb_value;
    double scale , *pcumul;
    ostringstream legend , title;


    nb_plot_set = 1;
    if (variance > 0.) {
      nb_plot_set += 3;
    }
    if (histo->variance == 0.) {
      nb_plot_set--;
    }
 
    plot_set = new MultiPlotSet(nb_plot_set);
    MultiPlotSet &plot = *plot_set;

    title.str("");
    title << STAT_label[STATL_DISTRIBUTION];
    if (histo) {
      title << " " << STAT_label[STATL_FIT];
    }
    plot.title = title.str();

    plot.border = "15 lw 0";

    // 1ere vue : ajustement

    plot_nb_value = plot_nb_value_computation(histo);

    plot[0].xrange = Range(0 , MAX(plot_nb_value - 1 , 1));
    plot[0].yrange = Range(0 , ceil(MAX(histo->max ,
                                    max * histo->nb_element) * YSCALE));

    if (plot_nb_value - 1 < TIC_THRESHOLD) {
      plot[0].xtics = 1;
    }

    plot[0].resize(2);

    plot[0][0].legend = STAT_label[STATL_FREQUENCY_DISTRIBUTION];

    plot[0][0].style = "impulses";

    histo->plotable_frequency_write(plot[0][0]);

    legend.str("");
    legend << STAT_label[STATL_DISTRIBUTION];
    plot_title_print(legend);
    plot[0][1].legend = legend.str();

    plot[0][1].style = "linespoints";

    plotable_mass_write(plot[0][1] , histo->nb_element);

    if (variance > 0.) {

      // calcul de la fonction de repartition de la loi empirique

      scale = histo->nb_element / (1. - complement);
      pcumul = histo->cumul_computation(scale);

      // 2eme vue : fonctions de repartition

      plot[1].xrange = Range(0 , plot_nb_value - 1);
      plot[1].yrange = Range(0. , 1. - complement);

      if (plot_nb_value - 1 < TIC_THRESHOLD) {
        plot[1].xtics = 1;
      }

      plot[1].resize(histo->variance > 0. ? 2 : 1);

      if (histo->variance > 0.) {
        legend.str("");
        legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
               << STAT_label[STATL_FUNCTION];
        plot[1][0].legend = legend.str();

        plot[1][0].style = "linespoints";

        histo->plotable_cumul_write(plot[1][0] , pcumul);
        i = 1;
      }

      else {
        i = 0;
      }

      legend.str("");
      legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION];
      plot[1][i].legend = legend.str();

      plot[1][i].style = "linespoints";

      plotable_cumul_write(plot[1][i]);

      // 3eme vue : mise en correspondance des fonctions de repartition en prenant
      // comme reference celle de la loi

      if (histo->variance > 0.) {
        title.str("");
        title << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
              << STAT_label[STATL_FUNCTION] << " " << STAT_label[STATL_MATCHING];
        plot[2].title = title.str();

        plot[2].xrange = Range(0. , 1. - complement);
        plot[2].yrange = Range(0. , 1. - complement);

        plot[2].grid = true;

        plot[2].xtics = 0.1;
        plot[2].ytics = 0.1;

        plot[2].resize(2);

/*        legend.str("");
        legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
               << STAT_label[STATL_FUNCTION];
               plot[2][0].legend = legend.str(); */

        plot[2][0].style = "linespoints";

        histo->plotable_cumul_matching_write(plot[2][0] , offset , nb_value , cumul , pcumul);

/*        legend.str("");
        legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
               << STAT_label[STATL_FUNCTION];
        plot[2][1].legend = legend.str();

        plot[2][1].style = "linespoints";

        plotable_cumul_matching_write(plot[2][1] , *this); */

        plot[2][1].style = "lines";

        plot[2][1].add_point(0. , 0.);
        plot[2][1].add_point(1. - complement , 1. - complement);
      }

      // courbes de concentration

      if (histo->variance == 0.) {
        i = 2;
      }
      else {
        i = 3;
      }

      plot[i].xrange = Range(0. , 1. - complement);
      plot[i].yrange = Range(0. , 1. - complement);

      plot[i].grid = true;

      plot[i].xtics = 0.1;
      plot[i].ytics = 0.1;

      plot[i].resize(histo->variance > 0. ? 3 : 2);

      if (histo->variance > 0.) {
        legend.str("");
        legend << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_CONCENTRATION] << " "
               << STAT_label[STATL_CURVE];
        plot[i][0].legend = legend.str();

        plot[i][0].style = "linespoints";

        histo->plotable_concentration_write(plot[i][0] , pcumul , scale);
        j = 1;
      }

      else {
        j = 0;
      }

      legend.str("");
      legend << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_CONCENTRATION] << " "
             << STAT_label[STATL_CURVE];
      plot[i][j].legend = legend.str();

      plot[i][j].style = "linespoints";

      plotable_concentration_write(plot[i][j]);
      j++;

      plot[i][j].style = "lines";

      plot[i][j].add_point(0. , 0.);
      plot[i][j].add_point(1. - complement , 1. - complement);
    }
  }

  else {
    StatError error;

    plot_set = Distribution::get_plotable_distributions(error , 0 , NULL);
  }

  return plot_set;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet DiscreteParametricModel.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* DiscreteParametricModel::get_plotable() const

{
  return get_plotable(frequency_distribution);
}


};  // namespace stat_tool
