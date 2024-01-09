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
 *       $Id: vectors.cpp 18024 2015-04-23 07:08:49Z guedon $
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
#include "markovian.h"
#include "vectors.h"
#include "stat_label.h"

// #include "quantile_computation.h"   probleme compilateur C++ Windows

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Vectors.
 *
 *--------------------------------------------------------------*/

Vectors::Vectors()

{
  nb_vector = 0;
  identifier = NULL;

  nb_variable = 0;

  type = NULL;
  min_value = NULL;
  max_value = NULL;
  min_interval = NULL;
  marginal_distribution = NULL;
  marginal_histogram = NULL;

  mean = NULL;
  covariance = NULL;

  int_vector = NULL;
  real_vector = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Initialisation d'un objet Vectors.
 *
 *  arguments : nombre de vecteurs, identificateurs des vecteurs,
 *              nombre de variables, type de chaque variable, flag initialisation.
 *
 *--------------------------------------------------------------*/

void Vectors::init(int inb_vector , int *iidentifier , int inb_variable ,
                   int *itype , bool init_flag)

{
  register int i , j;


  nb_vector = inb_vector;

  identifier = new int[nb_vector];

  if (iidentifier) {
    for (i = 0;i < nb_vector;i++) {
      identifier[i] = iidentifier[i];
    }
  }

  else {
    for (i = 0;i < nb_vector;i++) {
      identifier[i] = i + 1;
    }
  }

  nb_variable = inb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  min_interval = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = itype[i];
    min_value[i] = 0.;
    max_value[i] = 0.;
    min_interval[i] = 0.;
    marginal_distribution[i] = NULL;
    marginal_histogram[i] = NULL;
  }

  mean = new double[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    mean[i] = D_INF;
  }

  covariance = new double*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    covariance[i] = new double[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      covariance[i][j] = D_DEFAULT;
    }
  }

  int_vector = new int*[nb_vector];
  real_vector = new double*[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    int_vector[i] = new int[nb_variable];
    real_vector[i] = new double[nb_variable];

    if (init_flag) {
      for (j = 0;j < nb_variable;j++) {
        int_vector[i][j] = 0;
        real_vector[i][j] = 0.;
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Vectors.
 *
 *  arguments : nombre de vecteurs, identificateurs des vecteurs,
 *              nombre de variables, vecteurs entiers.
 *
 *--------------------------------------------------------------*/

Vectors::Vectors(int inb_vector , int *iidentifier , int inb_variable ,
                 int **iint_vector)

{
  register int i , j;
  int *itype;


  itype = new int[inb_variable];
  for (i = 0;i < inb_variable;i++) {
    itype[i] = INT_VALUE;
  }

  init(inb_vector , iidentifier , inb_variable , itype , false);
  delete [] itype;

  for (i = 0;i < nb_vector;i++) {
    for (j = 0;j < nb_variable;j++) {
      int_vector[i][j] = iint_vector[i][j];
    }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value_computation(i);
    max_value_computation(i);

    build_marginal_frequency_distribution(i);
    min_interval_computation(i);
  }

  covariance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Vectors.
 *
 *  arguments : nombre de vecteurs, identificateurs des vecteurs,
 *              nombre de variables, vecteurs reels.
 *
 *--------------------------------------------------------------*/

Vectors::Vectors(int inb_vector , int *iidentifier , int inb_variable ,
                 double **ireal_vector)

{
  register int i , j;
  int *itype;


  itype = new int[inb_variable];
  for (i = 0;i < inb_variable;i++) {
    itype[i] = REAL_VALUE;
  }

  init(inb_vector , iidentifier , inb_variable , itype , false);
  delete [] itype;

  for (i = 0;i < nb_vector;i++) {
    for (j = 0;j < nb_variable;j++) {
      real_vector[i][j] = ireal_vector[i][j];
    }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value_computation(i);
    max_value_computation(i);

    build_marginal_histogram(i);
    min_interval_computation(i);

    mean_computation(i);
    variance_computation(i);
  }

  covariance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Vectors.
 *
 *  arguments : nombre de vecteurs, identificateurs des vecteurs,
 *              nombre de variables, type de chaque variable,
 *              variables entieres, variables reelles, indexation variable.
 *
 *--------------------------------------------------------------*/

Vectors::Vectors(int inb_vector , int *iidentifier , int inb_variable , int *itype ,
                 int **iint_vector , double **ireal_vector , bool variable_index)

{
  register int i , j , k , m;


  init(inb_vector , iidentifier , inb_variable , itype , false);

  switch (variable_index) {

  case false : {
    i = 0;
    j = 0;
    for (k = 0;k < nb_variable;k++) {
      if (type[k] != REAL_VALUE) {
        for (m = 0;m < nb_vector;m++) {
          int_vector[m][k] = iint_vector[m][i];
        }
        i++;
      }

      else {
        for (m = 0;m < nb_vector;m++) {
          real_vector[m][k] = ireal_vector[m][j];
        }
        j++;
      }
    }
    break;
  }

  case true : {
    for (i = 0;i < nb_variable;i++) {
      if (type[i] != REAL_VALUE) {
        for (j = 0;j < nb_vector;j++) {
          int_vector[j][i] = iint_vector[j][i];
        }
      }

      else {
        for (j = 0;j < nb_vector;j++) {
          real_vector[j][i] = ireal_vector[j][i];
        }
      }
    }
    break;
  }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value_computation(i);
    max_value_computation(i);

    if (type[i] != REAL_VALUE) {
      build_marginal_frequency_distribution(i);
    }

    else {
      build_marginal_histogram(i);

      mean_computation(i);
      variance_computation(i);
    }

    min_interval_computation(i);
  }

  covariance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Vectors.
 *
 *  arguments : reference sur un objet Vectors, indice de la variable,
 *              type de la variable selectionnee.
 *
 *--------------------------------------------------------------*/

Vectors::Vectors(const Vectors &vec , int variable , int itype)

{
  register int i , j;


  nb_vector = vec.nb_vector;

  identifier = new int[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    identifier[i] = vec.identifier[i];
  }

  nb_variable = vec.nb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  min_interval = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    if (i != variable) {
      type[i] = vec.type[i];
      min_value[i] = vec.min_value[i];
      max_value[i] = vec.max_value[i];
      min_interval[i] = vec.min_interval[i];

      if (vec.marginal_distribution[i]) {
        marginal_distribution[i] = new FrequencyDistribution(*(vec.marginal_distribution[i]));
      }
      else {
        marginal_distribution[i] = NULL;
      }

      if (vec.marginal_histogram[i]) {
        marginal_histogram[i] = new Histogram(*(vec.marginal_histogram[i]));
      }
      else {
        marginal_histogram[i] = NULL;
      }
    }

    else {
      type[i] = itype;
      min_value[i] = 0.;
      max_value[i] = 0.;
      min_interval[i] = 0.;
      marginal_distribution[i] = NULL;
      marginal_histogram[i] = NULL;
    }
  }

  mean = new double[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    if (i != variable) {
      mean[i] = vec.mean[i];
    }
    else {
      mean[i] = D_INF;
    }
  }

  covariance = new double*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    covariance[i] = new double[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      if ((i != variable) && (j != variable)) {
        covariance[i][j] = vec.covariance[i][j];
      }
      else {
        covariance[i][j] = D_DEFAULT;
      }
    }
  }

  int_vector = new int*[nb_vector];
  real_vector = new double*[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    int_vector[i] = new int[nb_variable];
    real_vector[i] = new double[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      if (j != variable) {
        int_vector[i][j] = vec.int_vector[i][j];
        real_vector[i][j] = vec.real_vector[i][j];
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Vectors.
 *
 *  arguments : reference sur un objet Vectors, nombre de vecteurs,
 *              index des vecteurs selectionnes.
 *
 *--------------------------------------------------------------*/

Vectors::Vectors(const Vectors &vec , int inb_vector , int *index)

{
  register int i , j;
  int *pivector , *civector;


  nb_vector = inb_vector;

  identifier = new int[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    identifier[i] = vec.identifier[index[i]];
  }

  nb_variable = vec.nb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  min_interval = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = vec.type[i];
    marginal_distribution[i] = NULL;
    marginal_histogram[i] = NULL;
  }

  mean = new double[nb_variable];

  covariance = new double*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    covariance[i] = new double[nb_variable];
  }

  int_vector = new int*[nb_vector];
  real_vector = new double*[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    int_vector[i] = new int[nb_variable];
    real_vector[i] = new double[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      int_vector[i][j] = vec.int_vector[index[i]][j];
      real_vector[i][j] = vec.real_vector[index[i]][j];
    }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value_computation(i);
    max_value_computation(i);

    if (type[i] != REAL_VALUE) {
      build_marginal_frequency_distribution(i);
    }

    else  {
      build_marginal_histogram(i , vec.marginal_histogram[i]->step);

      mean_computation(i);
      variance_computation(i);
    }

    min_interval_computation(i);
  }

  covariance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Vectors.
 *
 *  argument : reference sur un objet Vectors.
 *
 *--------------------------------------------------------------*/

void Vectors::copy(const Vectors &vec)

{
  register int i , j;


  nb_vector = vec.nb_vector;

  identifier = new int[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    identifier[i] = vec.identifier[i];
  }

  nb_variable = vec.nb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  min_interval = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = vec.type[i];
    min_value[i] = vec.min_value[i];
    max_value[i] = vec.max_value[i];
    min_interval[i] = vec.min_interval[i];

    if (vec.marginal_distribution[i]) {
      marginal_distribution[i] = new FrequencyDistribution(*(vec.marginal_distribution[i]));
    }
    else {
      marginal_distribution[i] = NULL;
    }

    if (vec.marginal_histogram[i]) {
      marginal_histogram[i] = new Histogram(*(vec.marginal_histogram[i]));
    }
    else {
      marginal_histogram[i] = NULL;
    }
  }

  mean = new double[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    mean[i] = vec.mean[i];
  }

  covariance = new double*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    covariance[i] = new double[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      covariance[i][j] = vec.covariance[i][j];
    }
  }

  int_vector = new int*[nb_vector];
  real_vector = new double*[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    int_vector[i] = new int[nb_variable];
    real_vector[i] = new double[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      int_vector[i][j] = vec.int_vector[i][j];
      real_vector[i][j] = vec.real_vector[i][j];
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Vectors avec ajout d'une variable d'etat.
 *
 *  arguments : reference sur un objet Vectors.
 *
 *--------------------------------------------------------------*/

void Vectors::add_state_variable(const Vectors &vec)

{
  register int i , j;


  nb_vector = vec.nb_vector;

  identifier = new int[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    identifier[i] = vec.identifier[i];
  }

  nb_variable = vec.nb_variable + 1;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  min_interval = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  type[0] = STATE;
  min_value[0] = 0.;
  max_value[0] = 0.;
  min_interval[0] = 0.; 
  marginal_distribution[0] = NULL;
  marginal_histogram[0] = NULL;

  for (i = 0;i < vec.nb_variable;i++) {
    type[i + 1] = vec.type[i];
    min_value[i + 1] = vec.min_value[i];
    max_value[i + 1] = vec.max_value[i];
    min_interval[i + 1] = vec.min_interval[i];

    if (vec.marginal_distribution[i]) {
      marginal_distribution[i + 1] = new FrequencyDistribution(*(vec.marginal_distribution[i]));
    }
    else {
      marginal_distribution[i + 1] = NULL;
    }

    if (vec.marginal_histogram[i]) {
      marginal_histogram[i + 1] = new Histogram(*(vec.marginal_histogram[i]));
    }
    else {
      marginal_histogram[i + 1] = NULL;
    }
  }

  mean = new double[nb_variable];
  for (i = 0;i < vec.nb_variable;i++) {
    mean[i + 1] = vec.mean[i];
  }

  covariance = new double*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    covariance[i] = new double[nb_variable];
  }

  for (i = 0;i < vec.nb_variable;i++) {
    for (j = 0;j < vec.nb_variable;j++) {
      covariance[i + 1][j + 1] = vec.covariance[i][j];
    }
  }

  int_vector = new int*[nb_vector];
  real_vector = new double*[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    int_vector[i] = new int[nb_variable];
    real_vector[i] = new double[nb_variable];
    for (j = 0;j < vec.nb_variable;j++) {
      int_vector[i][j + 1] = vec.int_vector[i][j];
      real_vector[i][j + 1] = vec.real_vector[i][j];
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe Vectors.
 *
 *  arguments : reference sur un objet Vectors, type de transformation ('c' : copie,
 *              'a' : ajout d'une variable d'etat).
 *
 *--------------------------------------------------------------*/

Vectors::Vectors(const Vectors &vec , char transform)

{
  switch (transform) {
  case 'a' :
    Vectors::add_state_variable(vec);
    break;
  default :
    Vectors::copy(vec);
    break;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Vectors.
 *
 *--------------------------------------------------------------*/

void Vectors::remove()

{
  register int i;


  delete [] identifier;

  delete [] type;

  delete [] min_value;
  delete [] max_value;
  delete [] min_interval;

  if (marginal_distribution) {
    for (i = 0;i < nb_variable;i++) {
      delete marginal_distribution[i];
    }
    delete [] marginal_distribution;
  }

  if (marginal_histogram) {
    for (i = 0;i < nb_variable;i++) {
      delete marginal_histogram[i];
    }
    delete [] marginal_histogram;
  }

  delete [] mean;

  if (covariance) {
    for (i = 0;i < nb_variable;i++) {
      delete [] covariance[i];
    }
    delete [] covariance;
  }

  if (int_vector) {
    for (i = 0;i < nb_vector;i++) {
      delete [] int_vector[i];
    }
    delete [] int_vector;
  }

  if (real_vector) {
    for (i = 0;i < nb_vector;i++) {
      delete [] real_vector[i];
    }
    delete [] real_vector;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Vectors.
 *
 *--------------------------------------------------------------*/

Vectors::~Vectors()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Vectors.
 *
 *  argument : reference sur un objet Vectors.
 *
 *--------------------------------------------------------------*/

Vectors& Vectors::operator=(const Vectors &vec)

{
  if (&vec != this) {
    remove();
    copy(vec);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Construction des valeurs reelles pour les variables entieres.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Vectors::build_real_vector(int variable)

{
  register int i , j;


  for (i = 0;i < nb_variable;i++) {
    if (((variable == I_DEFAULT) || (variable == i)) && (type[i] == INT_VALUE)) {
      for (j = 0;j < nb_vector;j++) {
        real_vector[j][i] = int_vector[j][i];
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Extraction de la loi marginale empirique pour une variable entiere.
 *
 *  arguments : reference sur un objet StatError, variable.
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData* Vectors::extract(StatError &error , int variable) const

{
  bool status = true;
  DiscreteDistributionData *histo;


  histo = NULL;
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
      error.update(STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION]);
    }
  }

  if (status) {
    histo = new DiscreteDistributionData(*marginal_distribution[variable]);
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Verification de la liste complete des identificateurs.
 *
 *  arguments : reference sur un objet StatError, nombre d'individus,
 *              identificateurs des individus.
 *
 *--------------------------------------------------------------*/

bool identifier_checking(StatError &error , int nb_individual , int *identifier)

{
  bool status = true , *selected_identifier;
  register int i;
  int max_identifier;


  max_identifier = identifier[0];
  for (i = 1;i < nb_individual;i++) {
    if (identifier[i] > max_identifier) {
      max_identifier = identifier[i];
    }
  }

  selected_identifier = new bool[max_identifier + 1];
  for (i = 0;i <= max_identifier;i++) {
    selected_identifier[i] = false;
  }

  for (i = 0;i < nb_individual;i++) {
    if (selected_identifier[identifier[i]]) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_IDENTIFIER] << " " << identifier[i] << " "
                    << STAT_error[STATR_ALREADY_USED];
      error.update((error_message.str()).c_str());
    }
    else {
      selected_identifier[identifier[i]] = true;
    }
  }

  delete [] selected_identifier;

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Verification d'un objet Vectors.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

bool Vectors::check(StatError &error)

{
  bool status = true , lstatus;


  error.init();

  if (nb_variable > VECTOR_NB_VARIABLE) {
    status = false;
    error.update(STAT_error[STATR_NB_VARIABLE]);
  }

  lstatus = identifier_checking(error , nb_vector , identifier);
  if (!lstatus) {
    status = false;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Fusion d'objets Vectors.
 *
 *  argument : reference sur un objet StatError, nombre d'objets Vectors,
 *             pointeurs sur les objets Vectors.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::merge(StatError &error , int nb_sample , const Vectors **ivec) const

{
  bool status = true;
  register int i , j , k , m;
  int inb_vector , cumul_nb_vector , *iidentifier;
  const FrequencyDistribution **phisto;
  Vectors *vec;
  const Vectors **pvec;


  vec = NULL;
  error.init();

  for (i = 0;i < nb_sample;i++) {
    if (ivec[i]->nb_variable != nb_variable) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << STAT_error[STATR_NB_VARIABLE];
      error.correction_update((error_message.str()).c_str() , nb_variable);
    }

    else {
      for (j = 0;j < nb_variable;j++) {
        if (ivec[i]->type[j] != type[j]) {
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
    pvec = new const Vectors*[nb_sample];

    pvec[0] = this;
    for (i = 1;i < nb_sample;i++) {
      pvec[i] = ivec[i - 1];
    }

    // calcul du nombre de vecteurs

    inb_vector = 0;
    for (i = 0;i < nb_sample;i++) {
      inb_vector += pvec[i]->nb_vector;
    }

    // comparaison des identificateurs des vecteurs

    iidentifier = new int[inb_vector];

    cumul_nb_vector = 0;
    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < pvec[j]->nb_vector;k++) {
        iidentifier[i] = pvec[j]->identifier[k];

        for (m = 0;m < cumul_nb_vector;m++) {
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
      cumul_nb_vector += pvec[j]->nb_vector;
    }

    vec = new Vectors(inb_vector , iidentifier , nb_variable , type);
    delete [] iidentifier;

    // copie des vecteurs

    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < pvec[j]->nb_vector;k++) {
        for (m = 0;m < pvec[j]->nb_variable;m++) {
          vec->int_vector[i][m] = pvec[j]->int_vector[k][m];
          vec->real_vector[i][m] = pvec[j]->real_vector[k][m];
        }
        i++;
      }
    }

    phisto = new const FrequencyDistribution*[nb_sample];

    for (i = 0;i < vec->nb_variable;i++) {
      vec->min_value[i] = pvec[0]->min_value[i];
      vec->max_value[i] = pvec[0]->max_value[i];
      for (j = 1;j < nb_sample;j++) {
        if (pvec[j]->min_value[i] < vec->min_value[i]) {
          vec->min_value[i] = pvec[j]->min_value[i];
        }
        if (pvec[j]->max_value[i] > vec->max_value[i]) {
          vec->max_value[i] = pvec[j]->max_value[i];
        }
      }

      vec->min_interval[i] = pvec[0]->min_interval[i];
      for (j = 1;j < nb_sample;j++) {
        if (pvec[j]->min_interval[i] < vec->min_interval[i]) {
          vec->min_interval[i] = pvec[j]->min_interval[i];
        }
      }

      for (j = 0;j < nb_sample;j++) {
        if (pvec[j]->marginal_distribution[i]) {
          phisto[j] = pvec[j]->marginal_distribution[i];
        }
        else {
          break;
        }
      }

      if (j == nb_sample) {
        vec->marginal_distribution[i] = new FrequencyDistribution(nb_sample , phisto);

        vec->mean[i] = vec->marginal_distribution[i]->mean;
        vec->covariance[i][i] = vec->marginal_distribution[i]->variance;
      }

      else {
        vec->build_marginal_histogram(i);

        vec->mean_computation(i);
        vec->variance_computation(i);
      }
    }

    vec->covariance_computation();

    delete [] phisto;
    delete [] pvec;
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Translation des valeurs d'une variable.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              parametre de translation.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::shift(StatError &error , int variable , int shift_param) const

{
  bool status = true;
  register int i;
  Vectors *vec;


  vec = NULL;
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
    vec = new Vectors(*this , variable , type[variable]);

    switch (vec->type[variable]) {

    // translation des valeurs entieres

    case INT_VALUE : {
      for (i = 0;i < vec->nb_vector;i++) {
        vec->int_vector[i][variable] = int_vector[i][variable] + shift_param;
      }
      break;
    }

    // translation des valeurs reelles

    case REAL_VALUE : {
      for (i = 0;i < vec->nb_vector;i++) {
        vec->real_vector[i][variable] = real_vector[i][variable] + shift_param;
      }
      break;
    }
    }

    vec->min_value[variable] = min_value[variable] + shift_param;
    vec->max_value[variable] = max_value[variable] + shift_param;
    vec->min_interval[variable] = min_interval[variable];

    if ((vec->type[variable] == INT_VALUE) && (vec->min_value[variable] >= 0) &&
        (vec->max_value[variable] <= MARGINAL_DISTRIBUTION_MAX_VALUE)) {
      if (marginal_distribution[variable]) {
        vec->marginal_distribution[variable] = new FrequencyDistribution(*marginal_distribution[variable] , 's' , shift_param);

        vec->mean[variable] = vec->marginal_distribution[variable]->mean;
      }

      else {
        vec->build_marginal_frequency_distribution(variable);
      }
    }

    else {
      vec->build_marginal_histogram(variable);

      vec->mean[variable] = mean[variable] + shift_param;
    }

    for (i = 0;i < vec->nb_variable;i++) {
      vec->covariance[i][variable] = covariance[i][variable];
      vec->covariance[variable][i] = covariance[variable][i];
    }
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Translation des valeurs d'une variable reelle.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              parametre de translation.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::shift(StatError &error , int variable , double shift_param) const

{
  bool status = true;
  register int i;
  Vectors *vec;


  vec = NULL;
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
    vec = new Vectors(*this , variable , type[variable]);

    // translation des valeurs reelles

    for (i = 0;i < vec->nb_vector;i++) {
      vec->real_vector[i][variable] = real_vector[i][variable] + shift_param;
    }

    vec->min_value[variable] = min_value[variable] + shift_param;
    vec->max_value[variable] = max_value[variable] + shift_param;
    vec->min_interval[variable] = min_interval[variable];

    vec->build_marginal_histogram(variable);

    vec->mean[variable] = mean[variable] + shift_param;

    for (i = 0;i < vec->nb_variable;i++) {
      vec->covariance[i][variable] = covariance[i][variable];
      vec->covariance[variable][i] = covariance[variable][i];
    }
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Seuillage des valeurs d'une variable.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              seuil, mode (ABOVE/BELOW).
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::thresholding(StatError &error , int variable , int threshold , int mode) const

{
  bool status = true;
  register int i;
  Vectors *vec;


  vec = NULL;
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
    vec = new Vectors(*this , variable , type[variable]);

    switch (vec->type[variable]) {

    // seuillage des valeurs entieres

    case INT_VALUE : {
      switch (mode) {

      case ABOVE : {
        for (i = 0;i < vec->nb_vector;i++) {
          if (int_vector[i][variable] > threshold) {
            vec->int_vector[i][variable] = threshold;
          }
          else {
            vec->int_vector[i][variable] = int_vector[i][variable];
          }
        }
        break;
      }

      case BELOW : {
        for (i = 0;i < vec->nb_vector;i++) {
          if (int_vector[i][variable] < threshold) {
            vec->int_vector[i][variable] = threshold;
          }
          else {
            vec->int_vector[i][variable] = int_vector[i][variable];
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
        for (i = 0;i < vec->nb_vector;i++) {
          if (real_vector[i][variable] > threshold) {
            vec->real_vector[i][variable] = threshold;
          }
          else {
            vec->real_vector[i][variable] = real_vector[i][variable];
          }
        }
        break;
      }

      case BELOW : {
        for (i = 0;i < vec->nb_vector;i++) {
          if (real_vector[i][variable] < threshold) {
            vec->real_vector[i][variable] = threshold;
          }
          else {
            vec->real_vector[i][variable] = real_vector[i][variable];
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
      vec->min_value[variable] = min_value[variable];
      vec->max_value[variable] = threshold;
      break;
    case BELOW :
      vec->min_value[variable] = threshold;
      vec->max_value[variable] = max_value[variable];
      break;
    }

    if ((vec->type[variable] == INT_VALUE) && (vec->min_value[variable] >= 0) &&
        (vec->max_value[variable] <= MARGINAL_DISTRIBUTION_MAX_VALUE)) {
      vec->build_marginal_frequency_distribution(variable);
    }
    else {
      vec->build_marginal_histogram(variable);
      vec->mean_computation(variable);
    }

    vec->min_interval_computation(variable);

    vec->covariance_computation(variable);
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Seuillage des valeurs d'une variable reelle.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              seuil, mode (ABOVE/BELOW).
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::thresholding(StatError &error , int variable , double threshold , int mode) const

{
  bool status = true;
  register int i;
  Vectors *vec;


  vec = NULL;
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
    vec = new Vectors(*this , variable , type[variable]);

    // seuillage des valeurs reelles

    switch (mode) {

    case ABOVE : {
      for (i = 0;i < vec->nb_vector;i++) {
        if (real_vector[i][variable] > threshold) {
          vec->real_vector[i][variable] = threshold;
        }
        else {
          vec->real_vector[i][variable] = real_vector[i][variable];
        }
      }

      vec->min_value[variable] = min_value[variable];
      vec->max_value[variable] = threshold;
      break;
    }

    case BELOW : {
      for (i = 0;i < vec->nb_vector;i++) {
        if (real_vector[i][variable] < threshold) {
          vec->real_vector[i][variable] = threshold;
        }
        else {
          vec->real_vector[i][variable] = real_vector[i][variable];
        }
      }

      vec->min_value[variable] = threshold;
      vec->max_value[variable] = max_value[variable];
      break;
    }
    }

    vec->build_marginal_histogram(variable);

    vec->min_interval_computation(variable);

    vec->mean_computation(variable);
    vec->covariance_computation(variable);
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une variable.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              pas de regroupement, mode (FLOOR/ROUND/CEIL).
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::cluster(StatError &error , int variable , int step , int mode) const

{
  bool status = true;
  register int i;
  Vectors *vec;


  vec = NULL;
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

  if (step < 1) {
    status = false;
    error.update(STAT_error[STATR_CLUSTERING_STEP]);
  }

  if (status) {
    vec = new Vectors(*this , variable , type[variable]);

    switch (vec->type[variable]) {

    // regroupement des valeurs entieres

    case INT_VALUE : {
      switch (mode) {

      case FLOOR : {
        for (i = 0;i < vec->nb_vector;i++) {
          vec->int_vector[i][variable] = int_vector[i][variable] / step;
        }

        vec->min_value[variable] = (int)min_value[variable] / step;
        vec->max_value[variable] = (int)max_value[variable] / step;
//        vec->min_value[variable] = floor(min_value[variable] / step);
//        vec->max_value[variable] = floor(max_value[variable] / step);
        break;
      }

      case ROUND : {
        for (i = 0;i < vec->nb_vector;i++) {
          vec->int_vector[i][variable] = (int_vector[i][variable] + step / 2) / step;
//          vec->int_vector[i][variable] = (int)::round((double)int_vector[i][variable] / (double)step);
        }

        vec->min_value[variable] = ((int)min_value[variable] + step / 2) / step;
        vec->max_value[variable] = ((int)max_value[variable] + step / 2) / step;
//        vec->min_value[variable] = ::round(min_value[variable] / step);
//        vec->max_value[variable] = ::round(max_value[variable] / step);
        break;
      }

      case CEIL : {
        for (i = 0;i < vec->nb_vector;i++) {
          vec->int_vector[i][variable] = (int_vector[i][variable] + step - 1) / step;
//          vec->int_vector[i][variable] = (int)ceil((double)int_vector[i][variable] / (double)step);
        }

        vec->min_value[variable] = ((int)min_value[variable] + step - 1) / step;
        vec->max_value[variable] = ((int)max_value[variable] + step - 1) / step;
//        vec->min_value[variable] = ceil(min_value[variable] / step);
//        vec->max_value[variable] = ceil(max_value[variable] / step);
        break;
      }
      }

      if (marginal_distribution[variable]) {
        vec->marginal_distribution[variable] = new FrequencyDistribution(*marginal_distribution[variable] , 'c' , step , mode);

        vec->mean[variable] = vec->marginal_distribution[variable]->mean;
        vec->covariance[variable][variable] = vec->marginal_distribution[variable]->variance;
      }

      else {
        vec->build_marginal_frequency_distribution(variable);
      }

      vec->min_interval_computation(variable);
      break;
    }

    // regroupement des valeurs reelles

    case REAL_VALUE : {
      for (i = 0;i < vec->nb_vector;i++) {
        vec->real_vector[i][variable] = real_vector[i][variable] / step;
      }

      vec->min_value[variable] = min_value[variable] / step;
      vec->max_value[variable] = max_value[variable] / step;

      vec->build_marginal_histogram(variable , marginal_histogram[variable]->step / step);

      vec->min_interval_computation(variable);

      vec->mean[variable] = mean[variable] / step;
      vec->covariance[variable][variable] = covariance[variable][variable] / (step * step);
      break;
    }
    }

    vec->covariance_computation(variable);
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Transcodage des symboles d'une variable entiere.
 *
 *  arguments : reference sur un objet Vectors, indice de la variable,
 *              plus petit et plus grand symboles, table de transcodage des symboles.
 *
 *--------------------------------------------------------------*/

void Vectors::transcode(const Vectors &vec , int variable , int min_symbol ,
                        int max_symbol , int *symbol)

{
  register int i;


  for (i = 0;i < nb_vector;i++) {
    int_vector[i][variable] = symbol[vec.int_vector[i][variable] -
                                     (int)vec.min_value[variable]] + min_symbol;
  }

  min_value[variable] = min_symbol;
  max_value[variable] = max_symbol;

  build_marginal_frequency_distribution(variable);
  min_interval_computation(variable);

  covariance_computation(variable);
}


/*--------------------------------------------------------------*
 *
 *  Transcodage des symboles d'une variable entiere.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              table de transcodage des symboles.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::transcode(StatError &error , int variable , int *symbol) const

{
  bool status = true , *presence;
  register int i;
  int min_symbol , max_symbol;
  Vectors *vec;


  vec = NULL;
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

      vec = new Vectors(*this , variable , type[variable]);
      vec->transcode(*this , variable , min_symbol , max_symbol , symbol);
    }
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une variable entiere.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              nombre de classes, bornes pour regrouper les valeurs.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::cluster(StatError &error , int variable ,
                          int nb_class , int *ilimit) const

{
  bool status = true;
  register int i , j , k;
  int *int_limit , *symbol;
  double *real_limit;
  Vectors *vec;


  vec = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((nb_class < 2) || (nb_class >= (int)(max_value[variable] - min_value[variable]) + 1)) {
      status = false;
      error.update(STAT_error[STATR_NB_CLASS]);
    }
  }

  if (status) {
    if (type[variable] != REAL_VALUE) {
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

        vec = new Vectors(*this , variable , type[variable]);
        vec->transcode(*this , variable , 0 , nb_class - 1 , symbol);

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
        vec = new Vectors(*this , variable , INT_VALUE);
        vec->cluster(*this , variable , nb_class , real_limit);
      }

      delete [] real_limit;
    }
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une variable reelle.
 *
 *  arguments : reference sur un objet Vectors, indice de la variable,
 *              nombre de classes, bornes pour regrouper les valeurs.
 *
 *--------------------------------------------------------------*/

void Vectors::cluster(const Vectors &vec , int variable , int nb_class , double *limit)

{
  register int i , j;


  for (i = 0;i < nb_vector;i++) {
    for (j = 0;j < nb_class;j++) {
      if (vec.real_vector[i][variable] < limit[j + 1]) {
        int_vector[i][variable] = j;
        break;
      }
    }
  }

  min_value_computation(variable);
  max_value_computation(variable);

  build_marginal_frequency_distribution(variable);
  min_interval_computation(variable);

  covariance_computation(variable);
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une variable reelle.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              nombre de classes, bornes pour regrouper les valeurs.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::cluster(StatError &error , int variable ,
                          int nb_class , double *ilimit) const

{
  bool status = true;
  register int i;
  double *limit;
  Vectors *vec;


  vec = NULL;
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
      vec = new Vectors(*this , variable , INT_VALUE);
      vec->cluster(*this , variable , nb_class , limit);
    }

    delete [] limit;
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Changement d'unite d'une variable.
 *
 *  arguments : reference sur un objet StatError, variable, facteur d'echelle.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::scaling(StatError &error , int variable , int scaling_coeff) const

{
  bool status = true;
  register int i;
  Vectors *vec;


  vec = NULL;
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
    vec = new Vectors(*this , variable , type[variable]);

    switch (vec->type[variable]) {

    // mise a l'echelle des valeurs entieres

    case INT_VALUE : {
      for (i = 0;i < vec->nb_vector;i++) {
        vec->int_vector[i][variable] = int_vector[i][variable] * scaling_coeff;
      }
      break;
    }

    // mise a l'echelle des valeurs reelles

    case REAL_VALUE : {
      for (i = 0;i < vec->nb_vector;i++) {
        vec->real_vector[i][variable] = real_vector[i][variable] * scaling_coeff;
      }
      break;
    }
    }

    vec->min_value[variable] = min_value[variable] * scaling_coeff;
    vec->max_value[variable] = max_value[variable] * scaling_coeff;
    vec->min_interval[variable] = min_interval[variable] * scaling_coeff;

    vec->build_marginal_frequency_distribution(variable);

    vec->covariance_computation(variable);
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Changement d'unite d'une variable.
 *
 *  arguments : reference sur un objet StatError, variable, facteur d'echelle.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::scaling(StatError &error , int variable , double scaling_coeff) const

{
  bool status = true;
  register int i;
  Vectors *vec;


  vec = NULL;
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
    vec = new Vectors(*this , variable , REAL_VALUE);

    switch (type[variable]) {

    // mise a l'echelle des valeurs entieres

    case INT_VALUE : {
      for (i = 0;i < vec->nb_vector;i++) {
        vec->real_vector[i][variable] = int_vector[i][variable] * scaling_coeff;
      }
      break;
    }

    // mise a l'echelle des valeurs reelles

    case REAL_VALUE : {
      for (i = 0;i < vec->nb_vector;i++) {
        vec->real_vector[i][variable] = real_vector[i][variable] * scaling_coeff;
      }
      break;
    }
    }

    vec->min_value[variable] = min_value[variable] * scaling_coeff;
    vec->max_value[variable] = max_value[variable] * scaling_coeff;
    vec->min_interval[variable] = min_interval[variable] * scaling_coeff;

    vec->build_marginal_histogram(variable);

    vec->mean[variable] = mean[variable] * scaling_coeff;
    vec->covariance_computation(variable);
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Arrondi des valeurs d'une variable reelle.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              mode (FLOOR/ROUND/CEIL).
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::round(StatError &error , int variable , int mode) const

{
  bool status = true;
  register int i , j;
  int *itype;
  Vectors *vec;


  vec = NULL;
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
        itype[i] = INT_VALUE;
      }
    }

    else {
      for (i = 0;i < nb_variable;i++) {
        itype[i] = type[i];
      }
      itype[variable] = INT_VALUE;
    }

    vec = new Vectors(nb_vector , identifier , nb_variable , itype);
    delete [] itype;

    for (i = 0;i < vec->nb_vector;i++) {
      for (j = 0;j < vec->nb_variable;j++) {

        // copie des valeurs entieres

        if (type[j] != REAL_VALUE) {
          vec->int_vector[i][j] = int_vector[i][j];
        }

        else {

          // arrondi des valeurs reelles

          if ((variable == I_DEFAULT) || (variable == j)) {
            switch (mode) {
            case FLOOR :
              vec->int_vector[i][j] = (int)floor(real_vector[i][j]);
              break;
            case ROUND :
              vec->int_vector[i][j] = (int)::round(real_vector[i][j]);
              break;
            case CEIL :
              vec->int_vector[i][j] = (int)ceil(real_vector[i][j]);
              break;
            }
          }

          // copie des valeurs reelles

          else {
            vec->real_vector[i][j] = real_vector[i][j];
          }
        }
      }
    }

    for (i = 0;i < vec->nb_variable;i++) {
      if ((variable == I_DEFAULT) || (variable == i)) {
        switch (mode) {
        case FLOOR :
          vec->min_value[i] = floor(min_value[i]);
          vec->max_value[i] = floor(max_value[i]);
          break;
        case ROUND :
          vec->min_value[i] = ::round(min_value[i]);
          vec->max_value[i] = ::round(max_value[i]);
          break;
        case CEIL :
          vec->min_value[i] = ceil(min_value[i]);
          vec->max_value[i] = ceil(max_value[i]);
          break;
        }

        vec->build_marginal_frequency_distribution(i);
        vec->min_interval_computation(i);
      }

      else {
        vec->min_value[i] = min_value[i];
        vec->max_value[i] = max_value[i];
        vec->min_interval[i] = min_interval[i];

        if (marginal_distribution[i]) {
          vec->marginal_distribution[i] = new FrequencyDistribution(*marginal_distribution[i]);
        }
        if (marginal_histogram[i]) {
          vec->marginal_histogram[i] = new Histogram(*marginal_histogram[i]);
        }

        vec->mean[i] = mean[i];
        for (j = 0;j < vec->nb_variable;j++) {
          vec->covariance[i][j] = covariance[i][j];
        }
      }
    }

    vec->covariance_computation(variable);
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Selection de vecteurs sur les valeurs prises par une variable.
 *
 *  arguments : reference sur un objet StatError, stream, indice de la variable,
 *              bornes sur les valeurs, flag pour conserver ou rejeter
 *              les vecteurs selectionnes.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::value_select(StatError &error , ostream &os , int variable ,
                               int imin_value , int imax_value , bool keep) const

{
  bool status = true;
  register int i;
  int inb_vector , *index , *iidentifier;
  Vectors *vec;


  vec = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((imin_value > max_value[variable]) || (imin_value > imax_value)) {
      status = false;
      error.update(STAT_error[STATR_MIN_VALUE]);
    }
    if ((imax_value < min_value[variable]) || (imax_value < imin_value)) {
      status = false;
      error.update(STAT_error[STATR_MAX_VALUE]);
    }
  }

  if (status) {

    // selection des vecteurs

    iidentifier = new int[nb_vector];
    index = new int[nb_vector];
    inb_vector = 0;

    if (type[variable] != REAL_VALUE) {
      for (i = 0;i < nb_vector;i++) {
        if ((int_vector[i][variable] >= imin_value) && (int_vector[i][variable] <= imax_value)) {
          if (keep) {
            iidentifier[inb_vector] = identifier[i];
            index[inb_vector++] = i;
          }
        }

        else if (!keep) {
          iidentifier[inb_vector] = identifier[i];
          index[inb_vector++] = i;
        }
      }
    }

    else {
      for (i = 0;i < nb_vector;i++) {
        if ((real_vector[i][variable] >= imin_value) && (real_vector[i][variable] <= imax_value)) {
          if (keep) {
            iidentifier[inb_vector] = identifier[i];
            index[inb_vector++] = i;
          }
        }

        else if (!keep) {
          iidentifier[inb_vector] = identifier[i];
          index[inb_vector++] = i;
        }
      }
    }

    if (inb_vector == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    // copie des vecteurs

    if (status) {
      if (inb_vector <= DISPLAY_NB_INDIVIDUAL) {
        os << "\n" << STAT_label[inb_vector == 1 ? STATL_VECTOR : STATL_VECTORS] << ": ";
        for (i = 0;i < inb_vector;i++) {
          os << iidentifier[i] << ", ";
        }
        os << endl;
      }

      vec = new Vectors(*this , inb_vector , index);
    }

    delete [] iidentifier;
    delete [] index;
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Selection de vecteurs sur les valeurs prises par une variable reelle.
 *
 *  arguments : reference sur un objet StatError, stream, indice de la variable,
 *              bornes sur les valeurs, flag pour conserver ou rejeter
 *              les vecteurs selectionnes.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::value_select(StatError &error , ostream &os , int variable ,
                               double imin_value , double imax_value , bool keep) const

{
  bool status = true;
  register int i;
  int inb_vector , *index , *iidentifier;
  Vectors *vec;


  vec = NULL;
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
    
    if ((imin_value > max_value[variable]) || (imin_value > imax_value)) {
      status = false;
      error.update(STAT_error[STATR_MIN_VALUE]);
    }
    if ((imax_value < min_value[variable]) || (imax_value < imin_value)) {
      status = false;
      error.update(STAT_error[STATR_MAX_VALUE]);
    }
  }

  if (status) {

    // selection des vecteurs

    iidentifier = new int[nb_vector];
    index = new int[nb_vector];
    inb_vector = 0;

    for (i = 0;i < nb_vector;i++) {
      if ((real_vector[i][variable] >= imin_value) && (real_vector[i][variable] <= imax_value)) {
        if (keep) {
          iidentifier[inb_vector] = identifier[i];
          index[inb_vector++] = i;
        }
      }

      else if (!keep) {
        iidentifier[inb_vector] = identifier[i];
        index[inb_vector++] = i;
      }
    }

    if (inb_vector == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    // copie des vecteurs

    if (status) {
      if (inb_vector <= DISPLAY_NB_INDIVIDUAL) {
        os << "\n" << STAT_label[inb_vector == 1 ? STATL_VECTOR : STATL_VECTORS] << ": ";
        for (i = 0;i < inb_vector;i++) {
          os << iidentifier[i] << ", ";
        }
        os << endl;
      }

      vec = new Vectors(*this , inb_vector , index);
    }

    delete [] iidentifier;
    delete [] index;
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Verification d'une liste d'identificateurs.
 *
 *  arguments : reference sur un objet StatError, nombre d'individus,
 *              identificateurs des individus, nombre d'individus selectionnes,
 *              identificateurs des individus selectionnes, label du type de donnees.
 *
 *--------------------------------------------------------------*/

bool selected_identifier_checking(StatError &error , int nb_individual , int *identifier ,
                                  int nb_selected_individual , int *selected_identifier ,
                                  const char *data_label)

{
  bool status = true , *selected_individual;
  register int i , j;
  int max_identifier;


  max_identifier = selected_identifier[0];
  for (i = 1;i < nb_selected_individual;i++) {
    if (selected_identifier[i] > max_identifier) {
      max_identifier = selected_identifier[i];
    }
  }

  selected_individual = new bool[max_identifier + 1];
  for (i = 0;i <= max_identifier;i++) {
    selected_individual[i] = false;
  }

  for (i = 0;i < nb_selected_individual;i++) {
    for (j = 0;j < nb_individual;j++) {
      if (selected_identifier[i] == identifier[j]) {
        break;
      }
    }

    if (j == nb_individual) {
      status = false;
      ostringstream error_message;
      error_message << selected_identifier[i] << ": " << STAT_error[STATR_BAD] << " "
                    << data_label << " " << STAT_label[STATL_IDENTIFIER];
      error.update((error_message.str()).c_str());
    }

    else if (selected_individual[selected_identifier[i]]) {
      status = false;
      ostringstream error_message;
      error_message << data_label << " " << selected_identifier[i] << " "
                    << STAT_error[STATR_ALREADY_SELECTED];
      error.update((error_message.str()).c_str());
    }
    else {
      selected_individual[selected_identifier[i]] = true;
    }
  }

  delete [] selected_individual;

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Selection d'individus par l'identificateur.
 *
 *  arguments : nombre d'individus, identificateurs des individus,
 *              nombre d'individus selectionnes, identificateurs des individus selectionnes,
 *              flag pour conserver ou rejeter les individus selectionnes.
 *
 *--------------------------------------------------------------*/

int* identifier_select(int nb_individual , int *identifier , int nb_selected_individual ,
                       int *selected_identifier , bool keep)

{
  register int i , j , k;
  int *index;


  switch (keep) {

  case false : {
    index = new int[nb_individual - nb_selected_individual];

    i = 0;
    for (j = 0;j < nb_individual;j++) {
      for (k = 0;k < nb_selected_individual;k++) {
        if (selected_identifier[k] == identifier[j]) {
          break;
        }
      }

      if (k == nb_selected_individual) {
        index[i++] = j;
      }
    }
    break;
  }

  case true : {
    index = new int[nb_selected_individual];

    for (i = 0;i < nb_selected_individual;i++) {
      for (j = 0;j < nb_individual;j++) {
        if (selected_identifier[i] == identifier[j]) {
          index[i] = j;
          break;
        }
      }
    }
    break;
  }
  }

  return index;
}


/*--------------------------------------------------------------*
 *
 *  Selection de vecteurs par l'identificateur.
 *
 *  arguments : reference sur un objet StatError, nombre de vecteurs,
 *              identificateurs des vecteurs, flag pour conserver ou rejeter
 *              les vecteurs selectionnees.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::select_individual(StatError &error , int inb_vector ,
                                    int *iidentifier , bool keep) const

{
  bool status = true;
  int *index;
  Vectors *vec;


  vec = NULL;
  error.init();

  if ((inb_vector < 1) || (inb_vector > (keep ? nb_vector : nb_vector - 1))) {
    status = false;
    error.update(STAT_error[STATR_NB_VECTOR]);
  }

  else {
    status = selected_identifier_checking(error , nb_vector , identifier , inb_vector ,
                                          iidentifier , STAT_label[STATL_VECTOR]);
  }

  if (status) {
    index = identifier_select(nb_vector , identifier , inb_vector , iidentifier , keep);

    vec = new Vectors(*this , (keep ? inb_vector : nb_vector - inb_vector) , index);

    delete [] index;
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Selection de variables.
 *
 *  arguments : reference sur un objet Vectors, indices des variables selectionnees.
 *
 *--------------------------------------------------------------*/

void Vectors::select_variable(const Vectors &vec , int *variable)

{
  register int i , j;


  for (i = 0;i < nb_vector;i++) {
    for (j = 0;j < nb_variable;j++) {
      int_vector[i][j] = vec.int_vector[i][variable[j]];
      real_vector[i][j] = vec.real_vector[i][variable[j]];
    }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value[i] = vec.min_value[variable[i]];
    max_value[i] = vec.max_value[variable[i]];
    min_interval[i] = vec.min_interval[variable[i]];

    if (vec.marginal_distribution[variable[i]]) {
      marginal_distribution[i] = new FrequencyDistribution(*(vec.marginal_distribution[variable[i]]));
    }
    if (vec.marginal_histogram[variable[i]]) {
      marginal_histogram[i] = new Histogram(*(vec.marginal_histogram[variable[i]]));
    }

    mean[i] = vec.mean[variable[i]];
    for (j = 0;j < nb_variable;j++) {
      covariance[i][j] = vec.covariance[variable[i]][variable[j]];
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Selection de 2 variables.
 *
 *  arguments : variable explicative, variable expliquee.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::select_variable(int explanatory_variable , int response_variable) const

{
  int variable[2] , itype[2];
  Vectors *vec;


  variable[0] = explanatory_variable;
  variable[1] = response_variable;

  itype[0] = type[explanatory_variable];
  itype[1] = type[response_variable];

  vec = new Vectors(nb_vector , identifier , 2 , itype);

  vec->select_variable(*this , variable);

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Selection de variables.
 *
 *  arguments : nombre de variables, nombre de variables selectionnees,
 *              indices des variables selectionnees, flag pour conserver ou
 *              rejeter les variables selectionnees.
 *
 *--------------------------------------------------------------*/

int* select_variable(int nb_variable , int nb_selected_variable ,
                     int *selected_variable , bool keep)

{
  register int i , j , k;
  int *variable;


  for (i = 0;i < nb_selected_variable;i++) {
    selected_variable[i]--;
  }

  switch (keep) {

  case false : {
    variable = new int[nb_variable - nb_selected_variable];

    i = 0;
    for (j = 0;j < nb_variable;j++) {
      for (k = 0;k < nb_selected_variable;k++) {
        if (selected_variable[k] == j) {
          break;
        }
      }

      if (k == nb_selected_variable) {
        variable[i++] = j;
      }
    }
    break;
  }

  case true : {
    variable = new int[nb_selected_variable];

    for (i = 0;i < nb_selected_variable;i++) {
      variable[i] = selected_variable[i];
    }
    break;
  }
  }

  return variable;
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

Vectors* Vectors::select_variable(StatError &error , int inb_variable ,
                                  int *ivariable , bool keep) const

{
  bool status = true , *selected_variable;
  register int i;
  int bnb_variable , *variable , *itype;
  Vectors *vec;


  vec = NULL;
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
    variable = stat_tool::select_variable(nb_variable , inb_variable , ivariable , keep);

    bnb_variable = (keep ? inb_variable : nb_variable - inb_variable);

    itype = new int[bnb_variable];
    for (i = 0;i < bnb_variable;i++) {
      itype[i] = type[variable[i]];
    }

    vec = new Vectors(nb_vector , identifier , bnb_variable , itype);
    vec->select_variable(*this , variable);

    delete [] variable;
    delete [] itype;
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Suppression de la 1ere variable.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::remove_variable_1() const

{
  register int i;
  int *variable , *itype;
  Vectors *vec;


  variable = new int[nb_variable - 1];
  itype = new int[nb_variable - 1];
  for (i = 0;i < nb_variable - 1;i++) {
    variable[i] = i + 1;
    itype[i] = type[i + 1];
  }

  vec = new Vectors(nb_vector , identifier , nb_variable - 1 , itype);
  vec->select_variable(*this , variable);

  delete [] variable;
  delete [] itype;

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Concatenation des variables d'objets Vectors.
 *
 *  arguments : reference sur un objet StatError, nombre d'objets Vectors,
 *              pointeurs sur les objets Vectors, echantillon de reference pour les identificateurs.
 *
 *--------------------------------------------------------------*/

Vectors* Vectors::merge_variable(StatError &error , int nb_sample ,
                                 const Vectors **ivec , int ref_sample) const

{
  bool status = true;
  register int i , j , k;
  int inb_variable , *iidentifier , *itype;
  Vectors *vec;
  const Vectors **pvec;


  vec = NULL;
  error.init();

  for (i = 0;i < nb_sample;i++) {
    if (ivec[i]->nb_vector != nb_vector) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << STAT_error[STATR_NB_VECTOR];
      error.update((error_message.str()).c_str());
    }
  }

  if ((ref_sample != I_DEFAULT) && ((ref_sample < 1) || (ref_sample > nb_sample + 1))) {
    status = false;
    error.update(STAT_error[STATR_SAMPLE_INDEX]);
  }

  if (status) {
    nb_sample++;
    pvec = new const Vectors*[nb_sample];

    pvec[0] = this;
    inb_variable = nb_variable;
    for (i = 1;i < nb_sample;i++) {
      pvec[i] = ivec[i - 1];
      inb_variable += ivec[i - 1]->nb_variable;
    }

    // comparaison des identificateurs des vecteurs

    if (ref_sample == I_DEFAULT) {
      for (i = 0;i < nb_vector;i++) {
        for (j = 1;j < nb_sample;j++) {
          if (pvec[j]->identifier[i] != pvec[0]->identifier[i]) {
            break;
          }
        }
        if (j < nb_sample) {
          break;
        }
      }

      if (i < nb_vector) {
        iidentifier = NULL;
      }
      else {
        iidentifier = pvec[0]->identifier;
      }
    }

    else {
      ref_sample--;
      iidentifier = pvec[ref_sample]->identifier;
    }

    itype = new int[inb_variable];
    inb_variable = 0;
    for (i = 0;i < nb_sample;i++) {
      for (j = 0;j < pvec[i]->nb_variable;j++) {
        itype[inb_variable] = pvec[i]->type[j];
        if ((inb_variable > 0) && (itype[inb_variable] == STATE)) {
          itype[inb_variable] = INT_VALUE;
        }
        inb_variable++;
      }
    }

    vec = new Vectors(nb_vector , iidentifier , inb_variable , itype);
    delete [] itype;

    // copie des vecteurs

    for (i = 0;i < nb_vector;i++) {
      inb_variable = 0;
      for (j = 0;j < nb_sample;j++) {
        for (k = 0;k < pvec[j]->nb_variable;k++) {
          vec->int_vector[i][inb_variable] = pvec[j]->int_vector[i][k];
          vec->real_vector[i][inb_variable++] = pvec[j]->real_vector[i][k];
        }
      }
    }

    inb_variable = 0;
    for (i = 0;i < nb_sample;i++) {
      for (j = 0;j < pvec[i]->nb_variable;j++) {
        vec->min_value[inb_variable] = pvec[i]->min_value[j];
        vec->max_value[inb_variable] = pvec[i]->max_value[j];
        vec->min_interval[inb_variable] = pvec[i]->min_interval[j];

        if (pvec[i]->marginal_distribution[j]) {
          vec->marginal_distribution[inb_variable] = new FrequencyDistribution(*(pvec[i]->marginal_distribution[j]));
        }
        if (pvec[i]->marginal_histogram[j]) {
          vec->marginal_histogram[inb_variable] = new Histogram(*(pvec[i]->marginal_histogram[j]));
        }

        vec->mean[inb_variable] = pvec[i]->mean[j];
        vec->covariance[inb_variable][inb_variable] = pvec[i]->covariance[j][j];

        inb_variable++;
      }
    }

    vec->covariance_computation();

    delete [] pvec;
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Vectors a partir d'un fichier.
 *  Format : chaque ligne represente un individu.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

Vectors* vectors_ascii_read(StatError &error , const char *path)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status , lstatus;
  register int i , j;
  int line , read_line , initial_nb_line , nb_variable = 0 , nb_vector , index , *type;
  long int_value;
  double real_value;
  Vectors *vec;
  ifstream in_file(path);


  vec = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {

    // 1ere passe : analyse de la ligne definissant le nombre de variables

    status = true;
    line = 0;
    read_line = 0;

    while (buffer.readLine(in_file , false)) {
      line++;

      position = buffer.first('#');
      if (position != RW_NPOS) {
        buffer.remove(position);
      }
      i = 0;

      RWCTokenizer next(buffer);

      while (!((token = next()).isNull())) {
        switch (i) {

        // test nombre de variables

        case 0 : {
          lstatus = locale.stringToNum(token , &int_value);
          if (lstatus) {
            if ((int_value < 1) || (int_value > VECTOR_NB_VARIABLE)) {
              lstatus = false;
            }
            else {
              nb_variable = int_value;
            }
          }

          if (!lstatus) {
            status = false;
            error.update(STAT_parsing[STATP_NB_VARIABLE] , line , i + 1);
          }
          break;
        }

        // test mot cle VARIABLE(S)

        case 1 : {
          if (token != STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEY_WORD] ,
                                    STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] , line , i + 1);
          }
          break;
        }
        }

        i++;
      }

      if (i > 0) {
        if (i != 2) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }

        read_line++;
        break;
      }
    }

    if (read_line < 1) {
      status = false;
      error.update(STAT_parsing[STATP_FORMAT]);
    }

    // analyse des lignes definissant le type de chaque variable

    if (status) {
      type = new int[nb_variable];
      for (i = 0;i < nb_variable;i++) {
        type[i] = I_DEFAULT;
      }

      read_line = 0;

      while (buffer.readLine(in_file , false)) {
        line++;

        position = buffer.first('#');
        if (position != RW_NPOS) {
          buffer.remove(position);
        }
        i = 0;

        RWCTokenizer next(buffer);

        while (!((token = next()).isNull())) {
          switch (i) {

          // test mot cle VARIABLE

          case 0 : {
            if (token != STAT_word[STATW_VARIABLE]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_VARIABLE] , line , i + 1);
            }
            break;
          }

          // test index de la variable

          case 1 : {
            lstatus = locale.stringToNum(token , &int_value);
            if ((lstatus) && (int_value != read_line + 1)) {
              lstatus = false;
            }

            if (!lstatus) {
              status = false;
              error.correction_update(STAT_parsing[STATP_VARIABLE_INDEX] , read_line + 1 , line , i + 1);
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

          // test mot cle correspondant au type de la variable

          case 3 : {
            for (j = INT_VALUE;j <= STATE;j++) {
              if (token == STAT_variable_word[j]) {
                type[read_line] = j;
                break;
              }
            }

            if (j == REAL_VALUE + 1) {
              status = false;
              error.update(STAT_parsing[STATP_KEY_WORD] , line , i + 1);
            }
            break;
          }
          }

          i++;
        }

        if (i > 0) {
          if (i != 4) {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line);
          }

          read_line++;
          if (read_line == nb_variable) {
            break;
          }
        }
      }

      if (read_line < nb_variable) {
        status = false;
        error.update(STAT_parsing[STATP_FORMAT]);
      }

      initial_nb_line = line;
    }

    if (status) {

      // analyse des lignes et recherche du nombre de vecteurs

      line = 0;
      nb_vector = 0;

      while (buffer.readLine(in_file , false)) {
        line++;

#       ifdef DEBUG
        cout << line << "  " << buffer << endl;
#       endif

        position = buffer.first('#');
        if (position != RW_NPOS) {
          buffer.remove(position);
        }
        i = 0;

        RWCTokenizer next(buffer);

        while (!((token = next()).isNull())) {
          if (i <= nb_variable) {
            if (type[i] != REAL_VALUE) {
              lstatus = locale.stringToNum(token , &int_value);
            }
            else {
              lstatus = locale.stringToNum(token , &real_value);
            }

            if (!lstatus) {
              status = false;
              error.update(STAT_parsing[STATP_DATA_TYPE] , line , i + 1);
            }
          }

          i++;
        }

        // test nombre de valeurs par ligne constant

        if (i > 0) {
          if (i != nb_variable) {
            status = false;
            error.correction_update(STAT_parsing[STATP_NB_TOKEN] , nb_variable , line);
          }

          nb_vector++;
        }
      }

      if (nb_vector == 0) {
        status = false;
        error.update(STAT_parsing[STATP_EMPTY_SAMPLE]);
      }
    }

    // 2eme passe : copie des vecteurs

    if (status) {
//      in_file.close();
//      in_file.open(path , ios::in);

      in_file.clear();
      in_file.seekg(0 , ios::beg);

      line = 0;

      do {
        buffer.readLine(in_file , false);
        line++;

#       ifdef DEBUG
        cout << line << "  " << buffer << endl;
#       endif

      }
      while (line < initial_nb_line);

      vec = new Vectors(nb_vector , NULL , nb_variable , type);

      index = 0;

      while (buffer.readLine(in_file , false)) {
        position = buffer.first('#');
        if (position != RW_NPOS) {
          buffer.remove(position);
        }
        i = 0;

        RWCTokenizer next(buffer);

        while (!((token = next()).isNull())) {
          if (type[i] != REAL_VALUE) {
            locale.stringToNum(token , &int_value);
            vec->int_vector[index][i++] = int_value;
          }
          else {
            locale.stringToNum(token , &real_value);
            vec->real_vector[index][i++] = real_value;
          }
        }

        if (i > 0) {
          index++;
        }
      }

      for (i = 0;i < vec->nb_variable;i++) {
        vec->min_value_computation(i);
        vec->max_value_computation(i);

        vec->build_marginal_frequency_distribution(i);
        vec->min_interval_computation(i);
      }

      vec->covariance_computation();
    }
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet Vectors.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Vectors::line_write(ostream &os) const

{
  os << nb_vector << " " << STAT_label[nb_vector == 1 ? STATL_VECTOR : STATL_VECTORS] << "   "
     << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES];

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Vectors.
 *
 *  arguments : stream, flag niveau de detail, flag commentaire.
 *
 *--------------------------------------------------------------*/

ostream& Vectors::ascii_write(ostream &os , bool exhaustive , bool comment_flag) const

{
  register int i , j;
  int buff , width[2];
  long old_adjust;
  double **correlation;
  Test *test;

  old_adjust = os.setf(ios::right , ios::adjustfield);

  if (comment_flag) {
    os << "# ";
  }
  os << nb_vector << " " << STAT_label[nb_vector == 1 ? STATL_VECTOR : STATL_VECTORS] << endl;

  os << "\n" << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;

  for (i = 0;i < nb_variable;i++) {
    os << "\n" << STAT_word[STATW_VARIABLE] << " " << i + 1 << " : "
       << STAT_variable_word[type[i]];

    os << "   ";
    if (comment_flag) {
      os << "# ";
    }
    os << "("  << STAT_label[STATL_MIN_VALUE] << ": " << min_value[i] << ", "
       << STAT_label[STATL_MAX_VALUE] << ": " << max_value[i] << ")" << endl;

    if (marginal_distribution[i]) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";

      marginal_distribution[i]->ascii_characteristic_print(os , exhaustive , comment_flag);

      if ((marginal_distribution[i]->nb_value <= ASCII_NB_VALUE) || (exhaustive)) {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << "   | " << STAT_label[STATL_FREQUENCY] << endl;
        marginal_distribution[i]->ascii_print(os , comment_flag);
      }
    }

    else {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_SAMPLE_SIZE] << ": " << nb_vector << endl;

      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_MEAN] << ": " << mean[i] << "   "
         << STAT_label[STATL_VARIANCE] << ": " << covariance[i][i] << "   "
         << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(covariance[i][i]) << endl;

      if ((exhaustive) && (covariance[i][i] > 0.)) {
        if (comment_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << skewness_computation(i) << "   "
           << STAT_label[STATL_KURTOSIS_COEFF] << ": " << kurtosis_computation(i) << endl;
      }

      if (exhaustive) {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM] << endl;

        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_VALUE] << " | " << STAT_label[STATL_FREQUENCY] << endl;
        marginal_histogram[i]->ascii_print(os , comment_flag);
      }
    }

/* #   ifdef DEBUG
    if (comment_flag) {
      os << "# ";
    }

    if (type[i] == INT_VALUE) {
      int *int_value;

      int_value = new int[nb_vector];
      for (j = 0;j < nb_vector;j++) {
        int_value[j] = int_vector[j][i];
      }
      os << "quartile: " << quantile_computation(nb_vector , int_value , 0.25 , false) << " | "
         << "median: " << quantile_computation(nb_vector , int_value , 0.5 , false) << " | "
         << "quartile: " << quantile_computation(nb_vector , int_value , 0.75 , false) << endl;

      os << "quartile: " << quantile_computation(nb_vector , int_value , 0.25 , true) << " | "
         << "median: " << quantile_computation(nb_vector , int_value , 0.5 , true) << " | "
         << "quartile: " << quantile_computation(nb_vector , int_value , 0.75 , true) << endl;

      delete [] int_value;
    }

    else {
      double *real_value;

      real_value = new double[nb_vector];
      for (j = 0;j < nb_vector;j++) {
        real_value[j] = real_vector[j][i];
      }
      os << "quartile: " << quantile_computation(nb_vector , real_value , 0.25 , false) << " | "
         << "median: " << quantile_computation(nb_vector , real_value , 0.5 , false) << " | "
         << "quartile: " << quantile_computation(nb_vector , real_value , 0.75 , false) << endl;

      delete [] real_value;
    }
#   endif */

  }

  width[0] = column_width(nb_variable);

  if (exhaustive) {

    // calcul des largeurs des colonnes

    width[1] = 0;
    for (i = 0;i < nb_variable;i++) {
      buff = column_width(nb_variable , covariance[i]);
      if (buff > width[1]) {
        width[1] = buff;
      }
    }
    width[1] += ASCII_SPACE;

    // ecriture de la matrice de variance-covariance

    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_VARIANCE_COVARIANCE_MATRIX] << endl;

    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << setw(width[0] + width[1]) << 1;
    for (i = 1;i < nb_variable;i++) {
      os << setw(width[1]) << i + 1;
    }
    for (i = 0;i < nb_variable;i++) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << setw(width[0]) << i + 1;
      for (j = 0;j < nb_variable;j++) {
        os << setw(width[1]) << covariance[i][j];
      }
    }
    os << endl;
  }

  correlation = correlation_computation();

  // calcul des largeurs des colonnes

  width[1] = 0;
  for (i = 0;i < nb_variable;i++) {
    buff = column_width(nb_variable , correlation[i]);
    if (buff > width[1]) {
      width[1] = buff;
    }
  }
  width[1] += ASCII_SPACE;

  // ecriture de la matrice des coefficients de correlation

  os << "\n";
  if (comment_flag) {
    os << "# ";
  }
  os << STAT_label[STATL_CORRELATION_MATRIX] << endl;

  os << "\n";
  if (comment_flag) {
    os << "# ";
  }
  os << setw(width[0] + width[1]) << 1;
  for (i = 1;i < nb_variable;i++) {
    os << setw(width[1]) << i + 1;
  }
  for (i = 0;i < nb_variable;i++) {
    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << setw(width[0]) << i + 1;
    for (j = 0;j < nb_variable;j++) {
      os << setw(width[1]) << correlation[i][j];
    }
  }
  os << endl;

  // test du caractere significatif des coefficients de correlation

  if (nb_vector > 2) {
#   ifdef DEBUG
    double value = fabs(correlation[0][1]) * sqrt(nb_vector - 2) / sqrt(1 - correlation[0][1] * correlation[0][1]);

    test = new Test(STUDENT , false , nb_vector - 2 , I_DEFAULT , value);
    test->t_critical_probability_computation();

    os << "\n" << test->value << "   " << STAT_label[STATL_CRITICAL_PROBABILITY]
       << ": " << test->critical_probability << endl;

    delete test;
#   endif

    test = new Test(STUDENT , false , nb_vector - 2 , I_DEFAULT , D_DEFAULT);

    for (i = 0;i < NB_CRITICAL_PROBABILITY;i++) {
      test->critical_probability = ref_critical_probability[i];
      test->t_value_computation();

      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_REFERENCE] << " " << STAT_label[STATL_T_VALUE] << ": "
         << test->value << "   " << STAT_label[STATL_REFERENCE] << " "
         << STAT_label[STATL_CRITICAL_PROBABILITY] << ": " << test->critical_probability << endl;

      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_LIMIT_CORRELATION_COEFF] << ": "
         << test->value / sqrt(test->value * test->value + nb_vector - 2) << endl;
    }

    delete test;
  }

  for (i = 0;i < nb_variable;i++) {
    delete [] correlation[i];
  }
  delete [] correlation;

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Vectors.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Vectors::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os, exhaustive, false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Vectors dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Vectors::ascii_write(StatError &error , const char *path ,
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
    ascii_write(out_file , exhaustive , false);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des vecteurs discrets.
 *
 *  arguments : stream, flag commentaire, probabilites a posteriori des affectactions optimales
 *              entropies des affectation.
 *
 *--------------------------------------------------------------*/

ostream& Vectors::ascii_print(ostream &os , bool comment_flag ,
                              double *posterior_probability , double *entropy) const

{
  register int i , j;
  int buff , width[2];
  long old_adjust;


  old_adjust = os.setf(ios::right , ios::adjustfield);

  width[0] = 0;
  for (i = 0;i < nb_variable;i++) {
    if (type[i] != REAL_VALUE) {
      buff = column_width((int)min_value[i] , (int)max_value[i]);
      if (buff > width[0]) {
        width[0] = buff;
      }
    }

    else {
/*      buff = column_width(min_value[i] , max_value[i]);
      if (buff > width[0]) {
        width[0] = buff;
      } */

      for (j = 0;j < nb_vector;j++) {
        buff = column_width(1 , real_vector[j] + i);
        if (buff > width[0]) {
          width[0] = buff;
        }
      }
    }
  }
  width[1] = width[0] + ASCII_SPACE;

  os << "\n";
  for (i = 0;i < nb_vector;i++) {
    if (type[0] != REAL_VALUE) {
      os << setw(width[0]) << int_vector[i][0];
    }
    else {
      os << setw(width[0]) << real_vector[i][0];
    }

    for (j = 1;j < nb_variable;j++) {
      if (type[j] != REAL_VALUE) {
        os << setw(width[1]) << int_vector[i][j];
      }
      else {
        os << setw(width[1]) << real_vector[i][j];
      }
    }

    os << "   ";
    if (comment_flag) {
      os << "# ";
    }
    os << "(" << identifier[i] << ")";

    if ((posterior_probability) && (entropy)) {
      os << "   " << STAT_label[STATL_POSTERIOR_ASSIGNMENT_PROBABILITY]
         << ": " << posterior_probability[i] << "   "
         << STAT_label[STATL_ASSIGNMENT_ENTROPY] << ": " << entropy[i];
    }
    os << endl;
  }

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Vectors.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Vectors::ascii_data_write(ostream &os , bool exhaustive) const

{
  ascii_write(os , exhaustive , false);
  ascii_print(os , false);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Vectors dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Vectors::ascii_data_write(StatError &error , const char *path , bool exhaustive) const

{
  bool status = false;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    ascii_write(out_file , exhaustive , true);
    ascii_print(out_file , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Vectors dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool Vectors::spreadsheet_write(StatError &error , const char *path) const

{
  bool status;
  register int i , j;
  double **correlation;
  Test *test;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    out_file << nb_vector << "\t" << STAT_label[nb_vector == 1 ? STATL_VECTOR : STATL_VECTORS] << endl;

    out_file << "\n" << nb_variable << "\t" << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;

    for (i = 0;i < nb_variable;i++) {
      out_file << "\n" << STAT_word[STATW_VARIABLE] << "\t" << i + 1<< "\t"
               << STAT_variable_word[type[i]];

      out_file << "\t\t" << STAT_label[STATL_MIN_VALUE] << "\t" << min_value[i]
               << "\t\t" << STAT_label[STATL_MAX_VALUE] << "\t" << max_value[i] << endl;

      if (marginal_distribution[i]) {
        out_file << "\n" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
        marginal_distribution[i]->spreadsheet_characteristic_print(out_file);

        out_file << "\n\t" << STAT_label[STATL_FREQUENCY] << endl;
        marginal_distribution[i]->spreadsheet_print(out_file);
      }

      else {
        out_file << "\n" << STAT_label[STATL_SAMPLE_SIZE] << "\t" << nb_vector << endl;

        out_file << STAT_label[STATL_MEAN] << "\t" << mean[i] << "\t\t"
                 << STAT_label[STATL_VARIANCE] << "\t" << covariance[i][i] << "\t\t"
                 << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(covariance[i][i]) << endl;

        if (covariance[i][i] > 0.) {
          out_file << STAT_label[STATL_SKEWNESS_COEFF] << "\t" << skewness_computation(i) << "\t\t"
                   << STAT_label[STATL_KURTOSIS_COEFF] << "\t" << kurtosis_computation(i) << endl;
        }

        out_file << "\n" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM] << endl;
        out_file << "\n" << STAT_label[STATL_VALUE] << "\t" << STAT_label[STATL_FREQUENCY] << endl;
        marginal_histogram[i]->spreadsheet_print(out_file);
      }
    }

    // ecriture de la matrice de variance-covariance

    out_file << "\n" << STAT_label[STATL_VARIANCE_COVARIANCE_MATRIX] << endl;

    out_file << "\n";
    for (i = 0;i < nb_variable;i++) {
      out_file << "\t" << i + 1;
    }
    for (i = 0;i < nb_variable;i++) {
      out_file << "\n" << i + 1;
      for (j = 0;j < nb_variable;j++) {
        out_file << "\t" << covariance[i][j];
      }
    }
    out_file << endl;

    // ecriture de la matrice des coefficients de correlation

    correlation = correlation_computation();

    out_file << "\n" << STAT_label[STATL_CORRELATION_MATRIX] << endl;

    out_file << "\n";
    for (i = 0;i < nb_variable;i++) {
      out_file << "\t" << i + 1;
    }
    for (i = 0;i < nb_variable;i++) {
      out_file << "\n" << i + 1;
      for (j = 0;j < nb_variable;j++) {
        out_file << "\t" << correlation[i][j];
      }
    }
    out_file << endl;

    // test du caractere significatif des coefficients de correlation

    if (nb_vector > 2) {
      test = new Test(STUDENT , false , nb_vector - 2 , I_DEFAULT , D_DEFAULT);

      for (i = 0;i < NB_CRITICAL_PROBABILITY;i++) {
        test->critical_probability = ref_critical_probability[i];
        test->t_value_computation();

        out_file << "\n" << STAT_label[STATL_REFERENCE] << " " << STAT_label[STATL_T_VALUE] << "\t"
                 << test->value << "\t" << STAT_label[STATL_REFERENCE] << " "
                 << STAT_label[STATL_CRITICAL_PROBABILITY] << "\t" << test->critical_probability << endl;

        out_file << STAT_label[STATL_LIMIT_CORRELATION_COEFF] << "\t"
                 << test->value / sqrt(test->value * test->value + nb_vector - 2) << endl;
      }

      delete test;
    }

    for (i = 0;i < nb_variable;i++) {
      delete [] correlation[i];
    }
    delete [] correlation;

    // ecriture des vecteurs

    for (i = 0;i < nb_vector;i++) {
      out_file << "\n";
      for (j = 0;j < nb_variable;j++) {
        if (type[j] != REAL_VALUE) {
          out_file << int_vector[i][j] << "\t";
        }
        else {
          out_file << real_vector[i][j] << "\t";
        }
      }
      out_file << "\t" << identifier[i];
    }
    out_file << endl;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des vecteurs au format Gnuplot.
 *
 *  arguments : path, residus reduits.
 *
 *--------------------------------------------------------------*/


bool Vectors::plot_print(const char *path , double *standard_residual) const

{
  bool status = false;
  register int i , j;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < nb_variable;j++) {
        if (type[j] != REAL_VALUE) {
          out_file << int_vector[i][j] << " ";
        }
        else {
          out_file << real_vector[i][j] << " ";
        }
      }

      if (standard_residual) {
        out_file << standard_residual[i];
      }
      out_file << endl;
    }
  }

  return status;
}



/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Vectors.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Vectors::plot_write(StatError &error , const char *prefix ,
                         const char *title) const

{
  bool status;
  register int i , j , k , m , n;
  int **frequency;
  ostringstream *data_file_name;


  error.init();

  // ecriture des fichiers de donnees

  data_file_name = new ostringstream[nb_variable + 1];

  data_file_name[0] << prefix << 0 << ".dat";
  status = plot_print((data_file_name[0].str()).c_str());

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  else {
    for (i = 0;i < nb_variable;i++) {
      data_file_name[i + 1] << prefix << i + 1 << ".dat";

      if (marginal_distribution[i]) {
        marginal_distribution[i]->plot_print((data_file_name[i + 1].str()).c_str());
      }
      else {
        marginal_histogram[i]->plot_print((data_file_name[i + 1].str()).c_str());
      }
    }

    // ecriture des fichiers de commandes et des fichier d'impression

    for (i = 0;i < 2;i++) {
      for (j = 0;j < nb_variable;j++) {
        ostringstream file_name[2];

        switch (i) {

        case 0 : {
          if (nb_variable == 1) {
            file_name[0] << prefix << ".plot";
          }
          else {
            file_name[0] << prefix << j + 1 << ".plot";
          }
          break;
        }

        case 1 : {
          if (nb_variable == 1) {
            file_name[0] << prefix << ".print";
          }
          else {
            file_name[0] << prefix << j + 1 << ".print";
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
            file_name[1] << label(prefix) << j + 1 << ".ps";
          }
          out_file << "set output \"" << file_name[1].str() << "\"\n\n";
        }

        out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                 << "set title";
        if (title) {
          out_file << " \"" << title << "\"";
        }
        out_file << "\n\n";

        if (marginal_distribution[j]) {
          if (marginal_distribution[j]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }
          if ((int)(marginal_distribution[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [0:" << MAX(marginal_distribution[j]->nb_value - 1 , 1) << "] [0:"
                   << (int)(marginal_distribution[j]->max * YSCALE) + 1 << "] \""
                   << label((data_file_name[j + 1].str()).c_str()) << "\" using 1 title \""
                   << STAT_label[STATL_VARIABLE] << " " << j + 1 << " "
                   << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                   << "\" with impulses" << endl;

          if (marginal_distribution[j]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
          if ((int)(marginal_distribution[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics autofreq" << endl;
          }
        }

        else {
          if ((int)(marginal_histogram[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [" << marginal_histogram[j]->min_value << ":"
                   << marginal_histogram[j]->max_value << "] [0:"
                   << (int)(marginal_histogram[j]->max * YSCALE) + 1 << "] \""
                   << label((data_file_name[j + 1].str()).c_str()) << "\" using 1:2 title \""
                   << STAT_label[STATL_VARIABLE] << " " << j + 1 << " "
                   << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM]
                   << "\" with histeps" << endl;

          if ((int)(marginal_histogram[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics autofreq" << endl;
          }
        }

        if ((i == 0) && (nb_variable > 1)) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        for (k = 0;k < nb_variable;k++) {
          if (k != j) {
            out_file << "set xlabel \"" << STAT_label[STATL_VARIABLE] << " " << j + 1 << "\"" << endl;
            out_file << "set ylabel \"" << STAT_label[STATL_VARIABLE] << " " << k + 1 << "\"" << endl;

            if (max_value[j] - min_value[j] < TIC_THRESHOLD) {
              out_file << "set xtics " << MIN(min_value[j] , 0) << ",1" << endl;
            }
            if (max_value[k] - min_value[k] < TIC_THRESHOLD) {
              out_file << "set ytics " << MIN(min_value[k] , 0) << ",1" << endl;
            }

            if (((marginal_distribution[j]) &&
                 (marginal_distribution[j]->nb_value <= PLOT_NB_VALUE)) &&
                ((marginal_distribution[k]) &&
                 (marginal_distribution[k]->nb_value <= PLOT_NB_VALUE))) {
              frequency = joint_frequency_computation(j , k);

              for (m = marginal_distribution[j]->offset;m < marginal_distribution[j]->nb_value;m++) {
                for (n = marginal_distribution[k]->offset;n < marginal_distribution[k]->nb_value;n++) {
                  if (frequency[m][n] > 0) {
                    out_file << "set label \"" << frequency[m][n] << "\" at " << m << ","
                             << n << endl;
                  }
                }
              }

              for (m = 0;m < marginal_distribution[j]->nb_value;m++) {
                delete [] frequency[m];
              }
              delete [] frequency;
            }

            if ((min_value[j] >= 0.) && (max_value[j] - min_value[j] > min_value[j] * PLOT_RANGE_RATIO)) {
              out_file << "plot [" << 0;
            }
            else {
              out_file << "plot [" << min_value[j];
            }
            out_file << ":" << MAX(max_value[j] , min_value[j] + 1) << "] [";
            if ((min_value[k] >= 0.) && (max_value[k] - min_value[k] > min_value[k] * PLOT_RANGE_RATIO)) {
              out_file << 0;
            }
            else {
              out_file << min_value[k];
            }
            out_file << ":" << MAX(max_value[k] , min_value[k] + 1) << "] \""
                     << label((data_file_name[0].str()).c_str()) << "\" using "
                     << j + 1 << ":" << k + 1 << " notitle with points" << endl;

            out_file << "unset label" << endl;

            out_file << "set xlabel" << endl;
            out_file << "set ylabel" << endl;

            if (max_value[j] - min_value[j] < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
            if (max_value[k] - min_value[k] < TIC_THRESHOLD) {
              out_file << "set ytics autofreq" << endl;
            }

            if ((i == 0) && (((j < nb_variable - 1) && (k < nb_variable - 1)) ||
                 ((j == nb_variable - 1) && (k < nb_variable - 2)))) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;
          }
        }

        if (i == 1) {
          out_file << "\nset terminal x11" << endl;
        }

        out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
    }
  }

  delete [] data_file_name;

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des vecteurs.
 *
 *  arguments : reference sur un objet SinglePlot, indices des 2 variables.
 *
 *--------------------------------------------------------------*/

void Vectors::plotable_write(SinglePlot &plot , int variable1 , int variable2) const

{
  register int i;


  if (type[variable1] != REAL_VALUE) {
    if (type[variable2] != REAL_VALUE) {
      for (i = 0;i < nb_vector;i++) {
        plot.add_point(int_vector[i][variable1] , int_vector[i][variable2]);
      }
    }
    else {
      for (i = 0;i < nb_vector;i++) {
        plot.add_point(int_vector[i][variable1] , real_vector[i][variable2]);
      }
    }
  }

  else {
    if (type[variable2] != REAL_VALUE) {
      for (i = 0;i < nb_vector;i++) {
        plot.add_point(real_vector[i][variable1] , int_vector[i][variable2]);
      }
    }
    else {
      for (i = 0;i < nb_vector;i++) {
        plot.add_point(real_vector[i][variable1] , real_vector[i][variable2]);
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des frequences associees aux vecteurs.
 *
 *  arguments : reference sur un objet SinglePlot, indice des 2 variables.
 *
 *--------------------------------------------------------------*/

void Vectors::plotable_frequency_write(SinglePlot &plot , int variable1 , int variable2) const

{
  register int i , j;
  int **frequency;
  ostringstream label;


  frequency = joint_frequency_computation(variable1 , variable2);

  for (i = marginal_distribution[variable1]->offset;i < marginal_distribution[variable1]->nb_value;i++) {
    for (j = marginal_distribution[variable2]->offset;j < marginal_distribution[variable2]->nb_value;j++) {
      if (frequency[i][j] > 0) {
        label.str("");
        label << frequency[i][j];

        plot.add_text(i , j , label.str());
      }
    }
  }

  for (i = 0;i < marginal_distribution[variable1]->nb_value;i++) {
    delete [] frequency[i];
  }
  delete [] frequency;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet Vectors.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Vectors::get_plotable() const

{
  register int i , j , k;
  int nb_plot_set;
  double xmin , ymin;
  ostringstream legend , label;
  MultiPlotSet *plot_set;


  nb_plot_set = nb_variable * nb_variable;

  plot_set = new MultiPlotSet(nb_plot_set , nb_variable);
  MultiPlotSet &plot = *plot_set;

  plot.border = "15 lw 0";

  i = 0;
  for (j = 0;j < nb_variable;j++) {
    plot.variable_nb_viewpoint[j] = 1;
    plot.variable[i] = j;

    if (marginal_distribution[j]) {

      // vue : loi marginale empirique

      plot[i].xrange = Range(0 , MAX(marginal_distribution[j]->nb_value - 1 , 1));
      plot[i].yrange = Range(0 , ceil(marginal_distribution[j]->max * YSCALE));

      if (marginal_distribution[j]->nb_value - 1 < TIC_THRESHOLD) {
        plot[i].xtics = 1;
      }
      if (ceil(marginal_distribution[j]->max * YSCALE) < TIC_THRESHOLD) {
        plot[i].ytics = 1;
      }

      plot[i].resize(1);

      legend.str("");
      legend << STAT_label[STATL_VARIABLE] << " " << j + 1 << " "
             << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      plot[i][0].legend = legend.str();

      plot[i][0].style = "impulses";

      marginal_distribution[j]->plotable_frequency_write(plot[i][0]);
    }

    else {

      // vue : histogramme marginal

      plot[i].xrange = Range(marginal_histogram[j]->min_value , marginal_histogram[j]->max_value);
      plot[i].yrange = Range(0 , ceil(marginal_histogram[j]->max * YSCALE));

      if (ceil(marginal_histogram[j]->max * YSCALE) < TIC_THRESHOLD) {
        plot[i].ytics = 1;
      }

      plot[i].resize(1);

      legend.str("");
      legend << STAT_label[STATL_VARIABLE] << " " << j + 1 << " "
             << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM];
      plot[i][0].legend = legend.str();

      plot[i][0].style = "histeps";

      marginal_histogram[j]->plotable_write(plot[i][0]);
    }
    i++;

    for (k = 0;k < nb_variable;k++) {
      if (k != j) {
        plot.variable[i] = j;

         // vue : loi jointe empirique de 2 variables

        if ((min_value[j] >= 0.) && (max_value[j] - min_value[j] > min_value[j] * PLOT_RANGE_RATIO)) {
          xmin = 0.;
        }
        else {
          xmin = min_value[j];
        }
        plot[i].xrange = Range(xmin , MAX(max_value[j] , min_value[j] + 1));

        if ((min_value[k] >= 0.) && (max_value[k] - min_value[k] > min_value[k] * PLOT_RANGE_RATIO)) {
          ymin = 0.;
        }
        else {
          ymin = min_value[k];
        }
        plot[i].yrange = Range(ymin , MAX(max_value[k] , min_value[k] + 1));

        if (max_value[j] - min_value[j] < TIC_THRESHOLD) {
          plot[i].xtics = 1;
        }
        if (max_value[k] - min_value[k] < TIC_THRESHOLD) {
          plot[i].ytics = 1;
        }

        label.str("");
        label << STAT_label[STATL_VARIABLE] << " " << j + 1;
        plot[i].xlabel = label.str();

        label.str("");
        label << STAT_label[STATL_VARIABLE] << " " << k + 1;
        plot[i].ylabel = label.str();

        if (((marginal_distribution[j]) &&
             (marginal_distribution[j]->nb_value <= PLOT_NB_VALUE)) &&
            ((marginal_distribution[k]) &&
             (marginal_distribution[k]->nb_value <= PLOT_NB_VALUE))) {
          plot[i].resize(2);
        }
        else {
          plot[i].resize(1);
        }

        plot[i][0].style = "points";

        plotable_write(plot[i][0] , j , k);

        if (((marginal_distribution[j]) &&
             (marginal_distribution[j]->nb_value <= PLOT_NB_VALUE)) &&
            ((marginal_distribution[k]) &&
             (marginal_distribution[k]->nb_value <= PLOT_NB_VALUE))) {
          plot[i][1].label = "true";

          plotable_frequency_write(plot[i][1] , j , k);
        }

        i++;
      }
    }
  }

  return plot_set;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur minimum prise par une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Vectors::min_value_computation(int variable)

{
  register int i;


  if (type[variable] != REAL_VALUE) {
    min_value[variable] = int_vector[0][variable];

    for (i = 1;i < nb_vector;i++) {
      if (int_vector[i][variable] < min_value[variable]) {
        min_value[variable] = int_vector[i][variable];
      }
    }
  }

  else {
    min_value[variable] = real_vector[0][variable];

    for (i = 1;i < nb_vector;i++) {
      if (real_vector[i][variable] < min_value[variable]) {
        min_value[variable] = real_vector[i][variable];
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

void Vectors::max_value_computation(int variable)

{
  register int i;


  if (type[variable] != REAL_VALUE) {
    max_value[variable] = int_vector[0][variable];

    for (i = 1;i < nb_vector;i++) {
      if (int_vector[i][variable] > max_value[variable]) {
        max_value[variable] = int_vector[i][variable];
      }
    }
  }

  else {
    max_value[variable] = real_vector[0][variable];

    for (i = 1;i < nb_vector;i++) {
      if (real_vector[i][variable] > max_value[variable]) {
        max_value[variable] = real_vector[i][variable];
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction de la loi marginale empirique pour une variable entiere positive.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Vectors::build_marginal_frequency_distribution(int variable)

{
  if ((type[variable] != REAL_VALUE) && (min_value[variable] >= 0) &&
      (max_value[variable] <= MARGINAL_DISTRIBUTION_MAX_VALUE)) {
    register int i;


    marginal_distribution[variable] = new FrequencyDistribution((int)max_value[variable] + 1);

    for (i = 0;i < nb_vector;i++) {
      (marginal_distribution[variable]->frequency[int_vector[i][variable]])++;
    }

    marginal_distribution[variable]->offset = (int)min_value[variable];
    marginal_distribution[variable]->nb_element = nb_vector;
    marginal_distribution[variable]->max_computation();
    marginal_distribution[variable]->mean_computation();
    marginal_distribution[variable]->variance_computation();

    mean[variable] = marginal_distribution[variable]->mean;
    covariance[variable][variable] = marginal_distribution[variable]->variance;
  }

  else {
    build_marginal_histogram(variable);

    mean_computation(variable);
    variance_computation(variable);
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction de l'histogramme marginal pour une variable.
 *
 *  arguments : indice de la variable, pas de regroupement, valeur minimum.
 *
 *--------------------------------------------------------------*/

void Vectors::build_marginal_histogram(int variable , double step , double imin_value)

{
  if ((!marginal_histogram[variable]) || (step != marginal_histogram[variable]->step) ||
      (imin_value != D_INF)) {
    register int i;


    // construction de l'histogramme

    if (step == D_DEFAULT) {
      step = MAX(::round((max_value[variable] - min_value[variable]) * HISTOGRAM_FREQUENCY / nb_vector) , 1);

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

      marginal_histogram[variable]->nb_element = nb_vector;
      marginal_histogram[variable]->type = type[variable];
    }

    marginal_histogram[variable]->step = step;
    marginal_histogram[variable]->min_value = imin_value;
    marginal_histogram[variable]->max_value = ceil(max_value[variable] / step) * step;

#   ifdef DEBUG
    cout << "\nTEST: " << marginal_histogram[variable]->min_value << " " << marginal_histogram[variable]->max_value
         << " | " << marginal_histogram[variable]->nb_category 
        << " " << (marginal_histogram[variable]->max_value - marginal_histogram[variable]->min_value) / marginal_histogram[variable]->step << endl;
#    endif

    // calcul des frequences

    for (i = 0;i < marginal_histogram[variable]->nb_category;i++) {
      marginal_histogram[variable]->frequency[i] = 0;
    }

    if (type[variable] != REAL_VALUE) {
      for (i = 0;i < nb_vector;i++) {
//        (marginal_histogram[variable]->frequency[(int)((int_vector[i][variable] - imin_value) / step)])++;
        (marginal_histogram[variable]->frequency[(int)floor((int_vector[i][variable] - imin_value) / step)])++;
      }
    }

    else {
      for (i = 0;i < nb_vector;i++) {
//        (marginal_histogram[variable]->frequency[(int)((real_vector[i][variable] - imin_value) / step)])++;
        (marginal_histogram[variable]->frequency[(int)floor((real_vector[i][variable] - imin_value) / step)])++;
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

bool Vectors::select_step(StatError &error , int variable , double step ,
                          double imin_value)

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
 *  Calcul d'un ordre sur les vecteurs a partir des valeurs prises
 *  par une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

int* Vectors::order_computation(int variable) const

{
  register int i , j;
  int int_min , int_value , *index;
  double real_min , real_value;


  index = new int[nb_vector];

  i = 0;

  if (type[variable] != REAL_VALUE) {
    do {

      // recherche de la valeur minimum courante

      if (i == 0) {
        int_value = (int)min_value[variable];
      }

      else {
        int_min = (int)max_value[variable] + 1;
        for (j = 0;j < nb_vector;j++) {
          if ((int_vector[j][variable] > int_value) && (int_vector[j][variable] < int_min)) {
            int_min = int_vector[j][variable];
          }
        }
        int_value = int_min;
      }

      // recherche des vecteurs prenant pour la variable selectionnee
      // la valeur minimum courante

      for (j = 0;j < nb_vector;j++) {
        if (int_vector[j][variable] == int_value) {
          index[i++] = j;
        }
      }
    }
    while (i < nb_vector);
  }

  else {
    do {

      // recherche de la valeur minimum courante

      if (i == 0) {
        real_value = min_value[variable];
      }

      else {
        real_min = max_value[variable] + 1;
        for (j = 0;j < nb_vector;j++) {
          if ((real_vector[j][variable] > real_value) && (real_vector[j][variable] < real_min)) {
            real_min = real_vector[j][variable];
          }
        }
        real_value = real_min;
      }

      // recherche des vecteurs prenant pour la variable selectionnee
      // la valeur minimum courante

      for (j = 0;j < nb_vector;j++) {
        if (real_vector[j][variable] == real_value) {
          index[i++] = j;
        }
      }
    }
    while (i < nb_vector);
  }

  return index;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la fonction de repartition empirique pour une variable.
 *
 *  arguments : indice de la variable, (valeurs, fonction de repartition).
 *
 *--------------------------------------------------------------*/

int Vectors::cumulative_distribution_function_computation(int variable , double **cdf) const

{
  register int i , j;
  int cumul , int_min , int_value , frequency;
  double real_min , real_value;


  if (marginal_distribution[variable]) {
    i = marginal_distribution[variable]->cumulative_distribution_function_computation(cdf);
  }

  else {
    cdf[0] = new double[nb_vector];
    cdf[1] = new double[nb_vector];

    cumul = 0;
    i = 0;

    if (type[variable] != REAL_VALUE) {

      do {

        // recherche de la valeur minimum courante

        if (cumul == 0) {
          int_value = (int)min_value[variable];
        }

        else {
          int_min = (int)max_value[variable] + 1;
          for (j = 0;j < nb_vector;j++) {
            if ((int_vector[j][variable] > int_value) &&
                (int_vector[j][variable] < int_min)) {
              int_min = int_vector[j][variable];
            }
          }
          int_value = int_min;
        }

        // recherche du nombre de vecteurs prenant pour la variable selectionnee
        // la valeur minimum courante

        frequency = 0;
        for (j = 0;j < nb_vector;j++) {
          if (int_vector[j][variable] == int_value) {
            frequency++;
          }
        }

        cdf[0][i] = int_value;
        cdf[1][i] = (cumul + (double)(frequency + 1) / 2.) / (double)nb_vector;
        cumul += frequency;
        i++;
      }
      while (cumul < nb_vector);
    }

    else {
      do {

        // recherche de la valeur minimum courante

        if (cumul == 0) {
          real_value = min_value[variable];
        }

        else {
          real_min = max_value[variable] + 1;
          for (j = 0;j < nb_vector;j++) {
            if ((real_vector[j][variable] > real_value) &&
                (real_vector[j][variable] < real_min)) {
              real_min = real_vector[j][variable];
            }
          }
          real_value = real_min;
        }

        // recherche du nombre de vecteurs prenant pour la variable selectionnee
        // la valeur minimum courante

        frequency = 0;
        for (j = 0;j < nb_vector;j++) {
          if (real_vector[j][variable] == real_value) {
            frequency++;
          }
        }

        cdf[0][i] = real_value;
        cdf[1][i] = (cumul + (double)(frequency + 1) / 2.) / (double)nb_vector;
        cumul += frequency;
        i++;
      }
      while (cumul < nb_vector);
    }
  }

  return i;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'intervalle minimum entre 2 valeurs pour une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Vectors::min_interval_computation(int variable)

{
  if (marginal_distribution[variable]) {
    min_interval[variable] = marginal_distribution[variable]->min_interval_computation();
  }

  else {
    register int i , j;
    int int_min , int_value;
    double real_min , real_value;


    min_interval[variable] = max_value[variable] - min_value[variable];
    i = 0;

    if (type[variable] != REAL_VALUE) {
      do {

        // recherche de la valeur minimum courante

        if (i == 0) {
          int_value = (int)min_value[variable];
        }

        else {
          int_min = (int)max_value[variable] + 1;
          for (j = 0;j < nb_vector;j++) {
            if ((int_vector[j][variable] > int_value) &&
                (int_vector[j][variable] < int_min)) {
              int_min = int_vector[j][variable];
            }
          }

          if (int_min - int_value < min_interval[variable]) {
            min_interval[variable] = int_min - int_value;
          }
          int_value = int_min;
        }

        // recherche du nombre de vecteurs prenant pour la variable selectionnee
        // la valeur minimum courante

        for (j = 0;j < nb_vector;j++) {
          if (int_vector[j][variable] == int_value) {
            i++;
          }
        }
      }
      while (i < nb_vector);
    }

    else {
      do {

        // recherche de la valeur minimum courante

        if (i == 0) {
          real_value = min_value[variable];
        }

        else {
          real_min = max_value[variable] + 1;
          for (j = 0;j < nb_vector;j++) {
            if ((real_vector[j][variable] > real_value) &&
                (real_vector[j][variable] < real_min)) {
              real_min = real_vector[j][variable];
            }
          }

          if (real_min - real_value < min_interval[variable]) {
            min_interval[variable] = real_min - real_value;
          }
          real_value = real_min;
        }

        // recherche du nombre de vecteurs prenant pour la variable selectionnee
        // la valeur minimum courante

        for (j = 0;j < nb_vector;j++) {
          if (real_vector[j][variable] == real_value) {
            i++;
          }
        }
      }
      while (i < nb_vector);
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la moyenne d'une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Vectors::mean_computation(int variable)

{
  register int i;


  mean[variable] = 0.;

  if (type[variable] != REAL_VALUE) {
    for (i = 0;i < nb_vector;i++) {
      mean[variable] += int_vector[i][variable];
    }
  }

  else {
    for (i = 0;i < nb_vector;i++) {
      mean[variable] += real_vector[i][variable];
    }
  }

  mean[variable] /= nb_vector;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la variance d'une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Vectors::variance_computation(int variable)

{
  if (mean[variable] != D_INF) {
    covariance[variable][variable] = 0.;

    if (nb_vector > 1) {
      register int i;
      double diff;


      if (type[variable] != REAL_VALUE) {
        for (i = 0;i < nb_vector;i++) {
          diff = int_vector[i][variable] - mean[variable];
          covariance[variable][variable] += diff * diff;
        }
      }

      else{
        for (i = 0;i < nb_vector;i++) {
          diff = real_vector[i][variable] - mean[variable];
          covariance[variable][variable] += diff * diff;
        }
      }

      covariance[variable][variable] /= (nb_vector - 1);
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la matrice de variance-covariance.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Vectors::covariance_computation(int variable)

{
  if (mean[0] != D_INF) {
    register int i , j , k;


    for (i = 0;i < nb_variable;i++) {
      if ((variable == I_DEFAULT) || (i == variable)) {
        for (j = (variable == I_DEFAULT ? i + 1 : 0);j < nb_variable;j++) {
          covariance[i][j] = 0.;

          if (nb_vector > 1) {
            if ((type[i] != REAL_VALUE) && (type[j] != REAL_VALUE)) {
              for (k = 0;k < nb_vector;k++) {
                covariance[i][j] += (int_vector[k][i] - mean[i]) * (int_vector[k][j] - mean[j]);
              }
            }
            else if ((type[i] != REAL_VALUE) && (type[j] == REAL_VALUE)) {
              for (k = 0;k < nb_vector;k++) {
                covariance[i][j] += (int_vector[k][i] - mean[i]) * (real_vector[k][j] - mean[j]);
              }
            }
            else if ((type[i] == REAL_VALUE) && (type[j] != REAL_VALUE)) {
              for (k = 0;k < nb_vector;k++) {
                covariance[i][j] += (real_vector[k][i] - mean[i]) * (int_vector[k][j] - mean[j]);
              }
            }
//            else if ((type[i] == REAL_VALUE) && (type[j] == REAL_VALUE)) {
            else {
              for (k = 0;k < nb_vector;k++) {
                covariance[i][j] += (real_vector[k][i] - mean[i]) * (real_vector[k][j] - mean[j]);
              }
            }

            covariance[i][j] /= (nb_vector - 1);
          }

          covariance[j][i] = covariance[i][j];
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'ecart absolu moyen d'une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

double Vectors::mean_absolute_deviation_computation(int variable) const

{
  register int i;
  double mean_absolute_deviation = D_DEFAULT;


  if (mean[variable] != D_INF) {
    mean_absolute_deviation = 0.;

    if (type[variable] != REAL_VALUE) {
      for (i = 0;i < nb_vector;i++) {
        mean_absolute_deviation += fabs(int_vector[i][variable] - mean[variable]);
      }
    }

    else {
      for (i = 0;i < nb_vector;i++) {
        mean_absolute_deviation += fabs(real_vector[i][variable] - mean[variable]);
      }
    }

    mean_absolute_deviation /= nb_vector;
  }

  return mean_absolute_deviation;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la difference absolue moyenne pour une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

double Vectors::mean_absolute_difference_computation(int variable) const

{
  register int i , j;
  double mean_absolute_difference;


  mean_absolute_difference = 0.;

  if (nb_vector > 1) {
    if (type[variable] != REAL_VALUE) {
      for (i = 0;i < nb_vector;i++) {
        for (j = i + 1;j < nb_vector;j++) {
          mean_absolute_difference += abs(int_vector[i][variable] - int_vector[j][variable]);
        }
      }
    }

    else {
      for (i = 0;i < nb_vector;i++) {
        for (j = i + 1;j < nb_vector;j++) {
          mean_absolute_difference += fabs(real_vector[i][variable] - real_vector[j][variable]);
        }
      }
    }

    mean_absolute_difference = 2 * mean_absolute_difference / (nb_vector * (nb_vector - 1));
  }

  return mean_absolute_difference;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du coefficient d'asymetrie d'une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

double Vectors::skewness_computation(int variable) const

{
  register int i;
  double skewness = D_INF , diff;


  if ((mean[variable] != D_INF) && (covariance[variable][variable] != D_DEFAULT)) {
    skewness = 0.;

    if ((nb_vector > 2) && (covariance[variable][variable] > 0.)) {
      if (type[variable] != REAL_VALUE) {
        for (i = 0;i < nb_vector;i++) {
          diff = int_vector[i][variable] - mean[variable];
          skewness += diff * diff * diff;
        }
      }

      else {
        for (i = 0;i < nb_vector;i++) {
          diff = real_vector[i][variable] - mean[variable];
          skewness += diff * diff * diff;
        }
      }

      skewness = skewness * nb_vector / ((nb_vector - 1) * (nb_vector - 2) *
                  pow(covariance[variable][variable] , 1.5));
    }
  }

  return skewness;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'exces d'applatissement d'une variable :
 *  exces d'applatissement = coefficient d'applatissement - 3.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

double Vectors::kurtosis_computation(int variable) const

{
  register int i;
  double kurtosis = D_INF , diff;


  if ((mean[variable] != D_INF) && (covariance[variable][variable] != D_DEFAULT)) {
    if (covariance[variable][variable] == 0.) {
      kurtosis = -2.;
    }

    else {
      kurtosis = 0.;

      if (type[variable] != REAL_VALUE) {
        for (i = 0;i < nb_vector;i++) {
          diff = int_vector[i][variable] - mean[variable];
          kurtosis += diff * diff * diff * diff;
        }
      }

      else {
        for (i = 0;i < nb_vector;i++) {
          diff = real_vector[i][variable] - mean[variable];
          kurtosis += diff * diff * diff * diff;
        }
      }

      kurtosis = kurtosis / ((nb_vector - 1) * covariance[variable][variable] *
                  covariance[variable][variable]) - 3.;
    }
  }

  return kurtosis;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la matrice des coefficients de correlation a partir
 *  de la matrice de variance-covariance.
 *
 *--------------------------------------------------------------*/

double** Vectors::correlation_computation() const

{
  register int i , j;
  double norm , **correlation = NULL;


  if (covariance[0][0] != D_DEFAULT) {
    correlation = new double*[nb_variable];
    for (i = 0;i < nb_variable;i++) {
      correlation[i] = new double[nb_variable];
    }

    for (i = 0;i < nb_variable;i++) {
      correlation[i][i] = 1.;
      for (j = i + 1;j < nb_variable;j++) {
        correlation[i][j] = covariance[i][j];
        norm = covariance[i][i] * covariance[j][j];
        if (norm > 0.) {
          correlation[i][j] /= sqrt(norm);
        }
        correlation[j][i] = correlation[i][j];
      }
    }
  }

  return correlation;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la direction moyenne d'une variable circulaire.
 *
 *  arguments : indice de la variable, unite (DEGREE/RADIAN).
 *
 *--------------------------------------------------------------*/

double* Vectors::mean_direction_computation(int variable , int unit) const

{
  register int i;
  double *mean_direction;


  mean_direction = new double[4];
//  mean_direction = new double[2];

  mean_direction[0] = 0.;
  mean_direction[1] = 0.;

  switch (type[variable]) {

  case INT_VALUE : {
    switch (unit) {

    case DEGREE : {
      for (i = 0;i < nb_vector;i++) {
        mean_direction[0] += cos(int_vector[i][variable] * M_PI / 180);
        mean_direction[1] += sin(int_vector[i][variable] * M_PI / 180);
      }
      break;
    }

    case RADIAN : {
      for (i = 0;i < nb_vector;i++) {
        mean_direction[0] += cos(int_vector[i][variable]);
        mean_direction[1] += sin(int_vector[i][variable]);
      }
      break;
    }
    }
    break;
  }

  case REAL_VALUE : {
    switch (unit) {

    case DEGREE : {
      for (i = 0;i < nb_vector;i++) {
        mean_direction[0] += cos(real_vector[i][variable] * M_PI / 180);
        mean_direction[1] += sin(real_vector[i][variable] * M_PI / 180);
      }
      break;
    }

    case RADIAN : {
      for (i = 0;i < nb_vector;i++) {
        mean_direction[0] += cos(real_vector[i][variable]);
        mean_direction[1] += sin(real_vector[i][variable]);
      }
      break;
    }
    }
    break;
  }
  }

  mean_direction[0] /= nb_vector;
  mean_direction[1] /= nb_vector;

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


};  // namespace stat_tool
