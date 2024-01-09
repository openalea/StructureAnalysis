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
 *       $Id: sequences1.cpp 18073 2015-04-23 10:56:03Z guedon $
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



#include <sstream>

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

#include "renewal.h"
#include "sequences.h"
// #include "tops.h"
#include "sequence_label.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Sequences.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences()

{
  nb_sequence = 0;
  identifier = NULL;

  max_length = 0;
  cumul_length = 0;
  length = NULL;
  length_distribution = NULL;

  vertex_identifier = NULL;

  index_parameter_type = IMPLICIT_TYPE;
  index_parameter_distribution = NULL;
  index_interval = NULL;
  index_parameter = NULL;

  nb_variable = 0;

  type = NULL;
  min_value = NULL;
  max_value = NULL;
  marginal_distribution = NULL;
  marginal_histogram = NULL;

  int_sequence = NULL;
  real_sequence = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Sequences.
 *
 *  arguments : nombre de sequences, nombre de variables.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences(int inb_sequence , int inb_variable)

{
  register int i , j;


  nb_sequence = inb_sequence;

  identifier = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    identifier[i] = i + 1;
  }

  length = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    length[i] = 0;
  }

  max_length = 0;
  cumul_length = 0;
  length_distribution = NULL;

  vertex_identifier = NULL;

  index_parameter_type = IMPLICIT_TYPE;
  index_parameter_distribution = NULL;
  index_interval = NULL;
  index_parameter = NULL;

  nb_variable = inb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = INT_VALUE;
    min_value[i] = 0.;
    max_value[i] = 0.;
    marginal_distribution[i] = NULL;
    marginal_histogram[i] = NULL;
  }

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      int_sequence[i][j] = NULL;
      real_sequence[i][j] = NULL;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Initialisation d'un objet Sequences.
 *
 *  arguments : nombre de sequences, identificateurs des sequences,
 *              longueurs des sequences, identificateurs des vertex d'un MTG associe,
 *              type de parametre d'index, nombre de variables, type de chaque variable,
 *              flag copy of vertex identifiers, flag initialization.
 *
 *--------------------------------------------------------------*/

void Sequences::init(int inb_sequence , int *iidentifier , int *ilength ,
                     int **ivertex_identifier , int iindex_parameter_type ,
                     int inb_variable , int *itype , bool vertex_identifier_copy ,
                     bool init_flag)

{
  register int i , j , k;
  int blength;


  nb_sequence = inb_sequence;

  identifier = new int[nb_sequence];

  if (iidentifier) {
    for (i = 0;i < nb_sequence;i++) {
      identifier[i] = iidentifier[i];
    }
  }

  else {
    for (i = 0;i < nb_sequence;i++) {
      identifier[i] = i + 1;
    }
  }

  length = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    length[i] = ilength[i];
  }

  max_length_computation();
  cumul_length_computation();
  build_length_frequency_distribution();

  if (ivertex_identifier) {
    vertex_identifier = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      vertex_identifier[i] = new int[length[i]];

      if (vertex_identifier_copy) {
        for (j = 0;j < length[i];j++) {
          vertex_identifier[i][j] = ivertex_identifier[i][j];
        }
      }
    }
  }

  else {
    vertex_identifier = NULL;
  }

  index_parameter_type = iindex_parameter_type;
  index_parameter_distribution = NULL;
  index_interval = NULL;

  if (index_parameter_type != IMPLICIT_TYPE) {
    index_parameter = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      blength = ((index_parameter_type == POSITION) || (index_parameter_type == POSITION_INTERVAL) ? length[i] + 1 : length[i]);
      index_parameter[i] = new int[blength];

      if (init_flag) {
        for (j = 0;j < blength;j++) {
          index_parameter[i][j] = 0;
        }
      }
    }
  }

  else {
    index_parameter = NULL;
  }

  nb_variable = inb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = itype[i];
    min_value[i] = 0.;
    max_value[i] = 0.;
    marginal_distribution[i] = NULL;
    marginal_histogram[i] = NULL;
  }

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      if ((type[j] != REAL_VALUE) && (type[j] != AUXILIARY)) {
        int_sequence[i][j] = new int[length[i]];
        real_sequence[i][j] = NULL;

        if (init_flag) {
          for (k = 0;k < length[i];k++) {
            int_sequence[i][j][k] = 0;
          }
        }
      }

      else {
        int_sequence[i][j] = NULL;
        real_sequence[i][j] = new double[length[i]];

        if (init_flag) {
          for (k = 0;k < length[i];k++) {
            real_sequence[i][j][k] = 0.;
          }
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Initialisation d'un objet Sequences.
 *
 *  arguments : nombre de sequences, identificateurs des sequences,
 *              longueurs des sequences, nombre de variables, flag initialisation.
 *
 *--------------------------------------------------------------*/

void Sequences::init(int inb_sequence , int *iidentifier , int *ilength ,
                     int inb_variable , bool init_flag)

{
  register int i , j , k;


  nb_sequence = inb_sequence;

  identifier = new int[nb_sequence];
  if (iidentifier) {
    for (i = 0;i < nb_sequence;i++) {
      identifier[i] = iidentifier[i];
    }
  }
  else {
    for (i = 0;i < nb_sequence;i++) {
      identifier[i] = i + 1;
    }
  }

  length = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    length[i] = ilength[i];
  }

  max_length_computation();
  cumul_length_computation();
  build_length_frequency_distribution();

  vertex_identifier = NULL;

  index_parameter_type = IMPLICIT_TYPE;
  index_parameter_distribution = NULL;
  index_interval = NULL;
  index_parameter = NULL;

  nb_variable = inb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = INT_VALUE;
    min_value[i] = 0.;
    max_value[i] = 0.;
    marginal_distribution[i] = NULL;
    marginal_histogram[i] = NULL;
  }

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      int_sequence[i][j] = new int[length[i]];
      real_sequence[i][j] = NULL;

      if (init_flag) {
        for (k = 0;k < length[i];k++) {
          int_sequence[i][j][k] = 0;
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Sequences.
 *
 *  arguments : nombre de sequences, identificateurs des sequences,
 *              longueurs des sequences, type de parametre d'index (TIME/POSITION),
 *              nombre de variables, type des variables,
 *              (parametres d'index et) sequences entieres.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences(int inb_sequence , int *iidentifier , int *ilength ,
                     int iindex_parameter_type , int inb_variable , int itype ,
                     int ***iint_sequence)

{
  register int i , j , k;
  int *btype , *pisequence , *cisequence;


  btype = new int[inb_variable];
  for (i = 0;i < inb_variable;i++) {
    btype[i] = itype;
  }

  init(inb_sequence , iidentifier , ilength , NULL , iindex_parameter_type ,
       inb_variable , btype , false , false);
  delete [] btype;

//  if (index_parameter_type != IMPLICIT_TYPE) {
  if (index_parameter) {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < (index_parameter_type == POSITION ? length[i] + 1 : length[i]);j++) {
        index_parameter[i][j] = iint_sequence[i][0][j];
      }
    }

    build_index_parameter_frequency_distribution();

//    if ((index_parameter_type == TIME) || ((index_parameter_type == POSITION) &&
//        (type[0] != NB_INTERNODE))) {
    if ((index_parameter_type == TIME) || (index_parameter_type == POSITION)) {
      index_interval_computation();
    }
  }

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < nb_variable;j++) {
      pisequence = int_sequence[i][j];
      if (index_parameter) {
        cisequence = iint_sequence[i][j + 1];
      }
      else {
        cisequence = iint_sequence[i][j];
      }

      for (k = 0;k < length[i];k++) {
        *pisequence++ = *cisequence++;
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value_computation(i);
    max_value_computation(i);

    build_marginal_frequency_distribution(i);
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Sequences.
 *
 *  arguments : nombre de sequences, identificateurs des sequences, longueurs des sequences,
 *              nombre de variables, sequences reelles.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences(int inb_sequence , int *iidentifier , int *ilength ,
                     int inb_variable , double ***ireal_sequence)

{
  register int i , j , k;
  int *itype;


  itype = new int[inb_variable];
  for (i = 0;i < inb_variable;i++) {
    itype[i] = REAL_VALUE;
  }

  init(inb_sequence , iidentifier , ilength , NULL , IMPLICIT_TYPE ,
       inb_variable , itype , false , false);
  delete [] itype;

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < nb_variable;j++) {
      for (k = 0;k < length[i];k++) {
        real_sequence[i][j][k] = ireal_sequence[i][j][k];
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value_computation(i);
    max_value_computation(i);

    build_marginal_histogram(i);
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Sequences.
 *
 *  arguments : nombre de sequences, identificateurs des sequences,
 *              longueurs des sequences, identificateurs des vertex d'un MTG associe,
 *              type du parametre d'index (TIME/POSITION),
 *              parametres d'index, nombre de variables, type des variables,
 *              sequences entieres, sequences reelles.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences(int inb_sequence , int *iidentifier , int *ilength ,
                     int **ivertex_identifier , int iindex_parameter_type ,
                     int **iindex_parameter , int inb_variable , int *itype ,
                     int ***iint_sequence , double ***ireal_sequence)

{
  register int i , j , k , m , n;


  init(inb_sequence , iidentifier , ilength , ivertex_identifier ,
       iindex_parameter_type , inb_variable , itype , true , false);

//  if (index_parameter_type != IMPLICIT_TYPE) {
  if (index_parameter) {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < (index_parameter_type == POSITION ? length[i] + 1 : length[i]);j++) {
        index_parameter[i][j] = iindex_parameter[i][j];
      }
    }

    build_index_parameter_frequency_distribution();

    if ((index_parameter_type == TIME) || (index_parameter_type == POSITION)) {
      index_interval_computation();
    }
  }

  i = 0;
  j = 0;
  for (k = 0;k < nb_variable;k++) {
    switch (type[k]) {

    case INT_VALUE : {
      for (m = 0;m < nb_sequence;m++) {
        for (n = 0;n < length[m];n++) {
          int_sequence[m][k][n] = iint_sequence[m][i][n];
        }
      }
      i++;
      break;
    }

    case REAL_VALUE : {
      for (m = 0;m < nb_sequence;m++) {
        for (n = 0;n < length[m];n++) {
          real_sequence[m][k][n] = ireal_sequence[m][j][n];
        }
      }
      j++;
      break;
    }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value_computation(i);
    max_value_computation(i);

    switch (type[i]) {
    case INT_VALUE :
      build_marginal_frequency_distribution(i);
      break;
    case REAL_VALUE :
      build_marginal_histogram(i);
      break;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Sequences.
 *
 *  arguments : loi empirique des longueurs des sequences, nombre de variables,
 *              type de chaque variable, flag initialisation.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences(const FrequencyDistribution &ilength_distribution ,
                     int inb_variable , int *itype , bool init_flag)

{
  register int i , j , k;
  int *plength;


  nb_sequence = ilength_distribution.nb_element;

  identifier = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    identifier[i] = i + 1;
  }

  length = new int[nb_sequence];
  plength = length;
  for (i = ilength_distribution.offset;i < ilength_distribution.nb_value;i++) {
    for (j = 0;j < ilength_distribution.frequency[i];j++) {
      *plength++ = i;
    }
  }

  max_length = ilength_distribution.nb_value - 1;
  cumul_length_computation();
  length_distribution = new FrequencyDistribution(ilength_distribution);

  vertex_identifier = NULL;

  index_parameter_type = IMPLICIT_TYPE;
  index_parameter_distribution = NULL;
  index_interval = NULL;
  index_parameter = NULL;

  nb_variable = inb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  if (itype) {
    for (i = 0;i < nb_variable;i++) {
      type[i] = itype[i];
    }
  }

  else {
    type[0] = STATE;
    for (i = 1;i < nb_variable;i++) {
      type[i] = INT_VALUE;
    }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value[i] = 0.;
    max_value[i] = 0.;
    marginal_distribution[i] = NULL;
    marginal_histogram[i] = NULL;
  }

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      if (type[j] != REAL_VALUE) {
        int_sequence[i][j] = new int[length[i]];
        real_sequence[i][j] = NULL;

        if (init_flag) {
          for (k = 0;k < length[i];k++) {
            int_sequence[i][j][k] = 0;
          }
        }
      }

      else {
        int_sequence[i][j] = NULL;
        real_sequence[i][j] = new double[length[i]];

        if (init_flag) {
          for (k = 0;k < length[i];k++) {
            real_sequence[i][j][k] = 0.;
          }
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Sequences a partir d'un objet RenewalData.
 *
 *  argument : reference sur un objet RenewalData.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences(const RenewalData &timev)

{
  register int i , j;


  nb_sequence = timev.nb_element;

  identifier = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    identifier[i] = i + 1;
  }

  length = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    length[i] = timev.length[i];
  }

  max_length_computation();
  cumul_length_computation();
  build_length_frequency_distribution();

  vertex_identifier = NULL;

  index_parameter_type = IMPLICIT_TYPE;
  index_parameter_distribution = NULL;
  index_interval = NULL;
  index_parameter = NULL;

  nb_variable = 1;

  type = new int[nb_variable];
  type[0] = INT_VALUE;

  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];

    int_sequence[i][0] = new int[length[i]];
    real_sequence[i][0] = NULL;

    for (j = 0;j < length[i];j++) {
      int_sequence[i][0][j] = timev.sequence[i][j];
    }
  }

  min_value_computation(0);
  max_value_computation(0);
  build_marginal_frequency_distribution(0);
  marginal_histogram[0] = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Sequences.
 *
 *  argument : reference sur un objet Sequences, indice de la variable,
 *             type de la variable selectionnee.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences(const Sequences &seq , int variable , int itype)

{
  register int i , j , k;
  int blength;


  nb_sequence = seq.nb_sequence;

  identifier = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    identifier[i] = seq.identifier[i];
  }

  max_length = seq.max_length;
  cumul_length = seq.cumul_length;

  length = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    length[i] = seq.length[i];
  }

  length_distribution = new FrequencyDistribution(*(seq.length_distribution));

  if (seq.vertex_identifier) {
    vertex_identifier = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      vertex_identifier[i] = new int[length[i]];
      for (j = 0;j < length[i];j++) {
        vertex_identifier[i][j] = seq.vertex_identifier[i][j];
      }
    }
  }

  else {
    vertex_identifier = NULL;
  }

  index_parameter_type = seq.index_parameter_type;

  if (seq.index_parameter_distribution) {
    index_parameter_distribution = new FrequencyDistribution(*(seq.index_parameter_distribution));
  }
  else {
    index_parameter_distribution = NULL;
  }

  if (seq.index_interval) {
    index_interval = new FrequencyDistribution(*(seq.index_interval));
  }
  else {
    index_interval = NULL;
  }

  if (seq.index_parameter) {
    index_parameter = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      blength = (index_parameter_type == POSITION ? length[i] + 1 : length[i]);
      index_parameter[i] = new int[blength];
      for (j = 0;j < blength;j++) {
        index_parameter[i][j] = seq.index_parameter[i][j];
      }
    }
  }

  else {
    index_parameter = NULL;
  }

  nb_variable = seq.nb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    if (i != variable) {
      type[i] = seq.type[i];
      min_value[i] = seq.min_value[i];
      max_value[i] = seq.max_value[i];

      if (seq.marginal_distribution[i]) {
        marginal_distribution[i] = new FrequencyDistribution(*(seq.marginal_distribution[i]));
      }
      else {
        marginal_distribution[i] = NULL;
      }

      if (seq.marginal_histogram[i]) {
        marginal_histogram[i] = new Histogram(*(seq.marginal_histogram[i]));
      }
      else {
        marginal_histogram[i] = NULL;
      }
    }

    else {
      type[i] = itype;
      min_value[i] = 0.;
      max_value[i] = 0.;
      marginal_distribution[i] = NULL;
      marginal_histogram[i] = NULL;
    }
  }

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      if ((type[j] != REAL_VALUE) && (type[j] != AUXILIARY)) {
        int_sequence[i][j] = new int[length[i]];
        real_sequence[i][j] = NULL;

        if (j != variable) {
          for (k = 0;k < length[i];k++) {
            int_sequence[i][j][k] = seq.int_sequence[i][j][k];
          }
        }
      }

      else {
        int_sequence[i][j] = NULL;
        real_sequence[i][j] = new double[length[i]];

        if (j != variable) {
          for (k = 0;k < length[i];k++) {
            real_sequence[i][j][k] = seq.real_sequence[i][j][k];
          }
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Sequences.
 *
 *  arguments : reference sur un objet Sequences, nombre de sequences,
 *              indices des sequences selectionnees.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences(const Sequences &seq , int inb_sequence , int *index)

{
  register int i , j , k;
  int blength;


  nb_sequence = inb_sequence;

  identifier = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    identifier[i] = seq.identifier[index[i]];
  }

  length = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    length[i] = seq.length[index[i]];
  }

  max_length_computation();
  cumul_length_computation();
  build_length_frequency_distribution();

  if (seq.vertex_identifier) {
    vertex_identifier = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      vertex_identifier[i] = new int[length[i]];
      for (j = 0;j < length[i];j++) {
        vertex_identifier[i][j] = seq.vertex_identifier[index[i]][j];
      }
    }
  }

  else {
    vertex_identifier = NULL;
  }

  index_parameter_type = seq.index_parameter_type;

//  if (index_parameter_type != IMPLICIT_TYPE) {
  if (seq.index_parameter) {
    index_parameter = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      blength = (index_parameter_type == POSITION ? length[i] + 1 : length[i]);
      index_parameter[i] = new int[blength];
      for (j = 0;j < blength;j++) {
        index_parameter[i][j] = seq.index_parameter[index[i]][j];
      }
    }

    build_index_parameter_frequency_distribution();
  }

  else {
    index_parameter_distribution = NULL;
    index_interval = NULL;
    index_parameter = NULL;
  }

  nb_variable = seq.nb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = seq.type[i];
    marginal_distribution[i] = NULL;
    marginal_histogram[i] = NULL;
  }

  if (index_parameter) {
    index_interval_computation();
  }

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      if ((type[j] != REAL_VALUE) && (type[j] != AUXILIARY)) {
        int_sequence[i][j] = new int[length[i]];
        real_sequence[i][j] = NULL;

        for (k = 0;k < length[i];k++) {
          int_sequence[i][j][k] = seq.int_sequence[index[i]][j][k];
        }
      }

      else {
        int_sequence[i][j] = NULL;
        real_sequence[i][j] = new double[length[i]];

        for (k = 0;k < length[i];k++) {
          real_sequence[i][j][k] = seq.real_sequence[index[i]][j][k];
        }
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value_computation(i);
    max_value_computation(i);

    if (type[i] != AUXILIARY) {
      if (type[i] != REAL_VALUE) {
        build_marginal_frequency_distribution(i);
      }
      else {
        build_marginal_histogram(i , seq.marginal_histogram[i]->step);
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Sequences avec ajout de variables auxilliaires.
 *
 *  arguments : reference sur un objet Sequences, flags sur l'ajout
 *              des variables auxilliaires.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences(const Sequences &seq , bool *auxiliary)

{
  register int i , j , k , m;
  int blength;


  nb_sequence = seq.nb_sequence;

  identifier = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    identifier[i] = seq.identifier[i];
  }

  max_length = seq.max_length;
  cumul_length = seq.cumul_length;

  length = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    length[i] = seq.length[i];
  }

  length_distribution = new FrequencyDistribution(*(seq.length_distribution));

  if (seq.vertex_identifier) {
    vertex_identifier = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      vertex_identifier[i] = new int[length[i]];
      for (j = 0;j < length[i];j++) {
        vertex_identifier[i][j] = seq.vertex_identifier[i][j];
      }
    }
  }

  else {
    vertex_identifier = NULL;
  }

  index_parameter_type = seq.index_parameter_type;

  if (seq.index_parameter_distribution) {
    index_parameter_distribution = new FrequencyDistribution(*(seq.index_parameter_distribution));
  }
  else {
    index_parameter_distribution = NULL;
  }

  if (seq.index_interval) {
    index_interval = new FrequencyDistribution(*(seq.index_interval));
  }
  else {
    index_interval = NULL;
  }

  if (seq.index_parameter) {
    index_parameter = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      blength = (index_parameter_type == POSITION ? length[i] + 1 : length[i]);
      index_parameter[i] = new int[blength];
      for (j = 0;j < blength;j++) {
        index_parameter[i][j] = seq.index_parameter[i][j];
      }
    }
  }

  else {
    index_parameter = NULL;
  }

  nb_variable = seq.nb_variable;
  for (i = 0;i < seq.nb_variable;i++) {
    if (auxiliary[i]) {
      nb_variable++;
    }
  }

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  i = 0;
  for (j = 0;j < seq.nb_variable;j++) {
    type[i] = seq.type[j];
    min_value[i] = seq.min_value[j];
    max_value[i] = seq.max_value[j];

    if (seq.marginal_distribution[j]) {
      marginal_distribution[i] = new FrequencyDistribution(*(seq.marginal_distribution[j]));
    }
    else {
      marginal_distribution[i] = NULL;
    }

    if (seq.marginal_histogram[j]) {
      marginal_histogram[i] = new Histogram(*(seq.marginal_histogram[j]));
    }
    else {
      marginal_histogram[i] = NULL;
    }
    i++;

    if (auxiliary[j]) {
      type[i] = AUXILIARY;
      min_value[i] = 0.;
      max_value[i] = 0.;
      marginal_distribution[i] = NULL;
      marginal_histogram[i] = NULL;
      i++;
    }
  }

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];
    j = 0;

    for (k = 0;k < seq.nb_variable;k++) {
      if (seq.type[k] != REAL_VALUE) {
        int_sequence[i][j] = new int[length[i]];
        real_sequence[i][j] = NULL;

        for (m = 0;m < length[i];m++) {
          int_sequence[i][j][m] = seq.int_sequence[i][k][m];
        }
      }

      else {
        int_sequence[i][j] = NULL;
        real_sequence[i][j] = new double[length[i]];

        for (m = 0;m < length[i];m++) {
          real_sequence[i][j][m] = seq.real_sequence[i][k][m];
        }
      }

      j++;

      if (auxiliary[k]) {
        int_sequence[i][j] = NULL;
        real_sequence[i][j] = new double[length[i]];
        j++;
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Sequences.
 *
 *  argument : reference sur un objet Sequences.
 *
 *--------------------------------------------------------------*/

void Sequences::copy(const Sequences &seq)

{
  register int i , j , k;
  int blength;


  nb_sequence = seq.nb_sequence;

  identifier = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    identifier[i] = seq.identifier[i];
  }

  max_length = seq.max_length;
  cumul_length = seq.cumul_length;

  length = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    length[i] = seq.length[i];
  }

  length_distribution = new FrequencyDistribution(*(seq.length_distribution));

  if (seq.vertex_identifier) {
    vertex_identifier = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      vertex_identifier[i] = new int[length[i]];
      for (j = 0;j < length[i];j++) {
        vertex_identifier[i][j] = seq.vertex_identifier[i][j];
      }
    }
  }

  else {
    vertex_identifier = NULL;
  }

  index_parameter_type = seq.index_parameter_type;

  if (seq.index_parameter_distribution) {
    index_parameter_distribution = new FrequencyDistribution(*(seq.index_parameter_distribution));
  }
  else {
    index_parameter_distribution = NULL;
  }

  if (seq.index_interval) {
    index_interval = new FrequencyDistribution(*(seq.index_interval));
  }
  else {
    index_interval = NULL;
  }

  if (seq.index_parameter) {
    index_parameter = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      blength = (index_parameter_type == POSITION ? length[i] + 1 : length[i]);
      index_parameter[i] = new int[blength];
      for (j = 0;j < blength;j++) {
        index_parameter[i][j] = seq.index_parameter[i][j];
      }
    }
  }

  else {
    index_parameter = NULL;
  }

  nb_variable = seq.nb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = seq.type[i];
    min_value[i] = seq.min_value[i];
    max_value[i] = seq.max_value[i];

    if (seq.marginal_distribution[i]) {
      marginal_distribution[i] = new FrequencyDistribution(*(seq.marginal_distribution[i]));
    }
    else {
      marginal_distribution[i] = NULL;
    }

    if (seq.marginal_histogram[i]) {
      marginal_histogram[i] = new Histogram(*(seq.marginal_histogram[i]));
    }
    else {
      marginal_histogram[i] = NULL;
    }
  }

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      if ((type[j] != REAL_VALUE) && (type[j] != AUXILIARY)) {
        int_sequence[i][j] = new int[length[i]];
        real_sequence[i][j] = NULL;

        for (k = 0;k < length[i];k++) {
          int_sequence[i][j][k] = seq.int_sequence[i][j][k];
        }
      }

      else {
        int_sequence[i][j] = NULL;
        real_sequence[i][j] = new double[length[i]];

        for (k = 0;k < length[i];k++) {
          real_sequence[i][j][k] = seq.real_sequence[i][j][k];
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Sequences avec inversion du sens de parcours des sequences.
 *
 *  arguments : reference sur un objet Sequences.
 *
 *--------------------------------------------------------------*/

void Sequences::reverse(const Sequences &seq)

{
  register int i , j , k;
  int blength , end_position , *pidentifier , *cidentifier , *pindex_param ,
      *cindex_param , *pisequence , *cisequence;
  double *prsequence , *crsequence;


  nb_sequence = seq.nb_sequence;

  identifier = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    identifier[i] = seq.identifier[i];
  }

  max_length = seq.max_length;
  cumul_length = seq.cumul_length;

  length = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    length[i] = seq.length[i];
  }

  length_distribution = new FrequencyDistribution(*(seq.length_distribution));

  if (seq.vertex_identifier) {
    vertex_identifier = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      vertex_identifier[i] = new int[length[i]];

      pidentifier = vertex_identifier[i];
      cidentifier = seq.vertex_identifier[i] + length[i] - 1;
      for (j = 0;j < length[i];j++) {
        *pidentifier++ = *cidentifier--;
      }
    }
  }

  else {
    vertex_identifier = NULL;
  }

  index_parameter_type = seq.index_parameter_type;

  if (seq.index_parameter_distribution) {
    index_parameter_distribution = new FrequencyDistribution(*(seq.index_parameter_distribution));
  }
  else {
    index_parameter_distribution = NULL;
  }

  if (seq.index_interval) {
    index_interval = new FrequencyDistribution(*(seq.index_interval));
  }
  else {
    index_interval = NULL;
  }

  if (seq.index_parameter) {
    index_parameter = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      blength = (index_parameter_type == POSITION ? length[i] + 1 : length[i]);
      index_parameter[i] = new int[blength];
      pindex_param = index_parameter[i];

      if (index_parameter_type == POSITION) {
        cindex_param = seq.index_parameter[i] + length[i];
        end_position = *cindex_param--;
        for (j = 0;j < length[i];j++) {
          *pindex_param++ = end_position - *cindex_param--;
        }
        *pindex_param = end_position;
      }

      else if (index_parameter_type == TIME) {
        cindex_param = seq.index_parameter[i] + length[i] - 1;
        for (j = 0;j < length[i];j++) {
          *pindex_param++ = index_parameter_distribution->nb_value - *cindex_param--;
        }
      }

      else {
        cindex_param = seq.index_parameter[i] + length[i] - 1;
        for (j = 0;j < length[i];j++) {
          *pindex_param++ = *cindex_param--;
        }
      }
    }
  }

  else {
    index_parameter = NULL;
  }

  nb_variable = seq.nb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = seq.type[i];
    min_value[i] = seq.min_value[i];
    max_value[i] = seq.max_value[i];

    if (seq.marginal_distribution[i]) {
      marginal_distribution[i] = new FrequencyDistribution(*(seq.marginal_distribution[i]));
    }
    else {
      marginal_distribution[i] = NULL;
    }

    if (seq.marginal_histogram[i]) {
      marginal_histogram[i] = new Histogram(*(seq.marginal_histogram[i]));
    }
    else {
      marginal_histogram[i] = NULL;
    }
  }

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      if ((type[j] != REAL_VALUE) && (type[j] != AUXILIARY)) {
        int_sequence[i][j] = new int[length[i]];
        real_sequence[i][j] = NULL;

        pisequence = int_sequence[i][j];
        cisequence = seq.int_sequence[i][j] + length[i] - 1;
        for (k = 0;k < length[i];k++) {
          *pisequence++ = *cisequence--;
        }
      }

      else {
        int_sequence[i][j] = NULL;
        real_sequence[i][j] = new double[length[i]];

        prsequence = real_sequence[i][j];
        crsequence = seq.real_sequence[i][j] + length[i] - 1;
        for (k = 0;k < length[i];k++) {
          *prsequence++ = *crsequence--;
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Sequences avec ajout d'une variable d'etat.
 *
 *  argument : reference sur un objet Sequences.
 *
 *--------------------------------------------------------------*/

void Sequences::add_state_variable(const Sequences &seq)

{
  register int i , j , k;
  int blength;


  nb_sequence = seq.nb_sequence;

  identifier = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    identifier[i] = seq.identifier[i];
  }

  max_length = seq.max_length;
  cumul_length = seq.cumul_length;

  length = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    length[i] = seq.length[i];
  }

  length_distribution = new FrequencyDistribution(*(seq.length_distribution));

  if (seq.vertex_identifier) {
    vertex_identifier = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      vertex_identifier[i] = new int[length[i]];
      for (j = 0;j < length[i];j++) {
        vertex_identifier[i][j] = seq.vertex_identifier[i][j];
      }
    }
  }

  else {
    vertex_identifier = NULL;
  }

  index_parameter_type = seq.index_parameter_type;

  if (seq.index_parameter_distribution) {
    index_parameter_distribution = new FrequencyDistribution(*(seq.index_parameter_distribution));
  }
  else {
    index_parameter_distribution = NULL;
  }

  if (seq.index_interval) {
    index_interval = new FrequencyDistribution(*(seq.index_interval));
  }
  else {
    index_interval = NULL;
  }

  if (seq.index_parameter) {
    index_parameter = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      blength = (index_parameter_type == POSITION ? length[i] + 1 : length[i]);
      index_parameter[i] = new int[blength];
      for (j = 0;j < blength;j++) {
        index_parameter[i][j] = seq.index_parameter[i][j];
      }
    }
  }

  else {
    index_parameter = NULL;
  }

  nb_variable = seq.nb_variable + 1;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  type[0] = STATE;
  min_value[0] = 0.;
  max_value[0] = 0.;
  marginal_distribution[0] = NULL;
  marginal_histogram[0] = NULL;

  for (i = 0;i < seq.nb_variable;i++) {
    type[i + 1] = (seq.type[i] == STATE ? INT_VALUE : seq.type[i]);
    min_value[i + 1] = seq.min_value[i];
    max_value[i + 1] = seq.max_value[i];

    if (seq.marginal_distribution[i]) {
      marginal_distribution[i + 1] = new FrequencyDistribution(*(seq.marginal_distribution[i]));
    }
    else {
      marginal_distribution[i + 1] = NULL;
    }

    if (seq.marginal_histogram[i]) {
      marginal_histogram[i + 1] = new Histogram(*(seq.marginal_histogram[i]));
    }
    else {
      marginal_histogram[i + 1] = NULL;
    }
  }

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];

    int_sequence[i][0] = new int[length[i]];
    real_sequence[i][0] = NULL;

    for (j = 0;j < length[i];j++) {
      int_sequence[i][0][j] = 0;
    }

    for (j = 0;j < seq.nb_variable;j++) {
      if (seq.type[j] != REAL_VALUE) {
        int_sequence[i][j + 1] = new int[length[i]];
        real_sequence[i][j + 1] = NULL;

        for (k = 0;k < length[i];k++) {
          int_sequence[i][j + 1][k] = seq.int_sequence[i][j][k];
        }
      }

      else {
        int_sequence[i][j + 1] = NULL;
        real_sequence[i][j + 1] = new double[length[i]];

        for (k = 0;k < length[i];k++) {
          real_sequence[i][j + 1][k] = seq.real_sequence[i][j][k];
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Sequences avec suppression du parametre d'index.
 *
 *  argument : reference sur un objet Sequences.
 *
 *--------------------------------------------------------------*/

void Sequences::remove_index_parameter(const Sequences &seq)

{
  register int i , j , k;


  nb_sequence = seq.nb_sequence;

  identifier = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    identifier[i] = seq.identifier[i];
  }

  max_length = seq.max_length;
  cumul_length = seq.cumul_length;

  length = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    length[i] = seq.length[i];
  }

  length_distribution = new FrequencyDistribution(*(seq.length_distribution));

  if (seq.vertex_identifier) {
    vertex_identifier = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      vertex_identifier[i] = new int[length[i]];
      for (j = 0;j < length[i];j++) {
        vertex_identifier[i][j] = seq.vertex_identifier[i][j];
      }
    }
  }

  else {
    vertex_identifier = NULL;
  }

  index_parameter_type = IMPLICIT_TYPE;
  index_parameter_distribution = NULL;
  index_interval = NULL;
  index_parameter = NULL;

  nb_variable = seq.nb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal_distribution = new FrequencyDistribution*[nb_variable];
  marginal_histogram = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = seq.type[i];
    min_value[i] = seq.min_value[i];
    max_value[i] = seq.max_value[i];

    if (seq.marginal_distribution[i]) {
      marginal_distribution[i] = new FrequencyDistribution(*(seq.marginal_distribution[i]));
    }
    else {
      marginal_distribution[i] = NULL;
    }

    if (seq.marginal_histogram[i]) {
      marginal_histogram[i] = new Histogram(*(seq.marginal_histogram[i]));
    }
    else {
      marginal_histogram[i] = NULL;
    }
  }

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      if ((type[j] != REAL_VALUE) && (type[j] != AUXILIARY)) {
        int_sequence[i][j] = new int[length[i]];
        real_sequence[i][j] = NULL;

        for (k = 0;k < length[i];k++) {
          int_sequence[i][j][k] = seq.int_sequence[i][j][k];
        }
      }

      else {
        int_sequence[i][j] = NULL;
        real_sequence[i][j] = new double[length[i]];

        for (k = 0;k < length[i];k++) {
          real_sequence[i][j][k] = seq.real_sequence[i][j][k];
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Sequences avec transformation du parametre d'index implicite
 *  en parametre d'index explicite.
 *
 *  argument : reference sur un objet Sequences.
 *
 *--------------------------------------------------------------*/

void Sequences::explicit_index_parameter(const Sequences &seq)

{
  register int i , j;


  Sequences::copy(seq);

  index_parameter_type = TIME;

  index_parameter = new int*[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    index_parameter[i] = new int[length[i]];
    for (j = 0;j < length[i];j++) {
      index_parameter[i][j] = j;
    }
  }

  build_index_parameter_frequency_distribution();
  index_interval_computation();
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe Sequences.
 *
 *  arguments : reference sur un objet Sequences, type de transformation
 *             ('r' : inversion du sens de parcours, 'a' : addition d'une variable,
 *              'm' : suppression du parametre d'index, 'e' : parametre d'index rendu explicite).
 *
 *--------------------------------------------------------------*/

Sequences::Sequences(const Sequences &seq , char transform , int param)

{
  switch (transform) {
  case 'r' :
    Sequences::reverse(seq);
    break;
  case 'a' :
    Sequences::add_state_variable(seq);
    break;
  case 'm' :
    Sequences::remove_index_parameter(seq);
    break;
  case 'e' :
    Sequences::explicit_index_parameter(seq);
    break;
  default :
    Sequences::copy(seq);
    break;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Sequences.
 *
 *--------------------------------------------------------------*/

void Sequences::remove()

{
  register int i , j;


  delete [] identifier;

  delete [] length;
  delete length_distribution;

  if (vertex_identifier) {
    for (i = 0;i < nb_sequence;i++) {
      delete [] vertex_identifier[i];
    }
    delete [] vertex_identifier;
  }

  delete index_parameter_distribution;
  delete index_interval;

  if (index_parameter) {
    for (i = 0;i < nb_sequence;i++) {
      delete [] index_parameter[i];
    }
    delete [] index_parameter;
  }

  delete [] type;
  delete [] min_value;
  delete [] max_value;

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

  if (int_sequence) {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < nb_variable;j++) {
        delete [] int_sequence[i][j];
      }
      delete [] int_sequence[i];
    }
    delete [] int_sequence;
  }

  if (real_sequence) {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < nb_variable;j++) {
        delete [] real_sequence[i][j];
      }
      delete [] real_sequence[i];
    }
    delete [] real_sequence;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Sequences.
 *
 *--------------------------------------------------------------*/

Sequences::~Sequences()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Sequences.
 *
 *  argument : reference sur un objet Sequences.
 *
 *--------------------------------------------------------------*/

Sequences& Sequences::operator=(const Sequences &seq)

{
  if (&seq != this) {
    remove();
    copy(seq);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de la loi marginale empirique pour une variable entiere.
 *
 *  arguments : reference sur un objet StatError, variable.
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData* Sequences::extract(StatError &error , int variable) const

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
 *  Construction d'un object Vectors a partir d'un objet Sequences.
 *
 *  argument : flag variable d'index.
 *
 *--------------------------------------------------------------*/

Vectors* Sequences::build_vectors(bool index_variable) const

{
  register int i , j , k , m;
  int offset , *itype , **int_vector;
  double **real_vector;
  Vectors *vec;


  if (index_parameter) {
    index_variable = true;
  }

  switch (index_variable) {
  case false :
    offset = 0;
    break;
  case true :
    offset = 1;
    break;
  }

  itype = new int[nb_variable + offset];

  if (index_variable) {
    itype[0] = INT_VALUE;
  }

  for (i = 0;i < nb_variable;i++) {
    switch (type[i]) {
    case STATE :
      itype[i + offset] = INT_VALUE;
      break;
    case AUXILIARY :
      itype[i + offset] = REAL_VALUE;
      break;
    default :
      itype[i + offset] = type[i];
      break;
    }
  }

  int_vector = new int*[cumul_length];
  for (i = 0;i < cumul_length;i++) {
    int_vector[i] = new int[nb_variable + offset];
  }

  real_vector = new double*[cumul_length];
  for (i = 0;i < cumul_length;i++) {
    real_vector[i] = new double[nb_variable + offset];
  }

  i = 0;
  for (j = 0;j < nb_sequence;j++) {
    for (k = 0;k < length[j];k++) {
      if (index_variable) {
        if (index_parameter) {
          int_vector[i][0] = index_parameter[j][k];
        }
        else {
          int_vector[i][0] = k;
        }
      }

      for (m = 0;m < nb_variable;m++) {
        if ((type[m] != REAL_VALUE) && (type[m] != AUXILIARY)) {
          int_vector[i][m + offset] = int_sequence[j][m][k];
        }
        else {
          real_vector[i][m + offset] = real_sequence[j][m][k];
        }
      }

      i++;
    }
  }

  vec = new Vectors(cumul_length , NULL , nb_variable + offset , itype , int_vector , real_vector);
  delete [] itype;

  for (i = 0;i < cumul_length;i++) {
    delete [] int_vector[i];
  }
  delete [] int_vector;

  for (i = 0;i < cumul_length;i++) {
    delete [] real_vector[i];
  }
  delete [] real_vector;

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de mesures globales (longueur, temps avant la 1ere occurrence
 *  d'une valeur, nombre de series/d'occurrences d'une valeur) par sequence.
 *
 *  arguments : reference sur un objet StatError, type, variable, valeur.
 *
 *--------------------------------------------------------------*/

Vectors* Sequences::extract_vectors(StatError &error , int feature_type ,
                                    int variable , int value) const

{
  bool status = true;
  register int i , j;
  int begin_run , count , itype[1] , **int_vector;
  double **real_vector;
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

      if ((feature_type == SEQUENCE_CUMUL) || (feature_type == SEQUENCE_MEAN)) {
        if ((type[variable] != INT_VALUE) && (type[variable] != STATE) &&
            (type[variable] != REAL_VALUE)) {
          status = false;
          ostringstream error_message , correction_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE]
                             << " or " << STAT_variable_word[REAL_VALUE];
          error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
        }
      }

//      else if ((feature_type == FIRST_OCCURRENCE) || (feature_type == NB_RUN) ||
//               (feature_type == NB_OCCURRENCE)) {
      else {
        if ((type[variable] != INT_VALUE) && (type[variable] != STATE)) {
          status = false;
          ostringstream error_message , correction_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
          error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
        }

        if ((value < min_value[variable]) || (value > max_value[variable]) ||
            ((marginal_distribution[variable]) && (marginal_distribution[variable]->frequency[value] == 0))) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VALUE] << " " << value << " "
                        << STAT_error[STATR_NOT_PRESENT];
          error.update((error_message.str()).c_str());
        }
      }
    }
  }

  if (status) {
    switch (feature_type) {
    case SEQUENCE_CUMUL :
      itype[0] = type[variable];
      break;
    case SEQUENCE_MEAN :
      itype[0] = REAL_VALUE;
      break;
    default :
      itype[0] = INT_VALUE;
      break;
    }

    int_vector = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      int_vector[i] = new int[1];
    }

    real_vector = new double*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      real_vector[i] = new double[1];
    }

    switch (feature_type) {

    case LENGTH : {
      if (index_parameter_type == POSITION) {
        for (i = 0;i < nb_sequence;i++) {
          int_vector[i][0] = index_parameter[i][length[i]];
        }
      }

      else if ((index_parameter_type == TIME) && (index_interval->variance > 0.)) {
        for (i = 0;i < nb_sequence;i++) {
          int_vector[i][0] = index_parameter[i][length[i] - 1];
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          int_vector[i][0] = length[i];
        }
      }
      break;
    }

    case SEQUENCE_CUMUL : {
      if (type[variable] != REAL_VALUE) {
        for (i = 0;i < nb_sequence;i++) {
          int_vector[i][0] = 0;
          for (j = 0;j < length[i];j++) {
            int_vector[i][0] += int_sequence[i][variable][j];
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          real_vector[i][0] = 0.;
          for (j = 0;j < length[i];j++) {
            real_vector[i][0] += real_sequence[i][variable][j];
          }
        }
      }
      break;
    }

    case SEQUENCE_MEAN : {
      if (type[variable] != REAL_VALUE) {
        for (i = 0;i < nb_sequence;i++) {
          real_vector[i][0] = 0.;
          for (j = 0;j < length[i];j++) {
            real_vector[i][0] += int_sequence[i][variable][j];
          }
          real_vector[i][0] /= length[i];
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          real_vector[i][0] = 0.;
          for (j = 0;j < length[i];j++) {
            real_vector[i][0] += real_sequence[i][variable][j];
          }
          real_vector[i][0] /= length[i];
        }
      }
      break;
    }

    case FIRST_OCCURRENCE : {
      if (index_parameter_type != IMPLICIT_TYPE) {
        for (i = 0;i < nb_sequence;i++) {
          int_vector[i][0] = -1;
          for (j = 0;j < length[i];j++) {
            if (int_sequence[i][variable][j] == value) {
              int_vector[i][0] = index_parameter[i][j];
              break;
            }
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          int_vector[i][0] = -1;
          for (j = 0;j < length[i];j++) {
            if (int_sequence[i][variable][j] == value) {
              int_vector[i][0] = j;
              break;
            }
          }
        }
      }
      break;
    }

    case SOJOURN_TIME : {
      if ((index_parameter_type == TIME) && (index_interval->variance > 0.)) {  // pour les suivis de croissance manguier
        for (i = 0;i < nb_sequence;i++) {
          int_vector[i][0] = -1;
          if (int_sequence[i][variable][0] == value) {
            begin_run = 0;
          }

          for (j = 0;j < length[i] - 1;j++) {
            if (int_sequence[i][variable][j + 1] != int_sequence[i][variable][j]) {
              if (int_sequence[i][variable][j + 1] == value) {
                begin_run = index_parameter[i][j + 1];
              }
              else if (int_sequence[i][variable][j] == value) {
                int_vector[i][0] = index_parameter[i][j + 1] - begin_run - 1;
                break;
              }
            }
          }

          if ((j == length[i] - 1) && (int_sequence[i][variable][length[i] - 1] == value)) {
            int_vector[i][0] = index_parameter[i][j] - begin_run;
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
//          int_vector[i][0] = -1;
          int_vector[i][0] = 0;
          if (int_sequence[i][variable][0] == value) {
            begin_run = 0;
          }

          for (j = 0;j < length[i] - 1;j++) {
            if (int_sequence[i][variable][j + 1] != int_sequence[i][variable][j]) {
              if (int_sequence[i][variable][j + 1] == value) {
                begin_run = j + 1;
              }
              else if (int_sequence[i][variable][j] == value) {
                int_vector[i][0] = j + 1 - begin_run;
                break;
              }
            }
          }

          if ((j == length[i] - 1) && (int_sequence[i][variable][length[i] - 1] == value)) {
            int_vector[i][0] = length[i] - begin_run;
          }
        }
      }
      break;
    }

    case NB_RUN : {
      for (i = 0;i < nb_sequence;i++) {
        count = 0;
        if (int_sequence[i][variable][0] == value) {
          count++;
        }
        for (j = 1;j < length[i];j++) {
          if ((int_sequence[i][variable][j] != int_sequence[i][variable][j - 1]) &&
              (int_sequence[i][variable][j] == value)) {
            count++;
          }
        }

        int_vector[i][0] = count;
      }
      break;
    }

    case NB_OCCURRENCE : {
      for (i = 0;i < nb_sequence;i++) {
        count = 0;
        for (j = 0;j < length[i];j++) {
          if (int_sequence[i][variable][j] == value) {
            count++;
          }
        }

        int_vector[i][0] = count;
      }
      break;
    }
    }

    vec = new Vectors(nb_sequence , identifier , 1 , itype , int_vector , real_vector);

    for (i = 0;i < nb_sequence;i++) {
     delete [] int_vector[i];
    }
    delete [] int_vector;

    for (i = 0;i < nb_sequence;i++) {
      delete [] real_vector[i];
    }
    delete [] real_vector;
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet MarkovianSequences a partir d'un objet Sequences.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

MarkovianSequences* Sequences::markovian_sequences(StatError &error) const

{
  bool status = true;
  register int i;
  MarkovianSequences *seq;


  seq = NULL;
  error.init();

//  if (((index_parameter_type == TIME) && (index_interval->variance > 0.)) ||
//      (index_parameter_type == POSITION)) {
  if (index_parameter_type == POSITION) {
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

    else if (max_value[i] == min_value[i]) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_NB_VALUE];
      error.update((error_message.str()).c_str());
    }

    if ((type[i] == INT_VALUE) || (type[i] == STATE)) {
      if (min_value[i] < 0) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                      << STAT_error[STATR_POSITIVE_MIN_VALUE];
        error.update((error_message.str()).c_str());
      }

      if (!marginal_distribution[i]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                      << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {
    seq = new MarkovianSequences(*this);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Verification du caractere (strictement) croissant des parametres d'index.
 *
 *  arguments : reference sur un objet StatError, flag croissance strict ou non,
 *              label de l'objet.
 *
 *--------------------------------------------------------------*/

bool Sequences::increasing_index_parameter_checking(StatError &error , bool strict ,
                                                    const char *pattern_label) const

{
  bool status = true;
  register int i , j;


  for (i = 0;i < nb_sequence;i++) {
    for (j = 1;j < (index_parameter_type == POSITION ? length[i] + 1 : length[i]);j++) {
      if ((((!strict) || (j == length[i])) && (index_parameter[i][j] < index_parameter[i][j - 1])) ||
          ((strict) && (j < length[i]) && (index_parameter[i][j] <= index_parameter[i][j - 1]))) {
        status = false;
        ostringstream error_message;
        error_message << pattern_label << " " << i + 1 << ": "
                      << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME] : SEQ_label[SEQL_POSITION]) << " "
                      << index_parameter[i][j] << " " << STAT_error[STATR_NOT_ALLOWED];
        error.update((error_message.str()).c_str());
      }
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Verification du caractere (strictement) croissant des sequences
 *  pour une variable entiere.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              flag croissance strict ou non, labels de l'objet et de la variable.
 *
 *--------------------------------------------------------------*/

bool Sequences::increasing_sequence_checking(StatError &error , int variable , bool strict ,
                                             const char *pattern_label , const char *variable_label) const

{
  bool status = true;
  register int i , j;


  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 1;j < length[i];j++) {
        if (((!strict) && (int_sequence[i][variable][j] < int_sequence[i][variable][j - 1])) ||
            ((strict) && (int_sequence[i][variable][j] <= int_sequence[i][variable][j - 1]))) {
          status = false;
          ostringstream error_message;
          error_message << pattern_label << " " << i + 1 << ": " << STAT_label[STATL_VARIABLE] << " "
                        << variable + 1 << ": " << variable_label << " "
                        << int_sequence[i][variable][j] << " " << STAT_error[STATR_NOT_ALLOWED];
          error.update((error_message.str()).c_str());
        }
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 1;j < length[i];j++) {
        if (((!strict) && (real_sequence[i][variable][j] < real_sequence[i][variable][j - 1])) ||
            ((strict) && (real_sequence[i][variable][j] <= real_sequence[i][variable][j - 1]))) {
          status = false;
          ostringstream error_message;
          error_message << pattern_label << " " << i + 1 << ": " << STAT_label[STATL_VARIABLE] << " "
                        << variable + 1 << ": " << variable_label << " "
                        << real_sequence[i][variable][j] << " " << STAT_error[STATR_NOT_ALLOWED];
          error.update((error_message.str()).c_str());
        }
      }
    }
    break;
  }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Verification d'un objet Sequences.
 *
 *  arguments : reference sur un objet StatError, label.
 *
 *--------------------------------------------------------------*/

bool Sequences::check(StatError &error , const char *pattern_label)

{
  bool status = true , lstatus;


  error.init();

  if (nb_variable > SEQUENCE_NB_VARIABLE) {
    status = false;
    error.update(STAT_error[STATR_NB_VARIABLE]);
  }

  if (max_length == 1) {
    status = false;
    error.update(SEQ_parsing[SEQP_MAX_SEQUENCE_LENGTH]);
  }

  lstatus = identifier_checking(error , nb_sequence , identifier);
  if (!lstatus) {
    status = false;
  }

  if (index_parameter_type != IMPLICIT_TYPE) {
    lstatus = increasing_index_parameter_checking(error , (index_parameter_type == POSITION ? false : true) ,
                                                  pattern_label);

    if (!lstatus) {
      status = false;
    }
  }

  if (status) {
    if (index_parameter) {
      build_index_parameter_frequency_distribution();
    }
//    if ((index_parameter_type == TIME) || ((index_parameter_type == POSITION) &&
//         (type[0] != NB_INTERNODE))) {
    if ((index_parameter_type == TIME) || (index_parameter_type == POSITION)) {
      index_interval_computation();
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'un objet TimeEvents a partir d'un objet Sequences.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              nombre d'evenements, dates de debut et de fin,
 *              dates precedente et suivante.
 *
 *--------------------------------------------------------------*/

TimeEvents* Sequences::extract_time_events(StatError &error , int variable ,
                                           int begin_date , int end_date ,
                                           int previous_date , int next_date) const

{
  bool status = true , lstatus;
  register int i , j;
  int nb_element , previous , begin , end , next , *time , *nb_event , *pdate;
  TimeEvents *timev;


  timev = NULL;
  error.init();

  if (index_parameter_type != TIME) {
    status = false;
    error.correction_update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE] , SEQ_index_parameter_word[TIME]);
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
      lstatus = increasing_sequence_checking(error , variable , false , SEQ_label[SEQL_SEQUENCE] ,
                                             STAT_label[STATL_VALUE]);
      if (!lstatus) {
        status = false;
      }
    }
  }

  if (begin_date >= end_date) {
    status = false;
    error.update(SEQ_error[SEQR_DATE_ORDER]);
  }

  if (previous_date != I_DEFAULT) {
    if (previous_date > begin_date) {
      status = false;
      error.update(SEQ_error[SEQR_DATE_ORDER]);
    }
  }
  else {
    previous_date = begin_date;
  }

  if (next_date != I_DEFAULT) {
    if (next_date < end_date) {
      status = false;
      error.update(SEQ_error[SEQR_DATE_ORDER]);
    }
  }
  else {
    next_date = end_date;
  }

  if (status) {
    time = new int[nb_sequence];
    nb_event = new int[nb_sequence];
    nb_element = 0;

    for (i = 0;i < nb_sequence;i++) {
      pdate = index_parameter[i];
      previous = I_DEFAULT;
      begin = I_DEFAULT;
      end = I_DEFAULT;
      next = I_DEFAULT;

      for (j = 0;j < length[i];j++) {
        if (*pdate == previous_date) {
          previous = j;
        }
        if (*pdate == begin_date) {
          begin = j;
        }
        if (*pdate == end_date) {
          end = j;
        }
        if (*pdate == next_date) {
          next = j;
          break;
        }
        pdate++;
      }

      if ((previous != I_DEFAULT) && (begin != I_DEFAULT) && (end != I_DEFAULT) &&
          (next != I_DEFAULT) && ((previous == begin) || ((previous < begin) &&
            (int_sequence[i][variable][previous] < int_sequence[i][variable][begin]))) && ((end == next) ||
           ((end < next) && (int_sequence[i][variable][end] < int_sequence[i][variable][next])))) {
        time[nb_element] = end_date - begin_date;
        nb_event[nb_element++] = int_sequence[i][variable][end] - int_sequence[i][variable][begin];
      }
    }

    if (nb_element == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    else {
      timev = new TimeEvents(nb_element , time , nb_event);
    }

    delete [] time;
    delete [] nb_event;
  }

  return timev;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'un objet RenewalData a partir d'un objet Sequences.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              nombre d'evenements, index de debut et de fin.
 *
 *--------------------------------------------------------------*/

RenewalData* Sequences::extract_renewal_data(StatError &error , int variable ,
                                             int begin_index_parameter , int end_index_parameter) const

{
  bool status = true , lstatus;
  register int i , j;
  int nb_element , index , *ptime , *pnb_event , *pisequence , *cisequence;
  RenewalData *timev;


  timev = NULL;
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
      lstatus = increasing_sequence_checking(error , variable , false , SEQ_label[SEQL_SEQUENCE] ,
                                             STAT_label[STATL_VALUE]);
      if (!lstatus) {
        status = false;
      }
    }
  }

  if ((begin_index_parameter < 0) || (begin_index_parameter + 1 >= max_length) ||
      (begin_index_parameter > end_index_parameter)) {
    status = false;
    error.update(SEQ_error[SEQR_BEGIN_INDEX_PARAMETER]);
  }
  if ((end_index_parameter < 0) || (end_index_parameter + 1 >= max_length) ||
      (end_index_parameter < begin_index_parameter)) {
    status = false;
    error.update(SEQ_error[SEQR_END_INDEX_PARAMETER]);
  }

  if (status) {
    timev = new RenewalData(nb_sequence , end_index_parameter + 1 - begin_index_parameter);

    ptime = new int[nb_sequence];
    pnb_event = new int[nb_sequence];

    nb_element = 0;
    for (i = 0;i < nb_sequence;i++) {
      if (end_index_parameter + 1 < length[i]) {
        *ptime++ = end_index_parameter + 1 - begin_index_parameter;
        *pnb_event++ = int_sequence[i][variable][end_index_parameter + 1] -
                       int_sequence[i][variable][begin_index_parameter];

        timev->length[nb_element] = end_index_parameter + 1 - begin_index_parameter;
        timev->sequence[nb_element] = new int[timev->length[nb_element]];

        pisequence = timev->sequence[nb_element++];
        cisequence = int_sequence[i][variable] + begin_index_parameter;
        index = begin_index_parameter;
        for (j = begin_index_parameter + 1;j <= end_index_parameter + 1;j++) {
          *pisequence = *(cisequence + 1) - *cisequence;
          cisequence++;

          if (*pisequence > 0) {
            if (index == begin_index_parameter) {
              (timev->forward->frequency[j - index])++;
             }
             else {
               (timev->within->frequency[j - index])++;
             }
             index = j;
          }
          pisequence++;
        }

        if (index > begin_index_parameter) {
          (timev->backward->frequency[end_index_parameter + 1 - index])++;
        }
      }
    }

    // construction des echantillons {temps, nombre d'evenements, frequence} et
    // des lois empiriques du temps d'observation et du nombre d'evenements

    ptime -= nb_element;
    pnb_event -= nb_element;

    timev->build(nb_element , ptime , pnb_event);
    delete [] ptime;
    delete [] pnb_event;

    // extraction des caracteristiques des lois empiriques des intervalles de temps entre 2 evenements,
    // des intervalles de temps entre 2 evenements compris entre les 2 dates d'observation,
    // des intervalles de temps apres le dernier evenement, des intervalles de temps residuel

    timev->within->nb_value_computation();
    timev->within->offset_computation();
    timev->within->nb_element_computation();
    timev->within->max_computation();
    timev->within->mean_computation();
    timev->within->variance_computation();

    timev->backward->nb_value_computation();
    timev->backward->offset_computation();
    timev->backward->nb_element_computation();
    timev->backward->max_computation();
    timev->backward->mean_computation();
    timev->backward->variance_computation();

    timev->forward->nb_value_computation();
    timev->forward->offset_computation();
    timev->forward->nb_element_computation();
    timev->forward->max_computation();
    timev->forward->mean_computation();
    timev->forward->variance_computation();

    timev->build_index_event(1);

    if ((timev->backward->nb_element == 0) && (timev->forward->nb_element == 0)) {
      delete timev;
      timev = NULL;
      error.update(SEQ_error[SEQR_BOTH_END_CENSORED_INTERVAL]);
    }
  }

  return timev;
}


};  // namespace sequence_analysis
