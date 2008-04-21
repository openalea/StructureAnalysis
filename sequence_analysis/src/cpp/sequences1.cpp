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
#include "tool/rw_locale.h"
#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"
#include "stat_tool/vectors.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "renewal.h"
#include "sequences.h"
#include "tops.h"
#include "sequence_label.h"

using namespace std;


extern bool identifier_checking(Format_error &error , int nb_individual , int *individual_identifier);
extern bool selected_identifier_checking(Format_error &error , int nb_individual , int *identifier ,
                                         int nb_selected_individual , int *selected_identifier ,
                                         const char *data_label);
extern int* identifier_select(int nb_individual , int *identifier , int nb_selected_individual ,
                              int *selected_identifier , bool keep);
extern int* select_variable(int nb_variable , int nb_selected_variable ,
                            int *selected_variable , bool keep);

extern int column_width(int value);
extern int column_width(int nb_value , const double *value , double scale = 1.);

extern double standard_normal_value_computation(double critical_probability);
extern double t_value_computation(bool one_side , int df , double critical_probability);


/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Sequences.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences()

{
  nb_sequence = 0;
  identifier = 0;

  max_length = 0;
  cumul_length = 0;
  length = 0;
  hlength = 0;

  index_parameter_type = IMPLICIT_TYPE;
  hindex_parameter = 0;
  index_interval = 0;
  index_parameter = 0;

  nb_variable = 0;

  type = 0;
  min_value = 0;
  max_value = 0;
  marginal = 0;

  int_sequence = 0;
  real_sequence = 0;
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
  hlength = 0;

  index_parameter_type = IMPLICIT_TYPE;
  hindex_parameter = 0;
  index_interval = 0;
  index_parameter = 0;

  nb_variable = inb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = INT_VALUE;
    min_value[i] = 0.;
    max_value[i] = 0.;
    marginal[i] = 0;
  }

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      int_sequence[i][j] = 0;
      real_sequence[i][j] = 0;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Initialisation d'un objet Sequences.
 *
 *  arguments : nombre de sequences, identificateurs des sequences,
 *              longueurs des sequences, type de parametre d'index, nombre de variables,
 *              type de chaque variable, flag initialisation.
 *
 *--------------------------------------------------------------*/

void Sequences::init(int inb_sequence , int *iidentifier , int *ilength ,
                     int iindex_parameter_type , int inb_variable , int *itype ,
                     bool init_flag)

{
  register int i , j , k;
  int blength , *pisequence;
  double *prsequence;


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
  build_length_histogram();

  index_parameter_type = iindex_parameter_type;
  hindex_parameter = 0;
  index_interval = 0;

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
    index_parameter = 0;
  }

  nb_variable = inb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = itype[i];
    min_value[i] = 0.;
    max_value[i] = 0.;
    marginal[i] = 0;
  }

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      if ((type[j] != REAL_VALUE) && (type[j] != AUXILIARY)) {
        int_sequence[i][j] = new int[length[i]];
        real_sequence[i][j] = 0;

        if (init_flag) {
          pisequence = int_sequence[i][j];
          for (k = 0;k < length[i];k++) {
            *pisequence++ = 0;
          }
        }
      }

      else {
        int_sequence[i][j] = 0;
        real_sequence[i][j] = new double[length[i]];

        if (init_flag) {
          prsequence = real_sequence[i][j];
          for (k = 0;k < length[i];k++) {
            *prsequence++ = 0.;
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
 *  arguments : nombre de sequences, identificateurs des sequences,
 *              longueurs des sequences, type de parametre d'index (INT_VALUE/NB_INTERNODE),
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

  init(inb_sequence , iidentifier , ilength , iindex_parameter_type , inb_variable ,
       btype , false);
  delete [] btype;

//  if (index_parameter_type != IMPLICIT_TYPE) {
  if (index_parameter) {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < (index_parameter_type == POSITION ? length[i] + 1 : length[i]);j++) {
        index_parameter[i][j] = iint_sequence[i][0][j];
      }
    }

    build_index_parameter_histogram();

    if ((index_parameter_type == TIME) || ((index_parameter_type == POSITION) &&
        (type[0] != NB_INTERNODE))) {
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
    build_marginal_histogram(i);
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
  double *prsequence , *crsequence;


  itype = new int[inb_variable];
  for (i = 0;i < inb_variable;i++) {
    itype[i] = REAL_VALUE;
  }

  init(inb_sequence , iidentifier , ilength , IMPLICIT_TYPE , inb_variable ,
       itype , false);
  delete [] itype;

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < nb_variable;j++) {
      prsequence = real_sequence[i][j];
      crsequence = ireal_sequence[i][j];
      for (k = 0;k < length[i];k++) {
        *prsequence++ = *crsequence++;
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value_computation(i);
    max_value_computation(i);
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
  int *pisequence;


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
  build_length_histogram();

  index_parameter_type = IMPLICIT_TYPE;
  hindex_parameter = 0;
  index_interval = 0;
  index_parameter = 0;

  nb_variable = inb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = INT_VALUE;
    min_value[i] = 0.;
    max_value[i] = 0.;
    marginal[i] = 0;
  }

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      int_sequence[i][j] = new int[length[i]];
      real_sequence[i][j] = 0;

      if (init_flag) {
        pisequence = int_sequence[i][j];
        for (k = 0;k < length[i];k++) {
          *pisequence++ = 0;
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Sequences.
 *
 *  arguments : histogramme des longueurs des sequences, nombre de variables,
 *              flag initialisation.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences(const Histogram &ihlength , int inb_variable , bool init_flag)

{
  register int i , j , k;
  int *plength , *pisequence;


  nb_sequence = ihlength.nb_element;

  identifier = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    identifier[i] = i + 1;
  }

  length = new int[nb_sequence];
  plength = length;
  for (i = ihlength.offset;i < ihlength.nb_value;i++) {
    for (j = 0;j < ihlength.frequency[i];j++) {
      *plength++ = i;
    }
  }

  max_length = ihlength.nb_value - 1;
  cumul_length_computation();
  hlength = new Histogram(ihlength);

  index_parameter_type = IMPLICIT_TYPE;
  hindex_parameter = 0;
  index_interval = 0;
  index_parameter = 0;

  nb_variable = inb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = INT_VALUE;
    min_value[i] = 0.;
    max_value[i] = 0.;
    marginal[i] = 0;
  }

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      int_sequence[i][j] = new int[length[i]];
      real_sequence[i][j] = 0;

      if (init_flag) {
        pisequence = int_sequence[i][j];
        for (k = 0;k < length[i];k++) {
          *pisequence++ = 0;
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Sequences a partir d'un objet Renewal_data.
 *
 *  argument : reference sur un objet Renewal_data.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences(const Renewal_data &timev)

{
  register int i , j;
  int *pisequence , *cisequence;


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
  build_length_histogram();

  index_parameter_type = IMPLICIT_TYPE;
  hindex_parameter = 0;
  index_interval = 0;
  index_parameter = 0;

  nb_variable = 1;

  type = new int[nb_variable];
  type[0] = INT_VALUE;

  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal = new Histogram*[nb_variable];

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];

    int_sequence[i][0] = new int[length[i]];
    real_sequence[i][0] = 0;

    pisequence = int_sequence[i][0];
    cisequence = timev.sequence[i];
    for (j = 0;j < length[i];j++) {
      *pisequence++ = *cisequence++;
    }
  }

  min_value_computation(0);
  max_value_computation(0);
  build_marginal_histogram(0);
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
  int blength , *pindex_param , *cindex_param , *pisequence , *cisequence;
  double *prsequence , *crsequence;


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
  build_length_histogram();

  index_parameter_type = seq.index_parameter_type;

  nb_variable = seq.nb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = seq.type[i];
    marginal[i] = 0;
  }

//  if (index_parameter_type != IMPLICIT_TYPE) {
  if (seq.index_parameter) {
    index_parameter = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      blength = (index_parameter_type == POSITION ? length[i] + 1 : length[i]);
      index_parameter[i] = new int[blength];

      pindex_param = index_parameter[i];
      cindex_param = seq.index_parameter[index[i]];
      for (j = 0;j < blength;j++) {
        *pindex_param++ = *cindex_param++;
      }
    }

    build_index_parameter_histogram();
    index_interval_computation();
  }

  else {
    hindex_parameter = 0;
    index_interval = 0;
    index_parameter = 0;
  }

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      if ((type[j] != REAL_VALUE) && (type[j] != AUXILIARY)) {
        int_sequence[i][j] = new int[length[i]];
        real_sequence[i][j] = 0;

        pisequence = int_sequence[i][j];
        cisequence = seq.int_sequence[index[i]][j];
        for (k = 0;k < length[i];k++) {
          *pisequence++ = *cisequence++;
        }
      }

      else {
        int_sequence[i][j] = 0;
        real_sequence[i][j] = new double[length[i]];

        prsequence = real_sequence[i][j];
        crsequence = seq.real_sequence[index[i]][j];
        for (k = 0;k < length[i];k++) {
          *prsequence++ = *crsequence++;
        }
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
 *  arguments : reference sur un objet Sequences, flags sur l'ajout de variables
 *              de moyennes de segments.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences(const Sequences &seq , bool *segment_mean)

{
  register int i , j , k , m;
  int *pisequence , *cisequence;
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

  hlength = new Histogram(*(seq.hlength));

  index_parameter_type = seq.index_parameter_type;

  if (seq.hindex_parameter) {
    hindex_parameter = new Histogram(*(seq.hindex_parameter));
  }
  else {
    hindex_parameter = 0;
  }

  if (seq.index_interval) {
    index_interval = new Histogram(*(seq.index_interval));
  }
  else {
    index_interval = 0;
  }

  if (seq.index_parameter) {
    index_parameter = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      index_parameter[i] = new int[length[i]];
      for (j = 0;j < length[i];j++) {
        index_parameter[i][j] = seq.index_parameter[i][j];
      }
    }
  }
  else {
    index_parameter = 0;
  }

  nb_variable = seq.nb_variable;
  for (i = 0;i < seq.nb_variable;i++) {
    if (segment_mean[i]) {
      nb_variable++;
    }
  }

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal = new Histogram*[nb_variable];

  i = 0;
  for (j = 0;j < seq.nb_variable;j++) {
    type[i] = seq.type[j];
    min_value[i] = seq.min_value[j];
    max_value[i] = seq.max_value[j];

    if (seq.marginal[j]) {
      marginal[i] = new Histogram(*(seq.marginal[j]));
    }
    else {
      marginal[i] = 0;
    }
    i++;

    if (segment_mean[j]) {
      type[i] = AUXILIARY;
      min_value[i] = 0.;
      max_value[i] = 0.;
      marginal[i] = 0;
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
        real_sequence[i][j] = 0;

        pisequence = int_sequence[i][j];
        cisequence = seq.int_sequence[i][k];
        for (m = 0;m < length[i];m++) {
          *pisequence++ = *cisequence++;
        }
      }

      else {
        int_sequence[i][j] = 0;
        real_sequence[i][j] = new double[length[i]];

        prsequence = real_sequence[i][j];
        crsequence = seq.real_sequence[i][k];
        for (m = 0;m < length[i];m++) {
          *prsequence++ = *crsequence++;
        }
      }

      j++;

      if (segment_mean[k]) {
        int_sequence[i][j] = 0;
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
 *  arguments : reference sur un objet Sequences, flag inversion.
 *
 *--------------------------------------------------------------*/

void Sequences::copy(const Sequences &seq , bool reverse_flag)

{
  register int i , j , k;
  int blength , end_position , *pindex_param , *cindex_param , *pisequence , *cisequence;
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

  hlength = new Histogram(*(seq.hlength));

  index_parameter_type = seq.index_parameter_type;

  if (seq.hindex_parameter) {
    hindex_parameter = new Histogram(*(seq.hindex_parameter));
  }
  else {
    hindex_parameter = 0;
  }

  if (seq.index_interval) {
    index_interval = new Histogram(*(seq.index_interval));
  }
  else {
    index_interval = 0;
  }

  if (seq.index_parameter) {
    index_parameter = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      blength = (index_parameter_type == POSITION ? length[i] + 1 : length[i]);
      index_parameter[i] = new int[blength];
      pindex_param = index_parameter[i];

      switch (reverse_flag) {

      case false : {
        cindex_param = seq.index_parameter[i];
        for (j = 0;j < blength;j++) {
          *pindex_param++ = *cindex_param++;
        }
        break;
      }

      case true : {
        if (index_parameter_type == POSITION) {
          cindex_param = seq.index_parameter[i] + length[i];
          end_position = *cindex_param--;
          for (j = 0;j < length[i];j++) {
            *pindex_param++ = end_position - *cindex_param--;
          }
          *pindex_param = end_position;
        }

        else {
          cindex_param = seq.index_parameter[i] + length[i] - 1;
          for (j = 0;j < length[i];j++) {
            *pindex_param++ = *cindex_param--;
          }
        }
        break;
      }
      }
    }
  }

  else {
    index_parameter = 0;
  }

  nb_variable = seq.nb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = seq.type[i];
    min_value[i] = seq.min_value[i];
    max_value[i] = seq.max_value[i];

    if (seq.marginal[i]) {
      marginal[i] = new Histogram(*(seq.marginal[i]));
    }
    else {
      marginal[i] = 0;
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
        real_sequence[i][j] = 0;

        pisequence = int_sequence[i][j];

        switch (reverse_flag) {

        case false : {
          cisequence = seq.int_sequence[i][j];
          for (k = 0;k < length[i];k++) {
            *pisequence++ = *cisequence++;
          }
          break;
        }

        case true : {
          cisequence = seq.int_sequence[i][j] + length[i] - 1;
          for (k = 0;k < length[i];k++) {
            *pisequence++ = *cisequence--;
          }
          break;
        }
        }
      }

      else {
        int_sequence[i][j] = 0;
        real_sequence[i][j] = new double[length[i]];

        prsequence = real_sequence[i][j];

        switch (reverse_flag) {

        case false : {
          crsequence = seq.real_sequence[i][j];
          for (k = 0;k < length[i];k++) {
            *prsequence++ = *crsequence++;
          }
          break;
        }

        case true : {
          crsequence = seq.real_sequence[i][j] + length[i] - 1;
          for (k = 0;k < length[i];k++) {
            *prsequence++ = *crsequence--;
          }
          break;
        }
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
  int *pisequence , *cisequence;
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

  hlength = new Histogram(*(seq.hlength));

  index_parameter_type = seq.index_parameter_type;

  if (seq.hindex_parameter) {
    hindex_parameter = new Histogram(*(seq.hindex_parameter));
  }
  else {
    hindex_parameter = 0;
  }

  if (seq.index_interval) {
    index_interval = new Histogram(*(seq.index_interval));
  }
  else {
    index_interval = 0;
  }

  if (seq.index_parameter) {
    index_parameter = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      index_parameter[i] = new int[length[i]];
      for (j = 0;j < length[i];j++) {
        index_parameter[i][j] = seq.index_parameter[i][j];
      }
    }
  }
  else {
    index_parameter = 0;
  }

  nb_variable = seq.nb_variable + 1;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal = new Histogram*[nb_variable];

  type[0] = STATE;
  min_value[0] = 0.;
  max_value[0] = 0.;
  marginal[0] = 0;

  for (i = 0;i < seq.nb_variable;i++) {
    type[i + 1] = (seq.type[i] == STATE ? INT_VALUE : seq.type[i]);
    min_value[i + 1] = seq.min_value[i];
    max_value[i + 1] = seq.max_value[i];

    if (seq.marginal[i]) {
      marginal[i + 1] = new Histogram(*(seq.marginal[i]));
    }
    else {
      marginal[i + 1] = 0;
    }
  }

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];

    int_sequence[i][0] = new int[length[i]];
    real_sequence[i][0] = 0;

    for (j = 0;j < seq.nb_variable;j++) {
      if (seq.type[j] != REAL_VALUE) {
        int_sequence[i][j + 1] = new int[length[i]];
        real_sequence[i][j + 1] = 0;

        pisequence = int_sequence[i][j + 1];
        cisequence = seq.int_sequence[i][j];
        for (k = 0;k < length[i];k++) {
          *pisequence++ = *cisequence++;
        }
      }

      else {
        int_sequence[i][j + 1] = 0;
        real_sequence[i][j + 1] = new double[length[i]];

        prsequence = real_sequence[i][j + 1];
        crsequence = seq.real_sequence[i][j];
        for (k = 0;k < length[i];k++) {
          *prsequence++ = *crsequence++;
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
  int *pisequence , *cisequence;
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

  hlength = new Histogram(*(seq.hlength));

  index_parameter_type = IMPLICIT_TYPE;
  hindex_parameter = 0;
  index_interval = 0;
  index_parameter = 0;

  nb_variable = seq.nb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = seq.type[i];
    min_value[i] = seq.min_value[i];
    max_value[i] = seq.max_value[i];

    if (seq.marginal[i]) {
      marginal[i] = new Histogram(*(seq.marginal[i]));
    }
    else {
      marginal[i] = 0;
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
        real_sequence[i][j] = 0;

        pisequence = int_sequence[i][j];
        cisequence = seq.int_sequence[i][j];
        for (k = 0;k < length[i];k++) {
          *pisequence++ = *cisequence++;
        }
      }

      else {
        int_sequence[i][j] = 0;
        real_sequence[i][j] = new double[length[i]];

        prsequence = real_sequence[i][j];
        crsequence = seq.real_sequence[i][j];
        for (k = 0;k < length[i];k++) {
          *prsequence++ = *crsequence++;
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe Sequences.
 *
 *  arguments : reference sur un objet Sequences, type de transformation ('c' : copie,
 *              'a' : addition d'une variable, 'r' : suppression du parametre d'index),
 *              inversion.
 *
 *--------------------------------------------------------------*/

Sequences::Sequences(const Sequences &seq , char transform , int param)

{
  switch (transform) {
  case 'c' :
    Sequences::copy(seq , (param == REVERSE ? true : false));
    break;
  case 'a' :
    Sequences::add_state_variable(seq);
    break;
  case 'r' :
    Sequences::remove_index_parameter(seq);
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
  delete hlength;

  delete hindex_parameter;
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

  if (marginal) {
    for (i = 0;i < nb_variable;i++) {
      delete marginal[i];
    }
    delete [] marginal;
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
 *  Construction des sequences reelles pour les variables entieres.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Sequences::build_real_sequence(int variable)

{
  register int i , j , k;
  int *pisequence;
  double *prsequence;


  for (i = 0;i < nb_variable;i++) {
    if (((variable == I_DEFAULT) || (variable == i)) &&
        (type[i] == INT_VALUE) && (!real_sequence[0][i])) {
      for (j = 0;j < nb_sequence;j++) {
        real_sequence[j][i] = new double[length[j]];

        prsequence = real_sequence[j][i];
        pisequence = int_sequence[j][i];
        for (k = 0;k < length[j];k++) {
          *prsequence++ = *pisequence++;
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des sequences reelles pour les variables entieres.
 *
 *--------------------------------------------------------------*/

void Sequences::remove_real_sequence()

{
  register int i , j;

  for (i = 0;i < nb_variable;i++) {
    if (type[i] == INT_VALUE) {
      for (j = 0;j < nb_sequence;j++) {
        delete [] real_sequence[j][i];
        real_sequence[j][i] = 0;
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Extraction de la loi marginale empirique pour une variable entiere.
 *
 *  arguments : reference sur un objet Format_error, variable.
 *
 *--------------------------------------------------------------*/

Distribution_data* Sequences::extract(Format_error &error , int variable) const

{
  bool status = true;
  Distribution_data *histo;


  histo = 0;
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

    else if (!marginal[variable]) {
      status = false;
      error.update(STAT_error[STATR_MARGINAL_HISTOGRAM]);
    }
  }

  if (status) {
    histo = new Distribution_data(*(marginal[variable]));
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
  int offset , *itype;
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
    if (type[i] == STATE) {
      itype[i + offset] = INT_VALUE;
    }
    else if (type[i] == AUXILIARY) {
      itype[i + offset] = REAL_VALUE;
    }
    else {
      itype[i + offset] = type[i];
    }
  }

  vec = new Vectors(cumul_length , 0 , nb_variable + offset , itype);
  delete [] itype;

  i = 0;
  for (j = 0;j < nb_sequence;j++) {
    for (k = 0;k < length[j];k++) {
      if (index_variable) {
        if (index_parameter) {
          vec->int_vector[i][0] = index_parameter[j][k];
        }
        else {
          vec->int_vector[i][0] = k;
        }
      }

      for (m = 0;m < nb_variable;m++) {
        if ((type[m] != REAL_VALUE) && (type[m] != AUXILIARY)) {
          vec->int_vector[i][m + offset] = int_sequence[j][m][k];
        }
        else {
          vec->real_vector[i][m + offset] = real_sequence[j][m][k];
        }
      }

      i++;
    }
  }

  if (index_variable) {
    if (index_parameter) {
      vec->min_value[0] = hindex_parameter->offset;
      vec->max_value[0] = hindex_parameter->nb_value - 1;
    }
    else {
      vec->min_value[0] = 0;
      vec->max_value[0] = max_length - 1;
    }

    vec->build_marginal_histogram(0);
  }

  for (i = 0;i < nb_variable;i++) {
    vec->min_value[i + offset] = min_value[i];
    vec->max_value[i + offset] = max_value[i];

    if (marginal[i]) {
      vec->marginal[i + offset] = new Histogram(*marginal[i]);
      vec->mean[i + offset] = vec->marginal[i + offset]->mean;
      vec->covariance[i + offset][i + offset] = vec->marginal[i + offset]->variance;
    }
    else {
      vec->mean_computation(i + offset);
      vec->variance_computation(i + offset);
    }
  }

  vec->covariance_computation();

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de mesures globales (longueur, temps avant la 1ere occurrence
 *  d'une valeur, nombre de series/d'occurrences d'une valeur) par sequence.
 *
 *  arguments : reference sur un objet Format_error, type, variable, valeur.
 *
 *--------------------------------------------------------------*/

Vectors* Sequences::extract_vectors(Format_error &error , int feature_type ,
                                    int variable , int value) const

{
  bool status = true;
  register int i , j;
  int count , *pisequence , itype[1];
  double *prsequence;
  Vectors *vec;


  vec = 0;
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
            ((marginal[variable]) && (marginal[variable]->frequency[value] == 0))) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VALUE] << " " << value << " "
                        << SEQ_error[SEQR_NOT_PRESENT];
          error.update((error_message.str()).c_str());
        }
      }
    }
  }

  if (status) {
    if (feature_type == SEQUENCE_CUMUL) {
      itype[0] = type[variable];
    }
    else if (feature_type == SEQUENCE_MEAN) {
      itype[0] = REAL_VALUE;
    }
    else {
      itype[0] = INT_VALUE;
    }

    vec = new Vectors(nb_sequence , identifier , 1 , itype);

    switch (feature_type) {

    case LENGTH : {
      for (i = 0;i < nb_sequence;i++) {
        vec->int_vector[i][0] = length[i];
      }
      break;
    }

    case SEQUENCE_CUMUL : {
      if (type[variable] != REAL_VALUE) {
        for (i = 0;i < nb_sequence;i++) {
          vec->int_vector[i][0] = 0;
          pisequence = int_sequence[i][variable];
          for (j = 0;j < length[i];j++) {
            vec->int_vector[i][0] += *pisequence++;
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          vec->real_vector[i][0] = 0.;
          prsequence = real_sequence[i][variable];
          for (j = 0;j < length[i];j++) {
            vec->real_vector[i][0] += *prsequence++;
          }
        }
      }
      break;
    }

    case SEQUENCE_MEAN : {
      if (type[variable] != REAL_VALUE) {
        for (i = 0;i < nb_sequence;i++) {
          vec->real_vector[i][0] = 0.;
          pisequence = int_sequence[i][variable];
          for (j = 0;j < length[i];j++) {
            vec->real_vector[i][0] += *pisequence++;
          }
          vec->real_vector[i][0] /= length[i];
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          vec->real_vector[i][0] = 0.;
          prsequence = real_sequence[i][variable];
          for (j = 0;j < length[i];j++) {
            vec->real_vector[i][0] += *prsequence++;
          }
          vec->real_vector[i][0] /= length[i];
        }
      }
      break;
    }

    case FIRST_OCCURRENCE : {
      for (i = 0;i < nb_sequence;i++) {
        vec->int_vector[i][0] = -1;
        pisequence = int_sequence[i][variable];
        for (j = 0;j < length[i];j++) {
          if (*pisequence++ == value) {
            vec->int_vector[i][0] = j;
            break;
          }
        }
      }
      break;
    }

    case NB_RUN : {
      for (i = 0;i < nb_sequence;i++) {
        pisequence = int_sequence[i][variable];
        count = 0;
        if (*pisequence++ == value) {
          count++;
        }
        for (j = 1;j < length[i];j++) {
          if ((*pisequence != *(pisequence - 1)) && (*pisequence == value)) {
            count++;
          }
          pisequence++;
        }

        vec->int_vector[i][0] = count;
      }
      break;
    }

    case NB_OCCURRENCE : {
      for (i = 0;i < nb_sequence;i++) {
        pisequence = int_sequence[i][variable];
        count = 0;
        for (j = 0;j < length[i];j++) {
          if (*pisequence++ == value) {
            count++;
          }
        }

        vec->int_vector[i][0] = count;
      }
      break;
    }
    }

    vec->min_value_computation(0);
    vec->max_value_computation(0);
    vec->build_marginal_histogram(0);
  }

  return vec;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Markovian_sequences a partir d'un objet Sequences.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Markovian_sequences* Sequences::markovian_sequences(Format_error &error) const

{
  bool status = true;
  register int i;
  Markovian_sequences *seq;


  seq = 0;
  error.init();

  if (((index_parameter_type == TIME) && (index_interval->variance > 0.)) ||
      (index_parameter_type == POSITION)) {
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
      if (min_value[i] < 0.) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_POSITIVE_MIN_VALUE];
        error.update((error_message.str()).c_str());
      }

      if (!marginal[i]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                      << STAT_error[STATR_MARGINAL_HISTOGRAM];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {
    seq = new Markovian_sequences(*this);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Tops a partir d'un objet Sequences.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Tops* Sequences::tops(Format_error &error) const

{
  bool status = true;
  register int i;
  Tops *tops;


  tops = 0;
  error.init();

  if (index_parameter_type != POSITION) {
    status = false;
    error.correction_update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE] , SEQ_index_parameter_word[POSITION]);
  }

  if (nb_variable != 1) {
    status = false;
    error.correction_update(STAT_error[STATR_NB_VARIABLE] , 1);
  }

  else {
    if (type[0] != NB_INTERNODE) {
      status = false;
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , STAT_variable_word[NB_INTERNODE]);
    }

    else if (!marginal[0]) {
      status = false;
      error.update(STAT_error[STATR_MARGINAL_HISTOGRAM]);
    }
  }

  if (status) {
    for (i = 0;i < nb_sequence;i++) {
      if (index_parameter[i][0] == 0) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEQUENCE] << " " << i + 1 << ": "
                      << SEQ_parsing[SEQP_POSITION];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {
    tops = new Tops(*this);
  }

  return tops;
}


/*--------------------------------------------------------------*
 *
 *  Verification du caractere (strictement) croissant des parametres d'index.
 *
 *  arguments : reference sur un objet Format_error, flag croissance strict ou non,
 *              label de l'objet.
 *
 *--------------------------------------------------------------*/

bool Sequences::increasing_index_parameter_checking(Format_error &error , bool strict ,
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
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              flag croissance strict ou non, labels de l'objet et de la variable.
 *
 *--------------------------------------------------------------*/

bool Sequences::increasing_sequence_checking(Format_error &error , int variable , bool strict ,
                                             const char *pattern_label , const char *variable_label) const

{
  bool status = true;
  register int i , j;


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

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Verification d'un objet Sequences.
 *
 *  arguments : reference sur un objet Format_error, label.
 *
 *--------------------------------------------------------------*/

bool Sequences::check(Format_error &error , const char *pattern_label)

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
      build_index_parameter_histogram();
    }
    if ((index_parameter_type == TIME) || ((index_parameter_type == POSITION) &&
         (type[0] != NB_INTERNODE))) {
      index_interval_computation();
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'un objet Time_events a partir d'un objet Sequences.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              nombre d'evenements, dates de debut et de fin,
 *              dates precedente et suivante.
 *
 *--------------------------------------------------------------*/

Time_events* Sequences::extract_time_events(Format_error &error , int variable ,
                                            int begin_date , int end_date ,
                                            int previous_date , int next_date) const

{
  bool status = true , lstatus;
  register int i , j;
  int nb_element , previous , begin , end , next , *time , *nb_event , *pdate;
  Time_events *timev;


  timev = 0;
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
      timev = new Time_events(nb_element , time , nb_event);
    }

    delete [] time;
    delete [] nb_event;
  }

  return timev;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'un objet Renewal_data a partir d'un objet Sequences.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              nombre d'evenements, index de debut et de fin.
 *
 *--------------------------------------------------------------*/

Renewal_data* Sequences::extract_renewal_data(Format_error &error , int variable ,
                                              int begin_index_parameter , int end_index_parameter) const

{
  bool status = true , lstatus;
  register int i , j;
  int nb_element , index , *ptime , *pnb_event , *pisequence , *cisequence;
  Renewal_data *timev;


  timev = 0;
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
    timev = new Renewal_data(nb_sequence , end_index_parameter + 1 - begin_index_parameter);

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
    // des histogrammes du temps d'observation et du nombre d'evenements

    ptime -= nb_element;
    pnb_event -= nb_element;

    timev->build(nb_element , ptime , pnb_event);
    delete [] ptime;
    delete [] pnb_event;

    // extraction des caracteristiques des histogrammes des intervalles de temps entre 2 evenements,
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
      timev = 0;
      error.update(SEQ_error[SEQR_BOTH_END_CENSORED_INTERVAL]);
    }
  }

  return timev;
}


/*--------------------------------------------------------------*
 *
 *  Fusion d'objets Sequences.
 *
 *  arguments : reference sur un objet Format_error, nombre d'objets Sequences,
 *              pointeurs sur les objets Sequences.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::merge(Format_error &error , int nb_sample , const Sequences **iseq) const

{
  bool status = true;
  register int i , j , k , m , n;
  int inb_sequence , *ilength , *pindex_param , *cindex_param , *pisequence , *cisequence;
  double *prsequence , *crsequence;
  const Histogram **phisto;
  Sequences *seq;
  const Sequences **pseq;


  seq = 0;
  error.init();

  for (i = 0;i < nb_sample;i++) {
    if (iseq[i]->index_parameter_type != index_parameter_type) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_INDEX_PARAMETER_TYPE];

      if (index_parameter_type == IMPLICIT_TYPE) {
        error.update((error_message.str()).c_str());
      }
      else {
        error.correction_update((error_message.str()).c_str() , SEQ_index_parameter_word[index_parameter_type]);
      }
    }
  }

  for (i = 0;i < nb_sample;i++) {
    if (iseq[i]->nb_variable != nb_variable) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << STAT_error[STATR_NB_VARIABLE];
      error.correction_update((error_message.str()).c_str() , nb_variable);
    }

    else {
      for (j = 0;j < nb_variable;j++) {
        if (iseq[i]->type[j] != type[j]) {
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
    pseq = new const Sequences*[nb_sample];

    pseq[0] = this;
    for (i = 1;i < nb_sample;i++) {
      pseq[i] = iseq[i - 1];
    }

    // calcul du nombre et de la longueur des sequences

    inb_sequence = 0;
    for (i = 0;i < nb_sample;i++) {
      inb_sequence += pseq[i]->nb_sequence;
    }

    ilength = new int[inb_sequence];

    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < pseq[j]->nb_sequence;k++) {
        ilength[i++] = pseq[j]->length[k];
      }
    }

    seq = new Sequences(inb_sequence , 0 , ilength , index_parameter_type ,
                        nb_variable , type);
    delete [] ilength;

    phisto = new const Histogram*[nb_sample];

    // copie des sequences

    if (seq->index_parameter) {
      i = 0;
      for (j = 0;j < nb_sample;j++) {
        for (k = 0;k < pseq[j]->nb_sequence;k++) {
          pindex_param = seq->index_parameter[i];
          cindex_param = pseq[j]->index_parameter[k];
          for (m = 0;m < (pseq[j]->index_parameter_type == POSITION ? pseq[j]->length[k] + 1 : pseq[j]->length[k]);m++) {
            *pindex_param++ = *cindex_param++;
          }
          i++;
        }
      }

      for (i = 0;i < nb_sample;i++) {
        phisto[i] = pseq[i]->hindex_parameter;
      }
      seq->hindex_parameter = new Histogram(nb_sample , phisto);
    }

    if ((seq->index_parameter_type == TIME) || ((seq->index_parameter_type == POSITION) &&
         (seq->type[0] != NB_INTERNODE))) {
      for (i = 0;i < nb_sample;i++) {
        phisto[i] = pseq[i]->index_interval;
      }
      seq->index_interval = new Histogram(nb_sample , phisto);
    }

    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < pseq[j]->nb_sequence;k++) {
        for (m = 0;m < pseq[j]->nb_variable;m++) {
          if ((pseq[j]->type[m] != REAL_VALUE) && (pseq[j]->type[m] != AUXILIARY)) {
            pisequence = seq->int_sequence[i][m];
            cisequence = pseq[j]->int_sequence[k][m];
            for (n = 0;n < pseq[j]->length[k];n++) {
              *pisequence++ = *cisequence++;
            }
          }

          else {
            prsequence = seq->real_sequence[i][m];
            crsequence = pseq[j]->real_sequence[k][m];
            for (n = 0;n < pseq[j]->length[k];n++) {
              *prsequence++ = *crsequence++;
            }
          }
        }
        i++;
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value[i] = pseq[0]->min_value[i];
      seq->max_value[i] = pseq[0]->max_value[i];
      for (j = 1;j < nb_sample;j++) {
        if (pseq[j]->min_value[i] < seq->min_value[i]) {
          seq->min_value[i] = pseq[j]->min_value[i];
        }
        if (pseq[j]->max_value[i] > seq->max_value[i]) {
          seq->max_value[i] = pseq[j]->max_value[i];
        }
      }

      for (j = 0;j < nb_sample;j++) {
        if (pseq[j]->marginal[i]) {
          phisto[j] = pseq[j]->marginal[i];
        }
        else {
          break;
        }
      }

      if (j == nb_sample) {
        seq->marginal[i] = new Histogram(nb_sample , phisto);
      }
    }

    delete [] pseq;
    delete [] phisto;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Translation des valeurs d'une variable.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              parametre de translation.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::shift(Format_error &error , int variable , int shift_param) const

{
  bool status = true;
  register int i , j , k;
  int *pisequence , *cisequence;
  double *prsequence , *crsequence;
  Sequences *seq;


  seq = 0;
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
    seq = new Sequences(nb_sequence , identifier , length , index_parameter_type ,
                        nb_variable , type);

    // copie des parametres d'index

    if (hindex_parameter) {
      seq->hindex_parameter = new Histogram(*hindex_parameter);
    }
    if (index_interval) {
      seq->index_interval = new Histogram(*index_interval);
    }

    if (index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < (seq->index_parameter_type == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }
      }
    }

    for (i = 0;i < seq->nb_sequence;i++) {
      for (j = 0;j < seq->nb_variable;j++) {
        if ((seq->type[j] != REAL_VALUE) && (seq->type[j] != AUXILIARY)) {
          pisequence = seq->int_sequence[i][j];
          cisequence = int_sequence[i][j];

          // translation des valeurs entieres

          if (j == variable) {
            for (k = 0;k < seq->length[i];k++) {
              *pisequence++ = *cisequence++ + shift_param;
            }
          }

          // copie des valeurs entieres

          else {
            for (k = 0;k < seq->length[i];k++) {
              *pisequence++ = *cisequence++;
            }
          }
        }

        else {
          prsequence = seq->real_sequence[i][j];
          crsequence = real_sequence[i][j];

          // translation des valeurs reelles

          if ((j == variable) || ((j == variable + 1) && (seq->type[j] == AUXILIARY))) {
            for (k = 0;k < seq->length[i];k++) {
              *prsequence++ = *crsequence++ + shift_param;
            }
          }

          // copie des valeurs reelles

          else {
            for (k = 0;k < seq->length[i];k++) {
              *prsequence++ = *crsequence++;
            }
          }
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      if (i == variable) {
        seq->min_value[i] = min_value[i] + shift_param;
        seq->max_value[i] = max_value[i] + shift_param;

        if ((seq->type[i] == INT_VALUE) && (seq->min_value[i] >= 0) &&
            (seq->max_value[i] <= MARGINAL_MAX_VALUE)) {
          if (marginal[i]) {
            seq->marginal[i] = new Histogram(*marginal[i] , 's' , shift_param);
          }
          else {
            seq->build_marginal_histogram(i);
          }
        }
      }

      else {
        seq->min_value[i] = min_value[i];
        seq->max_value[i] = max_value[i];
        if (marginal[i]) {
          seq->marginal[i] = new Histogram(*marginal[i]);
        }
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Translation des valeurs d'une variable reelle.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              parametre de translation.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::shift(Format_error &error , int variable , double shift_param) const

{
  bool status = true;
  register int i , j , k;
  int *pisequence , *cisequence;
  double *prsequence , *crsequence;
  Sequences *seq;


  seq = 0;
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
    seq = new Sequences(nb_sequence , identifier , length , index_parameter_type ,
                        nb_variable , type);

    // copie des parametres d'index

    if (hindex_parameter) {
      seq->hindex_parameter = new Histogram(*hindex_parameter);
    }
    if (index_interval) {
      seq->index_interval = new Histogram(*index_interval);
    }

    if (index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < (seq->index_parameter_type == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }
      }
    }

    for (i = 0;i < seq->nb_sequence;i++) {
      for (j = 0;j < seq->nb_variable;j++) {

        // copie des valeurs entieres

        if ((seq->type[j] != REAL_VALUE) && (seq->type[j] != AUXILIARY)) {
          pisequence = seq->int_sequence[i][j];
          cisequence = int_sequence[i][j];
          for (k = 0;k < seq->length[i];k++) {
            *pisequence++ = *cisequence++;
          }
        }

        else {
          prsequence = seq->real_sequence[i][j];
          crsequence = real_sequence[i][j];

          // translation des valeurs reelles

          if ((j == variable) || ((j == variable + 1) && (seq->type[j] == AUXILIARY))) {
            for (k = 0;k < seq->length[i];k++) {
              *prsequence++ = *crsequence++ + shift_param;
            }
          }

          // copie des valeurs reelles

          else {
            for (k = 0;k < seq->length[i];k++) {
              *prsequence++ = *crsequence++;
            }
          }
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      if (i == variable) {
        seq->min_value[i] = min_value[i] + shift_param;
        seq->max_value[i] = max_value[i] + shift_param;
      }

      else {
        seq->min_value[i] = min_value[i];
        seq->max_value[i] = max_value[i];
        if (marginal[i]) {
          seq->marginal[i] = new Histogram(*marginal[i]);
        }
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une variable.
 *
 *  arguments : reference sur un objet Sequences, indice de la variable,
 *              pas de regroupement.
 *
 *--------------------------------------------------------------*/

void Sequences::cluster(const Sequences &seq , int variable , int step)

{
  register int i , j , k;
  int *pisequence , *cisequence;
  double *prsequence , *crsequence;


  // copie des parametres d'index

  if (seq.hindex_parameter) {
    hindex_parameter = new Histogram(*(seq.hindex_parameter));
  }
  if (seq.index_interval) {
    index_interval = new Histogram(*(seq.index_interval));
  }

  if (seq.index_parameter) {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < (index_parameter_type == POSITION ? length[i] + 1 : length[i]);j++) {
        index_parameter[i][j] = seq.index_parameter[i][j];
      }
    }
  }

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < nb_variable;j++) {
      if ((type[j] != REAL_VALUE) && (type[j] != AUXILIARY)) {
        pisequence = int_sequence[i][j];
        cisequence = seq.int_sequence[i][j];

        // regroupement des valeurs entieres

        if (j == variable) {
          for (k = 0;k < length[i];k++) {
            *pisequence++ = *cisequence++ / step;
          }
        }

        // copie des valeurs entieres

        else {
          for (k = 0;k < length[i];k++) {
            *pisequence++ = *cisequence++;
          }
        }
      }

      else {
        prsequence = real_sequence[i][j];
        crsequence = seq.real_sequence[i][j];

        // regroupement des valeurs reelles

        if ((j == variable) || ((j == variable + 1) && (type[j] == AUXILIARY))) {
          for (k = 0;k < length[i];k++) {
            *prsequence++ = *crsequence++ / step;
          }
        }

        // copie des valeurs reelles

        else {
          for (k = 0;k < length[i];k++) {
            *prsequence++ = *crsequence++;
          }
        }
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    if (i == variable) {
      if (type[i] != REAL_VALUE) {
        min_value[i] = (int)seq.min_value[i] / step;
        max_value[i] = (int)seq.max_value[i] / step;
      }
      else {
        min_value[i] = seq.min_value[i] / step;
        max_value[i] = seq.max_value[i] / step;
      }

      if (seq.marginal[i]) {
        marginal[i] = new Histogram(*(seq.marginal[i]) , 'c' , step);
      }
      else {
        build_marginal_histogram(i);
      }
    }

    else {
      min_value[i] = seq.min_value[i];
      max_value[i] = seq.max_value[i];
      if (seq.marginal[i]) {
        marginal[i] = new Histogram(*(seq.marginal[i]));
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une variable.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              pas de regroupement.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::cluster(Format_error &error , int variable , int step) const

{
  bool status = true;
  Sequences *seq;


  seq = 0;
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

    if ((type[variable] == INT_VALUE) && (variable + 1 < nb_variable) && ((type[variable + 1] != INT_VALUE) &&
         (type[variable + 1] != STATE) && (type[variable + 1] != REAL_VALUE))) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (step < 1) {
    status = false;
    error.update(STAT_error[STATR_CLUSTERING_STEP]);
  }

  if (status) {
    seq = new Sequences(nb_sequence , identifier , length , index_parameter_type ,
                        nb_variable , type);
    seq->cluster(*this , variable , step);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Transcodage des symboles d'une variable entiere.
 *
 *  arguments : reference sur un objet Sequences, indice de la variable,
 *              plus petit et plus grand symboles, table de transcodage des symboles,
 *              flag pour ajouter une variable.
 *
 *--------------------------------------------------------------*/

void Sequences::transcode(const Sequences &seq , int ivariable , int min_symbol ,
                          int max_symbol , int *symbol , bool add_flag)

{
  register int i , j , k;
  int variable , offset , *pisequence , *cisequence;
  double *prsequence , *crsequence;


  // copie des parametres d'index

  if (seq.hindex_parameter) {
    hindex_parameter = new Histogram(*(seq.hindex_parameter));
  }
  if (seq.index_interval) {
    index_interval = new Histogram(*(seq.index_interval));
  }

  if (seq.index_parameter) {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < (index_parameter_type == POSITION ? length[i] + 1 : length[i]);j++) {
        index_parameter[i][j] = seq.index_parameter[i][j];
      }
    }
  }

  switch (add_flag) {
  case false :
    variable = ivariable;
    offset = 0;
    break;
  case true :
    variable = 0;
    offset = 1;
    break;
  }

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < nb_variable;j++) {
      if ((type[j] != REAL_VALUE) && (type[j] != AUXILIARY)) {
        pisequence = int_sequence[i][j];

        // transcodage des symboles

        if (j == variable) {
          cisequence = seq.int_sequence[i][ivariable];
          for (k = 0;k < length[i];k++) {
            *pisequence++ = symbol[*cisequence++ - (int)seq.min_value[variable]] + min_symbol;
          }
        }

        // copie des valeurs entieres

        else {
          cisequence = seq.int_sequence[i][j - offset];
          for (k = 0;k < length[i];k++) {
            *pisequence++ = *cisequence++;
          }
        }
      }

      // copie des valeurs reelles

      else {
        prsequence = real_sequence[i][j];
        crsequence = seq.real_sequence[i][j - offset];
        for (k = 0;k < length[i];k++) {
          *prsequence++ = *crsequence++;
        }
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    if (i == variable) {
      min_value[i] = min_symbol;
      max_value[i] = max_symbol;
      build_marginal_histogram(i);
    }

    else {
      min_value[i] = seq.min_value[i - offset];
      max_value[i] = seq.max_value[i - offset];
      if (seq.marginal[i - offset]) {
        marginal[i] = new Histogram(*(seq.marginal[i - offset]));
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Transcodage des symboles d'une variable entiere.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              table de transcodage des symboles.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::transcode(Format_error &error , int variable , int *symbol) const

{
  bool status = true , *presence;
  register int i;
  int min_symbol , max_symbol;
  Sequences *seq;


  seq = 0;
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

    if ((variable + 1 < nb_variable) && ((type[variable + 1] != INT_VALUE) &&
         (type[variable + 1] != STATE) && (type[variable + 1] != REAL_VALUE))) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }

    if (status) {
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

      seq = new Sequences(nb_sequence , identifier , length , index_parameter_type ,
                          nb_variable , type);
      seq->transcode(*this , variable , min_symbol , max_symbol , symbol);
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une variable entiere.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              nombre de classes, bornes pour regrouper les valeurs.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::cluster(Format_error &error , int variable ,
                              int nb_class , int *ilimit) const

{
  bool status = true;
  register int i , j , k;
  int *int_limit , *symbol , *itype;
  double *real_limit;
  Sequences *seq;


  seq = 0;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != STATE) && (type[variable] != REAL_VALUE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    else if ((nb_class < 2) || (nb_class >= (int)(max_value[variable] - min_value[variable]) + 1)) {
      status = false;
      error.update(STAT_error[STATR_NB_CLASS]);
    }

    if ((variable + 1 < nb_variable) && ((type[variable + 1] != INT_VALUE) &&
         (type[variable + 1] != STATE) && (type[variable + 1] != REAL_VALUE))) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (status) {
    if ((type[variable] == INT_VALUE) || (type[variable] == STATE)) {
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

        seq = new Sequences(nb_sequence , identifier , length , index_parameter_type ,
                            nb_variable , type);
        seq->transcode(*this , variable , 0 , nb_class - 1 , symbol);

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
        itype = new int[nb_variable];

        for (i = 0;i < nb_variable;i++) {
          itype[i] = type[i];
        }
        itype[variable] = INT_VALUE;

        seq = new Sequences(nb_sequence , identifier , length , index_parameter_type ,
                            nb_variable , itype);
        delete [] itype;

        seq->cluster(*this , variable , nb_class , real_limit);
      }

      delete [] real_limit;
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une variable reelle.
 *
 *  arguments : reference sur un objet Sequences, indice de la variable,
 *              nombre de classes, bornes pour regrouper les valeurs.
 *
 *--------------------------------------------------------------*/

void Sequences::cluster(const Sequences &seq , int variable , int nb_class , double *limit)

{
  register int i , j , k , m;
  int *pisequence , *cisequence;
  double *prsequence , *crsequence;


  // copie des parametres d'index

  if (seq.hindex_parameter) {
    hindex_parameter = new Histogram(*(seq.hindex_parameter));
  }
  if (seq.index_interval) {
    index_interval = new Histogram(*(seq.index_interval));
  }

  if (seq.index_parameter) {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < (index_parameter_type == POSITION ? length[i] + 1 : length[i]);j++) {
        index_parameter[i][j] = seq.index_parameter[i][j];
      }
    }
  }

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < nb_variable;j++) {

      // copie des valeurs entieres

      if ((seq.type[j] != REAL_VALUE) && (seq.type[j] != AUXILIARY)) {
        pisequence = int_sequence[i][j];
        cisequence = seq.int_sequence[i][j];
        for (k = 0;k < length[i];k++) {
          *pisequence++ = *cisequence++;
        }
      }

      else {
        crsequence = seq.real_sequence[i][j];

        // regroupement des valeurs reelles

        if (j == variable) {
          pisequence = int_sequence[i][j];
          for (k = 0;k < length[i];k++) {
            for (m = 0;m < nb_class;m++) {
              if (*crsequence < limit[m + 1]) {
                *pisequence++ = m;
                break;
              }
            }
            crsequence++;
          }
        }

        // copie des valeurs reelles

        else {
          prsequence = real_sequence[i][j];
          for (k = 0;k < length[i];k++) {
            *prsequence++ = *crsequence++;
          }
        }
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    if (i == variable) {
      min_value_computation(i);
      max_value_computation(i);
      build_marginal_histogram(i);
    }

    else {
      min_value[i] = seq.min_value[i];
      max_value[i] = seq.max_value[i];
      if (seq.marginal[i]) {
        marginal[i] = new Histogram(*(seq.marginal[i]));
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une variable reelle.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              nombre de classes, bornes pour regrouper les valeurs.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::cluster(Format_error &error , int variable ,
                              int nb_class , double *ilimit) const

{
  bool status = true;
  register int i;
  int *itype;
  double *limit;
  Sequences *seq;


  seq = 0;
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

    if ((variable + 1 < nb_variable) && ((type[variable + 1] != INT_VALUE) &&
         (type[variable + 1] != STATE) && (type[variable + 1] != REAL_VALUE))) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
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
      itype = new int[nb_variable];

      for (i = 0;i < nb_variable;i++) {
        itype[i] = type[i];
      }
      itype[variable] = INT_VALUE;

      seq = new Sequences(nb_sequence , identifier , length , index_parameter_type ,
                          nb_variable , itype);
      delete [] itype;

      seq->cluster(*this , variable , nb_class , limit);
    }

    delete [] limit;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Changement d'unite d'une variable.
 *
 *  arguments : reference sur un objet Format_error, variable, facteur d'echelle.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::scaling(Format_error &error , int variable , int scaling_coeff) const

{
  bool status = true;
  register int i , j , k;
  int *pisequence , *cisequence;
  double *prsequence , *crsequence;
  Sequences *seq;


  seq = 0;
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
    seq = new Sequences(nb_sequence , identifier , length , index_parameter_type ,
                        nb_variable , type);

    if (hindex_parameter) {
      seq->hindex_parameter = new Histogram(*hindex_parameter);
    }
    if (index_interval) {
      seq->index_interval = new Histogram(*index_interval);
    }

    if (index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < (seq->index_parameter_type == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }
      }
    }

    for (i = 0;i < seq->nb_sequence;i++) {
      for (j = 0;j < seq->nb_variable;j++) {
        if ((seq->type[j] != REAL_VALUE) && (seq->type[j] != AUXILIARY)) {
          pisequence = seq->int_sequence[i][j];
          cisequence = int_sequence[i][j];

          // mise a l'echelle des valeurs entieres

          if (j == variable) {
            for (k = 0;k < seq->length[i];k++) {
              *pisequence++ = *cisequence++ * scaling_coeff;
            }
          }

          // copie des valeurs entieres

          else {
            for (k = 0;k < seq->length[i];k++) {
              *pisequence++ = *cisequence++;
            }
          }
        }

        else {
          prsequence = seq->real_sequence[i][j];
          crsequence = real_sequence[i][j];

          // mise a l'echelle des valeurs reelles

          if ((j == variable) || ((j == variable + 1) && (seq->type[j] == AUXILIARY))) {
            for (k = 0;k < seq->length[i];k++) {
              *prsequence++ = *crsequence++ * scaling_coeff;
            }
          }

          // copie des valeurs reelles

          else {
            for (k = 0;k < seq->length[i];k++) {
              *prsequence++ = *crsequence++;
            }
          }
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      if (i == variable) {
        seq->min_value[i] = min_value[i] * scaling_coeff;
        seq->max_value[i] = max_value[i] * scaling_coeff;
        seq->build_marginal_histogram(i);
      }

      else {
        seq->min_value[i] = min_value[i];
        seq->max_value[i] = max_value[i];
        if (marginal[i]) {
          seq->marginal[i] = new Histogram(*marginal[i]);
        }
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Arrondi des valeurs d'une variable reelle.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::round(Format_error &error , int variable) const

{
  bool status = true;
  register int i , j , k;
  int *itype , *pisequence , *cisequence;
  double *prsequence , *crsequence;
  Sequences *seq;


  seq = 0;
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
        if (type[i] == REAL_VALUE) {
          itype[i] = INT_VALUE;
        }
        else {
          itype[i] = type[i];
        }
      }
    }

    else {
      for (i = 0;i < nb_variable;i++) {
        itype[i] = type[i];
      }
      itype[variable] = INT_VALUE;
    }

    seq = new Sequences(nb_sequence , identifier , length , index_parameter_type ,
                        nb_variable , itype);
    delete [] itype;

    // copie des parametres d'index

    if (hindex_parameter) {
      seq->hindex_parameter = new Histogram(*hindex_parameter);
    }
    if (index_interval) {
      seq->index_interval = new Histogram(*index_interval);
    }

    if (index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < (seq->index_parameter_type == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }
      }
    }

    for (i = 0;i < seq->nb_sequence;i++) {
      for (j = 0;j < seq->nb_variable;j++) {

        // copie des valeurs entieres

        if ((type[j] != REAL_VALUE) && (type[j] != AUXILIARY)) {
          pisequence = seq->int_sequence[i][j];
          cisequence = int_sequence[i][j];
          for (k = 0;k < seq->length[i];k++) {
            *pisequence++ = *cisequence++;
          }
        }

        else {
          crsequence = real_sequence[i][j];

          // arrondi des valeurs reelles

          if (((variable == I_DEFAULT) && (type[j] == REAL_VALUE)) || (variable == j)) {
            pisequence = seq->int_sequence[i][j];
            for (k = 0;k < seq->length[i];k++) {
              *pisequence++ = (int)::round(*crsequence++);
            }
          }

          // copie des valeurs reelles

          else {
            prsequence = seq->real_sequence[i][j];
            for (k = 0;k < seq->length[i];k++) {
              *prsequence++ = *crsequence++;
            }
          }
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      if (((variable == I_DEFAULT) && (type[i] == REAL_VALUE)) || (variable == i)) {
        seq->min_value[i] = ::round(min_value[i]);
        seq->max_value[i] = ::round(max_value[i]);
        seq->build_marginal_histogram(i);
      }

      else {
        seq->min_value[i] = min_value[i];
        seq->max_value[i] = max_value[i];
        if (marginal[i]) {
          seq->marginal[i] = new Histogram(*marginal[i]);
        }
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Selection de sequences sur les valeurs prises par le parametre d'index.
 *
 *  arguments : reference sur un objet Format_error, bornes sur les parametres d'index,
 *              flag pour conserver ou rejeter les sequences selectionnees.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::index_parameter_select(Format_error &error , int min_index_parameter ,
                                             int max_index_parameter , bool keep) const

{
  bool status = true;
  register int i , j;
  int inb_sequence , *index;
  Sequences *seq;


  seq = 0;
  error.init();

  if ((index_parameter_type != TIME) && (index_parameter_type != POSITION)) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }

  else {
    if ((min_index_parameter < 0) || (min_index_parameter >= hindex_parameter->nb_value) ||
        (min_index_parameter > max_index_parameter)) {
      status = false;
      error.update(SEQ_error[SEQR_MIN_INDEX_PARAMETER]);
    }
    if ((max_index_parameter < hindex_parameter->offset) ||
        (max_index_parameter < min_index_parameter)) {
      status = false;
      error.update(SEQ_error[SEQR_MAX_INDEX_PARAMETER]);
    }
  }

  if (status) {

    // selection des sequences

    index = new int[nb_sequence];
    inb_sequence = 0;

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < (index_parameter_type == POSITION ? length[i] + 1 : length[i]);j++) {
        if ((index_parameter[i][j] >= min_index_parameter) &&
            (index_parameter[i][j] <= max_index_parameter)) {
          if (keep) {
            index[inb_sequence++] = i;
          }
          break;
        }
      }

      if ((!keep) && (j == (index_parameter_type == POSITION ? length[i] + 1 : length[i]))) {
        index[inb_sequence++] = i;
      }
    }

    if (inb_sequence == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    // copie des sequences

    if (status) {
      seq = new Sequences(*this , inb_sequence , index);
    }

    delete [] index;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Selection de sequences sur les valeurs prises par une variable.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              bornes sur les valeurs, flag pour conserver ou rejeter
 *              les sequences selectionnees.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::value_select(Format_error &error , int variable , int imin_value ,
                                   int imax_value , bool keep) const

{
  bool status = true;
  register int i , j;
  int inb_sequence , *index , *pisequence;
  double *prsequence;
  Sequences *seq;


  seq = 0;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != STATE) && (type[variable] != REAL_VALUE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    else {
      if ((imin_value > max_value[variable]) || (imin_value > imax_value)) {
        status = false;
        error.update(STAT_error[STATR_MIN_VALUE]);
      }
      if ((imax_value < min_value[variable]) || (imax_value < imin_value)) {
        status = false;
        error.update(STAT_error[STATR_MAX_VALUE]);
      }
    }
  }

  if (status) {

    // selection des sequences

    index = new int[nb_sequence];
    inb_sequence = 0;

    if (type[variable] != REAL_VALUE) {
      for (i = 0;i < nb_sequence;i++) {
        pisequence = int_sequence[i][variable];
        for (j = 0;j < length[i];j++) {
          if ((*pisequence >= imin_value) && (*pisequence <= imax_value)) {
            if (keep) {
              index[inb_sequence++] = i;
            }
            break;
          }

          pisequence++;
        }

        if ((!keep) && (j == length[i])) {
          index[inb_sequence++] = i;
        }
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        prsequence = real_sequence[i][variable];
        for (j = 0;j < length[i];j++) {
          if ((*prsequence >= imin_value) && (*prsequence <= imax_value)) {
            if (keep) {
              index[inb_sequence++] = i;
            }
            break;
          }

          prsequence++;
        }

        if ((!keep) && (j == length[i])) {
          index[inb_sequence++] = i;
        }
      }
    }

    if (inb_sequence == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    // copie des sequences

    if (status) {
      seq = new Sequences(*this , inb_sequence , index);
    }

    delete [] index;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Selection de sequences sur les valeurs prises par une variable reelle.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              bornes sur les valeurs, flag pour conserver ou rejeter
 *              les sequences selectionnees.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::value_select(Format_error &error , int variable , double imin_value ,
                                   double imax_value , bool keep) const

{
  bool status = true;
  register int i , j;
  int inb_sequence , *index;
  double *prsequence;
  Sequences *seq;


  seq = 0;
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

    else {
      if ((imin_value > max_value[variable]) || (imin_value > imax_value)) {
        status = false;
        error.update(STAT_error[STATR_MIN_VALUE]);
      }
      if ((imax_value < min_value[variable]) || (imax_value < imin_value)) {
        status = false;
        error.update(STAT_error[STATR_MAX_VALUE]);
      }
    }
  }

  if (status) {

    // selection des sequences

    index = new int[nb_sequence];
    inb_sequence = 0;

    for (i = 0;i < nb_sequence;i++) {
      prsequence = real_sequence[i][variable];
      for (j = 0;j < length[i];j++) {
        if ((*prsequence >= imin_value) && (*prsequence <= imax_value)) {
          if (keep) {
            index[inb_sequence++] = i;
          }
          break;
        }

        prsequence++;
      }

      if ((!keep) && (j == length[i])) {
        index[inb_sequence++] = i;
      }
    }

    if (inb_sequence == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    // copie des sequences

    if (status) {
      seq = new Sequences(*this , inb_sequence , index);
    }

    delete [] index;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Selection de sequences par l'identificateur.
 *
 *  arguments : reference sur un objet Format_error, nombre de sequences,
 *              identificateur des sequences, flag pour conserver ou rejeter
 *              les sequences selectionnees.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::select_individual(Format_error &error , int inb_sequence ,
                                        int *iidentifier , bool keep) const

{
  bool status = true;
  int *index;
  Sequences *seq;


  seq = 0;
  error.init();

  if ((inb_sequence < 1) || (inb_sequence > (keep ? nb_sequence : nb_sequence - 1))) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }

  else {
    status = selected_identifier_checking(error , nb_sequence , identifier , inb_sequence ,
                                          iidentifier , SEQ_label[SEQL_SEQUENCE]);
  }

  if (status) {
    index = identifier_select(nb_sequence , identifier , inb_sequence , iidentifier , keep);

    seq = new Sequences(*this , (keep ? inb_sequence : nb_sequence - inb_sequence) , index);

    delete [] index;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Suppression du parametre d'index.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::remove_index_parameter(Format_error &error) const

{
  Sequences *seq;


  error.init();

  if (!index_parameter) {
    seq = 0;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  else {
    seq = new Sequences(*this , 'r');
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Selection de variables.
 *
 *  arguments : reference sur un objet Sequences, indices des variables.
 *
 *--------------------------------------------------------------*/

void Sequences::select_variable(const Sequences &seq , int *variable)

{
  register int i , j , k;
  int *pisequence , *cisequence;
  double *prsequence , *crsequence;


  // copie des parametres d'index

  if (seq.hindex_parameter) {
    hindex_parameter = new Histogram(*(seq.hindex_parameter));
  }
  if (seq.index_interval) {
    index_interval = new Histogram(*(seq.index_interval));
  }

  if (seq.index_parameter) {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < (index_parameter_type == POSITION ? length[i] + 1 : length[i]);j++) {
        index_parameter[i][j] = seq.index_parameter[i][j];
      }
    }
  }

  // copie des valeurs

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < nb_variable;j++) {
      if ((type[j] != REAL_VALUE) && (type[j] != AUXILIARY)) {
        pisequence = int_sequence[i][j];
        cisequence = seq.int_sequence[i][variable[j]];
        for (k = 0;k < length[i];k++) {
          *pisequence++ = *cisequence++;
        }
      }

      else {
        prsequence = real_sequence[i][j];
        crsequence = seq.real_sequence[i][variable[j]];
        for (k = 0;k < length[i];k++) {
          *prsequence++ = *crsequence++;
        }
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    min_value[i] = seq.min_value[variable[i]];
    max_value[i] = seq.max_value[variable[i]];
    if (seq.marginal[variable[i]]) {
      marginal[i] = new Histogram(*(seq.marginal[variable[i]]));
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Selection de variables.
 *
 *  arguments : reference sur un objet Format_error, nombre de variables,
 *              indices des variables, flag pour conserver ou rejeter
 *              les variables selectionnees.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::select_variable(Format_error &error , int inb_variable ,
                                      int *ivariable , bool keep) const

{
  bool status = true , *selected_variable;
  register int i;
  int bnb_variable , *variable , *itype;
  Sequences *seq;


  seq = 0;
  error.init();

  if ((inb_variable < 1) || (inb_variable > nb_variable - 1)) {
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
    variable = ::select_variable(nb_variable , inb_variable , ivariable , keep);

    bnb_variable = (keep ? inb_variable : nb_variable - inb_variable);

    itype = new int[bnb_variable];
    for (i = 0;i < bnb_variable;i++) {
      itype[i] = type[variable[i]];
      if ((type[variable[i]] == AUXILIARY) && (variable[i - 1] != variable[i] - 1)) {
        itype[i] = REAL_VALUE;
      }
    }

    seq = new Sequences(nb_sequence , identifier , length , index_parameter_type ,
                        bnb_variable , itype);
    seq->select_variable(*this , variable);

    delete [] variable;
    delete [] itype;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Concatenation des variables d'objets Sequences.
 *
 *  arguments : reference sur un objet Format_error, nombre d'objets Sequences,
 *              pointeurs sur les objets Sequences, echantillon de reference pour les identificateurs.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::merge_variable(Format_error &error , int nb_sample ,
                                     const Sequences **iseq , int ref_sample) const

{
  bool status = true;
  register int i , j , k , m;
  int inb_variable , *iidentifier , *itype , *pisequence , *cisequence;
  double *prsequence , *crsequence;
  Sequences *seq;
  const Sequences **pseq;


  seq = 0;
  error.init();

  for (i = 0;i < nb_sample;i++) {
    if ((iseq[i]->index_parameter_type != IMPLICIT_TYPE) &&
        (iseq[i]->index_parameter_type != index_parameter_type)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_INDEX_PARAMETER_TYPE];

      if (index_parameter_type == IMPLICIT_TYPE) {
        error.update((error_message.str()).c_str());
      }
      else {
        error.correction_update((error_message.str()).c_str() , SEQ_index_parameter_word[index_parameter_type]);
      }
    }

    if (iseq[i]->nb_sequence != nb_sequence) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_NB_SEQUENCE];
      error.update((error_message.str()).c_str());
    }

    else {
      for (j = 0;j < nb_sequence;j++) {
        if (iseq[i]->length[j] != length[j]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                        << SEQ_label[SEQL_SEQUENCE] << " " << j + 1 << ": "
                        << SEQ_error[SEQR_SEQUENCE_LENGTH];
          error.update((error_message.str()).c_str());
        }

        else if ((iseq[i]->index_parameter_type != IMPLICIT_TYPE) &&
                 (iseq[i]->index_parameter_type == index_parameter_type)) {
          for (k = 0;k < (index_parameter_type == POSITION ? length[j] + 1 : length[j]);k++) {
            if (iseq[i]->index_parameter[j][k] != index_parameter[j][k]) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                            << SEQ_label[SEQL_SEQUENCE] << " " << j + 1 << ": "
                            << SEQ_label[SEQL_INDEX_PARAMETER] << " " << k << ": "
                            << SEQ_error[SEQR_INDEX_PARAMETER];
              error.update((error_message.str()).c_str());
            }
          }
        }
      }
    }
  }

  if ((ref_sample != I_DEFAULT) && ((ref_sample < 1) || (ref_sample > nb_sample + 1))) {
    status = false;
    error.update(STAT_error[STATR_SAMPLE_INDEX]);
  }

  if (status) {
    nb_sample++;
    pseq = new const Sequences*[nb_sample];

    pseq[0] = this;
    inb_variable = nb_variable;
    for (i = 1;i < nb_sample;i++) {
      pseq[i] = iseq[i - 1];
      inb_variable += iseq[i - 1]->nb_variable;
    }

    // calcul des identificateurs des sequences

    if (ref_sample == I_DEFAULT) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 1;j < nb_sample;j++) {
          if (pseq[j]->identifier[i] != pseq[0]->identifier[i]) {
            break;
          }
        }
        if (j < nb_sample) {
          break;
        }
      }

      if (i < nb_sequence) {
        iidentifier = 0;
      }
      else {
        iidentifier = pseq[0]->identifier;
      }
    }

    else {
      iidentifier = pseq[--ref_sample]->identifier;
    }

    itype = new int[inb_variable];
    inb_variable = 0;
    for (i = 0;i < nb_sample;i++) {
      for (j = 0;j < pseq[i]->nb_variable;j++) {
        itype[inb_variable] = pseq[i]->type[j];
        if ((inb_variable > 0) && (itype[inb_variable] == STATE)) {
          itype[inb_variable] = INT_VALUE;
        }
        inb_variable++;
      }
    }

    seq = new Sequences(nb_sequence , iidentifier , length , index_parameter_type ,
                        inb_variable , itype);
    delete [] itype;

    // copie des sequences

    if (hindex_parameter) {
      seq->hindex_parameter = new Histogram(*hindex_parameter);
    }
    if (index_interval) {
      seq->index_interval = new Histogram(*index_interval);
    }

    if (index_parameter) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < (index_parameter_type == POSITION ? length[i] + 1 : length[i]);j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }
      }
    }

    for (i = 0;i < nb_sequence;i++) {
      inb_variable = 0;
      for (j = 0;j < nb_sample;j++) {
        for (k = 0;k < pseq[j]->nb_variable;k++) {
          if ((seq->type[inb_variable] != REAL_VALUE) && (seq->type[inb_variable] != AUXILIARY)) {
            pisequence = seq->int_sequence[i][inb_variable++];
            cisequence = pseq[j]->int_sequence[i][k];
            for (m = 0;m < length[i];m++) {
              *pisequence++ = *cisequence++;
            }
          }

          else {
            prsequence = seq->real_sequence[i][inb_variable++];
            crsequence = pseq[j]->real_sequence[i][k];
            for (m = 0;m < length[i];m++) {
              *prsequence++ = *crsequence++;
            }
          }
        }
      }
    }

    inb_variable = 0;
    for (i = 0;i < nb_sample;i++) {
      for (j = 0;j < pseq[i]->nb_variable;j++) {
        seq->min_value[inb_variable] = pseq[i]->min_value[j];
        seq->max_value[inb_variable] = pseq[i]->max_value[j];
        if (pseq[i]->marginal[j]) {
          seq->marginal[inb_variable] = new Histogram(*(pseq[i]->marginal[j]));
        }
        inb_variable++;
      }
    }

    delete [] pseq;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Inversion du sens de parcours des sequences.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::reverse(Format_error &error) const

{
  Sequences *seq;


  error.init();

  if (index_parameter_type == TIME) {
    seq = 0;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER]);
  }

  else {
    seq = new Sequences(*this , 'c' , REVERSE);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Selection des sequences sur un critere de longueur.
 *
 *  arguments : reference sur un objet Format_error, bornes sur la longueur,
 *              flag pour conserver ou rejeter les sequences selectionnees.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::length_select(Format_error &error , int min_length ,
                                    int imax_length , bool keep) const

{
  bool status = true;
  register int i;
  int inb_sequence , *index;
  Sequences *seq;


  seq = 0;
  error.init();

  if ((min_length < 1) || (min_length >= hlength->nb_value) || (min_length > imax_length)) {
    status = false;
    error.update(SEQ_error[SEQR_MIN_SEQUENCE_LENGTH]);
  }
  if ((imax_length < hlength->offset) || (imax_length < min_length)) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_SEQUENCE_LENGTH]);
  }

  if (status) {

    // selection des sequences

    index = new int[nb_sequence];
    inb_sequence = 0;

    for (i = 0;i < nb_sequence;i++) {
      if ((length[i] >= min_length) && (length[i] <= imax_length)) {
        if (keep) {
          index[inb_sequence++] = i;
        }
      }

      else if (!keep) {
        index[inb_sequence++] = i;
      }
    }

    if (inb_sequence == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    // copie des sequences selectionnees

    if (status) {
      seq = new Sequences(*this , inb_sequence , index);
    }

    delete [] index;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Suppression des premieres/dernieres series d'une valeur donne.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              valeur, position ('b' : begin, 'e' : end),
 *              longueur maximum de la serie supprimee.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::remove_run(Format_error &error , int variable , int ivalue ,
                                 char position , int max_run_length) const

{
  bool status = true;
  register int i , j , k;
  int smax_length , inb_sequence , *iidentifier , *ilength , *index , *pindex_param ,
      *cindex_param , *pisequence , *cisequence;
  double *prsequence , *crsequence;
  Sequences *seq;


  seq = 0;
  error.init();

  if (index_parameter_type == POSITION) {
    error.update(SEQ_error[SEQR_INDEX_PARAMETER]);
    status = false;
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
      if (!marginal[variable]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << STAT_error[STATR_MARGINAL_HISTOGRAM];
        error.update((error_message.str()).c_str());
      }

      else if ((ivalue < marginal[variable]->offset) || (ivalue >= marginal[variable]->nb_value) ||
               (marginal[variable]->frequency[ivalue] == 0)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << STAT_label[STATL_VALUE] << " " << ivalue << " "
                      << SEQ_error[SEQR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if ((max_run_length < I_DEFAULT) || (max_run_length == 0)) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_RUN_LENGTH]);
  }

  if (status) {

    // calcul de la longueur des sequences

    iidentifier = new int[nb_sequence];
    ilength = new int[nb_sequence];
    index = new int[nb_sequence];
    inb_sequence = 0;

    for (i = 0;i < nb_sequence;i++) {
      if (max_run_length == I_DEFAULT) {
        smax_length = length[i];
      }
      else {
        smax_length = MIN(max_run_length , length[i]);
      }

      switch (position) {

      case 'b' : {
        cisequence = int_sequence[i][variable];
        for (j = 0;j < smax_length;j++) {
          if (*cisequence++ != ivalue) {
            break;
          }
        }
        break;
      }

      case 'e' : {
        cisequence = int_sequence[i][variable] + length[i];
        for (j = 0;j < smax_length;j++) {
          if (*--cisequence != ivalue) {
            break;
          }
        }
        break;
      }
      }

      if (j < length[i]) {
        iidentifier[inb_sequence] = identifier[i];
        ilength[inb_sequence] = length[i] - j;
        index[inb_sequence++] = i;
      }
    }

    // copie des sequences

    seq = new Sequences(inb_sequence , iidentifier , ilength , index_parameter_type ,
                        nb_variable , type);

    if (seq->index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        pindex_param = seq->index_parameter[i];

        switch (position) {
        case 'b' :
          cindex_param = index_parameter[index[i]] + length[index[i]] - seq->length[i];
          break;
        case 'e' :
          cindex_param = index_parameter[index[i]];
          break;
        }

        for (j = 0;j < seq->length[i];j++) {
          *pindex_param++ = *cindex_param++;
        }
      }

      seq->build_index_parameter_histogram();
      seq->index_interval_computation();
    }

    for (i = 0;i < seq->nb_sequence;i++) {
      for (j = 0;j < seq->nb_variable;j++) {
        if ((seq->type[j] != REAL_VALUE) && (seq->type[j] != AUXILIARY)) {
          pisequence = seq->int_sequence[i][j];

          switch (position) {
          case 'b' :
            cisequence = int_sequence[index[i]][j] + length[index[i]] - seq->length[i];
            break;
          case 'e' :
            cisequence = int_sequence[index[i]][j];
            break;
          }

          for (k = 0;k < seq->length[i];k++) {
            *pisequence++ = *cisequence++;
          }
        }

        else {
          prsequence = seq->real_sequence[i][j];

          switch (position) {
          case 'b' :
            crsequence = real_sequence[index[i]][j] + length[index[i]] - seq->length[i];
            break;
          case 'e' :
            crsequence = real_sequence[index[i]][j];
            break;
          }

          for (k = 0;k < seq->length[i];k++) {
            *prsequence++ = *crsequence++;
          }
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);
      seq->build_marginal_histogram(i);
    }

    delete [] iidentifier;
    delete [] ilength;
    delete [] index;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de sous-sequences.
 *
 *  arguments : reference sur un objet Format_error, parametres d'index minimum et
 *              maximum dans la sequence.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::index_parameter_extract(Format_error &error , int min_index_parameter ,
                                              int max_index_parameter) const

{
  bool status = true;
  register int i , j , k;
  int inb_sequence , *iidentifier , *ilength , *index , *first_index , *pindex_param ,
      *cindex_param , *pisequence , *cisequence;
  double *prsequence , *crsequence;
  Sequences *seq;


  seq = 0;
  error.init();

  if ((min_index_parameter < 0) || ((!index_parameter) && (min_index_parameter >= max_length)) ||
      ((max_index_parameter != I_DEFAULT) && (min_index_parameter > max_index_parameter))) {
    status = false;
    error.update(SEQ_error[SEQR_MIN_INDEX_PARAMETER]);
  }
  if ((max_index_parameter != I_DEFAULT) && ((max_index_parameter < 0) || ((!index_parameter) &&
        (max_index_parameter >= max_length)) || (max_index_parameter < min_index_parameter))) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_INDEX_PARAMETER]);
  }

  if (status) {

    // selection des sequences

    iidentifier = new int[nb_sequence];
    ilength = new int[nb_sequence];
    index = new int[nb_sequence];
    inb_sequence = 0;

    // parametre d'index explicite

    if (index_parameter) {
      first_index = new int[nb_sequence];

      if (max_index_parameter == I_DEFAULT) {
        for (i = 0;i < nb_sequence;i++) {
          if (index_parameter[i][length[i] - 1] >= min_index_parameter) {
            iidentifier[inb_sequence] = identifier[i];

            for (j = 0;j < length[i];j++) {
              if (index_parameter[i][j] >= min_index_parameter) {
                break;
              }
            }
            first_index[inb_sequence] = j;
            ilength[inb_sequence] = length[i] - j;

            index[inb_sequence++] = i;
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          if ((index_parameter[i][0] <= min_index_parameter) &&
              (index_parameter[i][length[i] - 1] >= max_index_parameter)) {
            iidentifier[inb_sequence] = identifier[i];

            for (j = 0;j < length[i];j++) {
              if (index_parameter[i][j] >= min_index_parameter) {
                break;
              }
            }
            first_index[inb_sequence] = j;
            ilength[inb_sequence] = length[i] - j;

            for (j = length[i] - 1;j >= 0;j--) {
              if (index_parameter[i][j] <= max_index_parameter) {
                break;
              }
            }
            ilength[inb_sequence] -= (length[i] - 1 - j);

            index[inb_sequence++] = i;
          }
        }
      }
    }

    // parametre d'index implicite

    else {
      if (max_index_parameter == I_DEFAULT) {
        for (i = 0;i < nb_sequence;i++) {
          if (length[i] > min_index_parameter) {
            iidentifier[inb_sequence] = identifier[i];
            ilength[inb_sequence] = length[i] - min_index_parameter;
            index[inb_sequence++] = i;
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          if (length[i] > max_index_parameter) {
            iidentifier[inb_sequence] = identifier[i];
            ilength[inb_sequence] = max_index_parameter - min_index_parameter + 1;
            index[inb_sequence++] = i;
          }
        }
      }
    }

    if (inb_sequence == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    // extraction des sous-sequences

    if (status) {
      seq = new Sequences(inb_sequence , iidentifier , ilength , index_parameter_type ,
                          nb_variable , type);

      // parametre d'index explicite

      if (seq->index_parameter) {
        for (i = 0;i < seq->nb_sequence;i++) {
          pindex_param = seq->index_parameter[i];
          cindex_param = index_parameter[index[i]] + first_index[i];
          for (k = 0;k < seq->length[i];k++) {
            *pindex_param++ = *cindex_param++;
          }

          if (seq->index_parameter_type == POSITION) {
            if (max_index_parameter == I_DEFAULT) {
              *pindex_param = *cindex_param;
            }
            else {
              *pindex_param = max_index_parameter;
            }
          }
        }
      }

      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < seq->nb_variable;j++) {
          if ((seq->type[j] != REAL_VALUE) && (seq->type[j] != AUXILIARY)) {
            pisequence = seq->int_sequence[i][j];
            if (seq->index_parameter) {
              cisequence = int_sequence[index[i]][j] + first_index[i];
            }
            else {
              cisequence = int_sequence[index[i]][j] + min_index_parameter;
            }

            for (k = 0;k < seq->length[i];k++) {
              *pisequence++ = *cisequence++;
            }
          }

          else {
            prsequence = seq->real_sequence[i][j];
            if (seq->index_parameter) {
              crsequence = real_sequence[index[i]][j] + first_index[i];
            }
            else {
              crsequence = real_sequence[index[i]][j] + min_index_parameter;
            }

            for (k = 0;k < seq->length[i];k++) {
              *prsequence++ = *crsequence++;
            }
          }
        }
      }

      if (hindex_parameter) {
        seq->build_index_parameter_histogram();
      }
      if (index_interval) {
        seq->index_interval_computation();
      }

      for (i = 0;i < seq->nb_variable;i++) {
        seq->min_value_computation(i);
        seq->max_value_computation(i);
        seq->build_marginal_histogram(i);
      }

      delete [] iidentifier;
      delete [] ilength;
      delete [] index;

      if (index_parameter) {
        delete [] first_index;
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Extraction par segmentation d'un objet Sequences.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              nombre de valeurs, valeurs, flag zones correspondant aux valeurs
 *              extraites/pas extraites.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::segmentation_extract(Format_error &error , int variable ,
                                           int nb_value , int *ivalue , bool keep) const

{
  bool status = true;
  register int i , j , k , m;
  int nb_present_value , nb_selected_value , *pfrequency , *selected_value , *itype ,
      *pindex_param , *cindex_param , *pisequence , *cisequence , *zone_length , **zone_begin;
  double *prsequence , *crsequence;
  Sequences *seq;


  seq = 0;
  error.init();

  if (index_parameter_type == POSITION) {
    error.update(SEQ_error[SEQR_INDEX_PARAMETER]);
    status = false;
  }
  if (nb_variable == 1) {
    status = false;
    error.update(STAT_error[STATR_NB_VARIABLE]);
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
      if (!marginal[variable]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << STAT_error[STATR_MARGINAL_HISTOGRAM];
        error.update((error_message.str()).c_str());
      }

      else {
        pfrequency = marginal[variable]->frequency + marginal[variable]->offset;
        nb_present_value = 0;
        for (i = marginal[variable]->offset;i < marginal[variable]->nb_value;i++) {
          if (*pfrequency > 0) {
            nb_present_value++;
          }
          pfrequency++;
        }

        if ((nb_value < 1) || (nb_value > (keep ? nb_present_value : nb_present_value - 1))) {
          status = false;
          error.update(SEQ_error[SEQR_NB_SELECTED_VALUE]);
        }

        else {
          selected_value = new int[marginal[variable]->nb_value];
          for (i = marginal[variable]->offset;i < marginal[variable]->nb_value;i++) {
            selected_value[i] = false;
          }

          for (i = 0;i < nb_value;i++) {
            if ((ivalue[i] < marginal[variable]->offset) || (ivalue[i] >= marginal[variable]->nb_value) ||
                (marginal[variable]->frequency[ivalue[i]] == 0)) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                            << STAT_label[STATL_VALUE] << " " << ivalue[i] << " "
                            << SEQ_error[SEQR_NOT_PRESENT];
              error.update((error_message.str()).c_str());
            }

            else if (selected_value[ivalue[i]]) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_VALUE] << " " << ivalue[i] << " "
                            << STAT_error[STATR_ALREADY_SELECTED];
              error.update((error_message.str()).c_str());
            }
            else {
              selected_value[ivalue[i]] = true;
            }
          }

          delete [] selected_value;
        }
      }
    }

    if ((variable + 1 < nb_variable) && ((type[variable + 1] != INT_VALUE) &&
         (type[variable + 1] != STATE) && (type[variable + 1] != REAL_VALUE))) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (status) {
    itype = new int[nb_variable - 1];
    for (i = 0;i < variable;i++) {
      itype[i] = type[i];
    }
    for (i = variable + 1;i < nb_variable;i++) {
      itype[i - 1] = type[i];
    }

    switch (keep) {

    case false : {
      nb_selected_value = nb_present_value - nb_value;
      selected_value = new int[nb_selected_value];

      pfrequency = marginal[variable]->frequency + marginal[variable]->offset;
      i = 0;

      for (j = marginal[variable]->offset;j < marginal[variable]->nb_value;j++) {
        if (*pfrequency > 0) {
          for (k = 0;k < nb_value;k++) {
            if (ivalue[k] == j) {
              break;
            }
          }

          if (k == nb_value) {
            selected_value[i++] = j;
          }
        }

        pfrequency++;
      }
      break;
    }

    case true : {
      nb_selected_value = nb_value;
      selected_value = ivalue;
      break;
    }
    }

    // initialisations

    zone_length = new int[cumul_length];

    zone_begin = new int*[cumul_length];
    for (i = 0;i < cumul_length;i++) {
      zone_begin[i] = 0;
    }

    // recherche des sequences a extraire

    i = -1;

    for (j = 0;j < nb_sequence;j++) {
      pisequence = int_sequence[j][variable];

      for (k = 0;k < nb_selected_value;k++) {
        if (*pisequence == selected_value[k]) {
          break;
        }
      }

      if (k < nb_selected_value) {
        i++;
        zone_begin[i] = new int[2];

        zone_begin[i][0] = j;
        zone_begin[i][1] = 0;
        zone_length[i] = 1;
      }

      for (k = 1;k < length[j];k++) {
        pisequence++;

        for (m = 0;m < nb_selected_value;m++) {
          if (*pisequence == selected_value[m]) {
            break;
          }
        }

        if (m  < nb_selected_value) {
          if (*pisequence != *(pisequence - 1)) {
            i++;
            zone_begin[i] = new int[2];

            zone_begin[i][0] = j;
            zone_begin[i][1] = k;
            zone_length[i] = 1;
          }
          else {
            zone_length[i]++;
          }
        }
      }
    }

    // creation de l'objet Sequences

    seq = new Sequences(i + 1 , 0 , zone_length , index_parameter_type ,
                        nb_variable - 1 , itype);

    // copie des sequences

    if (index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        pindex_param = seq->index_parameter[i];
        cindex_param = index_parameter[zone_begin[i][0]] + zone_begin[i][1];
        for (j = 0;j < seq->length[i];j++) {
          *pindex_param++ = *cindex_param++;
        }
      }
    }

    if (hindex_parameter) {
      seq->build_index_parameter_histogram();
    }
    if (index_interval) {
      seq->index_interval_computation();
    }

    for (i = 0;i < seq->nb_sequence;i++) {
      j = 0;
      for (k = 0;k < nb_variable;k++) {
        if (k != variable) {
          if ((type[k] != REAL_VALUE) && (type[k] != AUXILIARY)) {
            pisequence = seq->int_sequence[i][j];
            cisequence = int_sequence[zone_begin[i][0]][k] + zone_begin[i][1];
            for (m = 0;m < seq->length[i];m++) {
              *pisequence++ = *cisequence++;
            }
          }

          else {
            prsequence = seq->real_sequence[i][j];
            crsequence = real_sequence[zone_begin[i][0]][k] + zone_begin[i][1];
            for (m = 0;m < seq->length[i];m++) {
              *prsequence++ = *crsequence++;
            }
          }

          j++;
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);
      seq->build_marginal_histogram(i);
    }

    if (!keep) {
      delete [] selected_value;
    }

    delete [] zone_length;
    for (i = 0;i < seq->nb_sequence;i++) {
      delete [] zone_begin[i];
    }
    delete [] zone_begin;

    delete [] itype;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Cumul des valeurs successives des sequences.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::cumulate(Format_error &error , int variable) const

{
  bool status = true;
  register int i , j , k , m;
  int inb_variable , *itype , *pisequence , *cisequence;
  double *prsequence , *crsequence;
  Sequences *seq;


  seq = 0;
  error.init();

  if (variable != I_DEFAULT) {
    if ((variable < 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }
    else {
      variable--;
    }
  }

  for (i = 0;i < nb_variable;i++) {
    if (((variable == I_DEFAULT) || (variable == i)) &&
        ((type[i] != INT_VALUE) && (type[i] != STATE) && (type[i] != REAL_VALUE))) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (status) {
    if (variable == I_DEFAULT) {
      inb_variable = nb_variable;
      itype = type;
    }
    else {
      inb_variable = 1;
      itype = type + variable;
    }

    seq = new Sequences(nb_sequence , identifier , length , index_parameter_type ,
                        inb_variable , itype);

    // copie des parametres d'index

    if (hindex_parameter) {
      seq->hindex_parameter = new Histogram(*hindex_parameter);
    }
    if (index_interval) {
      seq->index_interval = new Histogram(*index_interval);
    }

    if (index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < (seq->index_parameter_type == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }
      }
    }

   // cumul des valeurs

    for (i = 0;i < nb_sequence;i++) {
      j = 0;
      for (k = 0;k < nb_variable;k++) {
        if ((variable == I_DEFAULT) || (variable == k)) {
          if (type[k] != REAL_VALUE) {
            pisequence = seq->int_sequence[i][j++];
            cisequence = int_sequence[i][k];

            *pisequence++ = *cisequence++;
            for (m = 1;m < length[i];m++) {
              *pisequence = *(pisequence - 1) + *cisequence++;
              pisequence++;
            }
          }

          else {
            prsequence = seq->real_sequence[i][j++];
            crsequence = real_sequence[i][k];

            *prsequence++ = *crsequence++;
            for (m = 1;m < length[i];m++) {
              *prsequence = *(prsequence - 1) + *crsequence++;
              prsequence++;
            }
          }
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);
      seq->build_marginal_histogram(i);
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Differenciation au 1er ordre des sequences.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              premier element de la sequence garde ou pas.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::difference(Format_error &error , int variable ,
                                 bool first_element) const

{
  bool status = true;
  register int i , j , k , m;
  int inb_variable , *ilength , *itype , *pindex_param , *cindex_param ,
      *pisequence , *cisequence;
  double *prsequence , *crsequence;
  Sequences *seq;


  seq = 0;
  error.init();

  if (index_parameter_type == POSITION) {
    status = increasing_index_parameter_checking(error , true , SEQ_label[SEQL_SEQUENCE]);
  }

  if (variable != I_DEFAULT) {
    if ((variable < 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }
    else {
      variable--;
    }
  }

  for (i = 0;i < nb_variable;i++) {
    if (((variable == I_DEFAULT) || (variable == i)) &&
        ((type[i] != INT_VALUE) && (type[i] != STATE) && (type[i] != REAL_VALUE))) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if ((!first_element) && (hlength->offset < 2)) {
    status = false;
    ostringstream correction_message;
    correction_message << STAT_error[STATR_GREATER_THAN] << " " << 1;
    error.correction_update(SEQ_error[SEQR_MIN_SEQUENCE_LENGTH] , (correction_message.str()).c_str());
  }

  if (status) {
    if (first_element) {
      ilength = length;
    }
    else {
      ilength = new int[nb_sequence];
      for (i = 0;i < nb_sequence;i++) {
        ilength[i] = length[i] - 1;
      }
    }

    if (variable == I_DEFAULT) {
      inb_variable = nb_variable;
    }
    else {
      inb_variable = 1;
    }

    if ((index_parameter_type != IMPLICIT_TYPE) && ((index_interval->mean != 1.) ||
         (index_interval->variance > 0.))) {
      itype = new int[inb_variable];
      for (i = 0;i < inb_variable;i++) {
        itype[i] = REAL_VALUE;
      }
    }

    else {
      if (variable == I_DEFAULT) {
        itype = type;
      }
      else {
        itype = type + variable;
      }
    }

    seq = new Sequences(nb_sequence , identifier , ilength , index_parameter_type ,
                        inb_variable , itype);

    if (!first_element) {
      delete [] ilength;
    }
    if ((index_parameter_type != IMPLICIT_TYPE) && ((index_interval->mean != 1.) ||
         (index_interval->variance > 0.))) {
      delete [] itype;
    }

    // copie des parametres index

    if (index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        pindex_param = seq->index_parameter[i];
        if (first_element) {
          cindex_param = index_parameter[i];
        }
        else {
          cindex_param = index_parameter[i] + 1;
        }

        for (j = 0;j < (seq->index_parameter_type == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
          *pindex_param++ = *cindex_param++;
        }
      }
    }

    if (first_element) {
      if (hindex_parameter) {
        seq->hindex_parameter = new Histogram(*hindex_parameter);
      }
      if (index_interval) {
        seq->index_interval = new Histogram(*index_interval);
      }
    }

    else {
      seq->build_index_parameter_histogram();
      if (index_interval) {
        seq->index_interval_computation();
      }
    }

    // differenciation des sequences

    if ((index_parameter_type != IMPLICIT_TYPE) && ((index_interval->mean != 1.) ||
         (index_interval->variance > 0.))) {
      for (i = 0;i < nb_sequence;i++) {
        j = 0;
        for (k = 0;k < nb_variable;k++) {
          if ((variable == I_DEFAULT) || (variable == k)) {
            prsequence = seq->real_sequence[i][j];
            cindex_param = index_parameter[i];

            if (type[k] != REAL_VALUE) {
              cisequence = int_sequence[i][k];

              if (first_element) {
                *prsequence++ = (double)*cisequence / (double)*cindex_param;
              }
              for (m = 0;m < length[i] - 1;m++) {
                *prsequence++ = (double)(*(cisequence + 1) - *cisequence) /
                                (double)(*(cindex_param + 1) - *cindex_param);
                cindex_param++;
                cisequence++;
              }
            }

            else {
              crsequence = real_sequence[i][k];

              if (first_element) {
                *prsequence++ = *crsequence / *cindex_param;
              }
              for (m = 0;m < length[i] - 1;m++) {
                *prsequence++ = (*(crsequence + 1) - *crsequence) /
                                (*(cindex_param + 1) - *cindex_param);
                cindex_param++;
                crsequence++;
              }
            }

            j++;
          }
        }
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        j = 0;
        for (k = 0;k < nb_variable;k++) {
          if ((variable == I_DEFAULT) || (variable == k)) {
            if (type[k] != REAL_VALUE) {
              pisequence = seq->int_sequence[i][j];
              cisequence = int_sequence[i][k];

              if (first_element) {
                *pisequence++ = *cisequence;
              }
              for (m = 0;m < length[i] - 1;m++) {
                *pisequence++ = *(cisequence + 1) - *cisequence;
                cisequence++;
              }
            }

            else {
              prsequence = seq->real_sequence[i][j];
              crsequence = real_sequence[i][k];

              if (first_element) {
                *prsequence++ = *crsequence;
              }
              for (m = 0;m < length[i] - 1;m++) {
                *prsequence++ = *(crsequence + 1) - *crsequence;
                crsequence++;
              }
            }

            j++;
          }
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);
      seq->build_marginal_histogram(i);
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Filtrage de type moyenne mobile des sequences.
 *
 *  arguments : reference sur un objet Format_error, demi-largeur du filtre,
 *              filtre, indice de la variable, debut/fin garde ou pas,
 *              tendance ou residus (par soustraction ou par division).
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::moving_average(Format_error &error , int nb_point , double *filter ,
                                     int variable , bool begin_end , int output) const

{
  bool status = true;
  register int i , j , k , m , n;
  int inb_variable , *ilength , *itype , *pindex_param , *cindex_param ,
      *pisequence , *cisequence;
  double *prsequence , *crsequence , *ppoint;
  Sequences *seq;


  seq = 0;
  error.init();

  if ((index_interval) && (index_interval->variance > 0.)) {
    status = false;
    error.update(SEQ_error[SEQR_UNEQUAL_INDEX_INTERVALS]);
  }

  if (variable != I_DEFAULT) {
    if ((variable < 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }
    else {
      variable--;
    }
  }

  for (i = 0;i < nb_variable;i++) {
    if (((variable == I_DEFAULT) || (variable == i)) &&
        ((type[i] != INT_VALUE) && (type[i] != STATE) && (type[i] != REAL_VALUE))) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if ((!begin_end) && (hlength->offset < 2 * nb_point + 1)) {
    status = false;
    ostringstream correction_message;
    correction_message << STAT_error[STATR_GREATER_THAN] << " " << 2 * nb_point;
    error.correction_update(SEQ_error[SEQR_MIN_SEQUENCE_LENGTH] , (correction_message.str()).c_str());
  }

  if (status) {
    if (begin_end) {
      ilength = length;
    }
    else {
      ilength = new int[nb_sequence];
      for (i = 0;i < nb_sequence;i++) {
        ilength[i] = length[i] - 2 * nb_point;
      }
    }

    if (variable == I_DEFAULT) {
      inb_variable = nb_variable;
    }
    else {
      inb_variable = 1;
    }

    if (output == SEQUENCE) {
      inb_variable *= 2;
    }

    itype = new int[inb_variable];

    if (output == SEQUENCE) {
      i = 0;
      for (j = 0;j < nb_variable;j++) {
        if ((variable == I_DEFAULT) || (variable == j)) {
          itype[i++] = type[j];
          itype[i++] = AUXILIARY;
        }
      }
    }

    else {
      for (i = 0;i < inb_variable;i++) {
        itype[i] = REAL_VALUE;
      }
    }

    seq = new Sequences(nb_sequence , identifier , ilength , index_parameter_type ,
                        inb_variable , itype);

    if (!begin_end) {
      delete [] ilength;
    }
    delete [] itype;

    // copie des parametres index

    if (index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        pindex_param = seq->index_parameter[i];
        if (begin_end) {
          cindex_param = index_parameter[i];
        }
        else {
          cindex_param = index_parameter[i] + nb_point;
        }

        for (j = 0;j < (seq->index_parameter_type == POSITION ? seq->length[i] + 1 : seq->length[i]);j++) {
          *pindex_param++ = *cindex_param++;
        }
      }
    }

    if (begin_end) {
      if (hindex_parameter) {
        seq->hindex_parameter = new Histogram(*hindex_parameter);
      }
      if (index_interval) {
        seq->index_interval = new Histogram(*index_interval);
      }
    }

    else {
      seq->build_index_parameter_histogram();
      if (index_interval) {
        seq->index_interval_computation();
      }
    }

    // filtrage de type moyenne mobile

    for (i = 0;i < nb_sequence;i++) {
      j = 0;
      for (k = 0;k < nb_variable;k++) {
        if ((variable == I_DEFAULT) || (variable == k)) {
          prsequence = seq->real_sequence[i][output == SEQUENCE ? j + 1 : j];

          if (type[k] != REAL_VALUE) {
            if (begin_end) {
              for (m = 0;m < MIN(nb_point , length[i]);m++) {
                cisequence = int_sequence[i][k];
                ppoint = filter;
                *prsequence = 0.;
                for (n = 0;n < 2 * nb_point + 1;n++) {
                  *prsequence += *cisequence * *ppoint++;
                  if ((m - nb_point + n >= 0) && (m - nb_point + n < length[i] - 1)) {
                    cisequence++;
                  }
                }
                prsequence++;
              }
            }

            for (m = nb_point;m < length[i] - nb_point;m++) {
              cisequence = int_sequence[i][k] + m - nb_point;
              ppoint = filter;
              *prsequence = 0.;
              for (n = 0;n < 2 * nb_point + 1;n++) {
                *prsequence += *cisequence++ * *ppoint++;
              }
              prsequence++;
            }

            if (begin_end) {
              for (m = MAX(length[i] - nb_point , nb_point);m < length[i];m++) {
                cisequence = int_sequence[i][k] + m - nb_point;
                ppoint = filter;
                *prsequence = 0.;
                for (n = 0;n < 2 * nb_point + 1;n++) {
                  *prsequence += *cisequence * *ppoint++;
                  if ((m - nb_point + n >= 0) && (m - nb_point + n < length[i] - 1)) {
                    cisequence++;
                  }
                }
                prsequence++;
              }
            }

            if (begin_end) {
              cisequence = int_sequence[i][k];
            }
            else {
              cisequence = int_sequence[i][k] + nb_point;
            }

            switch (output) {

            case SEQUENCE : {
              pisequence = seq->int_sequence[i][j];
              for (m = 0;m < seq->length[i];m++) {
                *pisequence++ = *cisequence++;
              }
              break;
            }

            case SUBTRACTION_RESIDUAL : {
              prsequence = seq->real_sequence[i][j];
              for (m = 0;m < seq->length[i];m++) {
                *prsequence = *cisequence++ - *prsequence;
                prsequence++;
              }
              break;
            }

            case DIVISION_RESIDUAL : {
              prsequence = seq->real_sequence[i][j];
              for (m = 0;m < seq->length[i];m++) {
                if (*prsequence != 0.) {
                  *prsequence = *cisequence / *prsequence;
                }
                prsequence++;
                cisequence++;
              }
              break;
            }
            }
          }

          else {
            if (begin_end) {
              for (m = 0;m < MIN(nb_point , length[i]);m++) {
                crsequence = real_sequence[i][k];
                ppoint = filter;
                *prsequence = 0.;
                for (n = 0;n < 2 * nb_point + 1;n++) {
                  *prsequence += *crsequence * *ppoint++;
                  if ((m - nb_point + n >= 0) && (m - nb_point + n < length[i] - 1)) {
                    crsequence++;
                  }
                }
                prsequence++;
              }
            }

            for (m = nb_point;m < length[i] - nb_point;m++) {
              crsequence = real_sequence[i][k] + m - nb_point;
              ppoint = filter;
              *prsequence = 0.;
              for (n = 0;n < 2 * nb_point + 1;n++) {
                *prsequence += *crsequence++ * *ppoint++;
              }
              prsequence++;
            }

            if (begin_end) {
              for (m = MAX(length[i] - nb_point , nb_point);m < length[i];m++) {
                crsequence = real_sequence[i][k] + m - nb_point;
                ppoint = filter;
                *prsequence = 0.;
                for (n = 0;n < 2 * nb_point + 1;n++) {
                  *prsequence += *crsequence * *ppoint++;
                  if ((m - nb_point + n >= 0) && (m - nb_point + n < length[i] - 1)) {
                    crsequence++;
                  }
                }
                prsequence++;
              }
            }

            prsequence = seq->real_sequence[i][j];
            if (begin_end) {
              crsequence = real_sequence[i][k];
            }
            else {
              crsequence = real_sequence[i][k] + nb_point;
            }

            switch (output) {

            case SEQUENCE : {
              for (m = 0;m < seq->length[i];m++) {
                *prsequence++ = *crsequence++;
              }
              break;
            }

            case SUBTRACTION_RESIDUAL : {
              for (m = 0;m < seq->length[i];m++) {
                *prsequence = *crsequence++ - *prsequence;
                prsequence++;
              }
              break;
            }

            case DIVISION_RESIDUAL : {
              for (m = 0;m < seq->length[i];m++) {
                if (*prsequence != 0.) {
                  *prsequence = *crsequence / *prsequence;
                }
                prsequence++;
                crsequence++;
              }
              break;
            }
            }
          }

          if (output == SEQUENCE) {
            j++;
          }
          j++;
        }
      }
    }

    if (output == SEQUENCE) {
      i = 0;

      if (begin_end) {
        for (j = 0;j < nb_variable;j++) {
          if ((variable == I_DEFAULT) || (variable == j)) {
            seq->min_value[i] = min_value[j];
            seq->max_value[i] = max_value[j];
            if (marginal[j]) {
              seq->marginal[i] = new Histogram(*marginal[j]);
            }
            i++;
            seq->min_value_computation(i);
            seq->max_value_computation(i);
            i++;
          }
        }
      }

      else {
        for (j = 0;j < nb_variable;j++) {
          if ((variable == I_DEFAULT) || (variable == j)) {
            seq->min_value_computation(i);
            seq->max_value_computation(i);
            seq->build_marginal_histogram(i);
            i++;
            seq->min_value_computation(i);
            seq->max_value_computation(i);
            i++;
          }
        }
      }
    }

    else {
      for (i = 0;i < seq->nb_variable;i++) {
        seq->min_value_computation(i);
        seq->max_value_computation(i);
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Filtrage de type moyenne mobile des sequences.
 *
 *  arguments : reference sur un objet Format_error, loi symmetrique,
 *              indice de la variable, debut/fin supprime ou pas,
 *              tendance ou residus (par soustraction ou par division).
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::moving_average(Format_error &error , const Distribution &dist ,
                                     int variable , bool begin_end , int output) const

{
  bool status = true;
  Sequences *seq;


  seq = 0;
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
    seq = moving_average(error , dist.nb_value / 2 , dist.mass , variable ,
                         begin_end , output);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des sequences de moyennes et d'ecart-types.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              tailles d'echantillons, flag calcul des ecart-types,
 *              sortie (sequences, residus ou residus standardisees).
 *
 *--------------------------------------------------------------*/

bool Sequences::pointwise_average_ascii_print(Format_error &error , const char *path ,
                                              int *size , bool standard_deviation ,
                                              int output) const

{
  bool status;
  register int i , j , k , m;
  int buff , inb_sequence , *width;
  long old_adjust;
  double standard_normal_value , half_confidence_interval , *t_value;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    old_adjust = out_file.flags(ios::adjustfield);

    if (standard_deviation) {

#     ifdef MESSAGE
      standard_normal_value = standard_normal_value_computation(0.025);
      cout << "\nTEST: " << standard_normal_value;
#     endif

      t_value = new double[length[0]];
      for (i = 0;i < length[0];i++) {
        if (size[i] > 1) {
          t_value[i] = t_value_computation(false , size[i] - 1 , 0.05);
        }
      }

#     ifdef MESSAGE
      cout << " | " << t_value[0] << " " << size[0] << endl;
#     endif

    }

    else {
      t_value = 0;
    }

    // calcul des largeurs des colonnes

    inb_sequence = nb_sequence;
    if (standard_deviation) {
      inb_sequence--;
    }

    width = new int[2 * nb_variable + 2];

    for (i = 0;i < nb_variable;i++) {
      width[i] = column_width(length[0] , real_sequence[0][i]);
      if (standard_deviation) {
        buff = column_width(length[nb_sequence - 1] , real_sequence[nb_sequence - 1][i]);
        if (buff > width[i]) {
          width[i] = buff;
        }
      }
      width[i] += ASCII_SPACE;
    }

    if (index_parameter) {
      width[nb_variable] = column_width(hindex_parameter->nb_value - 1);
    }
    else {
      width[nb_variable] = column_width(max_length);
    }

    width[nb_variable + 1] = column_width(nb_sequence) + ASCII_SPACE;

    if (nb_sequence < POINTWISE_AVERAGE_NB_SEQUENCE) {
      for (i = 0;i < nb_variable;i++) {
        width[nb_variable + i + 2] = 0;
        for (j = 1;j < inb_sequence;j++) {
          buff = column_width(length[j] , real_sequence[j][i]);
          if (buff > width[nb_variable + i + 2]) {
            width[nb_variable + i + 2] = buff;
          }
        }
        width[nb_variable + i + 2] += ASCII_SPACE;
      }
    }

    out_file.setf(ios::right , ios::adjustfield);

    switch (output) {
    case SUBTRACTION_RESIDUAL :
      out_file << STAT_label[STATL_RESIDUAL] << "\n" << endl;
      break;
    case STANDARDIZED_RESIDUAL :
      out_file << STAT_label[STATL_STANDARDIZED_RESIDUAL] << "\n" << endl;
      break;
    }

    for (i = 0;i < nb_variable;i++) {
      out_file << STAT_label[STATL_VARIABLE] << " " << i + 1 << endl;

      if (index_parameter_type == TIME) {
        out_file << "\n" << SEQ_label[SEQL_TIME];
      }
      else {
        out_file << "\n" << SEQ_label[SEQL_INDEX];
      }

      out_file << "   " << STAT_label[STATL_MEAN];
      if (standard_deviation) {
        out_file << "   " << STAT_label[STATL_MEAN_CONFIDENCE_INTERVAL]
                 << "   " << STAT_label[STATL_STANDARD_DEVIATION];
      }
      out_file << "   " << STAT_label[STATL_SAMPLE_SIZE];

      if (nb_sequence < POINTWISE_AVERAGE_NB_SEQUENCE) {
        out_file << "   ";
        for (j = 1;j < inb_sequence;j++) {
          out_file << "   " << SEQ_label[SEQL_SEQUENCE] << " " << identifier[j];
        }
      }
      out_file << endl;

      for (j = 0;j < length[0];j++) {
        out_file << setw(width[nb_variable]) << (index_parameter ? index_parameter[0][j] : j) << "  "
                 << setw(width[i]) << real_sequence[0][i][j];

        if (standard_deviation) {
          if (size[j] > 1) {
//            half_confidence_interval = standard_normal_value * real_sequence[nb_sequence - 1][i][j] / sqrt((double)size[j]);
            half_confidence_interval = t_value[j] * real_sequence[nb_sequence - 1][i][j] / sqrt((double)size[j]);
            out_file << setw(width[i]) << real_sequence[0][i][j] - half_confidence_interval
                     << setw(width[i]) << real_sequence[0][i][j] + half_confidence_interval;
          }

          else {
            out_file << setw(width[i]) << " "
                     << setw(width[i]) << " ";
          }

          out_file << setw(width[i]) << real_sequence[nb_sequence - 1][i][j];
        }

        out_file << setw(width[nb_variable + 1]) << size[j];
 
        if (inb_sequence - 1 < POINTWISE_AVERAGE_NB_SEQUENCE) {
          out_file << "   ";

          if (index_parameter) {
            for (k = 1;k < inb_sequence;k++) {
              for (m = 0;m < length[k];m++) {
                if (index_parameter[k][m] == index_parameter[0][j]) {
                  out_file << setw(width[nb_variable + i + 2]) << real_sequence[k][i][m];
                  break;
                }
              }

              if (m == length[k]) {
                out_file << setw(width[nb_variable + i + 2]) << " ";
              }
            }
          }

          else {
            for (k = 1;k < inb_sequence;k++) {
              if (j < length[k]) {
                out_file << setw(width[nb_variable + i + 2]) << real_sequence[k][i][j];
              }
              else {
                out_file << setw(width[nb_variable + i + 2]) << " ";
              }
            }
          }
        }

        out_file << endl;
      }

      out_file << endl;
    }

    delete [] width;
    delete [] t_value;

    out_file.setf((FMTFLAGS)old_adjust , ios::adjustfield);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des sequences de moyennes et d'ecart-types au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              tailles d'echantillons, flag calcul des ecart-types,
 *              sortie (sequences, residus ou residus standardisees).
 *
 *--------------------------------------------------------------*/

bool Sequences::pointwise_average_spreadsheet_print(Format_error &error , const char *path ,
                                                    int *size , bool standard_deviation ,
                                                    int output) const

{
  bool status;
  register int i , j , k , m;
  int inb_sequence;
  double standard_normal_value , half_confidence_interval , *t_value;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    if (standard_deviation) {
//      standard_normal_value = standard_normal_value_computation(0.025);

      t_value = new double[length[0]];
      for (i = 0;i < length[0];i++) {
        if (size[i] > 1) {
          t_value[i] = t_value_computation(false , size[i] - 1 , 0.05);
        }
      }
    }

    else {
      t_value = 0;
    }

    switch (output) {
    case SUBTRACTION_RESIDUAL :
      out_file << STAT_label[STATL_RESIDUAL] << "\n" << endl;
      break;
    case STANDARDIZED_RESIDUAL :
      out_file << STAT_label[STATL_STANDARDIZED_RESIDUAL] << "\n" << endl;
      break;
    }

    for (i = 0;i < nb_variable;i++) {
      out_file << STAT_label[STATL_VARIABLE] << "\t" << i + 1 << endl;

      if (index_parameter_type == TIME) {
        out_file << "\n" << SEQ_label[SEQL_TIME];
      }
      else {
        out_file << "\n" << SEQ_label[SEQL_INDEX];
      }

      out_file << "\t" << STAT_label[STATL_MEAN];
      if (standard_deviation) {
        out_file << "\t" << STAT_label[STATL_MEAN_CONFIDENCE_INTERVAL]
                 << "\t\t" << STAT_label[STATL_STANDARD_DEVIATION];
      }
      out_file << "\t" << STAT_label[STATL_SAMPLE_SIZE];

      if (nb_sequence < POINTWISE_AVERAGE_NB_SEQUENCE) {
        inb_sequence = nb_sequence;
        if (standard_deviation) {
          inb_sequence--;
        }

        out_file << "\t";
        for (j = 1;j < inb_sequence;j++) {
          out_file << "\t" << SEQ_label[SEQL_SEQUENCE] << " " << identifier[j];
        }
      }
      out_file << endl;

      for (j = 0;j < length[0];j++) {
        out_file << (index_parameter ? index_parameter[0][j] : j)
                 << "\t" << real_sequence[0][i][j];

        if (standard_deviation) {
          if (size[j] > 1) {
//            half_confidence_interval = standard_normal_value * real_sequence[nb_sequence - 1][i][j] / sqrt((double)size[j]);
            half_confidence_interval = t_value[j] * real_sequence[nb_sequence - 1][i][j] / sqrt((double)size[j]);
            out_file << "\t" << real_sequence[0][i][j] - half_confidence_interval
                     << "\t" << real_sequence[0][i][j] + half_confidence_interval;
          }

          else {
            out_file << "\t\t";
          }

          out_file << "\t" << real_sequence[nb_sequence - 1][i][j];
        }

        out_file << "\t" << size[j];
 
        if (inb_sequence - 1 < POINTWISE_AVERAGE_NB_SEQUENCE) {
          out_file << "\t";

          if (index_parameter) {
            for (k = 1;k < inb_sequence;k++) {
              out_file << "\t";
              for (m = 0;m < length[k];m++) {
                if (index_parameter[k][m] == index_parameter[0][j]) {
                  out_file << real_sequence[k][i][m];
                  break;
                }
              }
            }
          }

          else {
            for (k = 1;k < inb_sequence;k++) {
              out_file << "\t";
              if (j < length[k]) {
                out_file << real_sequence[k][i][j];
              }
            }
          }
        }

        out_file << endl;
      }

      out_file << endl;
    }

    delete [] t_value;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la sequence des moyennes (et des ecart-types).
 *
 *  arguments : reference sur un objet Format_error, flag calcul des ecart-types,
 *              sortie (sequences, residus ou residus standardisees),
 *              path, format ('a' : ASCII / 's' : Spreadsheet).
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::pointwise_average(Format_error &error , bool standard_deviation ,
                                        int output , const char *path , char format) const

{
  bool status = true;
  register int i , j , k;
  int inb_sequence , min_identifier , max_identifier , *iidentifier , *ilength , *itype ,
      *pindex_param , *cindex_param , *cisequence , *size;
  double diff , *prsequence , *crsequence , *pmean , *pstandard_deviation;
  Sequences *seq;


  seq = 0;
  error.init();

  if (nb_sequence == 1) {
    status = false;
    error.update(SEQ_error[SEQR_SINGLE_SEQUENCE]);
  }
  if ((index_parameter_type != IMPLICIT_TYPE) && (index_parameter_type != TIME)) {
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
  }

  if (status) {
    if ((output == STANDARDIZED_RESIDUAL) && (!standard_deviation)) {
      standard_deviation = true;
    }

    inb_sequence = nb_sequence + (standard_deviation ? 2 : 1);

    iidentifier = new int[inb_sequence];

    min_identifier = identifier[0];
    for (i = 0;i < nb_sequence;i++) {
      if (identifier[i] < min_identifier) {
        min_identifier = identifier[i];
      }

      iidentifier[i + 1] = identifier[i];
    }

    iidentifier[0] = min_identifier - 1;

    if (standard_deviation) {
      max_identifier = identifier[nb_sequence - 1];
      for (i = 0;i < nb_sequence - 1;i++) {
        if (identifier[i] > max_identifier) {
          max_identifier = identifier[i];
        }
      }

      iidentifier[inb_sequence - 1] = max_identifier + 1;
    }

    ilength = new int[inb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      ilength[i + 1] = length[i];
    }

    if (index_parameter) {
      ilength[0] = 0;
      for (i = hindex_parameter->offset;i < hindex_parameter->nb_value;i++) {
        if (hindex_parameter->frequency[i] > 0) {
          ilength[0]++;
        }
      }
    }

    else {
      ilength[0] = max_length;
    }

    if (standard_deviation) {
      ilength[inb_sequence - 1] = ilength[0];
    }

    itype = new int[nb_variable];
    for (i = 0;i < nb_variable;i++) {
      itype[i] = REAL_VALUE;
    }

    seq = new Sequences(inb_sequence , iidentifier , ilength , index_parameter_type ,
                        nb_variable , itype);

    delete [] iidentifier;
    delete [] ilength;
    delete [] itype;

    if (index_parameter) {
      pindex_param = seq->index_parameter[0];
      for (i = hindex_parameter->offset;i < hindex_parameter->nb_value;i++) {
        if (hindex_parameter->frequency[i] > 0) {
          *pindex_param++ = i;
        }
      }

      if (standard_deviation) {
        pindex_param = seq->index_parameter[seq->nb_sequence - 1];
        cindex_param = seq->index_parameter[0];
        for (i = 0;i < seq->length[seq->nb_sequence - 1];i++) {
          *pindex_param++ = *cindex_param++;
        }
      }

      // copie des parametres d'index

      for (i = 0;i < nb_sequence;i++) {
        pindex_param = seq->index_parameter[i + 1];
        cindex_param = index_parameter[i];
        for (j = 0;j < length[i];j++) {
          *pindex_param++ = *cindex_param++;
        }
      }

      seq->build_index_parameter_histogram();
      seq->index_interval_computation();
    }

    // calcul des tailles d'echantillons

    size = new int[seq->length[0]];

    if (index_parameter) {
      pindex_param = seq->index_parameter[0];
      i = 0;
      for (j = hindex_parameter->offset;j < hindex_parameter->nb_value;j++) {
        if (hindex_parameter->frequency[j] > 0) {
          size[i++] = hindex_parameter->frequency[j];
        }
      }
    }

    else {
      size[0] = nb_sequence;
      for (i = 1;i < max_length;i++) {
        size[i] = size[i - 1] - hlength->frequency[i];
      }
    }

    // calcul des moyennes

    for (i = 0;i < seq->nb_variable;i++) {
      prsequence = seq->real_sequence[0][i];
      for (j = 0;j < seq->length[0];j++) {
        *prsequence++ = 0.;
      }
    }

    if (index_parameter) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < nb_variable;j++) {
          pindex_param = seq->index_parameter[0];
          cindex_param = index_parameter[i];
          prsequence = seq->real_sequence[0][j];
          if (type[j] != REAL_VALUE) {
            cisequence = int_sequence[i][j];
            for (k = 0;k < seq->length[0];k++) {
              if (*pindex_param++ == *cindex_param) {
                cindex_param++;
                *prsequence += *cisequence++;
              }
              prsequence++;
            }
          }

          else {
            crsequence = real_sequence[i][j];
            for (k = 0;k < seq->length[0];k++) {
              if (*pindex_param++ == *cindex_param) {
                cindex_param++;
                *prsequence += *crsequence++;
              }
              prsequence++;
            }
          }
        }
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < nb_variable;j++) {
          prsequence = seq->real_sequence[0][j];
          if (type[j] != REAL_VALUE) {
            cisequence = int_sequence[i][j];
            for (k = 0;k < length[i];k++) {
              *prsequence++ += *cisequence++;
            }
          }

          else {
            crsequence = real_sequence[i][j];
            for (k = 0;k < length[i];k++) {
              *prsequence++ += *crsequence++;
            }
          }
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      for (j = 0;j < seq->length[0];j++) {
        seq->real_sequence[0][i][j] /= size[j];
      }
    }

    // calcul des ecart-types

    if (standard_deviation) {
      for (i = 0;i < seq->nb_variable;i++) {
        prsequence = seq->real_sequence[seq->nb_sequence - 1][i];
        for (j = 0;j < seq->length[seq->nb_sequence - 1];j++) {
          *prsequence++ = 0.;
        }
      }

      if (index_parameter) {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < nb_variable;j++) {
            pindex_param = seq->index_parameter[seq->nb_sequence - 1];
            cindex_param = index_parameter[i];
            prsequence = seq->real_sequence[seq->nb_sequence - 1][j];
            pmean = seq->real_sequence[0][j];
            if (type[j] != REAL_VALUE) {
              cisequence = int_sequence[i][j];
              for (k = 0;k < seq->length[seq->nb_sequence - 1];k++) {
                if (*pindex_param++ == *cindex_param) {
                  cindex_param++;
                  diff = *cisequence++ - *pmean;
                  *prsequence += diff * diff;
                }
                prsequence++;
                pmean++;
              }
            }

            else {
              crsequence = real_sequence[i][j];
              for (k = 0;k < seq->length[seq->nb_sequence - 1];k++) {
                if (*pindex_param++ == *cindex_param) {
                  cindex_param++;
                  diff = *crsequence++ - *pmean;
                  *prsequence += diff * diff;
                }
                prsequence++;
                pmean++;
              }
            }
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < nb_variable;j++) {
            prsequence = seq->real_sequence[seq->nb_sequence - 1][j];
            pmean = seq->real_sequence[0][j];
            if (type[j] != REAL_VALUE) {
              cisequence = int_sequence[i][j];
              for (k = 0;k < length[i];k++) {
                diff = *cisequence++ - *pmean++;
                *prsequence++ += diff * diff;
              }
            }

            else {
              crsequence = real_sequence[i][j];
              for (k = 0;k < length[i];k++) {
                diff = *crsequence++ - *pmean++;
                *prsequence++ += diff * diff;
              }
            }
          }
        }
      }

      for (i = 0;i < seq->nb_variable;i++) {
        for (j = 0;j < seq->length[seq->nb_sequence - 1];j++) {
          if (size[j] > 1) {
            seq->real_sequence[seq->nb_sequence - 1][i][j] = sqrt(seq->real_sequence[seq->nb_sequence - 1][i][j] /
                                                                  (size[j] - 1));
          }
        }
      }
    }

    switch (output) {

    // copie des sequences

    case SEQUENCE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < nb_variable;j++) {
          prsequence = seq->real_sequence[i + 1][j];
          if (type[j] != REAL_VALUE) {
            cisequence = int_sequence[i][j];
            for (k = 0;k < length[i];k++) {
              *prsequence++ = *cisequence++;
            }
          }

          else {
            crsequence = real_sequence[i][j];
            for (k = 0;k < length[i];k++) {
              *prsequence++ = *crsequence++;
            }
          }
        }
      }
      break;
    }

    // calcul des residus

    case SUBTRACTION_RESIDUAL : {
      if (index_parameter) {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < nb_variable;j++) {
            pindex_param = seq->index_parameter[0];
            cindex_param = index_parameter[i];
            prsequence = seq->real_sequence[i + 1][j];
            pmean = seq->real_sequence[0][j];
            if (type[j] != REAL_VALUE) {
              cisequence = int_sequence[i][j];
              for (k = 0;k < seq->length[0];k++) {
                if (*pindex_param++ == *cindex_param) {
                  cindex_param++;
                  *prsequence++ = *cisequence++ - *pmean;
                }
                pmean++;
              }
            }

            else {
              crsequence = real_sequence[i][j];
              for (k = 0;k < seq->length[0];k++) {
                if (*pindex_param++ == *cindex_param) {
                  cindex_param++;
                  *prsequence++ = *crsequence++ - *pmean;
                }
                pmean++;
              }
            }
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < nb_variable;j++) {
            prsequence = seq->real_sequence[i + 1][j];
            pmean = seq->real_sequence[0][j];
            if (type[j] != REAL_VALUE) {
              cisequence = int_sequence[i][j];
              for (k = 0;k < length[i];k++) {
                *prsequence++ = *cisequence++ - *pmean++;
              }
            }

            else {
              crsequence = real_sequence[i][j];
              for (k = 0;k < length[i];k++) {
                *prsequence++ = *crsequence++ - *pmean++;
              }
            }
          }
        }
      }
      break;
    }

    // calcul des residus standardises

    case STANDARDIZED_RESIDUAL : {
      if (index_parameter) {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < nb_variable;j++) {
            pindex_param = seq->index_parameter[0];
            cindex_param = index_parameter[i];
            prsequence = seq->real_sequence[i + 1][j];
            pmean = seq->real_sequence[0][j];
            pstandard_deviation = seq->real_sequence[seq->nb_sequence - 1][j];
            if (type[j] != REAL_VALUE) {
              cisequence = int_sequence[i][j];
              for (k = 0;k < seq->length[0];k++) {
                if (*pindex_param++ == *cindex_param) {
                  cindex_param++;
                  if (*pstandard_deviation > 0.) {
                    *prsequence++ = (*cisequence++ - *pmean) / *pstandard_deviation;
                  }
                  else {
                    *prsequence++ = 0.;
                    cisequence++;
                  }
                }
                pmean++;
                pstandard_deviation++;
              }
            }

            else {
              crsequence = real_sequence[i][j];
              for (k = 0;k < seq->length[0];k++) {
                if (*pindex_param++ == *cindex_param) {
                  cindex_param++;
                  if (*pstandard_deviation > 0.) {
                    *prsequence++ = (*crsequence++ - *pmean) / *pstandard_deviation;
                  }
                  else {
                    *prsequence++ = 0.;
                    crsequence++;
                  }
                }
                pmean++;
                pstandard_deviation++;
              }
            }
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < nb_variable;j++) {
            prsequence = seq->real_sequence[i + 1][j];
            pmean = seq->real_sequence[0][j];
            pstandard_deviation = seq->real_sequence[seq->nb_sequence - 1][j];
            if (type[j] != REAL_VALUE) {
              cisequence = int_sequence[i][j];
              for (k = 0;k < length[i];k++) {
                if (*pstandard_deviation > 0.) {
                  *prsequence++ = (*cisequence++ - *pmean++) / *pstandard_deviation;
                }
                else {
                  *prsequence++ = 0.;
                  cisequence++;
                }
                pstandard_deviation++;
              }
            }

            else {
              crsequence = real_sequence[i][j];
              for (k = 0;k < length[i];k++) {
                if (*pstandard_deviation > 0.) {
                  *prsequence++ = (*crsequence++ - *pmean++) / *pstandard_deviation;
                }
                else {
                  *prsequence++ = 0.;
                  crsequence++;
                }
                pstandard_deviation++;
              }
            }
          }
        }
      }
      break;
    }
    }

    if ((output == SEQUENCE) && (!standard_deviation)) {
      for (i = 0;i < seq->nb_variable;i++) {
        seq->min_value[i] = min_value[i];
        seq->max_value[i] = max_value[i];
      }
    }
    else {
      for (i = 0;i < seq->nb_variable;i++) {
        seq->min_value_computation(i);
        seq->max_value_computation(i);
      }
    }

    // ecriture des sequences de moyennes et d'ecart-types

    if (path) {
      switch (format) {
      case 'a' :
        status = seq->pointwise_average_ascii_print(error , path , size ,
                                                    standard_deviation , output);
        break;
      case 's' :
        status = seq->pointwise_average_spreadsheet_print(error , path , size ,
                                                          standard_deviation , output);
        break;
      }

      if (!status) {

#       ifdef MESSAGE
        cout << error;
#       endif

      }
    }

    delete [] size;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des sequences des temps de retour pour une valeur prise
 *  par une variable entiere.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable, valeur.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::recurrence_time_sequences(Format_error &error , int variable , int value) const

{
  bool status = true;
  register int i , j;
  int inb_sequence , ilength , previous_index , *pisequence , *prsequence;
  Sequences *seq;


  seq = 0;
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
      if (!marginal[variable]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << STAT_error[STATR_MARGINAL_HISTOGRAM];
        error.update((error_message.str()).c_str());
      }

      else if ((value < marginal[variable]->offset) || (value >= marginal[variable]->nb_value) ||
               (marginal[variable]->frequency[value] == 0)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << STAT_label[STATL_VALUE] << " " << value << " "
                      << SEQ_error[SEQR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {
    seq = new Sequences(nb_sequence , 1);

    // calcul des sequences des temps de retour

    inb_sequence = 0;

    for (i = 0;i < nb_sequence;i++) {
      pisequence = int_sequence[i][variable];
      previous_index = 0;
      ilength = 0;
      for (j = 0;j < length[i];j++) {
        if (*pisequence == value) {
          if (ilength == 0) {
            seq->int_sequence[inb_sequence][0] = new int[length[i]];
            prsequence = seq->int_sequence[inb_sequence][0];
          }

          *prsequence++ = j - previous_index;
          previous_index = j;
          ilength++;
        }

        pisequence++;
      }

      if (ilength > 0) {
        seq->length[inb_sequence] = ilength;
        seq->identifier[inb_sequence++] = identifier[i];
      }
    }

    seq->nb_sequence = inb_sequence;

    seq->max_length_computation();
    seq->cumul_length_computation();
    seq->build_length_histogram();

    seq->min_value_computation(0);
    seq->max_value_computation(0);
    seq->build_marginal_histogram(0);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des sequences des temps de sejour par une variable entiere.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::sojourn_time_sequences(Format_error &error , int variable) const

{
  bool status = true;
  register int i , j;
  int ilength , run_length , *pisequence , *pstate , *pssequence , itype[2];
  Sequences *seq;


  seq = 0;
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

    else if (!marginal[variable]) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_MARGINAL_HISTOGRAM];
      error.update((error_message.str()).c_str());
    }
  }

  if (status) {
    itype[0] = type[variable];
    itype[1] = INT_VALUE;

    seq = new Sequences(nb_sequence , identifier , length , IMPLICIT_TYPE , 2 , itype);

    // calcul des sequences de temps de sejour

    for (i = 0;i < nb_sequence;i++) {
      pstate = seq->int_sequence[i][0];
      pssequence = seq->int_sequence[i][1];
      pisequence = int_sequence[i][variable];
      run_length = 1;
      ilength = 0;

      for (j = 0;j < length[i] - 1;j++) {
        if (*(pisequence + 1) != *pisequence) {
          *pstate++ = *pisequence;
          *pssequence++ = run_length;
          run_length = 0;
          ilength++;
        }

        run_length++;
        pisequence++;
      }

      *pstate = *pisequence;
      *pssequence = run_length;

      seq->length[i] = ilength + 1;
    }

    seq->max_length_computation();
    seq->cumul_length_computation();
    delete seq->hlength;
    seq->build_length_histogram();

    for (i = 0;i < 2;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);
      seq->build_marginal_histogram(i);
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Discretisation des positions.
 *
 *  arguments : reference sur un objet Format_error,
 *              pas de discretisation.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::transform_position(Format_error &error , int step) const

{
  bool status = true;
  register int i , j , k , m;
  int inter_position , nb_unit , *ilength , **pisequence;
  Sequences *seq;


  seq = 0;
  error.init();

  if (index_parameter_type != POSITION) {
    status = false;
    error.correction_update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE] , SEQ_index_parameter_word[POSITION]);
  }

  else if ((step < 1) || ((index_interval) && (step > index_interval->mean))) {
    status = false;
    error.update(SEQ_error[SEQR_POSITION_STEP]);
  }

  for (i = 0;i < nb_variable;i++) {
    if ((type[i] != INT_VALUE) && (type[i] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (status) {
    ilength = new int[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      ilength[i] = index_parameter[i][length[i]] / step + 1 + length[i];
    }

    seq = new Sequences(nb_sequence , identifier , ilength , nb_variable);
    delete [] ilength;

    // extraction des sequences

    pisequence = new int*[nb_variable];

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < nb_variable;j++) {
        pisequence[j] = seq->int_sequence[i][j];
      }
      seq->length[i] = 0;

      for (j = 0;j <= length[i];j++) {
        if (j == 0) {
          inter_position = index_parameter[i][j];
        }
        else {
          inter_position = index_parameter[i][j] - index_parameter[i][j - 1];
        }

        nb_unit = (inter_position % step == 0 ? inter_position / step : inter_position / step + 1);
        if ((nb_unit > 0) && (j < length[i])) {
          nb_unit--;
        }

        for (k = 0;k < nb_variable;k++) {
          for (m = 0;m < nb_unit;m++) {
            *pisequence[k]++ = (int)min_value[k] - 1;
          }
          if (j < length[i]) {
            *pisequence[k]++ = int_sequence[i][k][j];
          }
        }

        if (j < length[i]) {
          seq->length[i] += nb_unit + 1;
        }
        else {
          seq->length[i] += nb_unit;
        }
      }
    }

    seq->max_length_computation();
    seq->cumul_length_computation();
    delete seq->hlength;
    seq->build_length_histogram();

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value[i] = min_value[i] - 1;
      seq->max_value[i] = max_value[i];
      seq->build_marginal_histogram(i);
    }

    delete [] pisequence;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Croisement des sequences.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::cross(Format_error &error) const

{
  bool status = true;
  register int i , j , k , m;
  int sense = 0 , *ilength;
  Sequences *seq;


  seq = 0;
  error.init();

  if (type[0] != INT_VALUE) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_VARIABLE] << " " << 1 << ": "
                  << STAT_error[STATR_VARIABLE_TYPE];
    error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
  }

  for (i = 1;i < nb_sequence;i++) {
    if (length[i] > length[i - 1]) {
      if (sense == 0) {
        sense++;
      }
      else if (sense == -1) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEQUENCE] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_LENGTH];
        error.update((error_message.str()).c_str());
      }
    }

    else if (length[i] < length[i - 1]) {
      if (sense == 0) {
        sense--;
      }
      else if (sense == 1) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEQUENCE] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_LENGTH];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {
    ilength = new int[max_length];
    for (i = 0;i < max_length;i++) {
      ilength[i] = nb_sequence;
    }

    seq = new Sequences(max_length , 0 , ilength , index_parameter_type ,
                        nb_variable , type);
    delete [] ilength;

    // constitution des sequences croisees

    for (i = 0;i < seq->nb_sequence;i++) {
      j = 0;
      while (length[j] <= i) {
        j++;
      }

      k = 0;
      do {
        for (m = 0;m < seq->nb_variable;m++) {
          if (seq->type[m] != REAL_VALUE) {
            seq->int_sequence[i][m][k] = int_sequence[j][m][i];
          }
          else {
            seq->real_sequence[i][m][k] = real_sequence[j][m][i];
          }
        }
        j++;
        k++;
      }
      while ((j < nb_sequence) && (length[j] > i));

      seq->length[i] = k;
    }

    seq->max_length = nb_sequence;
    seq->cumul_length = cumul_length;
    delete seq->hlength;
    seq->build_length_histogram();

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value[i] = min_value[i];
      seq->max_value[i] = max_value[i];
      if (marginal[i]) {
        seq->marginal[i] = new Histogram(*marginal[i]);
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la longueur maximum des sequences.
 *
 *--------------------------------------------------------------*/

void Sequences::max_length_computation()

{
  register int i;


  max_length = length[0];
  for (i = 1;i < nb_sequence;i++) {
    if (length[i] > max_length) {
      max_length = length[i];
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la longueur cumulee des sequences.
 *
 *--------------------------------------------------------------*/

void Sequences::cumul_length_computation()

{
  register int i;


  cumul_length = 0;
  for (i = 0;i < nb_sequence;i++) {
    cumul_length += length[i];
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction de l'histogramme des longueurs des sequences.
 *
 *--------------------------------------------------------------*/

void Sequences::build_length_histogram()

{
  register int i;


  hlength = new Histogram(max_length + 1);

  hlength->nb_element = nb_sequence;
  for (i = 0;i < nb_sequence;i++) {
    (hlength->frequency[length[i]])++;
  }

  hlength->nb_value_computation();
  hlength->offset_computation();
  hlength->max_computation();
  hlength->mean_computation();
  hlength->variance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Calcul des temps a partir des intervalles de temps /
 *  calcul des positions a partir des intervalles inter-positions.
 *
 *--------------------------------------------------------------*/

void Sequences::index_parameter_computation()

{
  if ((index_parameter_type == TIME_INTERVAL) || (index_parameter_type == POSITION_INTERVAL)) {
    register int i , j;


    switch (index_parameter_type) {
    case TIME_INTERVAL :
      index_parameter_type = TIME;
      break;
    case POSITION_INTERVAL :
      index_parameter_type = POSITION;
      break;
    }

    for (i = 0;i < nb_sequence;i++) {
      for (j = 1;j < (index_parameter_type == POSITION ? length[i] + 1 : length[i]);j++) {
        index_parameter[i][j] += index_parameter[i][j - 1];
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur minimum prise par le parametre d'index.
 *
 *--------------------------------------------------------------*/

int Sequences::min_index_parameter_computation() const

{
  register int i;
  int min_index_parameter = I_DEFAULT;


  if (index_parameter) {
    min_index_parameter = index_parameter[0][0];
    for (i = 1;i < nb_sequence;i++) {
      if (index_parameter[i][0] < min_index_parameter) {
        min_index_parameter = index_parameter[i][0];
      }
    }
  }

  return min_index_parameter;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur maximum prise par le parametre d'index.
 *
 *  argument : flag derniere position.
 *
 *--------------------------------------------------------------*/

int Sequences::max_index_parameter_computation(bool last_position) const

{
  register int i;
  int max_index_parameter = I_DEFAULT;


  if (index_parameter) {
    if ((index_parameter_type == TIME) || (last_position)) {
      max_index_parameter = index_parameter[0][length[0] - 1];
      for (i = 1;i < nb_sequence;i++) {
        if (index_parameter[i][length[i] - 1] > max_index_parameter) {
          max_index_parameter = index_parameter[i][length[i] - 1];
        }
      }
    }

    else {
      max_index_parameter = index_parameter[0][length[0]];
      for (i = 1;i < nb_sequence;i++) {
        if (index_parameter[i][length[i]] > max_index_parameter) {
          max_index_parameter = index_parameter[i][length[i]];
        }
      }
    }
  }

  return max_index_parameter;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'histogramme des parametres d'index.
 *
 *--------------------------------------------------------------*/

void Sequences::build_index_parameter_histogram()

{
  if (index_parameter) {
    register int i , j;


    hindex_parameter = new Histogram(max_index_parameter_computation() + 1);

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < (index_parameter_type == POSITION ? length[i] + 1 : length[i]);j++) {
        (hindex_parameter->frequency[index_parameter[i][j]])++;
      }
    }

    hindex_parameter->offset_computation();
    hindex_parameter->nb_element = cumul_length;
    if (index_parameter_type == POSITION) {
      hindex_parameter->nb_element += nb_sequence;
    }
    hindex_parameter->max_computation();
    hindex_parameter->mean_computation();
    hindex_parameter->variance_computation();
  }
}


/*--------------------------------------------------------------*
 *
 *  Extraction de l'histogramme des intervalles entre index successifs.
 *
 *--------------------------------------------------------------*/

void Sequences::index_interval_computation()

{
  if ((index_parameter_type == TIME) || ((index_parameter_type == POSITION) &&
       (type[0] != NB_INTERNODE))) {
    register int i , j;


    index_interval = new Histogram(max_index_parameter_computation(true) + 1);

    // constitution de l'histogramme des intervalles entre index successifs

    for (i = 0;i < nb_sequence;i++) {
      for (j = 1;j < length[i];j++) {
        (index_interval->frequency[index_parameter[i][j] - index_parameter[i][j - 1]])++;
      }
    }

    index_interval->nb_value_computation();
    index_interval->offset_computation();
    index_interval->nb_element = cumul_length - nb_sequence;
    index_interval->max_computation();
    index_interval->mean_computation();
    index_interval->variance_computation();
  }
}


/*--------------------------------------------------------------*
 *
 *  Extraction de l'histogramme des intervalles entre index successifs
 *  pour une valeur d'une variable entiere.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable, valeur.
 *
 *--------------------------------------------------------------*/

Histogram* Sequences::value_index_interval_computation(Format_error &error , int variable ,
                                                       int value) const

{
  bool status = true;
  register int i , j;
  int previous_index_param , *pindex_param , *pisequence;
  Histogram *value_index_interval;


  value_index_interval = 0;
  error.init();

  if ((index_parameter_type != TIME) && (index_parameter_type != POSITION)) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (!marginal[variable]) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_MARGINAL_HISTOGRAM];
      error.update((error_message.str()).c_str());
    }

    else if ((value < marginal[variable]->offset) || (value >= marginal[variable]->nb_value) ||
             (marginal[variable]->frequency[value] <= 1)) {
      status = false;
      error.update(SEQ_error[SEQR_VALUE]);
    }
  }

  if (status) {
    value_index_interval = new Histogram(max_index_parameter_computation(true) + 1);

    for (i = 0;i < nb_sequence;i++) {
      pindex_param = index_parameter[i];
      pisequence = int_sequence[i][variable];
      previous_index_param = I_DEFAULT;

      for (j = 0;j < length[i];j++) {
        if (*pisequence == value) {
          if (previous_index_param != I_DEFAULT) {
            (value_index_interval->frequency[*pindex_param - previous_index_param])++;
          }
          previous_index_param = *pindex_param;
        }

        pindex_param++;
        pisequence++;
      }
    }

    // extraction des caracteristiques de l'histogramme

    value_index_interval->nb_value_computation();
    value_index_interval->offset_computation();
    value_index_interval->nb_element_computation();

    if (value_index_interval->nb_element == 0) {
      delete value_index_interval;
      value_index_interval = 0;
      error.update(STAT_error[STATR_EMPTY_HISTOGRAM]);
    }

    else {
      value_index_interval->max_computation();
      value_index_interval->mean_computation();
      value_index_interval->variance_computation();
    }
  }

  return value_index_interval;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur minimum prise par une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Sequences::min_value_computation(int variable)

{
  register int i , j;
  int *pisequence;
  double *prsequence;


  if ((type[variable] != REAL_VALUE) && (type[variable] != AUXILIARY)) {
    min_value[variable] = int_sequence[0][variable][0];

    for (i = 0;i < nb_sequence;i++) {
      pisequence = int_sequence[i][variable];
      for (j = 0;j < length[i];j++) {
        if (*pisequence < min_value[variable]) {
          min_value[variable] = *pisequence;
        }
        pisequence++;
      }
    }
  }

  else {
    min_value[variable] = real_sequence[0][variable][0];

    for (i = 0;i < nb_sequence;i++) {
      prsequence = real_sequence[i][variable];
      for (j = 0;j < length[i];j++) {
        if (*prsequence < min_value[variable]) {
          min_value[variable] = *prsequence;
        }
        prsequence++;
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

void Sequences::max_value_computation(int variable)

{
  register int i , j;
  int *pisequence;
  double *prsequence;


  if ((type[variable] != REAL_VALUE) && (type[variable] != AUXILIARY)) {
    max_value[variable] = int_sequence[0][variable][0];

    for (i = 0;i < nb_sequence;i++) {
      pisequence = int_sequence[i][variable];
      for (j = 0;j < length[i];j++) {
        if (*pisequence > max_value[variable]) {
          max_value[variable] = *pisequence;
        }
        pisequence++;
      }
    }
  }

  else {
    max_value[variable] = real_sequence[0][variable][0];

    for (i = 0;i < nb_sequence;i++) {
      prsequence = real_sequence[i][variable];
      for (j = 0;j < length[i];j++) {
        if (*prsequence > max_value[variable]) {
          max_value[variable] = *prsequence;
        }
        prsequence++;
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la loi marginale empirique pour une variable entiere positive.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Sequences::marginal_histogram_computation(int variable)

{
  register int i , j;
  int *pfrequency , *pisequence;


  pfrequency = marginal[variable]->frequency;
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    *pfrequency++ = 0;
  }

  for (i = 0;i < nb_sequence;i++) {
    pisequence = int_sequence[i][variable];
    for (j = 0;j < length[i];j++) {
      (marginal[variable]->frequency[*pisequence++])++;
    }
  }

  marginal[variable]->offset = (int)min_value[variable];
  marginal[variable]->nb_element_computation();
  marginal[variable]->max_computation();
  marginal[variable]->mean_computation();
  marginal[variable]->variance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Construction de la loi marginale empirique pour une variable entiere positive.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Sequences::build_marginal_histogram(int variable)

{
  if ((type[variable] != REAL_VALUE) && (type[variable] != AUXILIARY) &&
      (min_value[variable] >= 0) && (max_value[variable] <= MARGINAL_MAX_VALUE)) {
    marginal[variable] = new Histogram((int)max_value[variable] + 1);
    marginal_histogram_computation(variable);
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la moyenne d'une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

double Sequences::mean_computation(int variable) const

{
  register int i , j;
  int *pisequence;
  double mean , *prsequence;


  if (marginal[variable]) {
    mean = marginal[variable]->mean;
  }

  else {
    mean = 0.;

    if (type[variable] != REAL_VALUE) {
      for (i = 0;i < nb_sequence;i++) {
        pisequence = int_sequence[i][variable];
        for (j = 0;j < length[i];j++) {
          mean += *pisequence++;
        }
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        prsequence = real_sequence[i][variable];
        for (j = 0;j < length[i];j++) {
          mean += *prsequence++;
        }
      }
    }

    mean /= cumul_length;
  }

  return mean;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la variance d'une variable.
 *
 *  arguments : indice de la variable, moyenne.
 *
 *--------------------------------------------------------------*/

double Sequences::variance_computation(int variable , double mean) const

{
  register int i , j;
  int *pisequence;
  double variance , diff , *prsequence;


  if (marginal[variable]) {
    variance = marginal[variable]->variance;
  }

  else {
    variance = 0.;

    if (cumul_length > 1) {
      if (type[variable] != REAL_VALUE) {
        for (i = 0;i < nb_sequence;i++) {
          pisequence = int_sequence[i][variable];
          for (j = 0;j < length[i];j++) {
            diff = *pisequence++ - mean;
            variance += diff * diff;
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          prsequence = real_sequence[i][variable];
          for (j = 0;j < length[i];j++) {
            diff = *prsequence++ - mean;
            variance += diff * diff;
          }
        }
      }

      variance /= (cumul_length - 1);
    }
  }

  return variance;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'ecart absolu moyen d'une variable.
 *
 *  arguments : indice de la variable, moyenne.
 *
 *--------------------------------------------------------------*/

double Sequences::mean_absolute_deviation_computation(int variable , double mean) const

{
  register int i , j;
  int *pisequence;
  double mean_absolute_deviation , *prsequence;


  if (marginal[variable]) {
    mean_absolute_deviation = marginal[variable]->mean_absolute_deviation_computation();
  }

  else {
    mean_absolute_deviation = 0.;

    if (type[variable] != REAL_VALUE) {
      for (i = 0;i < nb_sequence;i++) {
        pisequence = int_sequence[i][variable];
        for (j = 0;j < length[i];j++) {
          mean_absolute_deviation += fabs(*pisequence - mean);
          pisequence++;
        }
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        prsequence = real_sequence[i][variable];
        for (j = 0;j < length[i];j++) {
          mean_absolute_deviation += fabs(*prsequence - mean);
          prsequence++;
        }
      }
    }

    mean_absolute_deviation /= cumul_length;
  }

  return mean_absolute_deviation;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la difference absolue moyenne d'une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

double Sequences::mean_absolute_difference_computation(int variable) const

{
  register int i , j , k , m;
  int *pisequence1 , *pisequence2;
  double mean_absolute_difference , *prsequence1 , *prsequence2;


  mean_absolute_difference = 0.;

  if (cumul_length > 1) {
    if (type[variable] != REAL_VALUE) {
      for (i = 0;i < nb_sequence;i++) {
        pisequence1 = int_sequence[i][variable];
        for (j = 0;j < length[i];j++) {
          pisequence2 = int_sequence[i][variable] + j + 1;
          for (k = j + 1;k < length[i];k++) {
            mean_absolute_difference += abs(*pisequence1 - *pisequence2);
            pisequence2++;
          }
          pisequence1++;
        }

        for (j = i + 1;j < nb_sequence;j++) {
          pisequence1 = int_sequence[i][variable];
          for (k = 0;k < length[i];k++) {
            pisequence2 = int_sequence[j][variable];
            for (m = 0;m < length[j];m++) {
              mean_absolute_difference += abs(*pisequence1 - *pisequence2);
              pisequence2++;
            }
            pisequence1++;
          }
        }
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        prsequence1 = real_sequence[i][variable];
        for (j = 0;j < length[i];j++) {
          prsequence2 = real_sequence[i][variable] + j + 1;
          for (k = j + 1;k < length[i];k++) {
            mean_absolute_difference += fabs(*prsequence1 - *prsequence2);
            prsequence2++;
          }
          prsequence1++;
        }

        for (j = i + 1;j < nb_sequence;j++) {
          prsequence1 = real_sequence[i][variable];
          for (k = 0;k < length[i];k++) {
            prsequence2 = real_sequence[j][variable];
            for (m = 0;m < length[j];m++) {
              mean_absolute_difference += fabs(*prsequence1 - *prsequence2);
              prsequence2++;
            }
            prsequence1++;
          }
        }
      }
    }

    mean_absolute_difference = 2 * mean_absolute_difference /
                               (cumul_length * (double)(cumul_length - 1));
  }

  return mean_absolute_difference;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du coefficient d'asymetrie d'une variable.
 *
 *  arguments : indice de la variable, moyenne, variance.
 *
 *--------------------------------------------------------------*/

double Sequences::skewness_computation(int variable , double mean , double variance) const

{
  register int i , j;
  int *pisequence;
  double skewness , diff , *prsequence;


  if (marginal[variable]) {
    skewness = marginal[variable]->skewness_computation();
  }

  else {
    skewness = 0.;

    if ((cumul_length > 2) && (variance > 0.)) {
      if (type[variable] != REAL_VALUE) {
        for (i = 0;i < nb_sequence;i++) {
          pisequence = int_sequence[i][variable];
          for (j = 0;j < length[i];j++) {
            diff = *pisequence++ - mean;
            skewness += diff * diff * diff;
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          prsequence = real_sequence[i][variable];
          for (j = 0;j < length[i];j++) {
            diff = *prsequence++ - mean;
            skewness += diff * diff * diff;
          }
        }
      }

      skewness = skewness * cumul_length / ((cumul_length - 1) * (double)(cumul_length - 2) *
                  pow(variance , 1.5));
    }
  }

  return skewness;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'exces d'applatissement d'une variable :
 *  exces d'applatissement = coefficient d'applatissement - 3.
 *
 *  arguments : indice de la variable, moyenne, variance.
 *
 *--------------------------------------------------------------*/

double Sequences::kurtosis_computation(int variable , double mean , double variance) const

{
  register int i , j;
  int *pisequence;
  double kurtosis , diff , *prsequence;


  if (marginal[variable]) {
    kurtosis = marginal[variable]->kurtosis_computation();
  }

  else {
    if (variance == 0.) {
      kurtosis = -2.;
    }

    else {
      kurtosis = 0.;

      if (type[variable] != REAL_VALUE) {
        for (i = 0;i < nb_sequence;i++) {
          pisequence = int_sequence[i][variable];
          for (j = 0;j < length[i];j++) {
            diff = *pisequence++ - mean;
            kurtosis += diff * diff * diff * diff;
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          prsequence = real_sequence[i][variable];
          for (j = 0;j < length[i];j++) {
            diff = *prsequence++ - mean;
            kurtosis += diff * diff * diff * diff;
          }
        }
      }

      kurtosis = kurtosis / ((cumul_length - 1) * variance * variance) - 3.;
    }
  }

  return kurtosis;
}
