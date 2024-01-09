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
 *       $Id: distance_matrix.cpp 17995 2015-04-23 06:54:49Z guedon $
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
#include <iomanip>
#include <cstring>
#include <string.h>

#include "tool/config.h"

#include "stat_tools.h"
#include "distance_matrix.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe DistanceMatrix.
 *
 *--------------------------------------------------------------*/

DistanceMatrix::DistanceMatrix()

{
  nb_row = 0;
  nb_column = 0;

  row_identifier = NULL;
  column_identifier = NULL;

  distance = NULL;
  length = NULL;

  deletion_distance = NULL;
  nb_deletion = NULL;
  insertion_distance = NULL;
  nb_insertion = NULL;
  nb_match = NULL;
  substitution_distance = NULL;
  nb_substitution = NULL;
  transposition_distance = NULL;
  nb_transposition = NULL;

  label_size = 0;
  label = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe DistanceMatrix.
 *
 *  arguments : nombre de formes, label, identificateurs des formes,
 *
 *--------------------------------------------------------------*/

DistanceMatrix::DistanceMatrix(int nb_pattern , const char *ilabel , int *pattern_identifier)

{
  register int i , j;


  nb_row = nb_pattern;
  nb_column = nb_pattern;

  row_identifier = new int[nb_row];
  if (pattern_identifier) {
    for (i = 0;i < nb_row;i++) {
      row_identifier[i] = pattern_identifier[i];
    }
  }
  else {
    for (i = 0;i < nb_row;i++) {
      row_identifier[i] = i + 1;
    }
  }

  column_identifier = new int[nb_column];
  if (pattern_identifier) {
    for (i = 0;i < nb_column;i++) {
      column_identifier[i] = pattern_identifier[i];
    }
  }
  else {
    for (i = 0;i < nb_column;i++) {
      column_identifier[i] = i + 1;
    }
  }

  distance = new double*[nb_row];
  length = new int*[nb_row];
  for (i = 0;i < nb_row;i++) {
    distance[i] = new double[nb_column];
    length[i] = new int[nb_column];

    for (j = 0;j < nb_column;j++) {
      if (row_identifier[i] != column_identifier[j]) {
        distance[i][j] = -D_INF;
      }
      else {
        distance[i][j] = 0.;
      }

      length[i][j] = 0;
    }
  }

  deletion_distance = NULL;
  nb_deletion = NULL;
  insertion_distance = NULL;
  nb_insertion = NULL;
  nb_match = NULL;
  substitution_distance = NULL;
  nb_substitution = NULL;
  transposition_distance = NULL;
  nb_transposition = NULL;

  label_size = strlen(ilabel) + 1;
  label = new char[label_size];
  strcpy(label , ilabel);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe DistanceMatrix.
 *
 *  arguments : nombre de formes, identificateurs de ligne et de colonne, label,
 *              identificateurs des formes, flags operations de substitution et
 *              de transposition.
 *
 *--------------------------------------------------------------*/

DistanceMatrix::DistanceMatrix(int nb_pattern , int irow_identifier , int icolumn_identifier ,
                               const char *ilabel , int *pattern_identifier ,
                               bool substitution_flag , bool transposition_flag)

{
  register int i , j;


  nb_row = (irow_identifier == I_DEFAULT ? nb_pattern : 1);
  nb_column = (icolumn_identifier == I_DEFAULT ? nb_pattern : 1);

  row_identifier = new int[nb_row];

  if (irow_identifier == I_DEFAULT) {
    if (pattern_identifier) {
      for (i = 0;i < nb_row;i++) {
        row_identifier[i] = pattern_identifier[i];
      }
    }
    else {
      for (i = 0;i < nb_row;i++) {
        row_identifier[i] = i + 1;
      }
    }
  }
  else {
    row_identifier[0] = irow_identifier;
  }

  column_identifier = new int[nb_column];

  if (icolumn_identifier == I_DEFAULT) {
    if (pattern_identifier) {
      for (i = 0;i < nb_column;i++) {
        column_identifier[i] = pattern_identifier[i];
      }
    }
    else {
      for (i = 0;i < nb_column;i++) {
        column_identifier[i] = i + 1;
      }
    }
  }
  else {
    column_identifier[0] = icolumn_identifier;
  }

  distance = new double*[nb_row];
  length = new int*[nb_row];
  for (i = 0;i < nb_row;i++) {
    distance[i] = new double[nb_column];
    length[i] = new int[nb_column];

    for (j = 0;j < nb_column;j++) {
      if (row_identifier[i] != column_identifier[j]) {
        distance[i][j] = -D_INF;
      }
      else {
        distance[i][j] = 0.;
      }

      length[i][j] = 0;
    }
  }

  deletion_distance = new double*[nb_row];
  nb_deletion = new int*[nb_row];
  for (i = 0;i < nb_row;i++) {
    deletion_distance[i] = new double[nb_column];
    nb_deletion[i] = new int[nb_column];

    for (j = 0;j < nb_column;j++) {
      deletion_distance[i][j] = 0.;
      nb_deletion[i][j] = 0;
    }
  }

  insertion_distance = new double*[nb_row];
  nb_insertion = new int*[nb_row];
  for (i = 0;i < nb_row;i++) {
    insertion_distance[i] = new double[nb_column];
    nb_insertion[i] = new int[nb_column];

    for (j = 0;j < nb_column;j++) {
      insertion_distance[i][j] = 0.;
      nb_insertion[i][j] = 0;
    }
  }

  nb_match = new int*[nb_row];
  for (i = 0;i < nb_row;i++) {
    nb_match[i] = new int[nb_column];

    for (j = 0;j < nb_column;j++) {
      nb_match[i][j] = 0;
    }
  }

  if (substitution_flag) {
    substitution_distance = new double*[nb_row];
    nb_substitution = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      substitution_distance[i] = new double[nb_column];
      nb_substitution[i] = new int[nb_column];

      for (j = 0;j < nb_column;j++) {
        substitution_distance[i][j] = 0.;
        nb_substitution[i][j] = 0;
      }
    }
  }
  else {
    substitution_distance = NULL;
    nb_substitution = NULL;
  }

  if (transposition_flag) {
    transposition_distance = new double*[nb_row];
    nb_transposition = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      transposition_distance[i] = new double[nb_column];
      nb_transposition[i] = new int[nb_column];

      for (j = 0;j < nb_column;j++) {
        transposition_distance[i][j] = 0.;
        nb_transposition[i][j] = 0;
      }
    }
  }
  else {
    transposition_distance = NULL;
    nb_transposition = NULL;
  }

  label_size = strlen(ilabel) + 1;
  label = new char[label_size];
  strcpy(label , ilabel);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe DistanceMatrix.
 *
 *  arguments : reference sur un objet DistanceMatrix,
 *              nombre et identificateurs des formes selectionnees,
 *              flag pour conserver ou rejeter les formes selectionnees.
 *
 *--------------------------------------------------------------*/

DistanceMatrix::DistanceMatrix(const DistanceMatrix &dist_matrix , int inb_pattern ,
                               int *iidentifier , bool keep)

{
  register int i , j , k;
  int nb_pattern , dnb_pattern , *didentifier , *index , *rindex , *cindex;


  // calcul des identificateurs des lignes et des colonnes correspondant aux formes selectionnees

  dnb_pattern = MAX(dist_matrix.nb_row , dist_matrix.nb_column);
  if (dist_matrix.nb_row > 1) {
    didentifier = dist_matrix.row_identifier;
  }
  else {
    didentifier = dist_matrix.column_identifier;
  }

  switch (keep) {
  case false :
    nb_pattern = dnb_pattern - inb_pattern;
    break;
  case true :
    nb_pattern = inb_pattern;
    break;
  }

  index = new int[nb_pattern];

  switch (keep) {

  case false : {
    i = 0;
    for (j = 0;j < dnb_pattern;j++) {
      for (k = 0;k < inb_pattern;k++) {
        if (iidentifier[k] == didentifier[j]) {
          break;
        }
      }

      if (k == inb_pattern) {
        index[i++] = j;
      }
    }
    break;
  }

  case true : {
    for (i = 0;i < inb_pattern;i++) {
      for (j = 0;j < dnb_pattern;j++) {
        if (iidentifier[i] == didentifier[j]) {
          index[i] = j;
          break;
        }
      }
    }
    break;
  }
  }

  if (dist_matrix.nb_row == 1) {
    rindex = new int[1];
    rindex[0] = 1;
  }
  else {
    rindex = index;
  }

  if (dist_matrix.nb_column == 1) {
    cindex = new int[1];
    cindex[0] = 1;
  }
  else {
    cindex = index;
  }

  // copie des informations correspondants aux formes selectionnees

  nb_row = (dist_matrix.nb_row == 1 ? 1 : nb_pattern);
  nb_column = (dist_matrix.nb_column == 1 ? 1 : nb_pattern);

  row_identifier = new int[nb_row];
  for (i = 0;i < nb_row;i++) {
    row_identifier[i] = dist_matrix.row_identifier[rindex[i]];
  }

  column_identifier = new int[nb_column];
  for (i = 0;i < nb_column;i++) {
    column_identifier[i] = dist_matrix.column_identifier[cindex[i]];
  }

  distance = new double*[nb_row];
  length = new int*[nb_row];
  for (i = 0;i < nb_row;i++) {
    distance[i] = new double[nb_column];
    length[i] = new int[nb_column];
    for (j = 0;j < nb_column;j++) {
      distance[i][j] = dist_matrix.distance[rindex[i]][cindex[j]];
      length[i][j] = dist_matrix.length[rindex[i]][cindex[j]];
    }
  }

  if ((dist_matrix.deletion_distance) && (dist_matrix.nb_deletion)) {
    deletion_distance = new double*[nb_row];
    nb_deletion = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      deletion_distance[i] = new double[nb_column];
      nb_deletion[i] = new int[nb_column];
      for (j = 0;j < nb_column;j++) {
        deletion_distance[i][j] = dist_matrix.deletion_distance[rindex[i]][cindex[j]];
        nb_deletion[i][j] = dist_matrix.nb_deletion[rindex[i]][cindex[j]];
      }
    }
  }
  else {
    deletion_distance = NULL;
    nb_deletion = NULL;
  }

  if ((dist_matrix.insertion_distance) && (dist_matrix.nb_insertion)) {
    insertion_distance = new double*[nb_row];
    nb_insertion = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      insertion_distance[i] = new double[nb_column];
      nb_insertion[i] = new int[nb_column];
      for (j = 0;j < nb_column;j++) {
        insertion_distance[i][j] = dist_matrix.insertion_distance[rindex[i]][cindex[j]];
        nb_insertion[i][j] = dist_matrix.nb_insertion[rindex[i]][cindex[j]];
      }
    }
  }
  else {
    insertion_distance = NULL;
    nb_insertion = NULL;
  }

  if (dist_matrix.nb_match) {
    nb_match = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      nb_match[i] = new int[nb_column];
      for (j = 0;j < nb_column;j++) {
        nb_match[i][j] = dist_matrix.nb_match[rindex[i]][cindex[j]];
      }
    }
  }
  else {
    nb_match = NULL;
  }

  if ((dist_matrix.substitution_distance) && (dist_matrix.nb_substitution)) {
    substitution_distance = new double*[nb_row];
    nb_substitution = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      substitution_distance[i] = new double[nb_column];
      nb_substitution[i] = new int[nb_column];
      for (j = 0;j < nb_column;j++) {
        substitution_distance[i][j] = dist_matrix.substitution_distance[rindex[i]][cindex[j]];
        nb_substitution[i][j] = dist_matrix.nb_substitution[rindex[i]][cindex[j]];
      }
    }
  }
  else {
    substitution_distance = NULL;
    nb_substitution = NULL;
  }

  if ((dist_matrix.transposition_distance) && (dist_matrix.nb_transposition)) {
    transposition_distance = new double*[nb_row];
    nb_transposition = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      transposition_distance[i] = new double[nb_column];
      nb_transposition[i] = new int[nb_column];
      for (j = 0;j < nb_column;j++) {
        transposition_distance[i][j] = dist_matrix.transposition_distance[rindex[i]][cindex[j]];
        nb_transposition[i][j] = dist_matrix.nb_transposition[rindex[i]][cindex[j]];
      }
    }
  }
  else {
    transposition_distance = NULL;
    nb_transposition = NULL;
  }

  label_size = dist_matrix.label_size;
  label = new char[label_size];
  strcpy(label , dist_matrix.label);

  delete [] index;
  if (dist_matrix.nb_row == 1) {
    delete [] rindex;
  }
  if (dist_matrix.nb_column == 1) {
    delete [] cindex;
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe DistanceMatrix.
 *
 *  arguments : reference sur un objet DistanceMatrix, nombre de groupes, label.
 *
 *--------------------------------------------------------------*/

DistanceMatrix::DistanceMatrix(const DistanceMatrix &dist_matrix , int nb_cluster ,
                               const char *ilabel)

{
  register int i , j;


  nb_row = nb_cluster;
  nb_column = nb_cluster;

  row_identifier = new int[nb_row];
  for (i = 0;i < nb_row;i++) {
    row_identifier[i] = i + 1;
  }

  column_identifier = new int[nb_column];
  for (i = 0;i < nb_column;i++) {
    column_identifier[i] = i + 1;
  }

  distance = new double*[nb_row];
  length = new int*[nb_row];
  for (i = 0;i < nb_row;i++) {
    distance[i] = new double[nb_column];
    length[i] = new int[nb_column];
    for (j = 0;j < nb_column;j++) {
      distance[i][j] = 0.;
      length[i][j] = 0;
    }
  }

  if ((dist_matrix.deletion_distance) && (dist_matrix.nb_deletion)) {
    deletion_distance = new double*[nb_row];
    nb_deletion = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      deletion_distance[i] = new double[nb_column];
      nb_deletion[i] = new int[nb_column];
      for (j = 0;j < nb_column;j++) {
        deletion_distance[i][j] = 0.;
        nb_deletion[i][j] = 0;
      }
    }
  }

  else {
    deletion_distance = NULL;
    nb_deletion = NULL;
  }

  if ((dist_matrix.insertion_distance) && (dist_matrix.nb_insertion)) {
    insertion_distance = new double*[nb_row];
    nb_insertion = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      insertion_distance[i] = new double[nb_column];
      nb_insertion[i] = new int[nb_column];
      for (j = 0;j < nb_column;j++) {
        insertion_distance[i][j] = 0.;
        nb_insertion[i][j] = 0;
      }
    }
  }

  else {
    insertion_distance = NULL;
    nb_insertion = NULL;
  }

  if (dist_matrix.nb_match) {
    nb_match = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      nb_match[i] = new int[nb_column];
      for (j = 0;j < nb_column;j++) {
        nb_match[i][j] = 0;
      }
    }
  }

  else {
    nb_match = NULL;
  }

  if ((dist_matrix.substitution_distance) && (dist_matrix.nb_substitution)) {
    substitution_distance = new double*[nb_row];
    nb_substitution = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      substitution_distance[i] = new double[nb_column];
      nb_substitution[i] = new int[nb_column];
      for (j = 0;j < nb_column;j++) {
        substitution_distance[i][j] = 0.;
        nb_substitution[i][j] = 0;
      }
    }
  }

  else {
    substitution_distance = NULL;
    nb_substitution = NULL;
  }

  if ((dist_matrix.transposition_distance) && (dist_matrix.nb_transposition)) {
    transposition_distance = new double*[nb_row];
    nb_transposition = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      transposition_distance[i] = new double[nb_column];
      nb_transposition[i] = new int[nb_column];
      for (j = 0;j < nb_column;j++) {
        transposition_distance[i][j] = 0.;
        nb_transposition[i][j] = 0;
      }
    }
  }

  else {
    transposition_distance = NULL;
    nb_transposition = NULL;
  }

  label_size = strlen(ilabel) + 1;
  label = new char[label_size];
  strcpy(label , ilabel);
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet DistanceMatrix.
 *
 *  arguments : reference sur un objet DistanceMatrix,
 *              transformation de la matrice ('s' :  symetrisation, 'u' : denormalisation).
 *
 *--------------------------------------------------------------*/

void DistanceMatrix::copy(const DistanceMatrix &dist_matrix , char transform)

{
  register int i , j;


  nb_row = dist_matrix.nb_row;
  nb_column = dist_matrix.nb_column;

  row_identifier = new int[nb_row];
  for (i = 0;i < nb_row;i++) {
    row_identifier[i] = dist_matrix.row_identifier[i];
  }

  column_identifier = new int[nb_column];
  for (i = 0;i < nb_column;i++) {
    column_identifier[i] = dist_matrix.column_identifier[i];
  }

  distance = new double*[nb_row];
  length = new int*[nb_row];
  for (i = 0;i < nb_row;i++) {
    distance[i] = new double[nb_column];
    length[i] = new int[nb_column];
  }

  for (i = 0;i < nb_row;i++) {
    switch (transform) {

    case 's' : {
      distance[i][i] = dist_matrix.distance[i][i];
      length[i][i] = dist_matrix.length[i][i];
      for (j = i + 1;j < nb_column;j++) {
        if ((dist_matrix.distance[i][j] != -D_INF) && (dist_matrix.distance[j][i] != -D_INF)) {
          distance[i][j] = dist_matrix.distance[i][j] + dist_matrix.distance[j][i];
          length[i][j] = dist_matrix.length[i][j] + dist_matrix.length[j][i];
        }
        else {
          distance[i][j] = -D_INF;
          length[i][j] = 0;
        }
        distance[j][i] = distance[i][j];
        length[j][i] = length[i][j];
      }
      break;
    }

    case 'u' : {
      for (j = 0;j < nb_column;j++) {
        if ((dist_matrix.distance[i][j] != -D_INF) && (dist_matrix.length[i][j] > 0)) {
          distance[i][j] = dist_matrix.distance[i][j] / dist_matrix.length[i][j];
          length[i][j] = 1;
        }
        else {
          distance[i][j] = dist_matrix.distance[i][j];
          length[i][j] = dist_matrix.length[i][j];
        }
      }
      break;
    }

    default : {
      for (j = 0;j < nb_column;j++) {
        distance[i][j] = dist_matrix.distance[i][j];
        length[i][j] = dist_matrix.length[i][j];
      }
      break;
    }
    }
  }

  if (transform != 'u') {
    if ((dist_matrix.deletion_distance) && (dist_matrix.nb_deletion)) {
      deletion_distance = new double*[nb_row];
      nb_deletion = new int*[nb_row];
      for (i = 0;i < nb_row;i++) {
        deletion_distance[i] = new double[nb_column];
        nb_deletion[i] = new int[nb_column];
      }

      for (i = 0;i < nb_row;i++) {
        switch (transform) {

        case 's' : {
          deletion_distance[i][i] = dist_matrix.deletion_distance[i][i];
          nb_deletion[i][i] = dist_matrix.nb_deletion[i][i];
          for (j = i + 1;j < nb_column;j++) {
            if ((dist_matrix.distance[i][j] != -D_INF) && (dist_matrix.distance[j][i] != -D_INF)) {
              deletion_distance[i][j] = dist_matrix.deletion_distance[i][j] + dist_matrix.deletion_distance[j][i];
              nb_deletion[i][j] = dist_matrix.nb_deletion[i][j] + dist_matrix.nb_deletion[j][i];
            }
            else {
              deletion_distance[i][j] = 0.;
              nb_deletion[i][j] = 0;
            }
            deletion_distance[j][i] = deletion_distance[i][j];
            nb_deletion[j][i] = nb_deletion[i][j];
          }
          break;
        }

        default : {
          for (j = 0;j < nb_column;j++) {
            deletion_distance[i][j] = dist_matrix.deletion_distance[i][j];
            nb_deletion[i][j] = dist_matrix.nb_deletion[i][j];
          }
          break;
        }
        }
      }
    }

    else {
      deletion_distance = NULL;
      nb_deletion = NULL;
    }

    if ((dist_matrix.insertion_distance) && (dist_matrix.nb_insertion)) {
      insertion_distance = new double*[nb_row];
      nb_insertion = new int*[nb_row];
      for (i = 0;i < nb_row;i++) {
        insertion_distance[i] = new double[nb_column];
        nb_insertion[i] = new int[nb_column];
      }

      for (i = 0;i < nb_row;i++) {
        switch (transform) {

        case 's' : {
          insertion_distance[i][i] = dist_matrix.insertion_distance[i][i];
          nb_insertion[i][i] = dist_matrix.nb_insertion[i][i];
          for (j = i + 1;j < nb_column;j++) {
            if ((dist_matrix.distance[i][j] != -D_INF) && (dist_matrix.distance[j][i] != -D_INF)) {
              insertion_distance[i][j] = dist_matrix.insertion_distance[i][j] + dist_matrix.insertion_distance[j][i];
              nb_insertion[i][j] = dist_matrix.nb_insertion[i][j] + dist_matrix.nb_insertion[j][i];
            }
            else {
              insertion_distance[i][j] = 0.;
              nb_insertion[i][j] = 0;
            }
            insertion_distance[j][i] = insertion_distance[i][j];
            nb_insertion[j][i] = nb_insertion[i][j];
          }
          break;
        }

        default : {
          for (j = 0;j < nb_column;j++) {
            insertion_distance[i][j] = dist_matrix.insertion_distance[i][j];
            nb_insertion[i][j] = dist_matrix.nb_insertion[i][j];
          }
          break;
        }
        }
      }
    }

    else {
      insertion_distance = NULL;
      nb_insertion = NULL;
    }

    if (dist_matrix.nb_match) {
      nb_match = new int*[nb_row];
      for (i = 0;i < nb_row;i++) {
        nb_match[i] = new int[nb_column];
      }

      for (i = 0;i < nb_row;i++) {
        switch (transform) {

        case 's' : {
          nb_match[i][i] = dist_matrix.nb_match[i][i];
          for (j = i + 1;j < nb_column;j++) {
            if ((dist_matrix.distance[i][j] != -D_INF) && (dist_matrix.distance[j][i] != -D_INF)) {
              nb_match[i][j] = dist_matrix.nb_match[i][j] + dist_matrix.nb_match[j][i];
            }
            else {
              nb_match[i][j] = 0;
            }
            nb_match[j][i] = nb_match[i][j];
          }
          break;
        }

        default : {
          for (j = 0;j < nb_column;j++) {
            nb_match[i][j] = dist_matrix.nb_match[i][j];
          }
          break;
        }
        }
      }
    }

    else {
      nb_match = NULL;
    }

    if ((dist_matrix.substitution_distance) && (dist_matrix.nb_substitution)) {
      substitution_distance = new double*[nb_row];
      nb_substitution = new int*[nb_row];
      for (i = 0;i < nb_row;i++) {
        substitution_distance[i] = new double[nb_column];
        nb_substitution[i] = new int[nb_column];
      }

      for (i = 0;i < nb_row;i++) {
        switch (transform) {

        case 's' : {
          substitution_distance[i][i] = dist_matrix.substitution_distance[i][i];
          nb_substitution[i][i] = dist_matrix.nb_substitution[i][i];
          for (j = i + 1;j < nb_column;j++) {
            if ((dist_matrix.distance[i][j] != -D_INF) && (dist_matrix.distance[j][i] != -D_INF)) {
              substitution_distance[i][j] = dist_matrix.substitution_distance[i][j] + dist_matrix.substitution_distance[j][i];
              nb_substitution[i][j] = dist_matrix.nb_substitution[i][j] + dist_matrix.nb_substitution[j][i];
            }
            else {
              substitution_distance[i][j] = 0.;
              nb_substitution[i][j] = 0;
            }
            substitution_distance[j][i] = substitution_distance[i][j];
            nb_substitution[j][i] = nb_substitution[i][j];
          }
          break;
        }

        default : {
          for (j = 0;j < nb_column;j++) {
            substitution_distance[i][j] = dist_matrix.substitution_distance[i][j];
            nb_substitution[i][j] = dist_matrix.nb_substitution[i][j];
          }
          break;
        }
        }
      }
    }

    else {
      substitution_distance = NULL;
      nb_substitution = NULL;
    }

    if ((dist_matrix.transposition_distance) && (dist_matrix.nb_transposition)) {
      transposition_distance = new double*[nb_row];
      nb_transposition = new int*[nb_row];
      for (i = 0;i < nb_row;i++) {
        transposition_distance[i] = new double[nb_column];
        nb_transposition[i] = new int[nb_column];
      }

      for (i = 0;i < nb_row;i++) {
        switch (transform) {

        case 's' : {
          transposition_distance[i][i] = dist_matrix.transposition_distance[i][i];
          nb_transposition[i][i] = dist_matrix.nb_transposition[i][i];
          for (j = i + 1;j < nb_column;j++) {
            if ((dist_matrix.distance[i][j] != -D_INF) && (dist_matrix.distance[j][i] != -D_INF)) {
              transposition_distance[i][j] = dist_matrix.transposition_distance[i][j] + dist_matrix.transposition_distance[j][i];
              nb_transposition[i][j] = dist_matrix.nb_transposition[i][j] + dist_matrix.nb_transposition[j][i];
            }
            else {
              transposition_distance[i][j] = 0.;
              nb_transposition[i][j] = 0;
            }
            transposition_distance[j][i] = transposition_distance[i][j];
            nb_transposition[j][i] = nb_transposition[i][j];
          }
          break;
        }

        default : {
          for (j = 0;j < nb_column;j++) {
            transposition_distance[i][j] = dist_matrix.transposition_distance[i][j];
            nb_transposition[i][j] = dist_matrix.nb_transposition[i][j];
          }
          break;
        }
        }
      }
    }

    else {
      transposition_distance = NULL;
      nb_transposition = NULL;
    }
  }

  else {
    deletion_distance = NULL;
    nb_deletion = NULL;
    insertion_distance = NULL;
    nb_insertion = NULL;
    nb_match = NULL;
    substitution_distance = NULL;
    nb_substitution = NULL;
    transposition_distance = NULL;
    nb_transposition = NULL;
  }

  label_size = dist_matrix.label_size;
  label = new char[label_size];
  strcpy(label , dist_matrix.label);
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet DistanceMatrix.
 *
 *--------------------------------------------------------------*/

void DistanceMatrix::remove()

{
  register int i;


  delete [] row_identifier;
  delete [] column_identifier;

  if (distance) {
    for (i = 0;i < nb_row;i++) {
      delete [] distance[i];
    }
    delete [] distance;
  }

  if (length) {
    for (i = 0;i < nb_row;i++) {
      delete [] length[i];
    }
    delete [] length;
  }

  if (deletion_distance) {
    for (i = 0;i < nb_row;i++) {
      delete [] deletion_distance[i];
    }
    delete [] deletion_distance;
  }

  if (nb_deletion) {
    for (i = 0;i < nb_row;i++) {
      delete [] nb_deletion[i];
    }
    delete [] nb_deletion;
  }

  if (insertion_distance) {
    for (i = 0;i < nb_row;i++) {
      delete [] insertion_distance[i];
    }
    delete [] insertion_distance;
  }

  if (nb_insertion) {
    for (i = 0;i < nb_row;i++) {
      delete [] nb_insertion[i];
    }
    delete [] nb_insertion;
  }

  if (nb_match) {
    for (i = 0;i < nb_row;i++) {
      delete [] nb_match[i];
    }
    delete [] nb_match;
  }

  if (substitution_distance) {
    for (i = 0;i < nb_row;i++) {
      delete [] substitution_distance[i];
    }
    delete [] substitution_distance;
  }

  if (nb_substitution) {
    for (i = 0;i < nb_row;i++) {
      delete [] nb_substitution[i];
    }
    delete [] nb_substitution;
  }

  if (transposition_distance) {
    for (i = 0;i < nb_row;i++) {
      delete [] transposition_distance[i];
    }
    delete [] transposition_distance;
  }

  if (nb_transposition) {
    for (i = 0;i < nb_row;i++) {
      delete [] nb_transposition[i];
    }
    delete [] nb_transposition;
  }

  delete [] label;
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe DistanceMatrix.
 *
 *--------------------------------------------------------------*/

DistanceMatrix::~DistanceMatrix()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe DistanceMatrix.
 *
 *  argument : reference sur un objet DistanceMatrix.
 *
 *--------------------------------------------------------------*/

DistanceMatrix& DistanceMatrix::operator=(const DistanceMatrix &dist_matrix)

{
  if (&dist_matrix != this) {
    remove();
    copy(dist_matrix);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Selection de formes par l'identificateur.
 *
 *  arguments : reference sur un objet StatError, nombre de formes,
 *              identificateurs des formes, flag pour conserver ou rejeter
 *              les formes selectionnees.
 *
 *--------------------------------------------------------------*/

DistanceMatrix* DistanceMatrix::select_individual(StatError &error , int inb_pattern ,
                                                  int *iidentifier , bool keep) const

{
  bool status = true , *selected_pattern;
  register int i , j;
  int nb_pattern , max_identifier , *identifier;
  DistanceMatrix *dist_matrix;


  dist_matrix = NULL;
  error.init();

  nb_pattern = MAX(nb_row , nb_column);

  if (nb_pattern <= 2) {
    status = false;
    error.update(STAT_error[STATR_MATRIX_DIMENSIONS]);
  }

  if ((inb_pattern < 1) || (inb_pattern > (keep ? nb_pattern : nb_pattern - 1))) {
    status = false;
    ostringstream error_message;
    error_message << STAT_error[STATR_NUMBER] << " " << label << "s";
    error.update((error_message.str()).c_str());
  }

  if (status) {
    max_identifier = 1;
    for (i = 0;i < inb_pattern;i++) {
      if (iidentifier[i] > max_identifier) {
        max_identifier = iidentifier[i];
      }
    }

    selected_pattern = new bool[max_identifier + 1];
    for (i = 0;i <= max_identifier;i++) {
      selected_pattern[i] = false;
    }

    if (nb_row > 1) {
      identifier = row_identifier;
    }
    else {
      identifier = column_identifier;
    }

    for (i = 0;i < inb_pattern;i++) {
      for (j = 0;j < nb_pattern;j++) {
        if (iidentifier[i] == identifier[j]) {
          break;
        }
      }

      if (j == nb_pattern) {
        status = false;
        ostringstream error_message;
        error_message << iidentifier[i] << ": " << STAT_error[STATR_BAD] << " "
                      << label << " " << STAT_label[STATL_IDENTIFIER];
        error.update((error_message.str()).c_str());
      }

      else if (selected_pattern[iidentifier[i]]) {
        status = false;
        ostringstream error_message;
        error_message << label << " " << iidentifier[i] << " "
                      << STAT_error[STATR_ALREADY_SELECTED];
        error.update((error_message.str()).c_str());
      }
      else {
        selected_pattern[iidentifier[i]] = true;
      }
    }

    delete [] selected_pattern;
  }

  if (status) {
    dist_matrix = new DistanceMatrix(*this , inb_pattern , iidentifier , keep);
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*
 *
 *  Symetrisation d'une matrice des distances.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

DistanceMatrix* DistanceMatrix::symmetrize(StatError &error) const

{
  DistanceMatrix *dist_matrix;


  error.init();

  if (test_symmetry()) {
    dist_matrix = NULL;
    error.update(STAT_error[STATR_SYMMETRICAL_MATRIX]);
  }
  else {
    dist_matrix = new DistanceMatrix(*this , 's');
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*
 *
 *  Denormalisation d'une matrice des distances.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

DistanceMatrix* DistanceMatrix::unnormalize(StatError &error) const

{
  bool status = false;
  register int i , j;
  DistanceMatrix *dist_matrix;


  error.init();
  dist_matrix = NULL;

  for (i = 0;i < nb_row;i++) {
    for (j = 0;j < nb_column;j++) {
      if ((distance[i][j] != -D_INF) && (length[i][j] > 1)) {
        status = true;
        break;
      }
    }
    if (status == true) {
      break;
    }
  }

  if (!status) {
    error.update(STAT_error[STATR_UNNORMALIZED_DISSIMILARITY_MEASURES]);
  }

  if (!test_symmetry()) {
    status = false;
    error.update(STAT_error[STATR_UNSYMMETRICAL_MATRIX]);
  }

  if (status) {
    dist_matrix = new DistanceMatrix(*this , 'u');
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet DistanceMatrix.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& DistanceMatrix::line_write(ostream &os) const

{
  os << STAT_label[STATL_NB_ROW] << ": " << nb_row << "   "
     << STAT_label[STATL_NB_COLUMN] << ": " << nb_column;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des proprietes d'une matrice de distances.
 *
 *  arguments : pointeur sur les distances mormalisees, stream,
 *              format de sortie ('a' : ASCII, 's' : Spreadsheet).
 *
 *--------------------------------------------------------------*/

ostream& DistanceMatrix::property_print(double **normalized_distance ,
                                        ostream &os , char format) const

{
  bool status = true;
  register int i , j , k;
  int nb_combination , nb_symmetry , nb_triangle_inequality;


  if ((nb_row > 1) && (nb_row == nb_column)) {
    for (i = 0;i < nb_row;i++) {
      if (row_identifier[i] != column_identifier[i]) {
        status = false;
        break;
      }
    }
  }

  if (status) {

    // symetrie

    nb_combination = 0;
    nb_symmetry = 0;

    for (i = 0;i < nb_row;i++) {
      for (j = i + 1;j < nb_row;j++) {
        if ((normalized_distance[i][j] != -D_INF) && (normalized_distance[j][i] != -D_INF)) {
          nb_combination++;

          if ((normalized_distance[i][j] < normalized_distance[j][i] - DISTANCE_ROUNDNESS) ||
              (normalized_distance[i][j] > normalized_distance[j][i] + DISTANCE_ROUNDNESS)) {
            nb_symmetry++;

            switch (format) {
            case 'a' :
              os << label << " " << row_identifier[i] << ", "
                 << label << " " << row_identifier[j] << ": "
                 << STAT_error[STATR_SYMMETRY] << endl;
              break;
            case 's' :
              os << label << " " << row_identifier[i] << "\t"
                 << label << " " << row_identifier[j] << "\t"
                 << STAT_error[STATR_SYMMETRY] << endl;
              break;
            }
          }
        }
      }
    }

    switch (format) {
    case 'a' :
      os << "\n" << STAT_label[STATL_SYMMETRY_RATE] << ": "
         << (double)(nb_combination - nb_symmetry) / (double)nb_combination << endl;
      break;
    case 's' :
      os << "\n" << STAT_label[STATL_SYMMETRY_RATE] << "\t"
         << (double)(nb_combination - nb_symmetry) / (double)nb_combination << endl;
      break;
    }

    // inegalite triangulaire

    if (nb_row > 2) {
      nb_combination = 0;
      nb_triangle_inequality = 0;

      for (i = 0;i < nb_row;i++) {
        for (j = (nb_symmetry == 0 ? i + 1 : 0);j < nb_row;j++) {
          for (k = 0;k < nb_row;k++) {
            if ((i != j) && (k != i) && (k != j) && (normalized_distance[i][k] != -D_INF) &&
                (normalized_distance[k][j] != -D_INF) && (normalized_distance[i][j] != -D_INF)) {
//            if ((i != j) && (k != i) && (k != j) && (distance[i][k] != -D_INF) &&
//                (distance[k][j] != -D_INF) && (distance[i][j] != -D_INF)) {
              nb_combination++;

              if (normalized_distance[i][k] + normalized_distance[k][j] < normalized_distance[i][j]) {
//              if ((distance[i][k] + distance[k][j]) / (length[i][k] + length[k][j]) < distance[i][j] / length[i][j]) {
                nb_triangle_inequality++;

                switch (format) {
                case 'a' :
                  os << label << " " << row_identifier[i] << ", "
                     << label << " " << row_identifier[k] << ", "
                     << label << " " << row_identifier[j] << ": "
                     << STAT_error[STATR_TRIANGLE_INEQUALITY] << endl;
                  break;
                case 's' :
                  os << label << " " << row_identifier[i] << "\t"
                     << label << " " << row_identifier[k] << "\t"
                     << label << " " << row_identifier[j] << "\t"
                     << STAT_error[STATR_TRIANGLE_INEQUALITY] << endl;
                  break;
                }
              }
            }
          }
        }
      }

      switch (format) {
      case 'a' :
        os << "\n" << STAT_label[STATL_TRIANGLE_INEQUALITY_RATE] << ": "
           << (double)(nb_combination - nb_triangle_inequality) / (double)nb_combination << endl;
        break;
      case 's' :
        os << "\n" << STAT_label[STATL_TRIANGLE_INEQUALITY_RATE] << "\t"
           << (double)(nb_combination - nb_triangle_inequality) / (double)nb_combination << endl;
        break;
      }
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Cumul des termes d'une matrice entiere.
 *
 *  arguments : nombres de lignes et de colonnes, pointeur sur les valeurs.
 *
 *--------------------------------------------------------------*/

int cumul_computation(int nb_row , int nb_column , int **value)

{
  register int i , j;
  int cumul = 0;


  for (i = 0;i < nb_row;i++) {
    for (j = 0;j < nb_column;j++) {
      cumul += value[i][j];
    }
  }

  return cumul;
}


/*--------------------------------------------------------------*
 *
 *  Cumul des termes d'un vecteur de distances.
 *
 *  arguments : dimension, pointeur sur les distances.
 *
 *--------------------------------------------------------------*/

double cumul_distance_computation(int dim , double *distance)

{
  register int i;
  double cumul_distance = 0.;


  for (i = 0;i < dim;i++) {
    if (distance[i] != -D_INF) {
      cumul_distance += distance[i];
    }
    else {
      cumul_distance = -D_INF;
      break;
    }
  }

  return cumul_distance;
}


/*--------------------------------------------------------------*
 *
 *  Tri des formes par distance croissante.
 *
 *  arguments : nombre de formes, distances, nombre de formes triees.
 *
 *--------------------------------------------------------------*/

int* pattern_sort(int nb_pattern , double *distance , int nb_sorted_pattern)

{
  bool *selected_pattern;
  register int i , j;
  int *index;
  double min_distance;


  if (nb_sorted_pattern == I_DEFAULT) {
    nb_sorted_pattern = nb_pattern;
  }

  index = new int[nb_sorted_pattern];
  selected_pattern = new bool[nb_pattern];

  for (i = 0;i < nb_pattern;i++) {
    selected_pattern[i] = false;
  }

  for (i = 0;i < nb_sorted_pattern;i++) {
    min_distance = -D_INF * 10;
    for (j = 0;j < nb_pattern;j++) {
      if ((!selected_pattern[j]) && (distance[j] < min_distance)) {
        min_distance = distance[j];
        index[i] = j;
      }
    }

    selected_pattern[index[i]] = true;
  }

  delete [] selected_pattern;

  return index;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet DistanceMatrix.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& DistanceMatrix::ascii_write(ostream &os , bool exhaustive) const

{
  register int i , j;
  int dim , max_identifier , buff , cumul_nb_deletion , cumul_nb_insertion , cumul_nb_match ,
      cumul_nb_substitution , cumul_nb_transposition , cumul_nb_operation , *cumul_length ,
      *index , width[2];
  long old_adjust;
  double cumul_deletion_distance , cumul_insertion_distance , cumul_substitution_distance ,
         cumul_transposition_distance , *distance_vector , *cumul_distance , **normalized_distance;


  old_adjust = os.setf(ios::right , ios::adjustfield);

  normalized_distance = new double*[nb_row];
  for (i = 0;i < nb_row;i++) {
    normalized_distance[i] = new double[nb_column];
    for (j = 0;j < nb_column;j++) {
      normalized_distance[i][j] = distance[i][j];
      if ((distance[i][j] != -D_INF) && (length[i][j] > 1)) {
        normalized_distance[i][j] /= length[i][j];
      }
    }
  }

  // tri par distance croissante

  dim = MAX(nb_row , nb_column);
  distance_vector = new double[dim];

  if (nb_row == 1) {
    for (i = 0;i < nb_column;i++) {
      distance_vector[i] = normalized_distance[0][i];
    }
  }

  else if (nb_column == 1) {
    for (i = 0;i < nb_row;i++) {
      distance_vector[i] = normalized_distance[i][0];
    }
  }

  else if (nb_row == nb_column) {
    cumul_distance = new double[nb_row];
    cumul_length = new int[nb_row];

    for (i = 0;i < nb_row;i++) {
      cumul_distance[i] = stat_tool::cumul_distance_computation(nb_column , distance[i]);

      if (cumul_distance[i] == -D_INF) {
        distance_vector[i] = -D_INF;
      }
      else {
        cumul_length[i] = cumul_computation(1 , nb_column , length + i);
        distance_vector[i] = cumul_distance[i] / cumul_length[i];
      }
    }
  }

  index = pattern_sort(dim , distance_vector , I_DEFAULT);

  // calcul des largeurs des colonnes

  max_identifier = 0;
  for (i = 0;i < nb_row;i++) {
    if (row_identifier[i] > max_identifier) {
      max_identifier = row_identifier[i];
    }
  }
  width[0] = column_width(max_identifier);

  max_identifier = 0;
  for (i = 0;i < nb_column;i++) {
    if (column_identifier[i] > max_identifier) {
      max_identifier = column_identifier[i];
    }
  }
  width[1] = column_width(max_identifier);

  for (i = 0;i < nb_row;i++) {
    buff = column_width(nb_column , normalized_distance[i]);
    if (buff > width[1]) {
      width[1] = buff;
    }
  }
  width[1] += ASCII_SPACE;

  // ecriture de la matrice des distances

  os << STAT_label[STATL_DISTANCE] << " " << STAT_label[STATL_MATRIX] << endl;

  if (nb_row == nb_column) {
    os << "\n" << setw(width[0] + width[1]) << column_identifier[0];
    for (i = 1;i < nb_column;i++) {
      os << setw(width[1]) << column_identifier[i];
    }
    os << endl;

    for (i = 0;i < nb_row;i++) {
      os << setw(width[0]) << row_identifier[i];
      for (j = 0;j < nb_column;j++) {
        os << setw(width[1]) << normalized_distance[i][j];
      }
      os << endl;
    }
  }

  else {
    os << "\n" << setw(width[0] + width[1]) << column_identifier[nb_column == 1 ? 0 : index[0]];
    for (i = 1;i < nb_column;i++) {
      os << setw(width[1]) << column_identifier[index[i]];
    }
    os << endl;

    for (i = 0;i < nb_row;i++) {
      os << setw(width[0]) << row_identifier[nb_row == 1 ? i : index[i]];
      for (j = 0;j < nb_column;j++) {
        os << setw(width[1]) << normalized_distance[nb_row == 1 ? i : index[i]][nb_column == 1 ? j : index[j]];
      }
      os << endl;
    }
  }

  if ((nb_row > 1) && (nb_row == nb_column)) {
    if (exhaustive) {
      property_print(normalized_distance , os , 'a');
    }

    os << "\n" << STAT_label[STATL_CUMUL_DISTANCE] << " " << STAT_label[STATL_REFERENCE] << " " << label << endl;
    for (i = 0;i < nb_row;i++) {
      if (distance_vector[index[i]] == -D_INF) {
        os << setw(width[0]) << row_identifier[index[i]] << "  "
           << STAT_label[STATL_DISTANCE] << ": " << distance_vector[index[i]] << endl;
      }

      else {
        if ((nb_deletion) && (nb_insertion) && (nb_match)) {
          cumul_nb_deletion = cumul_computation(1 , nb_column , nb_deletion + index[i]);
          cumul_nb_insertion = cumul_computation(1 , nb_column , nb_insertion + index[i]);
          cumul_nb_match = cumul_computation(1 , nb_column , nb_match + index[i]);

          if (nb_substitution) {
            cumul_nb_substitution = cumul_computation(1 , nb_column , nb_substitution + index[i]);
          }
          else {
            cumul_nb_substitution = 0;
          }

          if (nb_transposition) {
            cumul_nb_transposition = cumul_computation(1 , nb_column , nb_transposition + index[i]);
          }
          else {
            cumul_nb_transposition = 0;
          }

          cumul_nb_operation = cumul_nb_deletion + cumul_nb_insertion + cumul_nb_match +
                               cumul_nb_substitution + cumul_nb_transposition;
        }

        os << setw(width[0]) << row_identifier[index[i]] << "  " << STAT_label[STATL_DISTANCE]
           << " (" << STAT_label[STATL_LENGTH] << "): " << distance_vector[index[i]] << " (";

        if ((nb_deletion) && (nb_insertion) && (nb_match)) {
          os << cumul_nb_operation;
        }
        else {
          os << cumul_length[index[i]];
        }
        os << ")";

        if ((deletion_distance) && (insertion_distance)) {
          cumul_deletion_distance = stat_tool::cumul_distance_computation(nb_column , deletion_distance[index[i]]);
          cumul_insertion_distance = stat_tool::cumul_distance_computation(nb_column , insertion_distance[index[i]]);
          os << " = " << cumul_deletion_distance / cumul_length[index[i]]
             << " (" << cumul_nb_deletion << " d) + "
             << cumul_insertion_distance / cumul_length[index[i]]
             << " (" << cumul_nb_insertion << " i) + 0 (" << cumul_nb_match << " m)";

          if (substitution_distance) {
            cumul_substitution_distance = stat_tool::cumul_distance_computation(nb_column , substitution_distance[index[i]]);
            os << " + " << cumul_substitution_distance / cumul_length[index[i]]
               << " (" << cumul_nb_substitution << " s)";
          }

          if (transposition_distance) {
            cumul_transposition_distance = stat_tool::cumul_distance_computation(nb_column , transposition_distance[index[i]]);
            os << " + " << cumul_transposition_distance / cumul_length[index[i]]
               << " (" << cumul_nb_transposition << " t)";
          }
        }

        os << endl;
      }
    }

    delete [] cumul_distance;
    delete [] cumul_length;
  }

  if ((exhaustive) && (nb_row > 1) && (nb_row == nb_column) && (nb_row <= ASCII_NB_INDIVIDUAL) &&
      (nb_deletion) && (nb_insertion) && (nb_match)) {
    os << "\n" << label << " " << STAT_label[STATL_DISTANCE] << endl;
    for (i = 0;i < nb_row;i++) {
//      for (j = 0;j <= i;j++) {
      for (j = 0;j < nb_column;j++) {
        os << setw(width[0]) << row_identifier[i] << " " << setw(width[0]) << column_identifier[j] << "  "
           << STAT_label[STATL_DISTANCE] << " (" << STAT_label[STATL_LENGTH] << "): "
           << normalized_distance[i][j];

        if ((normalized_distance[i][j] != 0.) && (normalized_distance[i][j] != -D_INF)) {
          cumul_nb_operation = nb_deletion[i][j] + nb_insertion[i][j] + nb_match[i][j];
          if (nb_substitution) {
            cumul_nb_operation += nb_substitution[i][j];
          }
          if (nb_transposition) {
            cumul_nb_operation += nb_transposition[i][j];
          }
          os << " (" << cumul_nb_operation << ")";

          os << " = " << deletion_distance[i][j] / length[i][j]
             << " (" << nb_deletion[i][j] << " d) + "
             << insertion_distance[i][j] / length[i][j]
             << " (" << nb_insertion[i][j] << " i) + 0 (" << nb_match[i][j] << " m)";

          if (substitution_distance) {
            os << " + " << substitution_distance[i][j] / length[i][j]
               << " (" << nb_substitution[i][j] << " s)";
          }

          if (transposition_distance) {
            os << " + " << transposition_distance[i][j] / length[i][j]
               << " (" << nb_transposition[i][j] << " t)";
          }
        }
        os << endl;
      }
    }
  }

  for (i = 0;i < nb_row;i++) {
    delete [] normalized_distance[i];
  }
  delete [] normalized_distance;

  delete [] distance_vector;
  delete [] index;

  if ((nb_deletion) && (nb_insertion) && (nb_match)) {
    cumul_nb_deletion = cumul_computation(nb_row , nb_column , nb_deletion);
    cumul_nb_insertion = cumul_computation(nb_row , nb_column , nb_insertion);
    cumul_nb_match = cumul_computation(nb_row , nb_column , nb_match);

    if (nb_substitution) {
      cumul_nb_substitution = cumul_computation(nb_row , nb_column , nb_substitution);
    }
    else {
      cumul_nb_substitution = 0;
    }

    if (nb_transposition) {
      cumul_nb_transposition = cumul_computation(nb_row , nb_column , nb_transposition);
    }
    else {
      cumul_nb_transposition = 0;
    }

    cumul_nb_operation = cumul_nb_deletion + cumul_nb_insertion + cumul_nb_match +
                         cumul_nb_substitution + cumul_nb_transposition;

    if (cumul_nb_operation > 0) {
      os << "\n" << STAT_label[STATL_DELETION] << " " << STAT_label[STATL_RATE] << ": "
         << (double)cumul_nb_deletion / (double)cumul_nb_operation << endl;
      os << STAT_label[STATL_INSERTION] << " " << STAT_label[STATL_RATE] << ": "
         << (double)cumul_nb_insertion / (double)cumul_nb_operation << endl;
      os << STAT_label[STATL_MATCH] << " " << STAT_label[STATL_RATE] << ": "
         << (double)cumul_nb_match / (double)cumul_nb_operation << endl;
      if (nb_substitution) {
        os << STAT_label[STATL_SUBSTITUTION] << " " << STAT_label[STATL_RATE] << ": "
           << (double)cumul_nb_substitution / (double)cumul_nb_operation << endl;
      }
      if (nb_transposition) {
        os << STAT_label[STATL_TRANSPOSITION] << " " << STAT_label[STATL_RATE] << ": "
           << (double)cumul_nb_transposition / (double)cumul_nb_operation << endl;
      }
    }
  }

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet DistanceMatrix dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool DistanceMatrix::ascii_write(StatError &error , const char *path ,
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
    ascii_write(out_file , exhaustive);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet DistanceMatrix au format tableur.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& DistanceMatrix::spreadsheet_write(ostream &os) const

{
  register int i , j;
  int dim , cumul_nb_deletion , cumul_nb_insertion , cumul_nb_match ,
      cumul_nb_substitution , cumul_nb_transposition , cumul_nb_operation ,
      *cumul_length , *index;
  double cumul_deletion_distance , cumul_insertion_distance , cumul_substitution_distance ,
         cumul_transposition_distance , *distance_vector , *cumul_distance , **normalized_distance;


  normalized_distance = new double*[nb_row];
  for (i = 0;i < nb_row;i++) {
    normalized_distance[i] = new double[nb_column];
    for (j = 0;j < nb_column;j++) {
      normalized_distance[i][j] = distance[i][j];
      if ((distance[i][j] != -D_INF) && (length[i][j] > 1)) {
        normalized_distance[i][j] /= length[i][j];
      }
    }
  }

  // tri par distance croissante

  dim = MAX(nb_row , nb_column);
  distance_vector = new double[dim];

  if (nb_row == 1) {
    for (i = 0;i < nb_column;i++) {
      distance_vector[i] = normalized_distance[0][i];
    }
  }

  else if (nb_column == 1) {
    for (i = 0;i < nb_row;i++) {
      distance_vector[i] = normalized_distance[i][0];
    }
  }

  else if (nb_row == nb_column) {
    cumul_distance = new double[nb_row];
    cumul_length = new int[nb_row];

    for (i = 0;i < nb_row;i++) {
      cumul_distance[i] = stat_tool::cumul_distance_computation(nb_column , distance[i]);

      if (cumul_distance[i] == -D_INF) {
        distance_vector[i] = -D_INF;
      }
      else {
        cumul_length[i] = cumul_computation(1 , nb_column , length + i);
        distance_vector[i] = cumul_distance[i] / cumul_length[i];
      }
    }
  }

  index = pattern_sort(dim , distance_vector , I_DEFAULT);

  // ecriture de la matrice des distances

  os << STAT_label[STATL_DISTANCE] << " " << STAT_label[STATL_MATRIX] << endl;

  if (nb_row == nb_column) {
    os << "\n";
    for (i = 0;i < nb_column;i++) {
      os << "\t" << column_identifier[i];
    }
    os << endl;

    for (i = 0;i < nb_row;i++) {
      os << row_identifier[i];
      for (j = 0;j < nb_column;j++) {
        os << "\t" << normalized_distance[i][j];
      }
      os << endl;
    }
  }

  else {
    os << "\n";
    for (i = 0;i < nb_column;i++) {
      os << "\t" << column_identifier[nb_column == 1 ? i : index[i]];
    }
    os << endl;

    for (i = 0;i < nb_row;i++) {
      os << row_identifier[nb_row == 1 ? i : index[i]];
      for (j = 0;j < nb_column;j++) {
        os << "\t" << normalized_distance[nb_row == 1 ? i : index[i]][nb_column == 1 ? j : index[j]];
      }
      os << endl;
    }
  }

  if ((nb_row > 1) && (nb_row == nb_column)) {
    os << "\n" << STAT_label[STATL_CUMUL_DISTANCE] << " " << STAT_label[STATL_REFERENCE] << " " << label << endl;
    for (i = 0;i < nb_row;i++) {
      if (distance_vector[index[i]] == -D_INF) {
        os << row_identifier[index[i]] << "\t"
           << STAT_label[STATL_DISTANCE] << "\t" << distance_vector[index[i]] << endl;
      }

      else {
        if ((nb_deletion) && (nb_insertion) && (nb_match)) {
          cumul_nb_deletion = cumul_computation(1 , nb_column , nb_deletion + index[i]);
          cumul_nb_insertion = cumul_computation(1 , nb_column , nb_insertion + index[i]);
          cumul_nb_match = cumul_computation(1 , nb_column , nb_match + index[i]);

          if (nb_substitution) {
            cumul_nb_substitution = cumul_computation(1 , nb_column , nb_substitution + index[i]);
          }
          else {
            cumul_nb_substitution = 0;
          }

          if (nb_transposition) {
            cumul_nb_transposition = cumul_computation(1 , nb_column , nb_transposition + index[i]);
          }
          else {
            cumul_nb_transposition = 0;
          }

          cumul_nb_operation = cumul_nb_deletion + cumul_nb_insertion + cumul_nb_match +
                               cumul_nb_substitution + cumul_nb_transposition;
        }

        os << row_identifier[index[i]] << "\t" << STAT_label[STATL_DISTANCE] << " ("
           << STAT_label[STATL_LENGTH] << ")\t" << distance_vector[index[i]] << " (";

        if ((nb_deletion) && (nb_insertion) && (nb_match)) {
          os << cumul_nb_operation;
        }
        else {
          os << cumul_length[index[i]];
        }
        os << ")";

        if ((deletion_distance) && (insertion_distance)) {
          cumul_deletion_distance = stat_tool::cumul_distance_computation(nb_column , deletion_distance[index[i]]);
          cumul_insertion_distance = stat_tool::cumul_distance_computation(nb_column , insertion_distance[index[i]]);
          os << "\t" << cumul_deletion_distance / cumul_length[index[i]]
             << " (" << cumul_nb_deletion << " d)\t"
             << cumul_insertion_distance / cumul_length[index[i]]
             << " (" << cumul_nb_insertion << " i)\t0 (" << cumul_nb_match << " m)";

          if (substitution_distance) {
            cumul_substitution_distance = stat_tool::cumul_distance_computation(nb_column , substitution_distance[index[i]]);
            os << "\t" << cumul_substitution_distance / cumul_length[index[i]]
               << " (" << cumul_nb_substitution << " s)";
          }

          if (transposition_distance) {
            cumul_transposition_distance = stat_tool::cumul_distance_computation(nb_column , transposition_distance[index[i]]);
            os << "\t" << cumul_transposition_distance / cumul_length[index[i]]
               << " (" << cumul_nb_transposition << " s)";
          }
        }

        os << endl;
      }
    }

    delete [] cumul_distance;
    delete [] cumul_length;
  }

  for (i = 0;i < nb_row;i++) {
    delete [] normalized_distance[i];
  }
  delete [] normalized_distance;

  delete [] distance_vector;
  delete [] index;

  if ((nb_deletion) && (nb_insertion) && (nb_match)) {
    cumul_nb_deletion = cumul_computation(nb_row , nb_column , nb_deletion);
    cumul_nb_insertion = cumul_computation(nb_row , nb_column , nb_insertion);
    cumul_nb_match = cumul_computation(nb_row , nb_column , nb_match);

    if (nb_substitution) {
      cumul_nb_substitution = cumul_computation(nb_row , nb_column , nb_substitution);
    }
    else {
      cumul_nb_substitution = 0;
    }

    if (nb_transposition) {
      cumul_nb_transposition = cumul_computation(nb_row , nb_column , nb_transposition);
    }
    else {
      cumul_nb_transposition = 0;
    }

    cumul_nb_operation = cumul_nb_deletion + cumul_nb_insertion + cumul_nb_match +
                         cumul_nb_substitution + cumul_nb_transposition;

    if (cumul_nb_operation > 0) {
      os << "\n" << STAT_label[STATL_DELETION] << " " << STAT_label[STATL_RATE] << "\t"
         << (double)cumul_nb_deletion / (double)cumul_nb_operation << endl;
      os << STAT_label[STATL_INSERTION] << " " << STAT_label[STATL_RATE] << "\t"
         << (double)cumul_nb_insertion / (double)cumul_nb_operation << endl;
      os << STAT_label[STATL_MATCH] << " " << STAT_label[STATL_RATE] << "\t"
         << (double)cumul_nb_match / (double)cumul_nb_operation << endl;
      if (nb_substitution) {
        os << STAT_label[STATL_SUBSTITUTION] << " " << STAT_label[STATL_RATE] << "\t"
           << (double)cumul_nb_substitution / (double)cumul_nb_operation << endl;
      }
      if (nb_transposition) {
        os << STAT_label[STATL_TRANSPOSITION] << " " << STAT_label[STATL_RATE] << "\t"
           << (double)cumul_nb_transposition / (double)cumul_nb_operation << endl;
      }
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet DistanceMatrix dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool DistanceMatrix::spreadsheet_write(StatError &error , const char *path) const

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
    spreadsheet_write(out_file);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet DistanceMatrix.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool DistanceMatrix::plot_write(StatError &error , const char *prefix ,
                                const char *title) const

{
  bool status = true;


  error.init();

  if ((nb_row == 1) && (nb_column == 1)) {
    status = false;
    error.update(STAT_error[STATR_MATRIX_DIMENSIONS]);
  }

  else {
    register int i , j;
    int plot_nb_pattern , *index;
    double *plot_distance;
    ostringstream data_file_name;


    plot_nb_pattern = MAX(nb_row , nb_column);
    plot_distance = new double[plot_nb_pattern];

    // tri par distance croissante

    if (nb_row == 1) {
      for (i = 0;i < nb_column;i++) {
        plot_distance[i] = distance[0][i];
        if ((distance[0][i] != -D_INF) && (length[0][i] > 1)) {
          plot_distance[i] /= length[0][i];
        }
      }
    }

    else if (nb_column == 1) {
      for (i = 0;i < nb_row;i++) {
        plot_distance[i] = distance[i][0];
        if ((distance[i][0] != -D_INF) && (length[i][0] > 1)) {
          plot_distance[i] /= length[i][0];
        }
      }
    }

    else if (nb_row == nb_column) {
      for (i = 0;i < nb_row;i++) {
        plot_distance[i] = stat_tool::cumul_distance_computation(nb_column , distance[i]);
        if (plot_distance[i] != -D_INF) {
          plot_distance[i] /= cumul_computation(1 , nb_column , length + i);
        }
      }
    }

    index = pattern_sort(plot_nb_pattern , plot_distance , I_DEFAULT);

    if (plot_distance[index[0]] == -D_INF) {
      status = false;
      error.update(STAT_error[STATR_INFINITE_DISTANCES]);
    }

    else {

      // ecriture du fichier de donnees

      data_file_name << prefix << ".dat";
      ofstream out_data_file((data_file_name.str()).c_str());

      if (!out_data_file) {
        status = false;
        error.update(STAT_error[STATR_FILE_PREFIX]);
      }

      else {
        for (i = 0;i < plot_nb_pattern;i++) {
          if (plot_distance[index[i]] == -D_INF) {
            plot_nb_pattern = i;
            break;
          }
          else {
            out_data_file << i + 1 << " " << plot_distance[index[i]] << endl;
          }
        }

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
            file_name[1] << stat_tool::label(prefix) << ".ps";
            out_file << "set output \"" << file_name[1].str() << "\"\n\n";
          }

          out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                   << "set title";
          if (title) {
            out_file << " \"" << title << "\"";
          }
          out_file << "\n\n";

          for (j = 0;j < plot_nb_pattern;j++) {
            if (nb_row > 1) {
              out_file << "set label \"" << row_identifier[index[j]] << "\" at "
                       << j + 1 << ", " << plot_distance[index[j]] << endl;
            }
            else {
              out_file << "set label \"" << column_identifier[index[j]] << "\" at "
                       << j + 1 << ", " << plot_distance[index[j]] << endl;
            }
          }
          out_file << endl;

          if (plot_nb_pattern < TIC_THRESHOLD) {
            out_file << "set xtics 1,1" << endl;
          }

          out_file << "plot [1:" << plot_nb_pattern << "] ["
                   << plot_distance[index[0]] * (1. - PLOT_YMARGIN) << ":"
                   << plot_distance[index[plot_nb_pattern - 1]] * (1. + PLOT_YMARGIN) << "] \""
                   << stat_tool::label((data_file_name.str()).c_str()) << "\" title \"";
          if (nb_row == 1) {
            out_file << STAT_label[STATL_REFERENCE] << " " << label << " " << row_identifier[0];
          }
          else if (nb_column == 1) {
            out_file << STAT_label[STATL_TEST] << " " << label << " " << column_identifier[0];
          }
          else if (nb_row == nb_column) {
            out_file << label;
          }
          out_file << "\" with linespoints" << endl;

          if (plot_nb_pattern < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }

          out_file << "\nunset label" << endl;

          if (i == 1) {
            out_file << "\nset terminal x11" << endl;
          }

          out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
        }
      }
    }

    delete [] plot_distance;
    delete [] index;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet DistanceMatrix.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* DistanceMatrix::get_plotable(StatError &error) const

{
  MultiPlotSet *plot_set;


  error.init();

  if ((nb_row == 1) && (nb_column == 1)) {
    plot_set = NULL;
    error.update(STAT_error[STATR_MATRIX_DIMENSIONS]);
  }

  else {
    register int i;
    int plot_nb_pattern , *index;
    double *plot_distance;
    ostringstream title , legend , identifier;


    // tri par distance croissante

    plot_nb_pattern = MAX(nb_row , nb_column);
    plot_distance = new double[plot_nb_pattern];

    if (nb_row == 1) {
      for (i = 0;i < nb_column;i++) {
        plot_distance[i] = distance[0][i];
        if ((distance[0][i] != -D_INF) && (length[0][i] > 1)) {
          plot_distance[i] /= length[0][i];
        }
      }
    }

    else if (nb_column == 1) {
      for (i = 0;i < nb_row;i++) {
        plot_distance[i] = distance[i][0];
        if ((distance[i][0] != -D_INF) && (length[i][0] > 1)) {
          plot_distance[i] /= length[i][0];
        }
      }
    }

    else if (nb_row == nb_column) {
      for (i = 0;i < nb_row;i++) {
        plot_distance[i] = stat_tool::cumul_distance_computation(nb_column , distance[i]);
        if (plot_distance[i] != -D_INF) {
          plot_distance[i] /= cumul_computation(1 , nb_column , length + i);
        }
      }
    }

    index = pattern_sort(plot_nb_pattern , plot_distance , I_DEFAULT);

    if (plot_distance[index[0]] == -D_INF) {
      plot_set = NULL;
      error.update(STAT_error[STATR_INFINITE_DISTANCES]);
    }

    else {
      plot_set = new MultiPlotSet(1);
      MultiPlotSet &plot = *plot_set;

      title.str("");
      title << STAT_label[STATL_DISTANCE] << " " << STAT_label[STATL_MATRIX];
      plot.title = title.str();

      plot.border = "15 lw 0";

      plot[0].resize(2);

      legend.str("");
      if (nb_row == 1) {
        legend << STAT_label[STATL_REFERENCE] << " " << label << " " << row_identifier[0];
      }
      else if (nb_column == 1) {
        legend << STAT_label[STATL_TEST] << " " << label << " " << column_identifier[0];
      }
      else if (nb_row == nb_column) {
        legend << label;
      }
      plot[0][0].legend = legend.str();

      plot[0][0].style = "linespoints";

      plot[0][1].label = "true";

      for (i = 0;i < plot_nb_pattern;i++) {
        if (plot_distance[index[i]] == -D_INF) {
          plot_nb_pattern = i;
          break;
        }
        else {
          plot[0][0].add_point(i + 1 , plot_distance[index[i]]);
        }

        identifier.str("");
        if (nb_row > 1) {
          identifier << row_identifier[index[i]];
        }
        else {
          identifier << column_identifier[index[i]];
        }

        plot[0][1].add_text(i + 1 , plot_distance[index[i]] , identifier.str());
      }

      plot[0].xrange = Range(1 , plot_nb_pattern);
      plot[0].yrange = Range(plot_distance[index[0]] * (1. - PLOT_YMARGIN) ,
                             plot_distance[index[plot_nb_pattern - 1]] * (1. + PLOT_YMARGIN));

      if (plot_nb_pattern < TIC_THRESHOLD) {
        plot[0].xtics = 1;
      }
    }

    delete [] plot_distance;
    delete [] index;
  }

  return plot_set;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet DistanceMatrix.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* DistanceMatrix::get_plotable() const

{
  MultiPlotSet *plot_set;
  StatError error;


  return get_plotable(error);
}


/*--------------------------------------------------------------*
 *
 *  Test de la symetrie de la matrice des distances.
 *
 *--------------------------------------------------------------*/

bool DistanceMatrix::test_symmetry() const

{
  bool symmetry = false;


  if ((nb_row > 1) && (nb_row == nb_column)) {
    register int i , j;
    double normalized_distance[2];


    symmetry = true;
    for (i = 0;i < nb_row;i++) {
      for (j = i + 1;j < nb_row;j++) {
        normalized_distance[0] = distance[i][j];
        if ((distance[i][j] != -D_INF) && (length[i][j] > 1)) {
          normalized_distance[0] /= length[i][j];
        }
        normalized_distance[1] = distance[j][i];
        if ((distance[j][i] != -D_INF) && (length[j][i] > 1)) {
          normalized_distance[1] /= length[j][i];
        }

        if (normalized_distance[0] != normalized_distance[1]) {
          symmetry = false;
          break;
        }
      }

      if (!symmetry) {
        break;
      }
    }
  }

  return symmetry;
}


/*--------------------------------------------------------------*
 *
 *  Mise a jour d'un objet DistanceMatrix avec le resultat d'un alignement.
 *
 *  arguments : identificateurs des 2 formes, distance, longueur de l'alignement,
 *              distances et nombres d'elisions, d'insertions, nombre de matchs,
 *              distances et nombres de substitutions et de transpositions.
 *
 *--------------------------------------------------------------*/

void DistanceMatrix::update(int irow_identifier , int icolumn_identifier , double idistance ,
                            int alignment_length , double ideletion_distance , int inb_deletion ,
                            double iinsertion_distance , int inb_insertion , int inb_match ,
                            double isubstitution_distance , int inb_substitution ,
                            double itransposition_distance , int inb_transposition)

{
  if (idistance != -D_INF) {
    register int i , j;


    for (i = 0;i < nb_row;i++) {
      if (irow_identifier == row_identifier[i]) {
        for (j = 0;j < nb_column;j++) {
          if (icolumn_identifier == column_identifier[j]) {
            distance[i][j] = idistance;
            length[i][j] = alignment_length;

            deletion_distance[i][j] = ideletion_distance;
            nb_deletion[i][j] = inb_deletion;
            insertion_distance[i][j] = iinsertion_distance;
            nb_insertion[i][j] = inb_insertion;
            nb_match[i][j] = inb_match;

            if ((substitution_distance) && (nb_substitution)) {
              substitution_distance[i][j] = isubstitution_distance;
              nb_substitution[i][j] = inb_substitution;
            }

            if ((transposition_distance) && (nb_transposition)) {
              transposition_distance[i][j] = itransposition_distance;
              nb_transposition[i][j] = inb_transposition;
            }
            break;
          }
        }

        break;
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Mise a jour d'un objet DistanceMatrix avec le resultat d'une comparaison.
 *
 *  arguments : identificateurs des 2 formes, distance, longueur.
 *
 *--------------------------------------------------------------*/

void DistanceMatrix::update(int irow_identifier , int icolumn_identifier ,
                            double idistance , int ilength)

{
  if (idistance != -D_INF) {
    register int i , j;


    for (i = 0;i < nb_row;i++) {
      if (irow_identifier == row_identifier[i]) {
        for (j = 0;j < nb_column;j++) {
          if (icolumn_identifier == column_identifier[j]) {
            distance[i][j] = idistance;
            length[i][j] = ilength;
            break;
          }
        }

        break;
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'une longueur cumulee.
 *
 *  arguments : flags sur les lignes et sur les colonnes.
 *
 *--------------------------------------------------------------*/

int DistanceMatrix::cumul_length_computation(bool *row_flag , bool *column_flag) const

{
  register int i , j;
  int cumul_length = 0;


  for (i = 0;i < nb_row;i++) {
    if (row_flag[i]) {
      for (j = 0;j < nb_column;j++) {
        if (column_flag[j]) {
          cumul_length += length[i][j];
        }
      }
    }
  }

  return cumul_length;
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'une distance cumulee.
 *
 *  arguments : flags sur les lignes et sur les colonnes.
 *
 *--------------------------------------------------------------*/

double DistanceMatrix::cumul_distance_computation(bool *row_flag , bool *column_flag) const

{
  register int i , j;
  double cumul_distance = 0.;


  for (i = 0;i < nb_row;i++) {
    if (row_flag[i]) {
      for (j = 0;j < nb_column;j++) {
        if (column_flag[j]) {
          if (distance[i][j] != -D_INF) {
            if (distance[i][j] > 0.) {
              cumul_distance += distance[i][j];
            }
          }
          else {
            cumul_distance = -D_INF;
            break;
          }
        }
      }

      if (cumul_distance == -D_INF) {
        break;
      }
    }
  }

  return cumul_distance;
}


};  // namespace stat_tool
