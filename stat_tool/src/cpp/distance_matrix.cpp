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



#include <sstream>
#include <iomanip>

#include "tool/config.h"

// #include <rw/vstream.h>
// #include <rw/rwfile.h>
#include "stat_tools.h"
#include "distance_matrix.h"
#include "stat_label.h"

using namespace std;


extern int column_width(int value);
extern int column_width(int nb_value , const double *value , double scale = 1.);
extern char* label(const char *file_name);



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Distance_matrix.
 *
 *--------------------------------------------------------------*/

Distance_matrix::Distance_matrix()

{
  nb_row = 0;
  nb_column = 0;

  row_identifier = 0;
  column_identifier = 0;

  distance = 0;
  length = 0;

  deletion_distance = 0;
  nb_deletion = 0;
  insertion_distance = 0;
  nb_insertion = 0;
  nb_match = 0;
  substitution_distance = 0;
  nb_substitution = 0;
  transposition_distance = 0;
  nb_transposition = 0;

  label_size = 0;
  label = 0;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Distance_matrix.
 *
 *  arguments : nombre de formes, label, identificateurs des formes,
 *
 *--------------------------------------------------------------*/

Distance_matrix::Distance_matrix(int nb_pattern , const char *ilabel , int *pattern_identifier)

{
  register int i , j;
  int *plength;
  double *pdistance;


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

    pdistance = distance[i];
    plength = length[i];
    for (j = 0;j < nb_column;j++) {
      if (row_identifier[i] != column_identifier[j]) {
        *pdistance++ = -D_INF;
      }
      else {
        *pdistance++ = 0.;
      }

      *plength++ = 0;
    }
  }

  deletion_distance = 0;
  nb_deletion = 0;
  insertion_distance = 0;
  nb_insertion = 0;
  nb_match = 0;
  substitution_distance = 0;
  nb_substitution = 0;
  transposition_distance = 0;
  nb_transposition = 0;

  label_size = strlen(ilabel) + 1;
  label = new char[label_size];
  strcpy(label , ilabel);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Distance_matrix.
 *
 *  arguments : nombre de formes, identificateurs de ligne et de colonne, label,
 *              identificateurs des formes, flags operations de substitution et
 *              de transposition.
 *
 *--------------------------------------------------------------*/

Distance_matrix::Distance_matrix(int nb_pattern , int irow_identifier , int icolumn_identifier ,
                                 const char *ilabel , int *pattern_identifier ,
                                 bool substitution_flag , bool transposition_flag)

{
  register int i , j;
  int *plength;
  double *pdistance;


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

    pdistance = distance[i];
    plength = length[i];
    for (j = 0;j < nb_column;j++) {
      if (row_identifier[i] != column_identifier[j]) {
        *pdistance++ = -D_INF;
      }
      else {
        *pdistance++ = 0.;
      }

      *plength++ = 0;
    }
  }

  deletion_distance = new double*[nb_row];
  nb_deletion = new int*[nb_row];
  for (i = 0;i < nb_row;i++) {
    deletion_distance[i] = new double[nb_column];
    nb_deletion[i] = new int[nb_column];

    pdistance = deletion_distance[i];
    plength = nb_deletion[i];
    for (j = 0;j < nb_column;j++) {
      *pdistance++ = 0.;
      *plength++ = 0;
    }
  }

  insertion_distance = new double*[nb_row];
  nb_insertion = new int*[nb_row];
  for (i = 0;i < nb_row;i++) {
    insertion_distance[i] = new double[nb_column];
    nb_insertion[i] = new int[nb_column];

    pdistance = insertion_distance[i];
    plength = nb_insertion[i];
    for (j = 0;j < nb_column;j++) {
      *pdistance++ = 0.;
      *plength++ = 0;
    }
  }

  nb_match = new int*[nb_row];
  for (i = 0;i < nb_row;i++) {
    nb_match[i] = new int[nb_column];
    plength = nb_match[i];
    for (j = 0;j < nb_column;j++) {
      *plength++ = 0;
    }
  }

  if (substitution_flag) {
    substitution_distance = new double*[nb_row];
    nb_substitution = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      substitution_distance[i] = new double[nb_column];
      nb_substitution[i] = new int[nb_column];

      pdistance = substitution_distance[i];
      plength = nb_substitution[i];
      for (j = 0;j < nb_column;j++) {
        *pdistance++ = 0.;
        *plength++ = 0;
      }
    }
  }
  else {
    substitution_distance = 0;
    nb_substitution = 0;
  }

  if (transposition_flag) {
    transposition_distance = new double*[nb_row];
    nb_transposition = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      transposition_distance[i] = new double[nb_column];
      nb_transposition[i] = new int[nb_column];

      pdistance = transposition_distance[i];
      plength = nb_transposition[i];
      for (j = 0;j < nb_column;j++) {
        *pdistance++ = 0.;
        *plength++ = 0;
      }
    }
  }
  else {
    transposition_distance = 0;
    nb_transposition = 0;
  }

  label_size = strlen(ilabel) + 1;
  label = new char[label_size];
  strcpy(label , ilabel);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Distance_matrix.
 *
 *  arguments : reference sur un objet Distance_matrix,
 *              nombre et identificateurs des formes selectionnees,
 *              flag pour conserver ou rejeter les formes selectionnees.
 *
 *--------------------------------------------------------------*/

Distance_matrix::Distance_matrix(const Distance_matrix &dist_matrix , int inb_pattern ,
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
    deletion_distance = 0;
    nb_deletion = 0;
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
    insertion_distance = 0;
    nb_insertion = 0;
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
    nb_match = 0;
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
    substitution_distance = 0;
    nb_substitution = 0;
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
    transposition_distance = 0;
    nb_transposition = 0;
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
 *  Constructeur de la classe Distance_matrix.
 *
 *  arguments : reference sur un objet Distance_matrix, nombre de groupes, label.
 *
 *--------------------------------------------------------------*/

Distance_matrix::Distance_matrix(const Distance_matrix &dist_matrix , int nb_cluster ,
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
    deletion_distance = 0;
    nb_deletion = 0;
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
    insertion_distance = 0;
    nb_insertion = 0;
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
    nb_match = 0;
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
    substitution_distance = 0;
    nb_substitution = 0;
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
    transposition_distance = 0;
    nb_transposition = 0;
  }

  label_size = strlen(ilabel) + 1;
  label = new char[label_size];
  strcpy(label , ilabel);
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Distance_matrix.
 *
 *  arguments : reference sur un objet Distance_matrix,
 *              transformation de la matrice ('s' :  symetrisation, 'u' : denormalisation).
 *
 *--------------------------------------------------------------*/

void Distance_matrix::copy(const Distance_matrix &dist_matrix , char transform)

{
  register int i , j;
  int *plength , *clength;
  double *pdistance , *cdistance;


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
      pdistance = distance[i];
      cdistance = dist_matrix.distance[i];
      plength = length[i];
      clength = dist_matrix.length[i];
      for (j = 0;j < nb_column;j++) {
        if ((*cdistance != -D_INF) && (*clength > 0)) {
          *pdistance++ = *cdistance++ / *clength++;
          *plength++ = 1;
        }
        else {
          *pdistance++ = *cdistance++;
          *plength++ = *clength++;
        }
      }
      break;
    }

    default : {
      pdistance = distance[i];
      cdistance = dist_matrix.distance[i];
      plength = length[i];
      clength = dist_matrix.length[i];
      for (j = 0;j < nb_column;j++) {
        *pdistance++ = *cdistance++;
        *plength++ = *clength++;
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
          pdistance = deletion_distance[i];
          cdistance = dist_matrix.deletion_distance[i];
          plength = nb_deletion[i];
          clength = dist_matrix.nb_deletion[i];
          for (j = 0;j < nb_column;j++) {
            *pdistance++ = *cdistance++;
            *plength++ = *clength++;
          }
          break;
        }
        }
      }
    }

    else {
      deletion_distance = 0;
      nb_deletion = 0;
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
            nb_insertion[j][i]= nb_insertion[i][j];
          }
          break;
        }

        default : {
          pdistance = insertion_distance[i];
          cdistance = dist_matrix.insertion_distance[i];
          plength = nb_insertion[i];
          clength = dist_matrix.nb_insertion[i];
          for (j = 0;j < nb_column;j++) {
            *pdistance++ = *cdistance++;
            *plength++ = *clength++;
          }
          break;
        }
        }
      }
    }

    else {
      insertion_distance = 0;
      nb_insertion = 0;
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
          plength = nb_match[i];
          clength = dist_matrix.nb_match[i];
          for (j = 0;j < nb_column;j++) {
            *plength++ = *clength++;
          }
          break;
        }
        }
      }
    }

    else {
      nb_match = 0;
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
          pdistance = substitution_distance[i];
          cdistance = dist_matrix.substitution_distance[i];
          plength = nb_substitution[i];
          clength = dist_matrix.nb_substitution[i];
          for (j = 0;j < nb_column;j++) {
            *pdistance++ = *cdistance++;
            *plength++ = *clength++;
          }
          break;
        }
        }
      }
    }

    else {
      substitution_distance = 0;
      nb_substitution = 0;
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
          pdistance = transposition_distance[i];
          cdistance = dist_matrix.transposition_distance[i];
          plength = nb_transposition[i];
          clength = dist_matrix.nb_transposition[i];
          for (j = 0;j < nb_column;j++) {
            *pdistance++ = *cdistance++;
            *plength++ = *clength++;
          }
          break;
        }
        }
      }
    }

    else {
      transposition_distance = 0;
      nb_transposition = 0;
    }
  }

  else {
    deletion_distance = 0;
    nb_deletion = 0;
    insertion_distance = 0;
    nb_insertion = 0;
    nb_match = 0;
    substitution_distance = 0;
    nb_substitution = 0;
    transposition_distance = 0;
    nb_transposition = 0;
  }

  label_size = dist_matrix.label_size;
  label = new char[label_size];
  strcpy(label , dist_matrix.label);
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Distance_matrix.
 *
 *--------------------------------------------------------------*/

void Distance_matrix::remove()

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
 *  Destructeur de la classe Distance_matrix.
 *
 *--------------------------------------------------------------*/

Distance_matrix::~Distance_matrix()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Distance_matrix.
 *
 *  argument : reference sur un objet Distance_matrix.
 *
 *--------------------------------------------------------------*/

Distance_matrix& Distance_matrix::operator=(const Distance_matrix &dist_matrix)

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
 *  arguments : reference sur un objet Format_error, nombre de formes,
 *              identificateurs des formes, flag pour conserver ou rejeter
 *              les formes selectionnees.
 *
 *--------------------------------------------------------------*/

Distance_matrix* Distance_matrix::select_individual(Format_error &error , int inb_pattern ,
                                                    int *iidentifier , bool keep) const

{
  bool status = true , *selected_pattern;
  register int i , j;
  int nb_pattern , max_identifier , *identifier;
  Distance_matrix *dist_matrix;


  dist_matrix = 0;
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
                      << label << " " << STAT_error[STATR_IDENTIFIER];
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
    dist_matrix = new Distance_matrix(*this , inb_pattern , iidentifier , keep);
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*
 *
 *  Symetrisation d'une matrice des distances.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Distance_matrix* Distance_matrix::symmetrize(Format_error &error) const

{
  Distance_matrix *dist_matrix;


  error.init();

  if (test_symmetry()) {
    dist_matrix = 0;
    error.update(STAT_error[STATR_SYMMETRICAL_MATRIX]);
  }
  else {
    dist_matrix = new Distance_matrix(*this , 's');
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*
 *
 *  Denormalisation d'une matrice des distances.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Distance_matrix* Distance_matrix::unnormalize(Format_error &error) const

{
  bool status = false;
  register int i , j;
  Distance_matrix *dist_matrix;


  error.init();
  dist_matrix = 0;

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
    dist_matrix = new Distance_matrix(*this , 'u');
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet Distance_matrix.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Distance_matrix::line_write(ostream &os) const

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

ostream& Distance_matrix::property_print(double **normalized_distance ,
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
  int cumul = 0 , *pvalue;


  for (i = 0;i < nb_row;i++) {
    pvalue = value[i];
    for (j = 0;j < nb_column;j++) {
      cumul += *pvalue++;
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
    if (*distance != -D_INF) {
      cumul_distance += *distance++;
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
 *  Ecriture d'un objet Distance_matrix.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Distance_matrix::ascii_write(ostream &os , bool exhaustive) const

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
      cumul_distance[i] = ::cumul_distance_computation(nb_column , distance[i]);

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
          cumul_deletion_distance = ::cumul_distance_computation(nb_column , deletion_distance[index[i]]);
          cumul_insertion_distance = ::cumul_distance_computation(nb_column , insertion_distance[index[i]]);
          os << " = " << cumul_deletion_distance / cumul_length[index[i]]
             << " (" << cumul_nb_deletion << " d) + "
             << cumul_insertion_distance / cumul_length[index[i]]
             << " (" << cumul_nb_insertion << " i) + 0 (" << cumul_nb_match << " m)";

          if (substitution_distance) {
            cumul_substitution_distance = ::cumul_distance_computation(nb_column , substitution_distance[index[i]]);
            os << " + " << cumul_substitution_distance / cumul_length[index[i]]
               << " (" << cumul_nb_substitution << " s)";
          }

          if (transposition_distance) {
            cumul_transposition_distance = ::cumul_distance_computation(nb_column , transposition_distance[index[i]]);
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

  if ((exhaustive) && (nb_row > 1) && (nb_row == nb_column) && (nb_row <= ASCII_NB_PATTERN) &&
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
 *  Ecriture d'un objet Distance_matrix dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Distance_matrix::ascii_write(Format_error &error , const char *path ,
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
 *  Ecriture d'un objet Distance_matrix au format tableur.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Distance_matrix::spreadsheet_write(ostream &os) const

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
      cumul_distance[i] = ::cumul_distance_computation(nb_column , distance[i]);

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
          cumul_deletion_distance = ::cumul_distance_computation(nb_column , deletion_distance[index[i]]);
          cumul_insertion_distance = ::cumul_distance_computation(nb_column , insertion_distance[index[i]]);
          os << "\t" << cumul_deletion_distance / cumul_length[index[i]]
             << " (" << cumul_nb_deletion << " d)\t"
             << cumul_insertion_distance / cumul_length[index[i]]
             << " (" << cumul_nb_insertion << " i)\t0 (" << cumul_nb_match << " m)";

          if (substitution_distance) {
            cumul_substitution_distance = ::cumul_distance_computation(nb_column , substitution_distance[index[i]]);
            os << "\t" << cumul_substitution_distance / cumul_length[index[i]]
               << " (" << cumul_nb_substitution << " s)";
          }

          if (transposition_distance) {
            cumul_transposition_distance = ::cumul_distance_computation(nb_column , transposition_distance[index[i]]);
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
 *  Ecriture d'un objet Distance_matrix dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Distance_matrix::spreadsheet_write(Format_error &error , const char *path) const

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
 *  Sortie Gnuplot d'un objet Distance_matrix.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Distance_matrix::plot_write(Format_error &error , const char *prefix ,
                                 const char *title) const

{
  bool status;
  register int i , j;
  int plot_nb_pattern , *index;
  double max , *plot_distance;
  ostringstream data_file_name;


  // ecriture du fichier de donnees

  data_file_name << prefix << ".dat";
  ofstream out_data_file((data_file_name.str()).c_str());

  if (!out_data_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  else {
    status = true;

    if ((nb_row == 1) && (nb_column == 1)) {
      plot_nb_pattern = 2;
      plot_distance = new double[1];

      plot_distance[0] = distance[0][0];
      if ((distance[0][0] != -D_INF) && (length[0][0] > 1)) {
        plot_distance[0] /= length[0][0];
      }

      out_data_file << 1 << " " << 0 << endl;
      out_data_file << 2 << " " << plot_distance[0] << endl;
    }

    else {
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
          plot_distance[i] = ::cumul_distance_computation(nb_column , distance[i]);
          if (plot_distance[i] != -D_INF) {
            plot_distance[i] /= cumul_computation(1 , nb_column , length + i);
          }
        }
      }

      index = pattern_sort(plot_nb_pattern , plot_distance , I_DEFAULT);

      for (i = 0;i < plot_nb_pattern;i++) {
        if (plot_distance[index[i]] == -D_INF) {
          plot_nb_pattern = i;
          break;
        }
        else {
          out_data_file << i + 1 << " " << plot_distance[index[i]] << endl;
        }
      }
    }

    // ecriture du fichier de commandes et du fichier d'impression

    if (plot_nb_pattern > 0) {
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
          file_name[1] << ::label(prefix) << ".ps";
          out_file << "set output \"" << file_name[1].str() << "\"\n\n";
        }

        out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                 << "set title";
        if (title) {
          out_file << " \"" << title << "\"";
        }
        out_file << "\n\n";

        if ((nb_row == 1) && (nb_column == 1)) {
          out_file << "set label \"" << row_identifier[0] << "\" at 1,0" << endl;
          out_file << "set label \"" << column_identifier[0] << "\" at 2,"
                   << plot_distance[0] << endl;

          max = plot_distance[0];
        }

        else {
          for (j = 0;j < plot_nb_pattern;j++) {
            if (nb_row > 1) {
              out_file << "set label \"" << row_identifier[index[j]] << "\" at " << j + 1 << ","
                       << plot_distance[index[j]] << endl;
            }
            else {
              out_file << "set label \"" << column_identifier[index[j]] << "\" at " << j + 1 << ","
                       << plot_distance[index[j]] << endl;
            }
          }

          max = plot_distance[index[plot_nb_pattern - 1]];
        }
        out_file << endl;

        if (plot_nb_pattern < TIC_THRESHOLD) {
          out_file << "set xtics 1,1" << endl;
        }

        out_file << "plot [1:" << plot_nb_pattern << "] [0:" << max * YSCALE
                 << "] \"" << ::label((data_file_name.str()).c_str()) << "\" title \"";
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

    delete [] plot_distance;
    if ((nb_row > 1) || (nb_column > 1)) {
      delete [] index;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWDEFINE_COLLECTABLE(Distance_matrix , STATI_DISTANCE_MATRIX);


RWspace Distance_matrix::binaryStoreSize() const

{
  RWspace size;


  size = sizeof(nb_row) + sizeof(nb_column) + sizeof(*row_identifier) * nb_row +
         sizeof(*column_identifier) * nb_column + sizeof(**distance) * nb_row * nb_column +
         sizeof(**length) * nb_row * nb_column;

  size += sizeof(true);
  if ((deletion_distance) && (nb_deletion)) {
    size += sizeof(**deletion_distance) * nb_row * nb_column + sizeof(**nb_deletion) * nb_row * nb_column;
  }

  size += sizeof(true);
  if ((insertion_distance) && (nb_insertion)) {
    size += sizeof(**insertion_distance) * nb_row * nb_column + sizeof(**nb_insertion) * nb_row * nb_column;
  }

  size += sizeof(true);
  if (nb_match) {
    size += sizeof(**nb_match) * nb_row * nb_column;
  }

  size += sizeof(true);
  if ((substitution_distance) && (nb_substitution)) {
    size += sizeof(**substitution_distance) * nb_row * nb_column + sizeof(**nb_substitution) * nb_row * nb_column;
  }

  size += sizeof(true);
  if ((transposition_distance) && (nb_transposition)) {
    size += sizeof(**transposition_distance) * nb_row * nb_column + sizeof(**nb_transposition) * nb_row * nb_column;
  }

  size += sizeof(label_size) + sizeof(*label) * label_size;

  return size;
}


void Distance_matrix::restoreGuts(RWvistream &is)

{
  bool status;
  register int i , j;


  remove();

  is >> nb_row >> nb_column;

  row_identifier = new int[nb_row];
  for (i = 0;i < nb_row;i++) {
    is >> row_identifier[i];
  }
  column_identifier = new int[nb_column];
  for (i = 0;i < nb_column;i++) {
    is >> column_identifier[i];
  }

  distance = new double*[nb_row];
  for (i = 0;i < nb_row;i++) {
    distance[i] = new double[nb_column];
    for (j = 0;j < nb_column;j++) {
      is >> distance[i][j];
    }
  }
  length = new int*[nb_row];
  for (i = 0;i < nb_row;i++) {
    length[i] = new int[nb_column];
    for (j = 0;j < nb_column;j++) {
      is >> length[i][j];
    }
  }

  is >> status;
  if (status) {
    deletion_distance = new double*[nb_row];
    for (i = 0;i < nb_row;i++) {
      deletion_distance[i] = new double[nb_column];
      for (j = 0;j < nb_column;j++) {
        is >> deletion_distance[i][j];
      }
    }
    nb_deletion = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      nb_deletion[i] = new int[nb_column];
      for (j = 0;j < nb_column;j++) {
        is >> nb_deletion[i][j];
      }
    }
  }
  else {
    deletion_distance = 0;
    nb_deletion = 0;
  }

  is >> status;
  if (status) {
    insertion_distance = new double*[nb_row];
    for (i = 0;i < nb_row;i++) {
      insertion_distance[i] = new double[nb_column];
      for (j = 0;j < nb_column;j++) {
        is >> insertion_distance[i][j];
      }
    }
    nb_insertion = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      nb_insertion[i] = new int[nb_column];
      for (j = 0;j < nb_column;j++) {
        is >> nb_insertion[i][j];
      }
    }
  }
  else {
    insertion_distance = 0;
    nb_insertion = 0;
  }

  is >> status;
  if (status) {
    nb_match = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      nb_match[i] = new int[nb_column];
      for (j = 0;j < nb_column;j++) {
        is >> nb_match[i][j];
      }
    }
  }
  else {
    nb_match = 0;
  }

  is >> status;
  if (status) {
    substitution_distance = new double*[nb_row];
    for (i = 0;i < nb_row;i++) {
      substitution_distance[i] = new double[nb_column];
      for (j = 0;j < nb_column;j++) {
        is >> substitution_distance[i][j];
      }
    }
    nb_substitution = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      nb_substitution[i] = new int[nb_column];
      for (j = 0;j < nb_column;j++) {
        is >> nb_substitution[i][j];
      }
    }
  }
  else {
    substitution_distance = 0;
    nb_substitution = 0;
  }

  is >> status;
  if (status) {
    transposition_distance = new double*[nb_row];
    for (i = 0;i < nb_row;i++) {
      transposition_distance[i] = new double[nb_column];
      for (j = 0;j < nb_column;j++) {
        is >> transposition_distance[i][j];
      }
    }
    nb_transposition = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      nb_transposition[i] = new int[nb_column];
      for (j = 0;j < nb_column;j++) {
        is >> nb_transposition[i][j];
      }
    }
  }
  else {
    transposition_distance = 0;
    nb_transposition = 0;
  }

  is >> label_size;
  label = new char[label_size];
  is.getString(label , label_size);
}


void Distance_matrix::restoreGuts(RWFile &file)

{
  bool status;
  register int i;


  remove();

  file.Read(nb_row);
  file.Read(nb_column);

  row_identifier = new int[nb_row];
  file.Read(row_identifier , nb_row);
  column_identifier = new int[nb_column];
  file.Read(column_identifier , nb_column);

  distance = new double*[nb_row];
  for (i = 0;i < nb_row;i++) {
    distance[i] = new double[nb_column];
    file.Read(distance[i] , nb_column);
  }
  length = new int*[nb_row];
  for (i = 0;i < nb_row;i++) {
    length[i] = new int[nb_column];
    file.Read(length[i] , nb_column);
  }

  file.Read(status);
  if (status) {
    deletion_distance = new double*[nb_row];
    for (i = 0;i < nb_row;i++) {
      deletion_distance[i] = new double[nb_column];
      file.Read(deletion_distance[i] , nb_column);
    }
    nb_deletion = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      nb_deletion[i] = new int[nb_column];
      file.Read(nb_deletion[i] , nb_column);
    }
  }
  else {
    deletion_distance = 0;
    nb_deletion = 0;
  }

  file.Read(status);
  if (status) {
    insertion_distance = new double*[nb_row];
    for (i = 0;i < nb_row;i++) {
      insertion_distance[i] = new double[nb_column];
      file.Read(insertion_distance[i] , nb_column);
    }
    nb_insertion = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      nb_insertion[i] = new int[nb_column];
      file.Read(nb_insertion[i] , nb_column);
    }
  }
  else {
    insertion_distance = 0;
    nb_insertion = 0;
  }

  file.Read(status);
  if (status) {
    nb_match = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      nb_match[i] = new int[nb_column];
      file.Read(nb_match[i] , nb_column);
    }
  }
  else {
    nb_match = 0;
  }

  file.Read(status);
  if (status) {
    substitution_distance = new double*[nb_row];
    for (i = 0;i < nb_row;i++) {
      substitution_distance[i] = new double[nb_column];
      file.Read(substitution_distance[i] , nb_column);
    }
    nb_substitution = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      nb_substitution[i] = new int[nb_column];
      file.Read(nb_substitution[i] , nb_column);
    }
  }
  else {
    substitution_distance = 0;
    nb_substitution = 0;
  }

  file.Read(status);
  if (status) {
    transposition_distance = new double*[nb_row];
    for (i = 0;i < nb_row;i++) {
      transposition_distance[i] = new double[nb_column];
      file.Read(transposition_distance[i] , nb_column);
    }
    nb_transposition = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      nb_transposition[i] = new int[nb_column];
      file.Read(nb_transposition[i] , nb_column);
    }
  }
  else {
    transposition_distance = 0;
    nb_transposition = 0;
  }

  file.Read(label_size);
  label = new char[label_size];
  file.Read(label , label_size);
}


void Distance_matrix::saveGuts(RWvostream &os) const

{
  register int i , j;


  os << nb_row << nb_column;

  for (i = 0;i < nb_row;i++) {
    os << row_identifier[i];
  }
  for (i = 0;i < nb_column;i++) {
    os << column_identifier[i];
  }

  for (i = 0;i < nb_row;i++) {
    for (j = 0;j < nb_column;j++) {
      os << distance[i][j];
    }
  }
  for (i = 0;i < nb_row;i++) {
    for (j = 0;j < nb_column;j++) {
      os << length[i][j];
    }
  }

  if ((deletion_distance) && (nb_deletion)) {
    os << true;
    for (i = 0;i < nb_row;i++) {
      for (j = 0;j < nb_column;j++) {
        os << deletion_distance[i][j];
      }
    }
    for (i = 0;i < nb_row;i++) {
      for (j = 0;j < nb_column;j++) {
        os << nb_deletion[i][j];
      }
    }
  }
  else {
    os << false;
  }

  if ((insertion_distance) && (nb_insertion)) {
    os << true;
    for (i = 0;i < nb_row;i++) {
      for (j = 0;j < nb_column;j++) {
        os << insertion_distance[i][j];
      }
    }
    for (i = 0;i < nb_row;i++) {
      for (j = 0;j < nb_column;j++) {
        os << nb_insertion[i][j];
      }
    }
  }
  else {
    os << false;
  }

  if (nb_match) {
    os << true;
    for (i = 0;i < nb_row;i++) {
      for (j = 0;j < nb_column;j++) {
        os << nb_match[i][j];
      }
    }
  }
  else {
    os << false;
  }

  if ((substitution_distance) && (nb_substitution)) {
    os << true;
    for (i = 0;i < nb_row;i++) {
      for (j = 0;j < nb_column;j++) {
        os << substitution_distance[i][j];
      }
    }
    for (i = 0;i < nb_row;i++) {
      for (j = 0;j < nb_column;j++) {
        os << nb_substitution[i][j];
      }
    }
  }
  else {
    os << false;
  }

  if ((transposition_distance) && (nb_transposition)) {
    os << true;
    for (i = 0;i < nb_row;i++) {
      for (j = 0;j < nb_column;j++) {
        os << transposition_distance[i][j];
      }
    }
    for (i = 0;i < nb_row;i++) {
      for (j = 0;j < nb_column;j++) {
        os << nb_transposition[i][j];
      }
    }
  }
  else {
    os << false;
  }

  os << label_size;
  os << label;
}


void Distance_matrix::saveGuts(RWFile &file) const

{
  register int i;


  file.Write(nb_row);
  file.Write(nb_column);

  file.Write(row_identifier , nb_row);
  file.Write(column_identifier , nb_column);

  for (i = 0;i < nb_row;i++) {
    file.Write(distance[i] , nb_column);
  }
  for (i = 0;i < nb_row;i++) {
    file.Write(length[i] , nb_column);
  }

  if ((deletion_distance) && (nb_deletion)) {
    file.Write(true);
    for (i = 0;i < nb_row;i++) {
      file.Write(deletion_distance[i] , nb_column);
    }
    for (i = 0;i < nb_row;i++) {
      file.Write(nb_deletion[i] , nb_column);
    }
  }
  else {
    file.Write(false);
  }

  if ((insertion_distance) && (nb_insertion)) {
    file.Write(true);
    for (i = 0;i < nb_row;i++) {
      file.Write(insertion_distance[i] , nb_column);
    }
    for (i = 0;i < nb_row;i++) {
      file.Write(nb_insertion[i] , nb_column);
    }
  }
  else {
    file.Write(false);
  }

  if (nb_match) {
    file.Write(true);
    for (i = 0;i < nb_row;i++) {
      file.Write(nb_match[i] , nb_column);
    }
  }
  else {
    file.Write(false);
  }

  if ((substitution_distance) && (nb_substitution)) {
    file.Write(true);
    for (i = 0;i < nb_row;i++) {
      file.Write(substitution_distance[i] , nb_column);
    }
    for (i = 0;i < nb_row;i++) {
      file.Write(nb_substitution[i] , nb_column);
    }
  }
  else {
    file.Write(false);
  }

  if ((transposition_distance) && (nb_transposition)) {
    file.Write(true);
    for (i = 0;i < nb_row;i++) {
      file.Write(transposition_distance[i] , nb_column);
    }
    for (i = 0;i < nb_row;i++) {
      file.Write(nb_transposition[i] , nb_column);
    }
  }
  else {
    file.Write(false);
  }

  file.Write(label_size);
  file.Write(label , label_size);
} */


/*--------------------------------------------------------------*
 *
 *  Test de la symetrie de la matrice des distances.
 *
 *--------------------------------------------------------------*/

bool Distance_matrix::test_symmetry() const

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
 *  Mise a jour d'un objet Distance_matrix avec le resultat d'un alignement.
 *
 *  arguments : identificateurs des 2 formes, distance, longueur de l'alignement,
 *              distances et nombres d'elisions, d'insertions, nombre de matchs,
 *              distances et nombres de substitutions et de transpositions.
 *
 *--------------------------------------------------------------*/

void Distance_matrix::update(int irow_identifier , int icolumn_identifier , double idistance ,
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
 *  Mise a jour d'un objet Distance_matrix avec le resultat d'une comparaison.
 *
 *  arguments : identificateurs des 2 formes, distance, longueur.
 *
 *--------------------------------------------------------------*/

void Distance_matrix::update(int irow_identifier , int icolumn_identifier , double idistance ,
                             int ilength)

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

int Distance_matrix::cumul_length_computation(bool *row_flag , bool *column_flag) const

{
  bool *pcolumn_flag;
  register int i , j;
  int cumul_length = 0 , *plength;


  for (i = 0;i < nb_row;i++) {
    if (*row_flag++) {
      pcolumn_flag = column_flag;
      plength = length[i];

      for (j = 0;j < nb_column;j++) {
        if (*pcolumn_flag++) {
          cumul_length += *plength;
        }
        plength++;
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

double Distance_matrix::cumul_distance_computation(bool *row_flag , bool *column_flag) const

{
  bool *pcolumn_flag;
  register int i , j;
  double cumul_distance = 0. , *pdistance;


  for (i = 0;i < nb_row;i++) {
    if (*row_flag++) {
      pcolumn_flag = column_flag;
      pdistance = distance[i];

      for (j = 0;j < nb_column;j++) {
        if (*pcolumn_flag++) {
          if (*pdistance != -D_INF) {
            if (*pdistance > 0.) {
              cumul_distance += *pdistance;
            }
          }
          else {
            cumul_distance = -D_INF;
            break;
          }
        }

        pdistance++;
      }

      if (cumul_distance == -D_INF) {
        break;
      }
    }
  }

  return cumul_distance;
}
