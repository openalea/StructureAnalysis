/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2017 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
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



#include <string>
#include <iomanip>

#include "distance_matrix.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {


extern int cumul_computation(int nb_row , int nb_column , int **value);
extern double cumul_distance_computation(int dim , double *distance);
extern int* pattern_sort(int nb_pattern , double *distance , int nb_sorted_pattern = I_DEFAULT);



/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the Dendrogram class.
 */
/*--------------------------------------------------------------*/

Dendrogram::Dendrogram()

{
  distance_matrix = NULL;

  scale = CHILD_CLUSTER_DISTANCE;

  nb_cluster = 0;
  cluster_nb_pattern = NULL;
  cluster_pattern = NULL;

  parent = NULL;
  child = NULL;

  child_distance = NULL;
  within_cluster_distance = NULL;
  between_cluster_distance = NULL;
  max_within_cluster_distance = NULL;
  min_between_cluster_distance = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Dendrogram class.
 *
 *  \param[in] dist_matrix reference on a DistanceMatrix object,
 *  \param[in] iscale      cluster distance scale.
 */
/*--------------------------------------------------------------*/

Dendrogram::Dendrogram(const DistanceMatrix &dist_matrix , cluster_scale iscale)

{
  int i;


  distance_matrix = new DistanceMatrix(dist_matrix);

  scale = iscale;

  nb_cluster = 2 * distance_matrix->nb_row - 1;

  cluster_nb_pattern = new int[nb_cluster];
  cluster_pattern = new int*[nb_cluster];

  for (i = 0;i < distance_matrix->nb_row;i++) {
    cluster_nb_pattern[i] = 1;
    cluster_pattern[i] = new int[1];
    cluster_pattern[i][0] = i;
  }
  for (i = distance_matrix->nb_row;i < 2 * distance_matrix->nb_row - 1;i++) {
    cluster_nb_pattern[i] = 0;
    cluster_pattern[i] = NULL;
  }

  parent = new int[nb_cluster];
  child = new int*[nb_cluster];
  for (i = 0;i < distance_matrix->nb_row;i++) {
    child[i] = NULL;
  }
  for (i = distance_matrix->nb_row;i < 2 * distance_matrix->nb_row - 1;i++) {
    child[i] = new int[2];
  }
  parent[2 * distance_matrix->nb_row - 2] = I_DEFAULT;

  child_distance = new double[nb_cluster];
  for (i = 0;i < distance_matrix->nb_row;i++) {
    child_distance[i] = D_DEFAULT;
  }

  within_cluster_distance = new double[nb_cluster];
  for (i = 0;i < distance_matrix->nb_row;i++) {
    within_cluster_distance[i] = 0.;
  }

  between_cluster_distance = new double[nb_cluster - 1];

  max_within_cluster_distance = new double[nb_cluster];
  for (i = 0;i < distance_matrix->nb_row;i++) {
    max_within_cluster_distance[i] = 0.;
  }

  min_between_cluster_distance = new double[nb_cluster - 1];
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a Dendrogram object.
 *
 *  \param[in] dendrogram reference on a Dendrogram object.
 */
/*--------------------------------------------------------------*/

void Dendrogram::copy(const Dendrogram &dendrogram)

{
  int i , j;


  distance_matrix = new DistanceMatrix(*(dendrogram.distance_matrix));

  scale = dendrogram.scale;

  nb_cluster = dendrogram.nb_cluster;

  cluster_nb_pattern = new int[nb_cluster];
  for (i = 0;i < nb_cluster;i++) {
    cluster_nb_pattern[i] = dendrogram.cluster_nb_pattern[i];
  }

  cluster_pattern = new int*[nb_cluster];
  for (i = 0;i < nb_cluster;i++) {
    cluster_pattern[i] = new int[cluster_nb_pattern[i]];
    for (j = 0;j < cluster_nb_pattern[i];j++) {
      cluster_pattern[i][j] = dendrogram.cluster_pattern[i][j];
    }
  }

  parent = new int[nb_cluster];
  for (i = 0;i < nb_cluster;i++) {
    parent[i] = dendrogram.parent[i];
  }

  child = new int*[nb_cluster];
  for (i = 0;i < distance_matrix->nb_row;i++) {
    child[i] = NULL;
  }
  for (i = distance_matrix->nb_row;i < 2 * distance_matrix->nb_row - 1;i++) {
    child[i] = new int[2];
    for (j = 0;j < 2;j++) {
      child[i][j] = dendrogram.child[i][j];
    }
  }

  child_distance = new double[nb_cluster];
  for (i = 0;i < nb_cluster;i++) {
    child_distance[i] = dendrogram.child_distance[i];
  }

  within_cluster_distance = new double[nb_cluster];
  for (i = 0;i < nb_cluster;i++) {
    within_cluster_distance[i] = dendrogram.within_cluster_distance[i];
  }

  between_cluster_distance = new double[nb_cluster - 1];
  for (i = 0;i < nb_cluster - 1;i++) {
    between_cluster_distance[i] = dendrogram.between_cluster_distance[i];
  }

  max_within_cluster_distance = new double[nb_cluster];
  for (i = 0;i < nb_cluster;i++) {
    max_within_cluster_distance[i] = dendrogram.max_within_cluster_distance[i];
  }

  min_between_cluster_distance = new double[nb_cluster - 1];
  for (i = 0;i < nb_cluster - 1;i++) {
    min_between_cluster_distance[i] = dendrogram.min_between_cluster_distance[i];
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the data members of a Dendrogram object.
 */
/*--------------------------------------------------------------*/

void Dendrogram::remove()

{
  int i;


  delete [] cluster_nb_pattern;

  if (cluster_pattern) {
    for (i = 0;i < nb_cluster;i++) {
      delete [] cluster_pattern[i];
    }
    delete [] cluster_pattern;
  }

  delete [] parent;

  if (child) {
    for (i = distance_matrix->nb_row;i < 2 * distance_matrix->nb_row - 1;i++) {
      delete [] child[i];
    }
    delete [] child;
  }

  delete distance_matrix;

  delete [] child_distance;
  delete [] within_cluster_distance;
  delete [] between_cluster_distance;
  delete [] max_within_cluster_distance;
  delete [] min_between_cluster_distance;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the Dendrogram class.
 */
/*--------------------------------------------------------------*/

Dendrogram::~Dendrogram()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the Dendrogram class.
 *
 *  \param[in] dendrogram reference on a Dendrogram object.
 *
 *  \return               Dendrogram object.
 */
/*--------------------------------------------------------------*/

Dendrogram& Dendrogram::operator=(const Dendrogram &dendrogram)

{
  if (&dendrogram != this) {
    remove();
    copy(dendrogram);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing on a single line of a Dendrogram object.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& Dendrogram::line_write(ostream &os) const

{
  os << distance_matrix->nb_row << " " << distance_matrix->label << "s";

  os << "   " << STAT_label[STATL_CHILD_CLUSTER_DISTANCE_COEFF] << ": "
     << coefficient_computation(CHILD_CLUSTER_DISTANCE);

  os << "   " << STAT_label[STATL_DIAMETER_COEFF] << ": "
     << coefficient_computation(DIAMETER);

/*  switch (scale) {
  case CHILD_CLUSTER_DISTANCE :
    os << "   " << STAT_label[STATL_CHILD_CLUSTER_DISTANCE_COEFF];
    break;
  case DIAMETER :
    os << "   " << STAT_label[STATL_DIAMETER_COEFF];
    break;
  }
  os << ": " << coefficient_computation(); */

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Sort of the between-cluster distances.
 *
 *  \return ordered distances.
 */
/*--------------------------------------------------------------*/

double* Dendrogram::distance_ordering() const

{
  int i , j;
  int offset;
  double *ordered_distance , *distance;


  ordered_distance = new double[distance_matrix->nb_row - 1];

  switch (scale) {
  case CHILD_CLUSTER_DISTANCE :
    distance = child_distance;
    break;
  case DIAMETER :
    distance = max_within_cluster_distance;
    break;
  }

  for (i = 2 * distance_matrix->nb_row - 2;i >= distance_matrix->nb_row;i--) {
    if (i == 2 * distance_matrix->nb_row - 2) {
      offset = 0;
    }
    else {
      for (j = 0;j < distance_matrix->nb_row - 1;j++) {
        if (cluster_pattern[i][0] == cluster_pattern[2 * distance_matrix->nb_row - 2][j]) {
          offset = j;
          break;
        }
      }
    }

    for (j = i - 1;j >= distance_matrix->nb_row;j--) {
      if (cluster_pattern[i][0] == cluster_pattern[j][0]) {
        ordered_distance[offset + cluster_nb_pattern[j] - 1] = distance[i];
        break;
      }
    }

    if (j < distance_matrix->nb_row) {
      ordered_distance[offset] = distance[i];
    }
  }

  return ordered_distance;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Dendrogram object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive detail level.
 */
/*--------------------------------------------------------------*/

ostream& Dendrogram::ascii_write(ostream &os , bool exhaustive) const

{
  int i , j , k;
  int buff , max_identifier , max_nb_character , previous_nb_character ,
      *nb_character , width[3];
  double min_distance , min_diff_distance , *ordered_distance , *distance;
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::right , ios::adjustfield);

  os << distance_matrix->nb_row << " " << distance_matrix->label << "s" << endl;

  // computation of the column widths

  width[0] = column_width(distance_matrix->nb_row) + 1;
  width[1] = column_width(distance_matrix->nb_row - 1 , child_distance + distance_matrix->nb_row);
  if (scale == CHILD_CLUSTER_DISTANCE) {
    width[2] = width[1];
  }
  buff = column_width(distance_matrix->nb_row - 1 , within_cluster_distance + distance_matrix->nb_row);
  if (buff > width[1]) {
    width[1] = buff;
  }
  buff = column_width(distance_matrix->nb_row - 2 , between_cluster_distance + distance_matrix->nb_row);
  if (buff > width[1]) {
    width[1] = buff;
  }
  buff = column_width(distance_matrix->nb_row - 1 , max_within_cluster_distance + distance_matrix->nb_row);
  if (scale == DIAMETER) {
    width[2] = buff;
  }
  if (buff > width[1]) {
    width[1] = buff;
  }
  buff = column_width(distance_matrix->nb_row - 2 , min_between_cluster_distance + distance_matrix->nb_row);
  if (buff > width[1]) {
    width[1] = buff;
  }
  width[1] += ASCII_SPACE;

  // writing of (i) the distances between child clusters, (ii) the wihin-cluster distances,
  // (iii) the between-cluster distances, (iv) the distances between the farthest individuals
  // within a cluster (diameter), (v) the distances between the closest individuals belonging to
  // two different clusters (separation) and (vi) the cluster compositions

  os << "\n        | " << STAT_label[STATL_CHILD] << " "  << STAT_label[STATL_CLUSTER] << " "
     << STAT_label[STATL_DISTANCE] << " | " << STAT_label[STATL_WITHIN] << "-" << STAT_label[STATL_CLUSTER] << " "
     << STAT_label[STATL_DISTANCE] << " | " << STAT_label[STATL_BETWEEN] << "-" << STAT_label[STATL_CLUSTER] << " "
     << STAT_label[STATL_DISTANCE] << " | " << STAT_label[STATL_DIAMETER] << " | " << STAT_label[STATL_SEPARATION]
     << " | " << STAT_label[STATL_COMPOSITION] << endl;
  for (i = distance_matrix->nb_row;i < 2 * distance_matrix->nb_row - 1;i++) {
    os << STAT_label[STATL_STEP] << setw(width[0]) << i + 1 - distance_matrix->nb_row
       << setw(width[1]) << child_distance[i]
       << setw(width[1]) << within_cluster_distance[i];
    if (i < 2 * distance_matrix->nb_row - 2) {
      os << setw(width[1]) << between_cluster_distance[i]
         << setw(width[1]) << max_within_cluster_distance[i]
         << setw(width[1]) << min_between_cluster_distance[i];
    }
    else {
      os << setw(width[1]) << " "
         << setw(width[1]) << max_within_cluster_distance[i]
         << setw(width[1]) << " ";
    }
    os << "   ";
    for (j = 0;j < cluster_nb_pattern[i];j++) {
      os << distance_matrix->row_identifier[cluster_pattern[i][j]];
      if (j < cluster_nb_pattern[i] - 1) {
        os << ", ";
      }
    }
    os << endl;
  }
  os << endl;

  ordered_distance = distance_ordering();

  // computation of the column width

  max_identifier = 0;
  for (i = 0;i < distance_matrix->nb_row;i++) {
    if (distance_matrix->row_identifier[i] > max_identifier) {
      max_identifier = distance_matrix->row_identifier[i];
    }
  }

  width[0] = column_width(max_identifier);

  // writing of ordered distances

  os << "\n" << STAT_label[STATL_DENDROGRAM_SCALE] << ": ";
  switch (scale) {
  case CHILD_CLUSTER_DISTANCE :
    os << STAT_label[STATL_CHILD] << " "  << STAT_label[STATL_CLUSTER] << " "
       << STAT_label[STATL_DISTANCE] << endl;
    break;
  case DIAMETER :
    os << STAT_label[STATL_DIAMETER] << endl;
    break;
  }

  for (i = 0;i < distance_matrix->nb_row - 1;i++) {
    os << setw(width[0]) << distance_matrix->row_identifier[cluster_pattern[nb_cluster - 1][i]]
       << setw(width[2]) << " ";
  }
  os << setw(width[0]) << distance_matrix->row_identifier[cluster_pattern[nb_cluster - 1][distance_matrix->nb_row - 1]] << endl;
  for (i = 0;i < distance_matrix->nb_row - 1;i++) {
    os << setw(width[0]) << " "
       << setw(width[2]) << ordered_distance[i];
  }
  os << endl;

  // writing of the dendrogram

  switch (scale) {
  case CHILD_CLUSTER_DISTANCE :
    distance = child_distance;
    break;
  case DIAMETER :
    distance = max_within_cluster_distance;
    break;
  }

  min_diff_distance = distance[2 * distance_matrix->nb_row - 2];
  for (i = distance_matrix->nb_row;i < 2 * distance_matrix->nb_row - 2;i++) {
    if ((distance[i] > distance[i - 1]) && (distance[i] - distance[i - 1] < min_diff_distance)) {
      min_diff_distance = distance[i] - distance[i - 1];
    }
  }

  max_nb_character = MIN(2 * (int)round(distance[2 * distance_matrix->nb_row - 2] / min_diff_distance) ,
                         LINE_NB_CHARACTER);
  if (max_nb_character == 0) {
    max_nb_character = LINE_NB_CHARACTER;
  }

  if (max_nb_character == LINE_NB_CHARACTER) {
    min_diff_distance = 2 * distance[2 * distance_matrix->nb_row - 2] / LINE_NB_CHARACTER;
  }

  nb_character = new int[distance_matrix->nb_row - 1];
  for (i = 0;i < distance_matrix->nb_row - 1;i++) {
    nb_character[i] = 2 * (int)round(ordered_distance[i] / min_diff_distance);
  }

  os << "\n";
  for (i = 0;i < distance_matrix->nb_row;i++) {
    os << setw(width[0]) << distance_matrix->row_identifier[cluster_pattern[nb_cluster - 1][i]]
       << " ";

    for (j = 0;j < (i == 0 ? max_nb_character : nb_character[i - 1]) - 1;j++) {
      os << "_";
    }

    if (i > 0) {
      os << "|";

      min_distance = ordered_distance[i - 1];
      previous_nb_character = (nb_character[i - 1] == 0 ? 1 : nb_character[i - 1]);

      for (j = i;j < distance_matrix->nb_row - 1;j++) {
        if (ordered_distance[j] > min_distance) {
          for (k = previous_nb_character;k < nb_character[j] - 1;k++) {
            os << " ";
          }
          if (nb_character[j] > previous_nb_character) {
            os << "|";
          }

          min_distance = ordered_distance[j];
          previous_nb_character = nb_character[j];
        }
      }
    }
    os << endl;

    os << setw(width[0]) << " " << " ";

    if (nb_character[i] == 0) {
      os << "|";
    }

    min_distance = 0.;
    previous_nb_character = (nb_character[i] == 0 ? 1 : 0);

    if (i < distance_matrix->nb_row - 1) {
      for (j = i;j < distance_matrix->nb_row - 1;j++) {
        if (ordered_distance[j] > min_distance) {
          for (k = previous_nb_character;k < nb_character[j] - 1;k++) {
            os << " ";
          }
          if (nb_character[j] > previous_nb_character) {
            os << "|";
          }

          min_distance = ordered_distance[j];
          previous_nb_character = nb_character[j];
        }
      }

      os << endl;
    }
  }

  delete [] ordered_distance;
  delete [] nb_character;

  os << "\n" << STAT_label[STATL_CHILD_CLUSTER_DISTANCE_COEFF] << ": "
     << coefficient_computation(CHILD_CLUSTER_DISTANCE) << endl;
  os << STAT_label[STATL_DIAMETER_COEFF] << ": "
     << coefficient_computation(DIAMETER) << endl;

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Dendrogram object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool Dendrogram::ascii_write(StatError &error , const string path ,
                             bool exhaustive) const

{
  bool status;
  ofstream out_file(path.c_str());


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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Dendrogram object in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool Dendrogram::spreadsheet_write(StatError &error , const string path) const

{
  bool status;
  int i , j;
  double *ordered_distance;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    out_file << distance_matrix->nb_row << " " << distance_matrix->label << "s" << endl;

    // writing of (i) the distances between child clusters, (ii) the wihin-cluster distances,
    // (iii) the between-cluster distances, (iv) the distances between the farthest individuals
    // within a cluster (diameter), (v) the distances between the closest individuals belonging to
    // two different clusters (separation) and (vi) the cluster compositions

    out_file << "\n\t" << STAT_label[STATL_CHILD] << " " << STAT_label[STATL_CLUSTER] << " "
             << STAT_label[STATL_DISTANCE] << "\t" << STAT_label[STATL_WITHIN] << "-" << STAT_label[STATL_CLUSTER] << " "
             << STAT_label[STATL_DISTANCE] << "\t" << STAT_label[STATL_BETWEEN] << "-" << STAT_label[STATL_CLUSTER] << " "
             << STAT_label[STATL_DISTANCE] << "\t" << STAT_label[STATL_DIAMETER] << "\t" << STAT_label[STATL_SEPARATION]
             << "\t" << STAT_label[STATL_COMPOSITION] << endl;
    for (i = distance_matrix->nb_row;i < 2 * distance_matrix->nb_row - 1;i++) {
      out_file << STAT_label[STATL_STEP] << " " << i + 1 - distance_matrix->nb_row << "\t"
               << child_distance[i] << "\t" << within_cluster_distance[i] << "\t";
      if (i < 2 * distance_matrix->nb_row - 2) {
        out_file << between_cluster_distance[i] << "\t" << max_within_cluster_distance[i] << "\t"
                 << min_between_cluster_distance[i] << "\t";
      }
      else {
        out_file << "\t" << max_within_cluster_distance[i] << "\t\t";
      }
      for (j = 0;j < cluster_nb_pattern[i];j++) {
        out_file  << "\t" << distance_matrix->row_identifier[cluster_pattern[i][j]];
      }
      out_file << endl;
    }
    out_file << endl;

    ordered_distance = distance_ordering();

    // writing of the ordered distances

    out_file << "\n" << STAT_label[STATL_DENDROGRAM_SCALE] << "\t";
    switch (scale) {
    case CHILD_CLUSTER_DISTANCE :
      out_file << STAT_label[STATL_CHILD] << " "  << STAT_label[STATL_CLUSTER] << " "
               << STAT_label[STATL_DISTANCE] << endl;
      break;
    case DIAMETER :
      out_file << STAT_label[STATL_DIAMETER] << endl;
      break;
    }

    for (i = 0;i < distance_matrix->nb_row - 1;i++) {
      out_file << distance_matrix->row_identifier[cluster_pattern[nb_cluster - 1][i]] << "\t\t";
    }
    out_file << distance_matrix->row_identifier[cluster_pattern[nb_cluster - 1][distance_matrix->nb_row - 1]] << endl;
    for (i = 0;i < distance_matrix->nb_row - 1;i++) {
      out_file << "\t" << ordered_distance[i] << "\t";
    }
    out_file << endl;

    delete [] ordered_distance;

    out_file << "\n" << STAT_label[STATL_CHILD_CLUSTER_DISTANCE_COEFF] << "\t"
             << coefficient_computation(CHILD_CLUSTER_DISTANCE) << endl;
    out_file << STAT_label[STATL_DIAMETER_COEFF] << "\t"
             << coefficient_computation(DIAMETER) << endl;
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the dendrogram on the basis of the cluster compositions.
 */
/*--------------------------------------------------------------*/

void Dendrogram::tree_computation()

{
  int i , j , k;


  for (i = 0;i < nb_cluster - 1;i++) {
    for (j = MAX(i + 1 , distance_matrix->nb_row);j < nb_cluster;j++) {
      if (cluster_pattern[i][0] == cluster_pattern[j][0]) {
        child[j][0] = i;
        parent[i] = j;
        break;
      }

      for (k = cluster_nb_pattern[j] - 1;k >= 0;k--) {
        if (cluster_pattern[i][cluster_nb_pattern[i] - 1] == cluster_pattern[j][k]) {
          child[j][1] = i;
          parent[i] = j;
          break;
        }
      }

      if (k >= 0) {
        break;
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the agglomerative/divisive coefficient (Kaufman & Rousseeuw, pp. 212, 263).
 *
 *  \param[in] iscale cluster distance scale,
 *
 *  \return agglomerative/divisive coefficient.
 */
/*--------------------------------------------------------------*/

double Dendrogram::coefficient_computation(cluster_scale iscale) const

{
  int i;
  double coeff , *distance;


/*  if (iscale == I_DEFAULT) {
    iscale = scale;
  } */

  switch (iscale) {
  case CHILD_CLUSTER_DISTANCE :
    distance = child_distance;
    break;
  case DIAMETER :
    distance = max_within_cluster_distance;
    break;
  }

  coeff = 0.;
  for (i = 0;i < distance_matrix->nb_row;i++) {
    coeff += distance[parent[i]];
  }

  coeff = 1. - coeff / (distance[nb_cluster - 1] * distance_matrix->nb_row);

  return coeff;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the agglomerative/divisive coefficient.
 *
 *  \return agglomerative/divisive coefficient.
 */
/*--------------------------------------------------------------*/

double Dendrogram::coefficient_computation() const

{
  coefficient_computation(scale);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Agglomerative hierarchical clustering algorithm (Kaufman & Rousseeuw, pp. 199-208).
 *
 *  \param[in] strategy  algorithm type,
 *  \param[in] criterion cluster merging criterion.
 *
 *  \return              Dendrogram object.
 */
/*--------------------------------------------------------------*/

Dendrogram* DistanceMatrix::agglomerative_hierarchical_clustering(hierarchical_strategy strategy ,
                                                                  linkage criterion) const

{
  int i , j , k;
  int index , index1 , index2 , icluster , *pattern_index , **cluster_pattern;
  double min_distance , *cumul_distance , **normalized_cluster_distance ,
         **normalized_pattern_distance , **max_between_cluster_distance;
  DistanceMatrix *dist_matrix;
  Clusters *clusters;
  Dendrogram *dendrogram;


  dist_matrix = new DistanceMatrix(*this , (test_symmetry() ? COPY : SYMMETRIZATION));

# ifdef DEBUG
  if (!(test_symmetry())) {
    cout << *dist_matrix << endl;
  }
# endif

  dendrogram = new Dendrogram(*dist_matrix , (criterion != FARTHEST_NEIGHBOR ? CHILD_CLUSTER_DISTANCE : DIAMETER));

  // initialization of the data structures of the algorithm

  clusters = new Clusters(*dist_matrix , nb_row);

  for (i = 0;i < nb_row;i++) {
    clusters->assignment[i] = i;
  }

  cluster_pattern = new int*[nb_row];
  for (i = 0;i < nb_row;i++) {
    clusters->cluster_nb_pattern[i] = 1;
    cluster_pattern[i] = new int[nb_row - i];
    cluster_pattern[i][0] = i;
  }

  for (i = 0;i < nb_row;i++) {
    for (j = 0;j < nb_row;j++) {
      clusters->distance[i][j] = dist_matrix->distance[i][j];
    }
  }

  for (i = 0;i < nb_row;i++) {
    for (j = 0;j < nb_row;j++) {
      clusters->length[i][j] = dist_matrix->length[i][j];
    }
  }

  normalized_pattern_distance = new double*[nb_row];
  for (i = 0;i < nb_row;i++) {
    normalized_pattern_distance[i] = new double[nb_row];
    for (j = 0;j < nb_row;j++) {
      normalized_pattern_distance[i][j] = clusters->distance[i][j];
      if ((clusters->distance[i][j] != -D_INF) && (clusters->length[i][j] > 1)) {
        normalized_pattern_distance[i][j] /= clusters->length[i][j];
      }
    }
  }

  normalized_cluster_distance = new double*[nb_row];
  for (i = 0;i < nb_row;i++) {
    normalized_cluster_distance[i] = new double[nb_row];
    for (j = 0;j < nb_row;j++) {
      normalized_cluster_distance[i][j] = normalized_pattern_distance[i][j];
    }
  }

  if (criterion == FARTHEST_NEIGHBOR) {
    max_between_cluster_distance = new double*[nb_row];
    for (i = 0;i < nb_row;i++) {
      max_between_cluster_distance[i] = new double[nb_row];
    }
  }

  if (strategy == ORDERING) {
    cumul_distance = new double[nb_row];

    // sort by increasing average distance to the other individuals

    for (i = 0;i < nb_row;i++) {
      cumul_distance[i] = stat_tool::cumul_distance_computation(nb_column , dist_matrix->distance[i]);
      if (cumul_distance[i] != -D_INF) {
        cumul_distance[i] /= cumul_computation(1 , nb_column , dist_matrix->length + i);
      }
    }

    pattern_index = pattern_sort(nb_row , cumul_distance);
  }

  for (i = 0;i < nb_row - 1;i++) {
    if (i < nb_row - 2) {

#     ifdef DEBUG
      for (j = 0;j < nb_row - i;j++) {
        for (k = 0;k < nb_row - i;k++) {
          cout << normalized_cluster_distance[j][k] << "  ";
          if ((j == k) && (normalized_cluster_distance[j][k] == 0.)) {
            cout << "       ";
          }
        }
        cout << endl;
      }
      cout << endl;
#     endif

      // determination of the current minimum distance

      switch (strategy) {

      case AGGLOMERATIVE : {
        min_distance = -D_INF;

        switch (criterion) {

        case NEAREST_NEIGHBOR : {
          for (j = 0;j < nb_row - 1;j++) {
            for (k = j + 1;k < nb_row;k++) {
              if((clusters->assignment[j] != clusters->assignment[k]) &&
                 (normalized_pattern_distance[j][k] < min_distance)) {
                min_distance = normalized_pattern_distance[j][k];
                index1 = clusters->assignment[j];
                index2 = clusters->assignment[k];
              }
            }
          }

          if (index1 > index2) {
            index = index1;
            index1 = index2;
            index2 = index;
          }
          break;
        }

        case FARTHEST_NEIGHBOR : {
          for (j = 0;j < nb_row - i;j++) {
            for (k = 0;k < nb_row - i;k++) {
              max_between_cluster_distance[j][k] = 0.;
            }
          }

          for (j = 0;j < nb_row - 1;j++) {
            for (k = j + 1;k < nb_row;k++) {
              if((clusters->assignment[j] != clusters->assignment[k]) &&
                 (normalized_pattern_distance[j][k] > max_between_cluster_distance[clusters->assignment[j]][clusters->assignment[k]])) {
                max_between_cluster_distance[clusters->assignment[j]][clusters->assignment[k]] = normalized_pattern_distance[j][k];
              }
            }
          }

          for (j = 0;j < nb_row - i - 1;j++) {
            for (k = j + 1;k < nb_row - i;k++) {
              max_between_cluster_distance[j][k] = MAX(max_between_cluster_distance[j][k] , max_between_cluster_distance[k][j]);
              if (max_between_cluster_distance[j][k] < min_distance) {
                min_distance = max_between_cluster_distance[j][k];
                index1 = j;
                index2 = k;
              }
            }
          }
          break;
        }

        case AVERAGE_NEIGHBOR : {
          for (j = 0;j < nb_row - i - 1;j++) {
            for (k = j + 1;k < nb_row - i;k++) {
              if (normalized_cluster_distance[j][k] < min_distance) {
                min_distance = normalized_cluster_distance[j][k];
                index1 = j;
                index2 = k;
              }
            }
          }
          break;
        }
        }

#       ifdef DEBUG
        cout << "\nTest: " << index1 << " | " << index2 << " | " << min_distance << endl;
#       endif

        break;
      }

      case ORDERING : {
        if (i == 0) {
          if (pattern_index[0] < pattern_index[1]) {
            index1 = pattern_index[0];
            index2 = pattern_index[1];
          }
          else {
            index1 = pattern_index[1];
            index2 = pattern_index[0];
          }
        }

        else {
          for (j = 0;j < nb_row - i;j++) {
            if (pattern_index[i + 1] == cluster_pattern[j][0]) {
              if (index1 < j) {
                index2 = j;
              }
              else {
                index2 = index1;
                index1 = j;
              }
              break;
            }
          }
        }
        break;
      }
      }
    }

    else {
      index1 = 0;
      index2 = 1;
    }

    dendrogram->child_distance[nb_row + i] = normalized_cluster_distance[index1][index2];

    (clusters->nb_cluster)--;

    // update of the between-cluster distance matrix

    for (j = 0;j < nb_row - i - 1;j++) {
      if (j < index2) {
        clusters->distance[index1][j] += clusters->distance[index2][j];
        clusters->length[index1][j] += clusters->length[index2][j];
        if (j == index1) {
          clusters->distance[index1][j] += clusters->distance[j][index2] + clusters->distance[index2][index2];
          clusters->length[index1][j] += clusters->length[j][index2] + clusters->length[index2][index2];
        }
      }

      else {
        clusters->distance[index1][j] = clusters->distance[index1][j + 1] + clusters->distance[index2][j + 1];
        clusters->length[index1][j] = clusters->length[index1][j + 1] + clusters->length[index2][j + 1];
      }

      normalized_cluster_distance[index1][j] = clusters->distance[index1][j];
      if ((clusters->distance[index1][j] != -D_INF) && (clusters->length[index1][j] > 1)) {
        normalized_cluster_distance[index1][j] /= clusters->length[index1][j];
      }

      if (j != index1) {
        clusters->distance[j][index1] = clusters->distance[index1][j];
        clusters->length[j][index1] = clusters->length[index1][j];
        normalized_cluster_distance[j][index1] = normalized_cluster_distance[index1][j];
      }
    }

    for (j = 0;j < index2;j++) {
      if (j != index1) {
        for (k = index2;k < nb_row - i - 1;k++) {
          clusters->distance[j][k] = clusters->distance[j][k + 1];
          clusters->length[j][k] = clusters->length[j][k + 1];
          normalized_cluster_distance[j][k] = normalized_cluster_distance[j][k + 1];

          clusters->distance[k][j] = clusters->distance[k + 1][j];
          clusters->length[k][j] = clusters->length[k + 1][j];
          normalized_cluster_distance[k][j] = normalized_cluster_distance[k + 1][j];
        }
      }
    }

    for (j = index2;j < nb_row - i - 1;j++) {
      for (k = j;k < nb_row - i - 1;k++) {
        clusters->distance[j][k] = clusters->distance[j + 1][k + 1];
        clusters->length[j][k] = clusters->length[j + 1][k + 1];
        normalized_cluster_distance[j][k] = normalized_cluster_distance[j + 1][k + 1];

        if (k > j) {
          clusters->distance[k][j] = clusters->distance[k + 1][j + 1];
          clusters->length[k][j] = clusters->length[k + 1][j + 1];
          normalized_cluster_distance[k][j] = normalized_cluster_distance[k + 1][j + 1];
        }
      }
    }

    dendrogram->within_cluster_distance[nb_row + i] = normalized_cluster_distance[index1][index1];
    if (i < nb_row - 2) {
      dendrogram->between_cluster_distance[nb_row + i] = clusters->between_cluster_distance_computation(index1);
    }

    icluster = clusters->assignment[cluster_pattern[index2][0]];
    for (j = 0;j < clusters->cluster_nb_pattern[index2];j++) {
      clusters->assignment[cluster_pattern[index2][j]] = clusters->assignment[cluster_pattern[index1][0]];
    }
    for (j = cluster_pattern[index2][0] + 1;j < nb_row;j++) {
      if (clusters->assignment[j] > icluster) {
        (clusters->assignment[j])--;
      }
    }

    // update of the cluster compositions

    for (j = 0;j < clusters->cluster_nb_pattern[index2];j++) {
      cluster_pattern[index1][clusters->cluster_nb_pattern[index1] + j] = cluster_pattern[index2][j];
    }
    clusters->cluster_nb_pattern[index1] += clusters->cluster_nb_pattern[index2];

    for (j = index2;j < nb_row - i - 1;j++) {
      clusters->cluster_nb_pattern[j] = clusters->cluster_nb_pattern[j + 1];
      for (k = 0;k < clusters->cluster_nb_pattern[j];k++) {
        cluster_pattern[j][k] = cluster_pattern[j + 1][k];
      }
    }

    dendrogram->cluster_nb_pattern[nb_row + i] = clusters->cluster_nb_pattern[index1];
    dendrogram->cluster_pattern[nb_row + i] = new int[dendrogram->cluster_nb_pattern[nb_row + i]];
    for (j = 0;j < dendrogram->cluster_nb_pattern[nb_row + i];j++) {
      dendrogram->cluster_pattern[nb_row + i][j] = cluster_pattern[index1][j];
    }

#   ifdef DEBUG
    cout << "\nTest: " << clusters->cluster_nb_pattern[index1] << " | "
         << clusters->cluster_nb_pattern[clusters->assignment[cluster_pattern[index1][0]]] << " | "
         << clusters->assignment[cluster_pattern[index1][0]] << endl;
    for (j = 0;j < nb_row;j++) {
      cout << clusters->assignment[j] << " ";
    }
    cout << endl;
#   endif

    // computation of the cluster diameters and separations

    dendrogram->max_within_cluster_distance[nb_row + i] = clusters->max_within_cluster_distance_computation(normalized_pattern_distance , clusters->assignment[cluster_pattern[index1][0]]);
    if (i < nb_row - 2) {
      dendrogram->min_between_cluster_distance[nb_row + i] = clusters->min_between_cluster_distance_computation(normalized_pattern_distance , clusters->assignment[cluster_pattern[index1][0]]);
    }
  }

  // computation of the dendrogram

  dendrogram->tree_computation();

  delete dist_matrix;
  clusters->nb_cluster = nb_row;
  delete clusters;

  for (i = 0;i < nb_row;i++) {
    delete [] cluster_pattern[i];
  }
  delete [] cluster_pattern;

  for (i = 0;i < nb_row;i++) {
    delete [] normalized_pattern_distance[i];
  }
  delete [] normalized_pattern_distance;

  for (i = 0;i < nb_row;i++) {
    delete [] normalized_cluster_distance[i];
  }
  delete [] normalized_cluster_distance;

  if (criterion == FARTHEST_NEIGHBOR) {
    for (i = 0;i < nb_row;i++) {
      delete [] max_between_cluster_distance[i];
    }
    delete [] max_between_cluster_distance;
  }

  if (strategy == ORDERING) {
    delete [] cumul_distance;
    delete [] pattern_index;
  }

  return dendrogram;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Divisive hierarchical clustering algorithm (Kaufman & Rousseeuw, pp. 253-259).
 *
 *  \return Dendrogram object.
 */
/*--------------------------------------------------------------*/

Dendrogram* DistanceMatrix::divisive_hierarchical_clustering() const

{
  int i , j , k;
  int bnb_pattern , icluster , new_cluster , outlying_pattern , *passignment1 , *passignment2 ,
      *cluster_identifier , *cluster_index;
  double distance , max_cumul_distance , *cumul_distance , **normalized_distance;
  DistanceMatrix *dist_matrix , *step_dist_matrix;
  Clusters *clusters , *step_clusters;
  Dendrogram *dendrogram;


  dist_matrix = new DistanceMatrix(*this , (test_symmetry() ? COPY : SYMMETRIZATION));

# ifdef DEBUG
  if (!(test_symmetry())) {
    cout << *dist_matrix << endl;
  }
# endif

  dendrogram = new Dendrogram(*dist_matrix , DIAMETER);

  // initialization of the data structures of the algorithm

  clusters = new Clusters(*dist_matrix , nb_row);
  clusters->nb_cluster = 1;

  clusters->cluster_nb_pattern[0] = nb_row;
  passignment1 = clusters->assignment;
  for (i = 0;i < clusters->nb_pattern;i++) {
    *passignment1++ = 0;
  }

  normalized_distance = new double*[nb_row];
  for (i = 0;i < nb_row;i++) {
    normalized_distance[i] = new double[nb_column];
    for (j = 0;j < nb_column;j++) {
      normalized_distance[i][j] = dist_matrix->distance[i][j];
      if ((dist_matrix->distance[i][j] != -D_INF) && (dist_matrix->length[i][j] > 1)) {
        normalized_distance[i][j] /= dist_matrix->length[i][j];
      }
    }
  }

  cluster_identifier = new int[nb_row];
  cluster_index = new int[nb_row];
  cumul_distance = new double[nb_row];

  for (i = nb_row - 2;i >= 0;i--) {
    if (i == nb_row - 2) {
      dendrogram->max_within_cluster_distance[nb_row + i] = clusters->max_within_cluster_distance_computation(normalized_distance , 0);
      icluster = 0;
      for (j = 0;j < clusters->nb_pattern;j++) {
        cluster_index[j] = j;
      }

      step_dist_matrix = new DistanceMatrix(*dist_matrix);
    }

    else {

      // determination of the cluster of maximum diameter

      dendrogram->max_within_cluster_distance[nb_row + i] = 0.;
      for (j = 0;j <= nb_row - i - 2;j++) {
        distance = clusters->max_within_cluster_distance_computation(normalized_distance , j);
        if (distance > dendrogram->max_within_cluster_distance[nb_row + i]) {
          dendrogram->max_within_cluster_distance[nb_row + i] = distance;
          icluster = j;
        }
      }

      if (dendrogram->max_within_cluster_distance[nb_row + i] == 0.) {
        bnb_pattern = 1;
        for (j = 0;j <= nb_row - i - 2;j++) {
          if (clusters->cluster_nb_pattern[j] > bnb_pattern) {
            bnb_pattern = clusters->cluster_nb_pattern[j];
            icluster = j;
          }
        }
      }

      dendrogram->min_between_cluster_distance[nb_row + i] = clusters->min_between_cluster_distance_computation(normalized_distance , icluster);

#     ifdef DEBUG
      cout << "\nCluster to be splitted: " << icluster << endl;
#     endif

      if (dendrogram->max_within_cluster_distance[nb_row + i] > 0.) {

        // extraction of the distance matrix corresponding to the selected cluster

        passignment1 = clusters->assignment;
        j = 0;
        for (k = 0;k < clusters->nb_pattern;k++) {
          if (*passignment1++ == icluster) {
            cluster_index[j] = k;
            cluster_identifier[j++] = row_identifier[k];
          }
        }

        step_dist_matrix = new DistanceMatrix(*dist_matrix , clusters->cluster_nb_pattern[icluster] ,
                                              cluster_identifier);
      }
    }

    clusters->cluster_distance_computation_1();

    dendrogram->within_cluster_distance[nb_row + i] = clusters->distance[icluster][icluster];
    if ((clusters->distance[icluster][icluster] != -D_INF) && (clusters->length[icluster][icluster] > 1)) {
      dendrogram->within_cluster_distance[nb_row + i] /= clusters->length[icluster][icluster];
    }
    if (i < nb_row - 2) {
      dendrogram->between_cluster_distance[nb_row + i] = clusters->between_cluster_distance_computation(icluster);
    }

    dendrogram->cluster_nb_pattern[nb_row + i] = clusters->cluster_nb_pattern[icluster];
    dendrogram->cluster_pattern[nb_row + i] = new int[dendrogram->cluster_nb_pattern[nb_row + i]];

    if (dendrogram->max_within_cluster_distance[nb_row + i] > 0.) {
      step_clusters = new Clusters(*step_dist_matrix , 2);

      // determination of the most outlying individual within the selected cluster

      max_cumul_distance = 0.;

      for (j = 0;j < step_dist_matrix->nb_row;j++) {
        cumul_distance[j] = stat_tool::cumul_distance_computation(step_dist_matrix->nb_column , step_dist_matrix->distance[j]);
        if (cumul_distance[j] != -D_INF) {
          cumul_distance[j] /= cumul_computation(1 , step_dist_matrix->nb_column , step_dist_matrix->length + j);
        }
        if (cumul_distance[j] > max_cumul_distance) {
          max_cumul_distance = cumul_distance[j];
          outlying_pattern = j;
        }
      }

#     ifdef DEBUG
      cout << "\nOutlying individual: " << outlying_pattern << endl;
#     endif

      passignment2 = step_clusters->assignment;
      for (j = 0;j < step_clusters->nb_pattern;j++) {
        *passignment2++ = 0;
      }
      step_clusters->assignment[outlying_pattern] = 1;

      // partition in 2 clusters of the selected cluster

      step_clusters->algorithmic_step_1();
      step_clusters->cluster_distance_computation_1();

      // update of the cluster compositions

      for (j = 0;j < dendrogram->cluster_nb_pattern[nb_row + i];j++) {
        dendrogram->cluster_pattern[nb_row + i][j] = cluster_index[j];
      }
      dendrogram->child_distance[nb_row + i] = step_clusters->distance[0][1];
      if ((step_clusters->distance[0][1] != -D_INF) && (step_clusters->length[0][1] > 1)) {
        dendrogram->child_distance[nb_row + i] /= step_clusters->length[0][1];
      }

      (clusters->nb_cluster)++;

      new_cluster = (step_clusters->assignment[0] == 0 ? 1 : 0);
      clusters->cluster_nb_pattern[nb_row - i - 1] = step_clusters->cluster_nb_pattern[new_cluster];
      clusters->cluster_nb_pattern[icluster] -= step_clusters->cluster_nb_pattern[new_cluster];

      passignment1 = clusters->assignment;
      passignment2 = step_clusters->assignment;
      for (j = 0;j < clusters->nb_pattern;j++) {
        if (*passignment1 == icluster) {
          if (*passignment2 == new_cluster) {
            *passignment1 = nb_row - i - 1;
          }
          passignment2++;
        }
        passignment1++;
      }

      delete step_dist_matrix;
      delete step_clusters;
    }

    else {

      // update of the cluster compositions

      passignment1 = clusters->assignment;
      j = 0;
      for (k = 0;k < clusters->nb_pattern;k++) {
        if (*passignment1++ == icluster) {
          dendrogram->cluster_pattern[nb_row + i][j++] = k;
          outlying_pattern = k;
        }
      }
      dendrogram->child_distance[nb_row + i] = 0.;

      (clusters->nb_cluster)++;

      clusters->cluster_nb_pattern[nb_row - i - 1] = 1;
      (clusters->cluster_nb_pattern[icluster])--;
      clusters->assignment[outlying_pattern] = nb_row - i - 1;
    }
  }

  // computation of the dendrogram

  dendrogram->tree_computation();

  // re-ordering of clusters

  for (i = nb_row;i < 2 * nb_row - 1;i++) {
    for (j = 0;j < dendrogram->cluster_nb_pattern[dendrogram->child[i][0]];j++) {
      dendrogram->cluster_pattern[i][j] = dendrogram->cluster_pattern[dendrogram->child[i][0]][j];
    }
    for (j = 0;j < dendrogram->cluster_nb_pattern[dendrogram->child[i][1]];j++) {
      dendrogram->cluster_pattern[i][dendrogram->cluster_nb_pattern[dendrogram->child[i][0]] + j] =
      dendrogram->cluster_pattern[dendrogram->child[i][1]][j];
    }
  }

  delete dist_matrix;
  delete clusters;

  for (i = 0;i < nb_row;i++) {
    delete [] normalized_distance[i];
  }
  delete [] normalized_distance;

  delete [] cluster_identifier;
  delete [] cluster_index;
  delete [] cumul_distance;

  return dendrogram;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Hierarchical clustering algorithms.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] display   flag for displaying the hierarchical clustering results,
 *  \param[in] strategy  algorithm type (AGGLOMERATIVE/DIVISIVE/ORDERING),
 *  \param[in] criterion cluster merging criterion (agglomerative algorithm),
 *  \param[in] path      file path,
 *  \param[in] format    output format (ASCII/SPREADSHEET).
 *
 *  \return              Dendrogram object.
 */
/*--------------------------------------------------------------*/

bool DistanceMatrix::hierarchical_clustering(StatError &error , bool display ,
                                             hierarchical_strategy strategy , linkage criterion ,
                                             const string path , output_format format) const

{
  bool status = true;
  Dendrogram *dendrogram;


  error.init();

  if (nb_row != nb_column) {
    status = false;
    error.correction_update(STAT_error[STATR_MATRIX_STRUCTURE] , STAT_error[STATR_SQUARE_MATRIX]);
  }

  else {
    if (strategy != DIVISIVE) {
      dendrogram = agglomerative_hierarchical_clustering(strategy , criterion);
    }
    else {
      dendrogram = divisive_hierarchical_clustering();
    }

    // writing of results

    if (display) {
      dendrogram->ascii_write(cout);
    }

    if (!path.empty()) {
      switch (format) {
      case ASCII :
        status = dendrogram->ascii_write(error , path);
        break;
      case SPREADSHEET :
        status = dendrogram->spreadsheet_write(error , path);
        break;
      }
    }

    delete dendrogram;
  }

  return status;
}


};  // namespace stat_tool
