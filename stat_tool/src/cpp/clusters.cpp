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
#include <sstream>
#include <iomanip>
#include <cstring>

#include "distance_matrix.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {


extern int cumul_computation(int nb_row , int nb_column , int **value);
extern double cumul_distance_computation(int dim , double *distance);
extern int* pattern_sort(int nb_pattern , double *distance , int nb_sorted_pattern = I_DEFAULT);



/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the Clusters class.
 */
/*--------------------------------------------------------------*/

Clusters::Clusters()

{
  distance_matrix = NULL;

  nb_pattern = 0;
  nb_cluster = 0;

  cluster_nb_pattern = NULL;
  assignment = NULL;

  pattern_distance = NULL;
  pattern_length = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Clusters class.
 *
 *  \param[in] dist_matrix reference on a DistanceMatrix object,
 *  \param[in] inb_cluster number of clusters.
 */
/*--------------------------------------------------------------*/

Clusters::Clusters(const DistanceMatrix &dist_matrix , int inb_cluster)
:DistanceMatrix(dist_matrix , inb_cluster , STAT_label[STATL_CLUSTER])

{
  int i , j;


  distance_matrix = new DistanceMatrix(dist_matrix);

  nb_pattern = distance_matrix->nb_row;
  nb_cluster = inb_cluster;

  cluster_nb_pattern = new int[nb_cluster];
  for (i = 0;i < nb_cluster;i++) {
    cluster_nb_pattern[i] = 0;
  }

  assignment = new int[nb_pattern];
  for (i = 0;i < nb_pattern;i++) {
    assignment[i] = I_DEFAULT;
  }

  pattern_distance = new double*[nb_pattern];
  for (i = 0;i < nb_pattern;i++) {
    pattern_distance[i] = new double[nb_cluster];
    for (j = 0;j < nb_cluster;j++) {
      pattern_distance[i][j] = 0.;
    }
  }

  pattern_length = new int*[nb_pattern];
  for (i = 0;i < nb_pattern;i++) {
    pattern_length[i] = new int[nb_cluster];
    for (j = 0;j < nb_cluster;j++) {
      pattern_length[i][j] = 0;
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Clusters class.
 *
 *  \param[in] dist_matrix         reference on a DistanceMatrix object,
 *  \param[in] inb_cluster         number of clusters,
 *  \param[in] icluster_nb_pattern pointer on the cluster sizes,
 *  \param[in] cluster_pattern     pointer on the cluster compositions.
 */
/*--------------------------------------------------------------*/

Clusters::Clusters(const DistanceMatrix &dist_matrix , int inb_cluster ,
                   int *icluster_nb_pattern , int **cluster_pattern)
:DistanceMatrix(dist_matrix , inb_cluster , STAT_label[STATL_CLUSTER])

{
  int i , j , k;


  distance_matrix = new DistanceMatrix(dist_matrix);

  nb_pattern = distance_matrix->nb_row;
  nb_cluster = inb_cluster;

  cluster_nb_pattern = new int[nb_cluster];
  for (i = 0;i < nb_cluster;i++) {
    cluster_nb_pattern[i] = icluster_nb_pattern[i];
  }

  assignment = new int[nb_pattern];
  for (i = 0;i < nb_pattern;i++) {
    for (j = 0;j < nb_cluster;j++) {
      for (k = 0;k < cluster_nb_pattern[j];k++) {
        if (cluster_pattern[j][k] == distance_matrix->row_identifier[i]) {
          assignment[i] = j;
          break;
        }
      }
      if (k < cluster_nb_pattern[j]) {
        break;
      }
    }
  }

  pattern_distance = new double*[nb_pattern];
  for (i = 0;i < nb_pattern;i++) {
    pattern_distance[i] = new double[nb_cluster + 1];
    for (j = 0;j <= nb_cluster;j++) {
      pattern_distance[i][j] = 0.;
    }
  }

  pattern_length = new int*[nb_pattern];
  for (i = 0;i < nb_pattern;i++) {
    pattern_length[i] = new int[nb_cluster + 1];
    for (j = 0;j <= nb_cluster;j++) {
      pattern_length[i][j] = 0;
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a Clusters object.
 *
 *  \param[in] clusters reference on a Clusters object.
 */
/*--------------------------------------------------------------*/

void Clusters::copy(const Clusters &clusters)

{
  int i , j;


  distance_matrix = new DistanceMatrix(*(clusters.distance_matrix));

  nb_pattern = clusters.nb_pattern;
  nb_cluster = clusters.nb_cluster;

  cluster_nb_pattern = new int[nb_cluster];
  for (i = 0;i < nb_cluster;i++) {
    cluster_nb_pattern[i] = clusters.cluster_nb_pattern[i];
  }

  assignment = new int[nb_pattern];
  for (i = 0;i < nb_pattern;i++) {
    assignment[i] = clusters.assignment[i];
  }

  pattern_distance = new double*[nb_pattern];
  for (i = 0;i < nb_pattern;i++) {
    pattern_distance[i] = new double[nb_cluster + 1];
    for (j = 0;j <= nb_cluster;j++) {
      pattern_distance[i][j] = clusters.pattern_distance[i][j];
    }
  }

  pattern_length = new int*[nb_pattern];
  for (i = 0;i < nb_pattern;i++) {
    pattern_length[i] = new int[nb_cluster + 1];
    for (j = 0;j <= nb_cluster;j++) {
      pattern_length[i][j] = clusters.pattern_length[i][j];
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the data members of a Clusters object.
 */
/*--------------------------------------------------------------*/

void Clusters::remove()

{
  int i;


  delete distance_matrix;

  delete [] cluster_nb_pattern;
  delete [] assignment;

  if (pattern_distance) {
    for (i = 0;i < nb_pattern;i++) {
      delete [] pattern_distance[i];
    }
    delete [] pattern_distance;
  }

  if (pattern_length) {
    for (i = 0;i < nb_pattern;i++) {
      delete [] pattern_length[i];
    }
    delete [] pattern_length;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the Clusters class.
 */
/*--------------------------------------------------------------*/

Clusters::~Clusters()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the Clusters class.
 *
 *  \param[in] clusters reference on a Clusters object.
 *
 *  \return             Clusters object.
 */
/*--------------------------------------------------------------*/

Clusters& Clusters::operator=(const Clusters &clusters)

{
  if (&clusters != this) {
    remove();
    DistanceMatrix::remove();

    DistanceMatrix::copy(clusters);
    copy(clusters);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing on a line of a Clusters object.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& Clusters::line_write(ostream &os) const

{
  os << nb_pattern << " " << distance_matrix->label << "s" << "   "
     << nb_cluster << " " << STAT_label[STATL_CLUSTERS];

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Sort of the individuals of a cluster by increasing average distance to
 *         another individual of the cluster.
 *
 *  \param[in] cluster cluster index.
 *
 *  \return            sorted individual indices.
 */
/*--------------------------------------------------------------*/

int* Clusters::pattern_sort(int cluster) const

{
  bool *selected_pattern;
  int i , j;
  int *index;
  double min_distance , *normalized_distance;


  index = new int[cluster_nb_pattern[cluster]];

  selected_pattern = new bool[nb_pattern];
  normalized_distance = new double[nb_pattern];

  for (i = 0;i < nb_pattern;i++) {
    if (assignment[i] == cluster) {
      selected_pattern[i] = false;

      normalized_distance[i] = pattern_distance[i][cluster];
      if ((pattern_distance[i][cluster] != -D_INF) && (pattern_length[i][cluster] > 1)) {
        normalized_distance[i] /= pattern_length[i][cluster];
      }
    }

    else {
      selected_pattern[i] = true;
    }
  }

  for (i = 0;i < cluster_nb_pattern[cluster];i++) {
    min_distance = -D_INF * 10;
    for (j = 0;j < nb_pattern;j++) {
      if ((!selected_pattern[j]) && (normalized_distance[j] < min_distance)) {
        min_distance = normalized_distance[j];
        index[i] = j;
      }
    }

    selected_pattern[index[i]] = true;
  }

  delete [] selected_pattern;
  delete [] normalized_distance;

  return index;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Clusters object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& Clusters::ascii_write(ostream &os , bool exhaustive) const

{
  int i , j;
  int buff , max_identifier , neighbor_cluster , most_distant_pattern , neighbor_pattern ,
      isolation , *index , *order , width[3];
  double **normalized_pattern_distance , **normalized_pattern_cluster_distance ,
         **normalized_cluster_distance;
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::right , ios::adjustfield);

  if (exhaustive) {
    order = new int[nb_pattern];
  }

  // writing of the cluster compositions

  for (i = 0;i < nb_cluster;i++) {
    index = pattern_sort(i);

    os << STAT_label[STATL_CLUSTER] << " " << i + 1 << " (" << cluster_nb_pattern[i] << " " << distance_matrix->label;
    if (cluster_nb_pattern[i] > 1) {
      os << "s";
    }
    os << "): ";

    for (j = 0;j < cluster_nb_pattern[i];j++) {
      os << distance_matrix->row_identifier[index[j]];
      if (j < cluster_nb_pattern[i] - 1) {
        os << ", ";
      }

      if (exhaustive) {
        *order++ = index[j];
      }
    }
    os << endl;

    delete [] index;
  }

  normalized_pattern_distance = new double*[nb_pattern];
  for (i = 0;i < nb_pattern;i++) {
    normalized_pattern_distance[i] = new double[nb_pattern];
    for (j = 0;j < nb_pattern;j++) {
      normalized_pattern_distance[i][j] = distance_matrix->distance[i][j];
      if ((distance_matrix->distance[i][j] != -D_INF) && (distance_matrix->length[i][j] > 1)) {
        normalized_pattern_distance[i][j] /= distance_matrix->length[i][j];
      }
    }
  }

  // computation of the column widths

  width[2] = 0;
  for (i = 0;i < nb_pattern;i++) {
    buff = column_width(nb_pattern , normalized_pattern_distance[i]);
    if (buff > width[2]) {
      width[2] = buff;
    }
  }
  width[2] += ASCII_SPACE;

  if (exhaustive) {
    order -= nb_pattern;

    normalized_pattern_cluster_distance = new double*[nb_pattern];
    for (i = 0;i < nb_pattern;i++) {
      normalized_pattern_cluster_distance[i] = new double[nb_cluster];
      for (j = 0;j < nb_cluster;j++) {
        normalized_pattern_cluster_distance[i][j] = pattern_distance[i][j];
        if ((pattern_distance[i][j] != -D_INF) && (pattern_length[i][j] > 1)) {
          normalized_pattern_cluster_distance[i][j] /= pattern_length[i][j];
        }
      }
    }

    // computation of the column widths

    max_identifier = 0;
    for (i = 0;i < nb_pattern;i++) {
      if (distance_matrix->row_identifier[i] > max_identifier) {
        max_identifier = distance_matrix->row_identifier[i];
      }
    }

    width[0] = column_width(max_identifier);

    width[1] = 0;
    for (i = 0;i < nb_pattern;i++) {
      buff = column_width(nb_cluster , normalized_pattern_cluster_distance[i]);
      if (buff > width[1]) {
        width[1] = buff;
      }
    }
    width[1] += ASCII_SPACE;

    // writing of the individual/cluster distance matrix

    os << "\n" << distance_matrix->label << "/" << STAT_label[STATL_CLUSTER] << " "
       << STAT_label[STATL_DISTANCE] << " " << STAT_label[STATL_MATRIX] << endl;

    os << "\n";
    for (i = 0;i < distance_matrix->label_size - 1;i++) {
      os << " ";
    }
    os << setw(width[0]) << " " << " ";
    for (i = 0;i < nb_cluster;i++) {
      os << " | " << STAT_label[STATL_CLUSTER] << " " << i + 1;
    }
    os << endl;

    for (i = 0;i < nb_pattern;i++) {
      os << distance_matrix->label << " " << setw(width[0]) << distance_matrix->row_identifier[order[i]] << " ";
      for (j = 0;j < nb_cluster;j++) {
        os << setw(width[1]) << normalized_pattern_cluster_distance[order[i]][j];
      }
      os << endl;
    }

    // writing of the distances between an individual and, (i) the cluster to which this individual is assigned,
    // (ii) the neighboring cluster, (iii) the farthest individual within the cluster,
    // (iv) the closest individual belonging to another cluster

    os << "\n";
    for (i = 0;i < distance_matrix->label_size - 1;i++) {
      os << " ";
    }
    os << setw(width[0]) << " " << "  | " << STAT_label[STATL_WITHIN] << "-" << STAT_label[STATL_CLUSTER] << " "
       << STAT_label[STATL_DISTANCE] << " | " << STAT_label[STATL_NEIGHBOR] << " " << STAT_label[STATL_CLUSTER]
       << " | " << STAT_label[STATL_NEIGHBOR] << " " << STAT_label[STATL_CLUSTER] << " " << STAT_label[STATL_DISTANCE]
       << " | " << STAT_label[STATL_MOST_DISTANT] << " " << distance_matrix->label << " | "
       << STAT_label[STATL_MOST_DISTANT] << " " << distance_matrix->label  << " " << STAT_label[STATL_DISTANCE]
       << " | " << STAT_label[STATL_NEIGHBOR] << " " << distance_matrix->label << " | "
       << STAT_label[STATL_NEIGHBOR] << " " << distance_matrix->label  << " " << STAT_label[STATL_DISTANCE] << endl;

    for (i = 0;i < nb_pattern;i++) {
      neighbor_cluster = neighbor_pattern_cluster_selection(normalized_pattern_cluster_distance , order[i]);
      most_distant_pattern = most_distant_pattern_selection(normalized_pattern_distance , order[i]);
      neighbor_pattern = neighbor_pattern_selection(normalized_pattern_distance , order[i]);

      os << distance_matrix->label << " " << setw(width[0]) << distance_matrix->row_identifier[order[i]] << " "
         << setw(width[1]) << normalized_pattern_cluster_distance[order[i]][assignment[order[i]]]
         << "  " << setw(width[0]) << neighbor_cluster + 1
         << setw(width[1]) << normalized_pattern_cluster_distance[order[i]][neighbor_cluster]
         << "  " << setw(width[0]) << distance_matrix->row_identifier[most_distant_pattern]
         << setw(width[2]) << normalized_pattern_distance[order[i]][most_distant_pattern]
         << "  " << setw(width[0]) << distance_matrix->row_identifier[neighbor_pattern]
         << setw(width[2]) << normalized_pattern_distance[order[i]][neighbor_pattern] << endl;
    }

    delete [] order;

    for (i = 0;i < nb_pattern;i++) {
      delete [] normalized_pattern_cluster_distance[i];
    }
    delete [] normalized_pattern_cluster_distance;
  }

  // computation of within-cluster and between-cluster distances

  normalized_cluster_distance = new double*[nb_cluster];
  for (i = 0;i < nb_cluster;i++) {
    normalized_cluster_distance[i] = new double[nb_cluster];
    for (j = 0;j < nb_cluster;j++) {
      normalized_cluster_distance[i][j] = distance[i][j];
      if ((distance[i][j] != -D_INF) && (length[i][j] > 1)) {
        normalized_cluster_distance[i][j] /= length[i][j];
      }
    }
  }

  // computation of the column widths

  width[0] = column_width(nb_cluster);

  width[1] = 0;
  for (i = 0;i < nb_cluster;i++) {
    buff = column_width(nb_cluster , normalized_cluster_distance[i]);
    if (buff > width[1]) {
      width[1] = buff;
    }
  }
  width[1] += ASCII_SPACE;

  // writing of the between-cluster distance matrix

  os << "\n" << STAT_label[STATL_CLUSTER] << " " << STAT_label[STATL_DISTANCE] << " "
     << STAT_label[STATL_MATRIX] << endl;

  os << "\n        " << setw(width[0]) << " ";
  for (i = 0;i < nb_cluster;i++) {
    os << " | " << STAT_label[STATL_CLUSTER] << " " << i + 1;
  }
  os << endl;

  for (i = 0;i < nb_cluster;i++) {
    os << STAT_label[STATL_CLUSTER] << " " << setw(width[0]) << i + 1 << " ";
    for (j = 0;j < nb_cluster;j++) {
      os << setw(width[1]) << normalized_cluster_distance[i][j];
    }
    os << endl;
  }

  // for each cluster, writing of (i) the within-cluster distance, (ii) the between-cluster distance,
  // (iii) the largest distance between individuals of a cluster (diameter),
  // (iv) the smallest distance between individuals, one belonging to the cluster and
  // the other to another cluster (separation)

  os << "\n        " << setw(width[0]) << " " << " | "
     << STAT_label[STATL_WITHIN] << "-" << STAT_label[STATL_CLUSTER] << " " << STAT_label[STATL_DISTANCE] << " | "
     << STAT_label[STATL_BETWEEN] << "-" << STAT_label[STATL_CLUSTER] << " " << STAT_label[STATL_DISTANCE] << " | "
     << STAT_label[STATL_DIAMETER] << " | " << STAT_label[STATL_SEPARATION] << endl;

  for (i = 0;i < nb_cluster;i++) {
    os << STAT_label[STATL_CLUSTER] << " " << setw(width[0]) << i + 1 << " "
       << setw(width[1]) << normalized_cluster_distance[i][i]
       << setw(width[1]) << between_cluster_distance_computation(i)
       << setw(width[2]) << max_within_cluster_distance_computation(normalized_pattern_distance , i)
       << setw(width[2]) << min_between_cluster_distance_computation(normalized_pattern_distance , i) << endl;
  }

  for (i = 0;i < nb_cluster;i++) {
    os << "\n" << STAT_label[STATL_CLUSTER] << " " << i + 1 << ": ";
    isolation = isolation_property(normalized_pattern_distance , i , CLUSTER_SCALE);
    if (isolation) {
      os << STAT_label[STATL_ISOLATED];
    }
    else {
      isolation = isolation_property(normalized_pattern_distance , i , INDIVIDUAL);
      if (isolation) {
        os << STAT_label[STATL_ISOLATED] << " (" << STAT_label[STATL_PATTERN_LEVEL] << ")";
      }
      else {
        os << STAT_label[STATL_NON_ISOLATED];
      }
    }
  }
  os << endl;

  for (i = 0;i < nb_pattern;i++) {
    delete [] normalized_pattern_distance[i];
  }
  delete [] normalized_pattern_distance;

  for (i = 0;i < nb_cluster;i++) {
    delete [] normalized_cluster_distance[i];
  }
  delete [] normalized_cluster_distance;

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Clusters object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool Clusters::ascii_write(StatError &error , const string path ,
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
 *  \brief Writing of a Clusters object in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool Clusters::spreadsheet_write(StatError &error , const string path) const

{
  bool status;
  int i , j;
  int isolation , neighbor_cluster , most_distant_pattern , neighbor_pattern , *index , *order;
  double **normalized_pattern_distance , **normalized_pattern_cluster_distance ,
         **normalized_cluster_distance;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    // writing of the cluster compositions

    order = new int[nb_pattern];

    for (i = 0;i < nb_cluster;i++) {
      index = pattern_sort(i);

      out_file << STAT_label[STATL_CLUSTER] << " " << i + 1 << " (" << cluster_nb_pattern[i] << " " << distance_matrix->label;
      if (cluster_nb_pattern[i] > 1) {
        out_file << "s";
      }
      out_file << ")";

      for (j = 0;j < cluster_nb_pattern[i];j++) {
        out_file << "\t" << distance_matrix->row_identifier[index[j]];

        *order++ = index[j];
      }
      out_file << endl;

      delete [] index;
    }

    order -= nb_pattern;

    normalized_pattern_distance = new double*[nb_pattern];
    for (i = 0;i < nb_pattern;i++) {
      normalized_pattern_distance[i] = new double[nb_pattern];
      for (j = 0;j < nb_pattern;j++) {
        normalized_pattern_distance[i][j] = distance_matrix->distance[i][j];
        if ((distance_matrix->distance[i][j] != -D_INF) && (distance_matrix->length[i][j] > 1)) {
          normalized_pattern_distance[i][j] /= distance_matrix->length[i][j];
        }
      }
    }

    normalized_pattern_cluster_distance = new double*[nb_pattern];
    for (i = 0;i < nb_pattern;i++) {
      normalized_pattern_cluster_distance[i] = new double[nb_cluster + 2];
      for (j = 0;j < nb_cluster;j++) {
        normalized_pattern_cluster_distance[i][j] = pattern_distance[i][j];
        if ((pattern_distance[i][j] != -D_INF) && (pattern_length[i][j] > 1)) {
          normalized_pattern_cluster_distance[i][j] /= pattern_length[i][j];
        }
      }
    }

    // writing of the individual/cluster distance matrix

    out_file << "\n" << distance_matrix->label << "/" << STAT_label[STATL_CLUSTER] << " "
             << STAT_label[STATL_DISTANCE] << " " << STAT_label[STATL_MATRIX] << endl;

    out_file << "\n";
    for (i = 0;i < nb_cluster;i++) {
      out_file << "\t" << STAT_label[STATL_CLUSTER] << " " << i + 1;
    }
    out_file << endl;

    for (i = 0;i < nb_pattern;i++) {
      out_file << distance_matrix->label << " " << distance_matrix->row_identifier[order[i]];
      for (j = 0;j < nb_cluster;j++) {
        out_file << "\t" << normalized_pattern_cluster_distance[order[i]][j];
      }
      out_file << endl;
    }

    // writing of the distances between an individual and, (i) the cluster to which this individual is assigned,
    // (ii) the neighboring cluster, (iii) the farthest individual within the cluster,
    // (iv) the closest individual belonging to another cluster

    out_file << "\n\t" << STAT_label[STATL_WITHIN] << "-" << STAT_label[STATL_CLUSTER] << " " << STAT_label[STATL_DISTANCE]
             << "\t" << STAT_label[STATL_NEIGHBOR] << " " << STAT_label[STATL_CLUSTER] << "\t"
             << STAT_label[STATL_NEIGHBOR] << " " << STAT_label[STATL_CLUSTER] << " " << STAT_label[STATL_DISTANCE]
             << "\t" << STAT_label[STATL_MOST_DISTANT] << " " << distance_matrix->label << "\t"
             << STAT_label[STATL_MOST_DISTANT] << " " << distance_matrix->label  << " " << STAT_label[STATL_DISTANCE]
             << "\t" << STAT_label[STATL_NEIGHBOR] << " " << distance_matrix->label << "\t"
             << STAT_label[STATL_NEIGHBOR] << " " << distance_matrix->label  << " " << STAT_label[STATL_DISTANCE] << endl;

    for (i = 0;i < nb_pattern;i++) {
      neighbor_cluster = neighbor_pattern_cluster_selection(normalized_pattern_cluster_distance , order[i]);
      most_distant_pattern = most_distant_pattern_selection(normalized_pattern_distance , order[i]);
      neighbor_pattern = neighbor_pattern_selection(normalized_pattern_distance , order[i]);

      out_file << distance_matrix->label << " " << distance_matrix->row_identifier[order[i]]
               << "\t" << normalized_pattern_distance[order[i]][assignment[order[i]]]
               << "\t" << neighbor_cluster + 1 << "\t" << normalized_pattern_distance[order[i]][neighbor_cluster]
               << "\t" << distance_matrix->row_identifier[most_distant_pattern]
               << "\t" << normalized_pattern_distance[order[i]][most_distant_pattern]
               << "\t" << distance_matrix->row_identifier[neighbor_pattern]
               << "\t" << normalized_pattern_distance[order[i]][neighbor_pattern] << endl;
    }

    delete [] order;

    for (i = 0;i < nb_pattern;i++) {
      delete [] normalized_pattern_cluster_distance[i];
    }
    delete [] normalized_pattern_cluster_distance;

    // computation of within-cluster and between-cluster distances

    normalized_cluster_distance = new double*[nb_cluster];
    for (i = 0;i < nb_cluster;i++) {
      normalized_cluster_distance[i] = new double[nb_cluster];
      for (j = 0;j < nb_cluster;j++) {
        normalized_cluster_distance[i][j] = distance[i][j];
        if ((distance[i][j] != -D_INF) && (length[i][j] > 1)) {
          normalized_cluster_distance[i][j] /= length[i][j];
        }
      }
    }

    // writing of the between-cluster distance matrix

    out_file << "\n" << STAT_label[STATL_CLUSTER] << " " << STAT_label[STATL_DISTANCE] << " "
             << STAT_label[STATL_MATRIX] << endl;

    out_file << "\n";
    for (i = 0;i < nb_cluster;i++) {
      out_file << "\t" << STAT_label[STATL_CLUSTER] << " " << i + 1;
    }
    out_file << endl;

    for (i = 0;i < nb_cluster;i++) {
      out_file << STAT_label[STATL_CLUSTER] << " " << i + 1;
      for (j = 0;j < nb_cluster;j++) {
        out_file << "\t" << normalized_cluster_distance[i][j];
      }
      out_file << endl;
    }

    // for each cluster, writing of (i) the within-cluster distance, (ii) the between-cluster distance,
    // (iii) the largest distance between individuals of the cluster (diameter),
    // (iv) the smallest distance between individuals, one belonging to the cluster and
    // the other to another cluster (separation)

    out_file << "\n\t" << STAT_label[STATL_WITHIN] << "-" << STAT_label[STATL_CLUSTER] << " " << STAT_label[STATL_DISTANCE]
             << "\t"  << STAT_label[STATL_BETWEEN] << "-" << STAT_label[STATL_CLUSTER] << " " << STAT_label[STATL_DISTANCE]
             << "\t" << STAT_label[STATL_DIAMETER] << "\t" << STAT_label[STATL_SEPARATION] << endl;

    for (i = 0;i < nb_cluster;i++) {
      out_file << STAT_label[STATL_CLUSTER] << " " << i + 1
               << "\t" << normalized_cluster_distance[i][i]
               << "\t" << between_cluster_distance_computation(i)
               << "\t" << max_within_cluster_distance_computation(normalized_pattern_distance , i)
               << "\t" << min_between_cluster_distance_computation(normalized_pattern_distance , i) << endl;
    }

    for (i = 0;i < nb_cluster;i++) {
      out_file << "\n" << STAT_label[STATL_CLUSTER] << " " << i + 1 << "\t";
      isolation = isolation_property(normalized_pattern_distance , i , CLUSTER_SCALE);
      if (isolation) {
        out_file << STAT_label[STATL_ISOLATED];
      }
      else {
        isolation = isolation_property(normalized_pattern_distance , i , INDIVIDUAL);
        if (isolation) {
          out_file << STAT_label[STATL_ISOLATED] << " (" << STAT_label[STATL_PATTERN_LEVEL] << ")";
        }
        else {
          out_file << STAT_label[STATL_NON_ISOLATED];
        }
      }
    }
    out_file << endl;

    for (i = 0;i < nb_pattern;i++) {
      delete [] normalized_pattern_distance[i];
    }
    delete [] normalized_pattern_distance;

    for (i = 0;i < nb_cluster;i++) {
      delete [] normalized_cluster_distance[i];
    }
    delete [] normalized_cluster_distance;
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a Clusters object using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool Clusters::plot_write(StatError &error , const char *prefix ,
                          const char *title) const

{
  bool status = true;
  int i , j , k;
  int max_nb_pattern , plot_nb_cluster , *plot_nb_pattern , **index;
  double min_distance , max_distance , **normalized_distance;
  ostringstream *data_file_name;
  ofstream *out_data_file;


  // writing of the data files

  data_file_name = new ostringstream[nb_cluster];

  data_file_name[0] << prefix << 1 << ".dat";
  out_data_file = new ofstream((data_file_name[0].str()).c_str());

  if (!out_data_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  else {
    index = new int*[nb_cluster];
    normalized_distance = new double*[nb_cluster];
    plot_nb_pattern = new int[nb_cluster];

    max_nb_pattern = 0;

    for (i = 0;i < nb_cluster;i++) {

      // sort of distances by increasing distance

      index[i] = pattern_sort(i);

      plot_nb_pattern[i] = cluster_nb_pattern[i];

      normalized_distance[i] = new double[cluster_nb_pattern[i]];
      for (j = 0;j < cluster_nb_pattern[i];j++) {
        normalized_distance[i][j] = pattern_distance[index[i][j]][i];
        if (normalized_distance[i][j] == -D_INF) {
          plot_nb_pattern[i] = j;
          break;
        }
        else if (pattern_length[index[i][j]][i] > 1) {
          normalized_distance[i][j] /= pattern_length[index[i][j]][i];
        }
      }

      if (plot_nb_pattern[i] > max_nb_pattern) {
        max_nb_pattern = plot_nb_pattern[i];
      }
    }

    if (max_nb_pattern == 1) {
      status = false;
      error.update(STAT_error[STATR_SINGLE_ELEMENT_CLUSTERS]);
    }

    else {
      plot_nb_cluster = 0;
      i = 0;

      for (j = 0;j < nb_cluster;j++) {
        if (plot_nb_pattern[j] > 1) {
          plot_nb_cluster++;

          if (i == 0) {
            min_distance = normalized_distance[j][0];
            max_distance = normalized_distance[j][plot_nb_pattern[j] - 1];
          }

          else {
            if (normalized_distance[j][0] < min_distance) {
              min_distance = normalized_distance[j][0];
            }
            if (normalized_distance[j][plot_nb_pattern[j] - 1] > max_distance) {
              max_distance = normalized_distance[j][plot_nb_pattern[j] - 1];
            }

            data_file_name[i] << prefix << i + 1 << ".dat";
            out_data_file = new ofstream((data_file_name[i].str()).c_str());
          }

          for (k = 0;k < plot_nb_pattern[j];k++) {
            *out_data_file << k + 1 << " " << normalized_distance[j][k] << endl;
          }
          out_data_file->close();
          delete out_data_file;

          i++;
        }
      }

      // writing of the script files

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

        out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n";

        out_file << "set title" << " \"";
        if (title) {
          out_file << title << " - ";
        }
        out_file << nb_cluster << " " << STAT_label[STATL_CLUSTERS] << "\"\n\n";

        for (j = 0;j < nb_cluster;j++) {
          if (plot_nb_pattern[j] > 1) {
            for (k = 0;k < plot_nb_pattern[j];k++) {
              out_file << "set label \"" << distance_matrix->row_identifier[index[j][k]] << "\" at "
                       << k + 1 << ", " << normalized_distance[j][k] << endl;
            }
            out_file << endl;
          }
        }

        if (max_nb_pattern < TIC_THRESHOLD) {
          out_file << "set xtics 1,1" << endl;
        }

        out_file << "plot [1:" << max_nb_pattern << "] ["
                 << min_distance * (1. - PLOT_YMARGIN)  << ":"
                 << max_distance * (1. + PLOT_YMARGIN) << "] ";

        j = 0;
        for (k = 0;k < nb_cluster;k++) {
          if (plot_nb_pattern[k] > 1) {
            out_file << "\"" << stat_tool::label((data_file_name[j].str()).c_str())
                     << "\" using 1:2 title \"" << STAT_label[STATL_CLUSTER] << " " << k + 1
                     << "\" with linespoints";

            if (j < plot_nb_cluster - 1) {
              out_file << ",\\";
            }
            out_file << endl;

            j++;
          }
        }

        if (max_nb_pattern < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }

        out_file << "\nunset label" << endl;

        if (i == 1) {
          out_file << "\nset terminal x11" << endl;
        }

        out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
    }

    for (i = 0;i < nb_cluster;i++) {
      delete [] index[i];
      delete [] normalized_distance[i];
    }
    delete [] index;
    delete [] normalized_distance;

    delete [] plot_nb_pattern;

  }
  delete [] data_file_name;

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a Clusters object.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* Clusters::get_plotable(StatError &error) const

{
  int i , j , k;
  int max_nb_pattern , plot_nb_cluster , *plot_nb_pattern , **index;
  double min_distance , max_distance , **normalized_distance;
  ostringstream title , legend , identifier;
  MultiPlotSet *plot_set;


  index = new int*[nb_cluster];
  normalized_distance = new double*[nb_cluster];
  plot_nb_pattern = new int[nb_cluster];

  max_nb_pattern = 0;

  for (i = 0;i < nb_cluster;i++) {

    // sort of distances by increasing distances

    index[i] = pattern_sort(i);

    plot_nb_pattern[i] = cluster_nb_pattern[i];

    normalized_distance[i] = new double[cluster_nb_pattern[i]];
    for (j = 0;j < cluster_nb_pattern[i];j++) {
      normalized_distance[i][j] = pattern_distance[index[i][j]][i];
      if (normalized_distance[i][j] == -D_INF) {
        plot_nb_pattern[i] = j;
        break;
      }
      else if (pattern_length[index[i][j]][i] > 1) {
        normalized_distance[i][j] /= pattern_length[index[i][j]][i];
      }
    }

    if (plot_nb_pattern[i] > max_nb_pattern) {
      max_nb_pattern = plot_nb_pattern[i];
    }
  }

  if (max_nb_pattern == 1) {
    plot_set = NULL;
    error.update(STAT_error[STATR_SINGLE_ELEMENT_CLUSTERS]);
  }

  else {
    plot_set = new MultiPlotSet(1);
    MultiPlotSet &plot = *plot_set;

    title.str("");
    title << nb_cluster << " " << STAT_label[STATL_CLUSTERS];
    plot.title = title.str();

    plot.border = "15 lw 0";

    plot_nb_cluster = 0;
    i = 0;

    for (j = 0;j < nb_cluster;j++) {
      if (plot_nb_pattern[j] > 1) {
        plot_nb_cluster++;

        if (i == 0) {
          min_distance = normalized_distance[j][0];
          max_distance = normalized_distance[j][plot_nb_pattern[j] - 1];
        }

        else {
          if (normalized_distance[j][0] < min_distance) {
            min_distance = normalized_distance[j][0];
          }
          if (normalized_distance[j][plot_nb_pattern[j] - 1] > max_distance) {
            max_distance = normalized_distance[j][plot_nb_pattern[j] - 1];
          }
        }

        i++;
      }
    }

    plot[0].xrange = Range(1 , max_nb_pattern);
    plot[0].yrange = Range(min_distance * (1. - PLOT_YMARGIN) ,
                           max_distance * (1. + PLOT_YMARGIN));

    if (max_nb_pattern < TIC_THRESHOLD) {
      plot[0].xtics = 1;
    }

    plot[0].resize(plot_nb_cluster * 2);

    i = 0;
    for (j = 0;j < nb_cluster;j++) {
      if (plot_nb_pattern[j] > 1) {
        legend.str("");
        legend << STAT_label[STATL_CLUSTER] << " " << j + 1;
        plot[0][i * 2].legend = legend.str();

        plot[0][i * 2].style = "linespoints";

        plot[0][i * 2 + 1].label = "true";

        for (k = 0;k < plot_nb_pattern[j];k++) {
          plot[0][i * 2].add_point(k + 1 , normalized_distance[j][k]);

          identifier.str("");
          identifier << distance_matrix->row_identifier[index[j][k]];
          plot[0][i * 2 + 1].add_text(k + 1 , normalized_distance[j][k] , identifier.str());
        }

        i++;
      }
    }
  }

  for (i = 0;i < nb_cluster;i++) {
    delete [] index[i];
    delete [] normalized_distance[i];
  }
  delete [] index;
  delete [] normalized_distance;

  delete [] plot_nb_pattern;

  return plot_set;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a Clusters object.
 *
 *  \return MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* Clusters::get_plotable() const

{
  MultiPlotSet *plot_set;
  StatError error;


  return get_plotable(error);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of cluster sizes.
 */
/*--------------------------------------------------------------*/

void Clusters::cluster_nb_pattern_computation()

{
  int i;


  for (i = 0;i < nb_cluster;i++) {
    cluster_nb_pattern[i] = 0;
  }

  for (i = 0;i < nb_pattern;i++) {
    cluster_nb_pattern[assignment[i]]++;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of individual-cluster distances and associated lengths
 *         on the basis of individual assignment to clusters.
 */
/*--------------------------------------------------------------*/

void Clusters::pattern_distance_computation()

{
  int i , j;


  for (i = 0;i < nb_pattern;i++) {
    for (j = 0;j < nb_cluster;j++) {
      pattern_distance[i][j] = 0.;
      pattern_length[i][j] = 0;
    }
  }

  for (i = 0;i < nb_pattern;i++) {
    for (j = 0;j < nb_pattern;j++) {
      if (pattern_distance[i][assignment[j]] != -D_INF) {
        if (distance_matrix->distance[i][j] != -D_INF) {
          pattern_distance[i][assignment[j]] += distance_matrix->distance[i][j];
          pattern_length[i][assignment[j]] += distance_matrix->length[i][j];
        }
        else {
          pattern_distance[i][assignment[j]] = -D_INF;
          pattern_length[i][assignment[j]] = 0;
        }
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of between-cluster distances and associated lengths.
 */
/*--------------------------------------------------------------*/

void Clusters::cluster_distance_computation_1()

{
  int i , j;


  for (i = 0;i < nb_cluster;i++) {
    for (j = 0;j < nb_cluster;j++) {
      distance[i][j] = 0.;
      length[i][j] = 0;
    }
  }

  for (i = 0;i < distance_matrix->nb_row;i++) {
    for (j = 0;j < distance_matrix->nb_column;j++) {
      if (distance[assignment[i]][assignment[j]] != -D_INF) {
        if (distance_matrix->distance[i][j] != -D_INF) {
          distance[assignment[i]][assignment[j]] += distance_matrix->distance[i][j];
          length[assignment[i]][assignment[j]] += distance_matrix->length[i][j];
        }
        else {
          distance[assignment[i]][assignment[j]] = -D_INF;
          length[assignment[i]][assignment[j]] = 0;
        }
      }
    }
  }

  if ((distance_matrix->deletion_distance) && (distance_matrix->nb_deletion)) {
    for (i = 0;i < nb_cluster;i++) {
      for (j = 0;j < nb_cluster;j++) {
        deletion_distance[i][j] = 0.;
        nb_deletion[i][j] = 0;
      }
    }

    for (i = 0;i < distance_matrix->nb_row;i++) {
      for (j = 0;j < distance_matrix->nb_column;j++) {
        if (deletion_distance[assignment[i]][assignment[j]] != -D_INF) {
          if (distance_matrix->deletion_distance[i][j] != -D_INF) {
            deletion_distance[assignment[i]][assignment[j]] += distance_matrix->deletion_distance[i][j];
            nb_deletion[assignment[i]][assignment[j]] += distance_matrix->nb_deletion[i][j];
          }
          else {
            deletion_distance[assignment[i]][assignment[j]] = -D_INF;
            nb_deletion[assignment[i]][assignment[j]] = 0;
          }
        }
      }
    }
  }

  if ((distance_matrix->insertion_distance) && (distance_matrix->nb_insertion)) {
    for (i = 0;i < nb_cluster;i++) {
      for (j = 0;j < nb_cluster;j++) {
        insertion_distance[i][j] = 0.;
        nb_insertion[i][j] = 0;
      }
    }

    for (i = 0;i < distance_matrix->nb_row;i++) {
      for (j = 0;j < distance_matrix->nb_column;j++) {
        if (insertion_distance[assignment[i]][assignment[j]] != -D_INF) {
          if (distance_matrix->insertion_distance[i][j] != -D_INF) {
            insertion_distance[assignment[i]][assignment[j]] += distance_matrix->insertion_distance[i][j];
            nb_insertion[assignment[i]][assignment[j]] += distance_matrix->nb_insertion[i][j];
          }
          else {
            insertion_distance[assignment[i]][assignment[j]] = -D_INF;
            nb_insertion[assignment[i]][assignment[j]] = 0;
          }
        }
      }
    }
  }

  if (distance_matrix->nb_match) {
    for (i = 0;i < nb_cluster;i++) {
      for (j = 0;j < nb_cluster;j++) {
        nb_match[i][j] = 0;
      }
    }

    for (i = 0;i < distance_matrix->nb_row;i++) {
      for (j = 0;j < distance_matrix->nb_column;j++) {
        nb_match[assignment[i]][assignment[j]] += distance_matrix->nb_match[i][j];
      }
    }
  }

  if ((distance_matrix->substitution_distance) && (distance_matrix->nb_substitution)) {
    for (i = 0;i < nb_cluster;i++) {
      for (j = 0;j < nb_cluster;j++) {
        substitution_distance[i][j] = 0.;
        nb_substitution[i][j] = 0;
      }
    }

    for (i = 0;i < distance_matrix->nb_row;i++) {
      for (j = 0;j < distance_matrix->nb_column;j++) {
        if (substitution_distance[assignment[i]][assignment[j]] != -D_INF) {
          if (distance_matrix->substitution_distance[i][j] != -D_INF) {
            substitution_distance[assignment[i]][assignment[j]] += distance_matrix->substitution_distance[i][j];
            nb_substitution[assignment[i]][assignment[j]] += distance_matrix->nb_substitution[i][j];
          }
          else {
            substitution_distance[assignment[i]][assignment[j]] = -D_INF;
            nb_substitution[assignment[i]][assignment[j]] = 0;
          }
        }
      }
    }
  }

  if ((distance_matrix->transposition_distance) && (distance_matrix->nb_transposition)) {
    for (i = 0;i < nb_cluster;i++) {
      for (j = 0;j < nb_cluster;j++) {
        transposition_distance[i][j] = 0.;
        nb_transposition[i][j] = 0;
      }
    }

    for (i = 0;i < distance_matrix->nb_row;i++) {
      for (j = 0;j < distance_matrix->nb_column;j++) {
        if (transposition_distance[assignment[i]][assignment[j]] != -D_INF) {
          if (distance_matrix->transposition_distance[i][j] != -D_INF) {
            transposition_distance[assignment[i]][assignment[j]] += distance_matrix->transposition_distance[i][j];
            nb_transposition[assignment[i]][assignment[j]] += distance_matrix->nb_transposition[i][j];
          }
          else {
            transposition_distance[assignment[i]][assignment[j]] = -D_INF;
            nb_transposition[assignment[i]][assignment[j]] = 0;
          }
        }
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of between-cluster distances and associated lengths
 *         on the basis of individual-cluster distances and associated lengths.
 */
/*--------------------------------------------------------------*/

void Clusters::cluster_distance_computation_2()

{
  int i , j;


  for (i = 0;i < nb_cluster;i++) {
    for (j = 0;j < nb_cluster;j++) {
      distance[i][j] = 0.;
      length[i][j] = 0;
    }
  }

  for (i = 0;i < nb_pattern;i++) {
    for (j = 0;j < nb_cluster;j++) {
      if (distance[assignment[i]][j] != -D_INF) {
        if (pattern_distance[i][j] != -D_INF) {
          distance[assignment[i]][j] += pattern_distance[i][j];
          length[assignment[i]][j] += pattern_length[i][j];
        }
        else {
          distance[assignment[i]][j] = -D_INF;
          length[assignment[i]][j] = 0;
        }
      }
    }
  }

# ifdef DEBUG
  {
    int k;
    double **normalized_distance;


    normalized_distance = new double*[nb_cluster];
    for (i = 0;i < nb_cluster;i++) {
      normalized_distance[i] = new double[nb_cluster];
      for (j = 0;j < nb_cluster;j++) {
        normalized_distance[i][j] = distance[i][j];
        if ((distance[i][j] != -D_INF) && (length[i][j] > 1)) {
          normalized_distance[i][j] /= length[i][j];
        }
      }
    }

    for (i = 0;i < nb_cluster;i++) {
      for (j = i + 1;j < nb_cluster;j++) {
        for (k = 0;k < nb_cluster;k++) {
          if ((i != j) && (k != i) && (k != j) && (normalized_distance[i][k] != -D_INF) &&
              (normalized_distance[k][j] != -D_INF) && (normalized_distance[i][j] != -D_INF)) {
            if (normalized_distance[i][k] + normalized_distance[k][j] < normalized_distance[i][j]) {
              cout << STAT_label[STATL_CLUSTER] << " " << i + 1 << ", "
                   << STAT_label[STATL_CLUSTER] << " " << k + 1 << ", "
                   << STAT_label[STATL_CLUSTER] << " " << j + 1 << ": "
                   << STAT_error[STATR_TRIANGLE_INEQUALITY] << endl;
            }
          }
        }
      }
    }

    for (i = 0;i < nb_cluster;i++) {
      delete [] normalized_distance[i];
    }
    delete [] normalized_distance;
  }
# endif

}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of the farthest individual assigned to the same cluster.
 *
 *  \param[in] normalized_distance normalized distance matrix,
 *  \param[in] ipattern            individual index.
 *
 *  \return                        farthest individual index.
 */
/*--------------------------------------------------------------*/

int Clusters::most_distant_pattern_selection(double **normalized_distance , int ipattern) const

{
  int i;
  int pattern = ipattern;
  double max_within_cluster_distance;


  if (cluster_nb_pattern[assignment[ipattern]] > 1) {
    max_within_cluster_distance = 0.;
    for (i = 0;i < nb_pattern;i++) {
      if ((assignment[i] == assignment[ipattern]) &&
          (normalized_distance[ipattern][i] > max_within_cluster_distance)) {
        max_within_cluster_distance = normalized_distance[ipattern][i];
        pattern = i;
      }
    }
  }

  return pattern;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of the closest individual assigned to another cluster.
 *
 *  \param[in] normalized_distance normalized distance matrix,
 *  \param[in] ipattern            individual index.
 *
 *  \return                        closest individual index.
 */
/*--------------------------------------------------------------*/

int Clusters::neighbor_pattern_selection(double **normalized_distance , int ipattern) const

{
  int i;
  int pattern;
  double min_between_cluster_distance;


  min_between_cluster_distance = -D_INF;
  for (i = 0;i < nb_pattern;i++) {
    if ((assignment[i] != assignment[ipattern]) &&
        (normalized_distance[ipattern][i] < min_between_cluster_distance)) {
      min_between_cluster_distance = normalized_distance[ipattern][i];
      pattern = i;
    }
  }

  return pattern;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of the closest cluster for a given individual.
 *
 *  \param[in] normalized_distance normalized distance matrix,
 *  \param[in] pattern             individual index.
 *
 *  \return                        closest cluster index.
 */
/*--------------------------------------------------------------*/

int Clusters::neighbor_pattern_cluster_selection(double **normalized_distance , int pattern) const

{
  int i;
  int cluster;
  double min_between_cluster_distance;


  min_between_cluster_distance = -D_INF;
  for (i = 0;i < nb_cluster;i++) {
    if ((i != assignment[pattern]) &&
        (normalized_distance[pattern][i] < min_between_cluster_distance)) {
      min_between_cluster_distance = normalized_distance[pattern][i];
      cluster = i;
    }
  }

  return cluster;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the cluster diameter (maximum within-cluster distance).
 *
 *  \param[in] normalized_distance normalized distance matrix,
 *  \param[in] cluster             cluster index.
 *
 *  \return                        cluster diameter.
 */
/*--------------------------------------------------------------*/

double Clusters::max_within_cluster_distance_computation(double **normalized_distance , int cluster) const

{
  int i , j;
  double max_within_cluster_distance;


  max_within_cluster_distance = 0.;

  if (cluster_nb_pattern[cluster] > 1) {
    for (i = 0;i < nb_pattern - 1;i++) {
      if (assignment[i] == cluster) {
        for (j = i + 1;j < nb_pattern;j++) {
          if ((assignment[j] == cluster) &&
              (normalized_distance[i][j] > max_within_cluster_distance)) {
            max_within_cluster_distance = normalized_distance[i][j];
          }
        }
      }
    }
  }

  return max_within_cluster_distance;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the cluster separation (minimum between-cluster distance).
 *
 *  \param[in] normalized_distance normalized distance matrix,
 *  \param[in] cluster             cluster index.
 *
 *  \return                        cluster separation.
 */
/*--------------------------------------------------------------*/

double Clusters::min_between_cluster_distance_computation(double **normalized_distance , int cluster) const

{
  int i , j;
  double min_between_cluster_distance;


  min_between_cluster_distance = -D_INF;

  for (i = 0;i < nb_pattern;i++) {
    if (assignment[i] == cluster) {
      for (j = 0;j < nb_pattern;j++) {
        if ((assignment[j] != cluster) &&
            (normalized_distance[i][j] < min_between_cluster_distance)) {
          min_between_cluster_distance = normalized_distance[i][j];
        }
      }
    }
  }

  return min_between_cluster_distance;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the cluster isolation property.
 *
 *  \param[in] normalized_distance normalized distance matrix,
 *  \param[in] cluster             cluster index,
 *  \param[in] scale               property at the individual (INDIVIDUAL) or cluster (CLUSTER_SCALE) level.
 *
 *  \return                        isolation property or not.
 */
/*--------------------------------------------------------------*/

bool Clusters::isolation_property(double **normalized_distance , int cluster ,
                                  isolation_scale scale) const

{
  bool isolation = true;
  int i , j;
  double max_within_distance , min_between_distance;


  if (scale == CLUSTER_SCALE) {
    max_within_distance = 0.;
    min_between_distance = -D_INF;
  }

  for (i = 0;i < nb_pattern;i++) {
    if (assignment[i] == cluster) {
      if (scale == INDIVIDUAL) {
        max_within_distance = 0.;
        min_between_distance = -D_INF;
      }

      for (j = 0;j < nb_pattern;j++) {
        if (j != i) {
          if (assignment[j] == cluster) {
            if (normalized_distance[i][j] > max_within_distance) {
              max_within_distance = normalized_distance[i][j];
            }
          }

          else {
            if (normalized_distance[i][j] < min_between_distance) {
              min_between_distance = normalized_distance[i][j];
            }
          }
        }
      }

      if ((scale == INDIVIDUAL) && (max_within_distance >= min_between_distance)) {
        isolation = false;
        break;
      }
    }
  }

  if ((scale == CLUSTER_SCALE) && (max_within_distance >= min_between_distance)) {
    isolation = false;
  }

  return isolation;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the average distance between a given cluster and
 *         all the other clusters.
 *
 *  \param[in] cluster cluster index.
 *
 *  \return            average distance.
 */
/*--------------------------------------------------------------*/

double Clusters::between_cluster_distance_computation(int cluster) const

{
  int i;
  int between_length;
  double between_distance;


  between_distance = 0.;
  between_length = 0;

  for (i = 0;i < nb_cluster;i++) {
    if (i != cluster) {
      if (distance[cluster][i] == -D_INF) {
        between_distance = -D_INF;
        break;
      }
      else {
        between_distance += distance[cluster][i];
        between_length += length[cluster][i];
      }
    }
  }

  if ((between_distance != -D_INF) && (between_length > 1)) {
    between_distance /= between_length;
  }

  return between_distance;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of within-cluster and between-cluster distances and
 *         writing of the results.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& Clusters::global_distance_ascii_print(ostream &os)

{
  int i , j;
  int within_length , between_length;
  double within_distance , between_distance;


  within_distance = 0.;
  within_length = 0;
  between_distance = 0.;
  between_length = 0;

  for (i = 0;i < nb_cluster;i++) {
    if (within_distance != -D_INF) {
      if (distance[i][i] == -D_INF) {
        within_distance = -D_INF;
      }
      else {
        within_distance += distance[i][i];
        within_length += length[i][i];
      }
    }

    if (between_distance != -D_INF) {
      for (j = 0;j < nb_cluster;j++) {
        if (j != i) {
          if (distance[i][j] == -D_INF) {
            between_distance = -D_INF;
            break;
          }
          else {
            between_distance += distance[i][j];
            between_length += length[i][j];
          }
        }
      }
    }
  }

  if ((within_distance != -D_INF) && (within_length > 1)) {
    within_distance /= within_length;
  }
  if ((between_distance != -D_INF) && (between_length > 1)) {
    between_distance /= between_length;
  }

  os << "\n" << STAT_label[STATL_WITHIN] << "-" << STAT_label[STATL_CLUSTER] << " " << STAT_label[STATL_DISTANCE]
     << ": " << within_distance << "   " << STAT_label[STATL_BETWEEN] << "-" << STAT_label[STATL_CLUSTER]
     << " " << STAT_label[STATL_DISTANCE] << ": " << between_distance;
  if ((within_distance != -D_INF) && (within_length > 0) && (between_distance != -D_INF) && (between_length > 0)) {
    os << "   " << STAT_label[STATL_RATIO] << ": " << within_distance / between_distance << "\n";
  }
  os << endl;

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Initialization of prototypes as a function of the average distance of an individual to
 *         the other individuals.
 */
/*--------------------------------------------------------------*/

void Clusters::prototype_initialization_1()

{
  int i;
  int *index;
  double *cumul_distance;


  cumul_distance = new double[distance_matrix->nb_row];

  for (i = 0;i < distance_matrix->nb_row;i++) {
    cumul_distance[i] = stat_tool::cumul_distance_computation(distance_matrix->nb_column , distance_matrix->distance[i]);
    if (cumul_distance[i] != -D_INF) {
      cumul_distance[i] /= cumul_computation(1 , distance_matrix->nb_column , distance_matrix->length + i);
    }
  }

  index = stat_tool::pattern_sort(distance_matrix->nb_row , cumul_distance);

  for (i = 0;i < nb_cluster;i++) {
    assignment[index[(int)round((double)(i * distance_matrix->nb_row) / (double)nb_cluster)]] = i;
  }

  delete [] cumul_distance;
  delete [] index;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Initialization of prototypes (Kaufman & Rousseeuw, pp. 102-103).
 */
/*--------------------------------------------------------------*/

void Clusters::prototype_initialization_2()

{
  int i , j , k , m;
  int *prototype , *index;
  double distance_diff , min_distance , max_distance , *cumul_distance , **normalized_distance;


  normalized_distance = new double*[distance_matrix->nb_row];
  for (i = 0;i < distance_matrix->nb_row;i++) {
    normalized_distance[i] = new double[distance_matrix->nb_column];
    for (j = 0;j < distance_matrix->nb_column;j++) {
      normalized_distance[i][j] = distance_matrix->distance[i][j];
      if ((distance_matrix->distance[i][j] != -D_INF) && (distance_matrix->length[i][j] > 1)) {
        normalized_distance[i][j] /= distance_matrix->length[i][j];
      }
    }
  }

  prototype = new int[nb_cluster];

  // determination of the first prototype (most central individual)

  cumul_distance = new double[distance_matrix->nb_row];

  for (i = 0;i < distance_matrix->nb_row;i++) {
    cumul_distance[i] = stat_tool::cumul_distance_computation(distance_matrix->nb_column , distance_matrix->distance[i]);
    if (cumul_distance[i] != -D_INF) {
      cumul_distance[i] /= cumul_computation(1 , distance_matrix->nb_column , distance_matrix->length + i);
    }
  }

  index = stat_tool::pattern_sort(distance_matrix->nb_row , cumul_distance , 1);

  prototype[0] = index[0];
  assignment[prototype[0]] = 0;

  delete [] cumul_distance;
  delete [] index;

  // determination of the other prototypes

  for (i = 1;i < nb_cluster;i++) {
    max_distance = -1.;
    for (j = 0;j < nb_pattern;j++) {
      if (assignment[j] == I_DEFAULT) {
        distance_diff = 0.;
        for (k = 0;k < nb_pattern;k++) {
          if (assignment[k] == I_DEFAULT) {
            min_distance = -D_INF * 10.;
            for (m = 0;m < i;m++) {
              if (normalized_distance[k][prototype[m]] < min_distance) {
                min_distance = normalized_distance[k][prototype[m]];
              }
            }

            if (min_distance > normalized_distance[k][j]) {
              distance_diff += min_distance - normalized_distance[k][j];
            }
          }
        }

        if (distance_diff > max_distance) {
          max_distance = distance_diff;
          prototype[i] = j;
        }
      }
    }

    assignment[prototype[i]] = i;
  }

  for (i = 0;i < distance_matrix->nb_row;i++) {
    delete [] normalized_distance[i];
  }
  delete [] normalized_distance;

  delete [] prototype;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Update of the distances between individuals of the old and new cluster and
 *         the associated lengths as a function of the change in assignment of an individual.
 *
 *  \param[in] pattern     individual index,
 *  \param[in] old_cluster old cluster index.
 *  \param[in] new_cluster new cluster index.
 */
/*--------------------------------------------------------------*/

void Clusters::pattern_distance_update(int pattern , int old_cluster , int new_cluster)

{
  bool *pattern_flag , *cluster_flag;
  int i , j;


  for (i = 0;i < nb_pattern;i++) {
    if (pattern_distance[i][old_cluster] != -D_INF) {
      pattern_distance[i][old_cluster] -= distance_matrix->distance[i][pattern];
      pattern_length[i][old_cluster] -= distance_matrix->length[i][pattern];
    }

    else {
      if (distance_matrix->distance[i][pattern] == -D_INF) {
        pattern_flag = new bool[nb_pattern];
        cluster_flag = new bool[nb_pattern];

        for (j = 0;j < nb_pattern;j++) {
          pattern_flag[j] = false;
          if (assignment[j] == old_cluster) {
            cluster_flag[j] = true;
          }
          else {
            cluster_flag[j] = false;
          }
        }

        pattern_flag[i] = true;
        cluster_flag[pattern] = false;

        pattern_distance[i][old_cluster] = distance_matrix->cumul_distance_computation(pattern_flag , cluster_flag);

        if (pattern_distance[i][old_cluster] != -D_INF) {
          pattern_length[i][old_cluster] = distance_matrix->cumul_length_computation(pattern_flag , cluster_flag);
        }
        else {
          pattern_length[i][old_cluster] = 0;
        }

        delete [] pattern_flag;
        delete [] cluster_flag;
      }
    }

    if (pattern_distance[i][new_cluster] != -D_INF) {
      if (distance_matrix->distance[i][pattern] != -D_INF) {
        pattern_distance[i][new_cluster] += distance_matrix->distance[i][pattern];
        pattern_length[i][new_cluster] += distance_matrix->length[i][pattern];
      }
      else {
        pattern_distance[i][new_cluster] = -D_INF;
        pattern_length[i][new_cluster] = 0;
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Step of a partitioning clustering algorithm.
 */
/*--------------------------------------------------------------*/

void Clusters::algorithmic_step_1()

{
  bool infinite_distance;
  int i , j;
  int iter , nb_change , *tmp_assignment;
  double min_distance , normalized_distance;


  infinite_distance = false;
  for (i = 0;i < distance_matrix->nb_row;i++) {
    for (j = 0;j < distance_matrix->nb_column;j++) {
      if (distance_matrix->distance[i][j] == -D_INF) {
        infinite_distance = true;
        break;
      }
    }
    if (infinite_distance) {
      break;
    }
  }

  tmp_assignment = new int[nb_pattern];

  // initial assignment of the individuals

  for (i = 0;i < nb_pattern;i++) {
    if (assignment[i] == I_DEFAULT) {
      min_distance = -D_INF * 10.;
      for (j = 0;j < nb_pattern;j++) {
        if (assignment[j] != I_DEFAULT) {
          normalized_distance = distance_matrix->distance[i][j];
          if ((distance_matrix->distance[i][j] != -D_INF) && (distance_matrix->length[i][j] > 1)) {
            normalized_distance /= distance_matrix->length[i][j];
          }

          if (normalized_distance < min_distance) {
            min_distance = normalized_distance;
            tmp_assignment[i] = assignment[j];
          }
        }
      }
    }
  }

  for (i = 0;i < nb_pattern;i++) {
    if (assignment[i] == I_DEFAULT) {
      assignment[i] = tmp_assignment[i];
    }
  }

  cluster_nb_pattern_computation();
  pattern_distance_computation();

# ifdef DEBUG
  cout << "\n";
  for (i = 0;i < nb_cluster;i++) {
    cout << STAT_label[STATL_CLUSTER] << " " << i + 1 << ":";
    for (j = 0;j < nb_pattern;j++) {
      if (assignment[j] == i) {
        cout << " " << distance_matrix->row_identifier[j];
      }
    }
    cout << endl;
  }
# endif

  iter = 0;
  do {
    iter++;
    nb_change = 0;

    // assignment of the individuals

    for (i = 0;i < nb_pattern;i++) {
      min_distance = -D_INF * 10.;

      for (j = 0;j < nb_cluster;j++) {
        normalized_distance = pattern_distance[i][j];
        if ((pattern_distance[i][j] != -D_INF) && (pattern_length[i][j] > 1)) {
          normalized_distance /= pattern_length[i][j];
        }

        if (normalized_distance < min_distance) {
          min_distance = normalized_distance;
          tmp_assignment[i] = j;
        }
      }

      // test old cluster empty and update of the size of the old and new clusters as a function
      // of the change in assignment of an individual

      if (tmp_assignment[i] != assignment[i]) {
        if (cluster_nb_pattern[assignment[i]] > 1) {
          cluster_nb_pattern[assignment[i]]--;
          cluster_nb_pattern[tmp_assignment[i]]++;
        }
        else {
          tmp_assignment[i] = assignment[i];
        }
      }

      // update of the distances between individuals of the old and new clusters
      // as a function of the change in assignment of an individual

      if ((iter > GLOBAL_NB_ITER) && (tmp_assignment[i] != assignment[i])) {
        nb_change++;
        pattern_distance_update(i , assignment[i] , tmp_assignment[i]);
        assignment[i] = tmp_assignment[i];
      }
    }

    // global update of the individual-cluster distances

    if (iter <= GLOBAL_NB_ITER) {
      for (i = 0;i < nb_pattern;i++) {
        if (tmp_assignment[i] != assignment[i]) {
          nb_change++;
          if (!infinite_distance) {
            pattern_distance_update(i , assignment[i] , tmp_assignment[i]);
          }
          assignment[i] = tmp_assignment[i];
        }
      }

      if (infinite_distance) {
        pattern_distance_computation();
      }
    }

#   ifdef DEBUG
    {
      cout << "\niteration " << iter << ": " << nb_change << endl;
      for (i = 0;i < nb_cluster;i++) {
        cout << STAT_label[STATL_CLUSTER] << " " << i + 1 << ":";
        for (j = 0;j < nb_pattern;j++) {
          if (assignment[j] == i) {
            cout << " " << distance_matrix->row_identifier[j];
          }
        }
        cout << endl;
      }

      cluster_distance_computation_2();
      global_distance_ascii_print(cout);
    }
#   endif

  }
  while ((iter < PARTITIONING_NB_ITER_1) && (nb_change > 0));

# ifdef DEBUG
  cout << "\n" << iter << " iterations" << endl;
# endif

  delete [] tmp_assignment;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Step of a partitioning clustering algorithm (Kaufman & Rousseeuw, pp. 103-104).
 */
/*--------------------------------------------------------------*/

void Clusters::algorithmic_step_2()

{
  int i , j , k , m;
  int iter , pattern , cluster , *prototype , *tmp_assignment , *index;
  double min_distance , min_swap_distance = D_INF , previous_min_swap_distance , *prototype_distance ,
         **normalized_distance , **swap_distance;


  // initializations

  normalized_distance = new double*[distance_matrix->nb_row];
  for (i = 0;i < distance_matrix->nb_row;i++) {
    normalized_distance[i] = new double[distance_matrix->nb_column];
    for (j = 0;j < distance_matrix->nb_column;j++) {
      normalized_distance[i][j] = distance_matrix->distance[i][j];
      if ((distance_matrix->distance[i][j] != -D_INF) && (distance_matrix->length[i][j] > 1)) {
        normalized_distance[i][j] /= distance_matrix->length[i][j];
      }
    }
  }

  prototype = new int[nb_cluster];
  tmp_assignment = new int[nb_pattern];

  for (i = 0;i < nb_pattern;i++) {
    if (assignment[i] != I_DEFAULT) {
      prototype[assignment[i]] = i;
    }
    tmp_assignment[i] = assignment[i];
  }

  // assignment of the non-selected individuals to clusters

  for (i = 0;i < nb_pattern;i++) {
    if (tmp_assignment[i] == I_DEFAULT) {
      min_distance = -D_INF * 10;
      for (j = 0;j < nb_cluster;j++) {
        if (normalized_distance[i][prototype[j]] < min_distance) {
          min_distance = normalized_distance[i][prototype[j]];
          assignment[i] = j;
        }
      }
    }
  }

  swap_distance = new double*[nb_pattern];
  for (i = 0;i < nb_pattern;i++) {
    swap_distance[i] = new double[nb_cluster];
  }

  prototype_distance = new double[nb_cluster + 1];

# ifdef DEBUG
  cout << "\n";
  for (i = 0;i < nb_cluster;i++) {
    cout << STAT_label[STATL_CLUSTER] << " " << i + 1 << " - " << distance_matrix->row_identifier[prototype[i]] << ":";
    for (j = 0;j < nb_pattern;j++) {
      if ((tmp_assignment[j] == I_DEFAULT) && (assignment[j] == i)) {
        cout << " " << distance_matrix->row_identifier[j];
      }
    }
    cout << endl;
  }
# endif

  iter = 0;
  do {
    iter++;
    previous_min_swap_distance = min_swap_distance;

    // computation of the costs of replacement of prototypes

    for (i = 0;i < nb_pattern;i++) {
      if (tmp_assignment[i] == I_DEFAULT) {
        for (j = 0;j < nb_cluster;j++) {
          swap_distance[i][j] = 0.;
          for (k = 0;k < nb_pattern;k++) {
            if (tmp_assignment[k] == I_DEFAULT) {
              for (m = 0;m < nb_cluster;m++) {
                prototype_distance[m] = normalized_distance[k][prototype[m]];
              }
              prototype_distance[nb_cluster] = normalized_distance[k][i];

              index = stat_tool::pattern_sort(nb_cluster + 1 , prototype_distance , 2);

              if (((index[0] == j) && (index[1] == nb_cluster)) || ((index[0] == nb_cluster) && (index[1] == j))) {
                swap_distance[i][j] += normalized_distance[k][i] - normalized_distance[k][prototype[j]];
              }
              else if (index[0] == j) {
                swap_distance[i][j] += normalized_distance[k][prototype[index[1]]] - normalized_distance[k][prototype[j]];
              }
              else if (index[0] == nb_cluster) {
                swap_distance[i][j] += normalized_distance[k][i] - normalized_distance[k][prototype[index[1]]];
              }

              delete [] index;
            }
          }
        }
      }
    }

    // choice of replacement of a prototype

    min_swap_distance = -D_INF * nb_pattern * nb_cluster;
    for (i = 0;i < nb_pattern;i++) {
      if (tmp_assignment[i] == I_DEFAULT) {
        for (j = 0;j < nb_cluster;j++) {
          if (swap_distance[i][j] < min_swap_distance) {
            min_swap_distance = swap_distance[i][j];
            pattern = i;
            cluster = j;
          }
        }
      }
    }

#   ifdef MESSAGE
    cout << "\nminimum swap distance : " << min_swap_distance
         << "   old " << STAT_label[STATL_PROTOTYPE] << ": " << distance_matrix->row_identifier[prototype[cluster]]
         << "   new " << STAT_label[STATL_PROTOTYPE] << ": " << distance_matrix->row_identifier[pattern] << endl;
#   endif

    if ((min_swap_distance < 0.) && (min_swap_distance > previous_min_swap_distance)) {
      tmp_assignment[prototype[cluster]] = I_DEFAULT;
      tmp_assignment[pattern] = cluster;
      assignment[pattern] = cluster;
      prototype[cluster] = pattern;

      // assignment of the non-selected individuals to clusters

      for (i = 0;i < nb_pattern;i++) {
        if (tmp_assignment[i] == I_DEFAULT) {
          min_distance = -D_INF * 10.;
          for (j = 0;j < nb_cluster;j++) {
            if (normalized_distance[i][prototype[j]] < min_distance) {
              min_distance = normalized_distance[i][prototype[j]];
              assignment[i] = j;
            }
          }
        }
      }

      cluster_nb_pattern_computation();
      pattern_distance_computation();

#     ifdef MESSAGE
      {
        cout << "\niteration " << iter << endl;
        for (i = 0;i < nb_cluster;i++) {
          cout << STAT_label[STATL_CLUSTER] << " " << i + 1 << " - " << distance_matrix->row_identifier[prototype[i]] << ":";
          for (j = 0;j < nb_pattern;j++) {
            if ((tmp_assignment[j] == I_DEFAULT) && (assignment[j] == i)) {
              cout << " " << distance_matrix->row_identifier[j];
            }
          }
          cout << endl;
        }

        cluster_distance_computation_2();
        global_distance_ascii_print(cout);
      }
#     endif

    }

    else {
      break;
    }
  }
  while (iter < PARTITIONING_NB_ITER_2);

# ifdef MESSAGE
  cout << "\n" << iter << " iterations" << endl;
# endif

  for (i = 0;i < distance_matrix->nb_row;i++) {
    delete [] normalized_distance[i];
  }
  delete [] normalized_distance;

  delete [] prototype;
  delete [] tmp_assignment;

  for (i = 0;i < nb_pattern;i++) {
    delete [] swap_distance[i];
  }
  delete [] swap_distance;

  delete [] prototype_distance;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Partitioning clustering algorithm.
 *
 *  \param[in] error          reference on a StatError object,
 *  \param[in] display        flag for displaying the partitioning results,
 *  \param[in] nb_cluster     number of clusters,
 *  \param[in] initialization initial prototypes,
 *  \param[in] algorithm      algorithm type.
 *
 *  \return                   Clusters object.
 */
/*--------------------------------------------------------------*/

Clusters* DistanceMatrix::partitioning(StatError &error , bool display , int nb_cluster ,
                                       int *prototype , int initialization , int algorithm) const

{
  bool status = true , *selected_prototype;
  int i , j;
  int max_identifier;
  DistanceMatrix *dist_matrix;
  Clusters *clusters;


  clusters = NULL;
  error.init();

  if (nb_row != nb_column) {
    status = false;
    error.correction_update(STAT_error[STATR_MATRIX_STRUCTURE] , STAT_error[STATR_SQUARE_MATRIX]);
  }
  else if (nb_row <= nb_cluster) {
    status = false;
    error.update(STAT_error[STATR_MATRIX_DIMENSIONS]);
  }

  if (nb_cluster < 2) {
    status = false;
    error.update(STAT_error[STATR_NB_CLUSTER]);
  }
  if (strcmp(label , STAT_label[STATL_CLUSTER]) == 0) {
    status = false;
    error.update(STAT_error[STATR_PATTERN_TYPE]);
  }

  if ((status) && (prototype)) {
    max_identifier = 0;
    for (i = 0;i < nb_cluster;i++) {
      if (prototype[i] > max_identifier) {
        max_identifier = prototype[i];
      }
    }

    selected_prototype = new bool[max_identifier + 1];
    for (i = 0;i <= max_identifier;i++) {
      selected_prototype[i] = false;
    }

    for (i = 0;i < nb_cluster;i++) {
      for (j = 0;j < nb_row;j++) {
        if (prototype[i] == row_identifier[j]) {
          break;
        }
      }

      if (j == nb_row) {
        status = false;
        ostringstream error_message;
        error_message << prototype[i] << ": " << STAT_error[STATR_PROTOTYPE_IDENTIFIER];
        error.update((error_message.str()).c_str());
      }

      else if (selected_prototype[prototype[i]]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_PROTOTYPE] << " " << prototype[i] << " "
                      << STAT_error[STATR_ALREADY_SELECTED];
        error.update((error_message.str()).c_str());
      }
      else {
        selected_prototype[prototype[i]] = true;
      }
    }

    delete [] selected_prototype;
  }

  if (status) {
    dist_matrix = new DistanceMatrix(*this , (test_symmetry() ? COPY : SYMMETRIZATION));

#   ifdef DEBUG
    if (!(test_symmetry())) {
      cout << *dist_matrix << endl;
    }
#   endif

    clusters = new Clusters(*dist_matrix , nb_cluster);

    // initialization of cluster prototypes

    if (prototype) {
      for (i = 0;i < clusters->nb_cluster;i++) {
        for (j = 0;j < clusters->nb_pattern;j++) {
          if (prototype[i] == clusters->distance_matrix->row_identifier[j]) {
            clusters->assignment[j] = i;
            break;
          }
        }
      }
    }

    else {
      switch (initialization) {
      case 1 :
        clusters->prototype_initialization_1();
        break;
      case 2 :
        clusters->prototype_initialization_2();
        break;
      }

#     ifdef MESSAGE
      for (i = 0;i < clusters->nb_cluster;i++) {
        for (j = 0;j < clusters->nb_pattern;j++) {
          if (clusters->assignment[j] == i) {
            cout << "\n" << STAT_label[STATL_CLUSTER] << " " << i + 1 << ": " << clusters->distance_matrix->row_identifier[j];
          }
        }
      }
      cout << endl;
#     endif
    }

    // iterative clustering algorithm

    switch (algorithm) {
    case 1 :
      clusters->algorithmic_step_1();
      break;
    case 2 :
      clusters->algorithmic_step_2();
      break;
    }

    clusters->cluster_distance_computation_1();

    if (display) {
      clusters->global_distance_ascii_print(cout);
    }

    delete dist_matrix;
  }

  return clusters;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Partitioning clustering algorithm.
 *
 *  \param[in] error              reference on a StatError object,
 *  \param[in] display            flag for displaying the partitioning results,
 *  \param[in] nb_cluster         number of clusters,
 *  \param[in] cluster_nb_pattern cluster sizes.
 *  \param[in] cluster_pattern    cluster compositions.
 *
 *  \return                       Clusters object.
 */
/*--------------------------------------------------------------*/

Clusters* DistanceMatrix::partitioning(StatError &error , bool display , int nb_cluster ,
                                       int *cluster_nb_pattern , int **cluster_pattern) const

{
  bool status = true , *selected_pattern;
  int i , j , k;
  int nb_pattern , max_identifier;
  DistanceMatrix *dist_matrix;
  Clusters *clusters;


  clusters = NULL;
  error.init();

  if (nb_row != nb_column) {
    status = false;
    error.correction_update(STAT_error[STATR_MATRIX_STRUCTURE] , STAT_error[STATR_SQUARE_MATRIX]);
  }
  else if (nb_row <= nb_cluster) {
    status = false;
    error.update(STAT_error[STATR_MATRIX_DIMENSIONS]);
  }

  if (nb_cluster < 2) {
    status = false;
    error.update(STAT_error[STATR_NB_CLUSTER]);
  }

  if (status) {
    nb_pattern = 0;
    for (i = 0;i < nb_cluster;i++) {
      nb_pattern += cluster_nb_pattern[i];
    }

    if (nb_pattern != nb_row) {
      status = false;
      ostringstream error_message;
      error_message << STAT_error[STATR_NUMBER] << " " << label << "s";
      error.update((error_message.str()).c_str());
    }

    max_identifier = 0;
    for (i = 0;i < nb_cluster;i++) {
      for (j = 0;j < cluster_nb_pattern[i];j++) {
        if (cluster_pattern[i][j] > max_identifier) {
          max_identifier = cluster_pattern[i][j];
        }
      }
    }

    selected_pattern = new bool[max_identifier + 1];
    for (i = 0;i <= max_identifier;i++) {
      selected_pattern[i] = false;
    }

    for (i = 0;i < nb_cluster;i++) {
      for (j = 0;j < cluster_nb_pattern[i];j++) {
        for (k = 0;k < nb_row;k++) {
          if (cluster_pattern[i][j] == row_identifier[k]) {
            break;
          }
        }

        if (k == nb_row) {
          status = false;
          ostringstream error_message;
          error_message << cluster_pattern[i][j] << ": " << STAT_error[STATR_BAD] << " "
                        << label << " " << STAT_label[STATL_IDENTIFIER];
          error.update((error_message.str()).c_str());
        }

        else if (selected_pattern[cluster_pattern[i][j]]) {
          status = false;
          ostringstream error_message;
          error_message << label << " " << cluster_pattern[i][j] << " "
                        << STAT_error[STATR_ALREADY_SELECTED];
          error.update((error_message.str()).c_str());
        }
        else {
          selected_pattern[cluster_pattern[i][j]] = true;
        }
      }
    }

    delete [] selected_pattern;
  }

  if (status) {
    dist_matrix = new DistanceMatrix(*this , (test_symmetry() ?  COPY : SYMMETRIZATION));

#   ifdef DEBUG
    if (!(test_symmetry())) {
      cout << *dist_matrix << endl;
    }
#   endif

    clusters = new Clusters(*dist_matrix , nb_cluster , cluster_nb_pattern , cluster_pattern);

    clusters->pattern_distance_computation();
    clusters->cluster_distance_computation_1();

    if (display) {
      clusters->global_distance_ascii_print(cout);
    }

    delete dist_matrix;
  }

  return clusters;
}


};  // namespace stat_tool
