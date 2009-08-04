/*------------------------------------------------------------------------------
 *
 *        VPlants.Stat_Tool : VPlants Statistics module
 *
 *        Copyright 2006-2007 INRIA - CIRAD - INRA
 *
 *        File author(s): Yann Guédon <yann.guedon@cirad.fr>
 *                        Jean-Baptiste Durand <Jean-Baptiste.Durand@imag.fr>
 *                        Samuel Dufour-Kowalski <samuel.dufour@sophia.inria.fr>
 *                        Christophe Pradal <christophe.prada@cirad.fr>
 *                        Thomas Cokelaer <Thomas.Cokelaer@inria.fr>
 *
 *        Distributed under the GPL 2.0 License.
 *        See accompanying file LICENSE.txt or copy at
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr
 *
 *        $Id$
 *
 *-----------------------------------------------------------------------------*/

#include "wrapper_util.h"
#include "export_base.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/distance_matrix.h"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/make_constructor.hpp>

#include "boost_python_aliases.h"

using namespace boost::python;

// DistanceMatrix

class DistanceMatrixWrap
{

public:

  WRAP_METHOD_FILE_ASCII_WRITE( Distance_matrix);

  WRAP_METHOD0(Distance_matrix, symmetrize, Distance_matrix);

  WRAP_METHOD0(Distance_matrix, unnormalize, Distance_matrix);

  static Clusters*
  partitioning_prototype(const Distance_matrix& dm, int nb_cluster,
      const boost::python::list& prototype, int initialization, int algorithm)
  {
    Clusters *ret;
    Format_error error;
    bool berror = false;

    std::stringstream output;
    int nb_proto = len(prototype);

    ostringstream error_message;

    if (nb_proto != 0)
      {
        stat_tool::wrap_util::auto_ptr_array<int> protos(new int[nb_proto]);
        for (int i = 0; i < nb_proto; i++)
          protos[i] = extract<int> (prototype[i]);
        ret = dm.partitioning(error, output, nb_cluster, protos.get(),
            initialization, algorithm);
      }
    else
      {
        int *protos = 0;
        ret = dm.partitioning(error, output, nb_cluster, protos,
            initialization, algorithm);
      }

    cout << output.str() << endl;
    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Clusters*
  partitioning_clusters(const Distance_matrix& dm,
      const boost::python::list& clusters)
  {
    Clusters *ret;
    Format_error error;

    std::stringstream output;
    int nb_cluster = len(clusters);
    int* cluster_nb_pattern;
    int** cluster_pattern;

    cluster_nb_pattern = new int[nb_cluster];
    cluster_pattern = new int*[nb_cluster];

    // Build dynamic 2D array
    try
      {
        for (int i = 0; i < nb_cluster; i++)
          {
            boost::python::list* l =
                extract<boost::python::list*> (clusters[i]);
            int nb_item = len(*l);
            cluster_nb_pattern[i] = nb_item;

            cluster_pattern[i] = new int[nb_item];

            for (int j = 0; j < cluster_nb_pattern[i]; j++)
              {
                cluster_pattern[i][j] = extract<int> ((*l)[j]);
              }
          }
      }
    catch (...)
      {
        // Free memory
        for (int i = 0; i < nb_cluster; i++)
          delete[] cluster_pattern[i];

        delete[] cluster_nb_pattern;
        delete[] cluster_pattern;
      }

    ret = dm.partitioning(error, output, nb_cluster, cluster_nb_pattern,
        cluster_pattern);

    // Free memory
    for (int i = 0; i < nb_cluster; i++)
      delete[] cluster_pattern[i];

    delete[] cluster_nb_pattern;
    delete[] cluster_pattern;

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static std::string
  hierarchical_clustering(const Distance_matrix& dm, int algorithm,
      int criterion, const char*path, char format)
  {
    Format_error error;
    std::stringstream output;
    bool ret;

    ret = dm.hierarchical_clustering(error, output, algorithm, criterion, path,
        format);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return output.str();

  }

  static Distance_matrix*
  read_from_file(char *filename)
  {
    Distance_matrix * data;
    Format_error error;

    // not implemented. see vectors.cpp for an exmaple.
    //data = distance_matrix_ascii_read(error, filename);

    if (data)
      {
        return data;
      }
    else
      {
        stat_tool::wrap_util::throw_error(error);
      }

  }
  ;

  static Distance_matrix*
  select_individual(const Distance_matrix& v,
      const boost::python::list& identifiers, bool keep)
  {
    Format_error error;
    Distance_matrix * ret = NULL;

    int nb_id = len(identifiers);

    stat_tool::wrap_util::auto_ptr_array<int> ids(new int[nb_id]);

    for (int i = 0; i < nb_id; i++)
      ids[i] = extract<int> (identifiers[i]);

    ret = v.select_individual(error, nb_id, ids.get(), keep);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static double
  get_length(Distance_matrix &input, int i, int j)
  {
    double ret;
    int row_max = input.get_nb_row();
    int column_max = input.get_nb_column();

    ostringstream error_message;
    error_message << "index not in valid range" << endl;

    CHECK(i, 0, row_max);
    CHECK(j, 0, column_max);

    ret = input.get_length(i, j);

    return ret;
  }

  static double
  get_distance(Distance_matrix &input, int i, int j)
  {
    double ret;
    int row_max = input.get_nb_row();
    int column_max = input.get_nb_column();

    ostringstream error_message;
    error_message << "index not in valid range" << "i must be less than "
        << row_max << " and j less than " << column_max << endl;

    if (i < row_max && j < column_max && i >= 0 && j >= 0)
      ret = input.get_distance(i, j);
    else
      {
        cout << error_message.str() << endl;
        ret = -1;
      }
    return ret;
  }

};

#define WRAP DistanceMatrixWrap
#define CLASS Distance_matrix

void
class_distance_matrix()
{


  class_< Distance_matrix, bases< STAT_interface> >
  ("_DistanceMatrix", "Distance Matrix", init< const Distance_matrix&>())

  DEF_INIT_MAKE_CONSTRUCTOR(WRAP::read_from_file, "todo")
  .def(self_ns::str(self)) // __str__
  .def("get_nb_row", &CLASS::get_nb_row, "get number of rows")
  .def("get_nb_column", &CLASS::get_nb_column, "get number of columns")
  .def("test_symmetry", &CLASS::test_symmetry, "returns True if symmetric")
  // test the validity of the arguments by using a wrapped function
  .def("get_distance", &WRAP::get_distance, ARGS("irow", "icolumn"), "todo")
  .def("get_length", WRAP::get_length,ARGS("irow", "icolumn"), "todo")
  .def("get_row_identifier", &CLASS::get_row_identifier,ARGS("index"), "todo")
  .def("get_column_identifier", &CLASS::get_column_identifier,ARGS("index"), "todo")

  // Clustering
  DEF_RETURN_VALUE_NO_ARGS("partitioning_prototype", WRAP::partitioning_prototype,"to be done")
  DEF_RETURN_VALUE_NO_ARGS("partitioning_clusters", WRAP::partitioning_clusters, "to be done")
  .def("hierarchical_clustering", WRAP::hierarchical_clustering, ARGS("algorithm","criterion","path","format"),"Clustering using hierarchical methods")
  DEF_RETURN_VALUE_NO_ARGS("symmetrize", WRAP::symmetrize, "symmetrize distance matrix")
  DEF_RETURN_VALUE_NO_ARGS("unnormalize", WRAP::unnormalize, "symmetrize distance matrix")
  DEF_RETURN_VALUE_NO_ARGS("file_ascii_write", WRAP::file_ascii_write, "Save vector summary into a file")
  DEF_RETURN_VALUE("select_individual", WRAP::select_individual, ARGS("identifiers", "keep"),"Select individuals given a list of identifiers")

;

/*
 Distance_matrix();
 Distance_matrix(int nb_pattern , const char *ilabel , int *pattern_identifier = 0);
 Distance_matrix(int nb_pattern , int irow_identifier , int icolumn_identifier , const char *ilabel , int *pattern_identifier = 0 ,                     bool substitution_flag = true , bool transposition_flag = false);
 Distance_matrix(const Distance_matrix &dist_matrix , int inb_pattern , int *iidentifier , bool keep = true);
 Distance_matrix(const Distance_matrix &dist_matrix , int nb_cluster , const char *ilabel);
 Distance_matrix(const Distance_matrix &dist_matrix , char transform = 'c')     { copy(dist_matrix , transform); }

 void update(...
 void update(int irow_identifier , int icolumn_identifier , double idistance , int ilength);

 double get_deletion_distance(int row , int column) const     { return deletion_distance[row][column]; }
 int get_nb_deletion(int row , int column) const     { return nb_deletion[row][column]; }
 double get_insertion_distance(int row , int column) const     { return insertion_distance[row][column]; }
 int get_nb_insertion(int row , int column) const     { return nb_insertion[row][column]; }

 int get_nb_match(int row , int column) const     { return nb_match[row][column]; }

 double get_substitution_distance(int row , int column) const     { return substitution_distance[row][column]; }
 int get_nb_substitution(int row , int clumn) const     { return nb_substitution[row][column]; }

 double get_transposition_distance(int row , int column) const     { return transposition_distance[row][column]; }
 int get_nb_transposition(int row , int column) const     { return nb_transposition[row][column]; }

 int get_label_size() const { return label_size; }     char* get_label() const { return label; }

 */

}
#undef WRAP
#undef CLASS

void
class_cluster()
{
  class_<Clusters, bases<Distance_matrix> > ("_Cluster", "Cluster", no_init)
  .def(init <const Distance_matrix &, int>())

  .def(self_ns::str(self)) // __str__

  ;

  /*
  Clusters();
     Clusters(const Distance_matrix &dist_matrix , int inb_cluster ,
              int *icluster_nb_pattern , int **cluster_pattern);
     Clusters(const Clusters &clusters):Distance_matrix(clusters) { copy(clusters); }


     std::ostream& line_write(std::ostream &os) const;

     std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
     bool ascii_write(Format_error &error , const char *path , bool exhaustive = false) const;
     bool spreadsheet_write(Format_error &error , const char *path) const;
     bool plot_write(Format_error &error , const char *prefix ,
                     const char *title = 0) const;
     MultiPlotSet* get_plotable() const;


     void cluster_nb_pattern_computation();
     void pattern_distance_computation();
     void cluster_distance_computation_1();
     void cluster_distance_computation_2();

     // acces membres de la classe

     Distance_matrix* get_distance_matrix() { return distance_matrix; }
     int get_nb_pattern() const { return nb_pattern; }
     int get_nb_cluster() const { return nb_cluster; }
     int get_cluster_nb_pattern(int cluster) const { return cluster_nb_pattern[cluster]; }
     int get_assignment(int pattern) const { return assignment[pattern]; }
     double get_pattern_distance(int pattern , int cluster) const
     { return pattern_distance[pattern][cluster]; }
     int get_pattern_length(int pattern , int cluster) const
     { return pattern_length[pattern][cluster]; }
     */
}

void
class_dendrogram()
{
  class_<Dendrogram, bases<STAT_interface> > ("_Dendrogram", "Dendrogram", no_init)
  .def(init <const Distance_matrix &, int>())

  .def(self_ns::str(self)) // __str__
;
}

/*
Dendrogram();
   Dendrogram(const Dendrogram &dendrogram) { copy(dendrogram); }

   std::ostream& line_write(std::ostream &os) const;

   std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
   bool ascii_write(Format_error &error , const char *path , bool exhaustive = false) const;
   bool spreadsheet_write(Format_error &error , const char *path) const;
   bool plot_write(Format_error &error , const char *prefix ,
                   const char *title = 0) const { return false; }
   // acces membres de la classe

   Distance_matrix* get_distance_matrix() { return distance_matrix; }
   int get_scale() const { return scale; }
   int get_nb_cluster() const { return nb_cluster; }
   int get_cluster_nb_pattern(int cluster) const { return cluster_nb_pattern[cluster]; }
   int get_cluster_pattern(int cluster , int index) const { return cluster_pattern[cluster][index]; }
   int get_parent(int cluster) const { return parent[cluster]; }
   int get_child(int cluster , int index) const { return child[cluster][index]; }
   double get_child_distance(int cluster) const { return child_distance[cluster]; }
   double get_intra_cluster_distance(int cluster) const { return intra_cluster_distance[cluster]; }
   double get_inter_cluster_distance(int cluster) const { return inter_cluster_distance[cluster]; }
   double get_max_intra_cluster_distance(int cluster) const { return max_intra_cluster_distance[cluster]; }
   double get_min_inter_cluster_distance(int cluster) const { return min_inter_cluster_distance[cluster]; }
*/
