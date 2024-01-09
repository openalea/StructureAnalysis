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
 *        $Id: export_distance_matrix.cpp 18032 2015-04-23 07:13:59Z guedon $
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
using namespace stat_tool;


// DistanceMatrix

class DistanceMatrixWrap
{

public:

  WRAP_METHOD_FILE_ASCII_WRITE( DistanceMatrix);

  WRAP_METHOD0(DistanceMatrix, symmetrize, DistanceMatrix);

  WRAP_METHOD0(DistanceMatrix, unnormalize, DistanceMatrix);

  static Clusters*
  partitioning_prototype(const DistanceMatrix& dm, int nb_cluster,
      const boost::python::list& prototype, int initialization, int algorithm)
  {
    Clusters *ret;
    StatError error;

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
  partitioning_clusters(const DistanceMatrix& dm,
      const boost::python::list& clusters)
  {
    Clusters *ret;
    StatError error;

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
  hierarchical_clustering(const DistanceMatrix& dm, int algorithm,
      int criterion, const char*path, char format)
  {
    StatError error;
    std::stringstream output;
    bool ret;

    ret = dm.hierarchical_clustering(error, output, algorithm, criterion, path,
        format);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return output.str();

  }

  static DistanceMatrix*
  read_from_file(char *filename)
  {
    DistanceMatrix * data;
    StatError error;

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

  static DistanceMatrix*
  select_individual(const DistanceMatrix& v,
      const boost::python::list& identifiers, bool keep)
  {
    StatError error;
    DistanceMatrix * ret = NULL;

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
  get_length(DistanceMatrix &input, int i, int j)
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
  get_deletion_distance(DistanceMatrix &input, int i, int j)
  {
    double ret;
    bool bret=false;
    int row_max = input.get_nb_row();
    int column_max = input.get_nb_column();

    ostringstream error_message;

    bret = input.is_deletion();
    if (bret == false)
      {
        error_message << "! no deletion data computed in this object.";
        cout << error_message.str() << endl;
        return -1;
    }

    if (i < row_max && j < column_max && i >= 0 && j >= 0)
      ret = input.get_deletion_distance(i, j);
    else
      {
        error_message << "index not in valid range" << "i must be less than "
        << row_max << " and j less than " << column_max << endl;
        cout << error_message.str() << endl;
        ret = -1;
      }
    return ret;
  }

  static double
  get_substitution_distance(DistanceMatrix &input, int i, int j)
  {
    double ret;
    bool bret=false;
    int row_max = input.get_nb_row();
    int column_max = input.get_nb_column();

    ostringstream error_message;

    bret = input.is_substitution();
    if (bret == false)
      {
        error_message << "! no substitution data computed in this object.";
        cout << error_message.str() << endl;
        return -1;
      }
    
    if (i < row_max && j < column_max && i >= 0 && j >= 0)
      ret = input.get_substitution_distance(i, j);
    else
      {
        error_message << "index not in valid range" << "i must be less than "
        << row_max << " and j less than " << column_max << endl;
        cout << error_message.str() << endl;
        ret = -1;
      }
    return ret;
  }
  
  
  
  static double
  get_transposition_distance(DistanceMatrix &input, int i, int j)
  {
    double ret;
    bool bret=false;
    int row_max = input.get_nb_row();
    int column_max = input.get_nb_column();

    ostringstream error_message;

    bret = input.is_transposition();
    if (bret == false)
      {
        error_message << "! no transposition data computed in this object.";
        cout << error_message.str() << endl;
        return -1;
    }

    if (i < row_max && j < column_max && i >= 0 && j >= 0)
      ret = input.get_transposition_distance(i, j);
    else
      {
        error_message << "index not in valid range" << "i must be less than "
        << row_max << " and j less than " << column_max << endl;
        cout << error_message.str() << endl;
        ret = -1;
      }
    return ret;
  }
  
  
  static double
  get_insertion_distance(DistanceMatrix &input, int i, int j)
  {
    double ret;
    bool bret=false;
    int row_max = input.get_nb_row();
    int column_max = input.get_nb_column();

    ostringstream error_message;

    bret = input.is_insertion();
    if (bret == false)
      {
        error_message << "! no insertion data computed in this object.";
        cout << error_message.str() << endl;
        return -1;
    }
    if (i < row_max && j < column_max && i >= 0 && j >= 0)
      ret = input.get_insertion_distance(i, j);
    else
      {
        error_message << "index not in valid range" << "i must be less than "
            << row_max << " and j less than " << column_max << endl;
        cout << error_message.str() << endl;
        ret = -1;
      }
    return ret;
  }
  
  static double
  get_distance(DistanceMatrix &input, int i, int j)
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

  static int
  get_nb_deletion(DistanceMatrix &input, int i, int j)
  {
    int ret;
    bool bret = false;
    int row_max = input.get_nb_row();
    int column_max = input.get_nb_column();

    ostringstream error_message;

    bret = input.is_deletion();
    if (bret == false)
      {
        error_message << "! no deletion data computed in this object.";
        cout << error_message.str() << endl;
        return -1;
      }
    if (i < row_max && j < column_max && i >= 0 && j >= 0)
          ret = input.get_nb_deletion(i, j);
    else
      {
        error_message << "index not in valid range" << "i must be less than "
            << row_max << " and j less than " << column_max << endl;
        cout << error_message.str() << endl;
        ret = -1;
      }
    return ret;
  }
  static int
  get_nb_insertion(DistanceMatrix &input, int i, int j)
  {
    int ret;
    bool bret = false;
    int row_max = input.get_nb_row();
    int column_max = input.get_nb_column();

    ostringstream error_message;

    bret = input.is_insertion();
    if (bret == false)
      {
        error_message << "! no insertion data computed in this object.";
        cout << error_message.str() << endl;
        return -1;
      }

    if (i < row_max && j < column_max && i >= 0 && j >= 0)
          ret = input.get_nb_insertion(i, j);
    else
      {
        error_message << "index not in valid range" << "i must be less than "
            << row_max << " and j less than " << column_max << endl;
        cout << error_message.str() << endl;
        ret = -1;
      }
    return ret;
  }
  static int
  get_nb_transposition(DistanceMatrix &input, int i, int j)
  {
    int ret;
    bool bret = false;
    int row_max = input.get_nb_row();
    int column_max = input.get_nb_column();

    ostringstream error_message;

    bret = input.is_transposition();
    if (bret == false)
      {
        error_message << "! no transposition data computed in this object.";
        cout << error_message.str() << endl;
        return -1;
      }

    if (i < row_max && j < column_max && i >= 0 && j >= 0)
          ret = input.get_nb_transposition(i, j);
    else
      {
        error_message << "index not in valid range" << "i must be less than "
            << row_max << " and j less than " << column_max << endl;
        cout << error_message.str() << endl;
        ret = -1;
      }
    return ret;
  }
  static int
  get_nb_substitution(DistanceMatrix &input, int i, int j)
  {
    int ret;
    bool bret = false;
    int row_max = input.get_nb_row();
    int column_max = input.get_nb_column();

    ostringstream error_message;

    bret = input.is_substitution();
    if (bret == false)
      {
        error_message << "! no substitution data computed in this object.";
        cout << error_message.str() << endl;
        return -1;
      }

    if (i < row_max && j < column_max && i >= 0 && j >= 0)
          ret = input.get_nb_substitution(i, j);
    else
      {
        error_message << "index not in valid range" << "i must be less than "
            << row_max << " and j less than " << column_max << endl;
        cout << error_message.str() << endl;
        ret = -1;
      }
    return ret;
  }

  static int
  get_nb_match(DistanceMatrix &input, int i, int j)
  {
    int ret;
    bool bret = false;
    int row_max = input.get_nb_row();
    int column_max = input.get_nb_column();

    ostringstream error_message;

    bret = input.is_match();
    if (bret == false)
      {
        error_message << "! no match data computed in this object.";
        cout << error_message.str() << endl;
        return -1;
      }

    if (i < row_max && j < column_max && i >= 0 && j >= 0)
          ret = input.get_nb_match(i, j);
    else
      {
        error_message << "index not in valid range" << "i must be less than "
            << row_max << " and j less than " << column_max << endl;
        cout << error_message.str() << endl;
        ret = -1;
      }
    return ret;
  }



};

#define WRAP DistanceMatrixWrap
#define CLASS DistanceMatrix

void
class_distance_matrix()
{


  class_< DistanceMatrix, bases< StatInterface> >
  ("_DistanceMatrix", "Distance Matrix", init< const DistanceMatrix&>())

  DEF_INIT_MAKE_CONSTRUCTOR(WRAP::read_from_file, "todo")

  .def(self_ns::str(self)) // __str__

  .add_property("nb_row", &CLASS::get_nb_row, "get number of rows")
  .add_property("nb_column", &CLASS::get_nb_column, "get number of columns")

  .def("test_symmetry", &CLASS::test_symmetry, "returns True if symmetric")
  // test the validity of the arguments by using a wrapped function
  .def("get_distance", &WRAP::get_distance,
      args("irow", "icolumn"), "returns distance between element i,j where i in [0, nbrow] and j in [0,nbcolum]")
  .def("get_length", WRAP::get_length,
      args("irow", "icolumn"), "returns length between element i,j where i in [0, nbrow] and j in [0,nbcolum]")
  .def("get_transposition_distance", WRAP::get_transposition_distance,
      args("irow", "icolumn"), "returns transposition between element i,j where i in [0, nbrow] and j in [0,nbcolum]")
  .def("get_insertion_distance", WRAP::get_insertion_distance,
      args("irow", "icolumn"), "returns insertion between element i,j where i in [0, nbrow] and j in [0,nbcolum]")
  .def("get_deletion_distance", WRAP::get_deletion_distance,
      args("irow", "icolumn"), "returns deletion between element i,j where i in [0, nbrow] and j in [0,nbcolum]")
  .def("get_substitution_distance", WRAP::get_substitution_distance,
      args("irow", "icolumn"), "returns substitution between element i,j where i in [0, nbrow] and j in [0,nbcolum]")

  .def("get_row_identifier", &CLASS::get_row_identifier,
      args("index"), "todo")
  .def("get_column_identifier", &CLASS::get_column_identifier,
      args("index"), "todo")
  
  .def("get_nb_substitution", WRAP::get_nb_substitution,
      args("irow", "icolumn"), "returns nb of substitution between element i,j where i in [0, nbrow] and j in [0,nbcolum]")
  .def("get_nb_deletion", WRAP::get_nb_deletion,
      args("irow", "icolumn"), "returns nb of deletion between element i,j where i in [0, nbrow] and j in [0,nbcolum]")
  .def("get_nb_transposition", WRAP::get_nb_transposition,
      args("irow", "icolumn"), "returns nb of transposition between element i,j where i in [0, nbrow] and j in [0,nbcolum]")
  .def("get_nb_match", WRAP::get_nb_match,
      args("irow", "icolumn"), "returns nb of match between element i,j where i in [0, nbrow] and j in [0,nbcolum]")
  .def("get_nb_insertion", WRAP::get_nb_insertion,
      args("irow", "icolumn"), "returns nb of insertion between element i,j where i in [0, nbrow] and j in [0,nbcolum]")

  // Clustering
  DEF_RETURN_VALUE_NO_ARGS("partitioning_prototype", WRAP::partitioning_prototype,
      "to be done")
  DEF_RETURN_VALUE_NO_ARGS("partitioning_clusters", WRAP::partitioning_clusters,
      "to be done")
  .def("hierarchical_clustering", WRAP::hierarchical_clustering,
      args("algorithm","criterion","path","format"),"Clustering using hierarchical methods")
  DEF_RETURN_VALUE_NO_ARGS("symmetrize", WRAP::symmetrize,
      "symmetrize distance matrix")
  DEF_RETURN_VALUE_NO_ARGS("unnormalize", WRAP::unnormalize,
      "symmetrize distance matrix")
  DEF_RETURN_VALUE_NO_ARGS("file_ascii_write", WRAP::file_ascii_write,
      "Save vector summary into a file")
  DEF_RETURN_VALUE("select_individual", WRAP::select_individual,
      args("identifiers", "keep"),"Select individuals given a list of identifiers")

;

/*
 void update(...
 void update(int irow_identifier , int icolumn_identifier , double idistance , int ilength);

 int get_label_size() const { return label_size; }     char* get_label() const { return label; }

 */

}
#undef WRAP
#undef CLASS

void
class_cluster()
{
  class_<Clusters, bases<DistanceMatrix> > ("_Cluster", "Cluster", no_init)
  .def(init <const DistanceMatrix &, int>())

  .def(self_ns::str(self)) // __str__

  ;

  /*
  Clusters();
     Clusters(const DistanceMatrix &dist_matrix , int inb_cluster ,
              int *icluster_nb_pattern , int **cluster_pattern);
     Clusters(const Clusters &clusters):DistanceMatrix(clusters) { copy(clusters); }


     std::ostream& line_write(std::ostream &os) const;

     std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
     bool ascii_write(StatError &error , const char *path , bool exhaustive = false) const;
     bool spreadsheet_write(StatError &error , const char *path) const;
     bool plot_write(StatError &error , const char *prefix ,
                     const char *title = 0) const;
     MultiPlotSet* get_plotable() const;


     void cluster_nb_pattern_computation();
     void pattern_distance_computation();
     void cluster_distance_computation_1();
     void cluster_distance_computation_2();

     // acces membres de la classe

     DistanceMatrix* get_distance_matrix() { return distance_matrix; }
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
  class_<Dendrogram, bases<StatInterface> > ("_Dendrogram", "Dendrogram", no_init)
  .def(init <const DistanceMatrix &, int>())

  .def(self_ns::str(self)) // __str__
;
}

/*
Dendrogram();
   Dendrogram(const Dendrogram &dendrogram) { copy(dendrogram); }

   std::ostream& line_write(std::ostream &os) const;

   std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
   bool ascii_write(StatError &error , const char *path , bool exhaustive = false) const;
   bool spreadsheet_write(StatError &error , const char *path) const;
   bool plot_write(StatError &error , const char *prefix ,
                   const char *title = 0) const { return false; }
   // acces membres de la classe

   DistanceMatrix* get_distance_matrix() { return distance_matrix; }
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
