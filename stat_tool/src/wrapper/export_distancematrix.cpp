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

using namespace boost::python;


// DistanceMatrix

class DistanceMatrixWrap
{

public:
  static Clusters* partitioning_prototype(const Distance_matrix& dm, int nb_cluster, 
					  const boost::python::list& prototype, 
					  int initialization, int algorithm)
  {
    Clusters *ret;
    Format_error error;
    bool berror = false;

    std::stringstream output;
    int nb_proto = len(prototype);

    ostringstream error_message;
    try{
      if (algorithm!=1 and algorithm!=2)
      {
        error_message << "Incorrect value for the algorithm argument (must be 1 or 2)"<<endl;
        PyErr_SetString(PyExc_TypeError, (error_message.str()).c_str());
        berror = true;
      }
    }
    catch(...)
    {
      berror = true;
    }   


    stat_tool::wrap_util::auto_ptr_array<int>
      protos(new int[nb_proto]);

    for (int i=0; i<nb_proto; i++)
      protos[i] = extract<int>(prototype[i]);

    ret = dm.partitioning(error, output, nb_cluster,
			  protos.get(), initialization, algorithm);

    if(!ret)
      stat_tool::wrap_util::throw_error(error);


    return ret;
  }


  static Clusters* partitioning_clusters(const Distance_matrix& dm, 
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
    try{
      for (int i=0; i<nb_cluster; i++)
	{
	  boost::python::list* l = extract<boost::python::list*>(clusters[i]);
	  int nb_item = len(*l);
	  cluster_nb_pattern[i] = nb_item;
	  
	  cluster_pattern[i] = new int[nb_item];

	  for (int j=0; j<cluster_nb_pattern[i]; j++)
	    cluster_pattern[i][j] = extract<int>((*l)[j]);
	}
    }
    catch(...)
      {
	// Free memory
	for (int i=0; i<nb_cluster; i++)
	  delete[] cluster_pattern[i];
	
	delete[] cluster_nb_pattern;
	delete[] cluster_pattern;
      }
    
    ret = dm.partitioning(error, output, nb_cluster,
			  cluster_nb_pattern, cluster_pattern);
      

    // Free memory
    for (int i=0; i<nb_cluster; i++)
      delete[] cluster_pattern[i];

    delete[] cluster_nb_pattern;
    delete[] cluster_pattern;

    if(!ret)
      stat_tool::wrap_util::throw_error(error);


    return ret;
  }


  static std::string hierarchical_clustering(const Distance_matrix& dm,
					     int algorithm, int criterion,
					     const string &path, char format)
  {
    Format_error error;
    std::stringstream output;
    bool ret;
    
    ret = dm.hierarchical_clustering(error , output,
				     algorithm, criterion,
				     path.c_str(), format);
    
    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return output.str();
    
  }



};



void class_distance_matrix()
{
  enum_<stat_tool::wrap_util::UniqueInt<3, 0> >("AlgoType")
      .value("AGGLOMERATIVE", AGGLOMERATIVE)
      .value("DIVISIVE", DIVISIVE)
      .value("ORDERING", ORDERING)
      .export_values()
      ;

    enum_<stat_tool::wrap_util::UniqueInt<2, 0> >("CriterionType")
      .value("FARTHEST_NEIGHBOR", FARTHEST_NEIGHBOR)
      .value("AVERAGING", AVERAGING)
      .export_values()
      ;

  class_< Distance_matrix, bases< STAT_interface > >
    ("_DistanceMatrix", "Distance Matrix", init<const Distance_matrix&>())


     .def(self_ns::str(self)) // __str__
    
    // Clustering
    .def("partitioning_prototype", DistanceMatrixWrap::partitioning_prototype,
	 return_value_policy< manage_new_object >() )
    .def("partitioning_clusters", DistanceMatrixWrap::partitioning_clusters,
	 return_value_policy< manage_new_object >() )
    .def("hierarchical_clustering", DistanceMatrixWrap::hierarchical_clustering, "Clustering using hierarchical methods")
    ;
}


void class_clusters()
{
  class_< Clusters, bases< Distance_matrix > >
    ("_Clusters", "Cluster")
    ;
}


void class_dendrogram()
{
  class_< Dendrogram, bases< STAT_interface > >
    ("_Dendrogram", "Dendrogram")
    ;
}


