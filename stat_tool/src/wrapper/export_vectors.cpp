/*------------------------------------------------------------------------------
 *
 *        VPlants.Stat_Tool : VPlants Statistics module
 *
 *        Copyright 2006-2007 INRIA - CIRAD - INRA
 *
 *        File author(s): Yann Guédon <yann.guedon@cirad.fr>
 *                        Jean-Baptiste Dur&& <Jean-Baptiste.Dur&&@imag.fr>
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

#include "stat_tool/vectors.h"
#include "stat_tool/regression.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/distribution.h"
#include "stat_tool/mv_mixture.h"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/make_constructor.hpp>

using namespace boost::python;
using namespace boost;




// Vectors

class VectorsWrap
{

public:
  static Vectors* read_from_file(char *filename)
  {
    Vectors *vec;
    Format_error error;

    vec = vectors_ascii_read(error, filename);

    if(vec) { return vec; }
    else { stat_tool::wrap_util::throw_error(error);   }

  };

  static Vectors* build_from_lists (boost::python::list& array)
  {
    int nb_vector = boost::python::len(array);
    int *identifier = 0;
    int nb_variable = -1;
    int **int_vector = 0;
    double **float_vector = 0;
    bool is_float = false;
    bool error = false;
    Vectors* vec = 0;

    try {
    // For each vector
      for(int vi=0; vi < nb_vector; vi++)
    {
      boost::python::list vec = boost::python::extract<boost::python::list>(array[vi]);

      // nb_variable is the length of the first vector
      if(vi == 0)
        {
          nb_variable = boost::python::len(vec);

          // Check the type
          boost::python::extract<int> get_int(vec[0]);
          if (get_int.check())  // Array of int
          is_float = false;
          else
        {
          boost::python::extract<double> get_double(vec[0]);
          if (get_double.check()) // Array of float
              is_float = true;
          else
            {
              // Error : type is neither double, nor int
              ostringstream error_message;
              error_message << "Incorrect type (expect int or float values)"<<vi<<endl;
              PyErr_SetString(PyExc_TypeError, (error_message.str()).c_str());
              error = true;
            }
        }
        }
      // Check the size of the other vectors
      else if(boost::python::len(vec) != nb_variable)
        {
          // Error : Bad size
          ostringstream error_message;
          error_message << "Incorrect size of vector "<<vi<<endl;
          PyErr_SetString(PyExc_ValueError, (error_message.str()).c_str());
          error = true;
        }


      // Build array if necessary
      if(!int_vector && !is_float)
        {
          int_vector = new int*[nb_vector];
          for(int e=0; e < nb_vector; e++)
        {
          int_vector[e] = new int[nb_variable];
        }

        }
      else if (!float_vector && is_float)
        {
          float_vector = new double*[nb_vector];
          for(int e=0; e < nb_vector; e++)
        {
          float_vector[e] = new double[nb_variable];
        }
        }

      // Extract each element of the vector
      for(int i=0; i<nb_variable; i++)
        {
          if(!is_float)
        {
          int v = boost::python::extract<int>(vec[i]);
          int_vector[vi][i] = v;
        }
          else
        {
          double v = boost::python::extract<double>(vec[i]);
          float_vector[vi][i] = v;
        }

        }
    } // end of for
    }
    catch(...)
      {
    error = true;
      }

    // Call constructor
    if(!error){
      if(!is_float) vec = new Vectors(nb_vector, identifier, nb_variable, int_vector);
      else vec = new Vectors(nb_vector, identifier, nb_variable, float_vector);
    }

    // Delete memory
    if (int_vector) {
      for (int i = 0; i < nb_vector; i++)
    {
      delete [] int_vector[i];
    }
      delete [] int_vector;
    }

    if (float_vector) {
      for (int i = 0; i < nb_vector; i++)
    {
      delete [] float_vector[i];
    }
      delete [] float_vector;
    }

    if(error)
      {
    throw_error_already_set();
      }

    return vec;

  };


  static boost::python::list get_item(const Vectors* vec, int index)
  {
    // Test index
    if(index<0 || index>=vec->get_nb_vector())
      {
    PyErr_SetString(PyExc_IndexError, "vector index out of bound");
    boost::python::throw_error_already_set();
      }

    boost::python::list l;

    int  nb_var = vec->get_nb_variable();
    for(int var=0; var<nb_var; var++)
      {
    if((vec->get_type(var) == INT_VALUE) || (vec->get_type(var) == STATE))
      l.append(vec->get_int_vector(index, var));
    else
      l.append(vec->get_real_vector(index, var));
      }

    return l;
  }

  static boost::python::list get_identifiers(const Vectors& vec)
  {
    boost::python::list l;

    int  nb_vec = vec.get_nb_vector();
    for(int v=0; v<nb_vec; v++)
    {
      l.append(vec.get_identifier(v));
    }

    return l;
  }

  static std::string ascii_data_write(const Vectors& d, bool exhaustive)
  {
    std::stringstream s;
    std::string res;


    d.ascii_data_write(s, exhaustive, true);
    res = s.str();
    return res;
  }

  static void file_ascii_write(const Vectors& d, const char* path, bool exhaustive)
  {
    bool result = true;
    Format_error error;

    result = d.ascii_write(error, path,exhaustive);
    if (!result)
       stat_tool::wrap_util::throw_error(error);

  }

 static void file_ascii_data_write(const Vectors& d, const char* path, bool exhaustive)
  {
    bool result = true;
    Format_error error;

    result = d.ascii_data_write(error, path,exhaustive);
    if (!result)
       stat_tool::wrap_util::throw_error(error);

  }

  static Vectors* value_select(const Vectors& v,  int variable,
                   const object& min, const object& max, bool keep)
  {
    Format_error error;
    Vectors * ret = NULL;

    boost::python::extract<int> get_min(min);
    boost::python::extract<int> get_max(max);

    if (get_min.check() && get_max.check())  // Array of int
      {
    int mi = get_min();
    int ma = get_max();
    ret = v.value_select(error, variable, mi, ma, keep);
      }
    else
      {
    double mi = extract<double>(min);
    double ma = extract<double>(max);
    ret = v.value_select(error, variable, mi, ma, keep);
      }

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }


  static Vectors* select_variable(const Vectors& v,
				  const boost::python::list& variables,
                  bool keep)
  {
    Format_error error;
    Vectors * ret = NULL;

    int nb_var = len(variables);
		stat_tool::wrap_util::auto_ptr_array<int> vars(new int[nb_var]);

    for (int i=0; i<nb_var; i++)
    	vars[i] = extract<int>(variables[i]);

    ret = v.select_variable(error, nb_var, vars.get(), keep);

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }


  static Vectors* select_individual(const Vectors& v,
		  const boost::python::list& identifiers,
		  bool keep)
  {
    Format_error error;
    Vectors * ret = NULL;

    int nb_id = len(identifiers);
    stat_tool::wrap_util::auto_ptr_array<int> ids(new int[nb_id]);

    for (int i=0; i<nb_id; i++)
    	ids[i] = extract<int>(identifiers[i]);

    ret = v.select_individual(error, nb_id, ids.get(), keep);

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }


  // Shift
  static Vectors* shift(const Vectors& v, int var, double param)
  {
    Format_error error;
    Vectors * ret = NULL;

    if(v.get_type(var-1) == REAL_VALUE)
      ret = v.shift(error, var, (double)param);
    else
      ret = v.shift(error, var, (int)param);

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }


  // Merge
  static Vectors* merge(const Vectors& v,  const boost::python::list& vecs)
  {
    Format_error error;
    Vectors * ret = NULL;

    int nb_vec = len(vecs);
    stat_tool::wrap_util::auto_ptr_array<const Vectors *>
      vects(new const Vectors*[nb_vec]);

    for (int i=0; i<nb_vec; i++)
      vects[i] = extract<Vectors*>(vecs[i]);

    ret = v.merge(error, nb_vec, vects.get());

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }


  static Vectors* merge_variable(const Vectors& v,
                 const boost::python::list& vecs, int ref_sample)
  {
    Format_error error;
    Vectors * ret = NULL;

    int nb_vec = len(vecs);
    stat_tool::wrap_util::auto_ptr_array<const Vectors *>
      vects(new const Vectors*[nb_vec]);

    for (int i=0; i<nb_vec; i++)
      vects[i] = extract<Vectors*>(vecs[i]);

    ret = v.merge_variable(error, nb_vec, vects.get(), ref_sample);

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }


  // Cluster
  static Vectors* cluster_step(const Vectors& v, int variable, int step)
  {
    Format_error error;
    Vectors* ret = v.cluster(error, variable, step);

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Vectors* cluster_limit(const Vectors& v, int variable,
                boost::python::list& limit
                )
  {

    Format_error error;

    int nb_limit = len(limit);
    bool is_float = true;
    int *lint = NULL;
    double *ldouble = NULL;
    Vectors* ret;

    // Test type
    boost::python::extract<int> get_int(limit[0]);
    if (get_int.check())
      {
    is_float = false;
    lint = new int[nb_limit];
      }
    else
      {
    ldouble = new double[nb_limit];
      }

    // Convert list
    for (int i=0; i<nb_limit; i++)
      {
    if(is_float)
      ldouble[i] = extract<int>(limit[i]);
    else
      lint[i] = extract<double>(limit[i]);
      }

    // Call correct function
    if(is_float)
      {
    ret = v.cluster(error, variable, nb_limit, ldouble);
    delete[] ldouble;
      }
    else
      {
    ret = v.cluster(error, variable, nb_limit, lint);
    delete[] lint;
      }

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }


  static Vectors* transcode(const Vectors& v, int variable,
                boost::python::list& symbol
                      )
  {

    Format_error error;

    int nb_symbol = len(symbol);
    stat_tool::wrap_util::auto_ptr_array<int>
      l(new int[nb_symbol]);

    int expected_nb_symbol = (int)(v.get_max_value(variable - 1)
                   - v.get_min_value(variable - 1)) + 1;

    if(nb_symbol != expected_nb_symbol)
      stat_tool::wrap_util::throw_error("Bad number of Symbol");

    for (int i=0; i<nb_symbol; i++)
      l[i] = extract<int>(symbol[i]);

    Vectors* ret = v.transcode(error, variable, l.get());


    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  // Mixture cluster
  static Mv_Mixture_data* mixture_cluster(const Vectors& v, const Mv_Mixture& m)
  {
    Format_error error;
    Mv_Mixture_data* ret = m.cluster(error, v);

    if (ret == NULL)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  // Compare
  static Distance_matrix* compare(const Vectors& v,
                  const Vector_distance& distance)
  {
    Format_error error;
    Distance_matrix* ret = NULL;

    ret = v.comparison(error, distance);

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  // Extract
  static Distribution_data* extract_histogram(const Vectors& v, int variable)
  {
    Format_error error;
    Distribution_data* ret = NULL;

    ret = v.extract(error, variable);

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }




  // Regression
  static Regression* linear_regression(const Vectors& v,
                       int explanatory_var, int response_var)
  {
    Format_error error;
    Regression * ret = NULL;

    ret = v.linear_regression(error, explanatory_var, response_var);

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }


  static Regression* moving_average_dist(const Vectors& v,
                     int explanatory_var, int response_var,
                     const Distribution &dist, char algo)
  {
    Format_error error;
    Regression * ret = NULL;

    if(algo != 'a' && algo != 's')
      {
    PyErr_SetString(PyExc_Exception, "Bad Algorithm");
    boost::python::throw_error_already_set();
      }


    ret = v.moving_average(error, explanatory_var, response_var, dist, algo);

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }


  static Regression* moving_average_list(const Vectors& v,
                     int explanatory_var, int response_var,
                     const boost::python::list& filter, char algo)
  {
    Format_error error;
    Regression * ret = NULL;

    int nb_var = len(filter);
    stat_tool::wrap_util::auto_ptr_array<double>
      vars(new double[nb_var]);

    for (int i=0; i<nb_var; i++)
      vars[i] = extract<double>(filter[i]);

    if(algo != 'a' && algo != 's')
      {
    PyErr_SetString(PyExc_Exception, "Bad Algorithm");
    boost::python::throw_error_already_set();
      }

    ret = v.moving_average(error, explanatory_var, response_var,
               nb_var, vars.get(), algo);

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }


  static Regression* nearest_neighbours(const Vectors& v,
                    int explanatory_var, int response_var,
                    double span , bool weighting)
  {
    Format_error error;
    Regression * ret = NULL;

    ret = v.nearest_neighbor_smoother(error, explanatory_var, response_var, span, weighting);

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Mv_Mixture* mixture_estimation_1(const Vectors& v, const Mv_Mixture& mixt,
                      int nb_iter,
                      boost::python::list force_param)
  {
    bool status = true, several_errors= false;
    Mv_Mixture* ret = NULL;
    bool *fparam= NULL;
    Format_error error;
    ostringstream error_message;
    int nb_fparam, p;
    const int nb_variables = v.get_nb_variable();
    object o;

    nb_fparam= boost::python::len(force_param);

    if (nb_fparam > 0) {
      if (nb_fparam != nb_variables) {
    status = false;
    error_message << "bad size of argument list: " << nb_fparam
              << ": should be the number of variables ("
              << nb_variables << ")";
    PyErr_SetString(PyExc_ValueError, (error_message.str()).c_str());
    throw_error_already_set();
      }
      else {
    fparam = new bool[nb_fparam];
    for (p = 0; p < nb_fparam; p++) {
      o = force_param[p];
      try {
        extract<bool> x(o);
        if (x.check())
          fparam[p]= x();
        else
          status=false;
      }
      catch (...) {
        status = false;
      }
      if (!status) {
        if (several_errors)
          error_message << endl;
        else
          several_errors = true;
        error_message << "incorrect type for element " << p
              << " of argument list: expecting a boolean";
      }
    }
    if (!status) {
      delete [] fparam;
      fparam = NULL;
      PyErr_SetString(PyExc_TypeError, (error_message.str()).c_str());
      throw_error_already_set();
    }
      }
    }
    else {
      nb_fparam = nb_variables;
      fparam = new bool[nb_fparam];
      for (p = 0; p < nb_fparam; p++)
    fparam[p] = false;
    }

    if (status) {

      ret = v.mixture_estimation(error, cout, mixt, nb_iter, fparam);
      if (fparam != NULL) {
    delete [] fparam;
    fparam = NULL;
      }

      if (ret == NULL)
    stat_tool::wrap_util::throw_error(error);

      if (error.get_nb_error() > 0) {
    ret->ascii_write(cout, true);
    delete ret;
    error_message << error << endl;
    PyErr_SetString(PyExc_UserWarning, (error_message.str()).c_str());
    throw_error_already_set();
      }
    }
    return ret;
  }

  static Mv_Mixture* mixture_estimation_2(const Vectors& v, int nb_component,
                      int nb_iter,
                      boost::python::list force_param)
  {
    bool status = true, several_errors= false;
    Mv_Mixture* ret = NULL;
    bool *fparam= NULL;
    Format_error error;
    ostringstream error_message;
    int nb_fparam, p;
    const int nb_variables = v.get_nb_variable();
    object o;

    nb_fparam= boost::python::len(force_param);

    if (nb_fparam > 0) {
      if (nb_fparam != nb_variables) {
    status = false;
    error_message << "bad size of argument list: " << nb_fparam
              << ": should be the number of variables ("
              << nb_variables << ")";
    PyErr_SetString(PyExc_ValueError, (error_message.str()).c_str());
    throw_error_already_set();
      }
      else {
    fparam = new bool[nb_fparam];
    for (p = 0; p < nb_fparam; p++) {
      o = force_param[p];
      try {
        extract<bool> x(o);
        if (x.check())
          fparam[p]= x();
        else
          status=false;
      }
      catch (...) {
        status = false;
      }
      if (!status) {
        if (several_errors)
          error_message << endl;
        else
          several_errors = true;
        error_message << "incorrect type for element " << p
              << " of argument list: expecting a boolean";
      }
    }
    if (!status) {
      delete [] fparam;
      fparam = NULL;
      PyErr_SetString(PyExc_TypeError, (error_message.str()).c_str());
      throw_error_already_set();
    }
      }
    }
    else {
      nb_fparam = nb_variables;
      fparam = new bool[nb_fparam];
      for (p = 0; p < nb_fparam; p++)
    fparam[p] = false;
    }

   if (status) {

      ret = v.mixture_estimation(error, cout, nb_component, nb_iter, fparam);
      if (fparam != NULL) {
    delete [] fparam;
    fparam = NULL;
      }

      if (ret == NULL)
    stat_tool::wrap_util::throw_error(error);

      if (error.get_nb_error() > 0) {
    ret->ascii_write(cout, true);
    delete ret;
    error_message << error << endl;
    PyErr_SetString(PyExc_UserWarning, (error_message.str()).c_str());
    throw_error_already_set();
      }
    }
    return ret;
  }



  static string contingency_table(const Vectors& v, int variable1, int variable2,
                   const string& filename, char format)
  {
     Format_error error;
     std::stringstream s;
     bool ret;
     ret = v.contingency_table(error, s, variable1, variable2, filename.c_str(), format);

     if(!ret)
       stat_tool::wrap_util::throw_error(error);

     return s.str();
  }

  static string variance_analysis(const Vectors& v, int class_variable,
                   int response_variable, int response_type,
                   const string& filename, char format)
  {
     Format_error error;
     std::stringstream s;
     bool ret;

     ret = v.variance_analysis(error, s, class_variable, response_variable, response_type,
             filename.c_str(), format);

     if(!ret)
       stat_tool::wrap_util::throw_error(error);

     return s.str();
  }
};



void class_vectors()
{

  class_< Vectors, bases< STAT_interface > >
    ("_Vectors", "Vectors (2 dimensions list)",
     init<>())
    // constructors
    .def("__init__", make_constructor(VectorsWrap::build_from_lists))
    .def("__init__", make_constructor(VectorsWrap::read_from_file))

    // Python Operators
    .def(self_ns::str(self)) // __str__
    .def("__len__", &Vectors::get_nb_vector,
     "Return the number of vectors")
    .def("__getitem__", VectorsWrap::get_item)

    .def("get_nb_variable", &Vectors::get_nb_variable,
     "Return the number of variables")
    .def("get_nb_vector", &Vectors::get_nb_vector,
     "Return the number of vectors")

    // Identifiers
    .def("get_identifiers", VectorsWrap::get_identifiers,
     "Return the list of identifiers")

    // Select
    .def("value_select", VectorsWrap::value_select,
     return_value_policy< manage_new_object >(),
     python::args("variable", "min", "max", "keep"),
     "Selection of individuals according to the values taken by a variable")

    .def("select_variable", VectorsWrap::select_variable,
     return_value_policy< manage_new_object >(),
     python::args("variables", "keep"),
     "select variable given a list of index")

    .def("select_individual", VectorsWrap::select_individual,
     return_value_policy< manage_new_object >(),
     python::args("identifiers", "keep"),
     "Select individuals given a list of identifiers")

    // Extract
    .def("extract", VectorsWrap::extract_histogram,
     return_value_policy< manage_new_object >(),
     python::args("variable"),
     "Extract histogram")

    // Regression
    .def("linear_regression", VectorsWrap::linear_regression,
     return_value_policy< manage_new_object >(),
     python::args("explanatory_var", "response_var"),
     "Linear regression")

    .def("moving_average_regression", VectorsWrap::moving_average_dist,
     return_value_policy< manage_new_object >(),
     python::args("explanatory_var", "response_var", "dist", "algo"),
     "Linear regression (moving average)")

    .def("moving_average_regression", VectorsWrap::moving_average_list,
     return_value_policy< manage_new_object >(),
     python::args("explanatory_var", "response_var", "filters", "algo"),
     "Linear regression (moving average)")

    .def("nearest_neighbours_regression", VectorsWrap::nearest_neighbours,
     return_value_policy< manage_new_object >(),
     python::args("explanatory_var", "response_var", "span", "weighting"),
     "Linear regression (nearest neighbours)")

    .def("mixture_estimation_wrap", VectorsWrap::mixture_estimation_1,
     return_value_policy< manage_new_object >(),
     python::args("initial_mixture", "nb_max_iteration", "force_param"),
     "Mixture estimation (EM algorithm with initial model)")

    .def("mixture_estimation_wrap", VectorsWrap::mixture_estimation_2,
     return_value_policy< manage_new_object >(),
     python::args("nb_component", "nb_max_iteration", "force_param"),
     "Mixture estimation (EM algorithm with fixed number of components)")

    // Merge
    .def("merge", VectorsWrap::merge,
     return_value_policy< manage_new_object >(),
     "Merge vectors")

    .def("merge_variable", VectorsWrap::merge_variable,
     return_value_policy< manage_new_object >(),
     "Merge variables" )

    // Cluster
    .def("cluster_step", VectorsWrap::cluster_step,
     return_value_policy< manage_new_object >(),
     python::args("variable", "step"),
     "Cluster Step"
     )

    .def("cluster_limit", VectorsWrap::cluster_limit,
     return_value_policy< manage_new_object >(),
     python::args("variable", "limits"),
     "Cluster limit")

    .def("transcode", VectorsWrap::transcode,
     return_value_policy< manage_new_object >(),
     python::args("variable", "symbols"),
     "Transcode")

    // Mixture clustering
    .def("mixture_cluster", VectorsWrap::mixture_cluster,
     return_value_policy< manage_new_object >(),
     python::args("model"),
     "Cluster individuals using mixture model"
     )

    // Compare
    .def("compare", VectorsWrap::compare,
     return_value_policy< manage_new_object >(),
     python::arg("distance"),
     "Compare Vectors given a VectorDistance")

    // Others
    .def("contingency_table", VectorsWrap::contingency_table,
     "Return a string with the contingency_table")

    .def("variance_analysis", VectorsWrap::variance_analysis,
     "Return a string with the variance analysis")


    // Output
    .def("ascii_data_write", VectorsWrap::ascii_data_write,
     "Return a string with the object representation")

    // save to file
    .def("file_ascii_write", VectorsWrap::file_ascii_write,
     "Save vector summary into a file")

   .def("file_ascii_data_write", VectorsWrap::file_ascii_data_write,
     "Save vector data into a file")

    // Shift
    .def("shift", VectorsWrap::shift,
     return_value_policy< manage_new_object >(),
     python::args("variable", "param"),
     "Shift")
    ;

}



// Vectors

class VectorDistanceWrap
{

public:

  static boost::shared_ptr<Vector_distance> build_from_types (boost::python::list& types,
                                  boost::python::list& weigths,
                                  int distance_type)
  {
    Vector_distance* dist;
    int nb_variable;

    nb_variable = boost::python::len(types);

    stat_tool::wrap_util::auto_ptr_array<double>
      variable_weight (new double[nb_variable]);

    stat_tool::wrap_util::auto_ptr_array<int>
      variable_type (new int[nb_variable]);

    // Extract each element of the vector
    for(int i=0; i<nb_variable; i++)
      {
    variable_type[i] = boost::python::extract<int>(types[i]);
    variable_weight[i] = boost::python::extract<double>(weigths[i]);
      }

    dist = new Vector_distance(nb_variable, variable_type.get(),
                   variable_weight.get(), distance_type);

    return boost::shared_ptr<Vector_distance>(dist);
  }


  static boost::shared_ptr<Vector_distance> read_from_file(char* filename)
  {
    Vector_distance* dist;
    Format_error error;
    dist = vector_distance_ascii_read(error, filename);

    if(!dist)
      stat_tool::wrap_util::throw_error(error);


    return boost::shared_ptr<Vector_distance>(dist);
  }

  static void file_ascii_write(const Vector_distance& d, const char* path, bool exhaustive)
  {
    bool result = true;
    Format_error error;

    result = d.ascii_write(error, path,exhaustive);
    if (!result)
       stat_tool::wrap_util::throw_error(error);

  }

};



void class_vectordistance()
{
  enum_<stat_tool::wrap_util::UniqueInt<2, 1> >("DistanceType")
    .value("ABSOLUTE_VALUE", ABSOLUTE_VALUE)
    .value("QUADRATIC", QUADRATIC)
    .export_values()
    ;


  class_< Vector_distance, bases< STAT_interface > >
    ("_VectorDistance", "Distance class", init<>())
    // constructors
    .def("__init__", make_constructor(VectorDistanceWrap::build_from_types))
    .def("__init__", make_constructor(VectorDistanceWrap::read_from_file))

    .def(self_ns::str(self))
    .def("__len__", &Vector_distance::get_nb_variable)

    .def("file_ascii_write", VectorDistanceWrap::file_ascii_write,
     "Save vector distance summary into a file")


    ;

}

