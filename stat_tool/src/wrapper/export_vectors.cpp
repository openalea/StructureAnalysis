/*------------------------------------------------------------------------------
 *
 *        VPlants.Stat_Tool : VPlants Statistics module
 *
 *        Copyright 2006-2014 CIRAD/INRA/Inria Virtual Plants
 *
 *        File author(s): Yann Guedon <yann.guedon@cirad.fr>
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
 *        $Id: export_vectors.cpp 18040 2015-04-23 07:17:16Z guedon $
 *
 *-----------------------------------------------------------------------------*/



#include "wrapper_util.h"
#include "export_base.h"

#include "stat_tool/stat_tools.h"

#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/regression.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/multivariate_mixture.h"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/make_constructor.hpp>

#include "boost_python_aliases.h"
#include "export_base.h"

using namespace boost::python;
using namespace boost;
using namespace stat_tool;


// Vectors
#define WRAP VectorsWrap
class WRAP
{

public:

  WRAP_METHOD2(Vectors,linear_regression, Regression, int, int);
  WRAP_METHOD2(Vectors,comparison, DistanceMatrix, VectorDistance, bool);
  WRAP_METHOD_SPREADSHEET_WRITE( Vectors);WRAP_METHOD2(Vectors,scaling,
       Vectors, int, int);
  WRAP_METHOD2(Vectors,round, Vectors, int, int);
  WRAP_METHOD_FILE_ASCII_WRITE( Vectors)

  static Vectors*
  read_from_file(char *filename)
  {
    Vectors *vec;
    StatError error;

    vec = vectors_ascii_read(error, filename);

    if (vec)
      {
        return vec;
      }
    else
      {
        stat_tool::wrap_util::throw_error(error);
      }

  }


  static Vectors*
  build_from_lists_and_types(boost::python::list& array, boost::python::list &ident, boost::python::list &user_types)
  {
    int nb_vector = boost::python::len(array);
    int nb_ident = boost::python::len(ident);
    int nb_types = boost::python::len(user_types);

    int dummy;

    int *types = 0; // 0 for int and 1 for float/double
    int *identifier = 0;
    int nb_variable = -1, this_nb_variable=-1;
    int **int_vector = 0;
    double **float_vector = 0;
    bool is_float = false;
    bool error = false;
    Vectors* vec = 0;

    // check that length are correct
    if (nb_types == 0 || nb_vector==0){
        return vec;
    }

    // get nb_variable
    boost::python::list vectemp = boost::python::extract<boost::python::list>(array[0]);
    nb_variable = boost::python::len(vectemp);

    // check hat nb_types == nb_vectors
    if (nb_types != nb_variable){
        error=true;
        ostringstream error_message;
        error_message << "Types size must be equal to number of variables."  << endl;
        PyErr_SetString(PyExc_ValueError, (error_message.str()).c_str());
        throw_error_already_set();
//        return vec;
    }

    // extract the identifiers if needed
    if (nb_ident > 0)
    {
      identifier = new int[nb_ident];
      for (int ii = 0; ii < nb_ident; ii++)
        {
          identifier[ii] = boost::python::extract<int>(ident[ii]);
        }
    }

    // check that nb_identifiers == nb_vector
    if (nb_ident != nb_vector){
        error=true;
        ostringstream error_message;
        error_message << "Idnetifiers length must be equal to number of vectors."  << endl;
        PyErr_SetString(PyExc_ValueError, (error_message.str()).c_str());
        throw_error_already_set();
    //    return vec;
    }

    // extract the types
    types = new int[nb_types];
    int nb_variable_int = 0;
    int nb_variable_float = 0;
    for (int ii = 0; ii < nb_types; ii++){
        dummy = boost::python::extract<int>(user_types[ii]);
        if (dummy==0){
            types[ii] = INT_VALUE;
            nb_variable_int++;
        }
        else {
            types[ii] = REAL_VALUE;
            nb_variable_float++;
        }
    }

    // count nb integer and float in types
    //
    int i,  index_int, index_float;

    // allocate only required memory according to types
    int_vector = new int*[nb_vector];
    float_vector = new double*[nb_vector];
    for (int e = 0; e < nb_vector; e++)
    {
       int_vector[e] = new int[nb_variable_int];
       float_vector[e] = new double[nb_variable_float];
    }

    // scan all vectors
    try {

    for (int vi = 0; vi < nb_vector; vi++)
    {
      // look at each vector
      boost::python::list vectemp = boost::python::extract<boost::python::list>(array[vi]);
      this_nb_variable = boost::python::len(vectemp);
      if (this_nb_variable!=nb_variable){
        ostringstream error_message;
        error_message << "vectors sizes must be equal"  << endl;
        PyErr_SetString(PyExc_ValueError, (error_message.str()).c_str());
        throw_error_already_set();
        return vec;
      }
      index_int = 0;
      index_float= 0;
      for (int i = 0; i < nb_variable; i++)
        {
          if (types[i]==INT_VALUE){
            int v = boost::python::extract<int>(vectemp[i]);
            int_vector[vi][index_int] = v;
            index_int++;
          }
          else{
              double v = boost::python::extract<double>(vectemp[i]);
              float_vector[vi][index_float] = v;
              index_float++;
          }
        }
    }
    }
    catch (...)
      {
        error = true;
      }
    if (!error){
      vec = new Vectors(nb_vector, identifier, nb_variable, types, int_vector, float_vector, false);
    }


    // Delete memory
    if (int_vector)
    {
      for (int i = 0; i < nb_vector; i++)
        {
          delete[] int_vector[i];
        }
      delete[] int_vector;
    }

    if (float_vector)
    {
      for (int i = 0; i < nb_vector; i++)
      {
        delete[] float_vector[i];
      }
      delete[] float_vector;
    }
    // delete identifier
    if (identifier)
    {
      delete[] identifier;
    }
    // delete type
    if (types)
    {
      delete[] types;
    }

    if (error)
      {
        throw_error_already_set();
      }

    return vec;
  }

  static Vectors*
  build_from_lists(boost::python::list& array, boost::python::list &ident)
  {
    int nb_vector = boost::python::len(array);
    int *identifier = 0;
    int nb_variable = -1;
    int **int_vector = 0;
    double **float_vector = 0;
    bool is_float = false;
    bool error = false;
    Vectors* vec = 0;

    // extract the identifiers if needed
    int nb_ident = boost::python::len(ident);
    if (nb_ident > 0)
      {
        identifier = new int[nb_ident];
        for (int ii = 0; ii < nb_ident; ii++)
          {
            identifier[ii] = boost::python::extract<int>(ident[ii]);
          }
      }

    // now the data itself.
    try
      {
        // For each vector
        for (int vi = 0; vi < nb_vector; vi++)
          {
            boost::python::list vec = boost::python::extract<
                boost::python::list>(array[vi]);

            // nb_variable is the length of the first vector
            if (vi == 0)
              {
                nb_variable = boost::python::len(vec);

                // Check the type
                boost::python::extract<int> get_int(vec[0]);
                if (get_int.check()) // Array of int
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
                        error_message
                            << "Incorrect type (expect int or float values)"
                            << vi << endl;
                        PyErr_SetString(PyExc_TypeError,
                            (error_message.str()).c_str());
                        error = true;
                      }
                  }
              }
            // Check the size of the other vectors

            else if (boost::python::len(vec) != nb_variable)
              {
                // Error : Bad size
                ostringstream error_message;
                error_message << "Incorrect size of vector " << vi << endl;
                PyErr_SetString(PyExc_ValueError, (error_message.str()).c_str());
                error = true;
              }

            // Build array if necessary
            if (!int_vector && !is_float)
              {
                int_vector = new int*[nb_vector];
                for (int e = 0; e < nb_vector; e++)
                  {
                    int_vector[e] = new int[nb_variable];
                  }

              }
            else if (!float_vector && is_float)
              {
                float_vector = new double*[nb_vector];
                for (int e = 0; e < nb_vector; e++)
                  {
                    float_vector[e] = new double[nb_variable];
                  }
              }

            // Extract each element of the vector
            for (int i = 0; i < nb_variable; i++)
              {
                if (!is_float)
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
    catch (...)
      {
        error = true;
      }

    // Call constructor
    if (!error)
      {
        if (!is_float)
          vec = new Vectors(nb_vector, identifier, nb_variable, int_vector);
        else
          vec = new Vectors(nb_vector, identifier, nb_variable, float_vector);
      }

    // Delete memory
    if (int_vector)
      {
        for (int i = 0; i < nb_vector; i++)
          {
            delete[] int_vector[i];
          }
        delete[] int_vector;
      }

    if (float_vector)
      {
        for (int i = 0; i < nb_vector; i++)
          {
            delete[] float_vector[i];
          }
        delete[] float_vector;
      }

    if (identifier)
      {
        delete[] identifier;
      }

    if (error)
      {
        throw_error_already_set();
      }

    return vec;

  }
  ;

  static boost::python::list
  get_item(const Vectors* vec, int index)
  {
    // Test index
    if (index < 0 || index >= vec->get_nb_vector())
      {
        PyErr_SetString(PyExc_IndexError, "vector index out of bound");
        boost::python::throw_error_already_set();
      }

    boost::python::list l;

    int nb_var = vec->get_nb_variable();
    for (int var = 0; var < nb_var; var++)
      {
        if ((vec->get_type(var) == INT_VALUE) || (vec->get_type(var) == STATE))
          l.append(vec->get_int_vector(index, var));
        else
          l.append(vec->get_real_vector(index, var));
      }

    return l;
  }

  static boost::python::list
  get_identifiers(const Vectors& vec)
  {
    boost::python::list l;

    int nb_vec = vec.get_nb_vector();
    for (int v = 0; v < nb_vec; v++)
      {
        l.append(vec.get_identifier(v));
      }

    return l;
  }

  static std::string
  ascii_data_write(const Vectors& d, bool exhaustive)
  {
    std::stringstream s;
    std::string res;

    d.ascii_data_write(s, exhaustive);
    res = s.str();

    return res;

  }

  static void
  file_ascii_data_write(const Vectors& d, const char* path, bool exhaustive)
  {
    bool result = true;
    StatError error;

    result = d.ascii_data_write(error, path, exhaustive);
    if (!result)
      stat_tool::wrap_util::throw_error(error);

  }

  static Vectors*
  value_select(const Vectors& v, int variable, const object& min,
      const object& max, bool keep)
  {
    StatError error;
    Vectors * ret = NULL;
    std::stringstream s;

    boost::python::extract<int> get_min(min);
    boost::python::extract<int> get_max(max);

    if (get_min.check() && get_max.check()) // Array of int

      {
        int mi = get_min();
        int ma = get_max();
        ret = v.value_select(error, s, variable, mi, ma, keep);
      }
    else
      {
        double mi = extract<double> (min);
        double ma = extract<double> (max);
        ret = v.value_select(error, s, variable, mi, ma, keep);
      }

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    cout << s.str() << endl;

    return ret;
  }

  static Vectors*
  select_variable(const Vectors& v, const boost::python::list& variables,
      bool keep)
  {
    StatError error;
    Vectors * ret = NULL;

    int nb_var = len(variables);
    stat_tool::wrap_util::auto_ptr_array<int> vars(new int[nb_var]);

    for (int i = 0; i < nb_var; i++)
      vars[i] = extract<int> (variables[i]);

    ret = v.select_variable(error, nb_var, vars.get(), keep);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Vectors*
  select_individual(const Vectors& v, const boost::python::list& identifiers,
      bool keep)
  {
    StatError error;
    Vectors * ret = NULL;

    int nb_id = len(identifiers);
    stat_tool::wrap_util::auto_ptr_array<int> ids(new int[nb_id]);

    for (int i = 0; i < nb_id; i++)
      ids[i] = extract<int> (identifiers[i]);

    ret = v.select_individual(error, nb_id, ids.get(), keep);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  // Shift
  static Vectors*
  shift(const Vectors& v, int var, double param)
  {
    StatError error;
    Vectors * ret = NULL;

    if (v.get_type(var - 1) == REAL_VALUE)
      ret = v.shift(error, var, (double) param);
    else
      ret = v.shift(error, var, (int) param);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  // Merge
  static Vectors*
  merge(const Vectors& v, const boost::python::list& vecs)
  {
    StatError error;
    Vectors * ret = NULL;

    int nb_vec = len(vecs);
    stat_tool::wrap_util::auto_ptr_array<const Vectors *> vects(
        new const Vectors*[nb_vec]);

    for (int i = 0; i < nb_vec; i++)
      vects[i] = extract<Vectors*> (vecs[i]);

    ret = v.merge(error, nb_vec, vects.get());

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Vectors*
  merge_variable(const Vectors& v, const boost::python::list& vecs,
      int ref_sample)
  {
    StatError error;
    Vectors * ret = NULL;

    int nb_vec = len(vecs);
    stat_tool::wrap_util::auto_ptr_array<const Vectors *> vects(
        new const Vectors*[nb_vec]);

    for (int i = 0; i < nb_vec; i++)
      vects[i] = extract<Vectors*> (vecs[i]);

    ret = v.merge_variable(error, nb_vec, vects.get(), ref_sample);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  // check
  static bool
  check(Vectors& v)
  {
    StatError error;
    bool ret = v.check(error);
    return ret;
  }

  // Cluster
  static Vectors*
  cluster_step(const Vectors& v, int variable, int step)
  {
    StatError error;
    Vectors* ret = v.cluster(error, variable, step);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Vectors*
  cluster_limit(const Vectors& v, int variable, boost::python::list& limit)
  {

    StatError error;

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
    for (int i = 0; i < nb_limit; i++)
      {
        if (is_float)
          ldouble[i] = extract<int> (limit[i]);
        else
          lint[i] = extract<double> (limit[i]);
      }

    // Call correct function
    if (is_float)
      {
        ret = v.cluster(error, variable, nb_limit, ldouble);
        delete[] ldouble;
      }
    else
      {
        ret = v.cluster(error, variable, nb_limit, lint);
        delete[] lint;
      }

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Vectors*
  transcode(const Vectors& v, int variable, boost::python::list& symbol)
  {

    StatError error;

    int nb_symbol = len(symbol);
    stat_tool::wrap_util::auto_ptr_array<int> l(new int[nb_symbol]);

    int expected_nb_symbol = (int) (v.get_max_value(variable - 1)
        - v.get_min_value(variable - 1)) + 1;

    if (nb_symbol != expected_nb_symbol)
      stat_tool::wrap_util::throw_error("Bad number of Symbol");

    for (int i = 0; i < nb_symbol; i++)
      l[i] = extract<int> (symbol[i]);

    Vectors* ret = v.transcode(error, variable, l.get());

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  // Mixture cluster
  static MultivariateMixtureData*
  mixture_cluster(const Vectors& v, const MultivariateMixture& m)
  {
    StatError error;
    MultivariateMixtureData* ret = NULL;
    ret = m.cluster(error, v);
    if (ret == NULL)
      stat_tool::wrap_util::throw_error(error);
    return ret;
  }

  // Extract
  static DiscreteDistributionData*
  extract_histogram(const Vectors& v, int variable)
  {
    StatError error;
    DiscreteDistributionData* ret = NULL;
    ret = v.extract(error, variable);
    if (!ret)
      stat_tool::wrap_util::throw_error(error);
    return ret;
  }

  static Regression*
  moving_average_dist(const Vectors& v, int explanatory_var, int response_var,
      const Distribution &dist, char algo)
  {
    StatError error;
    Regression * ret = NULL;
    ret = v.moving_average(error, explanatory_var, response_var, dist, algo);
    if (!ret)
      stat_tool::wrap_util::throw_error(error);
    return ret;
  }

  static Regression*
  moving_average_list(const Vectors& v, int explanatory_var, int response_var,
      const boost::python::list& filter, char algo)
  {
    StatError error;
    Regression * ret = NULL;
    double sum = 0;
    int nb_var = len(filter);
    stat_tool::wrap_util::auto_ptr_array<double> vars(
        new double[nb_var * 2 + 1]);

    // code from the former AML code
    int i = 0;
    nb_var--;

    for (i = 0; i < nb_var; i++)
      {
        vars[i] = extract<double> (filter[i]);
        vars[nb_var * 2 - i] = vars[i];
        sum += 2 * vars[i];
      }
    // i = n
    vars[i] = extract<double> (filter[i]);
    sum += vars[i];

    //normalization
    for (i = 0; i < 2 * nb_var + 1; i++)
      {
        vars[i] = vars[i] / sum;
      }

    ret = v.moving_average(error, explanatory_var, response_var, nb_var,
        vars.get(), algo);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Regression*
  nearest_neighbours(const Vectors& v, int explanatory_var, int response_var,
      double span, bool weighting)
  {
    StatError error;
    Regression * ret = NULL;

    ret = v.nearest_neighbor_smoother(error, explanatory_var, response_var,
        span, weighting);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static MultivariateMixture*
  mixture_estimation_model(const Vectors& v, const MultivariateMixture& mixt, int nb_iter,
      boost::python::list force_param)
  {
    bool status = true, several_errors = false;
    MultivariateMixture* ret = NULL;
    bool *fparam = NULL;
    StatError error;
    ostringstream error_message;
    int nb_fparam, p;
    const int nb_variables = v.get_nb_variable();
    object o;

    nb_fparam = boost::python::len(force_param);

    if (nb_fparam > 0)
      {
        if (nb_fparam != nb_variables)
          {
            status = false;
            error_message << "bad size of argument list: " << nb_fparam
                << ": should be the number of variables (" << nb_variables
                << ")";
            PyErr_SetString(PyExc_ValueError, (error_message.str()).c_str());
            throw_error_already_set();
          }
        else
          {
            fparam = new bool[nb_fparam];
            for (p = 0; p < nb_fparam; p++)
              {
                o = force_param[p];
                try
                  {
                    extract<bool> x(o);
                    if (x.check())
                      fparam[p] = x();
                    else
                      status = false;
                  }
                catch (...)
                  {
                    status = false;
                  }
                if (!status)
                  {
                    if (several_errors)
                      error_message << endl;
                    else
                      several_errors = true;
                    error_message << "incorrect type for element " << p
                        << " of argument list: expecting a boolean";
                  }
              }
            if (!status)
              {
                delete[] fparam;
                fparam = NULL;
                PyErr_SetString(PyExc_TypeError, (error_message.str()).c_str());
                throw_error_already_set();
              }
          }
      }
    else
      {
        nb_fparam = nb_variables;
        fparam = new bool[nb_fparam];
        for (p = 0; p < nb_fparam; p++)
          fparam[p] = false;
      }

    if (status)
      {

        ret = v.mixture_estimation(error, cout, mixt, nb_iter, fparam);
        if (fparam != NULL)
          {
            delete[] fparam;
            fparam = NULL;
          }

        if (ret == NULL)
          stat_tool::wrap_util::throw_error(error);

        if (error.get_nb_error() > 0)
          {
            ret->ascii_write(cout, true);
            delete ret;
            error_message << error << endl;
            PyErr_SetString(PyExc_UserWarning, (error_message.str()).c_str());
            throw_error_already_set();
          }
      }
    return ret;
  }

  static MultivariateMixture*
  mixture_estimation_nb_component(const Vectors& v, int nb_component, int nb_iter,
      boost::python::list force_param)
  {
    bool status = true, several_errors = false;
    MultivariateMixture* ret = NULL;
    bool *fparam = NULL;
    StatError error;
    ostringstream error_message;
    int nb_fparam, p;
    const int nb_variables = v.get_nb_variable();
    object o;

    nb_fparam = boost::python::len(force_param);

    if (nb_fparam > 0)
      {
        if (nb_fparam != nb_variables)
          {
            status = false;
            error_message << "bad size of argument list: " << nb_fparam
                << ": should be the number of variables (" << nb_variables
                << ")";
            PyErr_SetString(PyExc_ValueError, (error_message.str()).c_str());
            throw_error_already_set();
          }
        else
          {
            fparam = new bool[nb_fparam];
            for (p = 0; p < nb_fparam; p++)
              {
                o = force_param[p];
                try
                  {
                    extract<bool> x(o);
                    if (x.check())
                      fparam[p] = x();
                    else
                      status = false;
                  }
                catch (...)
                  {
                    status = false;
                  }
                if (!status)
                  {
                    if (several_errors)
                      error_message << endl;
                    else
                      several_errors = true;
                    error_message << "incorrect type for element " << p
                        << " of argument list: expecting a boolean";
                  }
              }
            if (!status)
              {
                delete[] fparam;
                fparam = NULL;
                PyErr_SetString(PyExc_TypeError, (error_message.str()).c_str());
                throw_error_already_set();
              }
          }
      }
    else
      {
        nb_fparam = nb_variables;
        fparam = new bool[nb_fparam];
        for (p = 0; p < nb_fparam; p++)
          fparam[p] = false;
      }

    if (status)
      {

        ret = v.mixture_estimation(error, cout, nb_component, nb_iter, fparam);
        if (fparam != NULL)
          {
            delete[] fparam;
            fparam = NULL;
          }

        if (ret == NULL)
          stat_tool::wrap_util::throw_error(error);

        if (error.get_nb_error() > 0)
          {
            ret->ascii_write(cout, true);
            delete ret;
            error_message << error << endl;
            PyErr_SetString(PyExc_UserWarning, (error_message.str()).c_str());
            throw_error_already_set();
          }
      }
    return ret;
  }

  static string
  contingency_table(const Vectors& v, int variable1, int variable2,
      const string& filename, char format)
  {
    StatError error;
    std::stringstream s;
    bool ret;
    ret = v.contingency_table(error, s, variable1, variable2, filename.c_str(),
        format);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return s.str();
  }

  static string
  variance_analysis(const Vectors& v, int class_variable,
      int response_variable, int response_type, const string& filename,
      char format)
  {
    StatError error;
    std::stringstream s;
    bool ret;

    ret = v.variance_analysis(error, s, class_variable, response_variable,
        response_type, filename.c_str(), format);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return s.str();
  }

  static string
  rank_correlation_computation(const Vectors& input, int type, const string &filename)
  {
    StatError error;
    std::stringstream os;
    bool ret;
    ret = input.rank_correlation_computation(error, os, type, filename.c_str());
    //std::cout << os.str()<<endl;
    return os.str();
  }

  static MultiPlotSet*
  get_plotable(const Vectors& p)
  {
    // to be checked
    MultiPlotSet* ret = p.get_plotable();
    return ret;
  }

  static bool
  select_step(Vectors &input, int variable, double step)
  {
    StatError error;
    bool ret;

    ret = input.select_step(error, variable, step);
    if (!ret)
      stat_tool::wrap_util::throw_error(error);
    return ret;
  }

  static Histogram*
  get_marginal_histogram(Vectors &input, int variable)
  {
    StatError error;
    Histogram *ret;
    if (variable<input.get_nb_variable() && variable>=0)
        {
            ret = input.get_marginal_histogram(variable);
            return ret;
        }
    else
        {
            cerr << "variable must be positive or zero and stricly below nb_variable"<< endl;
        }
    if (!ret)
      stat_tool::wrap_util::throw_error(error);

  }

  static void plot_write(const Vectors &input, const std::string& prefix,
  const std::string& title)
  {
      StatError error;
      input.plot_write(error, prefix.c_str(), title.c_str());
  }


};

void
class_vectors()
{

  class_< Vectors, bases< StatInterface > >
  ("_Vectors", "Vectors (2 dimensions list)", init<>())
  // constructors
  .def("__init__", make_constructor(VectorsWrap::build_from_lists))
  .def("__init__", make_constructor(VectorsWrap::build_from_lists_and_types))
  .def("__init__", make_constructor(VectorsWrap::read_from_file))

  // Python Operators
  .def(self_ns::str(self)) // __str__
  .def("__len__", &Vectors::get_nb_vector,
      "Return the number of vectors")
  .def("__getitem__", WRAP::get_item)

   // properties
  .add_property("nb_variable", &Vectors::get_nb_variable, "Return the number of variables")
  .add_property("nb_vector", &Vectors::get_nb_vector,"Return the number of vectors")

  // getter
  .def("get_identifiers", WRAP::get_identifiers, args("index"),
      "function that returns the list of identifiers")
  .def("get_min_value", &Vectors::get_min_value, args("index"),
      "Return the min value of a variable.  Example: v.get_min_value(4)")
  .def("get_max_value", &Vectors::get_max_value, args("index"),
      "Return the max value of a variable.  Example: v.get_max_value(4)")
  .def("get_mean", &Vectors::get_mean, args("index"),
      "Return the mean value of a vector. examepl v.get_mean(4)")
  .def("get_type", &Vectors::get_type, args("index"),
      "Return the type of a variable. example v.get_type(4)")
   //.def ("skewness", &Vectors::skewness_computation, args("variable"),
   //     "Returns skewness coefficient given variable")


  // python modules
  DEF_RETURN_VALUE("value_select", WRAP::value_select,
      args("variable", "min", "max", "keep"),"See ValueSelect")
  DEF_RETURN_VALUE("select_variable", WRAP::select_variable,
      args("variables", "keep"), "See SelectVariable")
  DEF_RETURN_VALUE("select_individual", WRAP::select_individual,
      args("identifiers", "keep"),"See SelectIndividual")
  DEF_RETURN_VALUE("extract", WRAP::extract_histogram,args("variable"),
      "Extract histogram")
  DEF_RETURN_VALUE("transcode", WRAP::transcode,
      args("variable", "symbols"), "See Transcode")
  DEF_RETURN_VALUE("cluster_step", WRAP::cluster_step,
      args("variable", "step"), "See Cluster" )
  DEF_RETURN_VALUE("cluster_limit", WRAP::cluster_limit,
      args("variable", "limits"), "See Cluster")
  DEF_RETURN_VALUE("shift", WRAP::shift,
      args("variable", "param"), "See Shift")
  DEF_RETURN_VALUE_NO_ARGS("merge", WRAP::merge,
      "See Merge")
  DEF_RETURN_VALUE("round", WRAP::scaling,
      args("variable", "scaling_coeff"), "See Round")
  DEF_RETURN_VALUE_NO_ARGS("merge_variable", WRAP::merge_variable,
      "See MergeVariable" )
  .def("contingency_table", WRAP::contingency_table,
      "See ContingencyTable")
  .def("variance_analysis", WRAP::variance_analysis,
      "See VarianceAnalysis")


  // Used in Python modules such as:
  // ---------------- Moving Average
  DEF_RETURN_VALUE("moving_average_regression_values", WRAP::moving_average_list,
      args("explanatory_var", "response_var", "dist", "algo"), "Linear regression (MovingAverage)")
  DEF_RETURN_VALUE("moving_average_regression_distribution", WRAP::moving_average_dist,
      args("explanatory_var", "response_var", "filters", "algo"), "Linear regression (See MovingAverage)")
  // ---------------------Estimation
  DEF_RETURN_VALUE("linear_regression", WRAP::linear_regression,
      args("explanatory_var", "response_var"), "TODO Linear regression")
  DEF_RETURN_VALUE("nearest_neighbours_regression", WRAP::nearest_neighbours,
      args("explanatory_var", "response_var", "span", "weighting"),
      "TODO Linear regression (nearest neighbours)")
  DEF_RETURN_VALUE("mixture_estimation_model", WRAP::mixture_estimation_model,
      args("initial_mixture", "nb_max_iteration", "force_param"),
      "TODO Mixture estimation (EM algorithm with initial model)")
  DEF_RETURN_VALUE("mixture_estimation_nb_component", WRAP::mixture_estimation_nb_component,
      args("nb_component", "nb_max_iteration", "force_param"),
      "TODO Mixture estimation (EM algorithm with fixed number of components)")



  DEF_RETURN_VALUE("mixture_cluster", WRAP::mixture_cluster,
      args("model"), "TODOCluster individuals using mixture model" )
  DEF_RETURN_VALUE("compare", WRAP::comparison,
      args("distance"), "TODOCompare Vectors given a VectorDistance")
  DEF_RETURN_VALUE("scaling", WRAP::scaling,
      args("variable", "scaling_coeff"), "TODOScales vectors")


  .def("variance_analysis", WRAP::variance_analysis,
      "Return a string with the variance analysis")
  .def("rank_correlation_computation", WRAP::rank_correlation_computation,
      args("type", "filename"), "Rank correlation computation")

  //others
  DEF_RETURN_VALUE_NO_ARGS("get_plotable", WRAP::get_plotable,
      "Return a plotable")
  .def("plot_write", WRAP::plot_write, args("prefix", "title"), "Write GNUPLOT files")

  .def("check", WRAP::check, args("todo"),
      "todo check vectors")
  .def("ascii_data_write", WRAP::ascii_data_write,
      "Return a string with the object representation")
  .def("file_ascii_write", WRAP::file_ascii_write,
      "Save vector summary into a file")
  .def("file_ascii_data_write", WRAP::file_ascii_data_write,
      "Save vector data into a file")
  .def("spreadsheet_write", WRAP::spreadsheet_write,
      "Save data into CSV file")
  .def("select_step", WRAP::select_step, 
    args("variable", "step"), "select_step(step) refedine the step of the histogram. Step must be >0")
  DEF_RETURN_VALUE("get_marginal_histogram", WRAP::get_marginal_histogram, 
    args("variable"), "get_marginal_histogram(nb_variable) construct marginal histogram of the vector given for the variable provided. The variable must be >=0 and less than nb_variable.")

  ;

  /*

    Vectors (int inb_vector , int *iidentifier , int inb_variable , int *itype ,  bool init_flag = false)
    Vectors(int inb_vector , int *iidentifier , int inb_variable , int *itype ,  int **iint_vector , double **ireal_vector);
    Vectors(const Vectors &vec , int inb_vector , int *index);
    





   Vectors();
   Vectors(const Vectors &vec , int inb_vector , int *index);

   bool check(StatError &error);

   -->bool plot_write(StatError &error , const char *prefix ,
   const char *title = 0) const;

   double mean_absolute_deviation_computation(int variable) const;
   double mean_absolute_difference_computation(int variable) const;
   double skewness_computation(int variable) const;
   double kurtosis_computation(int variable) const;

   // acces membres de la classe
   --     FrequencyDistribution* get_marginal(int variable) const { return marginal[variable]; }
   --double get_covariance(int variable1, int variable2) const { return covariance[variable1][variable2]; }
   */

}

// Vectors

class VectorDistanceWrap
{

public:

  static boost::shared_ptr<VectorDistance>
  build_from_types(boost::python::list& types, boost::python::list& weigths,
      int distance_type)
  {
    VectorDistance* dist;
    int nb_variable;

    nb_variable = boost::python::len(types);

    stat_tool::wrap_util::auto_ptr_array<double> variable_weight(
        new double[nb_variable]);

    stat_tool::wrap_util::auto_ptr_array<int> variable_type(
        new int[nb_variable]);

    // Extract each element of the vector
    for (int i = 0; i < nb_variable; i++)
      {
        variable_type[i] = boost::python::extract<int>(types[i]);
        variable_weight[i] = boost::python::extract<double>(weigths[i]);
      }

    dist = new VectorDistance(nb_variable, variable_type.get(),
        variable_weight.get(), distance_type);

    return boost::shared_ptr<VectorDistance>(dist);
  }

  static boost::shared_ptr<VectorDistance>
  read_from_file(char* filename)
  {
    VectorDistance* dist;
    StatError error;
    dist = vector_distance_ascii_read(error, filename);

    if (!dist)
      stat_tool::wrap_util::throw_error(error);

    return boost::shared_ptr<VectorDistance>(dist);
  }

  static void
  file_ascii_write(const VectorDistance& d, const char* path, bool exhaustive)
  {
    bool result = true;
    StatError error;

    result = d.ascii_write(error, path, exhaustive);
    if (!result)
      stat_tool::wrap_util::throw_error(error);

  }
  WRAP_METHOD_SPREADSHEET_WRITE(DistanceMatrix);

};



void class_vectordistance()
{

  class_< VectorDistance, bases< StatInterface > >
  ("_VectorDistance", "Distance class", init<>())
  // constructors
  .def("__init__", make_constructor(VectorDistanceWrap::build_from_types))
  .def("__init__", make_constructor(VectorDistanceWrap::read_from_file))

  .def(self_ns::str(self))
  .def("__len__", &VectorDistance::get_nb_variable)
  .def("file_ascii_write", VectorDistanceWrap::file_ascii_write, "Save vector distance summary into a file")

  .add_property("nb_variable", &VectorDistance::get_nb_variable, "returns number of variable")
  .add_property("distance_type", &VectorDistance::get_distance_type, "returns distance type")
  //   .def("spreadsheet_write", VectorDistanceWrap::spreadsheet_write, "Save data into CSV file")
    ;

  /*
  VectorDistance();
    VectorDistance(int inb_variable , int idistance_type , int *ivariable_type ,
                   double *iweight , int *inb_value , double ***isymbol_distance ,
                   int *iperiod);


    // fonctions pour la compatibilite avec la classe StatInterface

    double* max_symbol_distance_computation(int variable) const;

    void dispersion_computation(int variable , const FrequencyDistribution *marginal , double *rank = 0) const;

    // acces membres de la classe

    int get_variable_type(int variable) const { return variable_type[variable]; }
    double get_weight(int variable) const { return weight[variable]; }
    double get_dispersion(int variable) const { return dispersion[variable]; }
    int get_nb_value(int variable) const { return nb_value[variable]; }
    double get_symbol_distance(int variable , int symbol1 , int symbol2) const    { return symbol_distance[variable][symbol1][symbol2]; }
    int get_period(int variable) const { return period[variable]; }
    */

}

