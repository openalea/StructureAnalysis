/*------------------------------------------------------------------------------
 *
 *
 *        VPlants.Sequence_analysis : VPlants Statistics module
 *
 *        Copyright 2006-2007 INRIA - CIRAD - INRA
 *
 *        File author(s): Yann Guédon <yann.guedon@cirad.fr>
 *                        Thomas Cokelaer <Thomas.Cokelaer@inria.fr>
 *
 *        Distributed under the GPL 2.0 License.
 *        See accompanying file LICENSE.txt or copy at
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr
 *
 *        $Id:  $
 *
 *-----------------------------------------------------------------------------*/

#include "wrapper_util.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/vectors.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequence_analysis/renewal.h"
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/tops.h"
#include "sequence_analysis/sequence_label.h"

#include <boost/python.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/list.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>

#include "boost_python_aliases.h"

using namespace boost::python;
using namespace boost;
//using namespace stat_tool;

class SequencesWrap
{

public:

  static boost::shared_ptr<Sequences>
  sequences_from_file(char* filename)
  {
    Format_error error;
    Sequences *sequences = NULL;
    bool old_format = false;
    sequences = sequences_ascii_read(error, filename, old_format);
    return boost::shared_ptr<Sequences>(sequences);
  }

  static boost::shared_ptr<Sequences>
  sequences_from_file_old(char* filename, bool old_format)
  {
    Format_error error;
    Sequences *sequences = NULL;
    sequences = sequences_ascii_read(error, filename, old_format);
    return boost::shared_ptr<Sequences>(sequences);
  }

  static Sequences*
  build_from_lists(boost::python::list& input_list,
      boost::python::list& input_identifiers, int input_index_parameter_type)
  {

    // case 1 list of n lists of floats (only 1 variable and different sizes possible):
    //			seq = Sequences([ [1,1,1], [2,2,2,2,2]])
    //
    // case 2 list of n lists (different sizes) of variables (same size)
    //			seq = Sequences([ [ [1,1,1], [11,11,11] ],
    //							  [ [2,2,2,2,2], [22,22,22,22,22] ]
    //  						])

    // the length of the main lists to get the number of sequences
    int nb_sequences = boost::python::len(input_list);
    int nb_identifiers = boost::python::len(input_identifiers);

    //check in the python code anyway
    if (nb_sequences != nb_identifiers)
      {
        PyErr_SetString(PyExc_IndexError,
            "number of identifiers must be equal to number of sequences");
        boost::python::throw_error_already_set();
      }

    Sequences *sequences = 0;
    int nb_variables = 0;
    int nb_variables_check = 0;
    int *length;
    int *identifiers;
    //int ***int_sequence;
    double ***real_sequence;
    int ***int_sequence;
    bool is_float = false;
    bool is_int = false;
    bool is_sequence = false;
    boost::python::list sequence;
    int index_parameter_type = -1;
    int type = -1; //

    index_parameter_type = input_index_parameter_type;
    // length of each sequence will be store in this variable
    length = new int[nb_sequences];
    identifiers = new int[nb_sequences];

    //extract<boost::python::list> get_list(input_list[0]);
    boost::python::list seq = extract<boost::python::list> (input_list[0]);
    extract<boost::python::list> get_list(seq[0]);

    if (!get_list.check())
      {
        boost::python::list seq0 = extract<boost::python::list> (input_list[0]);
        object elt = seq0[0];
        extract<float> get_float(elt);
        extract<int> get_int(elt);

        if (get_int.check())
          is_int = true;
        else if (get_float.check())
          is_float = true;
      }
    else
      {

        boost::python::list sequence = extract<boost::python::list> (
            input_list[0]);
        boost::python::list variable = extract<boost::python::list> (
            sequence[0]);

        extract<float> get_float(variable[0]);
        extract<int> get_int(variable[0]);

        if (get_int.check())
          is_int = true;
        else if (get_float.check())
          is_float = true;
      }

    // allocate memory given the number of sequences and the type
    if (is_float)
      {
        real_sequence = new double**[nb_sequences];
        type = REAL_VALUE;

      }
    if (is_int)
      {
        int_sequence = new int**[nb_sequences];
        type = INT_VALUE;
      }

    // for each sequence, are we considering case 1 or 2 ?
    for (int seqi = 0; seqi < nb_sequences; seqi++)
      {
        // whatever case it is, the identifiers can be set here
        identifiers[seqi] = extract<int> (input_identifiers[seqi]);

        // try to get a single sequence
        boost::python::list sequence = extract<boost::python::list> (
            input_list[seqi]);

        // is it case 2 i.e. there is another nested list ?
        extract<boost::python::list> get_list(sequence[0]);

        // if not a list, we are in the case 1
        if (!get_list.check())
          {
            // the length of the sequence is simply:
            int N = len(sequence);
            length[seqi] = N;

            // and there is only 1 variable
            nb_variables = 1; // used by the constructor

            // so the data structure is as follows
            if (is_float)
              {
                real_sequence[seqi] = new double*[1];
                real_sequence[seqi][0] = new double[N];
              }
            else
              {
                int_sequence[seqi] = new int*[1];
                int_sequence[seqi][0] = new int[N];
              }

            // and we can populate the output using the sequence as vector
            for (int j = 0; j < N; j++)
              {
                if (is_float)
                  {
                    double v = extract<float> (sequence[j]);
                    real_sequence[seqi][0][j] = v;
                  }
                else
                  {
                    int v = extract<int> (sequence[j]);
                    int_sequence[seqi][0][j] = v;
                  }
              }
          }
        // if we can get the data, we have a list and therefore we are in case 2
        else
          {
            nb_variables = len(sequence); // used by the constructor

            // what are the number of variables in this constructor :
            int N_var = len(sequence);

            // loop over the variable to extract their contents
            for (int vari = 0; vari < N_var; vari++)
              {
                // get the variable of the current sequence
                boost::python::list variable = extract<boost::python::list> (
                    sequence[vari]);

                // and update the data structure
                int N = len(variable);
                length[seqi] = N;

                //real_sequence[seqi][vari] = new double[N];

                // before extracting the vector of the

                if (is_float)
                  {
                    if (vari == 0)
                      real_sequence[seqi] = new double*[N_var]; // only once
                    real_sequence[seqi][vari] = new double[N];
                  }
                if (is_int)
                  {
                    if (vari == 0)
                      int_sequence[seqi] = new int*[N_var]; // only once
                    int_sequence[seqi][vari] = new int[N];
                  }

                // get the vector corresponding to a sequence/variable pair
                for (int j = 0; j < N; j++)
                  {
                    if (is_float)
                      {
                        float v = extract<float> (variable[j]);
                        real_sequence[seqi][vari][j] = v;
                      }

                    if (is_int)
                      {
                        int v = extract<int> (variable[j]);
                        int_sequence[seqi][vari][j] = v;

                      }

                  }
              }
          }
      }

    if (is_float)
      {
        sequences = new Sequences(nb_sequences, identifiers, length,
            nb_variables, real_sequence);

        if (real_sequence)
          {
            for (int i = 0; i < nb_sequences; i++)
              {
                for (int j = 0; j < nb_variables; j++)
                  {
                    delete[] real_sequence[i][j];
                  }
                delete[] real_sequence[i];
              }
            delete[] real_sequence;
          }
      }

    if (is_int)
      {
        sequences = new Sequences(nb_sequences, identifiers, length,
            index_parameter_type, nb_variables, type, int_sequence);
        // freeing memory

        if (int_sequence)
          {
            for (int i = 0; i < nb_sequences; i++)
              {
                for (int j = 0; j < nb_variables; j++)
                  {
                    delete[] int_sequence[i][j];
                  }
                delete[] int_sequence[i];
              }

            delete[] int_sequence;
          }
      }

    delete[] length;
    delete[] identifiers;

    return sequences;

  }

  static Distribution_data*
  extract_value(const Sequences& input, int variable)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, extract, Distribution_data, variable);
  }

  static boost::python::list
  get_identifiers(const Sequences& seq)
  {
    boost::python::list l;

    int nb_seq = seq.get_nb_sequence();

    for (int s = 0; s < nb_seq; s++)
      {
        l.append(seq.get_identifier(s));
      }

    return l;
  }

  static std::string
  ascii_data_write(const Sequences& d, bool exhaustive)
  {
    std::stringstream os;
    std::string res;

    d.ascii_data_write(os, exhaustive, true);
    res = os.str();
    return res;
  }

  static void
  file_ascii_data_write(const Sequences& d, const char* path, bool exhaustive)
  {
    bool result = true;
    Format_error error;

    result = d.ascii_data_write(error, path, exhaustive);
    if (!result)
      sequence_analysis::wrap_util::throw_error(error);

  }

  static void
  file_ascii_write(const Sequences& d, const char* path, bool exhaustive)
  {
    bool result = true;
    Format_error error;

    result = d.ascii_write(error, path, exhaustive);
    if (!result)
      sequence_analysis::wrap_util::throw_error(error);

  }

  static boost::python::list
  get_item(const Sequences* seq, boost::python::tuple indexes)
  {
    int index_var = extract<int> (indexes[1]);
    int index_seq = extract<int> (indexes[0]);

    // Test index
    if (index_seq < 0 || index_seq >= seq->get_nb_sequence())
      {
        PyErr_SetString(PyExc_IndexError, "sequence index out of bound");
        boost::python::throw_error_already_set();
      }
    if (index_var < 0 || index_var >= seq->get_nb_variable())
      {
        PyErr_SetString(PyExc_IndexError, "variable index out of bound");
        boost::python::throw_error_already_set();
      }

    boost::python::list l;

    int nb_length = seq->get_length(index_seq);

    for (int index = 0; index < nb_length; index++)
      {
        if ((seq->get_type(index_var) == INT_VALUE)
            || (seq->get_type(index_var) == STATE))
          l.append(seq->get_int_sequence(index_seq, index_var, index));
        else
          l.append(seq->get_real_sequence(index_seq, index_var, index));
      }

    return l;
  }

  static Sequences*
  value_select(const Sequences& seq, int variable, const object& min,
      const object& max, bool keep)
  {
    HEADER(Sequences);
    std::stringstream s;
 
    boost::python::extract<int> get_min(min);
    boost::python::extract<int> get_max(max);

    if (get_min.check() && get_max.check()) // Array of int
      {
        int mi = get_min();
        int ma = get_max();
        ret = seq.value_select(error, s, variable, mi, ma, keep);
      }
    else
      {
        double mi = extract<double> (min);
        double ma = extract<double> (max);
        ret = seq.value_select(error, s, variable, mi, ma, keep);
      }

    cout << s.str() << endl;
    FOOTER;
  }

  static Sequences*
  select_variable(const Sequences& seq, const boost::python::list& variables,
      bool keep)
  {

    CREATE_ARRAY(variables, int, vars);

    SIMPLE_METHOD_TEMPLATE_1(seq, select_variable, Sequences, vars_size,
        vars.get(), keep);

  }

  static Sequences*
  select_individual(const Sequences& seq,
      const boost::python::list& identifiers, bool keep)
  {
    CREATE_ARRAY(identifiers, int, ids);

    SIMPLE_METHOD_TEMPLATE_1(seq, select_individual, Sequences, ids_size,
        ids.get(), keep);
  }

  // Shift
  static Sequences*
  shift(const Sequences& seq, int var, double param)
  {
    HEADER(Sequences);

    if (seq.get_type(var - 1) == REAL_VALUE)
      ret = seq.shift(error, var, (double) param);
    else
      ret = seq.shift(error, var, (int) param);

    FOOTER(ret);
  }

  // Merge
  static Sequences*
  merge(const Sequences& input_seq, const boost::python::list& seqs)
  {

    CREATE_ARRAY(seqs, const Sequences *, sequens)

    SIMPLE_METHOD_TEMPLATE_1(input_seq, merge, Sequences, sequens_size, sequens.get());
  }

  static Sequences*
  merge_variable(const Sequences& input_seq, const boost::python::list& seqs,
      int ref_sample)
  {

    CREATE_ARRAY(seqs, const Sequences *, data);

    SIMPLE_METHOD_TEMPLATE_1(input_seq, merge_variable, Sequences,data_size,
        data.get(), ref_sample);
  }

  // Cluster
  static Sequences*
  cluster_step(const Sequences& seq, int variable, int step)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, cluster, Sequences, variable, step);
  }

  static Sequences*
  cluster_limit(const Sequences& seq, int variable, boost::python::list& limit)
  {

    Format_error error;

    int nb_limit = len(limit);
    bool is_float = true;
    int *lint = NULL;
    double *ldouble = NULL;
    Sequences* ret;

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
          ldouble[i] = extract<double> (limit[i]);
        else
          lint[i] = extract<int> (limit[i]);
      }

    // Call correct function
    if (is_float)
      {
        ret = seq.cluster(error, variable, nb_limit, ldouble);
        delete[] ldouble;
      }
    else
      {
        ret = seq.cluster(error, variable, nb_limit, lint);
        delete[] lint;
      }

    FOOTER;
  }

  static Sequences*
  transcode(const Sequences& seq, int variable, boost::python::list& symbol)
  {

    int nb_symbol = len(symbol);
    sequence_analysis::wrap_util::auto_ptr_array<int> l(new int[nb_symbol]);

    int expected_nb_symbol = (int) (seq.get_max_value(variable - 1)
        - seq.get_min_value(variable - 1)) + 1;

    if (nb_symbol != expected_nb_symbol)
      sequence_analysis::wrap_util::throw_error("Bad number of Symbol");

    for (int i = 0; i < nb_symbol; i++)
      l[i] = extract<int> (symbol[i]);

    SIMPLE_METHOD_TEMPLATE_1(seq, transcode, Sequences, variable, l.get());
  }

  static Sequences*
  reverse(const Sequences& seq)
  {
    SIMPLE_METHOD_TEMPLATE_0(seq, reverse, Sequences);
  }

  static Sequences*
  length_select(const Sequences& input,  int min_length, int max_length, bool keep)
  {
     Format_error error; 
    Sequences* ret;
    std::ostringstream os;

    ret = input.length_select(error, os,
		min_length, max_length, keep);
    if (!ret) 
      sequence_analysis::wrap_util::throw_error(error);
    cout << os.str() << endl;
    return ret;
  }

  static Sequences*
  remove_run(const Sequences& seq, int variable, int ivalue, char position,
      int max_run_length)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, remove_run, Sequences, variable, ivalue,
        position, max_run_length);
  }

  static int
  get_type(const Sequences& seq, int index)
  {
    if (index < 0 || index >= seq.get_nb_variable())
      {
        PyErr_SetString(PyExc_IndexError,
            "index must be positive and less than number of variables");
        boost::python::throw_error_already_set();
      }

    return seq.get_type(index);
  }

  static int
  get_length(const Sequences& seq, int index)
  {
    if (index < 0 || index >= seq.get_nb_sequence())
      {
        PyErr_SetString(PyExc_IndexError,
            "index must be positive and less than number of sequences");
        boost::python::throw_error_already_set();
      }

    return seq.get_length(index);
  }

  static double
  get_min_value(const Sequences& seq, int variable)
  {
    if (variable < 0 || variable >= seq.get_nb_variable())
      {
        PyErr_SetString(PyExc_IndexError,
            "index must be positive and less than number of variables");
        boost::python::throw_error_already_set();
      }
    return seq.get_min_value(variable);
  }

  static double
  get_max_value(const Sequences& seq, int variable)
  {
    if (variable < 0 || variable >= seq.get_nb_variable())
      {
        PyErr_SetString(PyExc_IndexError,
            "index must be positive and less than number of variables");
        boost::python::throw_error_already_set();
      }
    return seq.get_max_value(variable);
  }

  static Sequences*
  scaling(const Sequences& seq, int variable, int scaling_coeff)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, scaling, Sequences, variable, scaling_coeff);
  }

  static Sequences*
  round(const Sequences& seq, int variable, int mode)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, round, Sequences, variable, mode);
  }

  //cumulate
  static Sequences*
  cumulate(const Sequences& seq, int variable)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, cumulate, Sequences, variable);
  }

  static Sequences*
  difference(const Sequences& seq, int variable, bool first)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, difference, Sequences, variable, first);
  }

  static Sequences*
  index_parameter_select(const Sequences& input, int min_index_parameter,
      int max_index_parameter, bool keep)
  {
    Format_error error; 
    Sequences* ret;
    std::ostringstream os;
    
    ret = input.index_parameter_select(error,  os,
		min_index_parameter, max_index_parameter, keep);
    if (!ret) 
      sequence_analysis::wrap_util::throw_error(error);
    cout << os.str() << endl;
    return ret;

  }

  static Sequences*
  remove_index_parameter(const Sequences& seq)
  {
    SIMPLE_METHOD_TEMPLATE_0(seq, remove_index_parameter, Sequences);
  }

  static Sequences*
  index_parameter_extract(const Sequences& seq, int min_index_parameter,
      int max_index_parameter)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, index_parameter_extract, Sequences,
        min_index_parameter, max_index_parameter);
  }

  static Sequences*
  segmentation_extract(const Sequences& seq, int variable,
      boost::python::list& input_values, bool keep)
  {
    CREATE_ARRAY(input_values, int, values);
    SIMPLE_METHOD_TEMPLATE_1(seq, segmentation_extract, Sequences,
        variable, values_size, values.get(), keep);
  }

  static Sequences*
  moving_average(const Sequences& seq, boost::python::list& input_values,
      int variable, bool begin_end, int output)
  {

    int nb_value = len(input_values);
    double *values;
    values = new double[nb_value * 2 + 1];
    for (int i = 0; i < nb_value; i++)
      {
        values[i] = extract<double> (input_values[i]);
        values[2 * nb_value - i] = values[i];
      }
    values[nb_value] = 0; //todo check that!!

    SIMPLE_METHOD_TEMPLATE_1(seq, moving_average, Sequences,
            nb_value, values, variable, begin_end,  output);

  }

  static Sequences*
  moving_average_from_distribution(const Sequences& seq,
      const Distribution& dist, int variable, bool begin_end, int output)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, moving_average, Sequences, dist, variable,
        begin_end, output);
  }

  static Sequences*
  pointwise_average(const Sequences& seq, bool standard_deviation, int output,
      const char *path, char format)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, pointwise_average, Sequences,
        standard_deviation, output, path, format);
  }

  static Sequences*
  recurrence_time_sequences(const Sequences& seq, int variable, int value)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, recurrence_time_sequences, Sequences,
        variable, value);
  }

  static Sequences*
  cross(const Sequences& seq)
  {
    SIMPLE_METHOD_TEMPLATE_0(seq, cross, Sequences);
  }

  static Sequences*
  transform_position(const Sequences& seq, int step)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, transform_position, Sequences, step);
  }

  static Sequences*
  sojourn_time_sequences(const Sequences& seq, int variable)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, sojourn_time_sequences, Sequences, variable);
  }

  static Correlation*
  partial_autocorrelation_computation(const Sequences& seq, int variable,
      int itype, int max_lag)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, partial_autocorrelation_computation,
        Correlation, variable, itype, max_lag);
  }

  static Vectors*
  extract_vectors(const Sequences& seq, int feature_type, int variable, int value)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, extract_vectors, Vectors, feature_type,
        variable, value);
  }

  static Distribution_data*
  extract_length(const Sequences& input)
  {
    Distribution_data *res;
    res = new Distribution_data(*(input.get_hlength()));
    return res;
  }

  static Markovian_sequences*
  markovian_sequences(const Sequences& input)
  {
    SIMPLE_METHOD_TEMPLATE_0(input, markovian_sequences, Markovian_sequences);
  }

  static Correlation*
  correlation_computation(const Sequences& input, int variable1, int variable2,
  int itype, int max_lag, int normalization)
  {
    SIMPLE_METHOD_TEMPLATE_1 (input, correlation_computation,
      Correlation, variable1, variable2, itype, max_lag, normalization)
  }

  static Distance_matrix*
  alignment_vector_distance(const Sequences& input,
      const Vector_distance &ivector_dist,  int ref_identifier,
      int test_identifier, bool begin_free, bool end_free, int indel_cost,
      double indel_factor, bool transposition_flag,
      double transposition_factor, const char *result_path, char result_format,
      const char *alignment_path, char alignment_format)
  {
    HEADER_OS(Distance_matrix);

    ret = input.alignment(error, &os, ivector_dist, ref_identifier, test_identifier, begin_free,
        end_free, indel_cost, indel_factor, transposition_flag,
        transposition_factor, result_path, result_format, alignment_path,
        alignment_format );

   FOOTER_OS;
  }

  static Distance_matrix*
  alignment(const Sequences& input,int ref_identifier,
  int test_identifier, bool begin_free, bool end_free,
  const char *result_path, char result_format,
  const char *alignment_path, char alignment_format)
  {
    HEADER_OS(Distance_matrix);

    ret = input.alignment(error, &os, ref_identifier, test_identifier, begin_free,
      end_free, result_path, result_format, alignment_path,
      alignment_format );
    FOOTER_OS;
  }

  static Renewal_data*
  extract_renewal_data(const Sequences & input,
    int variable , int begin_index , int end_index)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, extract_renewal_data, Renewal_data,
      variable, begin_index, end_index);
  }

  static Time_events*
  extract_time_events(const Sequences & input,
    int variable , int begin_date , int end_date ,
    int previous_date, int next_date)

  {
    SIMPLE_METHOD_TEMPLATE_1(input, extract_time_events, Time_events,
      variable, begin_date, end_date, previous_date, next_date);

  }

  static Vectors*
  build_vectors(const Sequences & input,
  bool index_variable)
  {
    Vectors* ret;

    ret = input.build_vectors(index_variable);

    return ret;
  }

  static Sequences*
  segmentation_change_point(const Sequences &input, int iidentifier,
    int nb_segment , boost::python::list input_change_point ,
    boost::python::list input_model_type , int output)
  {
    HEADER_OS(Sequences);

    int nb = len(input_change_point);
    int *change_point;
    change_point = new int[nb];
    for (int i = 0; i < nb; i++)
      change_point[i] = boost::python::extract<int> (input_change_point[i]);

    nb = len(input_model_type);
    int *model_type;
    model_type = new int[nb];
    for (int i = 0; i < nb; i++)
      model_type[i] = extract<int> (input_model_type[i]);

    ret = input.segmentation(error, os, iidentifier, nb_segment,
      change_point, model_type, output);

    FOOTER_OS;
  }


  static Sequences*
  segmentation_array(const Sequences &input, boost::python::list input_nb_segment,
      boost::python::list input_model_type , int iidentifier,  int output)
  {
    HEADER_OS(Sequences);

    int nb = len(input_nb_segment);
    int *nb_segment;
    nb_segment = new int[nb];
    for (int i = 0; i < nb; i++)
        nb_segment[i] = extract<int> (input_nb_segment[i]);

    nb = len(input_model_type);
    int *model_type;
    model_type = new int[nb];
    for (int i = 0; i < nb; i++)
        model_type[i] = extract<int> (input_model_type[i]);


    ret = input.segmentation(error, os, nb_segment, model_type, iidentifier,
        output);
    FOOTER_OS;
  }

  static Sequences*
  segmentation_model(const Sequences &input, int iidentifier,
    int max_nb_segment , boost::python::list input_model_type)
  {
    HEADER_OS(Sequences);

    int nb = len(input_model_type);
    int *model_type;
    model_type = new int[nb];
    for (int i = 0; i < nb; i++)
        model_type[i] = extract<int> (input_model_type[i]);

    ret = input.segmentation(error, os, iidentifier, max_nb_segment, model_type);
    FOOTER_OS;
  }

  static Sequences*
  segmentation_vector_distance(const Sequences &input,  int iidentifier,
    int nb_segment ,  const Vector_distance &ivector_dist, int output)
  {
    HEADER_OS(Sequences);
    ret = input.segmentation(error,iidentifier,  nb_segment,
      ivector_dist, os, output);
    FOOTER_OS;
  } 



};

// Boost declaration

void
class_sequences()
{

  enum_<sequence_analysis::wrap_util::UniqueInt<5, 100> > ("IndexParameterType")
    .value("IMPLICIT_TYPE", IMPLICIT_TYPE)
    .value("TIME", TIME)
    .value( "TIME_INTERVAL", TIME_INTERVAL)
    .value("POSITION", POSITION)
    .value("POSITION_INTERVAL", POSITION_INTERVAL)
    .export_values();

  enum_<sequence_analysis::wrap_util::UniqueInt<6, 101> > ("Type")
    .value("INT_VALUE", INT_VALUE)
    .value("REAL_VALUE", REAL_VALUE)
    .value("STATE", STATE)
    .value("OLD_INT_VALUE", OLD_INT_VALUE)
    .value("NB_INTERNODE", NB_INTERNODE)
    .value("AUXILIARY", AUXILIARY)
    .export_values();

  enum_<sequence_analysis::wrap_util::UniqueInt<4, 102> > ("Algorithm")
    .value( "CTM_BIC", CTM_BIC)
    .value("CTM_KT", CTM_KT)
    .value("LOCAL_BIC",LOCAL_BIC)
    .value("CONTEXT", CONTEXT)
    .export_values();

  enum_<sequence_analysis::wrap_util::UniqueInt<4, 103> > ("Estimator") .value(
      "MAXIMUM_LIKELIHOOD", MAXIMUM_LIKELIHOOD) .value("LAPLACE", LAPLACE) .value(
      "ADAPTATIVE_LAPLACE", ADAPTATIVE_LAPLACE) .value("UNIFORM_SUBSET",
      UNIFORM_SUBSET) .value("UNIFORM_CARDINALITY", UNIFORM_CARDINALITY) .export_values();

  enum_<sequence_analysis::wrap_util::UniqueInt<11, 104> > (
      "MarkovianSequenceType") .value("OBSERVATION", OBSERVATION) .value(
      "FIRST_OCCURRENCE", FIRST_OCCURRENCE) .value("RECURRENCE_TIME",
      RECURRENCE_TIME) .value("SOJOURN_TIME", SOJOURN_TIME) .value(
      "INITIAL_RUN", INITIAL_RUN) .value("FINAL_RUN", FINAL_RUN) .value(
      "NB_RUN", NB_RUN) .value("NB_OCCURRENCE", NB_OCCURRENCE) .value("LENGTH",
      LENGTH) .value("SEQUENCE_CUMUL", SEQUENCE_CUMUL) .value("SEQUENCE_MEAN",
      SEQUENCE_MEAN) .export_values();

  enum_<sequence_analysis::wrap_util::UniqueInt<2, 105> > ("IndelCost") .value(
      "ADAPTATIVE", ADAPTATIVE) .value("FIXED", FIXED) .export_values();

  class_<Sequences, bases<STAT_interface> > ("_Sequences", "Sequences")
    .def("__init__", make_constructor(SequencesWrap::sequences_from_file))
    .def("__init__", make_constructor(SequencesWrap::sequences_from_file_old))
    .def("__init__", make_constructor(SequencesWrap::build_from_lists))
    .def(init <const Renewal_data&>())

    // Python Operators

   .def(self_ns::str(self)) //__str__
   .def("__len__", &Sequences::get_nb_sequence,"Returns number of sequences")
   .def("__getitem__", SequencesWrap::get_item)

   //property
   .add_property("nb_sequence", &Sequences::get_nb_sequence, "Return the number of sequences")
   .add_property("nb_variable", &Sequences::get_nb_variable, "Return the number of variables")
   .add_property("max_length", &Sequences::get_max_length,"Return max length")
   .add_property("cumul_length", &Sequences::get_cumul_length,"Return cumul length")

   //index arguments, return single value
   .def("get_index_parameter_type", &Sequences::get_index_parameter_type,"return index parameter type")
   .def("get_min_value", &Sequences::get_min_value,args("index_var"), "return min value of variables")
   .def("get_max_value", &Sequences::get_max_value,args("index_var"), "return max value of variables")

   // index arguments, wrapping required
   .def("get_length", &SequencesWrap::get_length,args("index_seq"), "return length of a sequence")
   .def("get_type", &SequencesWrap::get_type,args("index_var"),"return  type")
   .def("get_identifiers", &SequencesWrap::get_identifiers,"returns list of identifiers")
   .def("ascii_data_write", SequencesWrap::ascii_data_write,"Return a string with the object representation")
   .def("file_ascii_write", SequencesWrap::file_ascii_write,"Save vector summary into a file")
   .def("file_ascii_data_write", SequencesWrap::file_ascii_data_write,"Save vector data into a file")


   DEF_RETURN_VALUE("alignment_vector_distance", SequencesWrap::alignment_vector_distance, args(""), "todo")
   DEF_RETURN_VALUE("alignment", SequencesWrap::alignment, args(""), "todo")
   DEF_RETURN_VALUE("build_vectors", SequencesWrap::build_vectors, args("index_variable"), "build a vector from sequence")
   DEF_RETURN_VALUE("cluster_step", SequencesWrap::cluster_step, args("variable", "step"),"Cluster Step")
   DEF_RETURN_VALUE("cluster_limit", SequencesWrap::cluster_limit, args("variable", "limits"),"Cluster limit")
   DEF_RETURN_VALUE("correlation_computation", SequencesWrap::correlation_computation, args("variable1", "variable2", "type","max_lag", "normalization"),"compute correlation")
   DEF_RETURN_VALUE("difference", SequencesWrap::difference,args("variable", "first_element"),"Difference")
   DEF_RETURN_VALUE("extract_value", SequencesWrap::extract_value, args("variable"),"Extract histogram")
   DEF_RETURN_VALUE("extract_renewal_data", SequencesWrap::extract_renewal_data, args("todo"),"Extract renewal_data")
   DEF_RETURN_VALUE("extract_time_events", SequencesWrap::extract_time_events, args("todo"),"Extract time events")
   DEF_RETURN_VALUE("index_parameter_select", SequencesWrap::index_parameter_select,args("min_index_parameter", "max_index_parameter", "keep"),"Select sequences in an index parameter range")
   DEF_RETURN_VALUE("index_parameter_extract", SequencesWrap::index_parameter_extract,args("min_index_parameter", "max_index_parameter"),"Select sequences in an index parameter range")
   DEF_RETURN_VALUE("length_select",SequencesWrap::length_select, args("min_length", "max_length", "keep"),"see LengthSelect")
   DEF_RETURN_VALUE("moving_average", SequencesWrap::moving_average,args("nb_point" ,"filter" , "variable" , "begin_end" , "output"),"Moving average from an array")
   DEF_RETURN_VALUE("moving_average_from_distribution", SequencesWrap::moving_average_from_distribution,args("nb_point" ,"filter" , "variable" , "begin_end" , "output"),"Moving average fron distribution ")
   DEF_RETURN_VALUE("partial_autocorrelation_computation", SequencesWrap::partial_autocorrelation_computation, args("variable", "itype", "max_lag"),"Transcode")
   DEF_RETURN_VALUE("pointwise_average", SequencesWrap::pointwise_average, args("standard_deviation", "output", "path", "format"), "Pointwise average")
   DEF_RETURN_VALUE("recurrence_time_sequences", SequencesWrap::recurrence_time_sequences,args("variable", "value"),"Recurrence time sequences")
   DEF_RETURN_VALUE("remove_run",SequencesWrap::remove_run, args("variable","ivalue", "position", "max_run_length"), "see RemoveRun")
   DEF_RETURN_VALUE("round", SequencesWrap::round, args("variable", "mode"),"round variable")
   DEF_RETURN_VALUE("scaling", SequencesWrap::scaling, args("variable", "scaling_coeff"),"scaling variable")
   DEF_RETURN_VALUE("select_variable", SequencesWrap::select_variable, args("variables", "keep"), "select variable given a list of index")
   DEF_RETURN_VALUE("select_individual", SequencesWrap::select_individual, args("identifiers", "keep"),"Select individuals given a list of identifiers")
   DEF_RETURN_VALUE("segmentation_extract",SequencesWrap::segmentation_extract, args("variable", "nb_values", "values", "keep"),"Segmentation extract")
   DEF_RETURN_VALUE("segmentation_model", SequencesWrap::segmentation_model, args("todo"), "todo")
   DEF_RETURN_VALUE("segmentation_change_point", SequencesWrap::segmentation_change_point, args("todo"), "todo")
   DEF_RETURN_VALUE("segmentation_array", SequencesWrap::segmentation_array, args("todo"), "todo")
   DEF_RETURN_VALUE("segmentation_vector_distance", SequencesWrap::segmentation_vector_distance, args("todo"), "todo")
   DEF_RETURN_VALUE("shift", SequencesWrap::shift, args("variable","param"),"Shift")
   DEF_RETURN_VALUE("sojourn_time_sequences", SequencesWrap::sojourn_time_sequences, args("variable"), "Sojourn time sequences")
   DEF_RETURN_VALUE("transcode", SequencesWrap::transcode, args("variable", "symbols"),"Transcode")
   DEF_RETURN_VALUE("transform_position", SequencesWrap::transform_position, args("step"), "Transform position")
   DEF_RETURN_VALUE("value_select", SequencesWrap::value_select,args("variable", "min", "max", "keep"),"Selection of individuals according to the values taken by a variable")

   DEF_RETURN_VALUE("extract_vectors", SequencesWrap::extract_vectors, args("type","variable","value"), "extract vectors")

   DEF_RETURN_VALUE_NO_ARGS("extract_length", SequencesWrap::extract_length, "extract length of the sequences and returns a vector")
   DEF_RETURN_VALUE_NO_ARGS("remove_index_parameter", SequencesWrap::remove_index_parameter,"Remove index parameter")
   DEF_RETURN_VALUE_NO_ARGS("reverse", SequencesWrap::reverse,"reverse")
   DEF_RETURN_VALUE_NO_ARGS("cross", SequencesWrap::cross, "Cross")
   DEF_RETURN_VALUE_NO_ARGS("cumulate", SequencesWrap::cumulate,"Cumulate")
   DEF_RETURN_VALUE_NO_ARGS("merge", SequencesWrap::merge, "Merge sequences")
   DEF_RETURN_VALUE_NO_ARGS("merge_variable", SequencesWrap::merge_variable, "Merge variables")
   DEF_RETURN_VALUE_NO_ARGS("markovian_sequences", SequencesWrap::markovian_sequences , "returns markovian sequence")
   DEF_RETURN_VALUE_NO_ARGS("extract_sequence_length", SequencesWrap::extract_length , "todo")

;

/*
 Sequences(int inb_sequence , int inb_variable);
 Sequences(int inb_sequence , int *iidentifier , int *ilength , int iindex_parameter_type ,	              int inb_variable , int *itype , bool init_flag = false)	    { init(inb_sequence , iidentifier , ilength , iindex_parameter_type , inb_variable ,	           itype , init_flag); }
 Sequences(int inb_sequence , int *iidentifier , int *ilength , int inb_variable , int *itype , bool init_flag = false){ init(inb_sequence , iidentifier , ilength , IMPLICIT_TYPE , inb_variable ,	           itype , init_flag); }
 Sequences(int inb_sequence , int *iidentifier , int *ilength , int inb_variable ,  bool init_flag = false)	    { init(inb_sequence , iidentifier , ilength , inb_variable , init_flag); }
 Sequences(const Histogram &ihlength , int inb_variable , bool init_flag = false);
 Sequences(const Sequences &seq , int inb_sequence , int *index);
 Sequences(const Sequences &seq , bool *segment_mean);
 Sequences(const Sequences &seq , char transform = 'c' , int param = DEFAULT);
 Sequences(const Histogram &ihlength , int inb_variable , bool init_flag = false);

 Tops* tops(Format_error &error) const;
 bool check(Format_error &error , const char *pattern_label);

 std::ostream& line_write(std::ostream &os) const;
 bool plot_data_write(Format_error &error , const char *prefix , const char *title = 0) const;
 bool spreadsheet_write(Format_error &error , const char *path) const;
 bool plot_write(Format_error &error , const char *prefix ,   const char *title = 0) const;

 int min_index_parameter_computation() const;
 int max_index_parameter_computation(bool last_position = false) const;

 void marginal_histogram_computation(int variable);
 double mean_computation(int variable) const;
 double variance_computation(int variable , double mean) const;
 double mean_absolute_deviation_computation(int variable , double mean) const;
 double mean_absolute_difference_computation(int variable) const;
 double skewness_computation(int variable , double mean , double variance) const;
 double kurtosis_computation(int variable , double mean , double variance) const;

 Histogram* value_index_interval_computation(Format_error &error , int variable , int value) const;
 Sequences* multiple_alignment(Format_error &error , std::ostream &os , const Vector_distance &ivector_dist ,   bool begin_free = false , bool end_free = false , int indel_cost = ADAPTATIVE ,double indel_factor = INDEL_FACTOR_N , int algorithm = AGGLOMERATIVE , const char *path = 0) const;
 Sequences* hierarchical_segmentation(Format_error &error , std::ostream &os , int iidentifier , int max_nb_segment , int *model_type) const;

 bool segment_profile_write(Format_error &error , std::ostream &os , int iidentifier ,int nb_segment , int *model_type , int output = SEGMENT , char format = 'a' , int segmentation = FORWARD_DYNAMIC_PROGRAMMING , int nb_segmentation = NB_SEGMENTATION) const;
 bool segment_profile_write(Format_error &error , const char *path , int iidentifier , int nb_segment , int *model_type , int output = SEGMENT , char format = 'a' , int segmentation = FORWARD_DYNAMIC_PROGRAMMING , int nb_segmentation = NB_SEGMENTATION) const;
 bool segment_profile_plot_write(Format_error &error , const char *prefix , int iidentifier , int nb_segment , int *model_type ,    int output = SEGMENT , const char *title = 0) const;

 Histogram* get_hlength() const { return hlength; }
 int get_index_parameter_type() const { return index_parameter_type; }
 Histogram* get_hindex_parameter() const { return hindex_parameter; }
 Histogram* get_index_interval() const { return index_interval; }
 int get_index_parameter(int iseq , int index) const  { return index_parameter[iseq][index]; }
 Histogram* get_marginal(int variable) const { return marginal[variable]; }
 int get_int_sequence(int iseq , int variable , int index) const    { return int_sequence[iseq][variable][index]; }
 double get_real_sequence(int iseq , int variable , int index) const    { return real_sequence[iseq][variable][index]; }

 */
}

#define WRAP SequenceCharacteristicsWrap
class WRAP
{
public:
  static boost::python::list
  get_initial_run(Sequence_characteristics& input)
  {

    boost::python::list list;

    for (int i = 0; i < input.get_nb_value(); i++)
      {
        //list.append(  boost::python::extract<Histogram>(input.get_initial_run(i)));
        list.append(i);

      }

    return list;
  }

  static Histogram*
  get_initial_run_from_index(Sequence_characteristics& input, int index)
  {
    return input.get_initial_run(index);
  }
};

void
class_sequence_characteristics()
{

class_<Sequence_characteristics> ("_SequenceCharacteristics", "SequenceCharacteristics")
.def(init<optional<int> >())
.def(init<Sequence_characteristics, optional<bool> >())
.def(init<Sequence_characteristics, optional<char> >())

.add_property("get_nb_value", &Sequence_characteristics::get_nb_value)

.def("get_initial_run", &WRAP::get_initial_run, "returns initial run")

DEF_RETURN_VALUE("get_initial_run_from_index", WRAP::get_initial_run_from_index,args("index"), "returns initial run")

DEF_RETURN_VALUE_NO_ARGS("get_index_value", &Sequence_characteristics::get_index_value, "get_index_value")
DEF_RETURN_VALUE_NO_ARGS("get_first_occurrence", &Sequence_characteristics::get_first_occurrence, "get first occurrence time")
DEF_RETURN_VALUE_NO_ARGS("get_recurrence_time", &Sequence_characteristics::get_recurrence_time, "get recurrence time")
DEF_RETURN_VALUE_NO_ARGS("get_sojourn_time", &Sequence_characteristics::get_sojourn_time,"returns sojourn time")
DEF_RETURN_VALUE_NO_ARGS("get_final_run", &Sequence_characteristics::get_final_run,"returns final run")
DEF_RETURN_VALUE_NO_ARGS("get_nb_run", &Sequence_characteristics::get_nb_run, "returns number of run")
DEF_RETURN_VALUE_NO_ARGS("get_nb_occurrence", &Sequence_characteristics::get_nb_occurrence, "returns number of ocurrences")
;

//DONE
/*      Histogram** get_initial_run() const { return initial_run; } */

}
#undef WRAP
