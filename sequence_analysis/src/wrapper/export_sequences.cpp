/*------------------------------------------------------------------------------
 *
 *
 *        VPlants.Sequence_analysis : VPlants Statistics module
 *
 *        Copyright 2006-2014 CIRAD/INRA/Inria Virtual Plants
 *
 *        File author(s): Yann Guedon <yann.guedon@cirad.fr>
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
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

#include "sequence_analysis/sequences.h"
#include "sequence_analysis/renewal.h"
// #include "sequence_analysis/tops.h"
#include "sequence_analysis/sequence_label.h"

#include <boost/python.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/list.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>

#include "boost_python_aliases.h"


using namespace boost::python;
using namespace boost;
using namespace stat_tool;
using namespace sequence_analysis;



class SequencesWrap
{

public:

  static boost::shared_ptr<Sequences>
  sequences_from_file(char *filename, bool old_format)
  {
    StatError error;
    Sequences *sequences = NULL;
    sequences = sequences_ascii_read(error, filename, old_format);
    return boost::shared_ptr<Sequences>(sequences);
  }

/*  static boost::shared_ptr<Sequences>
  sequences_from_file_old(char *filename, bool old_format)
  {
    StatError error;
    Sequences *sequences = NULL;
    sequences = sequences_ascii_read(error, filename, old_format);
    return boost::shared_ptr<Sequences>(sequences);
  }
*/

/*Sequences::Sequences(int inb_sequence , int *iidentifier , int *ilength ,
                     int **ivertex_identifier , int iindex_parameter_type ,
                     int **iindex_parameter , int inb_variable , int *itype ,
                     int ***iint_sequence , double ***ireal_sequence)

*/

//types is for sequences where all vectors have the same length
  static Sequences*
  build_from_lists_new(
        boost::python::list& input_sequences,
        boost::python::list& input_identifiers,
        boost::python::list& input_vertex_identifiers,
        boost::python::list& input_index_parameters,
        boost::python::list& input_types,
        int input_index_parameter_type
        )
  {
    int nb_sequences = boost::python::len(input_sequences);
    int nb_identifiers = boost::python::len(input_identifiers);
    int nb_types = boost::python::len(input_types);
    int nb_variables = boost::python::len(input_types);
    int nb_vectors = 0;

    Sequences *ret = NULL;

    int *lengths = NULL;
    int *identifiers = NULL;
    int *types = NULL;

    double ***real_sequences = NULL;
    int ***int_sequences = NULL;
    int **index_parameters = NULL;
    int **vertex_identifiers = NULL;
    int index_parameter_type = input_index_parameter_type;

    int nb_variable_float = 0;
    int nb_variable_int = 0;
    int index_int=0, index_float=0;
    boost::python::list sequences = extract<boost::python::list> (input_sequences);


    lengths = new int[nb_sequences];
    identifiers = new int[nb_identifiers];
    types = new int[nb_types];
    index_parameters = new int*[nb_sequences];



    //cout << "allocate memory and set values of index_parameters"<<endl;
    // allocate memory for index_parameters
    for (int iseq=0; iseq<nb_sequences; iseq++)
    {
        boost::python::list sequence = extract<boost::python::list> (input_sequences[iseq]);
        boost::python::list index_parameter = extract<boost::python::list> (input_index_parameters[iseq]);
        if (index_parameter_type == POSITION)
        {
            index_parameters[iseq] = new int[len(sequence)+1];
            for (int ilength=0; ilength<len(sequence)+1 ; ilength++)
            {
                index_parameters[iseq][ilength] = extract<int>(index_parameter[ilength]);
            }
        }
        else
        {
            index_parameters[iseq] = new int[len(sequence)];
            for (int ilength=0; ilength<len(sequence) ; ilength++)
            {
                index_parameters[iseq][ilength] = extract<int>(index_parameter[ilength]);
            }
        }
    }


    //cout << " nb identifiers="<< nb_identifiers<<endl;
    if (nb_identifiers > 0)
    {
      identifiers = new int[nb_identifiers];
      for (int ii = 0; ii < nb_identifiers; ii++)
        {
          identifiers[ii] = boost::python::extract<int>(input_identifiers[ii]);
        }
    }
    //cout << "identifier ok"<<endl;
    //types is for sequences where all vectors have the same lengths only and are homogeneous (same signature)
    // e.g., [ [ [1,1.],[2,2.] ] [ [3,3.],[4,4.] ] ] Type = [int, float]
    // but not [ [ [1.,1.],[2,2.] ] [ [3,3],[4,4.] ] ] Type = [int, float] sometimes is [float, int]
    // and not [ [ [1,1.],[2,2.] ] [ [3,3.,5],[4,4.,5] ] ] Type = [int, float] sometimes is [int, float, int]
    //cout << "nb_types="<<nb_types<<endl;
    if (nb_types!=0)
    {
      for (int ii = 0; ii < nb_types; ii++)
      {
        int dummy = boost::python::extract<int>(input_types[ii]);
        if (dummy==0){
            types[ii] = INT_VALUE;
            nb_variable_int++;
        }
        else {
            types[ii] = REAL_VALUE;
            nb_variable_float++;
        }
      }
    }
    //cout << "nb var int and float="<<nb_variable_int << " "<<nb_variable_float <<endl;
    // todo check that nb_sequences == nb_identifiers


    // get the list of sequences, assuming that input format is [ [ [1,1], [2,4] ], [ [1,4] ,[4,7] , [8,9]] ]
    // the most inner level being vectors of same length (2 variables here), the middle level being
    // the sequences of vectors, which may have different lengths (first sequence of length 2 and second of
    // length 3).

    // allocate memory
    if (nb_variable_int>0)
    {
        int_sequences = new int**[nb_sequences];
        for (int i=0; i<nb_sequences; i++)
        {
            int_sequences[i] = new int*[nb_variable_int];
        }
    }
    if (nb_variable_float>0)
    {
        real_sequences = new double**[nb_sequences];
        for (int i=0; i<nb_sequences; i++)
        {
            real_sequences[i] = new double*[nb_variable_float];
        }
    }


    vertex_identifiers = new int*[nb_sequences];

    //cout << "nb_sequences"<<len(input_sequences)<<endl;
    for (int iseq = 0; iseq < nb_sequences; iseq++)
    {
        boost::python::list sequence = extract<boost::python::list> (input_sequences[iseq]);
        boost::python::list vertex = extract<boost::python::list> (input_vertex_identifiers[iseq]);

        nb_vectors = len(sequence);
        //cout << "seq" << iseq << " and its length is " << len(sequence)<< endl;
        for (int kk=0; kk<nb_variable_int; kk++)
            int_sequences[iseq][kk] = new int[nb_vectors];
        for (int kk=0; kk<nb_variable_float; kk++)
            real_sequences[iseq][kk] = new double[nb_vectors];
        // fill the sequences lengths
        lengths[iseq] = len(sequence);
        //look at the first vector to get the length and therefore number of variables
        boost::python::list vector = extract<boost::python::list> (sequence[0]);
        nb_variables = len(vector);
        vertex_identifiers[iseq] = new int[nb_vectors];
        for (int ivec = 0; ivec < nb_vectors; ivec++)
        {

            boost::python::list vector = extract<boost::python::list> (sequence[ivec]);
            //cout << "length vector" << ivec << "="<< len(vector)<< " and expected is "<<  nb_variables <<endl;
            vertex_identifiers[iseq][ivec] = extract<int>(vertex[ivec]);


            index_int = 0;
            index_float= 0;

            for (int ivar=0; ivar < nb_variables; ivar++)
            {
                //cout << "------" << ivar<<endl;
                //TODO switch between int or real sequences
                if (types[ivar]==INT_VALUE)
                {
                    //cout <<"before affect int"<< extract<int>(vector[ivar])<<endl;
                    int_sequences[iseq][index_int][ivec] =  extract<int> (vector[ivar]);
                    //cout <<"aftere affect int"<<endl;
                    index_int++;
                }
                else
                {
                    //cout <<"before affect float"<<extract<double>(vector[ivar])<<endl;
                    real_sequences[iseq][index_float][ivec] =  extract<double> (vector[ivar]);
                    //cout <<"after affect float"<<endl;
                    index_float++;

                }
            }
        }
    }

/*
    for (int iseq = 0; iseq < nb_sequences; iseq++){
        boost::python::list sequence = extract<boost::python::list> (input_sequences[iseq]);
        nb_vectors = len(sequence);
        boost::python::list vector = extract<boost::python::list> (sequence[0]);
        for (int ivec = 0; ivec < nb_vectors; ivec++)
            cout << iseq <<","<<ivec<<"="<<vertex_identifiers[iseq][ivec] <<endl;
    }
*/

    //cout << "calling Sequences "<<endl;
    ret = new Sequences(nb_sequences, identifiers, lengths,
        vertex_identifiers, index_parameter_type, index_parameters,
        nb_variables, types, int_sequences, real_sequences);
    //cout << "calling Sequences done "<<endl;


    //cout << "freeing memory rela sequences "<<endl;
    if (real_sequences)
    {


      for (int i = 0; i < nb_sequences; i++)
      {
        if (nb_variable_float>0)
        {
          for (int j = 0; j < nb_variable_float; j++)
          {
            delete[] real_sequences[i][j];
          }
          delete[] real_sequences[i];
        }
      }
      delete[] real_sequences;
    }
    //cout << "freeing memory rela sequences done"<<endl;
    //cout << "freeing memory int sequences "<<endl;
    if (int_sequences)
    {
      for (int i = 0; i < nb_sequences; i++)
      {
        if (nb_variable_int>0)
        {
          for (int j = 0; j < nb_variable_int; j++)
          {
            delete[] int_sequences[i][j];
          }
          delete[] int_sequences[i];
        }
      }
      delete[] int_sequences;
    }

    if (vertex_identifiers)
    {
        for (int iseq = 0; iseq < nb_sequences; iseq++)
            delete [] vertex_identifiers[iseq];
      delete [] vertex_identifiers;
    }

    if (index_parameters)
    {
        for (int iseq = 0; iseq < nb_sequences; iseq++)
            delete [] index_parameters[iseq];
      delete [] index_parameters;
    }



    if (identifiers){
        delete[] identifiers;
    }
    if (lengths){
        delete[] lengths;
    }
    if (types){
        delete[] types;
    }

    return ret;
  }


  // very complicated implementation of Sequence, but seems to work for now
  static Sequences*
  build_from_lists(boost::python::list& input_list,
      boost::python::list& input_identifiers, int input_index_parameter_type)
  {

    // case 1 list of n lists of floats (only 1 variable and different sizes possible):
    //          seq = Sequences([ [1,1,1], [2,2,2,2,2]])
    //
    // case 2 list of n lists (different sizes) of variables (same size)
    //          seq = Sequences([ [ [1,1,1], [11,11,11] ],
    //                            [ [2,2,2,2,2], [22,22,22,22,22] ]
    //                          ])

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

    if (index_parameter_type != IMPLICIT_TYPE)
        nb_variables--;


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

  static DiscreteDistributionData*
  extract_value(const Sequences &input, int variable)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, extract, DiscreteDistributionData, variable);
  }

  static boost::python::list
  get_identifiers(const Sequences &seq)
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
  ascii_data_write(const Sequences &d, bool exhaustive)
  {
    std::stringstream os;
    std::string res;
    d.ascii_data_write(os, exhaustive, exhaustive);
    res = os.str();
    return res;
  }

  static void
  file_ascii_data_write(const Sequences &d, const char *path, bool exhaustive)
  {
    bool result = true;
    StatError error;

    result = d.ascii_data_write(error, path, exhaustive);
    if (!result)
      sequence_analysis::wrap_util::throw_error(error);

  }

  static void
  file_ascii_write(const Sequences &d, const char *path, bool exhaustive)
  {
    bool result = true;
    StatError error;

    result = d.ascii_write(error, path, exhaustive);
    if (!result)
      sequence_analysis::wrap_util::throw_error(error);

  }


  // used to return a vector
  //   seq = [[[1,2],[2,3]]]
  //  seq[0,1] -> [2,3] vector 1 of sequence 0
  static boost::python::list
  get_item_tuple(const Sequences *seq, boost::python::tuple indexes)
  {
    int index_vec = extract<int> (indexes[1]);
    int index_seq = extract<int> (indexes[0]);

    // Test index
    if (index_seq < 0 || index_seq >= seq->get_nb_sequence())
      {
        PyErr_SetString(PyExc_IndexError, "sequence index out of bound");
        boost::python::throw_error_already_set();
      }
    if (index_vec < 0 || index_vec >= seq->get_length(index_seq))
      {
        PyErr_SetString(PyExc_IndexError, "variable index out of bound");
        boost::python::throw_error_already_set();
      }

    boost::python::list l;

    int nb_length = seq->get_length(index_seq);
    int nb_variables = seq->get_nb_variable();
    for (int index = 0; index < nb_variables; index++)
      {
        if ((seq->get_type(index) == INT_VALUE)
            || (seq->get_type(index) == STATE))
          l.append(seq->get_int_sequence(index_seq, index, index_vec));
        else
          l.append(seq->get_real_sequence(index_seq, index, index_vec));
      }

    return l;
  }

  static boost::python::list
  get_item_int(const Sequences *seq, int index_seq)
  {
    // Test index
    if (index_seq < 0 || index_seq >= seq->get_nb_sequence())
      {
        PyErr_SetString(PyExc_IndexError, "sequence index out of bound");
        boost::python::throw_error_already_set();
      }

    boost::python::list l;
    int nb_length = seq->get_length(index_seq);
    int nb_var = seq->get_nb_variable();

    for (int indexvec = 0; indexvec<nb_length; indexvec++)
    {
      boost::python::list subl;


      for (int indexvar = 0; indexvar < nb_var; indexvar++)
      {

      //boost::python::list subl;
      //for (int index = 0; index < nb_length; index++)
      //{
        if ((seq->get_type(indexvar) == INT_VALUE)
            || (seq->get_type(indexvar) == STATE))
          subl.append(seq->get_int_sequence(index_seq, indexvar, indexvec));
        else
          subl.append(seq->get_real_sequence(index_seq, indexvar, indexvec));
      }
      l.append(subl);
    }
    return l;
  }




  static Sequences*
  value_select(const Sequences &seq, int variable, const object& min,
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
  select_variable(const Sequences &seq, const boost::python::list& variables,
      bool keep)
  {

    CREATE_ARRAY(variables, int, vars);

    SIMPLE_METHOD_TEMPLATE_1(seq, select_variable, Sequences, vars_size,
        vars.get(), keep);

  }

  static Sequences*
  select_individual(const Sequences &seq,
      const boost::python::list& identifiers, bool keep)
  {
    CREATE_ARRAY(identifiers, int, ids);

    SIMPLE_METHOD_TEMPLATE_1(seq, select_individual, Sequences, ids_size,
        ids.get(), keep);
  }

  // Shift
  static Sequences*
  shift(const Sequences &seq, int var, double param)
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
  merge(const Sequences &input_seq, const boost::python::list& seqs)
  {

    CREATE_ARRAY(seqs, const Sequences*, sequens)

    SIMPLE_METHOD_TEMPLATE_1(input_seq, merge, Sequences, sequens_size, sequens.get());
  }

  static Sequences*
  merge_variable(const Sequences &input_seq, const boost::python::list& seqs,
      int ref_sample)
  {

    CREATE_ARRAY(seqs, const Sequences*, data);

    SIMPLE_METHOD_TEMPLATE_1(input_seq, merge_variable, Sequences,data_size,
        data.get(), ref_sample);
  }

  // Cluster
  static Sequences*
  cluster_step(const Sequences &seq, int variable, int step)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, cluster, Sequences, variable, step);
  }

  static Sequences*
  cluster_limit(const Sequences &seq, int variable, boost::python::list& limit)
  {

    StatError error;

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
  transcode(const Sequences &seq, int variable, boost::python::list& input_symbols)
  {

    int nb_symbol =  boost::python::len(input_symbols);
    //seq.get_nb_sequence();
    sequence_analysis::wrap_util::auto_ptr_array<int> l(new int[nb_symbol]);

    int expected_nb_symbol = (int) (seq.get_max_value(variable - 1)
        - seq.get_min_value(variable - 1)) + 1;

    if (nb_symbol != expected_nb_symbol)
    {
        cout << "expected_nb_symbol="<<expected_nb_symbol <<endl;
        sequence_analysis::wrap_util::throw_error("Bad number of Symbol");
    }


    for (int i = 0; i < nb_symbol; i++)
      l[i] = extract<int> (input_symbols[i]);

    SIMPLE_METHOD_TEMPLATE_1(seq, transcode, Sequences, variable, l.get());

    //TODO in aml, if seq.transcode is correct, there is an additional call
    // to seq.markovian_sequences()
  }

  static Sequences*
  reverse(const Sequences &seq)
  {
    SIMPLE_METHOD_TEMPLATE_0(seq, reverse, Sequences);
  }

  static Sequences*
  length_select(const Sequences &input,  int min_length, int max_length, bool keep)
  {
     StatError error;
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
  remove_run(const Sequences &seq, int variable, int ivalue, char position,
      int max_run_length)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, remove_run, Sequences, variable, ivalue,
        position, max_run_length);
  }

  static int
  get_type(const Sequences &seq, int index)
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
  get_length(const Sequences &seq, int index)
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
  get_min_value(const Sequences &seq, int variable)
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
  get_index_parameter(const Sequences &seq, int iseq, int index)
  {
    if (iseq < 0 || iseq >= seq.get_nb_sequence())
      {
        PyErr_SetString(PyExc_IndexError,
            "id of sequence must be positive and less than number of sequences");
        boost::python::throw_error_already_set();
      }

    if (index < 0 || index >= seq.get_length(iseq))
      {
        PyErr_SetString(PyExc_IndexError,
            "index must be positive and less than sequence length");
        boost::python::throw_error_already_set();
      }

    return seq.get_index_parameter(iseq, index);
  }

  static double
  get_max_value(const Sequences &seq, int variable)
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
  scaling(const Sequences &seq, int variable, int scaling_coeff)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, scaling, Sequences, variable, scaling_coeff);
  }

  static Sequences*
  round(const Sequences &seq, int variable, int mode)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, round, Sequences, variable, mode);
  }

  //cumulate
  static Sequences*
  cumulate(const Sequences &seq, int variable)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, cumulate, Sequences, variable);
  }

  static Sequences*
  difference(const Sequences &seq, int variable, bool first)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, difference, Sequences, variable, first);
  }

  static Sequences*
  index_parameter_select(const Sequences &input, int min_index_parameter,
      int max_index_parameter, bool keep)
  {
    StatError error;
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
  remove_index_parameter(const Sequences &seq)
  {
    SIMPLE_METHOD_TEMPLATE_0(seq, remove_index_parameter, Sequences);
  }

  static Sequences*
  index_parameter_extract(const Sequences &seq, int min_index_parameter,
      int max_index_parameter)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, index_parameter_extract, Sequences,
        min_index_parameter, max_index_parameter);
  }

  static Sequences*
  segmentation_extract(const Sequences &seq, int variable,
      boost::python::list& input_values, bool keep)
  {
    CREATE_ARRAY(input_values, int, values);
    SIMPLE_METHOD_TEMPLATE_1(seq, segmentation_extract, Sequences,
        variable, values_size, values.get(), keep);
  }

  static Sequences*
  moving_average(const Sequences &seq, boost::python::list& input_values,
      int variable, bool begin_end, int output)
  {

    int nb_value = len(input_values);
    double sum = 0;
    int i = 0;

    sequence_analysis::wrap_util::auto_ptr_array<double> values(new double[nb_value * 2 + 1]);

    nb_value--;
    for (i = 0; i < nb_value; i++)
    {
      values[i] = extract<double> (input_values[i]);
      values[2 * nb_value - i] = values[i];
      sum += 2 * values[i];
    }
    // i = n
    values[i] = extract<double> (input_values[i]);
    sum += values[i];

    // check that sum is 1 ?

    //normalization
    for (i = 0; i < 2 * nb_value + 1; i++)
    {
      values[i] = values[i] / sum;
    }


    SIMPLE_METHOD_TEMPLATE_1(seq, moving_average, Sequences,
            nb_value, values.get(), variable, begin_end,  output);

  }

  static Sequences*
  moving_average_from_distribution(const Sequences &seq,
      const Distribution& dist, int variable, bool begin_end, int output)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, moving_average, Sequences, dist, variable,
        begin_end, output);
  }

  static Sequences*
  pointwise_average(const Sequences &seq, bool circular, bool standard_deviation,
      int output, const char *path, char format)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, pointwise_average, Sequences,
        circular, standard_deviation, output, path, format);
  }

  static Sequences*
  recurrence_time_sequences(const Sequences &seq, int variable, int value)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, recurrence_time_sequences, Sequences,
        variable, value);
  }

  static Sequences*
  cross(const Sequences &seq)
  {
    SIMPLE_METHOD_TEMPLATE_0(seq, cross, Sequences);
  }

  static Sequences*
  transform_position(const Sequences &seq, int step)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, transform_position, Sequences, step);
  }

  static Sequences*
  sojourn_time_sequences(const Sequences &seq, int variable)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, sojourn_time_sequences, Sequences, variable);
  }

  static Correlation*
  partial_autocorrelation_computation(const Sequences &seq, int variable,
      int itype, int max_lag)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, partial_autocorrelation_computation,
        Correlation, variable, itype, max_lag);
  }

  static Vectors*
  extract_vectors(const Sequences &seq, int feature_type, int variable, int value)
  {
    SIMPLE_METHOD_TEMPLATE_1(seq, extract_vectors, Vectors, feature_type,
        variable, value);
  }

  static DiscreteDistributionData*
  extract_length(const Sequences &input)
  {
    DiscreteDistributionData *res;
    res = new DiscreteDistributionData(*(input.get_length_distribution()));
    return res;
  }

  static MarkovianSequences*
  markovian_sequences(const Sequences &input)
  {
    SIMPLE_METHOD_TEMPLATE_0(input, markovian_sequences, MarkovianSequences);
  }

  static Correlation*
  correlation_computation(const Sequences &input, int variable1, int variable2,
  int itype, int max_lag, int normalization, bool individual_mean)
  {
    SIMPLE_METHOD_TEMPLATE_1 (input, correlation_computation,
      Correlation, variable1, variable2, itype, max_lag, normalization, individual_mean)
  }

  static Sequences*
  multiple_alignment(const Sequences &input, const VectorDistance &ivector_dist,
      bool begin_free, bool end_free, int indel_cost, double indel_factor, int algorithm, const char *path)
  {
    HEADER_OS(Sequences);
    ret = input.multiple_alignment(error, os, ivector_dist, begin_free,\
        end_free, indel_cost, indel_factor, algorithm, path);
    FOOTER_OS;

  }

  static DistanceMatrix*
  alignment_vector_distance(const Sequences &input,
      const VectorDistance &ivector_dist,  int ref_identifier,
      int test_identifier, bool begin_free, bool end_free, int indel_cost,
      double indel_factor, bool transposition_flag,
      double transposition_factor, const char *result_path, char result_format,
      const char *alignment_path, char alignment_format)
  {
    HEADER_OS(DistanceMatrix);

    ret = input.alignment(error, &os, ivector_dist, ref_identifier, test_identifier, begin_free,
        end_free, indel_cost, indel_factor, transposition_flag,
        transposition_factor, result_path, result_format, alignment_path,
        alignment_format );

   FOOTER_OS;
  }

  static DistanceMatrix*
  alignment(const Sequences &input,int ref_identifier,
  int test_identifier, bool begin_free, bool end_free,
  const char *result_path, char result_format,
  const char *alignment_path, char alignment_format)
  {
    HEADER_OS(DistanceMatrix);

    ret = input.alignment(error, &os, ref_identifier, test_identifier, begin_free,
      end_free, result_path, result_format, alignment_path,
      alignment_format );
    FOOTER_OS;
  }

  static RenewalData*
  extract_renewal_data(const Sequences &input,
    int variable , int begin_index , int end_index)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, extract_renewal_data, RenewalData,
      variable, begin_index, end_index);
  }

  static TimeEvents*
  extract_time_events(const Sequences &input,
    int variable , int begin_date , int end_date ,
    int previous_date, int next_date)

  {
    SIMPLE_METHOD_TEMPLATE_1(input, extract_time_events, TimeEvents,
      variable, begin_date, end_date, previous_date, next_date);

  }

  static Vectors*
  build_vectors(const Sequences &input,  bool index_variable)
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

  // fonction supprimee

/*  static Sequences*
  segmentation_vector_distance(const Sequences &input,  int iidentifier,
    int nb_segment ,  const VectorDistance &ivector_dist, int output)
  {
    HEADER_OS(Sequences);
    ret = input.segmentation(error,iidentifier,  nb_segment,
      ivector_dist, os, output);
    FOOTER_OS;
  } */

  static void plot_write(const Sequences &input, const std::string& prefix,
    const std::string& title)
  {
      StatError error;
      input.plot_write(error, prefix.c_str(), title.c_str());
  }
  static void plot_data_write(const Sequences &input, const std::string& prefix,
  const std::string& title)
  {
      StatError error;
      input.plot_data_write(error, prefix.c_str(), title.c_str());
  }

  /*static bool
  segment_profile_write(const Sequences &input, const char *prefix, int iidentifier, int nb_segment,
    boost::python::list model_type, int output, char title)
  {
    StatError error;
    bool ret;
    CREATE_ARRAY(model_type, int, models);

    ret = input.segment_profile_write(error, prefix, iidentifier, nb_segment,
         models.get(), output, title);
    FOOTER;
  }
  */
  static MultiPlotSet*
  get_plotable(const Sequences &p)
  {
    MultiPlotSet* ret = p.get_plotable();
    if (!ret)
      ERROR;
    return ret;
  }

  static MultiPlotSet*
  get_plotable_data(const Sequences &p)
  {
    StatError error;
    Sequences seq;
    seq = (Sequences )p;
    MultiPlotSet* ret = seq.get_plotable_data(error);
    if (!ret)
      ERROR;
    return ret;
  }

  static MultiPlotSet*
  segment_profile_plotable_write(const Sequences &input, int iidentifier ,
      int nb_segment, boost::python::list model_type , int output)
  {
    StatError error;
    CREATE_ARRAY(model_type, int, models);
    MultiPlotSet* ret = input.segment_profile_plotable_write(error, iidentifier, nb_segment, models.get(), output);

    if (!ret)
      ERROR;
    return ret;
  }

  static bool
  segment_profile_write(const Sequences &input,int iidentifier,
                               int nb_segment , boost::python::list& model_type , int output ,
                               char format, int segmentation ,
                               int nb_segmentation)
  {
    std::stringstream os;
    StatError error;
    bool ret;
    CREATE_ARRAY(model_type, int, models);
    ret = input.segment_profile_write(error, os, iidentifier,  nb_segment,
        models.get(), output, format, segmentation, nb_segmentation);
    if (!ret)
      sequence_analysis::wrap_util::throw_error(error);
    cout << os.str() << endl;
    return ret;
  }

  static bool
  select_step(Sequences &input, int variable, double step)
  {
    StatError error;
    bool ret;



    ret = input.select_step(error, variable, step);
    if (!ret)
      sequence_analysis::wrap_util::throw_error(error);
    return ret;
  }

  static Histogram*
  get_marginal_histogram(Sequences &input, int variable)
  {
    Histogram *ret;
    ret = input.get_marginal_histogram(variable);
    return ret;
  }


};

// Boost declaration

void
class_sequences()
{


  class_<Sequences, bases<StatInterface> > ("_Sequences", "Sequences")
    .def("__init__", make_constructor(SequencesWrap::sequences_from_file))
   // .def("__init__", make_constructor(SequencesWrap::sequences_from_file_old))
    //.def("__init__", make_constructor(SequencesWrap::build_from_lists))
    .def("__init__", make_constructor(SequencesWrap::build_from_lists_new))
    .def(init <const RenewalData&>())

    // Python Operators

   .def(self_ns::str(self)) //__str__
   .def("__len__", &Sequences::get_nb_sequence,"Returns number of sequences")
   .def("__getitem__", SequencesWrap::get_item_tuple)
   .def("__getitem__", SequencesWrap::get_item_int)

   //property
   .add_property("nb_sequence", &Sequences::get_nb_sequence, "Return the number of sequences")
   .add_property("nb_variable", &Sequences::get_nb_variable, "Return the number of variables")
   .add_property("max_length", &Sequences::get_max_length,"Return max length")
   .add_property("cumul_length", &Sequences::get_cumul_length,"Return cumul length")
   .add_property("index_parameter_type", &Sequences::get_index_parameter_type,"return index parameter type")

   .def("get_min_value", &Sequences::get_min_value,args("index_var"), "return min value of variables")
   .def("get_max_value", &Sequences::get_max_value,args("index_var"), "return max value of variables")
   // index arguments, wrapping required
   .def("get_length", &SequencesWrap::get_length,args("index_seq"), "return length of a sequence")
   .def("get_type", &SequencesWrap::get_type,args("index_var"),"return  type")
   .def("get_identifiers", &SequencesWrap::get_identifiers,"returns list of identifiers")
   .def("ascii_data_write", SequencesWrap::ascii_data_write,"Return a string with the object representation")
   .def("file_ascii_write", SequencesWrap::file_ascii_write,"Save vector summary into a file")
   .def("file_ascii_data_write", SequencesWrap::file_ascii_data_write,"Save vector data into a file")

   .def("get_index_parameter", SequencesWrap::get_index_parameter,args("iseq", "index"), "return index")

   DEF_RETURN_VALUE("alignment_vector_distance", SequencesWrap::alignment_vector_distance, args(""), "todo")
   DEF_RETURN_VALUE("alignment", SequencesWrap::alignment, args(""), "todo")
   DEF_RETURN_VALUE("multiple_alignment", SequencesWrap::multiple_alignment, args(""), "todo")
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
   DEF_RETURN_VALUE("pointwise_average", SequencesWrap::pointwise_average, args("circular", "standard_deviation", "output", "path", "format"), "Pointwise average")
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
//   DEF_RETURN_VALUE("segmentation_vector_distance", SequencesWrap::segmentation_vector_distance, args("todo"), "todo")
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
   DEF_RETURN_VALUE_NO_ARGS("get_plotable", SequencesWrap::get_plotable, "Return a plotable")
   DEF_RETURN_VALUE_NO_ARGS("get_plotable_data", SequencesWrap::get_plotable_data, "Return a plotable_data")

   .def("plot_write", SequencesWrap::plot_write, args("prefix", "title"), "Write GNUPLOT files")
   .def("plot_data_write", SequencesWrap::plot_data_write, args("prefix", "title"), "Write GNUPLOT files")
   DEF_RETURN_VALUE("segment_profile_plotable_write", SequencesWrap::segment_profile_plotable_write, args("identifier", "nb_segment", "model_type", "output"), "Write segment_profile")

    .def("segment_profile_write", SequencesWrap::segment_profile_write, args("sequences", "iidentifier","nb_segment", "model_type" , "output" ,"format","segmentation","nb_segmentation"), "segment profile write for Display")
    .def("select_step", SequencesWrap::select_step, args("variable", "step"), "select_step on sequences")
   DEF_RETURN_VALUE("get_marginal_histogram", SequencesWrap::get_marginal_histogram, args("variable"), "get_marginal_histogram wrapper")



;

/*
 Sequences(int inb_sequence , int inb_variable);
 Sequences(int inb_sequence , int *iidentifier , int *ilength , int iindex_parameter_type ,               int inb_variable , int *itype , bool init_flag = false)       { init(inb_sequence , iidentifier , ilength , iindex_parameter_type , inb_variable ,               itype , init_flag); }
 Sequences(int inb_sequence , int *iidentifier , int *ilength , int inb_variable , int *itype , bool init_flag = false){ init(inb_sequence , iidentifier , ilength , IMPLICIT_TYPE , inb_variable ,            itype , init_flag); }
 Sequences(int inb_sequence , int *iidentifier , int *ilength , int inb_variable ,  bool init_flag = false)     { init(inb_sequence , iidentifier , ilength , inb_variable , init_flag); }
 Sequences(const FrequencyDistribution &ilength_distribution , int inb_variable , bool init_flag = false);
 Sequences(const Sequences &seq , int inb_sequence , int *index);
 Sequences(const Sequences &seq , bool *segment_mean);
 Sequences(const Sequences &seq , char transform = 'c' , int param = DEFAULT);
 Sequences(const FrequencyDistribution &ilength_distribution , int inb_variable , bool init_flag = false);

 Tops* tops(StatError &error) const;
 bool check(StatError &error , const char *pattern_label);

 std::ostream& line_write(std::ostream &os) const;
 bool spreadsheet_write(StatError &error , const char *path) const;
 bool plot_write(StatError &error , const char *prefix ,   const char *title = 0) const;

 int min_index_parameter_computation() const;
 int max_index_parameter_computation(bool last_position = false) const;

 void marginal_histogram_computation(int variable);
 double mean_computation(int variable) const;
 double variance_computation(int variable , double mean) const;
 double mean_absolute_deviation_computation(int variable , double mean) const;
 double mean_absolute_difference_computation(int variable) const;
 double skewness_computation(int variable , double mean , double variance) const;
 double kurtosis_computation(int variable , double mean , double variance) const;

 FrequencyDistribution* value_index_interval_computation(StatError &error , int variable , int value) const;
 Sequences* hierarchical_segmentation(StatError &error , std::ostream &os , int iidentifier , int max_nb_segment , int *model_type) const;


    bool segment_profile_write(StatError &error , const char *path , int iidentifier ,
                               int nb_segment , int *model_type , int output = SEGMENT ,
                               char format = 'a' , int segmentation = FORWARD_DYNAMIC_PROGRAMMING ,
                               int nb_segmentation = NB_SEGMENTATION) const;
    bool segment_profile_plot_write(StatError &error , const char *prefix ,
                                    int iidentifier , int nb_segment , int *model_type ,
                                    int output = SEGMENT , const char *title = NULL) const;




 FrequencyDistribution* get_length_distribution() const { return length_distribution; }
 FrequencyDistribution* get_index_parameter_distribution() const { return index_parameter_distribution; }
 FrequencyDistribution* get_index_interval() const { return index_interval; }
 int get_index_parameter(int iseq , int index) const  { return index_parameter[iseq][index]; }
 FrequencyDistribution* get_marginal(int variable) const { return marginal[variable]; }
 int get_int_sequence(int iseq , int variable , int index) const    { return int_sequence[iseq][variable][index]; }
 double get_real_sequence(int iseq , int variable , int index) const    { return real_sequence[iseq][variable][index]; }

 */
}

#define WRAP SequenceCharacteristicsWrap
class WRAP
{
public:
  static boost::python::list
  get_initial_run(SequenceCharacteristics &input)
  {

    boost::python::list list;

/*    for (int i = 0; i < input.get_nb_value(); i++)
      {
        list.append(  boost::python::extract<FrequencyDistribution>(input.get_initial_run(i)));
        extract<FrequencyDistribution>(input.get_initial_run(i));
        list.append(i);

      } */

    return list;
  }

/*  static FrequencyDistribution*
  get_initial_run_from_index(SequenceCharacteristics &input, int index)
  {
    return input.get_initial_run(index);
  } */
};

void
class_sequence_characteristics()
{

  class_<SequenceCharacteristics> ("_SequenceCharacteristics", "SequenceCharacteristics")
//  .def(init<optional<int> >())
//  .def(init<SequenceCharacteristics, optional<bool> >())
//  .def(init<SequenceCharacteristics, optional<char> >())

/*  .add_property("get_nb_value", &SequenceCharacteristics::get_nb_value)

  .def("get_initial_run", &WRAP::get_initial_run, "returns initial run")

  DEF_RETURN_VALUE("get_initial_run_from_index", WRAP::get_initial_run_from_index,args("index"), "returns initial run")

  DEF_RETURN_VALUE_NO_ARGS("get_index_value", &SequenceCharacteristics::get_index_value, "get_index_value")
  DEF_RETURN_VALUE_NO_ARGS("get_first_occurrence", &SequenceCharacteristics::get_first_occurrence, "get first occurrence time")
  DEF_RETURN_VALUE_NO_ARGS("get_recurrence_time", &SequenceCharacteristics::get_recurrence_time, "get recurrence time")
  DEF_RETURN_VALUE_NO_ARGS("get_sojourn_time", &SequenceCharacteristics::get_sojourn_time,"returns sojourn time")
  DEF_RETURN_VALUE_NO_ARGS("get_final_run", &SequenceCharacteristics::get_final_run,"returns final run")
  DEF_RETURN_VALUE_NO_ARGS("get_nb_run", &SequenceCharacteristics::get_nb_run, "returns number of run")
  DEF_RETURN_VALUE_NO_ARGS("get_nb_occurrence", &SequenceCharacteristics::get_nb_occurrence, "returns number of ocurrences") */
  ;

//DONE
/*      FrequencyDistribution** get_initial_run() const { return initial_run; } */

}
#undef WRAP
