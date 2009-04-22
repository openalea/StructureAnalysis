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
 *        $Id: export_tops.cpp 6169 2009-04-01 16:42:59Z cokelaer $
 *
 *-----------------------------------------------------------------------------*/

#include "wrapper_util.h"


#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"
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

using namespace boost::python;
using namespace boost;
using namespace stat_tool;



/*
boost::python::list extract_list(object x)
{
	extract<list> get_list((x));

	// is it a list ?
	bool is_list1 get_list.check();
	bool is_list2 PyObject.IsInstance(x.ptr(), (PyObject*)&PyList_Type);
	if (is_list1 != is_list2){
		throw std::runtime_error("is_list1 == is_list2 failure!")
	}
	return get_list()

}*/
////////////////////// Export tops ////////////////////////////////////////

class SequencesWrap
{

public:

  static boost::shared_ptr<Sequences> sequences_from_file(char* filename)
  {
    Format_error error;
    Sequences *sequences= NULL;
    bool old_format=false;

    sequences = sequences_ascii_read(error, filename, old_format);

/*    if(!top_parameters)
    {
	  stat_tool::wrap_util::throw_error(error);
    }
*/
    return boost::shared_ptr<Sequences>(sequences);
  }


  static Sequences* build_from_lists(boost::python::list& input_list)
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
	int index_parameter_type;
	int type;

	// length of each sequence will be store in this variable
	length = new int[nb_sequences];
	identifiers = new int[nb_sequences];

	//extract<boost::python::list> get_list(input_list[0]);
	boost::python::list seq = extract<boost::python::list>(input_list[0]);
	extract<boost::python::list> get_list(seq[0]);

	if (!get_list.check())
	{
		boost::python::list seq0 = extract<boost::python::list>(input_list[0]);
		object elt= seq0[0];
		extract<float> get_float(elt);
		extract<int> get_int(elt);

		if (get_int.check())
			is_int = true;
		else
			if (get_float.check())
				is_float = true;
	}
	else
	{

		boost::python::list sequence = extract<boost::python::list>(input_list[0]);
		boost::python::list variable = extract<boost::python::list>(sequence[0]);

		extract<float> get_float(variable[0]);
		extract<int> get_int(variable[0]);

		if (get_int.check())
			is_int = true;
		else
			if (get_float.check())
				is_float = true;
	}

	// allocate memory given the number of sequences and the type
	if ( is_float )
		real_sequence = new double**[nb_sequences];
	if ( is_int )
		int_sequence = new int**[nb_sequences];

	//cerr << is_float << endl;
//	cerr << is_int << endl;

	// for each sequence, are we considering case 1 or 2 ?
	for (int seqi=0; seqi<nb_sequences; seqi++)
	{
		// whatever case it is, the identifiers can be set here
		identifiers[seqi] = seqi;

		// try to get a single sequence
		boost::python::list sequence = extract<boost::python::list>(input_list[seqi]);

		// is it case 2 i.e. there is another nested list ?
		extract<boost::python::list> get_list(sequence[0]);

		// if not a list, we are in the case 1
		if ( !get_list.check() )
		{
			// the length of the sequence is simply:
			int N = len(sequence);
			length[seqi] = N;

			// and there is only 1 variable
			nb_variables = 1;				// used by the constructor

			// so the data structure is as follows
			if ( is_float )
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
			for (int j=0; j<N; j++)
			{
				if ( is_float )
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
			nb_variables = len(sequence);			// used by the constructor

			// what are the number of variables in this constructor :
			int N_var = len(sequence);

			// loop over the variable to extract their contents
			for (int vari=0; vari<N_var; vari++)
			{
				// get the variable of the current sequence
				boost::python::list variable =
					extract<boost::python::list>(sequence[vari]);

				// and update the data structure
				int N = len(variable);
				length[seqi] = N;

				//real_sequence[seqi][vari] = new double[N];
		//		cerr << "---"<<N<< endl;

				// before extracting the vector of the

				if ( is_float )
				{
			//		cerr << " is float validated"<<endl;
					if (vari==0) real_sequence[seqi] = new double*[N_var]; // only once
					real_sequence[seqi][vari] = new double[N];
				}
				if ( is_int )
				{
			//		cerr << " is int chjosen"<<endl;
					if (vari==0) int_sequence[seqi] = new int*[N_var]; // only once
					int_sequence[seqi][vari] = new int[N];
				}

				// get the vector corresponding to a sequence/variable pair
				for (int j=0; j<N; j++)
				{
					if ( is_float )
					{
						float v = extract<float> (variable[j]);
						real_sequence[seqi][vari][j] = v;
					}

					if ( is_int )
					{
						int v = extract<int> (variable[j]);
						int_sequence[seqi][vari][j] = v;

					}

				}
			}
		}
	}

//				{
//					ostringstream error_message;
//					error_message << "Inconsistent number of variables between sequences" <<endl;
//					PyErr_SetString(PyExc_ValueError, (error_message.str().c_str()));
//					throw_error_already_set();
//					cerr <<"ERROR" <<endl;
//				}
	//index_parameter_type == TIME
//	index_parameter_type == POSITION


//type[0] != NB_INTERNODE


  //index_parameter_type == POSITION



	// now, we can call this constructor that returns a Sequences of REAL
	// and free memeory

//	cerr << "before constructor"<<endl;
	std::flush(cerr);
	if ( is_float )
    {
    	sequences = new Sequences(nb_sequences, identifiers, length, nb_variables,
    								real_sequence);

    	if (real_sequence){
			for (int i=0; i<nb_sequences; i++)
			{
				for (int j=0; j<nb_variables; j++)
				{
					delete [] real_sequence[i][j];
				}
				delete [] real_sequence[i];
			}
			delete [] real_sequence;
		}
    }

	if ( is_int )
	{
		sequences = new Sequences(nb_sequences, identifiers, length,
				index_parameter_type, nb_variables, type, int_sequence);
		// freeing memory

	    if (int_sequence){
			for (int i=0; i<nb_sequences; i++)
			{
				for (int j=0; j<nb_variables; j++)
				{
					delete [] int_sequence[i][j];
				}
				delete [] int_sequence[i];
			}

			delete [] int_sequence;
		}
	}

	delete [] length;
	delete [] identifiers;

	return sequences;

  }


  static Distribution_data* extract_histogram(const Sequences& seq, int variable)
  {
	Format_error error;
	Distribution_data *ret = NULL;

	ret = seq.extract(error, variable);

	if (!ret)
		stat_tool::wrap_util::throw_error(error);

  }

  static boost::python::list get_identifiers(const Sequences& seq)
    {
      boost::python::list l;

      int  nb_seq = seq.get_nb_sequence();

      for(int s=0; s<nb_seq; s++)
      {
        l.append(seq.get_identifier(s));
      }

      return l;
    }

    static std::string ascii_data_write(const Sequences& d, bool exhaustive)
    {
      std::stringstream s;
      std::string res;


      d.ascii_data_write(s, exhaustive, true);
      res = s.str();
      return res;
    }

    static void file_ascii_data_write(const Sequences& d, const char* path, bool exhaustive)
    {
      bool result = true;
      Format_error error;

      result = d.ascii_data_write(error, path,exhaustive);
      if (!result)
         stat_tool::wrap_util::throw_error(error);

    }


    static void file_ascii_write(const Sequences& d, const char* path, bool exhaustive)
    {
      bool result = true;
      Format_error error;

      result = d.ascii_write(error, path,exhaustive);
      if (!result)
         stat_tool::wrap_util::throw_error(error);

    }

    static boost::python::list get_item(const Sequences* seq,
    		boost::python::tuple indexes)
      {
    	int index_var = extract<int>(indexes[1]);
    	int index_seq = extract<int>(indexes[0]);

        // Test index
        if(index_seq<0 || index_seq>=seq->get_nb_sequence())
        {
        	PyErr_SetString(PyExc_IndexError, "sequence index out of bound");
        	boost::python::throw_error_already_set();
        }
        if(index_var<0 || index_var>=seq->get_nb_variable())
        {
             PyErr_SetString(PyExc_IndexError, "variable index out of bound");
             boost::python::throw_error_already_set();
        }

        boost::python::list l;

        int  nb_length = seq->get_length(index_seq);

        for(int index=0; index<nb_length; index++)
        {
        	if ((seq->get_type(index_var) == INT_VALUE) ||
        		(seq->get_type(index_var) == STATE))
        		l.append(seq->get_int_sequence(index_seq, index_var,index));
        	else
        		l.append(seq->get_real_sequence(index_seq, index_var, index));
        }

        return l;
      }


    static Sequences* value_select(const Sequences& seq,  int variable,
                      const object& min, const object& max, bool keep)
     {
       Format_error error;
       Sequences * ret = NULL;

       boost::python::extract<int> get_min(min);
       boost::python::extract<int> get_max(max);

       if (get_min.check() && get_max.check())  // Array of int
       {
    	   int mi = get_min();
    	   int ma = get_max();
    	   ret = seq.value_select(error, variable, mi, ma, keep);
       }
       else
       {
    	   double mi = extract<double>(min);
    	   double ma = extract<double>(max);
    	   ret = seq.value_select(error, variable, mi, ma, keep);
       }

       if(!ret)
         stat_tool::wrap_util::throw_error(error);
       return ret;
     }


    static Sequences* select_variable(const Sequences& seq,
    		const boost::python::list& variables,
            bool keep)
      {
        Format_error error;
        Sequences * ret = NULL;

        int nb_var = len(variables);
        stat_tool::wrap_util::auto_ptr_array<int> vars(new int[nb_var]);

        for (int i=0; i<nb_var; i++)
        	vars[i] = extract<int>(variables[i]);

        ret = seq.select_variable(error, nb_var, vars.get(), keep);

        if(!ret)
          stat_tool::wrap_util::throw_error(error);

        return ret;
      }

    static Sequences* select_individual(const Sequences& seq,
					  const boost::python::list& identifiers,
                      bool keep)
    {
      Format_error error;
      Sequences * ret = NULL;

      int nb_id = len(identifiers);
      stat_tool::wrap_util::auto_ptr_array<int> ids(new int[nb_id]);

      for (int i=0; i<nb_id; i++)
        ids[i] = extract<int>(identifiers[i]);

      ret = seq.select_individual(error, nb_id, ids.get(), keep);

      if(!ret)
        stat_tool::wrap_util::throw_error(error);

      return ret;
    }

    // Shift
    static Sequences* shift(const Sequences& seq, int var, double param)
      {
        Format_error error;
        Sequences * ret = NULL;

        if(seq.get_type(var-1) == REAL_VALUE)
          ret = seq.shift(error, var, (double)param);
        else
          ret = seq.shift(error, var, (int)param);

        if(!ret)
          stat_tool::wrap_util::throw_error(error);

        return ret;
      }


      // Merge
      static Sequences* merge(const Sequences& input_seq,
    		  const boost::python::list& seqs)
      {
        Format_error error;
        Sequences * ret = NULL;

        int nb_seq = len(seqs);
        stat_tool::wrap_util::auto_ptr_array<const Sequences *>
          sequens(new const Sequences*[nb_seq]);

        for (int i=0; i<nb_seq; i++)
          sequens[i] = extract<Sequences*>(seqs[i]);

        ret = input_seq.merge(error, nb_seq, sequens.get());

        if(!ret)
          stat_tool::wrap_util::throw_error(error);

        return ret;
      }

      static Sequences* merge_variable(const Sequences& input_seq,
                     const boost::python::list& seqs, int ref_sample)
      {
        Format_error error;
        Sequences * ret = NULL;

        int nb_seq = len(seqs);
        stat_tool::wrap_util::auto_ptr_array<const Sequences *>
          sequences(new const Sequences*[nb_seq]);

        for (int i=0; i<nb_seq; i++)
          sequences[i] = extract<Sequences*>(seqs[i]);

        ret = input_seq.merge_variable(error, nb_seq, sequences.get(), ref_sample);

        if(!ret)
          stat_tool::wrap_util::throw_error(error);

        return ret;
      }


      // Cluster
      static Sequences* cluster_step(const Sequences& seq, int variable, int step)
      {
        Format_error error;
        Sequences* ret = seq.cluster(error, variable, step);

        if(!ret)
          stat_tool::wrap_util::throw_error(error);

        return ret;
      }

      static Sequences* cluster_limit(const Sequences& seq, int variable,
                    boost::python::list& limit
                    )
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
            ret = seq.cluster(error, variable, nb_limit, ldouble);
            delete[] ldouble;
        }
        else
        {
        	ret = seq.cluster(error, variable, nb_limit, lint);
        	delete[] lint;
        }

        if(!ret)
          stat_tool::wrap_util::throw_error(error);

        return ret;
      }


      static Sequences* transcode(const Sequences& seq, int variable,
                    boost::python::list& symbol
                          )
      {

        Format_error error;

        int nb_symbol = len(symbol);
        stat_tool::wrap_util::auto_ptr_array<int>
          l(new int[nb_symbol]);

        int expected_nb_symbol = (int)(seq.get_max_value(variable - 1)
                       - seq.get_min_value(variable - 1)) + 1;

        if(nb_symbol != expected_nb_symbol)
          stat_tool::wrap_util::throw_error("Bad number of Symbol");

        for (int i=0; i<nb_symbol; i++)
          l[i] = extract<int>(symbol[i]);

        Sequences* ret = seq.transcode(error, variable, l.get());


        if(!ret)
          stat_tool::wrap_util::throw_error(error);

        return ret;
      }




};



// Boost declaration

void class_sequences()
{
  class_< Sequences, bases< STAT_interface > >
    ("_Sequences", "Sequences")
    .def("__init__", make_constructor(SequencesWrap::sequences_from_file))
    .def("__init__", make_constructor(SequencesWrap::build_from_lists))

    // Python Operators
    .def(self_ns::str(self)) //__str__
    .def("__len__", &Sequences::get_nb_sequence,
    		"Returns number of sequences")
    .def("__getitem__", SequencesWrap::get_item)


    .def("get_nb_sequence", &Sequences::get_nb_sequence,
		"Return the number of sequences")
    .def("get_nb_variable", &Sequences::get_nb_variable,
    		"Return the number of variables")

    // Identifiers
    .def("get_identifiers", &SequencesWrap::get_identifiers,
    		"Return list of identifiers")

    .def("extract", SequencesWrap::extract_histogram,
    		return_value_policy< manage_new_object >(),
    		python::args("variable"),
    		"Extract Histogram")

    // Select
    .def("value_select", SequencesWrap::value_select,
     return_value_policy< manage_new_object >(),
     python::args("variable", "min", "max", "keep"),
     "Selection of individuals according to the values taken by a variable")

    .def("select_variable", SequencesWrap::select_variable,
     return_value_policy< manage_new_object >(),
     python::args("variables", "keep"),
     "select variable given a list of index")

    .def("select_individual", SequencesWrap::select_individual,
     return_value_policy< manage_new_object >(),
     python::args("identifiers", "keep"),
     "Select individuals given a list of identifiers")

    // Merge
   .def("merge", SequencesWrap::merge,
    return_value_policy< manage_new_object >(),
    "Merge sequences")

   .def("merge_variable", SequencesWrap::merge_variable,
    return_value_policy< manage_new_object >(),
    "Merge variables" )


     // Cluster
    .def("cluster_step", SequencesWrap::cluster_step,
     return_value_policy< manage_new_object >(),
     python::args("variable", "step"),
     "Cluster Step"
     )

    .def("cluster_limit", SequencesWrap::cluster_limit,
     return_value_policy< manage_new_object >(),
     python::args("variable", "limits"),
     "Cluster limit")

    .def("transcode", SequencesWrap::transcode,
     return_value_policy< manage_new_object >(),
     python::args("variable", "symbols"),
     "Transcode")



    // Output
    .def("ascii_data_write", SequencesWrap::ascii_data_write,
    "Return a string with the object representation")

    // save to file
    .def("file_ascii_write", SequencesWrap::file_ascii_write,
    "Save vector summary into a file")

    .def("file_ascii_data_write", SequencesWrap::file_ascii_data_write,
    "Save vector data into a file")

    // Shift
    .def("shift", SequencesWrap::shift,
     return_value_policy< manage_new_object >(),
     python::args("variable", "param"),
     "Shift")
    ;
}

