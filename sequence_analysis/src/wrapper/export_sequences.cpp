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


#include "boost_python_aliases.h"

using namespace boost::python;
using namespace boost;
using namespace stat_tool;

class SequencesWrap {

public:

	static boost::shared_ptr<Sequences> sequences_from_file(char* filename) {
		Format_error error;
		Sequences *sequences = NULL;
		bool old_format = false;

		sequences = sequences_ascii_read(error, filename, old_format);

		/*    if(!top_parameters)
		 {
		 stat_tool::wrap_util::throw_error(error);
		 }
		 */
		return boost::shared_ptr<Sequences>(sequences);
	}


	static Sequences* build_from_lists(boost::python::list& input_list,
			boost::python::list& input_identifiers,
			int input_index_parameter_type) {

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
		if (nb_sequences!=nb_identifiers){
			PyErr_SetString(PyExc_IndexError, "number of identifiers must be equal to number of sequences");
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

		if (!get_list.check()) {
			boost::python::list seq0 = extract<boost::python::list> (
					input_list[0]);
			object elt = seq0[0];
			extract<float> get_float(elt);
			extract<int> get_int(elt);

			if (get_int.check())
				is_int = true;
			else if (get_float.check())
				is_float = true;
		} else {

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
		if (is_float){
			real_sequence = new double**[nb_sequences];
			type = REAL_VALUE;
		}
		if (is_int){
			int_sequence = new int**[nb_sequences];
			type = INT_VALUE;
		}

		//cerr << is_float << endl;
		//	cerr << is_int << endl;

		// for each sequence, are we considering case 1 or 2 ?
		for (int seqi = 0; seqi < nb_sequences; seqi++) {
			// whatever case it is, the identifiers can be set here
			identifiers[seqi] = extract<int>(input_identifiers[seqi]);

			// try to get a single sequence
			boost::python::list sequence = extract<boost::python::list> (
					input_list[seqi]);

			// is it case 2 i.e. there is another nested list ?
			extract<boost::python::list> get_list(sequence[0]);

			// if not a list, we are in the case 1
			if (!get_list.check()) {
				// the length of the sequence is simply:
				int N = len(sequence);
				length[seqi] = N;

				// and there is only 1 variable
				nb_variables = 1; // used by the constructor

				// so the data structure is as follows
				if (is_float) {
					real_sequence[seqi] = new double*[1];
					real_sequence[seqi][0] = new double[N];
				} else {
					int_sequence[seqi] = new int*[1];
					int_sequence[seqi][0] = new int[N];
				}

				// and we can populate the output using the sequence as vector
				for (int j = 0; j < N; j++) {
					if (is_float) {
						double v = extract<float> (sequence[j]);
						real_sequence[seqi][0][j] = v;
					} else {
						int v = extract<int> (sequence[j]);
						int_sequence[seqi][0][j] = v;
					}
				}
			}
			// if we can get the data, we have a list and therefore we are in case 2
			else {
				nb_variables = len(sequence); // used by the constructor

				// what are the number of variables in this constructor :
				int N_var = len(sequence);

				// loop over the variable to extract their contents
				for (int vari = 0; vari < N_var; vari++) {
					// get the variable of the current sequence
					boost::python::list variable =
							extract<boost::python::list> (sequence[vari]);

					// and update the data structure
					int N = len(variable);
					length[seqi] = N;

					//real_sequence[seqi][vari] = new double[N];
					//		cerr << "---"<<N<< endl;

					// before extracting the vector of the

					if (is_float) {
						//		cerr << " is float validated"<<endl;
						if (vari == 0)
							real_sequence[seqi] = new double*[N_var]; // only once
						real_sequence[seqi][vari] = new double[N];
					}
					if (is_int) {
						//		cerr << " is int chjosen"<<endl;
						if (vari == 0)
							int_sequence[seqi] = new int*[N_var]; // only once
						int_sequence[seqi][vari] = new int[N];
					}

					// get the vector corresponding to a sequence/variable pair
					for (int j = 0; j < N; j++) {
						if (is_float) {
							float v = extract<float> (variable[j]);
							real_sequence[seqi][vari][j] = v;
						}

						if (is_int) {
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

		if (is_float) {
			sequences = new Sequences(nb_sequences, identifiers, length,
					nb_variables, real_sequence);

			if (real_sequence) {
				for (int i = 0; i < nb_sequences; i++) {
					for (int j = 0; j < nb_variables; j++) {
						delete[] real_sequence[i][j];
					}
					delete[] real_sequence[i];
				}
				delete[] real_sequence;
			}
		}

		if (is_int) {
			sequences = new Sequences(nb_sequences, identifiers, length,
					index_parameter_type, nb_variables, type, int_sequence);
			// freeing memory

			if (int_sequence) {
				for (int i = 0; i < nb_sequences; i++) {
					for (int j = 0; j < nb_variables; j++) {
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

	static Distribution_data* extract_histogram(const Sequences& seq,
			int variable) {
		Format_error error;
		Distribution_data *ret = NULL;

		ret = seq.extract(error, variable);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;

	}

	static boost::python::list get_identifiers(const Sequences& seq) {
		boost::python::list l;

		int nb_seq = seq.get_nb_sequence();

		for (int s = 0; s < nb_seq; s++) {
			l.append(seq.get_identifier(s));
		}

		return l;
	}

	static std::string ascii_data_write(const Sequences& d, bool exhaustive) {
		std::stringstream s;
		std::string res;

		d.ascii_data_write(s, exhaustive, true);
		res = s.str();
		return res;
	}

	static void file_ascii_data_write(const Sequences& d, const char* path,
			bool exhaustive) {
		bool result = true;
		Format_error error;

		result = d.ascii_data_write(error, path, exhaustive);
		if (!result)
			stat_tool::wrap_util::throw_error(error);

	}

	static void file_ascii_write(const Sequences& d, const char* path,
			bool exhaustive) {
		bool result = true;
		Format_error error;

		result = d.ascii_write(error, path, exhaustive);
		if (!result)
			stat_tool::wrap_util::throw_error(error);

	}

	static boost::python::list get_item(const Sequences* seq,
			boost::python::tuple indexes) {
		int index_var = extract<int> (indexes[1]);
		int index_seq = extract<int> (indexes[0]);

		// Test index
		if (index_seq < 0 || index_seq >= seq->get_nb_sequence()) {
			PyErr_SetString(PyExc_IndexError, "sequence index out of bound");
			boost::python::throw_error_already_set();
		}
		if (index_var < 0 || index_var >= seq->get_nb_variable()) {
			PyErr_SetString(PyExc_IndexError, "variable index out of bound");
			boost::python::throw_error_already_set();
		}

		boost::python::list l;

		int nb_length = seq->get_length(index_seq);

		for (int index = 0; index < nb_length; index++) {
			if ((seq->get_type(index_var) == INT_VALUE) || (seq->get_type(
					index_var) == STATE))
				l.append(seq->get_int_sequence(index_seq, index_var, index));
			else
				l.append(seq->get_real_sequence(index_seq, index_var, index));
		}

		return l;
	}

	static Sequences* value_select(const Sequences& seq, int variable,
			const object& min, const object& max, bool keep) {
		Format_error error;
		Sequences * ret = NULL;

		boost::python::extract<int> get_min(min);
		boost::python::extract<int> get_max(max);

		if (get_min.check() && get_max.check()) // Array of int
		{
			int mi = get_min();
			int ma = get_max();
			ret = seq.value_select(error, variable, mi, ma, keep);
		} else {
			double mi = extract<double> (min);
			double ma = extract<double> (max);
			ret = seq.value_select(error, variable, mi, ma, keep);
		}

		if (!ret)
			stat_tool::wrap_util::throw_error(error);
		return ret;
	}

	static Sequences* select_variable(const Sequences& seq,
			const boost::python::list& variables, bool keep) {
		Format_error error;
		Sequences * ret = NULL;

		int nb_var = len(variables);
		stat_tool::wrap_util::auto_ptr_array<int> vars(new int[nb_var]);

		for (int i = 0; i < nb_var; i++)
			vars[i] = extract<int> (variables[i]);

		ret = seq.select_variable(error, nb_var, vars.get(), keep);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	static Sequences* select_individual(const Sequences& seq,
			const boost::python::list& identifiers, bool keep) {
		Format_error error;
		Sequences * ret = NULL;

		int nb_id = len(identifiers);
		stat_tool::wrap_util::auto_ptr_array<int> ids(new int[nb_id]);

		for (int i = 0; i < nb_id; i++)
			ids[i] = extract<int> (identifiers[i]);

		ret = seq.select_individual(error, nb_id, ids.get(), keep);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	// Shift
	static Sequences* shift(const Sequences& seq, int var, double param) {
		Format_error error;
		Sequences * ret = NULL;

		if (seq.get_type(var - 1) == REAL_VALUE)
			ret = seq.shift(error, var, (double) param);
		else
			ret = seq.shift(error, var, (int) param);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	// Merge
	static Sequences* merge(const Sequences& input_seq,
			const boost::python::list& seqs) {
		Format_error error;
		Sequences * ret = NULL;

		int nb_seq = len(seqs);
		stat_tool::wrap_util::auto_ptr_array<const Sequences *> sequens(
				new const Sequences*[nb_seq]);

		for (int i = 0; i < nb_seq; i++)
			sequens[i] = extract<Sequences*> (seqs[i]);

		ret = input_seq.merge(error, nb_seq, sequens.get());

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	static Sequences* merge_variable(const Sequences& input_seq,
			const boost::python::list& seqs, int ref_sample) {
		Format_error error;
		Sequences * ret = NULL;

		int nb_seq = len(seqs);
		stat_tool::wrap_util::auto_ptr_array<const Sequences *> sequences(
				new const Sequences*[nb_seq]);

		for (int i = 0; i < nb_seq; i++)
			sequences[i] = extract<Sequences*> (seqs[i]);

		ret = input_seq.merge_variable(error, nb_seq, sequences.get(),
				ref_sample);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	// Cluster
	static Sequences* cluster_step(const Sequences& seq, int variable, int step) {
		Format_error error;
		Sequences* ret = seq.cluster(error, variable, step);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	static Sequences* cluster_limit(const Sequences& seq, int variable,
			boost::python::list& limit) {

		Format_error error;

		int nb_limit = len(limit);
		bool is_float = true;
		int *lint = NULL;
		double *ldouble = NULL;
		Sequences* ret;

		// Test type
		boost::python::extract<int> get_int(limit[0]);
		if (get_int.check()) {
			is_float = false;
			lint = new int[nb_limit];
		} else {
			ldouble = new double[nb_limit];
		}

		// Convert list
		for (int i = 0; i < nb_limit; i++) {
			if (is_float)
				ldouble[i] = extract<int> (limit[i]);
			else
				lint[i] = extract<double> (limit[i]);
		}

		// Call correct function
		if (is_float) {
			ret = seq.cluster(error, variable, nb_limit, ldouble);
			delete[] ldouble;
		} else {
			ret = seq.cluster(error, variable, nb_limit, lint);
			delete[] lint;
		}

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	static Sequences* transcode(const Sequences& seq, int variable,
			boost::python::list& symbol) {

		Format_error error;

		int nb_symbol = len(symbol);
		stat_tool::wrap_util::auto_ptr_array<int> l(new int[nb_symbol]);

		int expected_nb_symbol = (int) (seq.get_max_value(variable - 1)
				- seq.get_min_value(variable - 1)) + 1;

		if (nb_symbol != expected_nb_symbol)
			stat_tool::wrap_util::throw_error("Bad number of Symbol");

		for (int i = 0; i < nb_symbol; i++)
			l[i] = extract<int> (symbol[i]);

		Sequences* ret = seq.transcode(error, variable, l.get());

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	// Reverse
	static Sequences* reverse(const Sequences& seq) {

		Format_error error;

		Sequences* ret = seq.reverse(error);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);
		return ret;
	}

	static Sequences* length_select(const Sequences& seq, int min_length,
			int max_length, bool keep) {

		Format_error error;

		Sequences* ret;

		ret = seq.length_select(error, min_length, max_length, keep);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}


	static Sequences* remove_run(const Sequences& seq, int variable,
			int ivalue, char position , int max_run_length) {
		Format_error error;

		Sequences* ret;

		ret = seq.remove_run(error, variable, ivalue,
			position, max_run_length);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	static int get_type(const Sequences& seq, int index) {
		if (index <0 or index>=seq.get_nb_variable()){
				PyErr_SetString(PyExc_IndexError,
					"index must be positive and less than number of variables");
					boost::python::throw_error_already_set();
		}

		return seq.get_type(index);
	}

	static int get_length(const Sequences& seq, int index) {
		if (index <0 or index>=seq.get_nb_sequence()){
				PyErr_SetString(PyExc_IndexError,
					"index must be positive and less than number of sequences");
					boost::python::throw_error_already_set();
		}

		return seq.get_length(index);
	}

	static double get_min_value(const Sequences& seq, int variable) {
		if (variable <0 or variable>=seq.get_nb_variable()){
				PyErr_SetString(PyExc_IndexError,
					"index must be positive and less than number of variables");
					boost::python::throw_error_already_set();
		}
		return seq.get_min_value(variable);
	}

	static double get_max_value(const Sequences& seq, int variable) {
		if (variable <0 or variable>=seq.get_nb_variable()){
				PyErr_SetString(PyExc_IndexError,
					"index must be positive and less than number of variables");
					boost::python::throw_error_already_set();
		}
		return seq.get_max_value(variable);
	}


	static Sequences* scaling(const Sequences& seq, int variable,
			int scaling_coeff) {

		Format_error error;
		Sequences* ret;

		ret = seq.scaling(error, variable, scaling_coeff);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	static Sequences* round(const Sequences& seq, int variable,
			int mode) {

		Format_error error;
		Sequences* ret;

		ret = seq.round(error, variable, mode);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	//cumulate
	static Sequences* cumulate(const Sequences& seq, int variable) {

		Format_error error;
		Sequences* ret;

		ret = seq.cumulate(error, variable);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	//differece
	static Sequences* difference(const Sequences& seq, int variable, bool first) {

		Format_error error;
		Sequences* ret;

		ret = seq.difference(error, variable, first);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}


	static Sequences* index_parameter_select(const Sequences& seq,
		int min_index_parameter, int max_index_parameter, bool keep) {

		Format_error error;
		Sequences* ret;

		ret = seq.index_parameter_select(error, min_index_parameter, max_index_parameter, keep);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	static Sequences* remove_index_parameter(const Sequences& seq){

		Format_error error;
		Sequences* ret;

		ret = seq.remove_index_parameter(error);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	static Sequences* index_parameter_extract(const Sequences& seq,
			int min_index_parameter, int max_index_parameter){

		Format_error error;
		Sequences* ret;

		ret = seq.index_parameter_extract(error, min_index_parameter, max_index_parameter);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	static Sequences* segmentation_extract(const Sequences& seq,
			int variable, int nb_value,
			boost::python::list& input_values, bool keep){

		Format_error error;
		Sequences* ret;
		int *values;

		values = new int[nb_value];


		for (int i = 0; i < nb_value; i++) {
			values[i] = extract<double> (input_values[i]);
		}

		ret = seq.segmentation_extract(error,
				variable, nb_value, values, keep);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		delete [] values;

		return ret;
	}

	static Sequences* moving_average(const Sequences& seq,
			 int nb_value , boost::python::list& input_values ,
			 int variable , bool begin_end , int output){

		Format_error error;
		Sequences* ret;

		double *values;
		values = new double[nb_value];

		for (int i = 0; i < nb_value; i++) {
			values[i] = extract<double> (input_values[i]);
		}

		ret = seq.moving_average(error,
				nb_value, values, variable, begin_end, output);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		delete [] values;

		return ret;
	}

	static Sequences* pointwise_average(const Sequences& seq,
		bool standard_deviation ,  int output , const char *path , char format){

		Format_error error;
		Sequences* ret;

		ret = seq.pointwise_average(error,
				standard_deviation, output, path, format);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	static Sequences* recurrence_time_sequences(const Sequences& seq,
			int variable, int value){

		Format_error error;
		Sequences* ret;

		ret = seq.recurrence_time_sequences(error, variable, value);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	static Sequences* cross(const Sequences& seq){

		Format_error error;
		Sequences* ret;

		ret = seq.cross(error);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}


	static Sequences* transform_position(const Sequences& seq, int step)
	{
		Format_error error;
		Sequences* ret;

		ret = seq.transform_position(error, step);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	static Sequences* sojourn_time_sequences(const Sequences& seq, int variable)
	{
		Format_error error;
		Sequences* ret;

		ret = seq.sojourn_time_sequences(error, variable);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

};

// Boost declaration

void class_sequences() {

	enum_<stat_tool::wrap_util::UniqueInt<5, 10> >("IndexParameterType")
	 .value("IMPLICIT_TYPE", IMPLICIT_TYPE)
	 .value("TIME", TIME)
	 .value("TIME_INTERVAL", TIME_INTERVAL)
	 .value("POSITION", POSITION)
	 .value("POSITION_INTERVAL", POSITION_INTERVAL)
	 .export_values();

        enum_<stat_tool::wrap_util::UniqueInt<6, 11> >("Type")
	 .value("INT_VALUE", INT_VALUE)
	 .value("REAL_VALUE", REAL_VALUE)
	 .value("STATE", STATE)
	 .value("OLD_INT_VALUE", OLD_INT_VALUE)
	 .value("NB_INTERNODE", NB_INTERNODE)
	 .value("AUXILIARY",   AUXILIARY)
	 .export_values();

	class_<Sequences, bases<STAT_interface> > ("_Sequences", "Sequences")
	.def(init <const Renewal_data&>())
	.def("__init__", make_constructor(SequencesWrap::sequences_from_file))
	.def("__init__", make_constructor(SequencesWrap::build_from_lists))

	// Python Operators
	.def(self_ns::str(self)) //__str__
	.def("__len__", &Sequences::get_nb_sequence,"Returns number of sequences")
	.def("__getitem__", SequencesWrap::get_item)

	.def("get_nb_sequence", &Sequences::get_nb_sequence, "Return the number of sequences")
	.def("get_nb_variable", &Sequences::get_nb_variable, "Return the number of variables")
	.def("get_max_length", &Sequences::get_max_length,"Return max length")
	.def("get_cumul_length", &Sequences::get_cumul_length,"Return cumul length")


	.def("get_length", &SequencesWrap::get_length,
		python::args("index_seq"), "return length of a sequence")

	.def("get_index_parameter_type", &Sequences::get_index_parameter_type,
		"return index parameter type")
	.def("get_type", &SequencesWrap::get_type,
		python::args("index_var"),
		"return  type")
	.def("get_min_value", &Sequences::get_min_value,
		python::args("index_var"), "return min value of variables")
	.def("get_max_value", &Sequences::get_max_value,
		python::args("index_var"), "return max value of variables")



	// Identifiers
	.def("get_identifiers", &SequencesWrap::get_identifiers,
			"Return list of identifiers")
	// Select
	.def("value_select", SequencesWrap::value_select,
		return_value_policy<manage_new_object> (),
		python::args("variable", "min", "max", "keep"),
		"Selection of individuals according to the values taken by a variable")
	.def("select_variable", SequencesWrap::select_variable,
		return_value_policy<manage_new_object> (),
              python::args("variables", "keep"),
                 "select variable given a list of index")
        .def("select_individual", SequencesWrap::select_individual,
		return_value_policy<manage_new_object> (),
		python::args("identifiers", "keep"),
		"Select individuals given a list of identifiers")
	// parameters
	.def("index_parameter_select", SequencesWrap::index_parameter_select,
		return_value_policy<manage_new_object> (),
		python::args("min_index_parameter", "max_index_parameter", "keep"),
		"Select sequences in an index parameter range")
	.def("index_parameter_extract", SequencesWrap::index_parameter_extract,
		return_value_policy<manage_new_object> (),
		python::args("min_index_parameter", "max_index_parameter"),
		"Select sequences in an index parameter range")
        .def("remove_index_parameter", SequencesWrap::remove_index_parameter,
		return_value_policy<manage_new_object> (),
		"Remove index parameter")
	// Reverse
	.def("reverse", SequencesWrap::reverse,
		return_value_policy<manage_new_object> (),
		"reverse")
	// Extract
	.def("extract", SequencesWrap::extract_histogram, return_value_policy<
		manage_new_object> (), python::args("variable"),
		"Extract histogram")
	// Merge
	.def("merge", SequencesWrap::merge,
			return_value_policy<manage_new_object> (), "Merge sequences")
	.def("merge_variable", SequencesWrap::merge_variable, return_value_policy<
			manage_new_object> (), "Merge variables")
	// Cluster
	.def("cluster_step", SequencesWrap::cluster_step, return_value_policy<
			manage_new_object> (), python::args("variable", "step"),
			"Cluster Step")
	.def("cluster_limit", SequencesWrap::cluster_limit, return_value_policy<
			manage_new_object> (), python::args("variable", "limits"),
			"Cluster limit")
	//Transcode
	.def("transcode", SequencesWrap::transcode, return_value_policy<
			manage_new_object> (), python::args("variable", "symbols"),
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
		return_value_policy<manage_new_object> (),
		python::args("variable","param"),
		"Shift")
	//remove run
	.def("remove_run",SequencesWrap::remove_run,
			return_value_policy<manage_new_object> (),
			python::args("variable","ivalue", "position", "max_run_length"),
			"RemoveRun")
	//length select
	.def("length_select",SequencesWrap::length_select,
		return_value_policy<manage_new_object> (),
		python::args("min_length", "max_length", "keep"),
		"RemoveRun")
	//scaling
	.def("scaling", SequencesWrap::scaling,
		return_value_policy<manage_new_object> (),
		python::args("variable", "scaling_coeff"),
		"scaling variable")
	//round
	.def("round", SequencesWrap::round,
		return_value_policy<manage_new_object> (),
		python::args("variable", "mode"),
		"round variable")
	//segmentation extract
	.def("segmentation_extract",SequencesWrap::segmentation_extract,
		return_value_policy<manage_new_object> (),
		python::args("variable", "nb_values", "values", "keep"),
		"Segmentation extract")
	//cumulate
	.def("cumulate", SequencesWrap::cumulate,
		return_value_policy<manage_new_object> (),
		python::args("variable"),
		"Cumulate")
	//difference
	.def("difference", SequencesWrap::difference,
		return_value_policy<manage_new_object> (),
		python::args("variable", "first_element"),
		"Difference")
	//moving average
	//Filtrage de type moyenne mobile des sequences.
	.def("moving average", SequencesWrap::moving_average,
		return_value_policy<manage_new_object> (),
		python::args("nb_point" ,"filter" , "variable" , "begin_end" , "output"),
		"Moving average ")
	//pointwise average
	.def("pointwise_average", SequencesWrap::pointwise_average,
		return_value_policy<manage_new_object> (),
		python::args("standard_deviation", "output", "path", "format"),
		"Pointwise average")

	DEF_RETURN_VALUE("recurrence_time_sequences", SequencesWrap::recurrence_time_sequences,ARGS("variable", "value"),"Recurrence time sequences")
	DEF_RETURN_VALUE("sojourn_time_sequences", SequencesWrap::sojourn_time_sequences, ARGS("variable"), "Sojourn time sequences")
	DEF_RETURN_VALUE("transform position", SequencesWrap::transform_position, ARGS("step"), "Transform position")
	DEF_RETURN_VALUE_NO_ARGS("cross", SequencesWrap::cross, "Cross")
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


	    Vectors* build_vectors(bool index_variable) const;
	    Vectors* extract_vectors(Format_error &error , int feature_type , int variable = I_DEFAULT ,int value = I_DEFAULT) const;

	    Markovian_sequences* markovian_sequences(Format_error &error) const;
	    Tops* tops(Format_error &error) const;

	    bool check(Format_error &error , const char *pattern_label);

	    Time_events* extract_time_events(Format_error &error , int variable , int begin_date , int end_date ,  int previous_date = I_DEFAULT , int next_date = I_DEFAULT) const;
	    Renewal_data* extract_renewal_data(Format_error &error , int variable , int begin_index , int end_index) const;

	    Sequences* cluster(Format_error &error , int variable , int step ,int mode = FLOOR) const;
	    Sequences* index_parameter_select(Format_error &error , int min_index_parameter ,    int max_index_parameter , bool keep) const;
	    Sequences* remove_index_parameter(Format_error &error) const;
	    Sequences* moving_average(Format_error &error , int nb_point , double *filter ,   int variable = I_DEFAULT , bool begin_end = false ,	                              int output = TREND) const;
	    Sequences* moving_average(Format_error &error , const Distribution &dist ,   int variable = I_DEFAULT , bool begin_end = false ,int output = TREND) const;


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

	    Correlation* correlation_computation(Format_error &error , int variable1 , int variable2 ,
	                                         int itype = PEARSON , int max_lag = I_DEFAULT ,
	                                         int normalization = EXACT) const;
	    Correlation* partial_autocorrelation_computation(Format_error &error , int variable ,
	                                                     int itype = PEARSON , int max_lag = I_DEFAULT) const;

	    Distance_matrix* alignment(Format_error &error , std::ostream *os , const Vector_distance &ivector_dist ,
	                               int ref_identifier = I_DEFAULT , int test_identifier = I_DEFAULT ,
	                               bool begin_free = false , bool end_free = false , int indel_cost = ADAPTATIVE ,
	                               double indel_factor = INDEL_FACTOR_1 , bool transposition_flag = false ,
	                               double transposition_factor = TRANSPOSITION_FACTOR ,
	                               const char *result_path = 0 , char result_format = 'a' ,
	                               const char *alignment_path = 0 , char alignment_format = 'a') const;
	    Distance_matrix* alignment(Format_error &error , std::ostream *os , int ref_identifier = I_DEFAULT ,
	                               int test_identifier = I_DEFAULT , bool begin_free = false , bool end_free = false ,
	                               const char *result_path = 0 , char result_format = 'a' ,
	                               const char *alignment_path = 0 , char alignment_format = 'a') const;

	    Sequences* multiple_alignment(Format_error &error , std::ostream &os , const Vector_distance &ivector_dist ,
	                                  bool begin_free = false , bool end_free = false , int indel_cost = ADAPTATIVE ,
	                                  double indel_factor = INDEL_FACTOR_N , int algorithm = AGGLOMERATIVE ,
	                                  const char *path = 0) const;

	    Sequences* segmentation(Format_error &error , std::ostream &os , int iidentifier ,
	                            int nb_segment , int *ichange_point , int *model_type ,
	                            int output = SEQUENCE) const;
	    Sequences* segmentation(Format_error &error , std::ostream &os , int *nb_segment ,
	                            int *model_type , int iidentifier = I_DEFAULT ,
	                            int output = SEQUENCE) const;
	    Sequences* segmentation(Format_error &error , std::ostream &os , int iidentifier ,
	                            int max_nb_segment , int *model_type) const;

	    Sequences* hierarchical_segmentation(Format_error &error , std::ostream &os , int iidentifier ,
	                                         int max_nb_segment , int *model_type) const;

	    Sequences* segmentation(Format_error &error , int iidentifier , int nb_segment ,
	                            const Vector_distance &ivector_dist , std::ostream &os ,
	                            int output = SEGMENT) const;

	    bool segment_profile_write(Format_error &error , std::ostream &os , int iidentifier ,
	                               int nb_segment , int *model_type , int output = SEGMENT ,
	                               char format = 'a' , int segmentation = FORWARD_DYNAMIC_PROGRAMMING ,
	                               int nb_segmentation = NB_SEGMENTATION) const;
	    bool segment_profile_write(Format_error &error , const char *path , int iidentifier ,
	                               int nb_segment , int *model_type , int output = SEGMENT ,
	                               char format = 'a' , int segmentation = FORWARD_DYNAMIC_PROGRAMMING ,
	                               int nb_segmentation = NB_SEGMENTATION) const;
	    bool segment_profile_plot_write(Format_error &error , const char *prefix ,
	                                    int iidentifier , int nb_segment , int *model_type ,
	                                    int output = SEGMENT , const char *title = 0) const;

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
class WRAP {
public:
  static boost::python::list get_initial_run(Sequence_characteristics& input)
  {

    boost::python::list list;

    for (int i=0; i< input.get_nb_value(); i++)
      {
        //list.append(  boost::python::extract<Histogram>(input.get_initial_run(i)));
        list.append(  i);

      }

    return list;
  }

  static Histogram* get_initial_run_from_index(Sequence_characteristics& input, int index)
  {
      return input.get_initial_run(index);
  }
};

void class_sequence_characteristics() {

    class_<Sequence_characteristics> ("_SequenceCharacteristics", "SequenceCharacteristics")
    .def(init<optional<int> >())
    .def(init<Sequence_characteristics, optional<bool> >())
    .def(init<Sequence_characteristics, optional<char> >())
    .def("get_nb_value", &Sequence_characteristics::get_nb_value)
    DEF_RETURN_VALUE_NO_ARGS("get_index_value", &Sequence_characteristics::get_index_value, "get_index_value")
    DEF_RETURN_VALUE_NO_ARGS("get_first_occurrence", &Sequence_characteristics::get_first_occurrence, "get first occurrence time")
    DEF_RETURN_VALUE_NO_ARGS("get_recurrence_time", &Sequence_characteristics::get_recurrence_time, "get recurrence time")
    DEF_RETURN_VALUE_NO_ARGS("get_sojourn_time", &Sequence_characteristics::get_sojourn_time,"returns sojourn time")
    .def("get_initial_run", &WRAP::get_initial_run, "returns initial run")
    DEF_RETURN_VALUE("get_initial_run_from_index", WRAP::get_initial_run_from_index,ARGS("index"), "returns initial run")
    DEF_RETURN_VALUE_NO_ARGS("get_final_run", &Sequence_characteristics::get_final_run,"returns final run")
    DEF_RETURN_VALUE_NO_ARGS("get_nb_run", &Sequence_characteristics::get_nb_run, "returns number of run")
    DEF_RETURN_VALUE_NO_ARGS("get_nb_occurrence", &Sequence_characteristics::get_nb_occurrence, "returns number of ocurrences")
    ;

    /*      Histogram** get_initial_run() const { return initial_run; } */
}
#undef WRAP
