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


  static Sequences* build_from_lists(boost::python::list& array)
  {
	int nb_sequences = boost::python::len(array);
	Sequences *sequences = 0;
	int nb_variables = 0;
	bool three_nested_lists;

	cerr <<"nb sequence="<< nb_sequences <<endl;
	std::flush(cerr);

	boost::python::list test = extract<boost::python::list>(extract<boost::python::list>(array));

	std::flush(cerr);

	// for each sequence
	for (int seqi=0; seqi<nb_sequences; seqi++)
	{
		boost::python::list variables =
			boost::python::extract<boost::python::list>(array[seqi]);

		nb_variables = boost::python::len(variables);

		cerr <<"nb variable="<< nb_variables <<endl;
		std::flush(cerr);

		// for each variable
		for (int vari=0; vari<nb_variables; vari++)
		{
			cerr <<"vari="<< vari <<endl;
			std::flush(cerr);

			try
			{
				boost::python::list vector = extract<boost::python::list>(variables[vari]);
				three_nested_lists = true;
				cerr <<"vector="<<len(vector)<<endl;
				std::flush(cerr);
			}
			catch(...)
			{
				three_nested_lists = false;
				cerr << "length="<<boost::python::len(variables[vari])<<endl;
				std::flush(cerr);
			}




		}
	}

	sequences = new Sequences(nb_sequences, nb_variables);
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

  static boost::python::list get_identifiers(const Sequences& seq, int iseq)
    {
      boost::python::list l;

      int  nb_seq = seq.get_nb_sequence();
      for(int s=0; s<nb_seq; s++)
      {
        l.append(seq.get_identifier(iseq));
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

};



// Boost declaration

void class_sequences()
{
  class_< Sequences, bases< STAT_interface > >
    ("_Sequences", "Sequences")
    .def("__init__", make_constructor(SequencesWrap::sequences_from_file))
    .def("__init__", make_constructor(SequencesWrap::build_from_lists))


    .def(self_ns::str(self)) //__str__
    .def("__len__", &Sequences::get_nb_sequence)
    .def("get_nb_sequence", &Sequences::get_nb_sequence)
    .def("get_nb_variable", &Sequences::get_nb_variable)

    .def("get_identifiers", &Sequences::get_identifier,
    "return list of identifiers")

    .def("extract", SequencesWrap::extract_histogram,
    return_value_policy< manage_new_object >(),
    python::args("variable"),
    "Extract Histogram")

    // Output
    .def("ascii_data_write", SequencesWrap::ascii_data_write,
    "Return a string with the object representation")

    // save to file
    .def("file_ascii_write", SequencesWrap::file_ascii_write,
    "Save vector summary into a file")

    .def("file_ascii_data_write", SequencesWrap::file_ascii_data_write,
    "Save vector data into a file")

    ;
}

