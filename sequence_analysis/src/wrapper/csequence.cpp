
// Includes ====================================================================

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/convolution.h"
#include "stat_tool/compound.h"
#include "stat_tool/discrete_mixture.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"

#include "sequence_analysis/sequences.h"
#include "sequence_analysis/renewal.h"

#include <boost/python.hpp>
// #include <boost/python/detail/api_placeholder.hpp>
// definition of boost::python::len
// #include <boost/python/make_constructor.hpp>
// definition of boost::python::make_constructor

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

template<int num, int id> struct UniqueInt
{
   int v;
   enum { value=num };

   UniqueInt(int _v) : v(_v) { }
   operator int() const { return v; }
};


struct Parametric_Wrapper: DiscreteParametric
{
    Parametric_Wrapper(PyObject* self_):
        DiscreteParametric(), self(self_) {}

    Parametric_Wrapper(PyObject* self_, int p0):
        DiscreteParametric(p0), self(self_) {}

    Parametric_Wrapper(PyObject* self_, int p0, int p1):
        DiscreteParametric(p0, p1), self(self_) {}

    Parametric_Wrapper(PyObject* self_, int p0, int p1, int p2):
        DiscreteParametric(p0, p1, p2), self(self_) {}

    Parametric_Wrapper(PyObject* self_, int p0, int p1, int p2, int p3):
        DiscreteParametric(p0, p1, p2, p3), self(self_) {}

    Parametric_Wrapper(PyObject* self_, int p0, int p1, int p2, int p3, double p4):
        DiscreteParametric(p0, p1, p2, p3, p4), self(self_) {}

    Parametric_Wrapper(PyObject* self_, int p0, int p1, int p2, int p3, double p4, double p5):
        DiscreteParametric(p0, p1, p2, p3, p4, p5), self(self_) {}

    Parametric_Wrapper(PyObject* self_, int p0, int p1, int p2, double p3, double p4):
        DiscreteParametric(p0, p1, p2, p3, p4), self(self_) {}

    Parametric_Wrapper(PyObject* self_, int p0, int p1, int p2, double p3, double p4, double p5):
        DiscreteParametric(p0, p1, p2, p3, p4, p5), self(self_) {}

    Parametric_Wrapper(PyObject* self_, const Distribution& p0):
        DiscreteParametric(p0), self(self_) {}

    Parametric_Wrapper(PyObject* self_, const Distribution& p0, int p1):
        DiscreteParametric(p0, p1), self(self_) {}

    Parametric_Wrapper(PyObject* self_, const Distribution& p0, double p1):
        DiscreteParametric(p0, p1), self(self_) {}

    Parametric_Wrapper(PyObject* self_, const DiscreteParametric& p0, double p1):
        DiscreteParametric(p0, p1), self(self_) {}

    Parametric_Wrapper(PyObject* self_, const Forward& p0):
        DiscreteParametric(p0), self(self_) {}

    Parametric_Wrapper(PyObject* self_, const DiscreteParametric& p0):
        DiscreteParametric(p0), self(self_) {}

    Parametric_Wrapper(PyObject* self_, const DiscreteParametric& p0, char p1):
        DiscreteParametric(p0, p1), self(self_) {}

    Parametric_Wrapper(PyObject* self_, const DiscreteParametric& p0, char p1, int p2):
        DiscreteParametric(p0, p1, p2), self(self_) {}

    std::basic_ostream<char,std::char_traits<char> >& plot_title_print(std::basic_ostream<char,std::char_traits<char> >& p0) const {
        return call_method< std::basic_ostream<char,std::char_traits<char> >& >(self, "plot_title_print", p0);
    }

    std::basic_ostream<char,std::char_traits<char> >& default_plot_title_print(std::basic_ostream<char,std::char_traits<char> >& p0) const {
        return DiscreteParametric::plot_title_print(p0);
    }

    PyObject* self;
};


BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Parametric_computation_overloads_0_2, computation, 0, 2)

/*************************************************************
 *
 *  Wrappers:
 */



}// namespace

// Module ======================================================================
BOOST_PYTHON_MODULE(csequence)
{
    class_< DiscreteParametric, bases< Distribution > , boost::noncopyable, Parametric_Wrapper >("_DiscreteParametric", init< const Forward& >())
        .def(init< optional< int, int, int, int, double, double > >())
        .def(init< int, int, int, double, double, optional< double > >())
        .def(init< const Distribution&, optional< int > >())
        .def(init< const Distribution&, double >())
        .def(init< const DiscreteParametric&, double >())
        .def(init< const DiscreteParametric&, optional< char, int > >())
        .def("parametric_mean_computation", &DiscreteParametric::parametric_mean_computation)
        .def("parametric_variance_computation", &DiscreteParametric::parametric_variance_computation)
        .def("parametric_skewness_computation", &DiscreteParametric::parametric_skewness_computation)
        .def("parametric_kurtosis_computation", &DiscreteParametric::parametric_kurtosis_computation)
        .def("computation", &DiscreteParametric::computation, Parametric_computation_overloads_0_2())
        .def("simulation", &DiscreteParametric::simulation)
        .def_readonly("ident", &DiscreteParametric::ident)
        .def_readonly("inf_bound", &DiscreteParametric::inf_bound)
        .def_readonly("sup_bound", &DiscreteParametric::sup_bound)
        .def_readonly("parameter", &DiscreteParametric::parameter)
        .def_readonly("probability", &DiscreteParametric::probability)
        .def(self_ns::str(self))
    ;

}
