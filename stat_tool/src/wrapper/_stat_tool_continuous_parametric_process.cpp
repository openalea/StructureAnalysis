#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>
#include <stat_tool/discrete_mixture.h>
#include <stat_tool/markovian.h>

void _stat_tool_continuous_parametric_process()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        void (::stat_tool::ContinuousParametricProcess::*method_pointer_a51fa7e43cd559d48e6193d1985ac289)(class ::stat_tool::ContinuousParametricProcess const &) = &::stat_tool::ContinuousParametricProcess::copy;
        void (::stat_tool::ContinuousParametricProcess::*method_pointer_c0ed3022bc4c5e1ba45f2928b3b5af9b)() = &::stat_tool::ContinuousParametricProcess::remove;
        int (::stat_tool::ContinuousParametricProcess::*method_pointer_d474df83ba195ee181fd24c1a5b56fc2)() const = &::stat_tool::ContinuousParametricProcess::nb_parameter_computation;
        double (::stat_tool::ContinuousParametricProcess::*method_pointer_0a554efc0f6c5d83a995f8fbf54cad2d)(class ::stat_tool::Distribution *) const = &::stat_tool::ContinuousParametricProcess::mean_computation;
        double (::stat_tool::ContinuousParametricProcess::*method_pointer_cebc3dc50e265dddb861748541b68e5b)(class ::stat_tool::Distribution *, double) const = &::stat_tool::ContinuousParametricProcess::variance_computation;
        void (::stat_tool::ContinuousParametricProcess::*method_pointer_e4038a834f085556ba7622988e069150)(enum ::stat_tool::angle_unit) = &::stat_tool::ContinuousParametricProcess::select_unit;
        void (::stat_tool::ContinuousParametricProcess::*method_pointer_f1e34b3cae8354bf88b0b69db8db7c67)(enum ::stat_tool::continuous_parametric, double, double, double, double) = &::stat_tool::ContinuousParametricProcess::init;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::ContinuousParametricProcess::*method_pointer_ca47b5e35f5153f0ae9795874e4324e4)(class ::std::basic_ostream<char, std::char_traits<char> > &) = &::stat_tool::ContinuousParametricProcess::interval_computation;
        boost::python::class_< class ::stat_tool::ContinuousParametricProcess, std::shared_ptr< class ::stat_tool::ContinuousParametricProcess > >("ContinuousParametricProcess", boost::python::no_init)
            .def(boost::python::init< int >())
            .def(boost::python::init< class ::stat_tool::ContinuousParametricProcess const & >())
            .def("copy", method_pointer_a51fa7e43cd559d48e6193d1985ac289)
            .def("remove", method_pointer_c0ed3022bc4c5e1ba45f2928b3b5af9b)
            .def("nb_parameter_computation", method_pointer_d474df83ba195ee181fd24c1a5b56fc2)
            .def("mean_computation", method_pointer_0a554efc0f6c5d83a995f8fbf54cad2d)
            .def("variance_computation", method_pointer_cebc3dc50e265dddb861748541b68e5b)
            .def("select_unit", method_pointer_e4038a834f085556ba7622988e069150)
            .def("init", method_pointer_f1e34b3cae8354bf88b0b69db8db7c67)
            .def("interval_computation", method_pointer_ca47b5e35f5153f0ae9795874e4324e4, boost::python::return_internal_reference<>())
            .def_readwrite("nb_state", &::stat_tool::ContinuousParametricProcess::nb_state)
            .def_readwrite("ident", &::stat_tool::ContinuousParametricProcess::ident)
            .def_readwrite("tied_location", &::stat_tool::ContinuousParametricProcess::tied_location)
            .def_readwrite("tied_dispersion", &::stat_tool::ContinuousParametricProcess::tied_dispersion)
            .def_readwrite("unit", &::stat_tool::ContinuousParametricProcess::unit)
            .def_readwrite("weight", &::stat_tool::ContinuousParametricProcess::weight)
            .def_readwrite("restoration_weight", &::stat_tool::ContinuousParametricProcess::restoration_weight);
}