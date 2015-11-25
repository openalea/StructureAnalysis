#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/discrete_mixture.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>

void _stat_tool_forward()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        void (::stat_tool::Forward::*method_pointer_bfc33d5dcf2453368e8c72f45ac32592)(class ::stat_tool::DiscreteParametric const &) = &::stat_tool::Forward::computation;
        boost::python::class_< class ::stat_tool::Forward, std::shared_ptr< class ::stat_tool::Forward >, boost::python::bases< class ::stat_tool::DiscreteParametric > >("Forward", boost::python::no_init)
            .def(boost::python::init< int, enum ::stat_tool::discrete_parametric, int, int, double, double >())
            .def(boost::python::init< class ::stat_tool::DiscreteParametric const &, int >())
            .def(boost::python::init< class ::stat_tool::Forward const &, int >())
            .def("computation", method_pointer_bfc33d5dcf2453368e8c72f45ac32592);
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::Forward >, std::shared_ptr< class ::stat_tool::DiscreteParametric > >();
}