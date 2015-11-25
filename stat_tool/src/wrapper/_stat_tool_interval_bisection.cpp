#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/discrete_mixture.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>

void _stat_tool_interval_bisection()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        double (*function_pointer_7b71d8cfeaa953759c1b5d47981c87de)(class ::stat_tool::Reestimation<double> *, class ::stat_tool::Reestimation<double> *) = ::stat_tool::interval_bisection;
        boost::python::def("interval_bisection", function_pointer_7b71d8cfeaa953759c1b5d47981c87de);
}