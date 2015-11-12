#include <boost/python.hpp>
#include <stat_tool/markovian.h>

void _stat_tool_von_mises_concentration_computation()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        double (*function_pointer_59079eb2bb22511da7ccea17665e1bb1)(double) = ::stat_tool::von_mises_concentration_computation;
        boost::python::def("von_mises_concentration_computation", function_pointer_59079eb2bb22511da7ccea17665e1bb1);
}