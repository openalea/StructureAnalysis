#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _stat_tool_variable_type()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::enum_< enum ::stat_tool::variable_type >("variable_type")
            .value("NOMINAL", ::stat_tool::variable_type::NOMINAL)
            .value("ORDINAL", ::stat_tool::variable_type::ORDINAL)
            .value("NUMERIC", ::stat_tool::variable_type::NUMERIC)
            .value("CIRCULAR", ::stat_tool::variable_type::CIRCULAR);
}