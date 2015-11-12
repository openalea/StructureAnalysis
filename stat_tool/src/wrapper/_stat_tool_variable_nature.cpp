#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _stat_tool_variable_nature()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::enum_< enum ::stat_tool::variable_nature >("variable_nature")
            .value("INT_VALUE", ::stat_tool::variable_nature::INT_VALUE)
            .value("REAL_VALUE", ::stat_tool::variable_nature::REAL_VALUE)
            .value("STATE", ::stat_tool::variable_nature::STATE)
            .value("OLD_INT_VALUE", ::stat_tool::variable_nature::OLD_INT_VALUE)
            .value("AUXILIARY", ::stat_tool::variable_nature::AUXILIARY);
}