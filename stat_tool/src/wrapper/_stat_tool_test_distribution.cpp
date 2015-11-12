#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _stat_tool_test_distribution()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::enum_< enum ::stat_tool::test_distribution >("test_distribution")
            .value("STANDARD_NORMAL", ::stat_tool::test_distribution::STANDARD_NORMAL)
            .value("CHI2", ::stat_tool::test_distribution::CHI2)
            .value("FISHER", ::stat_tool::test_distribution::FISHER)
            .value("STUDENT", ::stat_tool::test_distribution::STUDENT);
}