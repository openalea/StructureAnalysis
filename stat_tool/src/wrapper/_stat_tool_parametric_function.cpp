#include <boost/python.hpp>
#include <stat_tool/regression.h>

void _stat_tool_parametric_function()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::enum_< enum ::stat_tool::parametric_function >("parametric_function")
            .value("LINEAR_FUNCTION", ::stat_tool::parametric_function::LINEAR_FUNCTION)
            .value("LOGISTIC", ::stat_tool::parametric_function::LOGISTIC)
            .value("MONOMOLECULAR", ::stat_tool::parametric_function::MONOMOLECULAR)
            .value("NONPARAMETRIC_FUNCTION", ::stat_tool::parametric_function::NONPARAMETRIC_FUNCTION)
            .value("CONSTANT_FUNCTION", ::stat_tool::parametric_function::CONSTANT_FUNCTION);
}