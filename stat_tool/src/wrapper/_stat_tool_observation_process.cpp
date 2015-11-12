#include <boost/python.hpp>
#include <stat_tool/markovian.h>

void _stat_tool_observation_process()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::enum_< enum ::stat_tool::observation_process >("observation_process")
            .value("CATEGORICAL_PROCESS", ::stat_tool::observation_process::CATEGORICAL_PROCESS)
            .value("DISCRETE_PARAMETRIC", ::stat_tool::observation_process::DISCRETE_PARAMETRIC)
            .value("CONTINUOUS_PARAMETRIC", ::stat_tool::observation_process::CONTINUOUS_PARAMETRIC)
            .value("DEFAULT_PROCESS", ::stat_tool::observation_process::DEFAULT_PROCESS);
}