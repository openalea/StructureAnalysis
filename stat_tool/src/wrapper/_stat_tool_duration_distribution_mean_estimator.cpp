#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _stat_tool_duration_distribution_mean_estimator()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::enum_< enum ::stat_tool::duration_distribution_mean_estimator >("duration_distribution_mean_estimator")
            .value("COMPUTED", ::stat_tool::duration_distribution_mean_estimator::COMPUTED)
            .value("ESTIMATED", ::stat_tool::duration_distribution_mean_estimator::ESTIMATED)
            .value("ONE_STEP_LATE", ::stat_tool::duration_distribution_mean_estimator::ONE_STEP_LATE);
}