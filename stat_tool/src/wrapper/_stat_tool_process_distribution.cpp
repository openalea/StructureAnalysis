#include <boost/python.hpp>
#include <stat_tool/plotable.h>

void _stat_tool_process_distribution()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::enum_< enum ::stat_tool::process_distribution >("process_distribution")
            .value("SELF_TRANSITION", ::stat_tool::process_distribution::SELF_TRANSITION)
            .value("OBSERVATION", ::stat_tool::process_distribution::OBSERVATION)
            .value("INTENSITY", ::stat_tool::process_distribution::INTENSITY)
            .value("FIRST_OCCURRENCE", ::stat_tool::process_distribution::FIRST_OCCURRENCE)
            .value("RECURRENCE_TIME", ::stat_tool::process_distribution::RECURRENCE_TIME)
            .value("SOJOURN_TIME", ::stat_tool::process_distribution::SOJOURN_TIME)
            .value("INITIAL_RUN", ::stat_tool::process_distribution::INITIAL_RUN)
            .value("FINAL_RUN", ::stat_tool::process_distribution::FINAL_RUN)
            .value("NB_RUN", ::stat_tool::process_distribution::NB_RUN)
            .value("NB_OCCURRENCE", ::stat_tool::process_distribution::NB_OCCURRENCE)
            .value("COUNTING", ::stat_tool::process_distribution::COUNTING)
            .value("LENGTH", ::stat_tool::process_distribution::LENGTH);
}