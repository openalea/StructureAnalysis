#include <boost/python.hpp>
#include <stat_tool/stat_label.h>

void _stat_tool_stat_tool_parsing()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::enum_< enum ::stat_tool::stat_tool_parsing >("stat_tool_parsing")
            .value("STATP_KEY_WORD", ::stat_tool::stat_tool_parsing::STATP_KEY_WORD)
            .value("STATP_FORMAT", ::stat_tool::stat_tool_parsing::STATP_FORMAT)
            .value("STATP_DISTRIBUTION_NAME", ::stat_tool::stat_tool_parsing::STATP_DISTRIBUTION_NAME)
            .value("STATP_PARAMETER_NAME", ::stat_tool::stat_tool_parsing::STATP_PARAMETER_NAME)
            .value("STATP_PARAMETER_INDEX", ::stat_tool::stat_tool_parsing::STATP_PARAMETER_INDEX)
            .value("STATP_NB_PARAMETER", ::stat_tool::stat_tool_parsing::STATP_NB_PARAMETER)
            .value("STATP_SEPARATOR", ::stat_tool::stat_tool_parsing::STATP_SEPARATOR)
            .value("STATP_PARAMETER_VALUE", ::stat_tool::stat_tool_parsing::STATP_PARAMETER_VALUE)
            .value("STATP_PROBABILITY_SUM", ::stat_tool::stat_tool_parsing::STATP_PROBABILITY_SUM)
            .value("STATP_DATA_TYPE", ::stat_tool::stat_tool_parsing::STATP_DATA_TYPE)
            .value("STATP_NB_TOKEN", ::stat_tool::stat_tool_parsing::STATP_NB_TOKEN)
            .value("STATP_EMPTY_SAMPLE", ::stat_tool::stat_tool_parsing::STATP_EMPTY_SAMPLE)
            .value("STATP_VALUE_ORDER", ::stat_tool::stat_tool_parsing::STATP_VALUE_ORDER)
            .value("STATP_MAX_VALUE", ::stat_tool::stat_tool_parsing::STATP_MAX_VALUE)
            .value("STATP_NB_DISTRIBUTION", ::stat_tool::stat_tool_parsing::STATP_NB_DISTRIBUTION)
            .value("STATP_DISTRIBUTION_INDEX", ::stat_tool::stat_tool_parsing::STATP_DISTRIBUTION_INDEX)
            .value("STATP_WEIGHT_VALUE", ::stat_tool::stat_tool_parsing::STATP_WEIGHT_VALUE)
            .value("STATP_NB_COMPONENT", ::stat_tool::stat_tool_parsing::STATP_NB_COMPONENT)
            .value("STATP_COMPONENT_INDEX", ::stat_tool::stat_tool_parsing::STATP_COMPONENT_INDEX)
            .value("STATP_NB_VARIABLE", ::stat_tool::stat_tool_parsing::STATP_NB_VARIABLE)
            .value("STATP_VARIABLE_INDEX", ::stat_tool::stat_tool_parsing::STATP_VARIABLE_INDEX)
            .value("STATP_VARIABLE_TYPE", ::stat_tool::stat_tool_parsing::STATP_VARIABLE_TYPE)
            .value("STATP_NB_CATEGORY", ::stat_tool::stat_tool_parsing::STATP_NB_CATEGORY)
            .value("STATP_LOCAL_DISTANCE", ::stat_tool::stat_tool_parsing::STATP_LOCAL_DISTANCE)
            .value("STATP_TRIANGLE_INEQUALITY", ::stat_tool::stat_tool_parsing::STATP_TRIANGLE_INEQUALITY)
            .value("STATP_PERIOD_VALUE", ::stat_tool::stat_tool_parsing::STATP_PERIOD_VALUE)
            .value("STATP_NB_STATE", ::stat_tool::stat_tool_parsing::STATP_NB_STATE)
            .value("STATP_ORDER", ::stat_tool::stat_tool_parsing::STATP_ORDER)
            .value("STATP_INITIAL_PROBABILITY", ::stat_tool::stat_tool_parsing::STATP_INITIAL_PROBABILITY)
            .value("STATP_TRANSITION_PROBABILITY", ::stat_tool::stat_tool_parsing::STATP_TRANSITION_PROBABILITY)
            .value("STATP_CHAIN_STRUCTURE", ::stat_tool::stat_tool_parsing::STATP_CHAIN_STRUCTURE)
            .value("STATP_IRREDUCIBLE", ::stat_tool::stat_tool_parsing::STATP_IRREDUCIBLE)
            .value("STATP_STATE_INDEX", ::stat_tool::stat_tool_parsing::STATP_STATE_INDEX)
            .value("STATP_OUTPUT_INDEX", ::stat_tool::stat_tool_parsing::STATP_OUTPUT_INDEX)
            .value("STATP_OBSERVATION_PROBABILITY", ::stat_tool::stat_tool_parsing::STATP_OBSERVATION_PROBABILITY)
            .value("STATP_NB_OBSERVATION_DISTRIBUTION", ::stat_tool::stat_tool_parsing::STATP_NB_OBSERVATION_DISTRIBUTION)
            .value("STATP_NON_CONSECUTIVE_OUTPUTS", ::stat_tool::stat_tool_parsing::STATP_NON_CONSECUTIVE_OUTPUTS)
            .value("STATP_OBSERVATION_DISTRIBUTION_OVERLAP", ::stat_tool::stat_tool_parsing::STATP_OBSERVATION_DISTRIBUTION_OVERLAP)
            .value("STATP_OBSERVATION_DISTRIBUTION_NON_OVERLAP", ::stat_tool::stat_tool_parsing::STATP_OBSERVATION_DISTRIBUTION_NON_OVERLAP)
            .value("STATP_OBSERVATION_DISTRIBUTION_TYPE", ::stat_tool::stat_tool_parsing::STATP_OBSERVATION_DISTRIBUTION_TYPE)
            .value("STATP_NB_OUTPUT_PROCESS", ::stat_tool::stat_tool_parsing::STATP_NB_OUTPUT_PROCESS)
            .value("STATP_OUTPUT_PROCESS_INDEX", ::stat_tool::stat_tool_parsing::STATP_OUTPUT_PROCESS_INDEX);
}