#include <boost/python.hpp>
#include <stat_tool/markovian.h>

void _stat_tool_latent_structure_algorithm()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::enum_< enum ::stat_tool::latent_structure_algorithm >("latent_structure_algorithm")
            .value("NO_LATENT_STRUCTURE", ::stat_tool::latent_structure_algorithm::NO_LATENT_STRUCTURE)
            .value("FORWARD", ::stat_tool::latent_structure_algorithm::FORWARD)
            .value("VITERBI", ::stat_tool::latent_structure_algorithm::VITERBI)
            .value("GENERALIZED_VITERBI", ::stat_tool::latent_structure_algorithm::GENERALIZED_VITERBI)
            .value("FORWARD_BACKWARD_SAMPLING", ::stat_tool::latent_structure_algorithm::FORWARD_BACKWARD_SAMPLING)
            .value("FORWARD_DYNAMIC_PROGRAMMING", ::stat_tool::latent_structure_algorithm::FORWARD_DYNAMIC_PROGRAMMING);
}