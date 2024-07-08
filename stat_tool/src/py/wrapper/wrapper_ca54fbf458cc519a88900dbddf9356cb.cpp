#include "_stat_tool.h"


void wrapper_ca54fbf458cc519a88900dbddf9356cb()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::latent_structure_algorithm > enum_ca54fbf458cc519a88900dbddf9356cb("latent_structure_algorithm");
    enum_ca54fbf458cc519a88900dbddf9356cb.value("NO_LATENT_STRUCTURE", ::stat_tool::NO_LATENT_STRUCTURE);
    enum_ca54fbf458cc519a88900dbddf9356cb.value("FORWARD", ::stat_tool::FORWARD);
    enum_ca54fbf458cc519a88900dbddf9356cb.value("VITERBI", ::stat_tool::VITERBI);
    enum_ca54fbf458cc519a88900dbddf9356cb.value("GENERALIZED_VITERBI", ::stat_tool::GENERALIZED_VITERBI);
    enum_ca54fbf458cc519a88900dbddf9356cb.value("FORWARD_BACKWARD_SAMPLING", ::stat_tool::FORWARD_BACKWARD_SAMPLING);
    enum_ca54fbf458cc519a88900dbddf9356cb.value("FORWARD_DYNAMIC_PROGRAMMING", ::stat_tool::FORWARD_DYNAMIC_PROGRAMMING);

}