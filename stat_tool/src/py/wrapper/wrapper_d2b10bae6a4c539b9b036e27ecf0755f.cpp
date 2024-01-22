#include "_stat_tool.h"


void wrapper_d2b10bae6a4c539b9b036e27ecf0755f()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::frequency_distribution_transformation > enum_d2b10bae6a4c539b9b036e27ecf0755f("frequency_distribution_transformation");
    enum_d2b10bae6a4c539b9b036e27ecf0755f.value("FREQUENCY_DISTRIBUTION_COPY", ::stat_tool::FREQUENCY_DISTRIBUTION_COPY);
    enum_d2b10bae6a4c539b9b036e27ecf0755f.value("SHIFT", ::stat_tool::SHIFT);
    enum_d2b10bae6a4c539b9b036e27ecf0755f.value("CLUSTER", ::stat_tool::CLUSTER);

}