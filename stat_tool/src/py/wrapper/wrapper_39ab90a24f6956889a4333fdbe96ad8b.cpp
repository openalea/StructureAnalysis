#include "_stat_tool.h"


void wrapper_39ab90a24f6956889a4333fdbe96ad8b()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::state_type > enum_39ab90a24f6956889a4333fdbe96ad8b("state_type");
    enum_39ab90a24f6956889a4333fdbe96ad8b.value("TRANSIENT", ::stat_tool::TRANSIENT);
    enum_39ab90a24f6956889a4333fdbe96ad8b.value("RECURRENT", ::stat_tool::RECURRENT);
    enum_39ab90a24f6956889a4333fdbe96ad8b.value("ABSORBING", ::stat_tool::ABSORBING);

}