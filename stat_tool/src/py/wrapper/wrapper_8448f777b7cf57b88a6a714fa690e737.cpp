#include "_stat_tool.h"


void wrapper_8448f777b7cf57b88a6a714fa690e737()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::process_type > enum_8448f777b7cf57b88a6a714fa690e737("process_type");
    enum_8448f777b7cf57b88a6a714fa690e737.value("ORDINARY", ::stat_tool::ORDINARY);
    enum_8448f777b7cf57b88a6a714fa690e737.value("EQUILIBRIUM", ::stat_tool::EQUILIBRIUM);
    enum_8448f777b7cf57b88a6a714fa690e737.value("DEFAULT_TYPE", ::stat_tool::DEFAULT_TYPE);

}