#include "_stat_tool.h"


void wrapper_ee6a7db98b4a5be8b6b0ce867301d18c()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::variable_type > enum_ee6a7db98b4a5be8b6b0ce867301d18c("variable_type");
    enum_ee6a7db98b4a5be8b6b0ce867301d18c.value("NOMINAL", ::stat_tool::NOMINAL);
    enum_ee6a7db98b4a5be8b6b0ce867301d18c.value("ORDINAL", ::stat_tool::ORDINAL);
    enum_ee6a7db98b4a5be8b6b0ce867301d18c.value("NUMERIC", ::stat_tool::NUMERIC);
    enum_ee6a7db98b4a5be8b6b0ce867301d18c.value("CIRCULAR", ::stat_tool::CIRCULAR);

}