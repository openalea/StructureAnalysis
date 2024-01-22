#include "_stat_tool.h"


void wrapper_8590aefdf6b151df97896fc931a26e54()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::metric > enum_8590aefdf6b151df97896fc931a26e54("metric");
    enum_8590aefdf6b151df97896fc931a26e54.value("ABSOLUTE_VALUE", ::stat_tool::ABSOLUTE_VALUE);
    enum_8590aefdf6b151df97896fc931a26e54.value("QUADRATIC", ::stat_tool::QUADRATIC);

}