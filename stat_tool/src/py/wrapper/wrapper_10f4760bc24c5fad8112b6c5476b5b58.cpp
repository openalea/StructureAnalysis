#include "_stat_tool.h"


void wrapper_10f4760bc24c5fad8112b6c5476b5b58()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::distribution_computation > enum_10f4760bc24c5fad8112b6c5476b5b58("distribution_computation");
    enum_10f4760bc24c5fad8112b6c5476b5b58.value("STANDARD", ::stat_tool::STANDARD);
    enum_10f4760bc24c5fad8112b6c5476b5b58.value("RENEWAL", ::stat_tool::RENEWAL);

}