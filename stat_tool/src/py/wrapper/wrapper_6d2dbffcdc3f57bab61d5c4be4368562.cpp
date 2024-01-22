#include "_stat_tool.h"


void wrapper_6d2dbffcdc3f57bab61d5c4be4368562()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::threshold_direction > enum_6d2dbffcdc3f57bab61d5c4be4368562("threshold_direction");
    enum_6d2dbffcdc3f57bab61d5c4be4368562.value("ABOVE", ::stat_tool::ABOVE);
    enum_6d2dbffcdc3f57bab61d5c4be4368562.value("BELOW", ::stat_tool::BELOW);

}