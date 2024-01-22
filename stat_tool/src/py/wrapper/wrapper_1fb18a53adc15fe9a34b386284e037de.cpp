#include "_stat_tool.h"


void wrapper_1fb18a53adc15fe9a34b386284e037de()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::side_effect > enum_1fb18a53adc15fe9a34b386284e037de("side_effect");
    enum_1fb18a53adc15fe9a34b386284e037de.value("ZERO", ::stat_tool::ZERO);
    enum_1fb18a53adc15fe9a34b386284e037de.value("CONTINUATION", ::stat_tool::CONTINUATION);

}