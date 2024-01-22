#include "_stat_tool.h"


void wrapper_eabe4a0bbbd854b3a3634d28c155ddc1()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::compound_distribution > enum_eabe4a0bbbd854b3a3634d28c155ddc1("compound_distribution");
    enum_eabe4a0bbbd854b3a3634d28c155ddc1.value("SUM", ::stat_tool::SUM);
    enum_eabe4a0bbbd854b3a3634d28c155ddc1.value("ELEMENTARY", ::stat_tool::ELEMENTARY);
    enum_eabe4a0bbbd854b3a3634d28c155ddc1.value("COMPOUND", ::stat_tool::COMPOUND);

}