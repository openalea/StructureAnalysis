#include "_stat_tool.h"


void wrapper_96ea627501dc57cb8d2e50bd0c5dbad7()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::distribution_transformation > enum_96ea627501dc57cb8d2e50bd0c5dbad7("distribution_transformation");
    enum_96ea627501dc57cb8d2e50bd0c5dbad7.value("DISTRIBUTION_COPY", ::stat_tool::DISTRIBUTION_COPY);
    enum_96ea627501dc57cb8d2e50bd0c5dbad7.value("NORMALIZATION", ::stat_tool::NORMALIZATION);

}