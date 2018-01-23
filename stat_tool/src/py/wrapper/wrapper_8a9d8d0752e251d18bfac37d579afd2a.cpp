#include "_stat_tool.h"


void wrapper_8a9d8d0752e251d18bfac37d579afd2a()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::test_distribution > enum_8a9d8d0752e251d18bfac37d579afd2a("test_distribution");
    enum_8a9d8d0752e251d18bfac37d579afd2a.value("STANDARD_NORMAL", ::stat_tool::STANDARD_NORMAL);
    enum_8a9d8d0752e251d18bfac37d579afd2a.value("CHI2", ::stat_tool::CHI2);
    enum_8a9d8d0752e251d18bfac37d579afd2a.value("FISHER", ::stat_tool::FISHER);
    enum_8a9d8d0752e251d18bfac37d579afd2a.value("STUDENT", ::stat_tool::STUDENT);

}