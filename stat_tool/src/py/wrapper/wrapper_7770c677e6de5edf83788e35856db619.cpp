#include "_stat_tool.h"


void wrapper_7770c677e6de5edf83788e35856db619()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::penalty_type > enum_7770c677e6de5edf83788e35856db619("penalty_type");
    enum_7770c677e6de5edf83788e35856db619.value("FIRST_DIFFERENCE", ::stat_tool::FIRST_DIFFERENCE);
    enum_7770c677e6de5edf83788e35856db619.value("SECOND_DIFFERENCE", ::stat_tool::SECOND_DIFFERENCE);
    enum_7770c677e6de5edf83788e35856db619.value("ENTROPY", ::stat_tool::ENTROPY);

}