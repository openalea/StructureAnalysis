#include "_stat_tool.h"


void wrapper_487fbf2abc7d54409b96103aaef60d8d()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::count_pattern > enum_487fbf2abc7d54409b96103aaef60d8d("count_pattern");
    enum_487fbf2abc7d54409b96103aaef60d8d.value("RUN", ::stat_tool::RUN);
    enum_487fbf2abc7d54409b96103aaef60d8d.value("OCCURRENCE", ::stat_tool::OCCURRENCE);

}