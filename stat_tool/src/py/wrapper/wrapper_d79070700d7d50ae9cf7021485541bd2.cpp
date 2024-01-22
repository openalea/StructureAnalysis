#include "_stat_tool.h"


void wrapper_d79070700d7d50ae9cf7021485541bd2()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::error_type > enum_d79070700d7d50ae9cf7021485541bd2("error_type");
    enum_d79070700d7d50ae9cf7021485541bd2.value("ERROR", ::stat_tool::ERROR);
    enum_d79070700d7d50ae9cf7021485541bd2.value("WARNING", ::stat_tool::WARNING);

}