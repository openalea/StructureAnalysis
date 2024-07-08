#include "_stat_tool.h"


void wrapper_fd2600c41cee533c8c9f0e437b2ea032()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::rounding > enum_fd2600c41cee533c8c9f0e437b2ea032("rounding");
    enum_fd2600c41cee533c8c9f0e437b2ea032.value("FLOOR", ::stat_tool::FLOOR);
    enum_fd2600c41cee533c8c9f0e437b2ea032.value("ROUND", ::stat_tool::ROUND);
    enum_fd2600c41cee533c8c9f0e437b2ea032.value("CEIL", ::stat_tool::CEIL);

}