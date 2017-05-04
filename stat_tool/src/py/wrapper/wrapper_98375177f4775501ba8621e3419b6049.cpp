#include "_stat_tool.h"


void wrapper_98375177f4775501ba8621e3419b6049()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    int  (*function_pointer_98375177f4775501ba8621e3419b6049)(int ) = ::stat_tool::column_width;
    boost::python::def("column_width", function_pointer_98375177f4775501ba8621e3419b6049, "");
}