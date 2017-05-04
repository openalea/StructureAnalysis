#include "_stat_tool.h"


void wrapper_a4ec1929ca0853fbb1225d63f381e9e9()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    int  (*function_pointer_a4ec1929ca0853fbb1225d63f381e9e9)(double , double ) = ::stat_tool::column_width;
    boost::python::def("column_width", function_pointer_a4ec1929ca0853fbb1225d63f381e9e9, "");
}