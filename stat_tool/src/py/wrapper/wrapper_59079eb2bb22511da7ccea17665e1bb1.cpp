#include "_stat_tool.h"


void wrapper_59079eb2bb22511da7ccea17665e1bb1()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    double  (*function_pointer_59079eb2bb22511da7ccea17665e1bb1)(double ) = ::stat_tool::von_mises_concentration_computation;
    boost::python::def("von_mises_concentration_computation", function_pointer_59079eb2bb22511da7ccea17665e1bb1, "");
}