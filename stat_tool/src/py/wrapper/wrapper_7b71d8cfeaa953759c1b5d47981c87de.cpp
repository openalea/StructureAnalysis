#include "_stat_tool.h"


void wrapper_7b71d8cfeaa953759c1b5d47981c87de()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    double  (*function_pointer_7b71d8cfeaa953759c1b5d47981c87de)(class ::stat_tool::Reestimation< double > *, class ::stat_tool::Reestimation< double > *) = ::stat_tool::interval_bisection;
    boost::python::def("interval_bisection", function_pointer_7b71d8cfeaa953759c1b5d47981c87de, "");
}