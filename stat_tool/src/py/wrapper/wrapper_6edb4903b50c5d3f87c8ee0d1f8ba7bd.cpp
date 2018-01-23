#include "_stat_tool.h"


void wrapper_6edb4903b50c5d3f87c8ee0d1f8ba7bd()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::curve_transformation > enum_6edb4903b50c5d3f87c8ee0d1f8ba7bd("curve_transformation");
    enum_6edb4903b50c5d3f87c8ee0d1f8ba7bd.value("CURVE_COPY", ::stat_tool::CURVE_COPY);
    enum_6edb4903b50c5d3f87c8ee0d1f8ba7bd.value("SMOOTHING", ::stat_tool::SMOOTHING);

}