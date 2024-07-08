#include "_stat_tool.h"


void wrapper_fe202c1f62ab5724b611e7a0685a88c5()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::estimation_criterion > enum_fe202c1f62ab5724b611e7a0685a88c5("estimation_criterion");
    enum_fe202c1f62ab5724b611e7a0685a88c5.value("LIKELIHOOD", ::stat_tool::LIKELIHOOD);
    enum_fe202c1f62ab5724b611e7a0685a88c5.value("PENALIZED_LIKELIHOOD", ::stat_tool::PENALIZED_LIKELIHOOD);
    enum_fe202c1f62ab5724b611e7a0685a88c5.value("PARAMETRIC_REGULARIZATION", ::stat_tool::PARAMETRIC_REGULARIZATION);

}