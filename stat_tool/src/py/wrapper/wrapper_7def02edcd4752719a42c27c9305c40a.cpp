#include "_stat_tool.h"


void wrapper_7def02edcd4752719a42c27c9305c40a()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::parametric_function > enum_7def02edcd4752719a42c27c9305c40a("parametric_function");
    enum_7def02edcd4752719a42c27c9305c40a.value("LINEAR_FUNCTION", ::stat_tool::LINEAR_FUNCTION);
    enum_7def02edcd4752719a42c27c9305c40a.value("LOGISTIC", ::stat_tool::LOGISTIC);
    enum_7def02edcd4752719a42c27c9305c40a.value("MONOMOLECULAR", ::stat_tool::MONOMOLECULAR);
    enum_7def02edcd4752719a42c27c9305c40a.value("NONPARAMETRIC_FUNCTION", ::stat_tool::NONPARAMETRIC_FUNCTION);
    enum_7def02edcd4752719a42c27c9305c40a.value("CONSTANT_FUNCTION", ::stat_tool::CONSTANT_FUNCTION);

}