#include "_stat_tool.h"


void wrapper_a18fbbaaa69b5c948c8a228f2d6f4a61()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::observation_process > enum_a18fbbaaa69b5c948c8a228f2d6f4a61("observation_process");
    enum_a18fbbaaa69b5c948c8a228f2d6f4a61.value("CATEGORICAL_PROCESS", ::stat_tool::CATEGORICAL_PROCESS);
    enum_a18fbbaaa69b5c948c8a228f2d6f4a61.value("DISCRETE_PARAMETRIC", ::stat_tool::DISCRETE_PARAMETRIC);
    enum_a18fbbaaa69b5c948c8a228f2d6f4a61.value("CONTINUOUS_PARAMETRIC", ::stat_tool::CONTINUOUS_PARAMETRIC);
    enum_a18fbbaaa69b5c948c8a228f2d6f4a61.value("DEFAULT_PROCESS", ::stat_tool::DEFAULT_PROCESS);

}