#include "_stat_tool.h"


void wrapper_72b962b22d385cfc8c88fecd4d6b9b59()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::angle_unit > enum_72b962b22d385cfc8c88fecd4d6b9b59("angle_unit");
    enum_72b962b22d385cfc8c88fecd4d6b9b59.value("DEGREE", ::stat_tool::DEGREE);
    enum_72b962b22d385cfc8c88fecd4d6b9b59.value("RADIAN", ::stat_tool::RADIAN);

}