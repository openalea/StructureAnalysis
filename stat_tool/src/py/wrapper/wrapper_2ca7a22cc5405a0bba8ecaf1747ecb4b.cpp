#include "_stat_tool.h"


void wrapper_2ca7a22cc5405a0bba8ecaf1747ecb4b()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::linkage > enum_2ca7a22cc5405a0bba8ecaf1747ecb4b("linkage");
    enum_2ca7a22cc5405a0bba8ecaf1747ecb4b.value("NEAREST_NEIGHBOR", ::stat_tool::NEAREST_NEIGHBOR);
    enum_2ca7a22cc5405a0bba8ecaf1747ecb4b.value("FARTHEST_NEIGHBOR", ::stat_tool::FARTHEST_NEIGHBOR);
    enum_2ca7a22cc5405a0bba8ecaf1747ecb4b.value("AVERAGE_NEIGHBOR", ::stat_tool::AVERAGE_NEIGHBOR);

}