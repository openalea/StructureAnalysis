#include "_stat_tool.h"


void wrapper_a87ca98a1bfb5734926ff985d38c7479()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::isolation_scale > enum_a87ca98a1bfb5734926ff985d38c7479("isolation_scale");
    enum_a87ca98a1bfb5734926ff985d38c7479.value("INDIVIDUAL", ::stat_tool::INDIVIDUAL);
    enum_a87ca98a1bfb5734926ff985d38c7479.value("CLUSTER_SCALE", ::stat_tool::CLUSTER_SCALE);

}