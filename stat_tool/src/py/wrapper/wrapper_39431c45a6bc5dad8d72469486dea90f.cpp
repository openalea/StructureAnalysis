#include "_stat_tool.h"


void wrapper_39431c45a6bc5dad8d72469486dea90f()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::cluster_scale > enum_39431c45a6bc5dad8d72469486dea90f("cluster_scale");
    enum_39431c45a6bc5dad8d72469486dea90f.value("CHILD_CLUSTER_DISTANCE", ::stat_tool::CHILD_CLUSTER_DISTANCE);
    enum_39431c45a6bc5dad8d72469486dea90f.value("DIAMETER", ::stat_tool::DIAMETER);

}