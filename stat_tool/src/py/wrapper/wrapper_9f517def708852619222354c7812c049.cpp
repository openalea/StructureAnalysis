#include "_stat_tool.h"


void wrapper_9f517def708852619222354c7812c049()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::hierarchical_strategy > enum_9f517def708852619222354c7812c049("hierarchical_strategy");
    enum_9f517def708852619222354c7812c049.value("AGGLOMERATIVE", ::stat_tool::AGGLOMERATIVE);
    enum_9f517def708852619222354c7812c049.value("DIVISIVE", ::stat_tool::DIVISIVE);
    enum_9f517def708852619222354c7812c049.value("ORDERING", ::stat_tool::ORDERING);

}