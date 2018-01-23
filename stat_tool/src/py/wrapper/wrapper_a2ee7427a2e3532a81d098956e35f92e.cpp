#include "_stat_tool.h"


void wrapper_a2ee7427a2e3532a81d098956e35f92e()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::vector_transformation > enum_a2ee7427a2e3532a81d098956e35f92e("vector_transformation");
    enum_a2ee7427a2e3532a81d098956e35f92e.value("VECTOR_COPY", ::stat_tool::VECTOR_COPY);
    enum_a2ee7427a2e3532a81d098956e35f92e.value("ADD_COMPONENT_VARIABLE", ::stat_tool::ADD_COMPONENT_VARIABLE);

}