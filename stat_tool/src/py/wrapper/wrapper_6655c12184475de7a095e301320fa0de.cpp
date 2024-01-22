#include "_stat_tool.h"


void wrapper_6655c12184475de7a095e301320fa0de()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::variable_nature > enum_6655c12184475de7a095e301320fa0de("variable_nature");
    enum_6655c12184475de7a095e301320fa0de.value("INT_VALUE", ::stat_tool::INT_VALUE);
    enum_6655c12184475de7a095e301320fa0de.value("REAL_VALUE", ::stat_tool::REAL_VALUE);
    enum_6655c12184475de7a095e301320fa0de.value("STATE", ::stat_tool::STATE);
    enum_6655c12184475de7a095e301320fa0de.value("OLD_INT_VALUE", ::stat_tool::OLD_INT_VALUE);
    enum_6655c12184475de7a095e301320fa0de.value("AUXILIARY", ::stat_tool::AUXILIARY);

}