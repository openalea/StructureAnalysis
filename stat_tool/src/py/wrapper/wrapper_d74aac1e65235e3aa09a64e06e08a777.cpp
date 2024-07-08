#include "_stat_tool.h"


void wrapper_d74aac1e65235e3aa09a64e06e08a777()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::correlation_type > enum_d74aac1e65235e3aa09a64e06e08a777("correlation_type");
    enum_d74aac1e65235e3aa09a64e06e08a777.value("PEARSON", ::stat_tool::PEARSON);
    enum_d74aac1e65235e3aa09a64e06e08a777.value("SPEARMAN", ::stat_tool::SPEARMAN);
    enum_d74aac1e65235e3aa09a64e06e08a777.value("KENDALL", ::stat_tool::KENDALL);
    enum_d74aac1e65235e3aa09a64e06e08a777.value("SPEARMAN2", ::stat_tool::SPEARMAN2);

}