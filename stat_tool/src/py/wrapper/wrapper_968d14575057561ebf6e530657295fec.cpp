#include "_stat_tool.h"


void wrapper_968d14575057561ebf6e530657295fec()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::model_type > enum_968d14575057561ebf6e530657295fec("model_type");
    enum_968d14575057561ebf6e530657295fec.value("MIXTURE", ::stat_tool::MIXTURE);
    enum_968d14575057561ebf6e530657295fec.value("HIDDEN_MARKOV", ::stat_tool::HIDDEN_MARKOV);

}