#include "_stat_tool.h"


void wrapper_1fbcb2e114fb5e8489ea1e607e5b5208()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::tying_rule > enum_1fbcb2e114fb5e8489ea1e607e5b5208("tying_rule");
    enum_1fbcb2e114fb5e8489ea1e607e5b5208.value("INDEPENDENT", ::stat_tool::INDEPENDENT);
    enum_1fbcb2e114fb5e8489ea1e607e5b5208.value("CONVOLUTION_FACTOR", ::stat_tool::CONVOLUTION_FACTOR);
    enum_1fbcb2e114fb5e8489ea1e607e5b5208.value("SCALING_FACTOR", ::stat_tool::SCALING_FACTOR);

}