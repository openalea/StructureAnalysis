#include "_stat_tool.h"


void wrapper_bbcc92c9b55f53388df45016b45c90de()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::continuous_parametric > enum_bbcc92c9b55f53388df45016b45c90de("continuous_parametric");
    enum_bbcc92c9b55f53388df45016b45c90de.value("GAMMA", ::stat_tool::GAMMA);
    enum_bbcc92c9b55f53388df45016b45c90de.value("INVERSE_GAUSSIAN", ::stat_tool::INVERSE_GAUSSIAN);
    enum_bbcc92c9b55f53388df45016b45c90de.value("GAUSSIAN", ::stat_tool::GAUSSIAN);
    enum_bbcc92c9b55f53388df45016b45c90de.value("VON_MISES", ::stat_tool::VON_MISES);
    enum_bbcc92c9b55f53388df45016b45c90de.value("ZERO_INFLATED_GAMMA", ::stat_tool::ZERO_INFLATED_GAMMA);
    enum_bbcc92c9b55f53388df45016b45c90de.value("LINEAR_MODEL", ::stat_tool::LINEAR_MODEL);
    enum_bbcc92c9b55f53388df45016b45c90de.value("AUTOREGRESSIVE_MODEL", ::stat_tool::AUTOREGRESSIVE_MODEL);

}