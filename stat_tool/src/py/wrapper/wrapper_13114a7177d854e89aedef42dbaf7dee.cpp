#include "_stat_tool.h"


void wrapper_13114a7177d854e89aedef42dbaf7dee()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::censoring_estimator > enum_13114a7177d854e89aedef42dbaf7dee("censoring_estimator");
    enum_13114a7177d854e89aedef42dbaf7dee.value("PARTIAL_LIKELIHOOD", ::stat_tool::PARTIAL_LIKELIHOOD);
    enum_13114a7177d854e89aedef42dbaf7dee.value("COMPLETE_LIKELIHOOD", ::stat_tool::COMPLETE_LIKELIHOOD);
    enum_13114a7177d854e89aedef42dbaf7dee.value("KAPLAN_MEIER", ::stat_tool::KAPLAN_MEIER);

}