#include "_stat_tool.h"


void wrapper_f2037461cb4052b7bd89cc9ba087647c()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::duration_distribution_mean_estimator > enum_f2037461cb4052b7bd89cc9ba087647c("duration_distribution_mean_estimator");
    enum_f2037461cb4052b7bd89cc9ba087647c.value("COMPUTED", ::stat_tool::COMPUTED);
    enum_f2037461cb4052b7bd89cc9ba087647c.value("ESTIMATED", ::stat_tool::ESTIMATED);
    enum_f2037461cb4052b7bd89cc9ba087647c.value("ONE_STEP_LATE", ::stat_tool::ONE_STEP_LATE);

}