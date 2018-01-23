#include "_stat_tool.h"


void wrapper_918075f5d686506daddebdc8a4198d0f()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::process_distribution > enum_918075f5d686506daddebdc8a4198d0f("process_distribution");
    enum_918075f5d686506daddebdc8a4198d0f.value("SELF_TRANSITION", ::stat_tool::SELF_TRANSITION);
    enum_918075f5d686506daddebdc8a4198d0f.value("OBSERVATION", ::stat_tool::OBSERVATION);
    enum_918075f5d686506daddebdc8a4198d0f.value("INTENSITY", ::stat_tool::INTENSITY);
    enum_918075f5d686506daddebdc8a4198d0f.value("FIRST_OCCURRENCE", ::stat_tool::FIRST_OCCURRENCE);
    enum_918075f5d686506daddebdc8a4198d0f.value("RECURRENCE_TIME", ::stat_tool::RECURRENCE_TIME);
    enum_918075f5d686506daddebdc8a4198d0f.value("SOJOURN_TIME", ::stat_tool::SOJOURN_TIME);
    enum_918075f5d686506daddebdc8a4198d0f.value("INITIAL_RUN", ::stat_tool::INITIAL_RUN);
    enum_918075f5d686506daddebdc8a4198d0f.value("FINAL_RUN", ::stat_tool::FINAL_RUN);
    enum_918075f5d686506daddebdc8a4198d0f.value("NB_RUN", ::stat_tool::NB_RUN);
    enum_918075f5d686506daddebdc8a4198d0f.value("NB_OCCURRENCE", ::stat_tool::NB_OCCURRENCE);
    enum_918075f5d686506daddebdc8a4198d0f.value("COUNTING", ::stat_tool::COUNTING);
    enum_918075f5d686506daddebdc8a4198d0f.value("LENGTH", ::stat_tool::LENGTH);

}