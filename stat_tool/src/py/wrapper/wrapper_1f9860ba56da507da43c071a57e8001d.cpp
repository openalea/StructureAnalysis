#include "_stat_tool.h"


void wrapper_1f9860ba56da507da43c071a57e8001d()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::moving_average_method > enum_1f9860ba56da507da43c071a57e8001d("moving_average_method");
    enum_1f9860ba56da507da43c071a57e8001d.value("AVERAGING", ::stat_tool::AVERAGING);
    enum_1f9860ba56da507da43c071a57e8001d.value("LEAST_SQUARES", ::stat_tool::LEAST_SQUARES);

}