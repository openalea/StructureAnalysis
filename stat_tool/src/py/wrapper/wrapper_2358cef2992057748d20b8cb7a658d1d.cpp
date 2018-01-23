#include "_stat_tool.h"


void wrapper_2358cef2992057748d20b8cb7a658d1d()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::log_base > enum_2358cef2992057748d20b8cb7a658d1d("log_base");
    enum_2358cef2992057748d20b8cb7a658d1d.value("NATURAL", ::stat_tool::NATURAL);
    enum_2358cef2992057748d20b8cb7a658d1d.value("TWO", ::stat_tool::TWO);
    enum_2358cef2992057748d20b8cb7a658d1d.value("TEN", ::stat_tool::TEN);

}