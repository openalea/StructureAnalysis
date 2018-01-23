#include "_stat_tool.h"


void wrapper_ee282d2f51c8552ab1fe6c3d4c138e4b()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::enum_< enum ::std::_Ios_Openmode > enum_ee282d2f51c8552ab1fe6c3d4c138e4b("ios_openmode");
    enum_ee282d2f51c8552ab1fe6c3d4c138e4b.value("S__APP", ::std::_S_app);
    enum_ee282d2f51c8552ab1fe6c3d4c138e4b.value("S__ATE", ::std::_S_ate);
    enum_ee282d2f51c8552ab1fe6c3d4c138e4b.value("S__BIN", ::std::_S_bin);
    enum_ee282d2f51c8552ab1fe6c3d4c138e4b.value("S__IN", ::std::_S_in);
    enum_ee282d2f51c8552ab1fe6c3d4c138e4b.value("S__OUT", ::std::_S_out);
    enum_ee282d2f51c8552ab1fe6c3d4c138e4b.value("S__TRUNC", ::std::_S_trunc);
    enum_ee282d2f51c8552ab1fe6c3d4c138e4b.value("S__IOS_OPENMODE_END", ::std::_S_ios_openmode_end);
    enum_ee282d2f51c8552ab1fe6c3d4c138e4b.value("S__IOS_OPENMODE_MAX", ::std::_S_ios_openmode_max);
    enum_ee282d2f51c8552ab1fe6c3d4c138e4b.value("S__IOS_OPENMODE_MIN", ::std::_S_ios_openmode_min);

}