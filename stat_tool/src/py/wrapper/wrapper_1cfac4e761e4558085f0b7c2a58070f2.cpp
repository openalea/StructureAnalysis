#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> struct ::std::error_code const volatile * get_pointer<struct ::std::error_code const volatile >(struct ::std::error_code const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_1cfac4e761e4558085f0b7c2a58070f2()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    void  (::std::error_code::*method_pointer_983f6b690936572f9b22589850b87def)() = &::std::error_code::clear;
    int  (::std::error_code::*method_pointer_15ef1dc09b265385bc5d0bc13c03cd5f)() const = &::std::error_code::value;
    struct ::std::error_condition  (::std::error_code::*method_pointer_aa96e684628c562dbdbb67930ecd700f)() const = &::std::error_code::default_error_condition;
    ::std::string  (::std::error_code::*method_pointer_3d779aea45e7530a97991c2ae194a791)() const = &::std::error_code::message;
    boost::python::class_< struct ::std::error_code, autowig::Held< struct ::std::error_code >::Type > class_1cfac4e761e4558085f0b7c2a58070f2("ErrorCode", "", boost::python::no_init);
    class_1cfac4e761e4558085f0b7c2a58070f2.def(boost::python::init<  >(""));
    class_1cfac4e761e4558085f0b7c2a58070f2.def(boost::python::init< struct ::std::error_code const & >(""));
    class_1cfac4e761e4558085f0b7c2a58070f2.def("clear", method_pointer_983f6b690936572f9b22589850b87def, "");
    class_1cfac4e761e4558085f0b7c2a58070f2.def("value", method_pointer_15ef1dc09b265385bc5d0bc13c03cd5f, "");
    class_1cfac4e761e4558085f0b7c2a58070f2.def("default_error_condition", method_pointer_aa96e684628c562dbdbb67930ecd700f, "");
    class_1cfac4e761e4558085f0b7c2a58070f2.def("message", method_pointer_3d779aea45e7530a97991c2ae194a791, "");

}