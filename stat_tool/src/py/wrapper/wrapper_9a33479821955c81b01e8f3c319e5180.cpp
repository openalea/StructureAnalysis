#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> struct ::std::error_condition const volatile * get_pointer<struct ::std::error_condition const volatile >(struct ::std::error_condition const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_9a33479821955c81b01e8f3c319e5180()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    void  (::std::error_condition::*method_pointer_fba0ff58c5425aa799ffc6ee30f6d2c7)() = &::std::error_condition::clear;
    int  (::std::error_condition::*method_pointer_46f6bb84da0e5c45aed5f8ef3e6f7728)() const = &::std::error_condition::value;
    ::std::string  (::std::error_condition::*method_pointer_d5ecd8c37ae55c2d8c2e8e6f83d6c460)() const = &::std::error_condition::message;
    boost::python::class_< struct ::std::error_condition, autowig::Held< struct ::std::error_condition >::Type > class_9a33479821955c81b01e8f3c319e5180("ErrorCondition", "", boost::python::no_init);
    class_9a33479821955c81b01e8f3c319e5180.def(boost::python::init<  >(""));
    class_9a33479821955c81b01e8f3c319e5180.def(boost::python::init< struct ::std::error_condition const & >(""));
    class_9a33479821955c81b01e8f3c319e5180.def("clear", method_pointer_fba0ff58c5425aa799ffc6ee30f6d2c7, "");
    class_9a33479821955c81b01e8f3c319e5180.def("value", method_pointer_46f6bb84da0e5c45aed5f8ef3e6f7728, "");
    class_9a33479821955c81b01e8f3c319e5180.def("message", method_pointer_d5ecd8c37ae55c2d8c2e8e6f83d6c460, "");

}