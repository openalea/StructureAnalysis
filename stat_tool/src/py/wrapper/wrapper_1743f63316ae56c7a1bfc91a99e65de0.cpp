#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::error_condition const volatile * get_pointer<class ::std::error_condition const volatile >(class ::std::error_condition const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_1743f63316ae56c7a1bfc91a99e65de0()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    void  (::std::error_condition::*method_pointer_a8d37d62bd7d51128221f95783655119)() = &::std::error_condition::clear;
    int  (::std::error_condition::*method_pointer_d011ea043494523cbd290359b1df8d8e)() const = &::std::error_condition::value;
    ::std::string  (::std::error_condition::*method_pointer_8854cd90233c5637bb130aa921eedcf3)() const = &::std::error_condition::message;
    boost::python::class_< class ::std::error_condition, autowig::Held< class ::std::error_condition >::Type > class_1743f63316ae56c7a1bfc91a99e65de0("ErrorCondition", "", boost::python::no_init);
    class_1743f63316ae56c7a1bfc91a99e65de0.def(boost::python::init<  >(""));
    class_1743f63316ae56c7a1bfc91a99e65de0.def(boost::python::init< class ::std::error_condition const & >(""));
    class_1743f63316ae56c7a1bfc91a99e65de0.def("clear", method_pointer_a8d37d62bd7d51128221f95783655119, "");
    class_1743f63316ae56c7a1bfc91a99e65de0.def("value", method_pointer_d011ea043494523cbd290359b1df8d8e, "");
    class_1743f63316ae56c7a1bfc91a99e65de0.def("message", method_pointer_8854cd90233c5637bb130aa921eedcf3, "");

}