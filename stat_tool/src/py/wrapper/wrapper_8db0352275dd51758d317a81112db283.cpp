#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::error_code const volatile * get_pointer<class ::std::error_code const volatile >(class ::std::error_code const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_8db0352275dd51758d317a81112db283()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    void  (::std::error_code::*method_pointer_035dcb41157658e1ad48fd44f71ea579)() = &::std::error_code::clear;
    int  (::std::error_code::*method_pointer_b85dbc38f6d753868135a87acd658694)() const = &::std::error_code::value;
    class ::std::error_condition  (::std::error_code::*method_pointer_277754b735a054bd9e2836ce3cd4e26f)() const = &::std::error_code::default_error_condition;
    ::std::string  (::std::error_code::*method_pointer_017cf1f33e7558b2b4e76b618ea95ba8)() const = &::std::error_code::message;
    boost::python::class_< class ::std::error_code, autowig::Held< class ::std::error_code >::Type > class_8db0352275dd51758d317a81112db283("ErrorCode", "", boost::python::no_init);
    class_8db0352275dd51758d317a81112db283.def(boost::python::init<  >(""));
    class_8db0352275dd51758d317a81112db283.def(boost::python::init< class ::std::error_code const & >(""));
    class_8db0352275dd51758d317a81112db283.def("clear", method_pointer_035dcb41157658e1ad48fd44f71ea579, "");
    class_8db0352275dd51758d317a81112db283.def("value", method_pointer_b85dbc38f6d753868135a87acd658694, "");
    class_8db0352275dd51758d317a81112db283.def("default_error_condition", method_pointer_277754b735a054bd9e2836ce3cd4e26f, "");
    class_8db0352275dd51758d317a81112db283.def("message", method_pointer_017cf1f33e7558b2b4e76b618ea95ba8, "");

}