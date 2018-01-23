#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::initializer_list< enum ::stat_tool::process_distribution > const volatile * get_pointer<class ::std::initializer_list< enum ::stat_tool::process_distribution > const volatile >(class ::std::initializer_list< enum ::stat_tool::process_distribution > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_cc80f44cefab5e0f9c9fff48bf8aa217()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    ::std::initializer_list< enum ::stat_tool::process_distribution >::size_type  (::std::initializer_list< enum ::stat_tool::process_distribution >::*method_pointer_c748d79c46c35a699ee30f1117684b73)() const = &::std::initializer_list< enum ::stat_tool::process_distribution >::size;
    boost::python::class_< class ::std::initializer_list< enum ::stat_tool::process_distribution >, autowig::Held< class ::std::initializer_list< enum ::stat_tool::process_distribution > >::Type > class_cc80f44cefab5e0f9c9fff48bf8aa217("_InitializerList_cc80f44cefab5e0f9c9fff48bf8aa217", "", boost::python::no_init);
    class_cc80f44cefab5e0f9c9fff48bf8aa217.def("__len__", method_pointer_c748d79c46c35a699ee30f1117684b73, "");

}