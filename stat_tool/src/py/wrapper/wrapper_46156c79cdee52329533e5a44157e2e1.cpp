#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::initializer_list< enum ::stat_tool::discrete_parametric > const volatile * get_pointer<class ::std::initializer_list< enum ::stat_tool::discrete_parametric > const volatile >(class ::std::initializer_list< enum ::stat_tool::discrete_parametric > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_46156c79cdee52329533e5a44157e2e1()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    ::std::initializer_list< enum ::stat_tool::discrete_parametric >::size_type  (::std::initializer_list< enum ::stat_tool::discrete_parametric >::*method_pointer_1fee46c7672250a490e4eafe9fc16b56)() const = &::std::initializer_list< enum ::stat_tool::discrete_parametric >::size;
    boost::python::class_< class ::std::initializer_list< enum ::stat_tool::discrete_parametric >, autowig::Held< class ::std::initializer_list< enum ::stat_tool::discrete_parametric > >::Type > class_46156c79cdee52329533e5a44157e2e1("_InitializerList_46156c79cdee52329533e5a44157e2e1", "", boost::python::no_init);
    class_46156c79cdee52329533e5a44157e2e1.def("__len__", method_pointer_1fee46c7672250a490e4eafe9fc16b56, "");

}