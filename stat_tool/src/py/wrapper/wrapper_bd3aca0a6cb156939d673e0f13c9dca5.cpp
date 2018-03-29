#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::initializer_list< class ::stat_tool::DiscreteParametric > const volatile * get_pointer<class ::std::initializer_list< class ::stat_tool::DiscreteParametric > const volatile >(class ::std::initializer_list< class ::stat_tool::DiscreteParametric > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_bd3aca0a6cb156939d673e0f13c9dca5()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    ::std::initializer_list< class ::stat_tool::DiscreteParametric >::size_type  (::std::initializer_list< ::stat_tool::DiscreteParametric >::*method_pointer_6bfb0444a5295ae9b96748554a27e3b8)() const = &::std::initializer_list< class ::stat_tool::DiscreteParametric >::size;
    ::std::initializer_list< class ::stat_tool::DiscreteParametric >::const_iterator  (::std::initializer_list< ::stat_tool::DiscreteParametric >::*method_pointer_7c94877159bf5284a0356addd6b3b05c)() const = &::std::initializer_list< class ::stat_tool::DiscreteParametric >::begin;
    ::std::initializer_list< class ::stat_tool::DiscreteParametric >::const_iterator  (::std::initializer_list< ::stat_tool::DiscreteParametric >::*method_pointer_8a840a970901512fa23c5ee45013b5cb)() const = &::std::initializer_list< class ::stat_tool::DiscreteParametric >::end;
    boost::python::class_< class ::std::initializer_list< class ::stat_tool::DiscreteParametric >, autowig::Held< class ::std::initializer_list< class ::stat_tool::DiscreteParametric > >::Type > class_bd3aca0a6cb156939d673e0f13c9dca5("_InitializerList_bd3aca0a6cb156939d673e0f13c9dca5", "", boost::python::no_init);
    class_bd3aca0a6cb156939d673e0f13c9dca5.def("__len__", method_pointer_6bfb0444a5295ae9b96748554a27e3b8, "");
    class_bd3aca0a6cb156939d673e0f13c9dca5.def("begin", method_pointer_7c94877159bf5284a0356addd6b3b05c, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_bd3aca0a6cb156939d673e0f13c9dca5.def("end", method_pointer_8a840a970901512fa23c5ee45013b5cb, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");

}