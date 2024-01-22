#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_55736e7affca5b5fb21779515a8ee547(class ::std::move_iterator< class ::std::__wrap_iter< enum ::stat_tool::discrete_parametric * > > & instance, ::std::move_iterator< class ::std::__wrap_iter< enum ::stat_tool::discrete_parametric * > >::difference_type  param_in_0, const class ::std::move_iterator< class ::std::__wrap_iter< enum ::stat_tool::discrete_parametric * > > & param_out) { instance.operator+=(param_in_0) = param_out; }
    void method_decorator_7bc02e8a455857a8a3d9640af8c63dcb(class ::std::move_iterator< class ::std::__wrap_iter< enum ::stat_tool::discrete_parametric * > > & instance, const class ::std::move_iterator< class ::std::__wrap_iter< enum ::stat_tool::discrete_parametric * > > & param_out) { instance.operator++() = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::move_iterator< class ::std::__wrap_iter< enum ::stat_tool::discrete_parametric * > > const volatile * get_pointer<class ::std::move_iterator< class ::std::__wrap_iter< enum ::stat_tool::discrete_parametric * > > const volatile >(class ::std::move_iterator< class ::std::__wrap_iter< enum ::stat_tool::discrete_parametric * > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_846e84372f4157b69a2438c203b1a0f2()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    class ::std::move_iterator< class ::std::__wrap_iter< enum ::stat_tool::discrete_parametric * > > & (::std::move_iterator< ::std::__wrap_iter< enum ::stat_tool::discrete_parametric * > >::*method_pointer_55736e7affca5b5fb21779515a8ee547)(::std::move_iterator< class ::std::__wrap_iter< enum ::stat_tool::discrete_parametric * > >::difference_type ) = &::std::move_iterator< class ::std::__wrap_iter< enum ::stat_tool::discrete_parametric * > >::operator+=;
    class ::std::move_iterator< class ::std::__wrap_iter< enum ::stat_tool::discrete_parametric * > > & (::std::move_iterator< ::std::__wrap_iter< enum ::stat_tool::discrete_parametric * > >::*method_pointer_7bc02e8a455857a8a3d9640af8c63dcb)() = &::std::move_iterator< class ::std::__wrap_iter< enum ::stat_tool::discrete_parametric * > >::operator++;
    boost::python::class_< class ::std::move_iterator< class ::std::__wrap_iter< enum ::stat_tool::discrete_parametric * > >, autowig::Held< class ::std::move_iterator< class ::std::__wrap_iter< enum ::stat_tool::discrete_parametric * > > >::Type > class_846e84372f4157b69a2438c203b1a0f2("_MoveIterator_846e84372f4157b69a2438c203b1a0f2", "", boost::python::no_init);
    class_846e84372f4157b69a2438c203b1a0f2.def("__iadd__", method_pointer_55736e7affca5b5fb21779515a8ee547, boost::python::return_internal_reference<>(), "");
    class_846e84372f4157b69a2438c203b1a0f2.def("__iadd__", autowig::method_decorator_55736e7affca5b5fb21779515a8ee547);
    class_846e84372f4157b69a2438c203b1a0f2.def("__next__", method_pointer_7bc02e8a455857a8a3d9640af8c63dcb, boost::python::return_internal_reference<>(), "");
    class_846e84372f4157b69a2438c203b1a0f2.def("__next__", autowig::method_decorator_7bc02e8a455857a8a3d9640af8c63dcb);

}