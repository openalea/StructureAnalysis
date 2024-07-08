#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_039334763cb4576ea85791856bf8a255(class ::std::move_iterator< class ::std::__wrap_iter< double * > > & instance, ::std::move_iterator< class ::std::__wrap_iter< double * > >::difference_type  param_in_0, const class ::std::move_iterator< class ::std::__wrap_iter< double * > > & param_out) { instance.operator+=(param_in_0) = param_out; }
    void method_decorator_fbbd65b4ed735a87a619fe0ac9f59309(class ::std::move_iterator< class ::std::__wrap_iter< double * > > & instance, const class ::std::move_iterator< class ::std::__wrap_iter< double * > > & param_out) { instance.operator++() = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::move_iterator< class ::std::__wrap_iter< double * > > const volatile * get_pointer<class ::std::move_iterator< class ::std::__wrap_iter< double * > > const volatile >(class ::std::move_iterator< class ::std::__wrap_iter< double * > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_0e8aa013dbd85cbb828bb636ba93091e()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    class ::std::move_iterator< class ::std::__wrap_iter< double * > > & (::std::move_iterator< ::std::__wrap_iter< double * > >::*method_pointer_039334763cb4576ea85791856bf8a255)(::std::move_iterator< class ::std::__wrap_iter< double * > >::difference_type ) = &::std::move_iterator< class ::std::__wrap_iter< double * > >::operator+=;
    class ::std::move_iterator< class ::std::__wrap_iter< double * > > & (::std::move_iterator< ::std::__wrap_iter< double * > >::*method_pointer_fbbd65b4ed735a87a619fe0ac9f59309)() = &::std::move_iterator< class ::std::__wrap_iter< double * > >::operator++;
    boost::python::class_< class ::std::move_iterator< class ::std::__wrap_iter< double * > >, autowig::Held< class ::std::move_iterator< class ::std::__wrap_iter< double * > > >::Type > class_0e8aa013dbd85cbb828bb636ba93091e("_MoveIterator_0e8aa013dbd85cbb828bb636ba93091e", "", boost::python::no_init);
    class_0e8aa013dbd85cbb828bb636ba93091e.def("__iadd__", method_pointer_039334763cb4576ea85791856bf8a255, boost::python::return_internal_reference<>(), "");
    class_0e8aa013dbd85cbb828bb636ba93091e.def("__iadd__", autowig::method_decorator_039334763cb4576ea85791856bf8a255);
    class_0e8aa013dbd85cbb828bb636ba93091e.def("__next__", method_pointer_fbbd65b4ed735a87a619fe0ac9f59309, boost::python::return_internal_reference<>(), "");
    class_0e8aa013dbd85cbb828bb636ba93091e.def("__next__", autowig::method_decorator_fbbd65b4ed735a87a619fe0ac9f59309);

}