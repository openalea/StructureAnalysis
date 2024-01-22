#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_975a242cbfb4576bbfb2e69cc63a9ca2(class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< int, class ::std::allocator< int > > * > > & instance, ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< int, class ::std::allocator< int > > * > >::difference_type  param_in_0, const class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< int, class ::std::allocator< int > > * > > & param_out) { instance.operator+=(param_in_0) = param_out; }
    void method_decorator_ccec540d50e556c0bd6fc3195dd1fe02(class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< int, class ::std::allocator< int > > * > > & instance, const class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< int, class ::std::allocator< int > > * > > & param_out) { instance.operator++() = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< int, class ::std::allocator< int > > * > > const volatile * get_pointer<class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< int, class ::std::allocator< int > > * > > const volatile >(class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< int, class ::std::allocator< int > > * > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_72b4985b076f54f0a30d51f8f5f500be()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< int, class ::std::allocator< int > > * > > & (::std::move_iterator< ::std::__wrap_iter< ::std::vector< int, ::std::allocator< int > > * > >::*method_pointer_975a242cbfb4576bbfb2e69cc63a9ca2)(::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< int, class ::std::allocator< int > > * > >::difference_type ) = &::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< int, class ::std::allocator< int > > * > >::operator+=;
    class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< int, class ::std::allocator< int > > * > > & (::std::move_iterator< ::std::__wrap_iter< ::std::vector< int, ::std::allocator< int > > * > >::*method_pointer_ccec540d50e556c0bd6fc3195dd1fe02)() = &::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< int, class ::std::allocator< int > > * > >::operator++;
    boost::python::class_< class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< int, class ::std::allocator< int > > * > >, autowig::Held< class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< int, class ::std::allocator< int > > * > > >::Type > class_72b4985b076f54f0a30d51f8f5f500be("_MoveIterator_72b4985b076f54f0a30d51f8f5f500be", "", boost::python::no_init);
    class_72b4985b076f54f0a30d51f8f5f500be.def("__iadd__", method_pointer_975a242cbfb4576bbfb2e69cc63a9ca2, boost::python::return_internal_reference<>(), "");
    class_72b4985b076f54f0a30d51f8f5f500be.def("__iadd__", autowig::method_decorator_975a242cbfb4576bbfb2e69cc63a9ca2);
    class_72b4985b076f54f0a30d51f8f5f500be.def("__next__", method_pointer_ccec540d50e556c0bd6fc3195dd1fe02, boost::python::return_internal_reference<>(), "");
    class_72b4985b076f54f0a30d51f8f5f500be.def("__next__", autowig::method_decorator_ccec540d50e556c0bd6fc3195dd1fe02);

}