#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_21f4804dbdad52a395f37d9446658957(class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< double, class ::std::allocator< double > > * > > & instance, ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< double, class ::std::allocator< double > > * > >::difference_type  param_in_0, const class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< double, class ::std::allocator< double > > * > > & param_out) { instance.operator+=(param_in_0) = param_out; }
    void method_decorator_3c3424e1b9535d8e937eefc501e0a970(class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< double, class ::std::allocator< double > > * > > & instance, const class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< double, class ::std::allocator< double > > * > > & param_out) { instance.operator++() = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< double, class ::std::allocator< double > > * > > const volatile * get_pointer<class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< double, class ::std::allocator< double > > * > > const volatile >(class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< double, class ::std::allocator< double > > * > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_4bc3e7e222065c4eb16729fc2579f504()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< double, class ::std::allocator< double > > * > > & (::std::move_iterator< ::std::__wrap_iter< ::std::vector< double, ::std::allocator< double > > * > >::*method_pointer_21f4804dbdad52a395f37d9446658957)(::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< double, class ::std::allocator< double > > * > >::difference_type ) = &::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< double, class ::std::allocator< double > > * > >::operator+=;
    class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< double, class ::std::allocator< double > > * > > & (::std::move_iterator< ::std::__wrap_iter< ::std::vector< double, ::std::allocator< double > > * > >::*method_pointer_3c3424e1b9535d8e937eefc501e0a970)() = &::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< double, class ::std::allocator< double > > * > >::operator++;
    boost::python::class_< class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< double, class ::std::allocator< double > > * > >, autowig::Held< class ::std::move_iterator< class ::std::__wrap_iter< class ::std::vector< double, class ::std::allocator< double > > * > > >::Type > class_4bc3e7e222065c4eb16729fc2579f504("_MoveIterator_4bc3e7e222065c4eb16729fc2579f504", "", boost::python::no_init);
    class_4bc3e7e222065c4eb16729fc2579f504.def("__iadd__", method_pointer_21f4804dbdad52a395f37d9446658957, boost::python::return_internal_reference<>(), "");
    class_4bc3e7e222065c4eb16729fc2579f504.def("__iadd__", autowig::method_decorator_21f4804dbdad52a395f37d9446658957);
    class_4bc3e7e222065c4eb16729fc2579f504.def("__next__", method_pointer_3c3424e1b9535d8e937eefc501e0a970, boost::python::return_internal_reference<>(), "");
    class_4bc3e7e222065c4eb16729fc2579f504.def("__next__", autowig::method_decorator_3c3424e1b9535d8e937eefc501e0a970);

}