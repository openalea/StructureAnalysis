#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_7dc3925030d45f30b171d8363ea53ab0(class ::std::move_iterator< class ::std::__wrap_iter< int * > > & instance, ::std::move_iterator< class ::std::__wrap_iter< int * > >::difference_type  param_in_0, const class ::std::move_iterator< class ::std::__wrap_iter< int * > > & param_out) { instance.operator+=(param_in_0) = param_out; }
    void method_decorator_ebe1517f9177532ebcbcf29d0f213772(class ::std::move_iterator< class ::std::__wrap_iter< int * > > & instance, const class ::std::move_iterator< class ::std::__wrap_iter< int * > > & param_out) { instance.operator++() = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::move_iterator< class ::std::__wrap_iter< int * > > const volatile * get_pointer<class ::std::move_iterator< class ::std::__wrap_iter< int * > > const volatile >(class ::std::move_iterator< class ::std::__wrap_iter< int * > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_745e5d0d0b025a76b955a27c7f788b56()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    class ::std::move_iterator< class ::std::__wrap_iter< int * > > & (::std::move_iterator< ::std::__wrap_iter< int * > >::*method_pointer_7dc3925030d45f30b171d8363ea53ab0)(::std::move_iterator< class ::std::__wrap_iter< int * > >::difference_type ) = &::std::move_iterator< class ::std::__wrap_iter< int * > >::operator+=;
    class ::std::move_iterator< class ::std::__wrap_iter< int * > > & (::std::move_iterator< ::std::__wrap_iter< int * > >::*method_pointer_ebe1517f9177532ebcbcf29d0f213772)() = &::std::move_iterator< class ::std::__wrap_iter< int * > >::operator++;
    boost::python::class_< class ::std::move_iterator< class ::std::__wrap_iter< int * > >, autowig::Held< class ::std::move_iterator< class ::std::__wrap_iter< int * > > >::Type > class_745e5d0d0b025a76b955a27c7f788b56("_MoveIterator_745e5d0d0b025a76b955a27c7f788b56", "", boost::python::no_init);
    class_745e5d0d0b025a76b955a27c7f788b56.def("__iadd__", method_pointer_7dc3925030d45f30b171d8363ea53ab0, boost::python::return_internal_reference<>(), "");
    class_745e5d0d0b025a76b955a27c7f788b56.def("__iadd__", autowig::method_decorator_7dc3925030d45f30b171d8363ea53ab0);
    class_745e5d0d0b025a76b955a27c7f788b56.def("__next__", method_pointer_ebe1517f9177532ebcbcf29d0f213772, boost::python::return_internal_reference<>(), "");
    class_745e5d0d0b025a76b955a27c7f788b56.def("__next__", autowig::method_decorator_ebe1517f9177532ebcbcf29d0f213772);

}