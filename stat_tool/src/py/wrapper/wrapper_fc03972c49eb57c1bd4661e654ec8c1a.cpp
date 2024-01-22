#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_dde0d2e3ac7a5435a8d59edabe1a5393(class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::DiscreteParametric * > > & instance, ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::DiscreteParametric * > >::difference_type  param_in_0, const class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::DiscreteParametric * > > & param_out) { instance.operator+=(param_in_0) = param_out; }
    void method_decorator_c6c8c2961d1e58c8b05715e0e763b8f8(class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::DiscreteParametric * > > & instance, const class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::DiscreteParametric * > > & param_out) { instance.operator++() = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::DiscreteParametric * > > const volatile * get_pointer<class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::DiscreteParametric * > > const volatile >(class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::DiscreteParametric * > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_fc03972c49eb57c1bd4661e654ec8c1a()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::DiscreteParametric * > > & (::std::move_iterator< ::std::__wrap_iter< ::stat_tool::DiscreteParametric * > >::*method_pointer_dde0d2e3ac7a5435a8d59edabe1a5393)(::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::DiscreteParametric * > >::difference_type ) = &::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::DiscreteParametric * > >::operator+=;
    class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::DiscreteParametric * > > & (::std::move_iterator< ::std::__wrap_iter< ::stat_tool::DiscreteParametric * > >::*method_pointer_c6c8c2961d1e58c8b05715e0e763b8f8)() = &::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::DiscreteParametric * > >::operator++;
    boost::python::class_< class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::DiscreteParametric * > >, autowig::Held< class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::DiscreteParametric * > > >::Type > class_fc03972c49eb57c1bd4661e654ec8c1a("_MoveIterator_fc03972c49eb57c1bd4661e654ec8c1a", "", boost::python::no_init);
    class_fc03972c49eb57c1bd4661e654ec8c1a.def("__iadd__", method_pointer_dde0d2e3ac7a5435a8d59edabe1a5393, boost::python::return_internal_reference<>(), "");
    class_fc03972c49eb57c1bd4661e654ec8c1a.def("__iadd__", autowig::method_decorator_dde0d2e3ac7a5435a8d59edabe1a5393);
    class_fc03972c49eb57c1bd4661e654ec8c1a.def("__next__", method_pointer_c6c8c2961d1e58c8b05715e0e763b8f8, boost::python::return_internal_reference<>(), "");
    class_fc03972c49eb57c1bd4661e654ec8c1a.def("__next__", autowig::method_decorator_c6c8c2961d1e58c8b05715e0e763b8f8);

}