#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_f9d7e15c83925f859767232849add78d(class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::FrequencyDistribution * > > & instance, ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::FrequencyDistribution * > >::difference_type  param_in_0, const class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::FrequencyDistribution * > > & param_out) { instance.operator+=(param_in_0) = param_out; }
    void method_decorator_0416fe07542b56fab97d52a965a98e9d(class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::FrequencyDistribution * > > & instance, const class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::FrequencyDistribution * > > & param_out) { instance.operator++() = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::FrequencyDistribution * > > const volatile * get_pointer<class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::FrequencyDistribution * > > const volatile >(class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::FrequencyDistribution * > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_1c84c65db02e5441bd8e625042563a4d()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::FrequencyDistribution * > > & (::std::move_iterator< ::std::__wrap_iter< ::stat_tool::FrequencyDistribution * > >::*method_pointer_f9d7e15c83925f859767232849add78d)(::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::FrequencyDistribution * > >::difference_type ) = &::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::FrequencyDistribution * > >::operator+=;
    class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::FrequencyDistribution * > > & (::std::move_iterator< ::std::__wrap_iter< ::stat_tool::FrequencyDistribution * > >::*method_pointer_0416fe07542b56fab97d52a965a98e9d)() = &::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::FrequencyDistribution * > >::operator++;
    boost::python::class_< class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::FrequencyDistribution * > >, autowig::Held< class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::FrequencyDistribution * > > >::Type > class_1c84c65db02e5441bd8e625042563a4d("_MoveIterator_1c84c65db02e5441bd8e625042563a4d", "", boost::python::no_init);
    class_1c84c65db02e5441bd8e625042563a4d.def("__iadd__", method_pointer_f9d7e15c83925f859767232849add78d, boost::python::return_internal_reference<>(), "");
    class_1c84c65db02e5441bd8e625042563a4d.def("__iadd__", autowig::method_decorator_f9d7e15c83925f859767232849add78d);
    class_1c84c65db02e5441bd8e625042563a4d.def("__next__", method_pointer_0416fe07542b56fab97d52a965a98e9d, boost::python::return_internal_reference<>(), "");
    class_1c84c65db02e5441bd8e625042563a4d.def("__next__", autowig::method_decorator_0416fe07542b56fab97d52a965a98e9d);

}