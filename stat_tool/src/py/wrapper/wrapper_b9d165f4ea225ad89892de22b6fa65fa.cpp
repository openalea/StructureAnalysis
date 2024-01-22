#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_3dbf6d478adb54118a1bfc852d94f46d(class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::Vectors * > > & instance, ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::Vectors * > >::difference_type  param_in_0, const class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::Vectors * > > & param_out) { instance.operator+=(param_in_0) = param_out; }
    void method_decorator_9908e97917e35e7c8e22da9cccab260b(class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::Vectors * > > & instance, const class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::Vectors * > > & param_out) { instance.operator++() = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::Vectors * > > const volatile * get_pointer<class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::Vectors * > > const volatile >(class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::Vectors * > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_b9d165f4ea225ad89892de22b6fa65fa()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::Vectors * > > & (::std::move_iterator< ::std::__wrap_iter< ::stat_tool::Vectors * > >::*method_pointer_3dbf6d478adb54118a1bfc852d94f46d)(::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::Vectors * > >::difference_type ) = &::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::Vectors * > >::operator+=;
    class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::Vectors * > > & (::std::move_iterator< ::std::__wrap_iter< ::stat_tool::Vectors * > >::*method_pointer_9908e97917e35e7c8e22da9cccab260b)() = &::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::Vectors * > >::operator++;
    boost::python::class_< class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::Vectors * > >, autowig::Held< class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::Vectors * > > >::Type > class_b9d165f4ea225ad89892de22b6fa65fa("_MoveIterator_b9d165f4ea225ad89892de22b6fa65fa", "", boost::python::no_init);
    class_b9d165f4ea225ad89892de22b6fa65fa.def("__iadd__", method_pointer_3dbf6d478adb54118a1bfc852d94f46d, boost::python::return_internal_reference<>(), "");
    class_b9d165f4ea225ad89892de22b6fa65fa.def("__iadd__", autowig::method_decorator_3dbf6d478adb54118a1bfc852d94f46d);
    class_b9d165f4ea225ad89892de22b6fa65fa.def("__next__", method_pointer_9908e97917e35e7c8e22da9cccab260b, boost::python::return_internal_reference<>(), "");
    class_b9d165f4ea225ad89892de22b6fa65fa.def("__next__", autowig::method_decorator_9908e97917e35e7c8e22da9cccab260b);

}