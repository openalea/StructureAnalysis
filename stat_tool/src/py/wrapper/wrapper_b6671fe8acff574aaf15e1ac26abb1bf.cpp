#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::initializer_list< class ::stat_tool::Vectors > const volatile * get_pointer<class ::std::initializer_list< class ::stat_tool::Vectors > const volatile >(class ::std::initializer_list< class ::stat_tool::Vectors > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_b6671fe8acff574aaf15e1ac26abb1bf()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    ::std::initializer_list< class ::stat_tool::Vectors >::size_type  (::std::initializer_list< ::stat_tool::Vectors >::*method_pointer_f0c38db54eaf5cccb8683f4b773f53d8)() const = &::std::initializer_list< class ::stat_tool::Vectors >::size;
    ::std::initializer_list< class ::stat_tool::Vectors >::const_iterator  (::std::initializer_list< ::stat_tool::Vectors >::*method_pointer_91277965f1c15c0fbf0405a5658731be)() const = &::std::initializer_list< class ::stat_tool::Vectors >::begin;
    ::std::initializer_list< class ::stat_tool::Vectors >::const_iterator  (::std::initializer_list< ::stat_tool::Vectors >::*method_pointer_6863d2cb34e855a5aa8947a73244c876)() const = &::std::initializer_list< class ::stat_tool::Vectors >::end;
    boost::python::class_< class ::std::initializer_list< class ::stat_tool::Vectors >, autowig::Held< class ::std::initializer_list< class ::stat_tool::Vectors > >::Type > class_b6671fe8acff574aaf15e1ac26abb1bf("_InitializerList_b6671fe8acff574aaf15e1ac26abb1bf", "", boost::python::no_init);
    class_b6671fe8acff574aaf15e1ac26abb1bf.def("__len__", method_pointer_f0c38db54eaf5cccb8683f4b773f53d8, "");
    class_b6671fe8acff574aaf15e1ac26abb1bf.def("begin", method_pointer_91277965f1c15c0fbf0405a5658731be, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_b6671fe8acff574aaf15e1ac26abb1bf.def("end", method_pointer_6863d2cb34e855a5aa8947a73244c876, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");

}