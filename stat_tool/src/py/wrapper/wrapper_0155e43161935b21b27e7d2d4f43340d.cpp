#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::initializer_list< int > const volatile * get_pointer<class ::std::initializer_list< int > const volatile >(class ::std::initializer_list< int > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_0155e43161935b21b27e7d2d4f43340d()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    ::std::initializer_list< int >::size_type  (::std::initializer_list< int >::*method_pointer_4a9f00d671275a2fa22a34f7cf32edd1)() const = &::std::initializer_list< int >::size;
    boost::python::class_< class ::std::initializer_list< int >, autowig::Held< class ::std::initializer_list< int > >::Type > class_0155e43161935b21b27e7d2d4f43340d("_InitializerList_0155e43161935b21b27e7d2d4f43340d", "", boost::python::no_init);
    class_0155e43161935b21b27e7d2d4f43340d.def("__len__", method_pointer_4a9f00d671275a2fa22a34f7cf32edd1, "");

}