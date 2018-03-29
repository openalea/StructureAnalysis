#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::initializer_list< char > const volatile * get_pointer<class ::std::initializer_list< char > const volatile >(class ::std::initializer_list< char > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_821f52165ab1565093f9ea5b7cdd5401()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    ::std::initializer_list< char >::size_type  (::std::initializer_list< char >::*method_pointer_c1384da902e05873b9a5782e53b9e197)() const = &::std::initializer_list< char >::size;
    boost::python::class_< class ::std::initializer_list< char >, autowig::Held< class ::std::initializer_list< char > >::Type > class_821f52165ab1565093f9ea5b7cdd5401("_InitializerList_821f52165ab1565093f9ea5b7cdd5401", "", boost::python::no_init);
    class_821f52165ab1565093f9ea5b7cdd5401.def("__len__", method_pointer_c1384da902e05873b9a5782e53b9e197, "");

}