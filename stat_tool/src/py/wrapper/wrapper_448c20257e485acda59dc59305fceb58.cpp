#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const volatile * get_pointer<class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const volatile >(class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_448c20257e485acda59dc59305fceb58()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::class_< class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > >, autowig::Held< class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > >::Type > class_448c20257e485acda59dc59305fceb58("_BasicString_448c20257e485acda59dc59305fceb58", "", boost::python::no_init);
    class_448c20257e485acda59dc59305fceb58.def_readonly("npos", ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > >::npos, "");

}