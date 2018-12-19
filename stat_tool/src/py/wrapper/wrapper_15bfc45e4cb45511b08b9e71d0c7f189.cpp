#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> struct ::std::pair< struct ::std::pair< float, float >, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > > const volatile * get_pointer<struct ::std::pair< struct ::std::pair< float, float >, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > > const volatile >(struct ::std::pair< struct ::std::pair< float, float >, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_15bfc45e4cb45511b08b9e71d0c7f189()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    void  (::std::pair< ::std::pair< float, float >, ::std::basic_string< char, ::std::char_traits< char >, ::std::allocator< char > > >::*method_pointer_b07ed6ea233650e393937afc2a8b38f8)(struct ::std::pair< struct ::std::pair< float, float >, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > > &) = &::std::pair< struct ::std::pair< float, float >, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > >::swap;
    boost::python::class_< struct ::std::pair< struct ::std::pair< float, float >, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > >, autowig::Held< struct ::std::pair< struct ::std::pair< float, float >, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > > >::Type > class_15bfc45e4cb45511b08b9e71d0c7f189("_Pair_15bfc45e4cb45511b08b9e71d0c7f189", "", boost::python::no_init);
    class_15bfc45e4cb45511b08b9e71d0c7f189.def("swap", method_pointer_b07ed6ea233650e393937afc2a8b38f8, "");

}