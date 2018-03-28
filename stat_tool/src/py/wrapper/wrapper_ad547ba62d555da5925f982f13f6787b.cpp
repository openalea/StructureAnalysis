#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> struct ::std::pair< float, float > const volatile * get_pointer<struct ::std::pair< float, float > const volatile >(struct ::std::pair< float, float > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_ad547ba62d555da5925f982f13f6787b()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    void  (::std::pair< float, float >::*method_pointer_3af1e5c2286a53edb9924dabe24b751c)(struct ::std::pair< float, float > &) = &::std::pair< float, float >::swap;
    boost::python::class_< struct ::std::pair< float, float >, autowig::Held< struct ::std::pair< float, float > >::Type > class_ad547ba62d555da5925f982f13f6787b("_Pair_ad547ba62d555da5925f982f13f6787b", "", boost::python::no_init);
    class_ad547ba62d555da5925f982f13f6787b.def(boost::python::init< struct ::std::pair< float, float > const & >(""));
    class_ad547ba62d555da5925f982f13f6787b.def("swap", method_pointer_3af1e5c2286a53edb9924dabe24b751c, "");
    class_ad547ba62d555da5925f982f13f6787b.def_readwrite("first", &::std::pair< float, float >::first, "");
    class_ad547ba62d555da5925f982f13f6787b.def_readwrite("second", &::std::pair< float, float >::second, "");

}