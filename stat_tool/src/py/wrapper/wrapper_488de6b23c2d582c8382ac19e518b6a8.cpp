#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::ctype< char > const volatile * get_pointer<class ::std::ctype< char > const volatile >(class ::std::ctype< char > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_488de6b23c2d582c8382ac19e518b6a8()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    bool  (::std::ctype< char >::*method_pointer_9dd924e9316e5530bf52d195c4a12410)(::std::ctype_base::mask , char ) const = &::std::ctype< char >::is;
    ::std::ctype< char >::char_type  (::std::ctype< char >::*method_pointer_257942a509e35d0d9ebec35f44bdc3a9)(::std::ctype< char >::char_type ) const = &::std::ctype< char >::toupper;
    ::std::ctype< char >::char_type  (::std::ctype< char >::*method_pointer_de1db92fef075cff9ae73cb8c7cd384a)(::std::ctype< char >::char_type ) const = &::std::ctype< char >::tolower;
    ::std::ctype< char >::char_type  (::std::ctype< char >::*method_pointer_c07755c6029c5d788287a798b839abc5)(char ) const = &::std::ctype< char >::widen;
    char  (::std::ctype< char >::*method_pointer_625e9e2026b75458b11f4355faaa4269)(::std::ctype< char >::char_type , char ) const = &::std::ctype< char >::narrow;
    boost::python::class_< class ::std::ctype< char >, autowig::Held< class ::std::ctype< char > >::Type, boost::python::bases< class ::std::locale::facet, struct ::std::ctype_base >, boost::noncopyable > class_488de6b23c2d582c8382ac19e518b6a8("_Ctype_488de6b23c2d582c8382ac19e518b6a8", "", boost::python::no_init);
    class_488de6b23c2d582c8382ac19e518b6a8.def("is", method_pointer_9dd924e9316e5530bf52d195c4a12410, "");
    class_488de6b23c2d582c8382ac19e518b6a8.def("toupper", method_pointer_257942a509e35d0d9ebec35f44bdc3a9, "");
    class_488de6b23c2d582c8382ac19e518b6a8.def("tolower", method_pointer_de1db92fef075cff9ae73cb8c7cd384a, "");
    class_488de6b23c2d582c8382ac19e518b6a8.def("widen", method_pointer_c07755c6029c5d788287a798b839abc5, "");
    class_488de6b23c2d582c8382ac19e518b6a8.def("narrow", method_pointer_625e9e2026b75458b11f4355faaa4269, "");
    class_488de6b23c2d582c8382ac19e518b6a8.def_readonly("table_size", ::std::ctype< char >::table_size, "");

    if(autowig::Held< class ::std::ctype< char > >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::std::ctype< char > >::Type, autowig::Held< class ::std::locale::facet >::Type >();
        boost::python::implicitly_convertible< autowig::Held< class ::std::ctype< char > >::Type, autowig::Held< struct ::std::ctype_base >::Type >();
    }

}