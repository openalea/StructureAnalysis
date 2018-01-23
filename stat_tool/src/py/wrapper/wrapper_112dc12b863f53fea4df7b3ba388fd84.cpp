#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::basic_streambuf< char, struct ::std::char_traits< char > > const volatile * get_pointer<class ::std::basic_streambuf< char, struct ::std::char_traits< char > > const volatile >(class ::std::basic_streambuf< char, struct ::std::char_traits< char > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_112dc12b863f53fea4df7b3ba388fd84()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    class ::std::locale  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_73c494a2dfb4517985232b211533200c)(class ::std::locale const &) = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::pubimbue;
    class ::std::locale  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_5737de16043a5da78137948c088d3257)() const = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::getloc;
    int  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_93b4b685fee151f1be0336412d779ca5)() = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::pubsync;
    ::std::streamsize  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_4e292143959b5ff8b98ac3e995006f1e)() = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::in_avail;
    ::std::basic_streambuf< char, struct ::std::char_traits< char > >::int_type  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_ae703141ab07573fa22ed3f7cb36685b)() = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::snextc;
    ::std::basic_streambuf< char, struct ::std::char_traits< char > >::int_type  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_5f21eb74f7e25ac3b1d3250327b983f9)() = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::sbumpc;
    ::std::basic_streambuf< char, struct ::std::char_traits< char > >::int_type  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_270c3ce800bd51fb89b29421337e1cf9)() = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::sgetc;
    ::std::basic_streambuf< char, struct ::std::char_traits< char > >::int_type  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_c8f77cf365025e30b293b5402ef836c9)(::std::basic_streambuf< char, struct ::std::char_traits< char > >::char_type ) = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::sputbackc;
    ::std::basic_streambuf< char, struct ::std::char_traits< char > >::int_type  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_91712faddf415594898ed51802787bb7)() = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::sungetc;
    ::std::basic_streambuf< char, struct ::std::char_traits< char > >::int_type  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_c0ccc799bd54556c894d54f4f3b142e2)(::std::basic_streambuf< char, struct ::std::char_traits< char > >::char_type ) = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::sputc;
    void  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_e75b373786e35733a0a37732790ad698)() = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::stossc;
    void  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_3c5548ba706e54bd9826734a04fb8ac5)(::std::streamsize ) = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::__safe_gbump;
    void  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_f5a1f600482f58f5973cdcf5cd0e5aa3)(::std::streamsize ) = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::__safe_pbump;
    boost::python::class_< class ::std::basic_streambuf< char, struct ::std::char_traits< char > >, autowig::Held< class ::std::basic_streambuf< char, struct ::std::char_traits< char > > >::Type, boost::noncopyable > class_112dc12b863f53fea4df7b3ba388fd84("_BasicStreambuf_112dc12b863f53fea4df7b3ba388fd84", "", boost::python::no_init);
    class_112dc12b863f53fea4df7b3ba388fd84.def("pubimbue", method_pointer_73c494a2dfb4517985232b211533200c, "");
    class_112dc12b863f53fea4df7b3ba388fd84.def("getloc", method_pointer_5737de16043a5da78137948c088d3257, "");
    class_112dc12b863f53fea4df7b3ba388fd84.def("pubsync", method_pointer_93b4b685fee151f1be0336412d779ca5, "");
    class_112dc12b863f53fea4df7b3ba388fd84.def("in_avail", method_pointer_4e292143959b5ff8b98ac3e995006f1e, "");
    class_112dc12b863f53fea4df7b3ba388fd84.def("snextc", method_pointer_ae703141ab07573fa22ed3f7cb36685b, "");
    class_112dc12b863f53fea4df7b3ba388fd84.def("sbumpc", method_pointer_5f21eb74f7e25ac3b1d3250327b983f9, "");
    class_112dc12b863f53fea4df7b3ba388fd84.def("sgetc", method_pointer_270c3ce800bd51fb89b29421337e1cf9, "");
    class_112dc12b863f53fea4df7b3ba388fd84.def("sputbackc", method_pointer_c8f77cf365025e30b293b5402ef836c9, "");
    class_112dc12b863f53fea4df7b3ba388fd84.def("sungetc", method_pointer_91712faddf415594898ed51802787bb7, "");
    class_112dc12b863f53fea4df7b3ba388fd84.def("sputc", method_pointer_c0ccc799bd54556c894d54f4f3b142e2, "");
    class_112dc12b863f53fea4df7b3ba388fd84.def("stossc", method_pointer_e75b373786e35733a0a37732790ad698, "");
    class_112dc12b863f53fea4df7b3ba388fd84.def("safe_gbump", method_pointer_3c5548ba706e54bd9826734a04fb8ac5, "");
    class_112dc12b863f53fea4df7b3ba388fd84.def("safe_pbump", method_pointer_f5a1f600482f58f5973cdcf5cd0e5aa3, "");

}