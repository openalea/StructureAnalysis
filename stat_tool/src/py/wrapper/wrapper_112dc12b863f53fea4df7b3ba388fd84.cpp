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
    class ::std::locale  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_54276d56823054afb84691b4b124fde8)(class ::std::locale const &) = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::pubimbue;
    class ::std::locale  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_a99fbe59a5d357eeaf3fa986718bfe5c)() const = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::getloc;
    int  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_5ff19440f513591e8aba0bc6d4b01fa0)() = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::pubsync;
    ::std::streamsize  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_8ba9b8dc5c6d5e36bd89826496c613f2)() = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::in_avail;
    ::std::basic_streambuf< char, struct ::std::char_traits< char > >::int_type  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_a738a4ddd45a58478eb5f529e3f3ef07)() = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::snextc;
    ::std::basic_streambuf< char, struct ::std::char_traits< char > >::int_type  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_0a3e640f01895954b7cbb4614b093ab3)() = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::sbumpc;
    ::std::basic_streambuf< char, struct ::std::char_traits< char > >::int_type  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_de5234c42d2d552990d1ca0cf39e10f1)() = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::sgetc;
    ::std::basic_streambuf< char, struct ::std::char_traits< char > >::int_type  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_47441f24965d5428881f29c16b2da680)(::std::basic_streambuf< char, struct ::std::char_traits< char > >::char_type ) = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::sputbackc;
    ::std::basic_streambuf< char, struct ::std::char_traits< char > >::int_type  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_3df14853796b5377ba3d8488530b1139)() = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::sungetc;
    ::std::basic_streambuf< char, struct ::std::char_traits< char > >::int_type  (::std::basic_streambuf< char, ::std::char_traits< char > >::*method_pointer_223618bfb2b05a71a17b8305d93622c8)(::std::basic_streambuf< char, struct ::std::char_traits< char > >::char_type ) = &::std::basic_streambuf< char, struct ::std::char_traits< char > >::sputc;
    boost::python::class_< class ::std::basic_streambuf< char, struct ::std::char_traits< char > >, autowig::Held< class ::std::basic_streambuf< char, struct ::std::char_traits< char > > >::Type, boost::noncopyable > class_112dc12b863f53fea4df7b3ba388fd84("_BasicStreambuf_112dc12b863f53fea4df7b3ba388fd84", "", boost::python::no_init);
    class_112dc12b863f53fea4df7b3ba388fd84.def("pubimbue", method_pointer_54276d56823054afb84691b4b124fde8, "");
    class_112dc12b863f53fea4df7b3ba388fd84.def("getloc", method_pointer_a99fbe59a5d357eeaf3fa986718bfe5c, "");
    class_112dc12b863f53fea4df7b3ba388fd84.def("pubsync", method_pointer_5ff19440f513591e8aba0bc6d4b01fa0, "");
    class_112dc12b863f53fea4df7b3ba388fd84.def("in_avail", method_pointer_8ba9b8dc5c6d5e36bd89826496c613f2, "");
    class_112dc12b863f53fea4df7b3ba388fd84.def("snextc", method_pointer_a738a4ddd45a58478eb5f529e3f3ef07, "");
    class_112dc12b863f53fea4df7b3ba388fd84.def("sbumpc", method_pointer_0a3e640f01895954b7cbb4614b093ab3, "");
    class_112dc12b863f53fea4df7b3ba388fd84.def("sgetc", method_pointer_de5234c42d2d552990d1ca0cf39e10f1, "");
    class_112dc12b863f53fea4df7b3ba388fd84.def("sputbackc", method_pointer_47441f24965d5428881f29c16b2da680, "");
    class_112dc12b863f53fea4df7b3ba388fd84.def("sungetc", method_pointer_3df14853796b5377ba3d8488530b1139, "");
    class_112dc12b863f53fea4df7b3ba388fd84.def("sputc", method_pointer_223618bfb2b05a71a17b8305d93622c8, "");

}