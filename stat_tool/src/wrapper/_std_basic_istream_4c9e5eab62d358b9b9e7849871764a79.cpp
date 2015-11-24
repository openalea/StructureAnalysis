#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_basic_istream_4c9e5eab62d358b9b9e7849871764a79()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        class ::std::basic_istream<char, std::char_traits<char> > & (::std::basic_istream<char, std::char_traits<char> >::*method_pointer_0d4eece2efc25188b39238d091820df9)(class ::std::basic_streambuf<char, std::char_traits<char> > *) = &::std::basic_istream<char, std::char_traits<char> >::operator>>;
        long (::std::basic_istream<char, std::char_traits<char> >::*method_pointer_87fd1f4abb075a10a9dc5b72cb6c9ed0)() const = &::std::basic_istream<char, std::char_traits<char> >::gcount;
        int (::std::basic_istream<char, std::char_traits<char> >::*method_pointer_09ebe6dcfe2f5e12a2e7c0681fe5fd97)() = &::std::basic_istream<char, std::char_traits<char> >::get;
        class ::std::basic_istream<char, std::char_traits<char> > & (::std::basic_istream<char, std::char_traits<char> >::*method_pointer_db525be836955981a442b64d06de4f8a)(class ::std::basic_streambuf<char, std::char_traits<char> > &, char) = &::std::basic_istream<char, std::char_traits<char> >::get;
        class ::std::basic_istream<char, std::char_traits<char> > & (::std::basic_istream<char, std::char_traits<char> >::*method_pointer_db65491293e75417a7ed00a898476f50)(class ::std::basic_streambuf<char, std::char_traits<char> > &) = &::std::basic_istream<char, std::char_traits<char> >::get;
        class ::std::basic_istream<char, std::char_traits<char> > & (::std::basic_istream<char, std::char_traits<char> >::*method_pointer_444f2cab903a591a99a79cf9ce44acdb)(long, int) = &::std::basic_istream<char, std::char_traits<char> >::ignore;
        class ::std::basic_istream<char, std::char_traits<char> > & (::std::basic_istream<char, std::char_traits<char> >::*method_pointer_95cbdbee64b559579e6810cf0fb02240)(long) = &::std::basic_istream<char, std::char_traits<char> >::ignore;
        class ::std::basic_istream<char, std::char_traits<char> > & (::std::basic_istream<char, std::char_traits<char> >::*method_pointer_1ab68372ce895e18845e0fcdca579c9f)() = &::std::basic_istream<char, std::char_traits<char> >::ignore;
        int (::std::basic_istream<char, std::char_traits<char> >::*method_pointer_c455066a363a51948f128f19e530315e)() = &::std::basic_istream<char, std::char_traits<char> >::peek;
        class ::std::basic_istream<char, std::char_traits<char> > & (::std::basic_istream<char, std::char_traits<char> >::*method_pointer_195ad6eb26fa5983bb4312e5ffda8381)(char) = &::std::basic_istream<char, std::char_traits<char> >::putback;
        class ::std::basic_istream<char, std::char_traits<char> > & (::std::basic_istream<char, std::char_traits<char> >::*method_pointer_434888c1ac3259b4becb5f0662697391)() = &::std::basic_istream<char, std::char_traits<char> >::unget;
        int (::std::basic_istream<char, std::char_traits<char> >::*method_pointer_f6faf835a89255a98b31e3cd8820bc0a)() = &::std::basic_istream<char, std::char_traits<char> >::sync;
        class ::std::fpos<__mbstate_t> (::std::basic_istream<char, std::char_traits<char> >::*method_pointer_bb7929eca6f750b1974d57c75ef07923)() = &::std::basic_istream<char, std::char_traits<char> >::tellg;
        class ::std::basic_istream<char, std::char_traits<char> > & (::std::basic_istream<char, std::char_traits<char> >::*method_pointer_b2546caa87fc5fd38b044fcc514f80ff)(class ::std::fpos<__mbstate_t>) = &::std::basic_istream<char, std::char_traits<char> >::seekg;
        class ::std::basic_istream<char, std::char_traits<char> > & (::std::basic_istream<char, std::char_traits<char> >::*method_pointer_a7f4010e62b35daa963de03fbc8e2e25)(long, enum ::std::_Ios_Seekdir) = &::std::basic_istream<char, std::char_traits<char> >::seekg;
        boost::python::class_< class ::std::basic_istream<char, std::char_traits<char> >, std::shared_ptr< class ::std::basic_istream<char, std::char_traits<char> > >, boost::python::bases< class ::std::basic_ios<char, std::char_traits<char> > >, boost::noncopyable >("_BasicIstream_4c9e5eab62d358b9b9e7849871764a79", boost::python::no_init)
            .def(boost::python::init< class ::std::basic_streambuf<char, std::char_traits<char> > * >())
            .def("__rshift__", method_pointer_0d4eece2efc25188b39238d091820df9, boost::python::return_internal_reference<>())
            .def("gcount", method_pointer_87fd1f4abb075a10a9dc5b72cb6c9ed0)
            .def("get", method_pointer_09ebe6dcfe2f5e12a2e7c0681fe5fd97)
            .def("get", method_pointer_db525be836955981a442b64d06de4f8a, boost::python::return_internal_reference<>())
            .def("get", method_pointer_db65491293e75417a7ed00a898476f50, boost::python::return_internal_reference<>())
            .def("ignore", method_pointer_444f2cab903a591a99a79cf9ce44acdb, boost::python::return_internal_reference<>())
            .def("ignore", method_pointer_95cbdbee64b559579e6810cf0fb02240, boost::python::return_internal_reference<>())
            .def("ignore", method_pointer_1ab68372ce895e18845e0fcdca579c9f, boost::python::return_internal_reference<>())
            .def("peek", method_pointer_c455066a363a51948f128f19e530315e)
            .def("putback", method_pointer_195ad6eb26fa5983bb4312e5ffda8381, boost::python::return_internal_reference<>())
            .def("unget", method_pointer_434888c1ac3259b4becb5f0662697391, boost::python::return_internal_reference<>())
            .def("sync", method_pointer_f6faf835a89255a98b31e3cd8820bc0a)
            .def("tellg", method_pointer_bb7929eca6f750b1974d57c75ef07923)
            .def("seekg", method_pointer_b2546caa87fc5fd38b044fcc514f80ff, boost::python::return_internal_reference<>())
            .def("seekg", method_pointer_a7f4010e62b35daa963de03fbc8e2e25, boost::python::return_internal_reference<>());
        boost::python::implicitly_convertible< std::shared_ptr< class ::std::basic_istream<char, std::char_traits<char> > >, std::shared_ptr< class ::std::basic_ios<char, std::char_traits<char> > > >();
}