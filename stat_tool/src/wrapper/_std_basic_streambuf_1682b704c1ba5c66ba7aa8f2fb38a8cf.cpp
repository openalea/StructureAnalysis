#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_basic_streambuf_1682b704c1ba5c66ba7aa8f2fb38a8cf()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        class ::std::locale (::std::basic_streambuf<char, std::char_traits<char> >::*method_pointer_5f91820536bd56edaa6665a462e148db)(class ::std::locale const &) = &::std::basic_streambuf<char, std::char_traits<char> >::pubimbue;
        class ::std::locale (::std::basic_streambuf<char, std::char_traits<char> >::*method_pointer_c82af173b6a65520903dd8fbcd436978)() const = &::std::basic_streambuf<char, std::char_traits<char> >::getloc;
        int (::std::basic_streambuf<char, std::char_traits<char> >::*method_pointer_ba0b0ddea84950c183abbdace771a53b)() = &::std::basic_streambuf<char, std::char_traits<char> >::pubsync;
        long (::std::basic_streambuf<char, std::char_traits<char> >::*method_pointer_622b3e983c985691acc059e8a2648ec0)() = &::std::basic_streambuf<char, std::char_traits<char> >::in_avail;
        int (::std::basic_streambuf<char, std::char_traits<char> >::*method_pointer_bca0a9fc504f5e239c22df6761d6998e)() = &::std::basic_streambuf<char, std::char_traits<char> >::snextc;
        int (::std::basic_streambuf<char, std::char_traits<char> >::*method_pointer_d7573f6204a553e5a1f6150dad11ec95)() = &::std::basic_streambuf<char, std::char_traits<char> >::sbumpc;
        int (::std::basic_streambuf<char, std::char_traits<char> >::*method_pointer_0bf2d5f8826c5de1b429f6be002b6254)() = &::std::basic_streambuf<char, std::char_traits<char> >::sgetc;
        int (::std::basic_streambuf<char, std::char_traits<char> >::*method_pointer_6e8611294bf35fb8b5abb770f1d082aa)(char) = &::std::basic_streambuf<char, std::char_traits<char> >::sputbackc;
        int (::std::basic_streambuf<char, std::char_traits<char> >::*method_pointer_9eaa4a1508cf54e7a22544688cf0c840)() = &::std::basic_streambuf<char, std::char_traits<char> >::sungetc;
        int (::std::basic_streambuf<char, std::char_traits<char> >::*method_pointer_d4077a93f6f95daea17151c6f644a193)(char) = &::std::basic_streambuf<char, std::char_traits<char> >::sputc;
        void (::std::basic_streambuf<char, std::char_traits<char> >::*method_pointer_efac272f5d1f5c34b324f50800ed6b46)() = &::std::basic_streambuf<char, std::char_traits<char> >::stossc;
        void (::std::basic_streambuf<char, std::char_traits<char> >::*method_pointer_be15f1a0d14f55fa82edc4320a5a4b7e)(long) = &::std::basic_streambuf<char, std::char_traits<char> >::__safe_gbump;
        void (::std::basic_streambuf<char, std::char_traits<char> >::*method_pointer_938d42d1c3b85ff8b4b907137b25b4be)(long) = &::std::basic_streambuf<char, std::char_traits<char> >::__safe_pbump;
        boost::python::class_< class ::std::basic_streambuf<char, std::char_traits<char> >, std::shared_ptr< class ::std::basic_streambuf<char, std::char_traits<char> > >, boost::noncopyable >("_BasicStreambuf_1682b704c1ba5c66ba7aa8f2fb38a8cf", boost::python::no_init)
            .def("pubimbue", method_pointer_5f91820536bd56edaa6665a462e148db)
            .def("getloc", method_pointer_c82af173b6a65520903dd8fbcd436978)
            .def("pubsync", method_pointer_ba0b0ddea84950c183abbdace771a53b)
            .def("in_avail", method_pointer_622b3e983c985691acc059e8a2648ec0)
            .def("snextc", method_pointer_bca0a9fc504f5e239c22df6761d6998e)
            .def("sbumpc", method_pointer_d7573f6204a553e5a1f6150dad11ec95)
            .def("sgetc", method_pointer_0bf2d5f8826c5de1b429f6be002b6254)
            .def("sputbackc", method_pointer_6e8611294bf35fb8b5abb770f1d082aa)
            .def("sungetc", method_pointer_9eaa4a1508cf54e7a22544688cf0c840)
            .def("sputc", method_pointer_d4077a93f6f95daea17151c6f644a193)
            .def("stossc", method_pointer_efac272f5d1f5c34b324f50800ed6b46)
            .def("safe_gbump", method_pointer_be15f1a0d14f55fa82edc4320a5a4b7e)
            .def("safe_pbump", method_pointer_938d42d1c3b85ff8b4b907137b25b4be);
}