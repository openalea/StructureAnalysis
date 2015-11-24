#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_basic_filebuf_627f248b8746525980451ecc501709fb()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        bool (::std::basic_filebuf<char, std::char_traits<char> >::*method_pointer_c7c76158e1f75b3e9c5199c4eee96a67)() const = &::std::basic_filebuf<char, std::char_traits<char> >::is_open;
        class ::std::basic_filebuf<char, std::char_traits<char> > * (::std::basic_filebuf<char, std::char_traits<char> >::*method_pointer_2e0098da8bc353ce9402ebb05033fbd0)(class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, enum ::std::_Ios_Openmode) = &::std::basic_filebuf<char, std::char_traits<char> >::open;
        class ::std::basic_filebuf<char, std::char_traits<char> > * (::std::basic_filebuf<char, std::char_traits<char> >::*method_pointer_aeda7330bd1e5e11b20519759ae1ad72)() = &::std::basic_filebuf<char, std::char_traits<char> >::close;
        boost::python::class_< class ::std::basic_filebuf<char, std::char_traits<char> >, std::shared_ptr< class ::std::basic_filebuf<char, std::char_traits<char> > >, boost::python::bases< class ::std::basic_streambuf<char, std::char_traits<char> > >, boost::noncopyable >("_BasicFilebuf_627f248b8746525980451ecc501709fb", boost::python::no_init)
            .def(boost::python::init<  >())
            .def("is_open", method_pointer_c7c76158e1f75b3e9c5199c4eee96a67)
            .def("open", method_pointer_2e0098da8bc353ce9402ebb05033fbd0, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("close", method_pointer_aeda7330bd1e5e11b20519759ae1ad72, boost::python::return_value_policy< boost::python::reference_existing_object >());
        boost::python::implicitly_convertible< std::shared_ptr< class ::std::basic_filebuf<char, std::char_traits<char> > >, std::shared_ptr< class ::std::basic_streambuf<char, std::char_traits<char> > > >();
}