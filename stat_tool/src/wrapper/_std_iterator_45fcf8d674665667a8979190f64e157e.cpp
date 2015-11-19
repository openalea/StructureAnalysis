#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_iterator_45fcf8d674665667a8979190f64e157e()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::class_< struct ::std::iterator<std::random_access_iterator_tag, boost::io::detail::format_item<char, std::char_traits<char>, std::allocator<char> >, long, const boost::io::detail::format_item<char, std::char_traits<char>, std::allocator<char> > *, const boost::io::detail::format_item<char, std::char_traits<char>, std::allocator<char> > &>, std::shared_ptr< struct ::std::iterator<std::random_access_iterator_tag, boost::io::detail::format_item<char, std::char_traits<char>, std::allocator<char> >, long, const boost::io::detail::format_item<char, std::char_traits<char>, std::allocator<char> > *, const boost::io::detail::format_item<char, std::char_traits<char>, std::allocator<char> > &> > >("_Iterator_45fcf8d674665667a8979190f64e157e", boost::python::no_init);
}