#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_pair_78364d1e814957ee9af9561e93f33e69()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        void (::std::pair<float, float>::*method_pointer_292b4a6606725594a962461df676ed20)(struct ::std::pair<float, float> &) = &::std::pair<float, float>::swap;
        boost::python::class_< struct ::std::pair<float, float>, std::shared_ptr< struct ::std::pair<float, float> > >("_Pair_78364d1e814957ee9af9561e93f33e69", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< float const &, float const & >())
            .def(boost::python::init< struct ::std::pair<float, float> const & >())
            .def("swap", method_pointer_292b4a6606725594a962461df676ed20)
            .def_readwrite("first", &::std::pair<float, float>::first)
            .def_readwrite("second", &::std::pair<float, float>::second);
}