#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_ios_base()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        enum ::std::_Ios_Fmtflags (::std::ios_base::*method_pointer_da8f39db41115ef9a747a2936d8defee)() const = &::std::ios_base::flags;
        enum ::std::_Ios_Fmtflags (::std::ios_base::*method_pointer_2dea9cc757f65b719eaebb9610bc715c)(enum ::std::_Ios_Fmtflags) = &::std::ios_base::flags;
        enum ::std::_Ios_Fmtflags (::std::ios_base::*method_pointer_b7032ce59fcc58da92d0570cfa5813b0)(enum ::std::_Ios_Fmtflags) = &::std::ios_base::setf;
        enum ::std::_Ios_Fmtflags (::std::ios_base::*method_pointer_b43eb7f743dd5d72aab1fa5561ba401a)(enum ::std::_Ios_Fmtflags, enum ::std::_Ios_Fmtflags) = &::std::ios_base::setf;
        void (::std::ios_base::*method_pointer_206cd8d350ec5297bcff4dfb8ad141a0)(enum ::std::_Ios_Fmtflags) = &::std::ios_base::unsetf;
        long (::std::ios_base::*method_pointer_48b849583ebd5cffb48df9fffa5c39e4)() const = &::std::ios_base::precision;
        long (::std::ios_base::*method_pointer_7dce16c0188457d496873f7f8411c53d)(long) = &::std::ios_base::precision;
        long (::std::ios_base::*method_pointer_f9a9622f5d14595aad7b0ff13aac72ad)() const = &::std::ios_base::width;
        long (::std::ios_base::*method_pointer_eafd237d52e0502c9eae565718817c85)(long) = &::std::ios_base::width;
        bool (*method_pointer_642ee2614f1a5ca5a34ee6ccea5b0b39)(bool) = ::std::ios_base::sync_with_stdio;
        class ::std::locale (::std::ios_base::*method_pointer_5b6536b65c305eb8b11844e126b6799e)(class ::std::locale const &) = &::std::ios_base::imbue;
        class ::std::locale (::std::ios_base::*method_pointer_7abb1da4ca8651bfa91319be26c40f5d)() const = &::std::ios_base::getloc;
        class ::std::locale const & (::std::ios_base::*method_pointer_44e9c3dcf1c65e0c84afee91c1927f70)() const = &::std::ios_base::_M_getloc;
        int (*method_pointer_08afc6a17c1d570481821132c7e45e80)() = ::std::ios_base::xalloc;
        boost::python::class_< class ::std::ios_base, std::shared_ptr< class ::std::ios_base >, boost::noncopyable >("IosBase", boost::python::no_init)
            .def("flags", method_pointer_da8f39db41115ef9a747a2936d8defee)
            .def("flags", method_pointer_2dea9cc757f65b719eaebb9610bc715c)
            .def("setf", method_pointer_b7032ce59fcc58da92d0570cfa5813b0)
            .def("setf", method_pointer_b43eb7f743dd5d72aab1fa5561ba401a)
            .def("unsetf", method_pointer_206cd8d350ec5297bcff4dfb8ad141a0)
            .def("precision", method_pointer_48b849583ebd5cffb48df9fffa5c39e4)
            .def("precision", method_pointer_7dce16c0188457d496873f7f8411c53d)
            .def("width", method_pointer_f9a9622f5d14595aad7b0ff13aac72ad)
            .def("width", method_pointer_eafd237d52e0502c9eae565718817c85)
            .def("sync_with_stdio", method_pointer_642ee2614f1a5ca5a34ee6ccea5b0b39)
            .def("imbue", method_pointer_5b6536b65c305eb8b11844e126b6799e)
            .def("getloc", method_pointer_7abb1da4ca8651bfa91319be26c40f5d)
            .def("m__getloc", method_pointer_44e9c3dcf1c65e0c84afee91c1927f70, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("xalloc", method_pointer_08afc6a17c1d570481821132c7e45e80)
            .staticmethod("xalloc")
            .staticmethod("sync_with_stdio");
}