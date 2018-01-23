#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_4152bddb47d352a08f31cb74e24a1918(class ::std::ios_base & instance, int  param_in_0, long int param_out) { instance.iword(param_in_0) = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::ios_base const volatile * get_pointer<class ::std::ios_base const volatile >(class ::std::ios_base const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_5647113ef4105dfab0588ffcaf6c479b()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    ::std::streamsize  (::std::ios_base::*method_pointer_48b849583ebd5cffb48df9fffa5c39e4)() const = &::std::ios_base::precision;
    ::std::streamsize  (::std::ios_base::*method_pointer_7dce16c0188457d496873f7f8411c53d)(::std::streamsize ) = &::std::ios_base::precision;
    ::std::streamsize  (::std::ios_base::*method_pointer_f9a9622f5d14595aad7b0ff13aac72ad)() const = &::std::ios_base::width;
    ::std::streamsize  (::std::ios_base::*method_pointer_eafd237d52e0502c9eae565718817c85)(::std::streamsize ) = &::std::ios_base::width;
    bool  (*method_pointer_642ee2614f1a5ca5a34ee6ccea5b0b39)(bool ) = ::std::ios_base::sync_with_stdio;
    class ::std::locale  (::std::ios_base::*method_pointer_5b6536b65c305eb8b11844e126b6799e)(class ::std::locale const &) = &::std::ios_base::imbue;
    class ::std::locale  (::std::ios_base::*method_pointer_7abb1da4ca8651bfa91319be26c40f5d)() const = &::std::ios_base::getloc;
    class ::std::locale const & (::std::ios_base::*method_pointer_44e9c3dcf1c65e0c84afee91c1927f70)() const = &::std::ios_base::_M_getloc;
    int  (*method_pointer_08afc6a17c1d570481821132c7e45e80)() = ::std::ios_base::xalloc;
    long int & (::std::ios_base::*method_pointer_4152bddb47d352a08f31cb74e24a1918)(int ) = &::std::ios_base::iword;
    boost::python::class_< class ::std::ios_base, autowig::Held< class ::std::ios_base >::Type, boost::noncopyable > class_5647113ef4105dfab0588ffcaf6c479b("IosBase", "", boost::python::no_init);
    class_5647113ef4105dfab0588ffcaf6c479b.def("precision", method_pointer_48b849583ebd5cffb48df9fffa5c39e4, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("precision", method_pointer_7dce16c0188457d496873f7f8411c53d, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("width", method_pointer_f9a9622f5d14595aad7b0ff13aac72ad, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("width", method_pointer_eafd237d52e0502c9eae565718817c85, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("sync_with_stdio", method_pointer_642ee2614f1a5ca5a34ee6ccea5b0b39, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("imbue", method_pointer_5b6536b65c305eb8b11844e126b6799e, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("getloc", method_pointer_7abb1da4ca8651bfa91319be26c40f5d, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("m__getloc", method_pointer_44e9c3dcf1c65e0c84afee91c1927f70, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("xalloc", method_pointer_08afc6a17c1d570481821132c7e45e80, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("iword", method_pointer_4152bddb47d352a08f31cb74e24a1918, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("iword", autowig::method_decorator_4152bddb47d352a08f31cb74e24a1918);
    class_5647113ef4105dfab0588ffcaf6c479b.staticmethod("xalloc");
    class_5647113ef4105dfab0588ffcaf6c479b.staticmethod("sync_with_stdio");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("app", ::std::ios_base::app, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("ate", ::std::ios_base::ate, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("binary", ::std::ios_base::binary, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("in", ::std::ios_base::in, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("out", ::std::ios_base::out, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("trunc", ::std::ios_base::trunc, "");

}