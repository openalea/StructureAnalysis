#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::locale const volatile * get_pointer<class ::std::locale const volatile >(class ::std::locale const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_2f3439617e035c41b1282a03e900ef19()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    ::std::string  (::std::locale::*method_pointer_46df17c7345b5f31a17c34d545e61a9e)() const = &::std::locale::name;
    bool  (::std::locale::*method_pointer_ac4db59f68e6542297cf43a64fbe7d6d)(class ::std::locale const &) const = &::std::locale::operator==;
    bool  (::std::locale::*method_pointer_b101c6eca41f5d04bfa1d2df4905126c)(class ::std::locale const &) const = &::std::locale::operator!=;
    class ::std::locale  (*method_pointer_4af6317d8c7955feb0876f8da6a750ab)(class ::std::locale const &) = ::std::locale::global;
    class ::std::locale const & (*method_pointer_825f9221e5085f16a66c2a7f67516c72)() = ::std::locale::classic;
    boost::python::class_< class ::std::locale, autowig::Held< class ::std::locale >::Type > class_2f3439617e035c41b1282a03e900ef19("Locale", "", boost::python::no_init);
    class_2f3439617e035c41b1282a03e900ef19.def(boost::python::init<  >(""));
    class_2f3439617e035c41b1282a03e900ef19.def(boost::python::init< class ::std::locale const & >(""));
    class_2f3439617e035c41b1282a03e900ef19.def(boost::python::init< ::std::string const & >(""));
    class_2f3439617e035c41b1282a03e900ef19.def(boost::python::init< class ::std::locale const &, ::std::string const &, ::std::locale::category  >(""));
    class_2f3439617e035c41b1282a03e900ef19.def(boost::python::init< class ::std::locale const &, class ::std::locale const &, ::std::locale::category  >(""));
    class_2f3439617e035c41b1282a03e900ef19.def("name", method_pointer_46df17c7345b5f31a17c34d545e61a9e, "");
    class_2f3439617e035c41b1282a03e900ef19.def("__eq__", method_pointer_ac4db59f68e6542297cf43a64fbe7d6d, "");
    class_2f3439617e035c41b1282a03e900ef19.def("__neq__", method_pointer_b101c6eca41f5d04bfa1d2df4905126c, "");
    class_2f3439617e035c41b1282a03e900ef19.def("global", method_pointer_4af6317d8c7955feb0876f8da6a750ab, "");
    class_2f3439617e035c41b1282a03e900ef19.def("classic", method_pointer_825f9221e5085f16a66c2a7f67516c72, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_2f3439617e035c41b1282a03e900ef19.staticmethod("global");
    class_2f3439617e035c41b1282a03e900ef19.staticmethod("classic");
    class_2f3439617e035c41b1282a03e900ef19.def_readonly("none", ::std::locale::none, "");
    class_2f3439617e035c41b1282a03e900ef19.def_readonly("collate", ::std::locale::collate, "");
    class_2f3439617e035c41b1282a03e900ef19.def_readonly("ctype", ::std::locale::ctype, "");
    class_2f3439617e035c41b1282a03e900ef19.def_readonly("monetary", ::std::locale::monetary, "");
    class_2f3439617e035c41b1282a03e900ef19.def_readonly("numeric", ::std::locale::numeric, "");
    class_2f3439617e035c41b1282a03e900ef19.def_readonly("time", ::std::locale::time, "");
    class_2f3439617e035c41b1282a03e900ef19.def_readonly("messages", ::std::locale::messages, "");
    class_2f3439617e035c41b1282a03e900ef19.def_readonly("all", ::std::locale::all, "");

}