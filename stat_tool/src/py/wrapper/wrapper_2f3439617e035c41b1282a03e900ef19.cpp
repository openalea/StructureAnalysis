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
    ::std::string  (::std::locale::*method_pointer_8991966518a054baa9f3c4bf41f9c27d)() const = &::std::locale::name;
    bool  (::std::locale::*method_pointer_d6f72d3b43955c1bbf1e4531b69d14a9)(class ::std::locale const &) const = &::std::locale::operator==;
    bool  (::std::locale::*method_pointer_4bb8912bb35551caaf4fba37f7b53426)(class ::std::locale const &) const = &::std::locale::operator!=;
    class ::std::locale  (*method_pointer_5bf034e284795fd18ef049d2d37994ea)(class ::std::locale const &) = ::std::locale::global;
    class ::std::locale const & (*method_pointer_37beb60d04255aa09f819d3e2545c8ef)() = ::std::locale::classic;
    boost::python::class_< class ::std::locale, autowig::Held< class ::std::locale >::Type > class_2f3439617e035c41b1282a03e900ef19("Locale", "", boost::python::no_init);
    class_2f3439617e035c41b1282a03e900ef19.def(boost::python::init<  >(""));
    class_2f3439617e035c41b1282a03e900ef19.def(boost::python::init< class ::std::locale const & >(""));
    class_2f3439617e035c41b1282a03e900ef19.def(boost::python::init< class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const & >(""));
    class_2f3439617e035c41b1282a03e900ef19.def(boost::python::init< class ::std::locale const &, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const &, ::std::locale::category  >(""));
    class_2f3439617e035c41b1282a03e900ef19.def(boost::python::init< class ::std::locale const &, class ::std::locale const &, ::std::locale::category  >(""));
    class_2f3439617e035c41b1282a03e900ef19.def("name", method_pointer_8991966518a054baa9f3c4bf41f9c27d, "");
    class_2f3439617e035c41b1282a03e900ef19.def("__eq__", method_pointer_d6f72d3b43955c1bbf1e4531b69d14a9, "");
    class_2f3439617e035c41b1282a03e900ef19.def("__neq__", method_pointer_4bb8912bb35551caaf4fba37f7b53426, "");
    class_2f3439617e035c41b1282a03e900ef19.def("global", method_pointer_5bf034e284795fd18ef049d2d37994ea, "");
    class_2f3439617e035c41b1282a03e900ef19.def("classic", method_pointer_37beb60d04255aa09f819d3e2545c8ef, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_2f3439617e035c41b1282a03e900ef19.staticmethod("global");
    class_2f3439617e035c41b1282a03e900ef19.staticmethod("classic");
    class_2f3439617e035c41b1282a03e900ef19.def_readonly("none", ::std::locale::none, "");
    class_2f3439617e035c41b1282a03e900ef19.def_readonly("ctype", ::std::locale::ctype, "");
    class_2f3439617e035c41b1282a03e900ef19.def_readonly("numeric", ::std::locale::numeric, "");
    class_2f3439617e035c41b1282a03e900ef19.def_readonly("collate", ::std::locale::collate, "");
    class_2f3439617e035c41b1282a03e900ef19.def_readonly("time", ::std::locale::time, "");
    class_2f3439617e035c41b1282a03e900ef19.def_readonly("monetary", ::std::locale::monetary, "");
    class_2f3439617e035c41b1282a03e900ef19.def_readonly("messages", ::std::locale::messages, "");
    class_2f3439617e035c41b1282a03e900ef19.def_readonly("all", ::std::locale::all, "");

}