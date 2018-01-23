#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> struct ::std::ctype_base const volatile * get_pointer<struct ::std::ctype_base const volatile >(struct ::std::ctype_base const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_96902179d8e1527ab8396789f078e437()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::class_< struct ::std::ctype_base, autowig::Held< struct ::std::ctype_base >::Type > class_96902179d8e1527ab8396789f078e437("CtypeBase", "", boost::python::no_init);
    class_96902179d8e1527ab8396789f078e437.def_readonly("upper", ::std::ctype_base::upper, "");
    class_96902179d8e1527ab8396789f078e437.def_readonly("lower", ::std::ctype_base::lower, "");
    class_96902179d8e1527ab8396789f078e437.def_readonly("alpha", ::std::ctype_base::alpha, "");
    class_96902179d8e1527ab8396789f078e437.def_readonly("digit", ::std::ctype_base::digit, "");
    class_96902179d8e1527ab8396789f078e437.def_readonly("xdigit", ::std::ctype_base::xdigit, "");
    class_96902179d8e1527ab8396789f078e437.def_readonly("space", ::std::ctype_base::space, "");
    class_96902179d8e1527ab8396789f078e437.def_readonly("print", ::std::ctype_base::print, "");
    class_96902179d8e1527ab8396789f078e437.def_readonly("graph", ::std::ctype_base::graph, "");
    class_96902179d8e1527ab8396789f078e437.def_readonly("cntrl", ::std::ctype_base::cntrl, "");
    class_96902179d8e1527ab8396789f078e437.def_readonly("punct", ::std::ctype_base::punct, "");
    class_96902179d8e1527ab8396789f078e437.def_readonly("alnum", ::std::ctype_base::alnum, "");
    class_96902179d8e1527ab8396789f078e437.def_readonly("blank", ::std::ctype_base::blank, "");

}