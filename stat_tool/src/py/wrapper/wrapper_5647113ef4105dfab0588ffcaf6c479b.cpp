#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_64c8d61e29c35192899e6c44f68f830b(class ::std::ios_base & instance, int  param_in_0, long int param_out) { instance.iword(param_in_0) = param_out; }
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
    ::std::ios_base::fmtflags  (::std::ios_base::*method_pointer_b3fd7a5740895abd981922a0006d3840)() const = &::std::ios_base::flags;
    ::std::ios_base::fmtflags  (::std::ios_base::*method_pointer_f127286fc5d15f9ca0d358316fd56c6a)(::std::ios_base::fmtflags ) = &::std::ios_base::flags;
    ::std::ios_base::fmtflags  (::std::ios_base::*method_pointer_f1e19d8741cb5969a9224b901582d8ac)(::std::ios_base::fmtflags ) = &::std::ios_base::setf;
    ::std::ios_base::fmtflags  (::std::ios_base::*method_pointer_ea742a80e3aa5d33a6cf434cc7e90699)(::std::ios_base::fmtflags , ::std::ios_base::fmtflags ) = &::std::ios_base::setf;
    void  (::std::ios_base::*method_pointer_eaf94675a8b0542a95b772c09ac26f63)(::std::ios_base::fmtflags ) = &::std::ios_base::unsetf;
    ::std::streamsize  (::std::ios_base::*method_pointer_a9c0e49b7a9b56c79745590f4f7c5216)() const = &::std::ios_base::precision;
    ::std::streamsize  (::std::ios_base::*method_pointer_08b5b3b44a2a52be979d85ba932927cf)(::std::streamsize ) = &::std::ios_base::precision;
    ::std::streamsize  (::std::ios_base::*method_pointer_32a8e4c056145422b221fe45a9acd3a5)() const = &::std::ios_base::width;
    ::std::streamsize  (::std::ios_base::*method_pointer_9d1d138e1eed5486a6ef2dc82e67c1c1)(::std::streamsize ) = &::std::ios_base::width;
    class ::std::locale  (::std::ios_base::*method_pointer_5c2ca32ff79e5d719c7ec07cabcf9fc0)(class ::std::locale const &) = &::std::ios_base::imbue;
    class ::std::locale  (::std::ios_base::*method_pointer_d48e2a27590e51c0bbad409ca435ac42)() const = &::std::ios_base::getloc;
    int  (*method_pointer_0f09442e93195fe1a3ecec362dfe0587)() = ::std::ios_base::xalloc;
    long int & (::std::ios_base::*method_pointer_64c8d61e29c35192899e6c44f68f830b)(int ) = &::std::ios_base::iword;
    bool  (*method_pointer_b425de36779d59de90d603f38ba2af30)(bool ) = ::std::ios_base::sync_with_stdio;
    ::std::ios_base::iostate  (::std::ios_base::*method_pointer_c4b9769aeed058cfbdde32c60073a7f1)() const = &::std::ios_base::rdstate;
    void  (::std::ios_base::*method_pointer_7bad030138105c26950d350110a78ad2)(::std::ios_base::iostate ) = &::std::ios_base::clear;
    void  (::std::ios_base::*method_pointer_ec01e83573e6514a8aa1654ae0720a1c)(::std::ios_base::iostate ) = &::std::ios_base::setstate;
    bool  (::std::ios_base::*method_pointer_cf903d1778ce518693a2a0fd7310fad8)() const = &::std::ios_base::good;
    bool  (::std::ios_base::*method_pointer_e131459cba9b515fa2ce223cb41db0b9)() const = &::std::ios_base::eof;
    bool  (::std::ios_base::*method_pointer_996922c3581c542a9c6c6e8ffe580a7c)() const = &::std::ios_base::fail;
    bool  (::std::ios_base::*method_pointer_64b9fc18aed755848476d7ecb58e73c2)() const = &::std::ios_base::bad;
    ::std::ios_base::iostate  (::std::ios_base::*method_pointer_2f73140ee8565a7c82a57238874b6322)() const = &::std::ios_base::exceptions;
    void  (::std::ios_base::*method_pointer_1cf33f80e62d527db584c1c18be7ba52)(::std::ios_base::iostate ) = &::std::ios_base::exceptions;
    void  (::std::ios_base::*method_pointer_cbde59beae1c5519ad4f9474f830ecde)() = &::std::ios_base::__set_badbit_and_consider_rethrow;
    void  (::std::ios_base::*method_pointer_94a5ec0e584b54809c1d14d6a90bd55b)() = &::std::ios_base::__set_failbit_and_consider_rethrow;
    boost::python::class_< class ::std::ios_base, autowig::Held< class ::std::ios_base >::Type, boost::noncopyable > class_5647113ef4105dfab0588ffcaf6c479b("IosBase", "", boost::python::no_init);
    class_5647113ef4105dfab0588ffcaf6c479b.def("flags", method_pointer_b3fd7a5740895abd981922a0006d3840, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("flags", method_pointer_f127286fc5d15f9ca0d358316fd56c6a, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("setf", method_pointer_f1e19d8741cb5969a9224b901582d8ac, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("setf", method_pointer_ea742a80e3aa5d33a6cf434cc7e90699, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("unsetf", method_pointer_eaf94675a8b0542a95b772c09ac26f63, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("precision", method_pointer_a9c0e49b7a9b56c79745590f4f7c5216, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("precision", method_pointer_08b5b3b44a2a52be979d85ba932927cf, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("width", method_pointer_32a8e4c056145422b221fe45a9acd3a5, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("width", method_pointer_9d1d138e1eed5486a6ef2dc82e67c1c1, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("imbue", method_pointer_5c2ca32ff79e5d719c7ec07cabcf9fc0, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("getloc", method_pointer_d48e2a27590e51c0bbad409ca435ac42, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("xalloc", method_pointer_0f09442e93195fe1a3ecec362dfe0587, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("iword", method_pointer_64c8d61e29c35192899e6c44f68f830b, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("iword", autowig::method_decorator_64c8d61e29c35192899e6c44f68f830b);
    class_5647113ef4105dfab0588ffcaf6c479b.def("sync_with_stdio", method_pointer_b425de36779d59de90d603f38ba2af30, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("rdstate", method_pointer_c4b9769aeed058cfbdde32c60073a7f1, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("clear", method_pointer_7bad030138105c26950d350110a78ad2, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("setstate", method_pointer_ec01e83573e6514a8aa1654ae0720a1c, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("good", method_pointer_cf903d1778ce518693a2a0fd7310fad8, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("eof", method_pointer_e131459cba9b515fa2ce223cb41db0b9, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("fail", method_pointer_996922c3581c542a9c6c6e8ffe580a7c, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("bad", method_pointer_64b9fc18aed755848476d7ecb58e73c2, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("exceptions", method_pointer_2f73140ee8565a7c82a57238874b6322, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("exceptions", method_pointer_1cf33f80e62d527db584c1c18be7ba52, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("set_badbit_and_consider_rethrow", method_pointer_cbde59beae1c5519ad4f9474f830ecde, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def("set_failbit_and_consider_rethrow", method_pointer_94a5ec0e584b54809c1d14d6a90bd55b, "");
    class_5647113ef4105dfab0588ffcaf6c479b.staticmethod("xalloc");
    class_5647113ef4105dfab0588ffcaf6c479b.staticmethod("sync_with_stdio");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("boolalpha", ::std::ios_base::boolalpha, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("dec", ::std::ios_base::dec, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("fixed", ::std::ios_base::fixed, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("hex", ::std::ios_base::hex, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("internal", ::std::ios_base::internal, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("left", ::std::ios_base::left, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("oct", ::std::ios_base::oct, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("right", ::std::ios_base::right, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("scientific", ::std::ios_base::scientific, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("showbase", ::std::ios_base::showbase, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("showpoint", ::std::ios_base::showpoint, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("showpos", ::std::ios_base::showpos, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("skipws", ::std::ios_base::skipws, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("unitbuf", ::std::ios_base::unitbuf, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("uppercase", ::std::ios_base::uppercase, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("adjustfield", ::std::ios_base::adjustfield, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("basefield", ::std::ios_base::basefield, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("floatfield", ::std::ios_base::floatfield, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("badbit", ::std::ios_base::badbit, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("eofbit", ::std::ios_base::eofbit, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("failbit", ::std::ios_base::failbit, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("goodbit", ::std::ios_base::goodbit, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("app", ::std::ios_base::app, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("ate", ::std::ios_base::ate, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("binary", ::std::ios_base::binary, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("in", ::std::ios_base::in, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("out", ::std::ios_base::out, "");
    class_5647113ef4105dfab0588ffcaf6c479b.def_readonly("trunc", ::std::ios_base::trunc, "");

}