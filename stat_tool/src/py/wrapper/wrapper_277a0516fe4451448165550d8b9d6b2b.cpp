#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> struct ::std::char_traits< char > const volatile * get_pointer<struct ::std::char_traits< char > const volatile >(struct ::std::char_traits< char > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_277a0516fe4451448165550d8b9d6b2b()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    void  (*method_pointer_e76c133b786152b8bcf5d577c8cf59e0)(::std::char_traits< char >::char_type &, ::std::char_traits< char >::char_type const &) = ::std::char_traits< char >::assign;
//    bool  (*method_pointer_13ff000e892f5803b159660b8ca8e446)(::std::char_traits< char >::char_type const &, ::std::char_traits< char >::char_type const &) = ::std::char_traits< char >::eq;
//    bool  (*method_pointer_68da674b9adc5b81862624674da75a12)(::std::char_traits< char >::char_type const &, ::std::char_traits< char >::char_type const &) = ::std::char_traits< char >::lt;
//    ::std::char_traits< char >::char_type  (*method_pointer_12d856e62be0598da943a4031b8c7d0a)(::std::char_traits< char >::int_type const &) = ::std::char_traits< char >::to_char_type;
//    ::std::char_traits< char >::int_type  (*method_pointer_863af32f10485e3793a4327ab4dc03c7)(::std::char_traits< char >::char_type const &) = ::std::char_traits< char >::to_int_type;
//    bool  (*method_pointer_958820696dae5876b19e238cec2e0040)(::std::char_traits< char >::int_type const &, ::std::char_traits< char >::int_type const &) = ::std::char_traits< char >::eq_int_type;
    ::std::char_traits< char >::int_type  (*method_pointer_8e0245c906515b8bac9b3a032c10c3ec)() = ::std::char_traits< char >::eof;
//    ::std::char_traits< char >::int_type  (*method_pointer_5efa6ea3d7425d318a1f00b954dd4b73)(::std::char_traits< char >::int_type const &) = ::std::char_traits< char >::not_eof;
    boost::python::class_< struct ::std::char_traits< char >, autowig::Held< struct ::std::char_traits< char > >::Type > class_277a0516fe4451448165550d8b9d6b2b("_CharTraits_277a0516fe4451448165550d8b9d6b2b", "", boost::python::no_init);
    class_277a0516fe4451448165550d8b9d6b2b.def("assign", method_pointer_e76c133b786152b8bcf5d577c8cf59e0, "");
//    class_277a0516fe4451448165550d8b9d6b2b.def("eq", method_pointer_13ff000e892f5803b159660b8ca8e446, "");
//    class_277a0516fe4451448165550d8b9d6b2b.def("lt", method_pointer_68da674b9adc5b81862624674da75a12, "");
//    class_277a0516fe4451448165550d8b9d6b2b.def("to_char_type", method_pointer_12d856e62be0598da943a4031b8c7d0a, "");
//    class_277a0516fe4451448165550d8b9d6b2b.def("to_int_type", method_pointer_863af32f10485e3793a4327ab4dc03c7, "");
//    class_277a0516fe4451448165550d8b9d6b2b.def("eq_int_type", method_pointer_958820696dae5876b19e238cec2e0040, "");
    class_277a0516fe4451448165550d8b9d6b2b.def("eof", method_pointer_8e0245c906515b8bac9b3a032c10c3ec, "");
//    class_277a0516fe4451448165550d8b9d6b2b.def("not_eof", method_pointer_5efa6ea3d7425d318a1f00b954dd4b73, "");
    class_277a0516fe4451448165550d8b9d6b2b.staticmethod("eof");
//    class_277a0516fe4451448165550d8b9d6b2b.staticmethod("to_int_type");
//    class_277a0516fe4451448165550d8b9d6b2b.staticmethod("not_eof");
//    class_277a0516fe4451448165550d8b9d6b2b.staticmethod("eq_int_type");
//    class_277a0516fe4451448165550d8b9d6b2b.staticmethod("to_char_type");
//    class_277a0516fe4451448165550d8b9d6b2b.staticmethod("lt");
//    class_277a0516fe4451448165550d8b9d6b2b.staticmethod("eq");
    class_277a0516fe4451448165550d8b9d6b2b.staticmethod("assign");

}