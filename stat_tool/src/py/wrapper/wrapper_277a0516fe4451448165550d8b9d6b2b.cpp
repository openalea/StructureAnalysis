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
    void  (*method_pointer_eeca1afbc2655f9180b635ad664da8ef)(::std::char_traits< char >::char_type &, ::std::char_traits< char >::char_type const &) = ::std::char_traits< char >::assign;
    bool  (*method_pointer_fd19bd54cdbf58a1a43f8723fff837d2)(::std::char_traits< char >::char_type , ::std::char_traits< char >::char_type ) = ::std::char_traits< char >::eq;
    bool  (*method_pointer_259b431475d450bd86527ca37dcb4288)(::std::char_traits< char >::char_type , ::std::char_traits< char >::char_type ) = ::std::char_traits< char >::lt;
    ::std::char_traits< char >::int_type  (*method_pointer_25cec31515f45f5b8aabd26671025a08)(::std::char_traits< char >::int_type ) = ::std::char_traits< char >::not_eof;
    ::std::char_traits< char >::char_type  (*method_pointer_c47647219e8450888e83d417c33c6b31)(::std::char_traits< char >::int_type ) = ::std::char_traits< char >::to_char_type;
    ::std::char_traits< char >::int_type  (*method_pointer_03dec8203f335de186e878ea3a356957)(::std::char_traits< char >::char_type ) = ::std::char_traits< char >::to_int_type;
    bool  (*method_pointer_43fa1a2557245026869130b6b583bbad)(::std::char_traits< char >::int_type , ::std::char_traits< char >::int_type ) = ::std::char_traits< char >::eq_int_type;
    ::std::char_traits< char >::int_type  (*method_pointer_63db2f7e5aeb51888bb2feae4545ab13)() = ::std::char_traits< char >::eof;
    boost::python::class_< struct ::std::char_traits< char >, autowig::Held< struct ::std::char_traits< char > >::Type > class_277a0516fe4451448165550d8b9d6b2b("_CharTraits_277a0516fe4451448165550d8b9d6b2b", "", boost::python::no_init);
    class_277a0516fe4451448165550d8b9d6b2b.def("assign", method_pointer_eeca1afbc2655f9180b635ad664da8ef, "");
    class_277a0516fe4451448165550d8b9d6b2b.def("eq", method_pointer_fd19bd54cdbf58a1a43f8723fff837d2, "");
    class_277a0516fe4451448165550d8b9d6b2b.def("lt", method_pointer_259b431475d450bd86527ca37dcb4288, "");
    class_277a0516fe4451448165550d8b9d6b2b.def("not_eof", method_pointer_25cec31515f45f5b8aabd26671025a08, "");
    class_277a0516fe4451448165550d8b9d6b2b.def("to_char_type", method_pointer_c47647219e8450888e83d417c33c6b31, "");
    class_277a0516fe4451448165550d8b9d6b2b.def("to_int_type", method_pointer_03dec8203f335de186e878ea3a356957, "");
    class_277a0516fe4451448165550d8b9d6b2b.def("eq_int_type", method_pointer_43fa1a2557245026869130b6b583bbad, "");
    class_277a0516fe4451448165550d8b9d6b2b.def("eof", method_pointer_63db2f7e5aeb51888bb2feae4545ab13, "");
    class_277a0516fe4451448165550d8b9d6b2b.staticmethod("eof");
    class_277a0516fe4451448165550d8b9d6b2b.staticmethod("to_int_type");
    class_277a0516fe4451448165550d8b9d6b2b.staticmethod("not_eof");
    class_277a0516fe4451448165550d8b9d6b2b.staticmethod("eq_int_type");
    class_277a0516fe4451448165550d8b9d6b2b.staticmethod("assign");
    class_277a0516fe4451448165550d8b9d6b2b.staticmethod("lt");
    class_277a0516fe4451448165550d8b9d6b2b.staticmethod("eq");
    class_277a0516fe4451448165550d8b9d6b2b.staticmethod("to_char_type");

}