#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::basic_ios< char, struct ::std::char_traits< char > > const volatile * get_pointer<class ::std::basic_ios< char, struct ::std::char_traits< char > > const volatile >(class ::std::basic_ios< char, struct ::std::char_traits< char > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_f8b8546034205658b6e3e16175284f26()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    bool  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_c939bd3169d65779ae7c3655133ff114)() const = &::std::basic_ios< char, struct ::std::char_traits< char > >::operator!;
    ::std::ios_base::iostate  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_8a55ed3327ca5e22bf9102037ff39f6b)() const = &::std::basic_ios< char, struct ::std::char_traits< char > >::rdstate;
    void  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_aace051628e956a786dbda3dbe032018)(::std::ios_base::iostate ) = &::std::basic_ios< char, struct ::std::char_traits< char > >::clear;
    void  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_88cf2eea9f7f54329c35bfda4954f8c1)(::std::ios_base::iostate ) = &::std::basic_ios< char, struct ::std::char_traits< char > >::setstate;
    bool  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_5ba8868fa78e5e55acac642f651922d5)() const = &::std::basic_ios< char, struct ::std::char_traits< char > >::good;
    bool  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_b2ef8faca86b535a829f0285a46fb7d6)() const = &::std::basic_ios< char, struct ::std::char_traits< char > >::eof;
    bool  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_28f5c94a862f5c928437623c55385183)() const = &::std::basic_ios< char, struct ::std::char_traits< char > >::fail;
    bool  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_f2e31280c5c65b47b4c14be02a8fd976)() const = &::std::basic_ios< char, struct ::std::char_traits< char > >::bad;
    ::std::ios_base::iostate  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_53c9aed22b0f5390b578d325b23a1834)() const = &::std::basic_ios< char, struct ::std::char_traits< char > >::exceptions;
    void  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_bf8dafc283b35c06a54323456e0cc764)(::std::ios_base::iostate ) = &::std::basic_ios< char, struct ::std::char_traits< char > >::exceptions;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > * (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_98a6905451da555cabf4df07432d0f50)() const = &::std::basic_ios< char, struct ::std::char_traits< char > >::tie;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > * (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_135cbd441bf556e5a092bfc163339536)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > *) = &::std::basic_ios< char, struct ::std::char_traits< char > >::tie;
    class ::std::basic_streambuf< char, struct ::std::char_traits< char > > * (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_f8000cf4fad151899b6a7b17de9775d3)() const = &::std::basic_ios< char, struct ::std::char_traits< char > >::rdbuf;
    class ::std::basic_streambuf< char, struct ::std::char_traits< char > > * (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_872beaeb746f5b9483d9c56e83b7c653)(class ::std::basic_streambuf< char, struct ::std::char_traits< char > > *) = &::std::basic_ios< char, struct ::std::char_traits< char > >::rdbuf;
    ::std::basic_ios< char, struct ::std::char_traits< char > >::char_type  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_db54d3d01c5b542293efc1943a364ec1)() const = &::std::basic_ios< char, struct ::std::char_traits< char > >::fill;
    ::std::basic_ios< char, struct ::std::char_traits< char > >::char_type  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_df388ace61dd5a56b794deeba17cbe11)(::std::basic_ios< char, struct ::std::char_traits< char > >::char_type ) = &::std::basic_ios< char, struct ::std::char_traits< char > >::fill;
    class ::std::locale  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_f7a20494d562518b9c5d756f0fb0f150)(class ::std::locale const &) = &::std::basic_ios< char, struct ::std::char_traits< char > >::imbue;
    char  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_462930ec7e1651e2a01b7109256cf9dc)(::std::basic_ios< char, struct ::std::char_traits< char > >::char_type , char ) const = &::std::basic_ios< char, struct ::std::char_traits< char > >::narrow;
    ::std::basic_ios< char, struct ::std::char_traits< char > >::char_type  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_dbeffd7b3c7350e49e1b48601b5cc480)(char ) const = &::std::basic_ios< char, struct ::std::char_traits< char > >::widen;
    boost::python::class_< class ::std::basic_ios< char, struct ::std::char_traits< char > >, autowig::Held< class ::std::basic_ios< char, struct ::std::char_traits< char > > >::Type, boost::python::bases< class ::std::ios_base >, boost::noncopyable > class_f8b8546034205658b6e3e16175284f26("_BasicIos_f8b8546034205658b6e3e16175284f26", "", boost::python::no_init);
    class_f8b8546034205658b6e3e16175284f26.def(boost::python::init< class ::std::basic_streambuf< char, struct ::std::char_traits< char > > * >(""));
    class_f8b8546034205658b6e3e16175284f26.def("__not__", method_pointer_c939bd3169d65779ae7c3655133ff114, "");
    class_f8b8546034205658b6e3e16175284f26.def("rdstate", method_pointer_8a55ed3327ca5e22bf9102037ff39f6b, "");
    class_f8b8546034205658b6e3e16175284f26.def("clear", method_pointer_aace051628e956a786dbda3dbe032018, "");
    class_f8b8546034205658b6e3e16175284f26.def("setstate", method_pointer_88cf2eea9f7f54329c35bfda4954f8c1, "");
    class_f8b8546034205658b6e3e16175284f26.def("good", method_pointer_5ba8868fa78e5e55acac642f651922d5, "");
    class_f8b8546034205658b6e3e16175284f26.def("eof", method_pointer_b2ef8faca86b535a829f0285a46fb7d6, "");
    class_f8b8546034205658b6e3e16175284f26.def("fail", method_pointer_28f5c94a862f5c928437623c55385183, "");
    class_f8b8546034205658b6e3e16175284f26.def("bad", method_pointer_f2e31280c5c65b47b4c14be02a8fd976, "");
    class_f8b8546034205658b6e3e16175284f26.def("exceptions", method_pointer_53c9aed22b0f5390b578d325b23a1834, "");
    class_f8b8546034205658b6e3e16175284f26.def("exceptions", method_pointer_bf8dafc283b35c06a54323456e0cc764, "");
    class_f8b8546034205658b6e3e16175284f26.def("tie", method_pointer_98a6905451da555cabf4df07432d0f50, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f8b8546034205658b6e3e16175284f26.def("tie", method_pointer_135cbd441bf556e5a092bfc163339536, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f8b8546034205658b6e3e16175284f26.def("rdbuf", method_pointer_f8000cf4fad151899b6a7b17de9775d3, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f8b8546034205658b6e3e16175284f26.def("rdbuf", method_pointer_872beaeb746f5b9483d9c56e83b7c653, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f8b8546034205658b6e3e16175284f26.def("fill", method_pointer_db54d3d01c5b542293efc1943a364ec1, "");
    class_f8b8546034205658b6e3e16175284f26.def("fill", method_pointer_df388ace61dd5a56b794deeba17cbe11, "");
    class_f8b8546034205658b6e3e16175284f26.def("imbue", method_pointer_f7a20494d562518b9c5d756f0fb0f150, "");
    class_f8b8546034205658b6e3e16175284f26.def("narrow", method_pointer_462930ec7e1651e2a01b7109256cf9dc, "");
    class_f8b8546034205658b6e3e16175284f26.def("widen", method_pointer_dbeffd7b3c7350e49e1b48601b5cc480, "");

    if(autowig::Held< class ::std::basic_ios< char, struct ::std::char_traits< char > > >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::std::basic_ios< char, struct ::std::char_traits< char > > >::Type, autowig::Held< class ::std::ios_base >::Type >();
    }

}