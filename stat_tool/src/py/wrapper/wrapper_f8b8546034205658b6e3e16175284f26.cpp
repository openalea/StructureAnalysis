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
    bool  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_a0c0805d4b6b59ba88690cd3db015da7)() const = &::std::basic_ios< char, struct ::std::char_traits< char > >::operator!;
    bool  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_59f0426d0a165494be5b3531e0d18c9e)() const = &::std::basic_ios< char, struct ::std::char_traits< char > >::good;
    bool  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_c2f76483114053bf83e03aa784c56523)() const = &::std::basic_ios< char, struct ::std::char_traits< char > >::eof;
    bool  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_707433c165ba5532813dcef01c94f6a3)() const = &::std::basic_ios< char, struct ::std::char_traits< char > >::fail;
    bool  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_0302786b8fca56cba0fda223d7ca1019)() const = &::std::basic_ios< char, struct ::std::char_traits< char > >::bad;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > * (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_b8f82021bdb55d99ba890738ad574cd7)() const = &::std::basic_ios< char, struct ::std::char_traits< char > >::tie;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > * (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_3d881e27cd3d5e73bf07d33b5e9ef462)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > *) = &::std::basic_ios< char, struct ::std::char_traits< char > >::tie;
    class ::std::basic_streambuf< char, struct ::std::char_traits< char > > * (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_2389245b79115d4eb749dac1b5ed04c9)() const = &::std::basic_ios< char, struct ::std::char_traits< char > >::rdbuf;
    class ::std::basic_streambuf< char, struct ::std::char_traits< char > > * (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_5d783671cf0a5ed59c133d962c405171)(class ::std::basic_streambuf< char, struct ::std::char_traits< char > > *) = &::std::basic_ios< char, struct ::std::char_traits< char > >::rdbuf;
    ::std::basic_ios< char, struct ::std::char_traits< char > >::char_type  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_5a384e5d036f524d814f0b703e9bdbdc)() const = &::std::basic_ios< char, struct ::std::char_traits< char > >::fill;
    ::std::basic_ios< char, struct ::std::char_traits< char > >::char_type  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_4e3c6969231f525a9828caaba87d4fe7)(::std::basic_ios< char, struct ::std::char_traits< char > >::char_type ) = &::std::basic_ios< char, struct ::std::char_traits< char > >::fill;
    class ::std::locale  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_ee6b6443b5e25b858ac6a330c430783d)(class ::std::locale const &) = &::std::basic_ios< char, struct ::std::char_traits< char > >::imbue;
    char  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_76485232bb4956709eb82fac8224d994)(::std::basic_ios< char, struct ::std::char_traits< char > >::char_type , char ) const = &::std::basic_ios< char, struct ::std::char_traits< char > >::narrow;
    ::std::basic_ios< char, struct ::std::char_traits< char > >::char_type  (::std::basic_ios< char, ::std::char_traits< char > >::*method_pointer_852c2026975153f5ac272623962df45e)(char ) const = &::std::basic_ios< char, struct ::std::char_traits< char > >::widen;
    boost::python::class_< class ::std::basic_ios< char, struct ::std::char_traits< char > >, autowig::Held< class ::std::basic_ios< char, struct ::std::char_traits< char > > >::Type, boost::python::bases< class ::std::ios_base >, boost::noncopyable > class_f8b8546034205658b6e3e16175284f26("_BasicIos_f8b8546034205658b6e3e16175284f26", "", boost::python::no_init);
    class_f8b8546034205658b6e3e16175284f26.def(boost::python::init< class ::std::basic_streambuf< char, struct ::std::char_traits< char > > * >(""));
    class_f8b8546034205658b6e3e16175284f26.def("__not__", method_pointer_a0c0805d4b6b59ba88690cd3db015da7, "");
    class_f8b8546034205658b6e3e16175284f26.def("good", method_pointer_59f0426d0a165494be5b3531e0d18c9e, "");
    class_f8b8546034205658b6e3e16175284f26.def("eof", method_pointer_c2f76483114053bf83e03aa784c56523, "");
    class_f8b8546034205658b6e3e16175284f26.def("fail", method_pointer_707433c165ba5532813dcef01c94f6a3, "");
    class_f8b8546034205658b6e3e16175284f26.def("bad", method_pointer_0302786b8fca56cba0fda223d7ca1019, "");
    class_f8b8546034205658b6e3e16175284f26.def("tie", method_pointer_b8f82021bdb55d99ba890738ad574cd7, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f8b8546034205658b6e3e16175284f26.def("tie", method_pointer_3d881e27cd3d5e73bf07d33b5e9ef462, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f8b8546034205658b6e3e16175284f26.def("rdbuf", method_pointer_2389245b79115d4eb749dac1b5ed04c9, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f8b8546034205658b6e3e16175284f26.def("rdbuf", method_pointer_5d783671cf0a5ed59c133d962c405171, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f8b8546034205658b6e3e16175284f26.def("fill", method_pointer_5a384e5d036f524d814f0b703e9bdbdc, "");
    class_f8b8546034205658b6e3e16175284f26.def("fill", method_pointer_4e3c6969231f525a9828caaba87d4fe7, "");
    class_f8b8546034205658b6e3e16175284f26.def("imbue", method_pointer_ee6b6443b5e25b858ac6a330c430783d, "");
    class_f8b8546034205658b6e3e16175284f26.def("narrow", method_pointer_76485232bb4956709eb82fac8224d994, "");
    class_f8b8546034205658b6e3e16175284f26.def("widen", method_pointer_852c2026975153f5ac272623962df45e, "");

    if(autowig::Held< class ::std::basic_ios< char, struct ::std::char_traits< char > > >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::std::basic_ios< char, struct ::std::char_traits< char > > >::Type, autowig::Held< class ::std::ios_base >::Type >();
    }

}