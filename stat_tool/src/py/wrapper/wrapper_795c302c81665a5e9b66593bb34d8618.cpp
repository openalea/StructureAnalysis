#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::Convolution const volatile * get_pointer<class ::stat_tool::Convolution const volatile >(class ::stat_tool::Convolution const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_795c302c81665a5e9b66593bb34d8618()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    class ::stat_tool::DiscreteParametricModel * (::stat_tool::Convolution::*method_pointer_877bbe1644c65ee6ab38b0d65b2525e0)(class ::stat_tool::StatError &, int ) const = &::stat_tool::Convolution::extract;
    class ::stat_tool::ConvolutionData * (::stat_tool::Convolution::*method_pointer_22c6701b93ab560b966fd9b4b7acc860)(class ::stat_tool::StatError &) const = &::stat_tool::Convolution::extract_data;
    class ::stat_tool::Convolution * (*method_pointer_b4c0743fe62e55eaac78d2d7ee44569f)(class ::stat_tool::StatError &, class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > const &) = ::stat_tool::Convolution::build;
    class ::stat_tool::Convolution * (*method_pointer_31c5771ec4805559b05f69cd82804cf6)(class ::stat_tool::StatError &, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const, double ) = ::stat_tool::Convolution::ascii_read;
    class ::stat_tool::ConvolutionData * (::stat_tool::Convolution::*method_pointer_30d362d92bec5625af53132bae339399)(class ::stat_tool::StatError &, int ) const = &::stat_tool::Convolution::simulation;
    class ::stat_tool::ConvolutionData * (::stat_tool::Convolution::*method_pointer_7ef5bd984eb15ea2a7279ab780fc6836)() const = &::stat_tool::Convolution::get_convolution_data;
    int  (::stat_tool::Convolution::*method_pointer_b8c7b4200a91527cb38ad8f217089c7d)() const = &::stat_tool::Convolution::get_nb_distribution;
    class ::stat_tool::DiscreteParametric * (::stat_tool::Convolution::*method_pointer_35714d9d7c605e439c910b173beb37a6)(int ) const = &::stat_tool::Convolution::get_distribution;
    boost::python::class_< class ::stat_tool::Convolution, autowig::Held< class ::stat_tool::Convolution >::Type, boost::python::bases< class ::stat_tool::StatInterface, class ::stat_tool::Distribution > > class_795c302c81665a5e9b66593bb34d8618("Convolution", "Convolution of discrete distributions\n\n", boost::python::no_init);
    class_795c302c81665a5e9b66593bb34d8618.def(boost::python::init<  >(""));
    class_795c302c81665a5e9b66593bb34d8618.def(boost::python::init< int , class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > const & >(""));
    class_795c302c81665a5e9b66593bb34d8618.def(boost::python::init< class ::stat_tool::DiscreteParametric const &, class ::stat_tool::DiscreteParametric const & >(""));
    class_795c302c81665a5e9b66593bb34d8618.def(boost::python::init< class ::stat_tool::Convolution const &, bool  >(""));
    class_795c302c81665a5e9b66593bb34d8618.def("extract", method_pointer_877bbe1644c65ee6ab38b0d65b2525e0, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_795c302c81665a5e9b66593bb34d8618.def("extract_data", method_pointer_22c6701b93ab560b966fd9b4b7acc860, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_795c302c81665a5e9b66593bb34d8618.def("build", method_pointer_b4c0743fe62e55eaac78d2d7ee44569f, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_795c302c81665a5e9b66593bb34d8618.def("ascii_read", method_pointer_31c5771ec4805559b05f69cd82804cf6, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_795c302c81665a5e9b66593bb34d8618.def("simulation", method_pointer_30d362d92bec5625af53132bae339399, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_795c302c81665a5e9b66593bb34d8618.def("get_convolution_data", method_pointer_7ef5bd984eb15ea2a7279ab780fc6836, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_795c302c81665a5e9b66593bb34d8618.def("get_nb_distribution", method_pointer_b8c7b4200a91527cb38ad8f217089c7d, "");
    class_795c302c81665a5e9b66593bb34d8618.def("get_distribution", method_pointer_35714d9d7c605e439c910b173beb37a6, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_795c302c81665a5e9b66593bb34d8618.staticmethod("ascii_read");
    class_795c302c81665a5e9b66593bb34d8618.staticmethod("build");

    if(autowig::Held< class ::stat_tool::Convolution >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::Convolution >::Type, autowig::Held< class ::stat_tool::StatInterface >::Type >();
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::Convolution >::Type, autowig::Held< class ::stat_tool::Distribution >::Type >();
    }

}