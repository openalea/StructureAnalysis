#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::Curves const volatile * get_pointer<class ::stat_tool::Curves const volatile >(class ::stat_tool::Curves const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_42d9fc9d560e519a8cee703c630b1d4b()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    void  (::stat_tool::Curves::*method_pointer_8204b1f9428953b793ca6a11c508e6c0)(class ::stat_tool::Curves const &) = &::stat_tool::Curves::copy;
    void  (::stat_tool::Curves::*method_pointer_dbd677cf5bd1512fb46974b8e337b517)(class ::stat_tool::Curves const &, int ) = &::stat_tool::Curves::smooth;
    void  (::stat_tool::Curves::*method_pointer_d941510568995943b11f44a2dfd1e379)() = &::stat_tool::Curves::remove;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::Curves::*method_pointer_b7c4c0f6e1af57669ae27ed3d493dd4c)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &, bool , class ::stat_tool::Curves const *) const = &::stat_tool::Curves::ascii_print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::Curves::*method_pointer_c67075366e29526295826ca7a5cefbd7)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &, class ::stat_tool::Curves const *) const = &::stat_tool::Curves::spreadsheet_print;
    int  (::stat_tool::Curves::*method_pointer_0e94622c5d265fca996dfce88a750f4c)() const = &::stat_tool::Curves::plot_length_computation;
    void  (::stat_tool::Curves::*method_pointer_60d8c3bcd5b85350996fb9d25c4809f9)(int , class ::stat_tool::SinglePlot &) const = &::stat_tool::Curves::plotable_write;
    void  (::stat_tool::Curves::*method_pointer_0e6ca46a684b5037acf82119fb1cd1dc)(class ::stat_tool::MultiPlot &) const = &::stat_tool::Curves::plotable_write;
    void  (::stat_tool::Curves::*method_pointer_b388de4b0d8f56ce825d0ab68adee121)(class ::stat_tool::SinglePlot &) const = &::stat_tool::Curves::plotable_frequency_write;
    int  (::stat_tool::Curves::*method_pointer_6d359e59d0cf51fda7730f88e372280a)() const = &::stat_tool::Curves::max_frequency_computation;
    int  (::stat_tool::Curves::*method_pointer_e94311af09b856d6889b53a1f222abda)() const = &::stat_tool::Curves::nb_element_computation;
    double  (::stat_tool::Curves::*method_pointer_9a007c299165594499a96997b32dda94)(int ) const = &::stat_tool::Curves::mean_computation;
    double  (::stat_tool::Curves::*method_pointer_b8f4f2ed0ffe5884933926b6a43d42fa)(int , double ) const = &::stat_tool::Curves::total_square_sum_computation;
    boost::python::class_< class ::stat_tool::Curves, autowig::Held< class ::stat_tool::Curves >::Type > class_42d9fc9d560e519a8cee703c630b1d4b("Curves", "Family of curves with frequencies\n\n", boost::python::no_init);
    class_42d9fc9d560e519a8cee703c630b1d4b.def(boost::python::init<  >(""));
    class_42d9fc9d560e519a8cee703c630b1d4b.def(boost::python::init< int , int , bool , bool , bool  >(""));
    class_42d9fc9d560e519a8cee703c630b1d4b.def(boost::python::init< class ::stat_tool::Curves const &, enum ::stat_tool::curve_transformation , int  >(""));
    class_42d9fc9d560e519a8cee703c630b1d4b.def(boost::python::init< class ::stat_tool::Distribution const & >(""));
    class_42d9fc9d560e519a8cee703c630b1d4b.def(boost::python::init< class ::stat_tool::FrequencyDistribution const & >(""));
    class_42d9fc9d560e519a8cee703c630b1d4b.def("copy", method_pointer_8204b1f9428953b793ca6a11c508e6c0, "");
    class_42d9fc9d560e519a8cee703c630b1d4b.def("smooth", method_pointer_dbd677cf5bd1512fb46974b8e337b517, "");
    class_42d9fc9d560e519a8cee703c630b1d4b.def("remove", method_pointer_d941510568995943b11f44a2dfd1e379, "");
    class_42d9fc9d560e519a8cee703c630b1d4b.def("ascii_print", method_pointer_b7c4c0f6e1af57669ae27ed3d493dd4c, boost::python::return_internal_reference<>(), "");
    class_42d9fc9d560e519a8cee703c630b1d4b.def("spreadsheet_print", method_pointer_c67075366e29526295826ca7a5cefbd7, boost::python::return_internal_reference<>(), "");
    class_42d9fc9d560e519a8cee703c630b1d4b.def("plot_length_computation", method_pointer_0e94622c5d265fca996dfce88a750f4c, "");
    class_42d9fc9d560e519a8cee703c630b1d4b.def("plotable_write", method_pointer_60d8c3bcd5b85350996fb9d25c4809f9, "");
    class_42d9fc9d560e519a8cee703c630b1d4b.def("plotable_write", method_pointer_0e6ca46a684b5037acf82119fb1cd1dc, "");
    class_42d9fc9d560e519a8cee703c630b1d4b.def("plotable_frequency_write", method_pointer_b388de4b0d8f56ce825d0ab68adee121, "");
    class_42d9fc9d560e519a8cee703c630b1d4b.def("max_frequency_computation", method_pointer_6d359e59d0cf51fda7730f88e372280a, "");
    class_42d9fc9d560e519a8cee703c630b1d4b.def("nb_element_computation", method_pointer_e94311af09b856d6889b53a1f222abda, "");
    class_42d9fc9d560e519a8cee703c630b1d4b.def("mean_computation", method_pointer_9a007c299165594499a96997b32dda94, "");
    class_42d9fc9d560e519a8cee703c630b1d4b.def("total_square_sum_computation", method_pointer_b8f4f2ed0ffe5884933926b6a43d42fa, "");
    class_42d9fc9d560e519a8cee703c630b1d4b.def_readwrite("nb_curve", &::stat_tool::Curves::nb_curve, "");
    class_42d9fc9d560e519a8cee703c630b1d4b.def_readwrite("length", &::stat_tool::Curves::length, "");
    class_42d9fc9d560e519a8cee703c630b1d4b.def_readwrite("offset", &::stat_tool::Curves::offset, "");

}