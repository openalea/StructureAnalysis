#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::ContinuousParametric const volatile * get_pointer<class ::stat_tool::ContinuousParametric const volatile >(class ::stat_tool::ContinuousParametric const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_16995d5158735f999b84ef5efdd8439b()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    void  (::stat_tool::ContinuousParametric::*method_pointer_51d3f6aaff1c56f9958e931bf18034cc)(class ::stat_tool::ContinuousParametric const &) = &::stat_tool::ContinuousParametric::copy;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::ContinuousParametric::*method_pointer_6e54558dbcae5453ad8cd8699b675a7b)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &, bool ) const = &::stat_tool::ContinuousParametric::ascii_parameter_print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::ContinuousParametric::*method_pointer_c069a296b1645684914b340b3c6efce1)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &, bool ) const = &::stat_tool::ContinuousParametric::ascii_characteristic_print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::ContinuousParametric::*method_pointer_07c437edf4365f0fb66e46d4f89672de)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &, bool , bool , class ::stat_tool::Histogram const *, class ::stat_tool::FrequencyDistribution const *) = &::stat_tool::ContinuousParametric::ascii_print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::ContinuousParametric::*method_pointer_7cb49e652f1a5820b6c078dfbb8431bc)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &) const = &::stat_tool::ContinuousParametric::spreadsheet_parameter_print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::ContinuousParametric::*method_pointer_d5cb552b875c5d35bc2b8ac79f171990)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &, bool , class ::stat_tool::Histogram const *, class ::stat_tool::FrequencyDistribution const *) = &::stat_tool::ContinuousParametric::spreadsheet_print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::ContinuousParametric::*method_pointer_94e780767b535d7990b88b31fe8ce165)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &) const = &::stat_tool::ContinuousParametric::spreadsheet_characteristic_print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::ContinuousParametric::*method_pointer_efbc93c3a0a55739b657618faeadb96b)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &) const = &::stat_tool::ContinuousParametric::plot_title_print;
    void  (::stat_tool::ContinuousParametric::*method_pointer_e07a6588990458909e46971a75350dae)(class ::stat_tool::SinglePlot &, class ::stat_tool::Histogram const *, class ::stat_tool::FrequencyDistribution const *) = &::stat_tool::ContinuousParametric::plotable_write;
    int  (::stat_tool::ContinuousParametric::*method_pointer_eb05dd7eac87515c9bc1748a0004daa7)() const = &::stat_tool::ContinuousParametric::nb_parameter_computation;
    double  (::stat_tool::ContinuousParametric::*method_pointer_2912c7228ded5fe0ad4edc7b170b41fd)(double , double ) const = &::stat_tool::ContinuousParametric::von_mises_mass_computation;
    double  (::stat_tool::ContinuousParametric::*method_pointer_34ac7fb3a6e95c50b999ebc59d6f7bbe)(double , double ) const = &::stat_tool::ContinuousParametric::mass_computation;
    void  (::stat_tool::ContinuousParametric::*method_pointer_2753680aa4bf56e888b2b587bf22bfd9)() = &::stat_tool::ContinuousParametric::von_mises_cumul_computation;
    double  (::stat_tool::ContinuousParametric::*method_pointer_744741da83c85b21b16946091ba3b49c)(class ::stat_tool::ContinuousParametric &) = &::stat_tool::ContinuousParametric::sup_norm_distance_computation;
    double  (::stat_tool::ContinuousParametric::*method_pointer_2a908d5857785e3ba04ff72f253dd84b)(class ::stat_tool::FrequencyDistribution const &, int ) const = &::stat_tool::ContinuousParametric::likelihood_computation;
    double  (::stat_tool::ContinuousParametric::*method_pointer_b440105ffde95546868447c3762fa258)() = &::stat_tool::ContinuousParametric::simulation;
    boost::python::class_< class ::stat_tool::ContinuousParametric, autowig::Held< class ::stat_tool::ContinuousParametric >::Type > class_16995d5158735f999b84ef5efdd8439b("ContinuousParametric", "Continuous parametric distribution\n\n", boost::python::no_init);
    class_16995d5158735f999b84ef5efdd8439b.def(boost::python::init< enum ::stat_tool::continuous_parametric , double , double , double , enum ::stat_tool::angle_unit  >(""));
    class_16995d5158735f999b84ef5efdd8439b.def(boost::python::init< class ::stat_tool::ContinuousParametric const & >(""));
    class_16995d5158735f999b84ef5efdd8439b.def("copy", method_pointer_51d3f6aaff1c56f9958e931bf18034cc, "");
    class_16995d5158735f999b84ef5efdd8439b.def("ascii_parameter_print", method_pointer_6e54558dbcae5453ad8cd8699b675a7b, boost::python::return_internal_reference<>(), "");
    class_16995d5158735f999b84ef5efdd8439b.def("ascii_characteristic_print", method_pointer_c069a296b1645684914b340b3c6efce1, boost::python::return_internal_reference<>(), "");
    class_16995d5158735f999b84ef5efdd8439b.def("ascii_print", method_pointer_07c437edf4365f0fb66e46d4f89672de, boost::python::return_internal_reference<>(), "");
    class_16995d5158735f999b84ef5efdd8439b.def("spreadsheet_parameter_print", method_pointer_7cb49e652f1a5820b6c078dfbb8431bc, boost::python::return_internal_reference<>(), "");
    class_16995d5158735f999b84ef5efdd8439b.def("spreadsheet_print", method_pointer_d5cb552b875c5d35bc2b8ac79f171990, boost::python::return_internal_reference<>(), "");
    class_16995d5158735f999b84ef5efdd8439b.def("spreadsheet_characteristic_print", method_pointer_94e780767b535d7990b88b31fe8ce165, boost::python::return_internal_reference<>(), "");
    class_16995d5158735f999b84ef5efdd8439b.def("plot_title_print", method_pointer_efbc93c3a0a55739b657618faeadb96b, boost::python::return_internal_reference<>(), "");
    class_16995d5158735f999b84ef5efdd8439b.def("plotable_write", method_pointer_e07a6588990458909e46971a75350dae, "");
    class_16995d5158735f999b84ef5efdd8439b.def("nb_parameter_computation", method_pointer_eb05dd7eac87515c9bc1748a0004daa7, "");
    class_16995d5158735f999b84ef5efdd8439b.def("von_mises_mass_computation", method_pointer_2912c7228ded5fe0ad4edc7b170b41fd, "");
    class_16995d5158735f999b84ef5efdd8439b.def("mass_computation", method_pointer_34ac7fb3a6e95c50b999ebc59d6f7bbe, "");
    class_16995d5158735f999b84ef5efdd8439b.def("von_mises_cumul_computation", method_pointer_2753680aa4bf56e888b2b587bf22bfd9, "");
    class_16995d5158735f999b84ef5efdd8439b.def("sup_norm_distance_computation", method_pointer_744741da83c85b21b16946091ba3b49c, "");
    class_16995d5158735f999b84ef5efdd8439b.def("likelihood_computation", method_pointer_2a908d5857785e3ba04ff72f253dd84b, "");
    class_16995d5158735f999b84ef5efdd8439b.def("simulation", method_pointer_b440105ffde95546868447c3762fa258, "");
    class_16995d5158735f999b84ef5efdd8439b.def_readwrite("ident", &::stat_tool::ContinuousParametric::ident, "");
    class_16995d5158735f999b84ef5efdd8439b.def_readwrite("min_value", &::stat_tool::ContinuousParametric::min_value, "");
    class_16995d5158735f999b84ef5efdd8439b.def_readwrite("max_value", &::stat_tool::ContinuousParametric::max_value, "");
    class_16995d5158735f999b84ef5efdd8439b.def_readwrite("unit", &::stat_tool::ContinuousParametric::unit, "");
    class_16995d5158735f999b84ef5efdd8439b.def_readwrite("slope_standard_deviation", &::stat_tool::ContinuousParametric::slope_standard_deviation, "");
    class_16995d5158735f999b84ef5efdd8439b.def_readwrite("sample_size", &::stat_tool::ContinuousParametric::sample_size, "");

}