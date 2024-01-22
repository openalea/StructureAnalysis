#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::Chain const volatile * get_pointer<class ::stat_tool::Chain const volatile >(class ::stat_tool::Chain const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_d8c688ecc06b59f3b584e1e6735b6825()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    void  (::stat_tool::Chain::*method_pointer_a255fb643996510ca503e0599fd445c5)(class ::stat_tool::Chain const &) = &::stat_tool::Chain::parameter_copy;
    void  (::stat_tool::Chain::*method_pointer_3e1e9bbe45495cac8fe305ca8064a73d)(class ::stat_tool::Chain const &) = &::stat_tool::Chain::copy;
    void  (::stat_tool::Chain::*method_pointer_06eebefd01815b2cbf83a87c5b7f4349)() = &::stat_tool::Chain::remove;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::Chain::*method_pointer_6a05059fdf2d56138b0f7b3927533f71)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &, bool ) const = &::stat_tool::Chain::ascii_print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::Chain::*method_pointer_86303b0d2e9554b7896d67841d4850c9)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &) const = &::stat_tool::Chain::spreadsheet_print;
    void  (::stat_tool::Chain::*method_pointer_a7a1ab89213756b282b6ca868005bb8b)() = &::stat_tool::Chain::create_cumul;
    void  (::stat_tool::Chain::*method_pointer_4d2ce28bf1025608a64163e69ce7beb3)() = &::stat_tool::Chain::cumul_computation;
    void  (::stat_tool::Chain::*method_pointer_e4b964206af254b7aa6dc1c1df876b0e)() = &::stat_tool::Chain::remove_cumul;
    void  (::stat_tool::Chain::*method_pointer_479008557ee351fa94dd837219f04e24)() = &::stat_tool::Chain::log_computation;
    void  (::stat_tool::Chain::*method_pointer_a6a8144bc965575cb12cb1bcd1b31328)() = &::stat_tool::Chain::probability_accessibility_computation;
    bool  (::stat_tool::Chain::*method_pointer_0e422f44e4fd596992f861c8035d3ae2)() const = &::stat_tool::Chain::parallel_initial_state;
    void  (::stat_tool::Chain::*method_pointer_a1c3a26e57da572fbbe0a5fdf925e646)(double , bool ) = &::stat_tool::Chain::thresholding;
    int  (::stat_tool::Chain::*method_pointer_169d063bc4985ff49208c4e648dc7b11)(double ) const = &::stat_tool::Chain::nb_parameter_computation;
    double  (::stat_tool::Chain::*method_pointer_31c042fd128258f39ed5c70ff948be31)(class ::stat_tool::ChainData const &) const = &::stat_tool::Chain::chi2_value_computation;
    void  (::stat_tool::Chain::*method_pointer_ad88b036fac55065aae22812dc3a4b4e)(bool , double ) = &::stat_tool::Chain::init;
    double  (::stat_tool::Chain::*method_pointer_e303c8fae3345709a2ee6b3a84cd16d5)(class ::stat_tool::ChainData const &, bool ) const = &::stat_tool::Chain::likelihood_computation;
    void  (::stat_tool::Chain::*method_pointer_895ffb83a01d57c99eeb11a71fef736d)(class ::stat_tool::ChainData const &, class ::stat_tool::Test &) const = &::stat_tool::Chain::chi2_fit;
    boost::python::class_< class ::stat_tool::Chain, autowig::Held< class ::stat_tool::Chain >::Type > class_d8c688ecc06b59f3b584e1e6735b6825("Chain", "Markov chain\n\n", boost::python::no_init);
    class_d8c688ecc06b59f3b584e1e6735b6825.def(boost::python::init< enum ::stat_tool::process_type , int , bool  >(""));
    class_d8c688ecc06b59f3b584e1e6735b6825.def(boost::python::init< enum ::stat_tool::process_type , int , int , bool  >(""));
    class_d8c688ecc06b59f3b584e1e6735b6825.def(boost::python::init< class ::stat_tool::Chain const & >(""));
    class_d8c688ecc06b59f3b584e1e6735b6825.def("parameter_copy", method_pointer_a255fb643996510ca503e0599fd445c5, "");
    class_d8c688ecc06b59f3b584e1e6735b6825.def("copy", method_pointer_3e1e9bbe45495cac8fe305ca8064a73d, "");
    class_d8c688ecc06b59f3b584e1e6735b6825.def("remove", method_pointer_06eebefd01815b2cbf83a87c5b7f4349, "");
    class_d8c688ecc06b59f3b584e1e6735b6825.def("ascii_print", method_pointer_6a05059fdf2d56138b0f7b3927533f71, boost::python::return_internal_reference<>(), "");
    class_d8c688ecc06b59f3b584e1e6735b6825.def("spreadsheet_print", method_pointer_86303b0d2e9554b7896d67841d4850c9, boost::python::return_internal_reference<>(), "");
    class_d8c688ecc06b59f3b584e1e6735b6825.def("create_cumul", method_pointer_a7a1ab89213756b282b6ca868005bb8b, "");
    class_d8c688ecc06b59f3b584e1e6735b6825.def("cumul_computation", method_pointer_4d2ce28bf1025608a64163e69ce7beb3, "");
    class_d8c688ecc06b59f3b584e1e6735b6825.def("remove_cumul", method_pointer_e4b964206af254b7aa6dc1c1df876b0e, "");
    class_d8c688ecc06b59f3b584e1e6735b6825.def("log_computation", method_pointer_479008557ee351fa94dd837219f04e24, "");
    class_d8c688ecc06b59f3b584e1e6735b6825.def("probability_accessibility_computation", method_pointer_a6a8144bc965575cb12cb1bcd1b31328, "");
    class_d8c688ecc06b59f3b584e1e6735b6825.def("parallel_initial_state", method_pointer_0e422f44e4fd596992f861c8035d3ae2, "");
    class_d8c688ecc06b59f3b584e1e6735b6825.def("thresholding", method_pointer_a1c3a26e57da572fbbe0a5fdf925e646, "");
    class_d8c688ecc06b59f3b584e1e6735b6825.def("nb_parameter_computation", method_pointer_169d063bc4985ff49208c4e648dc7b11, "");
    class_d8c688ecc06b59f3b584e1e6735b6825.def("chi_2__value_computation", method_pointer_31c042fd128258f39ed5c70ff948be31, "");
    class_d8c688ecc06b59f3b584e1e6735b6825.def("init", method_pointer_ad88b036fac55065aae22812dc3a4b4e, "");
    class_d8c688ecc06b59f3b584e1e6735b6825.def("likelihood_computation", method_pointer_e303c8fae3345709a2ee6b3a84cd16d5, "");
    class_d8c688ecc06b59f3b584e1e6735b6825.def("chi_2__fit", method_pointer_895ffb83a01d57c99eeb11a71fef736d, "");
    class_d8c688ecc06b59f3b584e1e6735b6825.def_readwrite("type", &::stat_tool::Chain::type, "");
    class_d8c688ecc06b59f3b584e1e6735b6825.def_readwrite("nb_state", &::stat_tool::Chain::nb_state, "");
    class_d8c688ecc06b59f3b584e1e6735b6825.def_readwrite("nb_row", &::stat_tool::Chain::nb_row, "");
    class_d8c688ecc06b59f3b584e1e6735b6825.def_readwrite("nb_component", &::stat_tool::Chain::nb_component, "");

}