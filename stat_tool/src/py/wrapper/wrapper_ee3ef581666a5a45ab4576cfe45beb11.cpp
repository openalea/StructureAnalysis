#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::RegressionKernel const volatile * get_pointer<class ::stat_tool::RegressionKernel const volatile >(class ::stat_tool::RegressionKernel const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_ee3ef581666a5a45ab4576cfe45beb11()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    void  (::stat_tool::RegressionKernel::*method_pointer_b9f6269292575bc5927fad1a476f13df)(class ::stat_tool::RegressionKernel const &) = &::stat_tool::RegressionKernel::copy;
    void  (::stat_tool::RegressionKernel::*method_pointer_7e6ebd55b9275f0fb16ac80ccec42546)() = &::stat_tool::RegressionKernel::remove;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::RegressionKernel::*method_pointer_1f424571439a5615880e9280d5bf874f)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &) const = &::stat_tool::RegressionKernel::ascii_parameter_print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::RegressionKernel::*method_pointer_446035408a0e5771bff827367ceedecb)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &) const = &::stat_tool::RegressionKernel::ascii_formal_print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::RegressionKernel::*method_pointer_5a8f7741571c546ebc31fe0c1a4978b7)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &) const = &::stat_tool::RegressionKernel::ascii_print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::RegressionKernel::*method_pointer_1f16783679cd570fbf72708823f9761d)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &) const = &::stat_tool::RegressionKernel::spreadsheet_print;
    void  (::stat_tool::RegressionKernel::*method_pointer_c522407fe3b85da5a427c7119f1d1aa0)(class ::stat_tool::SinglePlot &) const = &::stat_tool::RegressionKernel::plotable_write;
    void  (::stat_tool::RegressionKernel::*method_pointer_386cc5e062f15a5d84da5c522d8c71b9)() = &::stat_tool::RegressionKernel::computation;
    double  (::stat_tool::RegressionKernel::*method_pointer_f5f8fb700b86503b877dbdf9e84fe7de)() const = &::stat_tool::RegressionKernel::min_computation;
    double  (::stat_tool::RegressionKernel::*method_pointer_1ff834e792c7558989512355ed6e11cb)() const = &::stat_tool::RegressionKernel::max_computation;
    boost::python::class_< class ::stat_tool::RegressionKernel, autowig::Held< class ::stat_tool::RegressionKernel >::Type > class_ee3ef581666a5a45ab4576cfe45beb11("RegressionKernel", "Regression kernel class\n\n", boost::python::no_init);
    class_ee3ef581666a5a45ab4576cfe45beb11.def(boost::python::init<  >(""));
    class_ee3ef581666a5a45ab4576cfe45beb11.def(boost::python::init< enum ::stat_tool::parametric_function , int , int  >(""));
    class_ee3ef581666a5a45ab4576cfe45beb11.def(boost::python::init< class ::stat_tool::RegressionKernel const & >(""));
    class_ee3ef581666a5a45ab4576cfe45beb11.def("copy", method_pointer_b9f6269292575bc5927fad1a476f13df, "");
    class_ee3ef581666a5a45ab4576cfe45beb11.def("remove", method_pointer_7e6ebd55b9275f0fb16ac80ccec42546, "");
    class_ee3ef581666a5a45ab4576cfe45beb11.def("ascii_parameter_print", method_pointer_1f424571439a5615880e9280d5bf874f, boost::python::return_internal_reference<>(), "");
    class_ee3ef581666a5a45ab4576cfe45beb11.def("ascii_formal_print", method_pointer_446035408a0e5771bff827367ceedecb, boost::python::return_internal_reference<>(), "");
    class_ee3ef581666a5a45ab4576cfe45beb11.def("ascii_print", method_pointer_5a8f7741571c546ebc31fe0c1a4978b7, boost::python::return_internal_reference<>(), "");
    class_ee3ef581666a5a45ab4576cfe45beb11.def("spreadsheet_print", method_pointer_1f16783679cd570fbf72708823f9761d, boost::python::return_internal_reference<>(), "");
    class_ee3ef581666a5a45ab4576cfe45beb11.def("plotable_write", method_pointer_c522407fe3b85da5a427c7119f1d1aa0, "");
    class_ee3ef581666a5a45ab4576cfe45beb11.def("computation", method_pointer_386cc5e062f15a5d84da5c522d8c71b9, "");
    class_ee3ef581666a5a45ab4576cfe45beb11.def("min_computation", method_pointer_f5f8fb700b86503b877dbdf9e84fe7de, "");
    class_ee3ef581666a5a45ab4576cfe45beb11.def("max_computation", method_pointer_1ff834e792c7558989512355ed6e11cb, "");
    class_ee3ef581666a5a45ab4576cfe45beb11.def_readwrite("ident", &::stat_tool::RegressionKernel::ident, "");
    class_ee3ef581666a5a45ab4576cfe45beb11.def_readwrite("min_value", &::stat_tool::RegressionKernel::min_value, "");
    class_ee3ef581666a5a45ab4576cfe45beb11.def_readwrite("max_value", &::stat_tool::RegressionKernel::max_value, "");
    class_ee3ef581666a5a45ab4576cfe45beb11.def_readwrite("regression_df", &::stat_tool::RegressionKernel::regression_df, "");
    class_ee3ef581666a5a45ab4576cfe45beb11.def_readwrite("residual_df", &::stat_tool::RegressionKernel::residual_df, "");
    class_ee3ef581666a5a45ab4576cfe45beb11.def_readwrite("nb_parameter", &::stat_tool::RegressionKernel::nb_parameter, "");

}