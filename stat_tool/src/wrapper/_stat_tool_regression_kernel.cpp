#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>
#include <stat_tool/regression.h>

void _stat_tool_regression_kernel()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        void (::stat_tool::RegressionKernel::*method_pointer_b9f6269292575bc5927fad1a476f13df)(class ::stat_tool::RegressionKernel const &) = &::stat_tool::RegressionKernel::copy;
        void (::stat_tool::RegressionKernel::*method_pointer_7e6ebd55b9275f0fb16ac80ccec42546)() = &::stat_tool::RegressionKernel::remove;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::RegressionKernel::*method_pointer_38793ad727d75fd99ae45c0bab92cb45)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::RegressionKernel::ascii_parameter_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::RegressionKernel::*method_pointer_e622661965515b658a029be2f8dacc2b)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::RegressionKernel::ascii_formal_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::RegressionKernel::*method_pointer_73826b86625f5048b5e6daf11a737321)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::RegressionKernel::ascii_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::RegressionKernel::*method_pointer_4f733e68331455deb67c59b49e50f996)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::RegressionKernel::spreadsheet_print;
        void (::stat_tool::RegressionKernel::*method_pointer_c522407fe3b85da5a427c7119f1d1aa0)(class ::stat_tool::SinglePlot &) const = &::stat_tool::RegressionKernel::plotable_write;
        void (::stat_tool::RegressionKernel::*method_pointer_386cc5e062f15a5d84da5c522d8c71b9)() = &::stat_tool::RegressionKernel::computation;
        double (::stat_tool::RegressionKernel::*method_pointer_f5f8fb700b86503b877dbdf9e84fe7de)() const = &::stat_tool::RegressionKernel::min_computation;
        double (::stat_tool::RegressionKernel::*method_pointer_1ff834e792c7558989512355ed6e11cb)() const = &::stat_tool::RegressionKernel::max_computation;
        boost::python::class_< class ::stat_tool::RegressionKernel, std::shared_ptr< class ::stat_tool::RegressionKernel > >("RegressionKernel", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< enum ::stat_tool::parametric_function, int, int >())
            .def(boost::python::init< class ::stat_tool::RegressionKernel const & >())
            .def("copy", method_pointer_b9f6269292575bc5927fad1a476f13df)
            .def("remove", method_pointer_7e6ebd55b9275f0fb16ac80ccec42546)
            .def("ascii_parameter_print", method_pointer_38793ad727d75fd99ae45c0bab92cb45, boost::python::return_internal_reference<>())
            .def("ascii_formal_print", method_pointer_e622661965515b658a029be2f8dacc2b, boost::python::return_internal_reference<>())
            .def("ascii_print", method_pointer_73826b86625f5048b5e6daf11a737321, boost::python::return_internal_reference<>())
            .def("spreadsheet_print", method_pointer_4f733e68331455deb67c59b49e50f996, boost::python::return_internal_reference<>())
            .def("plotable_write", method_pointer_c522407fe3b85da5a427c7119f1d1aa0)
            .def("computation", method_pointer_386cc5e062f15a5d84da5c522d8c71b9)
            .def("min_computation", method_pointer_f5f8fb700b86503b877dbdf9e84fe7de)
            .def("max_computation", method_pointer_1ff834e792c7558989512355ed6e11cb)
            .def_readwrite("ident", &::stat_tool::RegressionKernel::ident)
            .def_readwrite("min_value", &::stat_tool::RegressionKernel::min_value)
            .def_readwrite("max_value", &::stat_tool::RegressionKernel::max_value)
            .def_readwrite("regression_df", &::stat_tool::RegressionKernel::regression_df)
            .def_readwrite("residual_df", &::stat_tool::RegressionKernel::residual_df)
            .def_readwrite("nb_parameter", &::stat_tool::RegressionKernel::nb_parameter);
}