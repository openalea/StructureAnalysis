#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>
#include <stat_tool/discrete_mixture.h>

void _stat_tool_discrete_mixture()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        class ::stat_tool::DiscreteParametricModel * (::stat_tool::DiscreteMixture::*method_pointer_1485577473cd5114be525e738eeca625)(class ::stat_tool::StatError &, int) const = &::stat_tool::DiscreteMixture::extract;
        class ::stat_tool::DiscreteMixtureData * (::stat_tool::DiscreteMixture::*method_pointer_bc96663c51af53f0b4910aa28c6a6518)(class ::stat_tool::StatError &) const = &::stat_tool::DiscreteMixture::extract_data;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::DiscreteMixture::*method_pointer_1e3a103ec28d5f7cafedb0ec49f93683)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::DiscreteMixture::line_write;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::DiscreteMixture::*method_pointer_439404d771385f1399eff99ce154ef2b)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool) const = &::stat_tool::DiscreteMixture::ascii_write;
        class ::stat_tool::MultiPlotSet * (::stat_tool::DiscreteMixture::*method_pointer_9c4f8dd882725a3790263fb7a23fa1ee)() const = &::stat_tool::DiscreteMixture::get_plotable;
        void (::stat_tool::DiscreteMixture::*method_pointer_280860cefdb855b48e419113a29ed6f0)(int, double, bool) = &::stat_tool::DiscreteMixture::computation;
        double (::stat_tool::DiscreteMixture::*method_pointer_f2bb8971d02a546e9c6e18ba9ed7ae81)(class ::stat_tool::DiscreteMixtureData const &) const = &::stat_tool::DiscreteMixture::likelihood_computation;
        class ::stat_tool::DiscreteMixtureData * (::stat_tool::DiscreteMixture::*method_pointer_6dc1d470644e5ccd97db63585a6215f8)(class ::stat_tool::StatError &, int) const = &::stat_tool::DiscreteMixture::simulation;
        class ::stat_tool::DiscreteMixtureData * (::stat_tool::DiscreteMixture::*method_pointer_7d90346215f9587694ba87c1d315be33)() const = &::stat_tool::DiscreteMixture::get_mixture_data;
        int (::stat_tool::DiscreteMixture::*method_pointer_ec26191b00235d48be2d6b293afbe914)() const = &::stat_tool::DiscreteMixture::get_nb_component;
        class ::stat_tool::DiscreteParametric * (::stat_tool::DiscreteMixture::*method_pointer_3b21e929cddc54d08e2c0d4e46dcc2d9)() const = &::stat_tool::DiscreteMixture::get_weight;
        class ::stat_tool::DiscreteParametric * (::stat_tool::DiscreteMixture::*method_pointer_348c66d6dde35bc08b6f5e06a873d4d1)(int) const = &::stat_tool::DiscreteMixture::get_component;
        boost::python::class_< class ::stat_tool::DiscreteMixture, std::shared_ptr< class ::stat_tool::DiscreteMixture >, boost::python::bases< class ::stat_tool::StatInterface, class ::stat_tool::Distribution > >("DiscreteMixture", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::stat_tool::DiscreteMixture const &, bool >())            .def("extract", method_pointer_1485577473cd5114be525e738eeca625, boost::python::return_value_policy< boost::python::reference_existing_object >())            .def("extract_data", method_pointer_bc96663c51af53f0b4910aa28c6a6518, boost::python::return_value_policy< boost::python::reference_existing_object >())            .def("line_write", method_pointer_1e3a103ec28d5f7cafedb0ec49f93683, boost::python::return_internal_reference<>())            .def("ascii_write", method_pointer_439404d771385f1399eff99ce154ef2b, boost::python::return_internal_reference<>())            .def("get_plotable", method_pointer_9c4f8dd882725a3790263fb7a23fa1ee, boost::python::return_value_policy< boost::python::reference_existing_object >())            .def("computation", method_pointer_280860cefdb855b48e419113a29ed6f0)            .def("likelihood_computation", method_pointer_f2bb8971d02a546e9c6e18ba9ed7ae81)            .def("simulation", method_pointer_6dc1d470644e5ccd97db63585a6215f8, boost::python::return_value_policy< boost::python::reference_existing_object >())            .def("get_mixture_data", method_pointer_7d90346215f9587694ba87c1d315be33, boost::python::return_value_policy< boost::python::reference_existing_object >())            .def("get_nb_component", method_pointer_ec26191b00235d48be2d6b293afbe914)            .def("get_weight", method_pointer_3b21e929cddc54d08e2c0d4e46dcc2d9, boost::python::return_value_policy< boost::python::reference_existing_object >())            .def("get_component", method_pointer_348c66d6dde35bc08b6f5e06a873d4d1, boost::python::return_value_policy< boost::python::reference_existing_object >());
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::DiscreteMixture >, std::shared_ptr< class ::stat_tool::StatInterface > >();
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::DiscreteMixture >, std::shared_ptr< class ::stat_tool::Distribution > >();
}