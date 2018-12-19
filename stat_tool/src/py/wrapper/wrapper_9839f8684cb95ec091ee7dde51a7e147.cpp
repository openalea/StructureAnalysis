#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::DiscreteMixture const volatile * get_pointer<class ::stat_tool::DiscreteMixture const volatile >(class ::stat_tool::DiscreteMixture const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_9839f8684cb95ec091ee7dde51a7e147()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    class ::stat_tool::DiscreteParametricModel * (::stat_tool::DiscreteMixture::*method_pointer_1485577473cd5114be525e738eeca625)(class ::stat_tool::StatError &, int ) const = &::stat_tool::DiscreteMixture::extract;
    class ::stat_tool::DiscreteMixtureData * (::stat_tool::DiscreteMixture::*method_pointer_bc96663c51af53f0b4910aa28c6a6518)(class ::stat_tool::StatError &) const = &::stat_tool::DiscreteMixture::extract_data;
    class ::stat_tool::DiscreteMixture * (*method_pointer_5866548e0e425a42b514a33827a6fb99)(class ::stat_tool::StatError &, class ::std::vector< double, class ::std::allocator< double > > const &, class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > const &) = ::stat_tool::DiscreteMixture::build;
    class ::stat_tool::DiscreteMixture * (*method_pointer_cd140b7b74fe5dee9e735d20c4318030)(class ::stat_tool::StatError &, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const, double ) = ::stat_tool::DiscreteMixture::ascii_read;
    void  (::stat_tool::DiscreteMixture::*method_pointer_280860cefdb855b48e419113a29ed6f0)(int , double , bool ) = &::stat_tool::DiscreteMixture::computation;
    double  (::stat_tool::DiscreteMixture::*method_pointer_f2bb8971d02a546e9c6e18ba9ed7ae81)(class ::stat_tool::DiscreteMixtureData const &) const = &::stat_tool::DiscreteMixture::likelihood_computation;
    class ::stat_tool::DiscreteMixtureData * (::stat_tool::DiscreteMixture::*method_pointer_6dc1d470644e5ccd97db63585a6215f8)(class ::stat_tool::StatError &, int ) const = &::stat_tool::DiscreteMixture::simulation;
    class ::stat_tool::DiscreteMixtureData * (::stat_tool::DiscreteMixture::*method_pointer_7d90346215f9587694ba87c1d315be33)() const = &::stat_tool::DiscreteMixture::get_mixture_data;
    int  (::stat_tool::DiscreteMixture::*method_pointer_ec26191b00235d48be2d6b293afbe914)() const = &::stat_tool::DiscreteMixture::get_nb_component;
    class ::stat_tool::DiscreteParametric * (::stat_tool::DiscreteMixture::*method_pointer_3b21e929cddc54d08e2c0d4e46dcc2d9)() const = &::stat_tool::DiscreteMixture::get_weight;
    class ::stat_tool::DiscreteParametric * (::stat_tool::DiscreteMixture::*method_pointer_348c66d6dde35bc08b6f5e06a873d4d1)(int ) const = &::stat_tool::DiscreteMixture::get_component;
    boost::python::class_< class ::stat_tool::DiscreteMixture, autowig::Held< class ::stat_tool::DiscreteMixture >::Type, boost::python::bases< class ::stat_tool::StatInterface, class ::stat_tool::Distribution > > class_9839f8684cb95ec091ee7dde51a7e147("DiscreteMixture", "Mixture of discrete distributions\n\n", boost::python::no_init);
    class_9839f8684cb95ec091ee7dde51a7e147.def(boost::python::init<  >(""));
    class_9839f8684cb95ec091ee7dde51a7e147.def(boost::python::init< int , class ::std::vector< double, class ::std::allocator< double > > const &, class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > const & >(""));
    class_9839f8684cb95ec091ee7dde51a7e147.def(boost::python::init< class ::stat_tool::DiscreteMixture const &, bool  >(""));
    class_9839f8684cb95ec091ee7dde51a7e147.def("extract", method_pointer_1485577473cd5114be525e738eeca625, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_9839f8684cb95ec091ee7dde51a7e147.def("extract_data", method_pointer_bc96663c51af53f0b4910aa28c6a6518, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_9839f8684cb95ec091ee7dde51a7e147.def("build", method_pointer_5866548e0e425a42b514a33827a6fb99, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_9839f8684cb95ec091ee7dde51a7e147.def("ascii_read", method_pointer_cd140b7b74fe5dee9e735d20c4318030, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_9839f8684cb95ec091ee7dde51a7e147.def("computation", method_pointer_280860cefdb855b48e419113a29ed6f0, "");
    class_9839f8684cb95ec091ee7dde51a7e147.def("likelihood_computation", method_pointer_f2bb8971d02a546e9c6e18ba9ed7ae81, "");
    class_9839f8684cb95ec091ee7dde51a7e147.def("simulation", method_pointer_6dc1d470644e5ccd97db63585a6215f8, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_9839f8684cb95ec091ee7dde51a7e147.def("get_mixture_data", method_pointer_7d90346215f9587694ba87c1d315be33, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_9839f8684cb95ec091ee7dde51a7e147.def("get_nb_component", method_pointer_ec26191b00235d48be2d6b293afbe914, "");
    class_9839f8684cb95ec091ee7dde51a7e147.def("get_weight", method_pointer_3b21e929cddc54d08e2c0d4e46dcc2d9, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_9839f8684cb95ec091ee7dde51a7e147.def("get_component", method_pointer_348c66d6dde35bc08b6f5e06a873d4d1, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_9839f8684cb95ec091ee7dde51a7e147.staticmethod("ascii_read");
    class_9839f8684cb95ec091ee7dde51a7e147.staticmethod("build");

    if(autowig::Held< class ::stat_tool::DiscreteMixture >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::DiscreteMixture >::Type, autowig::Held< class ::stat_tool::StatInterface >::Type >();
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::DiscreteMixture >::Type, autowig::Held< class ::stat_tool::Distribution >::Type >();
    }

}