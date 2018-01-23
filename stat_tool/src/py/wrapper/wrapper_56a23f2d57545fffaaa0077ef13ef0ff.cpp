#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::DiscreteMixtureData const volatile * get_pointer<class ::stat_tool::DiscreteMixtureData const volatile >(class ::stat_tool::DiscreteMixtureData const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_56a23f2d57545fffaaa0077ef13ef0ff()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    class ::stat_tool::DiscreteDistributionData * (::stat_tool::DiscreteMixtureData::*method_pointer_7d5bf00793c45f8ca8afc677d8825fd7)(class ::stat_tool::StatError &, int ) const = &::stat_tool::DiscreteMixtureData::extract;
    double  (::stat_tool::DiscreteMixtureData::*method_pointer_4bb3b235714c596f88fa6dbda61c3583)() const = &::stat_tool::DiscreteMixtureData::information_computation;
    class ::stat_tool::DiscreteMixture * (::stat_tool::DiscreteMixtureData::*method_pointer_87c6ed1036fe58e8a09e4a6b3dea1ac9)() const = &::stat_tool::DiscreteMixtureData::get_mixture;
    int  (::stat_tool::DiscreteMixtureData::*method_pointer_ccc8e0ad33b95412bd8532b25f0424ff)() const = &::stat_tool::DiscreteMixtureData::get_nb_component;
    class ::stat_tool::FrequencyDistribution * (::stat_tool::DiscreteMixtureData::*method_pointer_709feca86c75592499578eefef265ea6)() const = &::stat_tool::DiscreteMixtureData::get_weight;
    class ::stat_tool::FrequencyDistribution * (::stat_tool::DiscreteMixtureData::*method_pointer_bc26366684395f8c85296fb6178556d0)(int ) const = &::stat_tool::DiscreteMixtureData::get_component;
    boost::python::class_< class ::stat_tool::DiscreteMixtureData, autowig::Held< class ::stat_tool::DiscreteMixtureData >::Type, boost::python::bases< class ::stat_tool::StatInterface, class ::stat_tool::FrequencyDistribution > > class_56a23f2d57545fffaaa0077ef13ef0ff("DiscreteMixtureData", "Data structure corresponding to a mixture of discrete distributions\n\n", boost::python::no_init);
    class_56a23f2d57545fffaaa0077ef13ef0ff.def(boost::python::init<  >(""));
    class_56a23f2d57545fffaaa0077ef13ef0ff.def(boost::python::init< class ::stat_tool::FrequencyDistribution const &, int  >(""));
    class_56a23f2d57545fffaaa0077ef13ef0ff.def(boost::python::init< class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::DiscreteMixture const * >(""));
    class_56a23f2d57545fffaaa0077ef13ef0ff.def(boost::python::init< class ::stat_tool::DiscreteMixture const & >(""));
    class_56a23f2d57545fffaaa0077ef13ef0ff.def(boost::python::init< class ::stat_tool::DiscreteMixtureData const &, bool  >(""));
    class_56a23f2d57545fffaaa0077ef13ef0ff.def("extract", method_pointer_7d5bf00793c45f8ca8afc677d8825fd7, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_56a23f2d57545fffaaa0077ef13ef0ff.def("information_computation", method_pointer_4bb3b235714c596f88fa6dbda61c3583, "");
    class_56a23f2d57545fffaaa0077ef13ef0ff.def("get_mixture", method_pointer_87c6ed1036fe58e8a09e4a6b3dea1ac9, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_56a23f2d57545fffaaa0077ef13ef0ff.def("get_nb_component", method_pointer_ccc8e0ad33b95412bd8532b25f0424ff, "");
    class_56a23f2d57545fffaaa0077ef13ef0ff.def("get_weight", method_pointer_709feca86c75592499578eefef265ea6, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_56a23f2d57545fffaaa0077ef13ef0ff.def("get_component", method_pointer_bc26366684395f8c85296fb6178556d0, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");

    if(autowig::Held< class ::stat_tool::DiscreteMixtureData >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::DiscreteMixtureData >::Type, autowig::Held< class ::stat_tool::StatInterface >::Type >();
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::DiscreteMixtureData >::Type, autowig::Held< class ::stat_tool::FrequencyDistribution >::Type >();
    }

}