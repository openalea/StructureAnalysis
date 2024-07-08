#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::Compound const volatile * get_pointer<class ::stat_tool::Compound const volatile >(class ::stat_tool::Compound const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_e95f705e6c0c58008e5a01fae73debe6()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    class ::stat_tool::CompoundData * (::stat_tool::Compound::*method_pointer_e2e33880d5e45709a3165f902811fef7)(class ::stat_tool::StatError &) const = &::stat_tool::Compound::extract_data;
    class ::stat_tool::Compound * (*method_pointer_b2fdb2488947512f8ba0092d7c5ea96d)(class ::stat_tool::StatError &, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const, double ) = ::stat_tool::Compound::ascii_read;
    void  (::stat_tool::Compound::*method_pointer_40b9d2106f415ac190bab7e1b8715f95)(int , double , bool , bool ) = &::stat_tool::Compound::computation;
    class ::stat_tool::CompoundData * (::stat_tool::Compound::*method_pointer_b0550c167c9d5c97a464129d6e723e37)(class ::stat_tool::StatError &, int ) const = &::stat_tool::Compound::simulation;
    class ::stat_tool::CompoundData * (::stat_tool::Compound::*method_pointer_b25b77fef5d259b79275ea462cfa8cea)() const = &::stat_tool::Compound::get_compound_data;
    class ::stat_tool::DiscreteParametric * (::stat_tool::Compound::*method_pointer_9c6954c5fb705e09ad736fb70b9ed568)() const = &::stat_tool::Compound::get_sum_distribution;
    class ::stat_tool::DiscreteParametric * (::stat_tool::Compound::*method_pointer_aef6ebd1e1e252a491bf74f73217bdf4)() const = &::stat_tool::Compound::get_distribution;
    boost::python::class_< class ::stat_tool::Compound, autowig::Held< class ::stat_tool::Compound >::Type, boost::python::bases< class ::stat_tool::StatInterface, class ::stat_tool::Distribution > > class_e95f705e6c0c58008e5a01fae73debe6("Compound", "Compound distribution\n\n", boost::python::no_init);
    class_e95f705e6c0c58008e5a01fae73debe6.def(boost::python::init<  >(""));
    class_e95f705e6c0c58008e5a01fae73debe6.def(boost::python::init< class ::stat_tool::DiscreteParametric const &, class ::stat_tool::DiscreteParametric const &, double  >(""));
    class_e95f705e6c0c58008e5a01fae73debe6.def(boost::python::init< class ::stat_tool::DiscreteParametric const &, class ::stat_tool::DiscreteParametric const &, enum ::stat_tool::compound_distribution  >(""));
    class_e95f705e6c0c58008e5a01fae73debe6.def(boost::python::init< class ::stat_tool::Compound const &, bool  >(""));
    class_e95f705e6c0c58008e5a01fae73debe6.def("extract_data", method_pointer_e2e33880d5e45709a3165f902811fef7, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e95f705e6c0c58008e5a01fae73debe6.def("ascii_read", method_pointer_b2fdb2488947512f8ba0092d7c5ea96d, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e95f705e6c0c58008e5a01fae73debe6.def("computation", method_pointer_40b9d2106f415ac190bab7e1b8715f95, "");
    class_e95f705e6c0c58008e5a01fae73debe6.def("simulation", method_pointer_b0550c167c9d5c97a464129d6e723e37, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e95f705e6c0c58008e5a01fae73debe6.def("get_compound_data", method_pointer_b25b77fef5d259b79275ea462cfa8cea, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e95f705e6c0c58008e5a01fae73debe6.def("get_sum_distribution", method_pointer_9c6954c5fb705e09ad736fb70b9ed568, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e95f705e6c0c58008e5a01fae73debe6.def("get_distribution", method_pointer_aef6ebd1e1e252a491bf74f73217bdf4, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e95f705e6c0c58008e5a01fae73debe6.staticmethod("ascii_read");

    if(autowig::Held< class ::stat_tool::Compound >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::Compound >::Type, autowig::Held< class ::stat_tool::StatInterface >::Type >();
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::Compound >::Type, autowig::Held< class ::stat_tool::Distribution >::Type >();
    }

}