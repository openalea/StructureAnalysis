#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::VectorDistance const volatile * get_pointer<class ::stat_tool::VectorDistance const volatile >(class ::stat_tool::VectorDistance const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_01172f25879c56bd9d3f985ccab12de5()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    class ::stat_tool::VectorDistance * (*method_pointer_8bbc9d8c19005cb09c1e01a07a0af378)(class ::stat_tool::StatError &, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const) = ::stat_tool::VectorDistance::ascii_read;
    void  (::stat_tool::VectorDistance::*method_pointer_c48cb58c82df5fd0845f4da056f382ed)(int , double ) const = &::stat_tool::VectorDistance::dispersion_update;
    int  (::stat_tool::VectorDistance::*method_pointer_f522981ef750530a886bbf3e448c387b)() const = &::stat_tool::VectorDistance::get_nb_variable;
    enum ::stat_tool::metric  (::stat_tool::VectorDistance::*method_pointer_8366a2e9ef805d4e98c43bcf3593ddfa)() const = &::stat_tool::VectorDistance::get_distance_type;
    enum ::stat_tool::variable_type  (::stat_tool::VectorDistance::*method_pointer_74ca394b8dd25ed5a5e109512c3edd88)(int ) const = &::stat_tool::VectorDistance::get_var_type;
    double  (::stat_tool::VectorDistance::*method_pointer_4a8076132c285fa6800f866ed97a0fd3)(int ) const = &::stat_tool::VectorDistance::get_weight;
    double  (::stat_tool::VectorDistance::*method_pointer_99c788dc725954bbba868c6f3b38ad3a)(int ) const = &::stat_tool::VectorDistance::get_dispersion;
    int  (::stat_tool::VectorDistance::*method_pointer_a1b01637456f5084b89fff4122e32f10)(int ) const = &::stat_tool::VectorDistance::get_nb_value;
    double  (::stat_tool::VectorDistance::*method_pointer_610596ad20295a5cac0e85a6352d7c9e)(int , int , int ) const = &::stat_tool::VectorDistance::get_category_distance;
    int  (::stat_tool::VectorDistance::*method_pointer_97cf13043c2f5e40babb6f83fc923071)(int ) const = &::stat_tool::VectorDistance::get_period;
    boost::python::class_< class ::stat_tool::VectorDistance, autowig::Held< class ::stat_tool::VectorDistance >::Type, boost::python::bases< class ::stat_tool::StatInterface > > class_01172f25879c56bd9d3f985ccab12de5("VectorDistance", "Parameterization of a distance between vectors with heterogeneous\nvariables\n\n", boost::python::no_init);
    class_01172f25879c56bd9d3f985ccab12de5.def(boost::python::init<  >(""));
    class_01172f25879c56bd9d3f985ccab12de5.def(boost::python::init< class ::stat_tool::VectorDistance const & >(""));
    class_01172f25879c56bd9d3f985ccab12de5.def("ascii_read", method_pointer_8bbc9d8c19005cb09c1e01a07a0af378, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_01172f25879c56bd9d3f985ccab12de5.def("dispersion_update", method_pointer_c48cb58c82df5fd0845f4da056f382ed, "");
    class_01172f25879c56bd9d3f985ccab12de5.def("get_nb_variable", method_pointer_f522981ef750530a886bbf3e448c387b, "");
    class_01172f25879c56bd9d3f985ccab12de5.def("get_distance_type", method_pointer_8366a2e9ef805d4e98c43bcf3593ddfa, "");
    class_01172f25879c56bd9d3f985ccab12de5.def("get_var_type", method_pointer_74ca394b8dd25ed5a5e109512c3edd88, "");
    class_01172f25879c56bd9d3f985ccab12de5.def("get_weight", method_pointer_4a8076132c285fa6800f866ed97a0fd3, "");
    class_01172f25879c56bd9d3f985ccab12de5.def("get_dispersion", method_pointer_99c788dc725954bbba868c6f3b38ad3a, "");
    class_01172f25879c56bd9d3f985ccab12de5.def("get_nb_value", method_pointer_a1b01637456f5084b89fff4122e32f10, "");
    class_01172f25879c56bd9d3f985ccab12de5.def("get_category_distance", method_pointer_610596ad20295a5cac0e85a6352d7c9e, "");
    class_01172f25879c56bd9d3f985ccab12de5.def("get_period", method_pointer_97cf13043c2f5e40babb6f83fc923071, "");
    class_01172f25879c56bd9d3f985ccab12de5.staticmethod("ascii_read");

    if(autowig::Held< class ::stat_tool::VectorDistance >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::VectorDistance >::Type, autowig::Held< class ::stat_tool::StatInterface >::Type >();
    }

}