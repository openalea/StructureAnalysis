#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::ContinuousParametricProcess const volatile * get_pointer<class ::stat_tool::ContinuousParametricProcess const volatile >(class ::stat_tool::ContinuousParametricProcess const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_9ba34b38b8a75c968e25c23f75118106()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    void  (::stat_tool::ContinuousParametricProcess::*method_pointer_a51fa7e43cd559d48e6193d1985ac289)(class ::stat_tool::ContinuousParametricProcess const &) = &::stat_tool::ContinuousParametricProcess::copy;
    void  (::stat_tool::ContinuousParametricProcess::*method_pointer_c0ed3022bc4c5e1ba45f2928b3b5af9b)() = &::stat_tool::ContinuousParametricProcess::remove;
    int  (::stat_tool::ContinuousParametricProcess::*method_pointer_d474df83ba195ee181fd24c1a5b56fc2)() const = &::stat_tool::ContinuousParametricProcess::nb_parameter_computation;
    double  (::stat_tool::ContinuousParametricProcess::*method_pointer_0a554efc0f6c5d83a995f8fbf54cad2d)(class ::stat_tool::Distribution *) const = &::stat_tool::ContinuousParametricProcess::mean_computation;
    double  (::stat_tool::ContinuousParametricProcess::*method_pointer_cebc3dc50e265dddb861748541b68e5b)(class ::stat_tool::Distribution *, double ) const = &::stat_tool::ContinuousParametricProcess::variance_computation;
    void  (::stat_tool::ContinuousParametricProcess::*method_pointer_e4038a834f085556ba7622988e069150)(enum ::stat_tool::angle_unit ) = &::stat_tool::ContinuousParametricProcess::select_unit;
    void  (::stat_tool::ContinuousParametricProcess::*method_pointer_f1e34b3cae8354bf88b0b69db8db7c67)(enum ::stat_tool::continuous_parametric , double , double , double , double ) = &::stat_tool::ContinuousParametricProcess::init;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::ContinuousParametricProcess::*method_pointer_6b33212305505b149ccad021a22673f9)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &) = &::stat_tool::ContinuousParametricProcess::interval_computation;
    boost::python::class_< class ::stat_tool::ContinuousParametricProcess, autowig::Held< class ::stat_tool::ContinuousParametricProcess >::Type > class_9ba34b38b8a75c968e25c23f75118106("ContinuousParametricProcess", "Continuous parametric observation process\n\n", boost::python::no_init);
    class_9ba34b38b8a75c968e25c23f75118106.def(boost::python::init< int  >(""));
    class_9ba34b38b8a75c968e25c23f75118106.def(boost::python::init< class ::stat_tool::ContinuousParametricProcess const & >(""));
    class_9ba34b38b8a75c968e25c23f75118106.def("copy", method_pointer_a51fa7e43cd559d48e6193d1985ac289, "");
    class_9ba34b38b8a75c968e25c23f75118106.def("remove", method_pointer_c0ed3022bc4c5e1ba45f2928b3b5af9b, "");
    class_9ba34b38b8a75c968e25c23f75118106.def("nb_parameter_computation", method_pointer_d474df83ba195ee181fd24c1a5b56fc2, "");
    class_9ba34b38b8a75c968e25c23f75118106.def("mean_computation", method_pointer_0a554efc0f6c5d83a995f8fbf54cad2d, "");
    class_9ba34b38b8a75c968e25c23f75118106.def("variance_computation", method_pointer_cebc3dc50e265dddb861748541b68e5b, "");
    class_9ba34b38b8a75c968e25c23f75118106.def("select_unit", method_pointer_e4038a834f085556ba7622988e069150, "");
    class_9ba34b38b8a75c968e25c23f75118106.def("init", method_pointer_f1e34b3cae8354bf88b0b69db8db7c67, "");
    class_9ba34b38b8a75c968e25c23f75118106.def("interval_computation", method_pointer_6b33212305505b149ccad021a22673f9, boost::python::return_internal_reference<>(), "");
    class_9ba34b38b8a75c968e25c23f75118106.def_readwrite("nb_state", &::stat_tool::ContinuousParametricProcess::nb_state, "");
    class_9ba34b38b8a75c968e25c23f75118106.def_readwrite("ident", &::stat_tool::ContinuousParametricProcess::ident, "");
    class_9ba34b38b8a75c968e25c23f75118106.def_readwrite("tied_location", &::stat_tool::ContinuousParametricProcess::tied_location, "");
    class_9ba34b38b8a75c968e25c23f75118106.def_readwrite("tied_dispersion", &::stat_tool::ContinuousParametricProcess::tied_dispersion, "");
    class_9ba34b38b8a75c968e25c23f75118106.def_readwrite("offset", &::stat_tool::ContinuousParametricProcess::offset, "");
    class_9ba34b38b8a75c968e25c23f75118106.def_readwrite("unit", &::stat_tool::ContinuousParametricProcess::unit, "");

}