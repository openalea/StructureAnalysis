#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::CategoricalProcess const volatile * get_pointer<class ::stat_tool::CategoricalProcess const volatile >(class ::stat_tool::CategoricalProcess const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_de7886a6095e5fde840400e992e24385()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    void  (::stat_tool::CategoricalProcess::*method_pointer_9f738e654a3153949fb7de2b20354bad)(class ::stat_tool::CategoricalProcess const &) = &::stat_tool::CategoricalProcess::copy;
    void  (::stat_tool::CategoricalProcess::*method_pointer_3ad9bf8558705dc4a1d8455af57c9823)() = &::stat_tool::CategoricalProcess::remove;
    bool  (::stat_tool::CategoricalProcess::*method_pointer_23690c6039895e51b89862cc91bda6ef)() const = &::stat_tool::CategoricalProcess::test_hidden;
    void  (::stat_tool::CategoricalProcess::*method_pointer_a1a2ac36602f5025be1030997df0ef30)(double ) = &::stat_tool::CategoricalProcess::thresholding;
    int  (::stat_tool::CategoricalProcess::*method_pointer_8951048518bd56cdb8bea63a4171760a)(double ) const = &::stat_tool::CategoricalProcess::nb_parameter_computation;
    class ::stat_tool::Distribution * (::stat_tool::CategoricalProcess::*method_pointer_69fc5792cee357a68e2fe90f027996ac)(class ::stat_tool::Distribution *) = &::stat_tool::CategoricalProcess::mixture_computation;
    void  (::stat_tool::CategoricalProcess::*method_pointer_907f57d2b6675b84b5584f60a77b5fcb)() = &::stat_tool::CategoricalProcess::init;
    boost::python::class_< class ::stat_tool::CategoricalProcess, autowig::Held< class ::stat_tool::CategoricalProcess >::Type > class_de7886a6095e5fde840400e992e24385("CategoricalProcess", "Categorical observation process\n\n", boost::python::no_init);
    class_de7886a6095e5fde840400e992e24385.def(boost::python::init< int , int , bool  >(""));
    class_de7886a6095e5fde840400e992e24385.def(boost::python::init< class ::stat_tool::CategoricalProcess const & >(""));
    class_de7886a6095e5fde840400e992e24385.def("copy", method_pointer_9f738e654a3153949fb7de2b20354bad, "");
    class_de7886a6095e5fde840400e992e24385.def("remove", method_pointer_3ad9bf8558705dc4a1d8455af57c9823, "");
    class_de7886a6095e5fde840400e992e24385.def("test_hidden", method_pointer_23690c6039895e51b89862cc91bda6ef, "");
    class_de7886a6095e5fde840400e992e24385.def("thresholding", method_pointer_a1a2ac36602f5025be1030997df0ef30, "");
    class_de7886a6095e5fde840400e992e24385.def("nb_parameter_computation", method_pointer_8951048518bd56cdb8bea63a4171760a, "");
    class_de7886a6095e5fde840400e992e24385.def("mixture_computation", method_pointer_69fc5792cee357a68e2fe90f027996ac, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_de7886a6095e5fde840400e992e24385.def("init", method_pointer_907f57d2b6675b84b5584f60a77b5fcb, "");
    class_de7886a6095e5fde840400e992e24385.def_readwrite("nb_state", &::stat_tool::CategoricalProcess::nb_state, "");
    class_de7886a6095e5fde840400e992e24385.def_readwrite("nb_value", &::stat_tool::CategoricalProcess::nb_value, "");

}