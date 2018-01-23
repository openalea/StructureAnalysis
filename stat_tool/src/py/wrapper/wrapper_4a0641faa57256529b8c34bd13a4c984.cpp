#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::ChainReestimation< int > const volatile * get_pointer<class ::stat_tool::ChainReestimation< int > const volatile >(class ::stat_tool::ChainReestimation< int > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_4a0641faa57256529b8c34bd13a4c984()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    void  (::stat_tool::ChainReestimation< int >::*method_pointer_3dbcd9b08f305ec1b045caa7f09f5049)() = &::stat_tool::ChainReestimation< int >::init;
    void  (::stat_tool::ChainReestimation< int >::*method_pointer_4f76e31b5f9a5f10a230cbd30f5b0b27)(class ::stat_tool::ChainReestimation< int > const &) = &::stat_tool::ChainReestimation< int >::copy;
    void  (::stat_tool::ChainReestimation< int >::*method_pointer_1d3bc5d3ddf85900be092a22c8302c66)() = &::stat_tool::ChainReestimation< int >::remove;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::ChainReestimation< int >::*method_pointer_2fc87c6b91e75363bb690841b1960298)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &) const = &::stat_tool::ChainReestimation< int >::print;
    boost::python::class_< class ::stat_tool::ChainReestimation< int >, autowig::Held< class ::stat_tool::ChainReestimation< int > >::Type > class_4a0641faa57256529b8c34bd13a4c984("_ChainReestimation_4a0641faa57256529b8c34bd13a4c984", "", boost::python::no_init);
    class_4a0641faa57256529b8c34bd13a4c984.def(boost::python::init<  >(""));
    class_4a0641faa57256529b8c34bd13a4c984.def(boost::python::init< enum ::stat_tool::process_type , int , int , bool  >(""));
    class_4a0641faa57256529b8c34bd13a4c984.def(boost::python::init< class ::stat_tool::ChainReestimation< int > const & >(""));
    class_4a0641faa57256529b8c34bd13a4c984.def("init", method_pointer_3dbcd9b08f305ec1b045caa7f09f5049, "");
    class_4a0641faa57256529b8c34bd13a4c984.def("copy", method_pointer_4f76e31b5f9a5f10a230cbd30f5b0b27, "");
    class_4a0641faa57256529b8c34bd13a4c984.def("remove", method_pointer_1d3bc5d3ddf85900be092a22c8302c66, "");
    class_4a0641faa57256529b8c34bd13a4c984.def("print", method_pointer_2fc87c6b91e75363bb690841b1960298, boost::python::return_internal_reference<>(), "");
    class_4a0641faa57256529b8c34bd13a4c984.def_readwrite("type", &::stat_tool::ChainReestimation< int >::type, "");
    class_4a0641faa57256529b8c34bd13a4c984.def_readwrite("nb_state", &::stat_tool::ChainReestimation< int >::nb_state, "");
    class_4a0641faa57256529b8c34bd13a4c984.def_readwrite("nb_row", &::stat_tool::ChainReestimation< int >::nb_row, "");

}