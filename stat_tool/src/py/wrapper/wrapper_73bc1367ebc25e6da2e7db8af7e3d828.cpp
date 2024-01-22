#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::Histogram const volatile * get_pointer<class ::stat_tool::Histogram const volatile >(class ::stat_tool::Histogram const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_73bc1367ebc25e6da2e7db8af7e3d828()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    void  (::stat_tool::Histogram::*method_pointer_c75330da123c5028b9eab03621ce51a6)(class ::stat_tool::Histogram const &) = &::stat_tool::Histogram::copy;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::Histogram::*method_pointer_086b952cde195856bb1d8d78ef82f59a)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &, bool ) const = &::stat_tool::Histogram::ascii_print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::Histogram::*method_pointer_c394d1a2975554e3b60b8790aeccdcce)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &) const = &::stat_tool::Histogram::spreadsheet_print;
    void  (::stat_tool::Histogram::*method_pointer_9fdb0ce985a755cfbd7d4382c1b686ab)(class ::stat_tool::SinglePlot &) const = &::stat_tool::Histogram::plotable_write;
    void  (::stat_tool::Histogram::*method_pointer_44692df90a8a558b9f0b9e588c6a1cc9)() = &::stat_tool::Histogram::max_computation;
    boost::python::class_< class ::stat_tool::Histogram, autowig::Held< class ::stat_tool::Histogram >::Type > class_73bc1367ebc25e6da2e7db8af7e3d828("Histogram", "Histogram\n\n", boost::python::no_init);
    class_73bc1367ebc25e6da2e7db8af7e3d828.def(boost::python::init< int , bool  >(""));
    class_73bc1367ebc25e6da2e7db8af7e3d828.def(boost::python::init< class ::stat_tool::FrequencyDistribution const & >(""));
    class_73bc1367ebc25e6da2e7db8af7e3d828.def(boost::python::init< class ::stat_tool::Histogram const & >(""));
    class_73bc1367ebc25e6da2e7db8af7e3d828.def("copy", method_pointer_c75330da123c5028b9eab03621ce51a6, "");
    class_73bc1367ebc25e6da2e7db8af7e3d828.def("ascii_print", method_pointer_086b952cde195856bb1d8d78ef82f59a, boost::python::return_internal_reference<>(), "");
    class_73bc1367ebc25e6da2e7db8af7e3d828.def("spreadsheet_print", method_pointer_c394d1a2975554e3b60b8790aeccdcce, boost::python::return_internal_reference<>(), "");
    class_73bc1367ebc25e6da2e7db8af7e3d828.def("plotable_write", method_pointer_9fdb0ce985a755cfbd7d4382c1b686ab, "");
    class_73bc1367ebc25e6da2e7db8af7e3d828.def("max_computation", method_pointer_44692df90a8a558b9f0b9e588c6a1cc9, "");
    class_73bc1367ebc25e6da2e7db8af7e3d828.def_readwrite("nb_element", &::stat_tool::Histogram::nb_element, "");
    class_73bc1367ebc25e6da2e7db8af7e3d828.def_readwrite("nb_bin", &::stat_tool::Histogram::nb_bin, "");
    class_73bc1367ebc25e6da2e7db8af7e3d828.def_readwrite("bin_width", &::stat_tool::Histogram::bin_width, "");
    class_73bc1367ebc25e6da2e7db8af7e3d828.def_readwrite("max", &::stat_tool::Histogram::max, "");
    class_73bc1367ebc25e6da2e7db8af7e3d828.def_readwrite("type", &::stat_tool::Histogram::type, "");
    class_73bc1367ebc25e6da2e7db8af7e3d828.def_readwrite("min_value", &::stat_tool::Histogram::min_value, "");
    class_73bc1367ebc25e6da2e7db8af7e3d828.def_readwrite("max_value", &::stat_tool::Histogram::max_value, "");

}