#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>
#include <stat_tool/discrete_mixture.h>
#include <stat_tool/markovian.h>

void _stat_tool_categorical_process()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        void (::stat_tool::CategoricalProcess::*method_pointer_9f738e654a3153949fb7de2b20354bad)(class ::stat_tool::CategoricalProcess const &) = &::stat_tool::CategoricalProcess::copy;
        void (::stat_tool::CategoricalProcess::*method_pointer_3ad9bf8558705dc4a1d8455af57c9823)() = &::stat_tool::CategoricalProcess::remove;
        bool (::stat_tool::CategoricalProcess::*method_pointer_23690c6039895e51b89862cc91bda6ef)() const = &::stat_tool::CategoricalProcess::test_hidden;
        void (::stat_tool::CategoricalProcess::*method_pointer_a1a2ac36602f5025be1030997df0ef30)(double) = &::stat_tool::CategoricalProcess::thresholding;
        int (::stat_tool::CategoricalProcess::*method_pointer_8951048518bd56cdb8bea63a4171760a)(double) const = &::stat_tool::CategoricalProcess::nb_parameter_computation;
        class ::stat_tool::Distribution * (::stat_tool::CategoricalProcess::*method_pointer_69fc5792cee357a68e2fe90f027996ac)(class ::stat_tool::Distribution *) = &::stat_tool::CategoricalProcess::mixture_computation;
        void (::stat_tool::CategoricalProcess::*method_pointer_907f57d2b6675b84b5584f60a77b5fcb)() = &::stat_tool::CategoricalProcess::init;
        boost::python::class_< class ::stat_tool::CategoricalProcess, std::shared_ptr< class ::stat_tool::CategoricalProcess > >("CategoricalProcess", boost::python::no_init)
            .def(boost::python::init< int, int, int >())
            .def(boost::python::init< class ::stat_tool::CategoricalProcess const & >())            .def("copy", method_pointer_9f738e654a3153949fb7de2b20354bad)            .def("remove", method_pointer_3ad9bf8558705dc4a1d8455af57c9823)            .def("test_hidden", method_pointer_23690c6039895e51b89862cc91bda6ef)            .def("thresholding", method_pointer_a1a2ac36602f5025be1030997df0ef30)            .def("nb_parameter_computation", method_pointer_8951048518bd56cdb8bea63a4171760a)            .def("mixture_computation", method_pointer_69fc5792cee357a68e2fe90f027996ac, boost::python::return_value_policy< boost::python::reference_existing_object >())            .def("init", method_pointer_907f57d2b6675b84b5584f60a77b5fcb)
            .def_readwrite("nb_state", &::stat_tool::CategoricalProcess::nb_state)
            .def_readwrite("nb_value", &::stat_tool::CategoricalProcess::nb_value)
            .def_readwrite("weight", &::stat_tool::CategoricalProcess::weight)
            .def_readwrite("mixture", &::stat_tool::CategoricalProcess::mixture)
            .def_readwrite("restoration_weight", &::stat_tool::CategoricalProcess::restoration_weight)
            .def_readwrite("restoration_mixture", &::stat_tool::CategoricalProcess::restoration_mixture);
}