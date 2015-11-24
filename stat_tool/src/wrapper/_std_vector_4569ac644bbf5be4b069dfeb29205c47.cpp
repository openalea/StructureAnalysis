#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>
#include <stat_tool/discrete_mixture.h>

void _std_vector_4569ac644bbf5be4b069dfeb29205c47()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        unsigned long (::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >::*method_pointer_93bca36f5d1553c890d2dcff1b30bb8d)() const = &::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >::size;
        unsigned long (::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >::*method_pointer_9eb729d97cac5423969d3f80919ed45d)() const = &::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >::max_size;
        void (::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >::*method_pointer_3b446d48ce1a5eed90e8bafc1d90b838)(unsigned long) = &::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >::resize;
        void (::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >::*method_pointer_7c97c8fd5c2c5a849fdc83e78d9617a7)() = &::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >::shrink_to_fit;
        unsigned long (::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >::*method_pointer_ff7745e4f3cc5116b389c07cf53c25a7)() const = &::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >::capacity;
        bool (::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >::*method_pointer_6a97351f40155d4983f2d4880d84dffd)() const = &::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >::empty;
        void (::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >::*method_pointer_c0da4672f6a258f7936fe39e55dc6b59)(unsigned long) = &::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >::reserve;
        void (::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >::*method_pointer_7bc1bd284ca05b8f99009c9709650572)() = &::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >::pop_back;
        void (::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >::*method_pointer_981386a66c1857b2b46856a7cff6dd3d)(class ::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> > &) = &::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >::swap;
        void (::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >::*method_pointer_382038a4d9ac54198930cd6c851a6ad1)() = &::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >::clear;
        boost::python::class_< class ::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >, std::shared_ptr< class ::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> > > >("_Vector_4569ac644bbf5be4b069dfeb29205c47", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> > const & >())
            .def("size", method_pointer_93bca36f5d1553c890d2dcff1b30bb8d)
            .def("max_size", method_pointer_9eb729d97cac5423969d3f80919ed45d)
            .def("resize", method_pointer_3b446d48ce1a5eed90e8bafc1d90b838)
            .def("shrink_to_fit", method_pointer_7c97c8fd5c2c5a849fdc83e78d9617a7)
            .def("capacity", method_pointer_ff7745e4f3cc5116b389c07cf53c25a7)
            .def("empty", method_pointer_6a97351f40155d4983f2d4880d84dffd)
            .def("reserve", method_pointer_c0da4672f6a258f7936fe39e55dc6b59)
            .def("pop_back", method_pointer_7bc1bd284ca05b8f99009c9709650572)
            .def("swap", method_pointer_981386a66c1857b2b46856a7cff6dd3d)
            .def("clear", method_pointer_382038a4d9ac54198930cd6c851a6ad1);
    struct vector_4569ac644bbf5be4b069dfeb29205c47_from_python
    {
        vector_4569ac644bbf5be4b069dfeb29205c47_from_python()
        {
            boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id< class ::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> > >());
        }

        static void* convertible(PyObject* obj_ptr)
        { return obj_ptr; }

        static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
            void* storage = ((boost::python::converter::rvalue_from_python_storage< class ::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> > >*)data)->storage.bytes;
            new (storage) class ::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >();
            data->convertible = storage;
            class ::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >& result = *((class ::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> >*)storage);
            unsigned int i = 0;
            for(;; i++)
            {
                boost::python::handle<> py_elem_hdl(boost::python::allow_null(PyIter_Next(obj_iter.get())));
                if(PyErr_Occurred())
                { boost::python::throw_error_already_set(); }
                if(!py_elem_hdl.get())
                { break; }
                boost::python::object py_elem_obj(py_elem_hdl);
                result.push_back((class ::stat_tool::Reestimation<double> *)(boost::python::extract< class ::stat_tool::Reestimation<double> * >(py_elem_obj)));
            }
        }
    };

    vector_4569ac644bbf5be4b069dfeb29205c47_from_python();
}