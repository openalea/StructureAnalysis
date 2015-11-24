#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>
#include <stat_tool/discrete_mixture.h>

void _std_vector_f226c5efe1e45c11a9d350384badc7a2()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        unsigned long (::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >::*method_pointer_ee47f8f66f2953898595767b019e6f82)() const = &::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >::size;
        unsigned long (::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >::*method_pointer_e0f8c8d6203f5c21aa08d6d09664b329)() const = &::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >::max_size;
        void (::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >::*method_pointer_0f1093db4df95c2b8f5f8ee838aba8d1)(unsigned long) = &::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >::resize;
        void (::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >::*method_pointer_35a77f1e6af7524d9261842225cef7c6)() = &::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >::shrink_to_fit;
        unsigned long (::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >::*method_pointer_fdc7cfe6e6e35b53aa4b87e087d7a7ef)() const = &::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >::capacity;
        bool (::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >::*method_pointer_2ca69eb3e0eb54f2a197d891d2351237)() const = &::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >::empty;
        void (::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >::*method_pointer_759393bcc9f15a8295df7eb2d9daafab)(unsigned long) = &::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >::reserve;
        void (::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >::*method_pointer_945e0450719d532e88308b9d56e9e9a8)() = &::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >::pop_back;
        void (::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >::*method_pointer_d5e196ba76155948861537a3cfb87585)(class ::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> > &) = &::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >::swap;
        void (::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >::*method_pointer_af46c102840851ec8826c78657463a9f)() = &::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >::clear;
        boost::python::class_< class ::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >, std::shared_ptr< class ::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> > > >("_Vector_f226c5efe1e45c11a9d350384badc7a2", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> > const & >())
            .def("size", method_pointer_ee47f8f66f2953898595767b019e6f82)
            .def("max_size", method_pointer_e0f8c8d6203f5c21aa08d6d09664b329)
            .def("resize", method_pointer_0f1093db4df95c2b8f5f8ee838aba8d1)
            .def("shrink_to_fit", method_pointer_35a77f1e6af7524d9261842225cef7c6)
            .def("capacity", method_pointer_fdc7cfe6e6e35b53aa4b87e087d7a7ef)
            .def("empty", method_pointer_2ca69eb3e0eb54f2a197d891d2351237)
            .def("reserve", method_pointer_759393bcc9f15a8295df7eb2d9daafab)
            .def("pop_back", method_pointer_945e0450719d532e88308b9d56e9e9a8)
            .def("swap", method_pointer_d5e196ba76155948861537a3cfb87585)
            .def("clear", method_pointer_af46c102840851ec8826c78657463a9f);
    struct vector_f226c5efe1e45c11a9d350384badc7a2_from_python
    {
        vector_f226c5efe1e45c11a9d350384badc7a2_from_python()
        {
            boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id< class ::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> > >());
        }

        static void* convertible(PyObject* obj_ptr)
        { return obj_ptr; }

        static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
            void* storage = ((boost::python::converter::rvalue_from_python_storage< class ::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> > >*)data)->storage.bytes;
            new (storage) class ::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >();
            data->convertible = storage;
            class ::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >& result = *((class ::std::vector<stat_tool::Reestimation<int> *, std::allocator<stat_tool::Reestimation<int> *> >*)storage);
            unsigned int i = 0;
            for(;; i++)
            {
                boost::python::handle<> py_elem_hdl(boost::python::allow_null(PyIter_Next(obj_iter.get())));
                if(PyErr_Occurred())
                { boost::python::throw_error_already_set(); }
                if(!py_elem_hdl.get())
                { break; }
                boost::python::object py_elem_obj(py_elem_hdl);
                result.push_back((class ::stat_tool::Reestimation<int> *)(boost::python::extract< class ::stat_tool::Reestimation<int> * >(py_elem_obj)));
            }
        }
    };

    vector_f226c5efe1e45c11a9d350384badc7a2_from_python();
}