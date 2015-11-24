#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_vector_a0d6e5da7dec5b598e07f046e4fb632f()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        void (::std::vector<int, std::allocator<int> >::*method_pointer_4ccc81500a895368a9e8d22d68286743)(unsigned long, int const &) = &::std::vector<int, std::allocator<int> >::assign;
        unsigned long (::std::vector<int, std::allocator<int> >::*method_pointer_e07ec54b5ecf51aba9ea64933865d66f)() const = &::std::vector<int, std::allocator<int> >::size;
        unsigned long (::std::vector<int, std::allocator<int> >::*method_pointer_71bf07dd40255055b8a324e995732a46)() const = &::std::vector<int, std::allocator<int> >::max_size;
        void (::std::vector<int, std::allocator<int> >::*method_pointer_33512a4a00db51a98e4d13123a3bad29)(unsigned long) = &::std::vector<int, std::allocator<int> >::resize;
        void (::std::vector<int, std::allocator<int> >::*method_pointer_178ab65394295db6a364a19c28967f04)(unsigned long, int const &) = &::std::vector<int, std::allocator<int> >::resize;
        void (::std::vector<int, std::allocator<int> >::*method_pointer_a7b649cf3de95b3ab2863b316b860dfc)() = &::std::vector<int, std::allocator<int> >::shrink_to_fit;
        unsigned long (::std::vector<int, std::allocator<int> >::*method_pointer_2449152b3f41528eb23d2589039de7c2)() const = &::std::vector<int, std::allocator<int> >::capacity;
        bool (::std::vector<int, std::allocator<int> >::*method_pointer_6763fbcbb2be52318e0f476bfa741210)() const = &::std::vector<int, std::allocator<int> >::empty;
        void (::std::vector<int, std::allocator<int> >::*method_pointer_b4515a904899566390be0b4ff0194a91)(unsigned long) = &::std::vector<int, std::allocator<int> >::reserve;
        int const & (::std::vector<int, std::allocator<int> >::*method_pointer_96e6e91e2dc15b239c742513ae954e9d)(unsigned long) const = &::std::vector<int, std::allocator<int> >::operator[];
        int const & (::std::vector<int, std::allocator<int> >::*method_pointer_206368d540105e6b838fbd180d49ba5c)(unsigned long) const = &::std::vector<int, std::allocator<int> >::at;
        int const & (::std::vector<int, std::allocator<int> >::*method_pointer_638ebf46a23d50798c1d1afc6006bcc9)() const = &::std::vector<int, std::allocator<int> >::front;
        int const & (::std::vector<int, std::allocator<int> >::*method_pointer_37735e7cc42857a58f3b52a6b27205b7)() const = &::std::vector<int, std::allocator<int> >::back;
        void (::std::vector<int, std::allocator<int> >::*method_pointer_3b9d531e6aa554fdaa7229f51472adc2)(int const &) = &::std::vector<int, std::allocator<int> >::push_back;
        void (::std::vector<int, std::allocator<int> >::*method_pointer_8d12d434085a5f54a8c38754cb8ce092)() = &::std::vector<int, std::allocator<int> >::pop_back;
        void (::std::vector<int, std::allocator<int> >::*method_pointer_730a17c027895110a75f98866282c7c5)(class ::std::vector<int, std::allocator<int> > &) = &::std::vector<int, std::allocator<int> >::swap;
        void (::std::vector<int, std::allocator<int> >::*method_pointer_9d23ebaeda6d5d49b49b9958e15b5e51)() = &::std::vector<int, std::allocator<int> >::clear;
        boost::python::class_< class ::std::vector<int, std::allocator<int> >, std::shared_ptr< class ::std::vector<int, std::allocator<int> > > >("_Vector_a0d6e5da7dec5b598e07f046e4fb632f", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< unsigned long, class ::std::allocator<int> const & >())
            .def(boost::python::init< unsigned long, int const &, class ::std::allocator<int> const & >())
            .def(boost::python::init< class ::std::vector<int, std::allocator<int> > const & >())
            .def(boost::python::init< class ::std::vector<int, std::allocator<int> > const &, class ::std::allocator<int> const & >())
            .def(boost::python::init< class ::std::initializer_list<int>, class ::std::allocator<int> const & >())
            .def("assign", method_pointer_4ccc81500a895368a9e8d22d68286743)
            .def("size", method_pointer_e07ec54b5ecf51aba9ea64933865d66f)
            .def("max_size", method_pointer_71bf07dd40255055b8a324e995732a46)
            .def("resize", method_pointer_33512a4a00db51a98e4d13123a3bad29)
            .def("resize", method_pointer_178ab65394295db6a364a19c28967f04)
            .def("shrink_to_fit", method_pointer_a7b649cf3de95b3ab2863b316b860dfc)
            .def("capacity", method_pointer_2449152b3f41528eb23d2589039de7c2)
            .def("empty", method_pointer_6763fbcbb2be52318e0f476bfa741210)
            .def("reserve", method_pointer_b4515a904899566390be0b4ff0194a91)
            .def("__getitem__", method_pointer_96e6e91e2dc15b239c742513ae954e9d, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("at", method_pointer_206368d540105e6b838fbd180d49ba5c, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("front", method_pointer_638ebf46a23d50798c1d1afc6006bcc9, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("back", method_pointer_37735e7cc42857a58f3b52a6b27205b7, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("push_back", method_pointer_3b9d531e6aa554fdaa7229f51472adc2)
            .def("pop_back", method_pointer_8d12d434085a5f54a8c38754cb8ce092)
            .def("swap", method_pointer_730a17c027895110a75f98866282c7c5)
            .def("clear", method_pointer_9d23ebaeda6d5d49b49b9958e15b5e51);
    struct vector_a0d6e5da7dec5b598e07f046e4fb632f_from_python
    {
        vector_a0d6e5da7dec5b598e07f046e4fb632f_from_python()
        {
            boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id< class ::std::vector<int, std::allocator<int> > >());
        }

        static void* convertible(PyObject* obj_ptr)
        { return obj_ptr; }

        static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
            void* storage = ((boost::python::converter::rvalue_from_python_storage< class ::std::vector<int, std::allocator<int> > >*)data)->storage.bytes;
            new (storage) class ::std::vector<int, std::allocator<int> >();
            data->convertible = storage;
            class ::std::vector<int, std::allocator<int> >& result = *((class ::std::vector<int, std::allocator<int> >*)storage);
            unsigned int i = 0;
            for(;; i++)
            {
                boost::python::handle<> py_elem_hdl(boost::python::allow_null(PyIter_Next(obj_iter.get())));
                if(PyErr_Occurred())
                { boost::python::throw_error_already_set(); }
                if(!py_elem_hdl.get())
                { break; }
                boost::python::object py_elem_obj(py_elem_hdl);
                result.push_back((int)(boost::python::extract< int >(py_elem_obj)));
            }
        }
    };

    vector_a0d6e5da7dec5b598e07f046e4fb632f_from_python();
}