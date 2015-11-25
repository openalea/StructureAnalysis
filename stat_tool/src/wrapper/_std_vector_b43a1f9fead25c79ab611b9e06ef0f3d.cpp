#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_vector_b43a1f9fead25c79ab611b9e06ef0f3d()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        void (::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::*method_pointer_abd86f72c4ac5033b84faac5503beaa8)(unsigned long, enum ::stat_tool::discrete_parametric const &) = &::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::assign;
        unsigned long (::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::*method_pointer_666a071997205752925f6f9c79933ecf)() const = &::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::size;
        unsigned long (::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::*method_pointer_1c5e0068a9505137a12724ebb63b376c)() const = &::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::max_size;
        void (::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::*method_pointer_fc6d4371e2d0536f805770b3972a26bc)(unsigned long) = &::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::resize;
        void (::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::*method_pointer_abe1cec4a19f54239bbcc86e2ab5f9a7)(unsigned long, enum ::stat_tool::discrete_parametric const &) = &::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::resize;
        void (::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::*method_pointer_186f81303b8f5457bb668fe70be845bf)() = &::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::shrink_to_fit;
        unsigned long (::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::*method_pointer_ff25a8eb9b115e418e7e48e95a5a81eb)() const = &::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::capacity;
        bool (::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::*method_pointer_a39fd1f637215e8ba66e9254ffe86b66)() const = &::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::empty;
        void (::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::*method_pointer_f09c9f0c1bb5595c973c5d486139b498)(unsigned long) = &::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::reserve;
        enum ::stat_tool::discrete_parametric const & (::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::*method_pointer_042123e773175ea9b593e3ad3de9300b)(unsigned long) const = &::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::operator[];
        enum ::stat_tool::discrete_parametric const & (::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::*method_pointer_af6d16cacc705b26a659fb9c730a3479)(unsigned long) const = &::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::at;
        enum ::stat_tool::discrete_parametric const & (::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::*method_pointer_28cd5cbfe85f5fc3b85d987051c17b00)() const = &::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::front;
        enum ::stat_tool::discrete_parametric const & (::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::*method_pointer_a06cca498e995dcf91e0714c9e2859c7)() const = &::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::back;
        void (::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::*method_pointer_fa61ebe74cca5352bcbde273b15590c3)(enum ::stat_tool::discrete_parametric const &) = &::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::push_back;
        void (::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::*method_pointer_779154f311fc5f72ad9b57295b1ae514)() = &::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::pop_back;
        void (::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::*method_pointer_61ea22237ed15ecdad44c1ec127e92b6)(class ::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> > &) = &::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::swap;
        void (::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::*method_pointer_0d494c1958fc5e2cadaeae7759721d20)() = &::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >::clear;
        boost::python::class_< class ::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >, std::shared_ptr< class ::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> > > >("_Vector_b43a1f9fead25c79ab611b9e06ef0f3d", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> > const & >())
            .def("assign", method_pointer_abd86f72c4ac5033b84faac5503beaa8)
            .def("size", method_pointer_666a071997205752925f6f9c79933ecf)
            .def("max_size", method_pointer_1c5e0068a9505137a12724ebb63b376c)
            .def("resize", method_pointer_fc6d4371e2d0536f805770b3972a26bc)
            .def("resize", method_pointer_abe1cec4a19f54239bbcc86e2ab5f9a7)
            .def("shrink_to_fit", method_pointer_186f81303b8f5457bb668fe70be845bf)
            .def("capacity", method_pointer_ff25a8eb9b115e418e7e48e95a5a81eb)
            .def("empty", method_pointer_a39fd1f637215e8ba66e9254ffe86b66)
            .def("reserve", method_pointer_f09c9f0c1bb5595c973c5d486139b498)
            .def("__getitem__", method_pointer_042123e773175ea9b593e3ad3de9300b, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("at", method_pointer_af6d16cacc705b26a659fb9c730a3479, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("front", method_pointer_28cd5cbfe85f5fc3b85d987051c17b00, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("back", method_pointer_a06cca498e995dcf91e0714c9e2859c7, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("push_back", method_pointer_fa61ebe74cca5352bcbde273b15590c3)
            .def("pop_back", method_pointer_779154f311fc5f72ad9b57295b1ae514)
            .def("swap", method_pointer_61ea22237ed15ecdad44c1ec127e92b6)
            .def("clear", method_pointer_0d494c1958fc5e2cadaeae7759721d20);
    struct vector_b43a1f9fead25c79ab611b9e06ef0f3d_from_python
    {
        vector_b43a1f9fead25c79ab611b9e06ef0f3d_from_python()
        {
            boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id< class ::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> > >());
        }

        static void* convertible(PyObject* obj_ptr)
        { return obj_ptr; }

        static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
            void* storage = ((boost::python::converter::rvalue_from_python_storage< class ::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> > >*)data)->storage.bytes;
            new (storage) class ::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >();
            data->convertible = storage;
            class ::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >& result = *((class ::std::vector<stat_tool::discrete_parametric, std::allocator<stat_tool::discrete_parametric> >*)storage);
            unsigned int i = 0;
            for(;; i++)
            {
                boost::python::handle<> py_elem_hdl(boost::python::allow_null(PyIter_Next(obj_iter.get())));
                if(PyErr_Occurred())
                { boost::python::throw_error_already_set(); }
                if(!py_elem_hdl.get())
                { break; }
                boost::python::object py_elem_obj(py_elem_hdl);
                result.push_back((enum ::stat_tool::discrete_parametric)(boost::python::extract< enum ::stat_tool::discrete_parametric >(py_elem_obj)));
            }
        }
    };

    vector_b43a1f9fead25c79ab611b9e06ef0f3d_from_python();
}