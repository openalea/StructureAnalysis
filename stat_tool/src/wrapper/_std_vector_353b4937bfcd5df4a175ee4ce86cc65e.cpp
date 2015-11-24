#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_vector_353b4937bfcd5df4a175ee4ce86cc65e()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        void (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_1e960f79d6585bc38acb49374308d0c7)(unsigned long, class ::stat_tool::MultiPlot const &) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::assign;
        unsigned long (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_b27a21b3a1b359ad9e1525ba95889172)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::size;
        unsigned long (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_bfbaf1f6ef5f5de0894d99e9cb08d1e8)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::max_size;
        void (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_d3d1d6bcd36552b3a9d00a12c7fb5914)(unsigned long) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::resize;
        void (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_2716e006579e5128b03fb3359850bb6f)(unsigned long, class ::stat_tool::MultiPlot const &) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::resize;
        void (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_787a19fa0c825310b116e3af10c0d1e1)() = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::shrink_to_fit;
        unsigned long (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_578e2742f84d5f72ab2df2174f633fca)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::capacity;
        bool (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_a608e3fca0bb505fb88bb9ecf9f98e83)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::empty;
        void (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_9cb2aa9a2b0a551a81a050b90ba90aba)(unsigned long) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::reserve;
        class ::stat_tool::MultiPlot & (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_b33e92976e59532b99563f931c194509)(unsigned long) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::operator[];
        class ::stat_tool::MultiPlot const & (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_fe27993b65a854c987b983ceff076a8a)(unsigned long) const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::operator[];
        class ::stat_tool::MultiPlot & (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_948799718e2b5295a98fef05cc84c6b9)(unsigned long) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::at;
        class ::stat_tool::MultiPlot const & (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_a2a8a118e903577ba230286a13e4a516)(unsigned long) const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::at;
        class ::stat_tool::MultiPlot & (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_f59ba0fe4c39568c8e96bc9199f317ca)() = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::front;
        class ::stat_tool::MultiPlot const & (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_312d887aa45e5a9186391111122428ca)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::front;
        class ::stat_tool::MultiPlot & (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_be9d3b9c32c85cf8b39f812f835ad2b4)() = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::back;
        class ::stat_tool::MultiPlot const & (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_42114a1d3c9d5702857cb80e08e89ade)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::back;
        class ::stat_tool::MultiPlot * (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_2c3088b9e89d521da686542b25d2069d)() = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::data;
        class ::stat_tool::MultiPlot const * (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_8f520dbcc16658198b3da9bbae7c71ad)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::data;
        void (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_42da52736f6b52c3a8874835f3062487)(class ::stat_tool::MultiPlot const &) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::push_back;
        void (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_7022ca3e7c5557f089da3a139313c76a)() = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::pop_back;
        void (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_3cd5f3b75e0f57bebbecbe69f20456a4)(class ::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > &) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::swap;
        void (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_3468843386c55add8cf2668640bb3df0)() = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::clear;
        boost::python::class_< class ::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >, std::shared_ptr< class ::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > > >("_Vector_353b4937bfcd5df4a175ee4ce86cc65e", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > const & >())
            .def("assign", method_pointer_1e960f79d6585bc38acb49374308d0c7)
            .def("size", method_pointer_b27a21b3a1b359ad9e1525ba95889172)
            .def("max_size", method_pointer_bfbaf1f6ef5f5de0894d99e9cb08d1e8)
            .def("resize", method_pointer_d3d1d6bcd36552b3a9d00a12c7fb5914)
            .def("resize", method_pointer_2716e006579e5128b03fb3359850bb6f)
            .def("shrink_to_fit", method_pointer_787a19fa0c825310b116e3af10c0d1e1)
            .def("capacity", method_pointer_578e2742f84d5f72ab2df2174f633fca)
            .def("empty", method_pointer_a608e3fca0bb505fb88bb9ecf9f98e83)
            .def("reserve", method_pointer_9cb2aa9a2b0a551a81a050b90ba90aba)
            .def("__setitem__", method_pointer_b33e92976e59532b99563f931c194509, boost::python::return_internal_reference<>())
            .def("__getitem__", method_pointer_fe27993b65a854c987b983ceff076a8a, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("at", method_pointer_948799718e2b5295a98fef05cc84c6b9, boost::python::return_internal_reference<>())
            .def("at", method_pointer_a2a8a118e903577ba230286a13e4a516, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("front", method_pointer_f59ba0fe4c39568c8e96bc9199f317ca, boost::python::return_internal_reference<>())
            .def("front", method_pointer_312d887aa45e5a9186391111122428ca, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("back", method_pointer_be9d3b9c32c85cf8b39f812f835ad2b4, boost::python::return_internal_reference<>())
            .def("back", method_pointer_42114a1d3c9d5702857cb80e08e89ade, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("data", method_pointer_2c3088b9e89d521da686542b25d2069d, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("data", method_pointer_8f520dbcc16658198b3da9bbae7c71ad, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("push_back", method_pointer_42da52736f6b52c3a8874835f3062487)
            .def("pop_back", method_pointer_7022ca3e7c5557f089da3a139313c76a)
            .def("swap", method_pointer_3cd5f3b75e0f57bebbecbe69f20456a4)
            .def("clear", method_pointer_3468843386c55add8cf2668640bb3df0);
    struct vector_353b4937bfcd5df4a175ee4ce86cc65e_from_python
    {
        vector_353b4937bfcd5df4a175ee4ce86cc65e_from_python()
        {
            boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id< class ::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > >());
        }

        static void* convertible(PyObject* obj_ptr)
        { return obj_ptr; }

        static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
            void* storage = ((boost::python::converter::rvalue_from_python_storage< class ::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > >*)data)->storage.bytes;
            new (storage) class ::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >();
            data->convertible = storage;
            class ::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >& result = *((class ::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >*)storage);
            unsigned int i = 0;
            for(;; i++)
            {
                boost::python::handle<> py_elem_hdl(boost::python::allow_null(PyIter_Next(obj_iter.get())));
                if(PyErr_Occurred())
                { boost::python::throw_error_already_set(); }
                if(!py_elem_hdl.get())
                { break; }
                boost::python::object py_elem_obj(py_elem_hdl);
                result.push_back((class ::stat_tool::MultiPlot)(boost::python::extract< class ::stat_tool::MultiPlot >(py_elem_obj)));
            }
        }
    };

    vector_353b4937bfcd5df4a175ee4ce86cc65e_from_python();
}