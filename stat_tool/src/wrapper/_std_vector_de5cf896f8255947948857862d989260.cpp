#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_vector_de5cf896f8255947948857862d989260()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        void (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_e5d231f897de572e95406ba1a8108d23)(unsigned long, class ::stat_tool::SinglePlot const &) = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::assign;
        unsigned long (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_46f905bd25fc504d9118b3382811bc21)() const = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::size;
        unsigned long (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_2362da34d1295f96b81f47de686d576e)() const = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::max_size;
        void (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_5505c5cdbb0b59b089b98148cad8dcdd)(unsigned long) = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::resize;
        void (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_a816a3a5bf2c5e30817920db01fcb4a0)(unsigned long, class ::stat_tool::SinglePlot const &) = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::resize;
        void (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_10f3d326100653ecbab3ed3aa7c8f1e5)() = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::shrink_to_fit;
        unsigned long (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_f32749e8968e59bf88fb7712cbe95dc6)() const = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::capacity;
        bool (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_f910c980fcb555b29a5f344821c26e5b)() const = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::empty;
        void (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_156545923095563894f270b3746a29a2)(unsigned long) = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::reserve;
        class ::stat_tool::SinglePlot & (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_719718c81f6a51b0a2b630aa51ae9a5c)(unsigned long) = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::operator[];
        class ::stat_tool::SinglePlot const & (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_3d2434fc12085bcc993c928ece231f73)(unsigned long) const = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::operator[];
        class ::stat_tool::SinglePlot & (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_dd459506874d5994bccef1c5e3cc84fe)(unsigned long) = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::at;
        class ::stat_tool::SinglePlot const & (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_1356eb1458e75eb4996fb44014aefdbf)(unsigned long) const = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::at;
        class ::stat_tool::SinglePlot & (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_cde854cf7ad65ba68fb842a2dea67a28)() = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::front;
        class ::stat_tool::SinglePlot const & (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_6aa2e3bac6c55e44b18ffcbac0513082)() const = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::front;
        class ::stat_tool::SinglePlot & (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_ef97643dd6ce563498a314e8f524b499)() = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::back;
        class ::stat_tool::SinglePlot const & (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_331bdcc41f045f2b993eacd03767e5ea)() const = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::back;
        class ::stat_tool::SinglePlot * (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_563562405af6516ead9daa982e4f96a3)() = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::data;
        class ::stat_tool::SinglePlot const * (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_227ea08fbaf75560b11ab36b07c31541)() const = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::data;
        void (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_8f85a7030a125725b8ee45d862b748c7)(class ::stat_tool::SinglePlot const &) = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::push_back;
        void (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_1c02d61525af53cf8fa914fd5a8ae7f4)() = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::pop_back;
        void (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_0058f91a44af5e22910f8f742735503f)(class ::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> > &) = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::swap;
        void (::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::*method_pointer_546b94b9a8f15aca9a6e00823a3bb3ac)() = &::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >::clear;
        boost::python::class_< class ::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >, std::shared_ptr< class ::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> > > >("_Vector_de5cf896f8255947948857862d989260", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> > const & >())
            .def("assign", method_pointer_e5d231f897de572e95406ba1a8108d23)
            .def("size", method_pointer_46f905bd25fc504d9118b3382811bc21)
            .def("max_size", method_pointer_2362da34d1295f96b81f47de686d576e)
            .def("resize", method_pointer_5505c5cdbb0b59b089b98148cad8dcdd)
            .def("resize", method_pointer_a816a3a5bf2c5e30817920db01fcb4a0)
            .def("shrink_to_fit", method_pointer_10f3d326100653ecbab3ed3aa7c8f1e5)
            .def("capacity", method_pointer_f32749e8968e59bf88fb7712cbe95dc6)
            .def("empty", method_pointer_f910c980fcb555b29a5f344821c26e5b)
            .def("reserve", method_pointer_156545923095563894f270b3746a29a2)
            .def("__setitem__", method_pointer_719718c81f6a51b0a2b630aa51ae9a5c, boost::python::return_internal_reference<>())
            .def("__getitem__", method_pointer_3d2434fc12085bcc993c928ece231f73, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("at", method_pointer_dd459506874d5994bccef1c5e3cc84fe, boost::python::return_internal_reference<>())
            .def("at", method_pointer_1356eb1458e75eb4996fb44014aefdbf, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("front", method_pointer_cde854cf7ad65ba68fb842a2dea67a28, boost::python::return_internal_reference<>())
            .def("front", method_pointer_6aa2e3bac6c55e44b18ffcbac0513082, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("back", method_pointer_ef97643dd6ce563498a314e8f524b499, boost::python::return_internal_reference<>())
            .def("back", method_pointer_331bdcc41f045f2b993eacd03767e5ea, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("data", method_pointer_563562405af6516ead9daa982e4f96a3, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("data", method_pointer_227ea08fbaf75560b11ab36b07c31541, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("push_back", method_pointer_8f85a7030a125725b8ee45d862b748c7)
            .def("pop_back", method_pointer_1c02d61525af53cf8fa914fd5a8ae7f4)
            .def("swap", method_pointer_0058f91a44af5e22910f8f742735503f)
            .def("clear", method_pointer_546b94b9a8f15aca9a6e00823a3bb3ac);
    struct vector_de5cf896f8255947948857862d989260_from_python
    {
        vector_de5cf896f8255947948857862d989260_from_python()
        {
            boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id< class ::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> > >());
        }

        static void* convertible(PyObject* obj_ptr)
        { return obj_ptr; }

        static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
            void* storage = ((boost::python::converter::rvalue_from_python_storage< class ::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> > >*)data)->storage.bytes;
            new (storage) class ::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >();
            data->convertible = storage;
            class ::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >& result = *((class ::std::vector<stat_tool::SinglePlot, std::allocator<stat_tool::SinglePlot> >*)storage);
            unsigned int i = 0;
            for(;; i++)
            {
                boost::python::handle<> py_elem_hdl(boost::python::allow_null(PyIter_Next(obj_iter.get())));
                if(PyErr_Occurred())
                { boost::python::throw_error_already_set(); }
                if(!py_elem_hdl.get())
                { break; }
                boost::python::object py_elem_obj(py_elem_hdl);
                result.push_back((class ::stat_tool::SinglePlot)(boost::python::extract< class ::stat_tool::SinglePlot >(py_elem_obj)));
            }
        }
    };

    vector_de5cf896f8255947948857862d989260_from_python();
}