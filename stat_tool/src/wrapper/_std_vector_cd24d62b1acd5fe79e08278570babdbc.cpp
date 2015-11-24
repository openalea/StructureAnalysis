#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_vector_cd24d62b1acd5fe79e08278570babdbc()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        void (::std::vector<bool, std::allocator<bool> >::*method_pointer_418f2db7f7e25ae38255d7947e05b835)(unsigned long, bool const &) = &::std::vector<bool, std::allocator<bool> >::assign;
        struct ::std::_Bit_iterator (::std::vector<bool, std::allocator<bool> >::*method_pointer_6d1dbb115fbf531188f7f299955bac1b)() = &::std::vector<bool, std::allocator<bool> >::begin;
        struct ::std::_Bit_const_iterator (::std::vector<bool, std::allocator<bool> >::*method_pointer_56dd11d5947f5454b00f314f5bf377ef)() const = &::std::vector<bool, std::allocator<bool> >::begin;
        struct ::std::_Bit_iterator (::std::vector<bool, std::allocator<bool> >::*method_pointer_6f375ced59ce52e0a221fa59b494a829)() = &::std::vector<bool, std::allocator<bool> >::end;
        struct ::std::_Bit_const_iterator (::std::vector<bool, std::allocator<bool> >::*method_pointer_052b9d8d84605330bd485fcbd5e6d68b)() const = &::std::vector<bool, std::allocator<bool> >::end;
        struct ::std::_Bit_const_iterator (::std::vector<bool, std::allocator<bool> >::*method_pointer_07cd1a353b1d53fd9efaf2699b0688e2)() const = &::std::vector<bool, std::allocator<bool> >::cbegin;
        struct ::std::_Bit_const_iterator (::std::vector<bool, std::allocator<bool> >::*method_pointer_0ec19f671b53592b86a221d861e82eee)() const = &::std::vector<bool, std::allocator<bool> >::cend;
        unsigned long (::std::vector<bool, std::allocator<bool> >::*method_pointer_da55176bbbf253f096d0db9fd2816e43)() const = &::std::vector<bool, std::allocator<bool> >::size;
        unsigned long (::std::vector<bool, std::allocator<bool> >::*method_pointer_25e9e3037ae85438aed888ee77d21021)() const = &::std::vector<bool, std::allocator<bool> >::max_size;
        unsigned long (::std::vector<bool, std::allocator<bool> >::*method_pointer_0e2cc425e52e59519a372d0d4a7462e2)() const = &::std::vector<bool, std::allocator<bool> >::capacity;
        bool (::std::vector<bool, std::allocator<bool> >::*method_pointer_25d2bc0a1f5e5b75a3df011b3b1d6698)() const = &::std::vector<bool, std::allocator<bool> >::empty;
        struct ::std::_Bit_reference (::std::vector<bool, std::allocator<bool> >::*method_pointer_701946b6ae585a25b89ef10c402f1d7d)(unsigned long) = &::std::vector<bool, std::allocator<bool> >::operator[];
        bool (::std::vector<bool, std::allocator<bool> >::*method_pointer_644ab4bbdc6e5ace957a14bcf82bfb95)(unsigned long) const = &::std::vector<bool, std::allocator<bool> >::operator[];
        struct ::std::_Bit_reference (::std::vector<bool, std::allocator<bool> >::*method_pointer_0269d8fcea55583eabdaf082d8fbb19d)(unsigned long) = &::std::vector<bool, std::allocator<bool> >::at;
        bool (::std::vector<bool, std::allocator<bool> >::*method_pointer_b4bc34a43bea5ce49e516807dbb814b2)(unsigned long) const = &::std::vector<bool, std::allocator<bool> >::at;
        void (::std::vector<bool, std::allocator<bool> >::*method_pointer_211a8c447c5b58599202a844e7373586)(unsigned long) = &::std::vector<bool, std::allocator<bool> >::reserve;
        struct ::std::_Bit_reference (::std::vector<bool, std::allocator<bool> >::*method_pointer_433b25775dac5d37b0ab79bf65744c42)() = &::std::vector<bool, std::allocator<bool> >::front;
        bool (::std::vector<bool, std::allocator<bool> >::*method_pointer_d1a95eb25aea533a9a43dc21c4f77261)() const = &::std::vector<bool, std::allocator<bool> >::front;
        struct ::std::_Bit_reference (::std::vector<bool, std::allocator<bool> >::*method_pointer_1997fed230575279aa0ccbca96ca3b10)() = &::std::vector<bool, std::allocator<bool> >::back;
        bool (::std::vector<bool, std::allocator<bool> >::*method_pointer_d99102e2601c53899429eb2d469db19b)() const = &::std::vector<bool, std::allocator<bool> >::back;
        void (::std::vector<bool, std::allocator<bool> >::*method_pointer_9d98ebb88127524dabd40f98e8b50e35)() = &::std::vector<bool, std::allocator<bool> >::data;
        void (::std::vector<bool, std::allocator<bool> >::*method_pointer_76b209658ea3598a8b3d0cbe047cb4fe)(bool) = &::std::vector<bool, std::allocator<bool> >::push_back;
        void (::std::vector<bool, std::allocator<bool> >::*method_pointer_a2b7bf46ecc05a4face42af845e3eeb0)(class ::std::vector<bool, std::allocator<bool> > &) = &::std::vector<bool, std::allocator<bool> >::swap;
        void (*method_pointer_838a405cc4b65da2aa8cec3ebac3a0f8)(struct ::std::_Bit_reference, struct ::std::_Bit_reference) = ::std::vector<bool, std::allocator<bool> >::swap;
        struct ::std::_Bit_iterator (::std::vector<bool, std::allocator<bool> >::*method_pointer_3a80f2beb62d522ba7b0c161eb36ca7a)(struct ::std::_Bit_const_iterator, bool const &) = &::std::vector<bool, std::allocator<bool> >::insert;
        struct ::std::_Bit_iterator (::std::vector<bool, std::allocator<bool> >::*method_pointer_a6f5e26a82e9579a80da1b0d3e425cad)(struct ::std::_Bit_const_iterator, unsigned long, bool const &) = &::std::vector<bool, std::allocator<bool> >::insert;
        void (::std::vector<bool, std::allocator<bool> >::*method_pointer_ca55e32cb6965320b7747d8885bcca90)() = &::std::vector<bool, std::allocator<bool> >::pop_back;
        struct ::std::_Bit_iterator (::std::vector<bool, std::allocator<bool> >::*method_pointer_a62c9adb539c5cbfb55ce8eafa9f6061)(struct ::std::_Bit_const_iterator) = &::std::vector<bool, std::allocator<bool> >::erase;
        struct ::std::_Bit_iterator (::std::vector<bool, std::allocator<bool> >::*method_pointer_4fced5b0c5065a43b6b43503dc8bce86)(struct ::std::_Bit_const_iterator, struct ::std::_Bit_const_iterator) = &::std::vector<bool, std::allocator<bool> >::erase;
        void (::std::vector<bool, std::allocator<bool> >::*method_pointer_9f55e2d993ad540a894ad9119c988ffa)(unsigned long, bool) = &::std::vector<bool, std::allocator<bool> >::resize;
        void (::std::vector<bool, std::allocator<bool> >::*method_pointer_0983d2884e59504e84be145d2a2a97c7)() = &::std::vector<bool, std::allocator<bool> >::shrink_to_fit;
        void (::std::vector<bool, std::allocator<bool> >::*method_pointer_a59861ce1c355ec387d453da3a3f3f3b)() = &::std::vector<bool, std::allocator<bool> >::flip;
        void (::std::vector<bool, std::allocator<bool> >::*method_pointer_d356cc281b5856f4a404f23732065fe4)() = &::std::vector<bool, std::allocator<bool> >::clear;
        boost::python::class_< class ::std::vector<bool, std::allocator<bool> >, std::shared_ptr< class ::std::vector<bool, std::allocator<bool> > > >("_Vector_cd24d62b1acd5fe79e08278570babdbc", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::std::allocator<bool> const & >())
            .def(boost::python::init< unsigned long, class ::std::allocator<bool> const & >())
            .def(boost::python::init< class ::std::vector<bool, std::allocator<bool> > const & >())
            .def("assign", method_pointer_418f2db7f7e25ae38255d7947e05b835)
            .def("begin", method_pointer_6d1dbb115fbf531188f7f299955bac1b)
            .def("begin", method_pointer_56dd11d5947f5454b00f314f5bf377ef)
            .def("end", method_pointer_6f375ced59ce52e0a221fa59b494a829)
            .def("end", method_pointer_052b9d8d84605330bd485fcbd5e6d68b)
            .def("cbegin", method_pointer_07cd1a353b1d53fd9efaf2699b0688e2)
            .def("cend", method_pointer_0ec19f671b53592b86a221d861e82eee)
            .def("size", method_pointer_da55176bbbf253f096d0db9fd2816e43)
            .def("max_size", method_pointer_25e9e3037ae85438aed888ee77d21021)
            .def("capacity", method_pointer_0e2cc425e52e59519a372d0d4a7462e2)
            .def("empty", method_pointer_25d2bc0a1f5e5b75a3df011b3b1d6698)
            .def("__setitem__", method_pointer_701946b6ae585a25b89ef10c402f1d7d)
            .def("__getitem__", method_pointer_644ab4bbdc6e5ace957a14bcf82bfb95)
            .def("at", method_pointer_0269d8fcea55583eabdaf082d8fbb19d)
            .def("at", method_pointer_b4bc34a43bea5ce49e516807dbb814b2)
            .def("reserve", method_pointer_211a8c447c5b58599202a844e7373586)
            .def("front", method_pointer_433b25775dac5d37b0ab79bf65744c42)
            .def("front", method_pointer_d1a95eb25aea533a9a43dc21c4f77261)
            .def("back", method_pointer_1997fed230575279aa0ccbca96ca3b10)
            .def("back", method_pointer_d99102e2601c53899429eb2d469db19b)
            .def("data", method_pointer_9d98ebb88127524dabd40f98e8b50e35)
            .def("push_back", method_pointer_76b209658ea3598a8b3d0cbe047cb4fe)
            .def("swap", method_pointer_a2b7bf46ecc05a4face42af845e3eeb0)
            .def("swap", method_pointer_838a405cc4b65da2aa8cec3ebac3a0f8)
            .def("insert", method_pointer_3a80f2beb62d522ba7b0c161eb36ca7a)
            .def("insert", method_pointer_a6f5e26a82e9579a80da1b0d3e425cad)
            .def("pop_back", method_pointer_ca55e32cb6965320b7747d8885bcca90)
            .def("erase", method_pointer_a62c9adb539c5cbfb55ce8eafa9f6061)
            .def("erase", method_pointer_4fced5b0c5065a43b6b43503dc8bce86)
            .def("resize", method_pointer_9f55e2d993ad540a894ad9119c988ffa)
            .def("shrink_to_fit", method_pointer_0983d2884e59504e84be145d2a2a97c7)
            .def("flip", method_pointer_a59861ce1c355ec387d453da3a3f3f3b)
            .def("clear", method_pointer_d356cc281b5856f4a404f23732065fe4)
            .staticmethod("swap");
    struct vector_cd24d62b1acd5fe79e08278570babdbc_from_python
    {
        vector_cd24d62b1acd5fe79e08278570babdbc_from_python()
        {
            boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id< class ::std::vector<bool, std::allocator<bool> > >());
        }

        static void* convertible(PyObject* obj_ptr)
        { return obj_ptr; }

        static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
            void* storage = ((boost::python::converter::rvalue_from_python_storage< class ::std::vector<bool, std::allocator<bool> > >*)data)->storage.bytes;
            new (storage) class ::std::vector<bool, std::allocator<bool> >();
            data->convertible = storage;
            class ::std::vector<bool, std::allocator<bool> >& result = *((class ::std::vector<bool, std::allocator<bool> >*)storage);
            unsigned int i = 0;
            for(;; i++)
            {
                boost::python::handle<> py_elem_hdl(boost::python::allow_null(PyIter_Next(obj_iter.get())));
                if(PyErr_Occurred())
                { boost::python::throw_error_already_set(); }
                if(!py_elem_hdl.get())
                { break; }
                boost::python::object py_elem_obj(py_elem_hdl);
                result.push_back((bool)(boost::python::extract< bool >(py_elem_obj)));
            }
        }
    };

    vector_cd24d62b1acd5fe79e08278570babdbc_from_python();
}