



#define ARGS boost::python::args




// calls METHOD_NAME of the INPUT_CLASS
// :param input_class: a class instance
// :param var1: a variable
// :returns: OUTPUT_TYPE
#define WRAP_METHOD1(INPUT_TYPE, METHOD_NAME, OUTPUT_TYPE, VARTYPE1) \
static OUTPUT_TYPE* METHOD_NAME(const INPUT_TYPE& input_class, VARTYPE1 var1) \
{\
    Format_error error; \
    OUTPUT_TYPE* ret = NULL; \
    ret = input_class.METHOD_NAME(error, var1); \
    if(!ret) stat_tool::wrap_util::throw_error(error); \
    return ret; \
}
// same but no arguments
#define WRAP_METHOD0(INPUT_TYPE, METHOD_NAME, OUTPUT_TYPE) \
static OUTPUT_TYPE* METHOD_NAME(const INPUT_TYPE& input_class) \
{\
    Format_error error; \
    OUTPUT_TYPE* ret = NULL; \
    ret = input_class.METHOD_NAME(error); \
    if(!ret) stat_tool::wrap_util::throw_error(error); \
    return ret; \
}



// Specific methods
// ---------------------- file_ascii_write method ---------------------------------------
#define WRAP_METHOD_FILE_ASCII_WRITE(INPUT_CLASS) \
  static void file_ascii_write(const INPUT_CLASS& m, const char* path, bool exhaustive)\
  {\
     bool result = true;\
     Format_error error;\
     result = m.ascii_write(error, path, exhaustive);\
     if (!result)\
        stat_tool::wrap_util::throw_error(error);\
   }\

// ---------------------- survival_get_plotable ---------------------------------------
#define WRAP_METHOD_SURVIVAL_GET_PLOTABLE(INPUT_CLASS) \
static MultiPlotSet* survival_get_plotable(const INPUT_CLASS& p) \
  { \
    Format_error error; \
    MultiPlotSet* ret = p.survival_get_plotable(error); \
    if(!ret)\
      stat_tool::wrap_util::throw_error(error);\
    return ret;\
  }

// ---------------------- survival_plot_write ---------------------------------------
#define WRAP_METHOD_SURVIVAL_PLOT_WRITE(INPUT_CLASS) \
  static void survival_plot_write(const INPUT_CLASS& p, const std::string& prefix, const std::string& title) \
  {\
    Format_error error;\
    if(!p.survival_plot_write(error, prefix.c_str(), title.c_str()))\
      stat_tool::wrap_util::throw_error(error);\
  }


// ---------------------- survival_spreadsheet_write ---------------------------------------
#define WRAP_METHOD_SURVIVAL_SPREADSHEET_WRITE(INPUT_CLASS) \
 static void survival_spreadsheet_write(const INPUT_CLASS& p, const std::string& filename) \
  {\
    Format_error error;\
    if(!p.survival_spreadsheet_write(error, filename.c_str()))\
      stat_tool::wrap_util::throw_error(error);\
  }

// ---------------------- survival_ascii_write ---------------------------------------
#define WRAP_METHOD_SURVIVAL_ASCII_WRITE(INPUT_CLASS) \
static std::string survival_ascii_write(const INPUT_CLASS& p) \
  {\
    std::stringstream s; \
    std::string res; \
    p.survival_ascii_write(s);\
    res = s.str();\
    return res;\
  }






// INPUT PARAMETER VALIDATION -------------



#define CHECK(VALUE, MIN, MAX) \
	try{\
		if (VALUE< MIN || VALUE>=MAX)\
         { PyErr_SetString(PyExc_TypeError,\
        		(error_message.str()).c_str());\
        	berror = true;\
		  }\
		  }\
		     catch(...)\
		  {\
			  berror = true;\
           }






// boost python declarations -------------------------------------------
//
// quick alias to avoid writting the return_value_policy...
#define DEF_RETURN_VALUE(NAME, REFERENCE, ARGS, DOCSTRING) \
    .def(NAME, REFERENCE, return_value_policy< manage_new_object >(), ARGS, DOCSTRING)

#define DEF_RETURN_VALUE_NO_ARGS(NAME, REFERENCE, DOCSTRING) \
    .def(NAME, REFERENCE, return_value_policy< manage_new_object >(), DOCSTRING)

#define DEF_INIT_MAKE_CONSTRUCTOR(WRAPPED_FUNCTION, DOCSTRING) \
	.def("__init__", make_constructor(WRAPPED_FUNCTION), DOCSTRING)



// def("__len__", &Class::method, "docstring")
#define DEF_LEN(CLASS, FUNCTION_NAME) \
    .def("__len__", &CLASS::FUNCTION_NAME, "Return the size of the Class instance")

#define DEF_STR() .def(self_ns::str(self)) // __str__

