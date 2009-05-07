


// convert a boost python list into a CPP array
// INPUT is the input variable name
// TYPE is the type of the CPP array to be returned into the variable data
// returns data, the list of object and size its length
#define CREATE_ARRAY(INPUT, TYPE)\
  int size = len(INPUT);\
  sequence_analysis::wrap_util::auto_ptr_array<TYPE> data(new TYPE[size]);\
  for (int i = 0; i < size; i++)\
      data[i] = boost::python::extract<TYPE>(INPUT[i]);


// prototype not included since the name cannot be overloaded
// INPUT the input argument
// METHOD_NAME: name of the method that will be called on INPUT
// OUTPUT_TYPE: type returned by METHOD_NAME
// optional list of arguments to be used by METHOD_NAME
#define SIMPLE_METHOD_TEMPLATE_1(INPUT, METHOD_NAME, OUTPUT_TYPE, ...)\
  Format_error error; \
    OUTPUT_TYPE* ret;\
    ret = INPUT.METHOD_NAME(error, __VA_ARGS__ );\
    if (!ret)\
      sequence_analysis::wrap_util::throw_error(error);\
    return ret;\

#define SIMPLE_METHOD_TEMPLATE_0(INPUT, METHOD_NAME, OUTPUT_TYPE)\
  Format_error error; \
    OUTPUT_TYPE* ret;\
    ret = INPUT.METHOD_NAME(error);\
    if (!ret)\
      sequence_analysis::wrap_util::throw_error(error);\
    return ret;\




// Specific methods
// ---------------------- file_ascii_write (ascii_write) method ----------------
#define WRAP_METHOD_FILE_ASCII_WRITE(INPUT_CLASS) \
  static void file_ascii_write(const INPUT_CLASS& m, const char* path, bool exhaustive)\
  {\
     bool result = true;\
     Format_error error;\
     result = m.ascii_write(error, path, exhaustive);\
     if (!result)\
        sequence_analysis::wrap_util::throw_error(error);\
   }

// ---------------------- ascii_write method -----------------------------------
#define WRAP_METHOD_ASCII_WRITE(INPUT_CLASS) \
static std::string ascii_write(const INPUT_CLASS& input, bool exhaustive)\
  {\
    std::stringstream s;\
    std::string res;\
    input.ascii_write(s, exhaustive);\
    res = s.str();\
    return res;\
  }

// ---------------------- plot_write method ------------------------------------
#define WRAP_METHOD_PLOT_WRITE(INPUT_CLASS) \
  static void plot_write(const INPUT_CLASS& input, const std::string& prefix, const std::string& title)\
  {\
	  Format_error error;\
      if(! input.plot_write(error, prefix.c_str(), title.c_str()))\
          sequence_analysis::wrap_util::throw_error(error);\
  }

// ---------------------- spreadshhet_write method -----------------------------
#define WRAP_METHOD_SPREADSHEET_WRITE(INPUT_CLASS) \
  static void spreadsheet_write(const INPUT_CLASS& input, const std::string& filename)\
  {\
    Format_error error;\
    if(! input.spreadsheet_write(error, filename.c_str()))\
       sequence_analysis::wrap_util::throw_error(error);\
  }


// ---------------------- survival_get_plotable --------------------------------
#define WRAP_METHOD_SURVIVAL_GET_PLOTABLE(INPUT_CLASS) \
static MultiPlotSet* survival_get_plotable(const INPUT_CLASS& p) \
  { \
    Format_error error; \
    MultiPlotSet* ret = p.survival_get_plotable(error); \
    if(!ret)\
      sequence_analysis::wrap_util::throw_error(error);\
    return ret;\
  }

// ---------------------- survival_plot_write ---------------------------------
#define WRAP_METHOD_SURVIVAL_PLOT_WRITE(INPUT_CLASS) \
  static void survival_plot_write(const INPUT_CLASS& p, const std::string& prefix, const std::string& title) \
  {\
    Format_error error;\
    if(!p.survival_plot_write(error, prefix.c_str(), title.c_str()))\
      sequence_analysis::wrap_util::throw_error(error);\
  }


// ---------------------- survival_spreadsheet_write ---------------------------------------
#define WRAP_METHOD_SURVIVAL_SPREADSHEET_WRITE(INPUT_CLASS) \
 static void survival_spreadsheet_write(const INPUT_CLASS& p, const std::string& filename) \
  {\
    Format_error error;\
    if(!p.survival_spreadsheet_write(error, filename.c_str()))\
      sequence_analysis::wrap_util::throw_error(error);\
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
	if (VALUE< MIN || VALUE>=MAX)\
     {\
    	 PyErr_SetString(PyExc_TypeError,(error_message.str()).c_str());\
    	 throw_error_already_set();\
     }\


// boost python declarations -------------------------------------------
//
// quick alias to avoid writting the return_value_policy...


#define DEF_RETURN_VALUE(NAME, REFERENCE, ARGS, DOCSTRING) \
    .def(NAME, REFERENCE, return_value_policy< manage_new_object >(), ARGS, DOCSTRING)

#define DEF_RETURN_VALUE_NO_ARGS(NAME, REFERENCE, DOCSTRING) \
    .def(NAME, REFERENCE, return_value_policy< manage_new_object >(), DOCSTRING)

#define DEF_INIT_MAKE_CONSTRUCTOR(WRAPPED_FUNCTION, DOCSTRING) \
	.def("__init__", make_constructor(WRAPPED_FUNCTION), DOCSTRING)

