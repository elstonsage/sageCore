#include "func/PythonInterface.h"

jmp_buf context;

namespace SAGE {
namespace FUNC {

//============================================================================
//
//  CONSTRUCTOR
//
//============================================================================
PythonInterface::PythonInterface(cerrorstream errors)
    : my_errors   (errors), 
      my_globals  (NULL), 
      my_obj_code (NULL)
{
  // These flags make Python act as a single unit executable without certain options.
  // Otherwise, python produces spurious warnings.
  Py_NoSiteFlag             = 1;
  Py_FrozenFlag             = 1;
  Py_UseClassExceptionsFlag = 0;

  // Initialize python:
  Py_Initialize();
  
  empty_python_path();

  // Import math, parser:

  PyRun_SimpleString("from math import *");
  PyRun_SimpleString("import parser");

  // Add the 'extract_names' function to the Python environment:

  define_extract_names();
  define_marker_models();
  
  PyObject*   main_module = PyImport_ImportModule("__main__");
  
  // Produces borrowed reference.  
  //
  my_globals = PyModule_GetDict(main_module);        
  Py_DECREF(main_module);
  
  restrict_python();
}

//=====================================================
//
//  DESTRUCTOR
//
//=====================================================
PythonInterface::~PythonInterface()
{
  Py_XDECREF(my_obj_code);
  
  unrestrict_python();
  Py_Finalize();
}

//====================================================
//
//  runPythonString(...)
//
//====================================================
bool 
PythonInterface::runPythonString(const std::string & source) const
{
  return PyRun_SimpleString(const_cast<char*>(source.c_str())) == 0;
}

//====================================================
//
//  isNameInEnvironment(...)
//
//====================================================
bool 
PythonInterface::isNameInEnvironment(const std::string & name) const
{
  return PyDict_GetItemString(my_globals, const_cast<char*>(name.c_str())) != NULL;
}

//=====================================================
//
//  names(...)
//
//=====================================================
std::vector<std::string>
PythonInterface::getNameList(const std::string& expr)
{
  if(!addStringToEnvironment("expression", expr))
    throw PythonException();

  if(!runPythonString("names = extract_names(parser.expr(expression).totuple()).keys()"))
  {
    my_errors << py_error_msg() << "\n";

    throw PythonException();
  }

  std::vector<std::string> names;

  PyObject * py_names = PyDict_GetItemString(my_globals, "names");

  if(py_names && PyList_Check(py_names))
  {
    for(int i = 0; i < PyList_Size(py_names); ++i)
    {
      PyObject * name = PyList_GetItem(py_names, i);

      if(!name)
        continue;

      if(!PyString_Check(name))
        continue;

      char * c_name = PyString_AsString(name);
              
      if(c_name)
        names.push_back(std::string(c_name));
    }
  }
  
  return names;
}

// - Create Python byte code from expression.  Return value indicates
//   whether operation succeded.
// 
bool                
PythonInterface::compile(string expr) 
{
  my_obj_code = Py_CompileString(const_cast<char*>(expr.c_str()), "", Py_eval_input);

  if(!my_obj_code)
  {
    my_errors << py_error_msg() << "\n";
    return false;
  }
  else
  {
    return true;
  }
}
    
//==================================================
//
//  addNoneToEnvironment(...)
//
//==================================================
bool
PythonInterface::addNoneToEnvironment(const std::string & name) const
{
  return my_globals ? PyDict_SetItemString(my_globals, const_cast<char*>(name.c_str()), Py_None) == 0 : false;
}

//==================================================
//
//  addDoubleToEnvironment(...)
//
//==================================================
bool
PythonInterface::addDoubleToEnvironment(const std::string & name, double value) const
{
  bool       succeeded = false;
  PyObject * p_name    = PyString_FromString(const_cast<char*>(name.c_str()));
  PyObject * p_value   = PyFloat_FromDouble(value);
  
  if(my_globals && p_name && p_value)
    succeeded = PyDict_SetItem(my_globals, p_name, p_value) == 0;
  
  Py_XDECREF(p_name);
  Py_XDECREF(p_value);

  return succeeded;
}

//===================================================
//
//  addStringToEnvironment(...)
//
//===================================================
bool
PythonInterface::addStringToEnvironment(const std::string & name, const std::string & value) const
{
  bool       succeeded = false;
  PyObject * p_name    = PyString_FromString(const_cast<char*>(name.c_str()));
  PyObject * p_value   = PyString_FromString(value.c_str());
  
  if(my_globals && p_name && p_value)
    succeeded = PyDict_SetItem(my_globals, p_name, p_value) == 0;
  
  Py_XDECREF(p_name);
  Py_XDECREF(p_value);

  return succeeded;
}

//======================================
//
//  addAlleleListToEnvironment(...)
//
//======================================
bool
PythonInterface::addAlleleListToEnvironment(const std::string & name, const std::string & allele1, const std::string & allele2) const
{
  if(!my_globals || (name == "") || (allele1 == "") || (allele2 == ""))
  {
    return false;
  }
  else
  {
    return runPythonString(
      name + " = []" + "\n" +
      name + ".append(\"" + allele1 + "\")" + "\n" +
      name + ".append(\"" + allele2 + "\")" + "\n");
  }
}

// - Handler for alarm signal.
//
void
PythonInterface::timeout(int return_value)
{
#ifndef NO_TIMER
#ifndef SUN
#ifndef NO_POSIX_SIGNALS
  longjmp(context, return_value);
#endif
#endif
#endif
}

//==============================================
//
//  calculateValueInEnvironment
//
//==============================================
bool
PythonInterface::calculateValueInEnvironment(
  const    std::string & name, 
  const    std::string & expr,
  unsigned int           time_limit) 
{
  bool  success = false;

  string  source = "";
  source += name + " = " + expr + "\n";
  char*  c_source = const_cast<char*>(source.c_str());
  
  // - Set a timer so that Python code won't run indefinitely.
  //
  // NOTE: for reasons that aren't clear, can't write to my_errors 
  //       after longjmp() has occured.
  //

#ifndef NO_TIMER
#ifndef SUN
#ifndef NO_POSIX_SIGNALS
  if(setjmp(context))
  { 
    exit_on_timeout(expr); 
  }
  signal(SIGALRM, timeout);
  alarm(time_limit);
#endif
#endif
#endif

  PyObject* result = PyRun_String(c_source, Py_file_input, my_globals, my_globals);

#ifndef NO_TIMER
#ifndef SUN
#ifndef NO_POSIX_SIGNALS
  alarm(DISABLE_TIMER);
#endif  
#endif
#endif

  if(result)
  {
    success = true;
    Py_DECREF(result);
  }
  else 
  {
    my_errors << py_error_msg() << "\n";
  }
  
  return success;
}

// - Empty the path which Python uses to search for modules
//   except for the item for builtin modules.
//
void
PythonInterface::empty_python_path() const
{
  string  source = "";
  
  source += "import sys\n";
  source += "   \n";
  source += "#new_list = []  \n";
  source += "#for item in sys.path:  \n";
  source += "#  if len(item) > 6:  \n";
  source += "#    if item[-7:] == 'dynload':  \n"; 
  source += "#      new_list.append(item)  \n";
  source += "   \n";
  source += "sys.path = []   \n";
  source += "   \n";
  source += "del sys\n";
  
  char*  c_source = const_cast<char*>(source.c_str());
  PyRun_SimpleString(c_source);
}

// - Limit the Python environment wh. is accessible.
//
void
PythonInterface::restrict_python() const
{
  std::string source = "";
  
  // Selected built in names.
  source += "allowed = [";
  source += "'float', 'filter', 'round', 'reduce', 'repr', 'str', 'len', ";
  source += "'cmp', 'chr', 'abs', 'divmod', 'slice', 'tuple', 'ord', 'pow', ";
  source += "'eval', 'min', 'max', 'map', 'hex', 'int', 'None', 'type'";  
  source += "]";
  
  char* c_source = const_cast<char*>(source.c_str());
  PyRun_SimpleString(c_source);
  
  source = "";
  source += "disallowed = {}\n";
  source += "for name in __builtins__.__dict__.keys():\n";
  source += "  if not name in allowed:\n";
  source += "    disallowed[name] = __builtins__.__dict__[name]\n";
  source += "    del __builtins__.__dict__[name]\n";
  
  c_source = const_cast<char*>(source.c_str());
  PyRun_SimpleString(c_source);
}

// - Remove restrictions on Python environment.
//
void
PythonInterface::unrestrict_python() const
{
  std::string  source = "";
  
  source += "for name in disallowed.keys():\n";
  source += "  __builtins__.__dict__[name] = disallowed[name]\n";
  
  char*  c_source = const_cast<char*>(source.c_str());
  PyRun_SimpleString(c_source);
}


// - Run my_obj_code and return the result.  Return quiet::NaN 
//   if operation fails.  
// 
double                
PythonInterface::run(const std::string & expr, unsigned int time_limit)
{
  double return_value = numeric_limits<double>::quiet_NaN();

  if(my_globals && my_obj_code)
  {

    // Set a timer so that Python code won't run indefinitely.
    // NOTE: for reasons that aren't clear, can't write to my_errors after longjmp() has occured.

#ifndef NO_TIMER
#ifndef SUN    
#ifndef NO_POSIX_SIGNALS
    if(setjmp(context))
    { 
      exit_on_timeout(expr);
    }
    signal(SIGALRM, timeout);
    alarm(time_limit);
#endif
#endif
#endif
    
    PyObject* result = PyEval_EvalCode((PyCodeObject*)my_obj_code, my_globals, my_globals);
                                        
#ifndef NO_TIMER
#ifndef SUN                                        
#ifndef NO_POSIX_SIGNALS
    alarm(DISABLE_TIMER);
#endif
#endif
#endif

    if(result)
    {
      if(PyInt_Check(result))
      {
        long long_result = PyInt_AsLong(result);
        return_value = static_cast<double>(long_result);
      }
      else if(PyFloat_Check(result))
      {
        return_value = PyFloat_AsDouble(result);
      }
    }
    else
    {
      my_errors << py_error_msg() << "\n";
    }
    
    Py_XDECREF(result);
  }

  return return_value;
}
      
// - Return Python exception information.
//
string 
PythonInterface::py_error_msg() const
{
  string message = "";

  PyObject* type;
  PyObject* value;
  PyObject* traceback;
  
  PyErr_Fetch(&type, &value, &traceback);
  if(type)
  {
    PyObject* py_str = PyObject_Str(type);
    if(py_str)
    {
      char* type_string = PyString_AsString(py_str);
      if(type_string)
      {
        message += string(type_string) + "  ";
      }
      Py_DECREF(py_str);
    }
    Py_DECREF(type);
  }
  
  if(value)
  {
    PyObject* py_str = PyObject_Str(value);
    if(py_str)
    {
      char* value_string = PyString_AsString(py_str);
      if(value_string)
      {
        message += string(value_string) + "  ";
      }
      Py_DECREF(py_str);
    }
    Py_DECREF(value);
  }
  
  if(traceback)
  {
    PyObject* py_str = PyObject_Str(traceback);
    if(py_str)
    {
      char* traceback_string = PyString_AsString(py_str);
      if(traceback_string)
      {
        message += string(traceback_string) + "  ";
      }
      Py_DECREF(py_str);
    }
    Py_DECREF(traceback);
  }
  
  if(message == "")
  {
    message += "Couldn't intercept python exception.";
  }
  
#ifdef PY_ERROR_MSG
  return message;
#else
  return "";
#endif
}

// - Define function in the Python environment to extract names from
//   a Python parse tree.
//
void
PythonInterface::define_extract_names() const
{
  std::ostringstream source_code;

  source_code << "# - Define NAME locally to avoid importing token." << std::endl
              << "#"                                                 << std::endl
              << "NAME = 1"                                          << std::endl
              << "def extract_names(tree, names = None):"            << std::endl
              << "  if names == None:"                               << std::endl
              << "    names = {}"                                    << std::endl
              << "  #if type(tree) is types.TupleType:"              << std::endl
              << "  if type(tree) is type(()):"                      << std::endl
              << "    if len(tree) == 2:"                            << std::endl
              << "      if tree[0] == NAME:"                         << std::endl
              << "        names[tree[1]] = 0"                        << std::endl
              << "        return names"                              << std::endl
              << "    for node in tree:"                             << std::endl
              << "      extract_names(node, names)"                  << std::endl
              << "  return names"                                    << std::endl;
  
  PyRun_SimpleString(const_cast<char*>(source_code.str().c_str()));
}

// - Define in Python environment functions to calculate a trait value
//   based on marker alleles using dominant, recessive, and additive
//   models.
//
void
PythonInterface::define_marker_models() const
{
  string source = "";
  
  source += "# - marker is a tuple of two strings representing the\n";
  source += "#  alleles for an individual at a particular locus.\n";
  source += "#  allele is a string representing the name of the allele\n";
  source += "#  that the model is based on.\n";
  source += "#      \n";
  source += "def dominant(marker, allele):\n";
  source += "  return (marker[0] == allele) or (marker[1] == allele)\n";
  source += "def recessive(marker, allele):\n";
  source += "  return (marker[0] == allele) and (marker[1] == allele)\n";
  source += "def additive(marker, allele):\n";
  source += "  return (marker[0] == allele) + (marker[1] == allele)\n";
  
  source += "# - alternative names for the preceding three functions.\n";
  source += "#      \n";
  source += "def dom(marker, allele):\n";
  source += "  return dominant(marker, allele)\n";
  source += "def rec(marker, allele):\n";
  source += "  return recessive(marker, allele)\n";
  source += "def add(marker, allele):\n";
  source += "  return additive(marker, allele)\n";
  
  source += "# - Does the value of the genotype match the two alleles?\n";
  source += "#      \n";
  source += "def genotype(marker, allele1, allele2):\n";
  source += "  return ((marker[0] == allele1 and marker[1] == allele2) or\n";
  source += "          (marker[0] == allele2 and marker[1] == allele1)   )\n";
  source += "def gen(marker, allele1, allele2):\n";
  source += "  return genotype(marker, allele1, allele2)\n";

  runPythonString(source);
}

//===============================================
//
//  extract_names(...)
//
//===============================================
bool                
PythonInterface::extract_names(const string & expr) 
{
  if(!addStringToEnvironment("expression", expr))
    return false;;

  std::string source_code = "";

  source_code += "ast   = parser.expr(expression)\n";
  source_code += "tree  = ast.totuple()\n";
  source_code += "names = extract_names(tree).keys()";
  
  PyObject * result = PyRun_String(const_cast<char*>(source_code.c_str()), Py_file_input, my_globals, my_globals);
  
  if(result)
  {
    Py_DECREF(result);

    return true;
  }
  else 
  {
    my_errors << py_error_msg() << "\n";

    return false;
  }
}

void
PythonInterface::exit_on_timeout(const string& expr)
{
  my_errors << priority(fatal) 
            << "Timed out evaluating expression, '" 
            << expr
            << "'.  Change expression, increase time limit attribute"
            << " of the function block parameter, or run program when your" 
            << " machine is not so heavily loaded.  Program exiting ..." 
            << endl;

  exit(TIMED_OUT);
}

} // End namespace FUNC
} // End namespace SAGE
