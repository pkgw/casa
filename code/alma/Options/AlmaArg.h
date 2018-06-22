

#include <stdcasa/optionparser.h>
#include <string>
#include <iostream>
#include <cstdlib>
#include <climits>
#include <alma/ASDM/Misc.h>

#if !defined(ALMAARG_H)

namespace alma {

  using namespace asdm;

  // extend the option::Arg struct to provide the necessary argument checking for the alma apps
  // for use with the optionparser suite

  struct AlmaArg: public option::Arg
  {
    static void printError(const char* msg1, const option::Option& opt, const char* msg2)
    {
      std::cerr << msg1 << std::string(opt.name,opt.namelen) << msg2 << std::endl;
    }
    
    // unknown argument
    static option::ArgStatus Unknown(const option::Option& option, bool msg)
    {
      if (msg) printError("Unknown option '", option, "'\n");
      return option::ARG_ILLEGAL;
    }
    
    // an associated value is required opt=arg, no constraints on type
    static option::ArgStatus Required(const option::Option& option, bool msg)
    {
      if (option.arg != 0)
	return option::ARG_OK;
      
      if (msg) printError("Option '", option, "' requires an argument\n");
      return option::ARG_ILLEGAL;
    }
    
    // an associated value is required and must be an a long integer
    static option::ArgStatus Long(const option::Option& option, bool msg)
    {
      char* endptr = 0;
      if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
      if (endptr != option.arg && *endptr == 0)
	return option::ARG_OK;
      
      if (msg) printError("Option '", option, "' requires a long integer argument\n");
      return option::ARG_ILLEGAL;
    }
    
    // an associated value is required and must be an unsigned int
    static option::ArgStatus uInt(const option::Option& option, bool msg)
    {
      char* endptr = 0;
      unsigned long ulval = 0;
      // uses strtoul, but must also check that converted value would fit in an unsigned int
      if (option.arg != 0 && (ulval=strtoul(option.arg, &endptr, 10))){};
      if (endptr != option.arg && *endptr == 0 && ulval <= UINT_MAX) 
	return option::ARG_OK;
      
      if (msg) printError("Option '", option, "' requires an unsigned integer argument\n");
      return option::ARG_ILLEGAL;	
    }
    
    // an associated value is required and must be a valid floating point (double)
    static option::ArgStatus Float(const option::Option& option, bool msg)
    {
      char* endptr = 0;
      if (option.arg != 0 && strtod(option.arg, &endptr)){};
      if (endptr != option.arg && *endptr == 0)
	return option::ARG_OK;
      
      if (msg) printError("Option '", option, "' requires a floating point argument\n");
      return option::ARG_ILLEGAL;
    }
    
    // an associated value is required and must be either "true" or "false"
    // ignores case
    static option::ArgStatus Bool(const option::Option& option, bool msg)
    {
      if (option.arg != 0) {
	std::string argVal(option.arg);
	trim(argVal);
	argVal = str_tolower(argVal);
	if (argVal == "true" || argVal == "false") 
	  return option::ARG_OK;
      }
      if (msg) printError("Option '", option, "' requires a value of 'true' or 'false'\n");
      return option::ARG_ILLEGAL;
    }
  };
  
}

#define ALMAARG_H
#endif
