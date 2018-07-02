#include <stdcasa/optionparser.h>
#include <string>
#include <iostream>
#include <cstdlib>
#include <climits>

namespace wvr {

    namespace arg {
        namespace index {
            enum optionIndex { UNKNOWN, HELP, MS, OUTPUT, TOFFSET, NSOL, SEGFIELD,
                               SEGSOURCE, REVERSE, REVERSESPW, DISPERSE, CONT, WVRFLAG,
                               SOURCEFLAG, STATFIELD, STATSOURCE, TIE, SMOOTH,
                               SCALE, MAXDISTM, MINNUMANTS, MINGOODFRAC, USEFIELDTAB,
                               SPW, WVRSPW, REFANT, OFFSETS };
        }

        struct check : public option::Arg {
            static void printError(const char* msg1, const option::Option& opt, const char* msg2) {
                std::cerr << msg1 << std::string(opt.name,opt.namelen) << msg2 << std::endl;
            }
            static option::ArgStatus Unknown(const option::Option& option, bool msg) {
                if (msg) printError("Unknown option '", option, "'\n");
                return option::ARG_ILLEGAL;
            }

            // an associated value is required and must be a valid floating point (double)
            static option::ArgStatus Float(const option::Option& option, bool msg) {
                char* endptr = 0;
                if (option.arg != 0 && strtod(option.arg, &endptr)){};
                if (endptr != option.arg && *endptr == 0)
                    return option::ARG_OK;
      
                if (msg) printError("Option '", option, "' requires a floating point argument\n");
                return option::ARG_ILLEGAL;
            }

            // an associated value is required and must be an integer
            static option::ArgStatus Int(const option::Option& option, bool msg) {
                char* endptr = 0;
                if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
                if (endptr != option.arg && *endptr == 0)
                    return option::ARG_OK;
      
                if (msg) printError("Option '", option, "' requires an integer argument\n");
                return option::ARG_ILLEGAL;
            }

            // an associated value is required opt=arg, no constraints on type
            static option::ArgStatus Required(const option::Option& option, bool msg) {
                if (option.arg != 0)
                    return option::ARG_OK;
      
                if (msg) printError("Option '", option, "' requires an argument\n");
                return option::ARG_ILLEGAL;
            }

            // option that is no longer supported
            static option::ArgStatus Deprecated(const option::Option& option, bool msg) {
                printError("Option '", option, "' is no longer valid, see options with --help\n");
                return option::ARG_ILLEGAL;	
            }
        };
    }
}
