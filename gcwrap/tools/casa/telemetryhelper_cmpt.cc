#include "crash/utils/CrashReportHelper.h"
#include <telemetryhelper_cmpt.h>
namespace casac {
    telemetryhelper::telemetryhelper()
    {

    }

    telemetryhelper::~telemetryhelper()
    {

    }
    //typedef std::vector<unsigned> tmac;

    std::string telemetryhelper::getUniqueId () {
        CrashReportHelper ch;
        return ch.getUniqueId();
    }
}
