casa_add_library( crash
    utils/CrashReportHelper.cc
)

#if (UseCrashReporter)#
#    casa_add_executable (NONE CrashReportPoster apps/reporter/CrashReportPoster.cc ${Breakpad_IncludeRoot}/common/linux/http_upload.cc utils/CrashReportHelper)
#    target_link_libraries (CrashReportPoster CrashReportHelper  -ldl)
#    target_include_directories (CrashReportPoster PUBLIC ${Breakpad_IncludeRoot})
#endif ()

if (UseCrashReporter)#
    casa_add_executable (NONE CrashReportPoster apps/reporter/CrashReportPoster.cc ${Breakpad_IncludeRoot}/common/linux/http_upload.cc utils/CrashReportHelper)
    target_link_libraries (CrashReportPoster -ldl)
    target_include_directories (CrashReportPoster PUBLIC ${Breakpad_IncludeRoot})
endif ()
