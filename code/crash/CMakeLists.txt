install (FILES
    utils/CrashReportHelper.h
    DESTINATION include/casacode/crash/utils
)

if (UseCrashReporter AND NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    casa_add_executable (NONE CrashReportPoster apps/reporter/CrashReportPoster.cc ${Breakpad_IncludeRoot}/common/linux/http_upload.cc)
    target_link_libraries (CrashReportPoster -ldl)
    target_include_directories (CrashReportPoster PUBLIC ${Breakpad_IncludeRoot} ${PROJECT_SOURCE_DIR}/crash/utils)
endif ()
