set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-underscoring -fPIC")

set(SOURCES
    "src/VVSS3_0.f"
    "src/VVSS2_0.f"
    "src/VVS2_3.f"
    "src/HelAmps_SMEFT_FF_2.cc"
    "src/Parameters_SMEFT_FF_2.cc"
    "SubProcesses/P1_Sigma_SMEFT_FF_2_gg_bbxbbx/P1_Sigma_SMEFT_FF_2_gg_bbxbbx.cc"
    )

add_library(me_pp_hh_bbbb OBJECT ${SOURCES}) 

target_include_directories(me_pp_hh_bbbb PRIVATE "include")
target_include_directories(me_pp_hh_bbbb PRIVATE "SubProcesses/P1_Sigma_SMEFT_FF_2_gg_bbxbbx/")
