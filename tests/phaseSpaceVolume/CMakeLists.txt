### BlockD

add_executable(test_psVol_blockD "blockD.cc")
target_link_libraries(test_psVol_blockD momemta)
set_target_properties(test_psVol_blockD PROPERTIES OUTPUT_NAME
  "phaseSpaceVolume_blockD.exe")

### BlockF

add_executable(test_psVol_blockF "blockF.cc")
target_link_libraries(test_psVol_blockF momemta)
set_target_properties(test_psVol_blockF PROPERTIES OUTPUT_NAME
  "phaseSpaceVolume_blockF.exe")

### BlockA

add_executable(test_psVol_blockA "blockA.cc")
target_link_libraries(test_psVol_blockA momemta)
set_target_properties(test_psVol_blockA PROPERTIES OUTPUT_NAME
  "phaseSpaceVolume_blockA.exe")
