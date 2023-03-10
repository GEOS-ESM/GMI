esma_set_this ()

set (src_directories
  .
  GMI_GridComp
  GMI_GridComp/GmiChemistry
  GMI_GridComp/GmiChemistry/AerosolDust
  GMI_GridComp/GmiChemistry/ioChemistry
  GMI_GridComp/GmiChemistry/photolysis/CloudJ
  GMI_GridComp/GmiChemistry/photolysis/fastJX53c_ref
  GMI_GridComp/GmiChemistry/photolysis/fast_JX
  GMI_GridComp/GmiChemistry/photolysis/fast_JX65
  GMI_GridComp/GmiChemistry/photolysis/fastj
  GMI_GridComp/GmiChemistry/photolysis/fast_JX53b
  GMI_GridComp/GmiChemistry/photolysis/fast_JX53c
  GMI_GridComp/GmiChemistry/photolysis/lookup
  GMI_GridComp/GmiChemistry/photolysis/utils
  GMI_GridComp/GmiChemistry/sad
  GMI_GridComp/GmiChemistry/smv2chem
  GMI_GridComp/GmiDeposition
  GMI_GridComp/GmiEmission
  GMI_GridComp/GmiEmission/GOCARTroutines
  GMI_GridComp/GmiEmission/gsfc
  GMI_GridComp/GmiEmission/harvard
  GMI_GridComp/GmiEmission/ioEmission
  GMI_GridComp/GmiEmission/lightning
  GMI_GridComp/GmiEmission/llnl
  GMI_GridComp/GmiShared/GmiESMF
  GMI_GridComp/GmiShared/GmiIOutilities
  GMI_GridComp/GmiShared/GmiSupportingModules
  GMI_GridComp/GmiSpeciesConcentration/spcConcentrationMethod
  )

# Set (i.e. hardcode in) the default GMI mechanism
set (GMI_MECHANISM "StratTrop_Orig" CACHE STRING "GMI Mechanism chosen by the user at CMake configure time (default: StratTrop_Orig)")
# Compile list of possible mechanisms
set (mechdir ${CMAKE_CURRENT_SOURCE_DIR}/GMI_GridComp/GmiChemistry)

# Glob through the files relative to the mechdir
file (GLOB_RECURSE mechfiles CONFIGURE_DEPENDS RELATIVE ${mechdir} setkin_chem_mech.txt)

# This will return a list like: StratTrop_HFC_S/setkin_chem_mech.txt;StratTrop_Orig/setkin_chem_mech.txt
# Now, loop through mechfiles and extract the directory names
set (ALLOWED_MECHANISM)
foreach (mechfile ${mechfiles})
  get_filename_component (mechdir ${mechfile} DIRECTORY)
  list (APPEND ALLOWED_MECHANISM ${mechdir})
endforeach ()

if (GMI_MECHANISM IN_LIST ALLOWED_MECHANISM)
  set (MECHANISM_DIR "${GMI_MECHANISM}")
else ()
  message(FATAL_ERROR "GMI_MECHANISM must be one of ${ALLOWED_MECHANISM}, is set to: ${GMI_MECHANISM}")
endif()
message(STATUS "GMI Mechanism chosen: ${MECHANISM_DIR}")
list (APPEND src_directories GMI_GridComp/GmiChemistry/${MECHANISM_DIR})

set (srcs)
foreach (dir ${src_directories})
  file (GLOB tmpsrcs CONFIGURE_DEPENDS ${dir}/*.[fF] ${dir}/*.[fF]90 ${dir}/*.c)
  list (APPEND srcs ${tmpsrcs})
endforeach()

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   if (CMAKE_BUILD_TYPE MATCHES Release)
      set (CMAKE_Fortran_FLAGS_RELEASE "")
      foreach (src ${srcs})
         if (${src} MATCHES GMIchem_GridCompMod.F90)
            set_source_files_properties (${src} PROPERTIES COMPILE_FLAGS "${FOPT2} ${GEOS_Fortran_Release_FPE_Flags} ${ALIGNCOM}")
         elseif (${src} MATCHES GmiChem_GridCompClassMod.F90)
            set_source_files_properties (${src} PROPERTIES COMPILE_FLAGS "${FOPT1} ${GEOS_Fortran_Release_FPE_Flags} ${ALIGNCOM}")
         else ()
            set_source_files_properties (${src} PROPERTIES COMPILE_FLAGS "${GEOS_Fortran_FLAGS_RELEASE} ${EXTENDED_SOURCE}")
         endif ()
      endforeach ()
   endif ()
endif ()
if (CMAKE_BUILD_TYPE MATCHES Release)
   set_source_files_properties( GMI_GridComp/GmiChem_GridCompMod.F90
      PROPERTIES COMPILE_FLAGS "-O2 ${common_Fortran_flags} ${GEOS_Fortran_Release_FPE_Flags} ${ALIGNCOM}")
   set_source_files_properties( GMI_GridComp/GmiChem_GridCompClassMod.F90
      PROPERTIES COMPILE_FLAGS "-O1 ${common_Fortran_flags} ${GEOS_Fortran_Release_FPE_Flags} ${ALIGNCOM}")
endif ()



esma_add_library (${this} SRCS ${srcs} DEPENDENCIES Chem_Shared GEOS_Shared esmf Process_Library)

target_include_directories (${this}
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/GMI_GridComp>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/GMI_GridComp/GmiChemistry/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/GMI_GridComp/GmiChemistry/photolysis/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/GMI_GridComp/GmiChemistry/${GMI_MECHANISM}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/GMI_GridComp/GmiShared/GmiInclude>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/GMI_GridComp/GmiChemistry/photolysis/CloudJ>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/GMI_GridComp/GmiChemistry/photolysis/fast_JX>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/GMI_GridComp/GmiChemistry/photolysis/fast_JX65>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/GMI_GridComp/GmiChemistry/photolysis/fast_JX53b>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/GMI_GridComp/GmiChemistry/photolysis/fastJX53c_ref>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/GMI_GridComp/GmiChemistry/photolysis/fastj>
  )

target_include_directories (${this} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/GMI_GridComp>)
target_include_directories (${this} PUBLIC ${INC_ESMF})



set (gmi_gc ${CMAKE_CURRENT_SOURCE_DIR}/GMI_GridComp)
set (registry ${gmi_gc}/GmiChemistry/${GMI_MECHANISM}/GMI_Registry.rc)
set (mechanism ${gmi_gc}/GmiChemistry/${GMI_MECHANISM}/setkin_lchem.h)
set (deposParam ${gmi_gc}/GmiChemistry/${GMI_MECHANISM}/setkin_depos.h)
set (mechParam ${gmi_gc}/GmiChemistry/${GMI_MECHANISM}/setkin_par.h)
set (mechRC ${gmi_gc}/GmiChemistry/${GMI_MECHANISM}/setkin_chem_mech.txt)
set (first_stage_files Deposition_Registry___.rc Reactions_Registry___.rc GMI_Tendency_Registry___.rc setkin_chem_mech.txt___.rc GMICHEM_History___.rc)

find_file (generator
  NAME mapl_acg.pl
  PATHS ${esma_include}/MAPL_Base ${esma_etc}/MAPL)

set (name GMICHEM)

add_custom_command (
  OUTPUT ${first_stage_files}
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/gmi_acg.pl -R -r ${registry} -m ${mechanism} -d ${deposParam} -p ${mechParam}
  COMMAND ${CMAKE_COMMAND} -E copy  ${mechRC} setkin_chem_mech.txt___.rc

  COMMAND ${generator} ${acg_flags} -N ${name} Deposition_Registry___.rc
  COMMAND ${CMAKE_COMMAND} -E copy ${name}_DeclarePointer___.h GMI_GridComp/Deposition_DeclarePointer___.h
  COMMAND ${CMAKE_COMMAND} -E copy ${name}_GetPointer___.h GMI_GridComp/Deposition_GetPointer___.h
  COMMAND ${CMAKE_COMMAND} -E copy ${name}_ExportSpec___.h Deposition_ExportSpec___.h

  COMMAND ${generator} ${acg_flags} -N ${name} Reactions_Registry___.rc
  COMMAND ${CMAKE_COMMAND} -E copy ${name}_DeclarePointer___.h GMI_GridComp/Reactions_DeclarePointer___.h
  COMMAND ${CMAKE_COMMAND} -E copy ${name}_GetPointer___.h GMI_GridComp/Reactions_GetPointer___.h
  COMMAND ${CMAKE_COMMAND} -E copy ${name}_ExportSpec___.h Reactions_ExportSpec___.h

  COMMAND ${generator} ${acg_flags} -N ${name} GMI_Tendency_Registry___.rc
  COMMAND ${CMAKE_COMMAND} -E copy ${name}_DeclarePointer___.h GMI_GridComp/Tendency_DeclarePointer___.h
  COMMAND ${CMAKE_COMMAND} -E copy ${name}_GetPointer___.h GMI_GridComp/Tendency_GetPointer___.h
  COMMAND ${CMAKE_COMMAND} -E copy ${name}_ExportSpec___.h Tendency_ExportSpec___.h

  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/gmi_acg.pl -v -r ${registry} -m ${mechanism} -d ${deposParam} -p ${mechParam}
  COMMAND ${CMAKE_COMMAND} -E copy GMI_GridComp/Deposition_GetPointer2___.h GMI_GridComp/Deposition_GetPointer___.h
  COMMAND ${CMAKE_COMMAND} -E copy GMI_GridComp/Reactions_GetPointer2___.h GMI_GridComp/Reactions_GetPointer___.h

  COMMAND ${generator} ${acg_flags} -N ${name} ${registry}
  COMMAND ${CMAKE_COMMAND} -E copy ${name}_DeclarePointer___.h GMI_GridComp/Tendency_DeclarePointer___.h
  COMMAND ${CMAKE_COMMAND} -E copy ${name}_History___.rc GMI_GridComp/${name}_History___.rc

  MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/gmi_acg.pl
  DEPENDS ${registry} ${mechanism} ${deposParam} ${mechParam} ${mechRC}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "First stage generation of automatic code in GMI"
  )

add_custom_target (phony_gmi DEPENDS ${first_stage_files})
add_dependencies (${this} phony_gmi)

#add_custom_command (
#  OUTPUT something
#  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/gmi_acg.pl -v
#  MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/gmi_acg.pl
#  DEPENDS
#  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#  COMMENT "2nd stage generation of automatic code in GMI"
#  )

set ( rc_files
   GMI_GridComp/GMI_ExtData.rc
   GMI_GridComp/GMI_ExtData.yaml
   GMI_GridComp/GmiChemistry/${GMI_MECHANISM}/GMI_GridComp.rc
   GMI_GridComp/GmiChemistry/${GMI_MECHANISM}/GMI_Registry.rc
   GMI_GridComp/GmiChemistry/${GMI_MECHANISM}/GMI_Mech_Registry.rc
)
install (
   FILES ${rc_files}
   DESTINATION etc
   )

foreach (file ${first_stage_files})
   install(
      FILES ${CMAKE_CURRENT_BINARY_DIR}/${file}
      DESTINATION etc
      )
endforeach ()

