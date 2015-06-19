# basic source code CMakeLists.txt

##include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

##include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{CGAL_DIR}/include )
find_library( CGAL_LIBRARY      NAMES CGAL      PATH $ENV{CGAL_DIR}/lib )
find_library( CGAL_CORE_LIBRARY NAMES CGAL_Core PATH $ENV{CGAL_DIR}/lib )
find_library( GMP_LIBRARY       NAMES gmp       PATH $ENV{CGAL_DIR}/lib )
find_library( MPFR_LIBRARY      NAMES mpfr      PATH $ENV{CGAL_DIR}/lib )

art_make( BASENAME_ONLY
          MODULE_LIBRARIES  RecoObjects
			    RecoBase
			    MCBase
			    Geometry
			    Geometry_service
                            LArProperties_service
			    DetectorProperties_service
			    Utilities
			    ${SIMULATIONBASE}
			    ${ART_FRAMEWORK_CORE}
			    ${ART_FRAMEWORK_PRINCIPAL}
			    ${ART_FRAMEWORK_SERVICES_REGISTRY}
			    ${ART_FRAMEWORK_SERVICES_BASIC}
			    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			    ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			    ${ART_PERSISTENCY_COMMON}
			    ${ART_UTILITIES}
			    ${MF_MESSAGELOGGER}
			    ${MF_UTILITIES}
			    ${FHICLCPP}
			    ${CETLIB}
			    ${CLHEP}
			    ${ROOT_BASIC_LIB_LIST}
			    ${G4_LIB_LIST}
			    ${CGAL_LIBRARY}
			    ${CGAL_CORE_LIBRARY}
			    ${GMP_LIBRARY}
			    ${MPFR_LIBRARY}
          )

install_headers()
install_source()
install_fhicl()
