# basic source code CMakeLists.txt
include_directories ( ${CMAKE_CURRENT_SOURCE_DIR} )

FILE( GLOB src_files *.cxx )

SET(CET_DICT_FUNCTIONS TRUE)


art_make( LIBRARY_NAME  EventDisplayLibs
          LIB_LIBRARIES ColorDrawingOptions_service
	                InfoTransfer_service
	                Filters
                        LArProperties_service
                        DetectorProperties_service
 			Simulation
			RawData
			RecoBase
			RecoObjects
			Utilities
			AnalysisBase
			AnalysisAlg
		 	CalibrationDBI_IOVData
			CalibrationDBI_WebDBI
	                ${SIMULATIONBASE}
			${EVENTDISPLAYBASE}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_PERSISTENCY_COMMON}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
                        ${ROOT_EVE_LIB_LIST}
                        ${ROOT_EVE}
                        ${ROOT_EG}
                        ${ROOT_RGL}
			${ROOT_X3d}
           SERVICE_LIBRARIES
			    EventDisplayLibs
    			    ${ART_FRAMEWORK_CORE}
   			    ${ART_UTILITIES}
    			    ${ART_FRAMEWORK_IO_ROOT_ROOTINPUT_SOURCE}
   			    ${ART_FRAMEWORK_IO_ROOT}
   			    ${ART_PERSISTENCY_PROVENANCE}
    			    ${MF_MESSAGELOGGER}
   			    ${ROOT_EVE}
	   MODULE_LIBRARIES EventDisplayLibs
	                    InfoTransfer_service
	                    Filters
                            LArProperties_service
                            DetectorProperties_service
 			    Simulation
			    RawData
			    RecoBase
			    RecoObjects
			    Utilities
			    AnalysisBase
			    AnalysisAlg
	                    ${SIMULATIONBASE}
			    ${EVENTDISPLAYBASE}
                            ${ART_FRAMEWORK_CORE}
			    ${ART_FRAMEWORK_PRINCIPAL}
			    ${ART_PERSISTENCY_COMMON}
			    ${ART_UTILITIES}
			    ${MF_MESSAGELOGGER}
			    ${MF_UTILITIES}
	                    ${FHICLCPP}
			    ${CETLIB}
                            ${ROOT_GUI}
                            ${ROOT_GEOM}
                            ${ROOT_EVE}
                            ${ROOT_EG}
                            ${ROOT_RGL}
                )

install_headers()
install_source()
install_fhicl()
