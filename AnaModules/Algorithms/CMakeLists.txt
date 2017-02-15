art_make( BASENAME_ONLY
          LIBRARY_NAME UBooneDevAlgs
          LIB_LIBRARIES 
			lardataobj_RawData
                        lardataobj_RecoBase
			lardata_RecoObjects
			lardata_ArtDataHelper
                        lardata_Utilities
                        larcore_Geometry
	  		larcore_Geometry_Geometry_service
			lardata_DetectorInfoServices_DetectorPropertiesServiceStandard_service
			larevt_Filters
			larevt_CalibrationDBI_IOVData
			larevt_CalibrationDBI_Providers
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
			${ART_UTILITIES}
                        ${FHICLCPP}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
                        ${FHICLCPP}
			${CETLIB}
                        ${Boost_SYSTEM_LIBRARY}
	)

install_headers()
install_fhicl()
install_source()
