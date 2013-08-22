cmake_minimum_required(VERSION 2.8.9)


#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
set(EXTENSION_NAME SurfaceReconstruction)
set(EXTENSION_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/SurfaceReconstruction")
set(EXTENSION_CATEGORY "Surface Models")
set(EXTENSION_CONTRIBUTORS "David-Doria Arnaud-Gelas, Jean-Christophe Fillion-Robin (Kitware)")
set(EXTENSION_DESCRIPTION "This is an implementation of the Poisson surface reconstruction algorithm in the VTK framework")
set(EXTENSION_ICONURL "http://viewvc.slicer.org/viewvc.cgi/Slicer4/trunk/Extensions/Testing/SurfaceReconstruction/SurfaceReconstruction.png?revision=21745&view=co")
set(EXTENSION_SCREENSHOTURLS "http://wiki.slicer.org/slicerWiki/images/a/ab/Slicer-r19441-SurfaceReconstruction-screenshot.png http://wiki.slicer.org/slicerWiki/images/1/1e/Slicer-r19441-SurfaceReconstruction-screenshot-2.png")
endif()

#-----------------------------------------------------------------------------
set(MODULE_NAME SurfaceReconstruction)


#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()

#-----------------------------------------------------------------------------
set(MODULE_TARGET_LIBRARIES
  ${ITK_LIBRARIES}
  ${VTK_LIBRARIES}
  )

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  ADDITIONAL_SRCS
	vtkPoissonReconstruction.cxx 
	Factor.cpp 
	Geometry.cpp 
	MarchingCubes.cpp 
	MultiGridOctest.cpp 
	MultiGridOctreeData.cpp
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  include(${Slicer_EXTENSION_CPACK})
endif()