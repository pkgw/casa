#
# CASA - Common Astronomy Software Applications
# Copyright (C) 2010 by ESO (in the framework of the ALMA collaboration)
#
# This file is part of CASA.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

casa_add_library( mstransform  
 MSTransform/MSTransform.cc
 MSTransform/MSTransformManager.cc
 MSTransform/MSTransformDataHandler.cc
 MSTransform/MSTransformRegridder.cc
 MSTransform/MSTransformBufferImpl.cc
 MSTransform/MSTransformIterator.cc
 MSTransform/MSTransformIteratorFactory.cc
 MSTransform/MSUVBin.cc
 MSTransform/StatWt.cc
 TVI/UtilsTVI.cc 
 TVI/DenoisingLib.cc 
 TVI/FreqAxisTVI.cc 
 TVI/ChannelAverageTVI.cc
 TVI/ConvolutionTVI.cc
 TVI/HanningSmoothTVI.cc	 
 TVI/PhaseShiftingTVI.cc 
 TVI/RegriddingTVI.cc  
 TVI/UVContSubTVI.cc
 TVI/PolAverageTVI.cc
 TVI/StatWtTVI.cc
 TVI/StatWtTVIFactory.cc
 TVI/StatWtTVILayerFactory.cc
)

install (FILES
    	MSTransform/MSTransform.h
	MSTransform/MSTransformManager.h
	MSTransform/MSTransformDataHandler.h
	MSTransform/MSTransformRegridder.h
	MSTransform/MSTransformBufferImpl.h
	MSTransform/MSTransformIterator.h
	MSTransform/MSTransformIteratorFactory.h
	MSTransform/MSUVBin.h
	MSTransform/StatWt.h
	DESTINATION include/casacode/mstransform/MSTransform
	)
	
install (FILES
 	MSTransform/test/tMSTransformIterator.h	
	DESTINATION include/casacode/mstransform/MSTransform/test
	)		
	
install (FILES
 	TVI/UtilsTVI.h	
 	TVI/DenoisingLib.h	
 	TVI/FreqAxisTVI.h	
	TVI/ChannelAverageTVI.h	
	TVI/ConvolutionTVI.h	
	TVI/HanningSmoothTVI.h
	TVI/PhaseShiftingTVI.h		
 	TVI/RegriddingTVI.h
 	TVI/UVContSubTVI.h 	
 	TVI/PolAverageTVI.h
    TVI/StatWtTVI.h
    TVI/StatWtTVIFactory.h
    TVI/StatWtTVILayerFactory.h
	DESTINATION include/casacode/mstransform/TVI
	)
	
install (FILES
 	TVI/test/TestUtilsTVI.h	
	DESTINATION include/casacode/mstransform/TVI/test
	)	

casa_add_assay( mstransform MSTransform/test/tMSBin.cc )
#casa_add_assay( mstransform MSTransform/test/tStatWt.cc )
casa_add_executable( mstransform msuvbin apps/msuvbin/msuvbin.cc )
casa_add_executable( mstransform fixspwbackport apps/fixspwbackport/fixspwbackport.cc )

casa_add_google_test( MODULES mstransform SOURCES TVI/test/tRecursiveVi2Layers_GT.cc )
casa_add_google_test( MODULES mstransform SOURCES TVI/test/PolAverageTVI_GTest.cc )
#casa_add_google_test( MODULES mstransform SOURCES TVI/test/tChannelAverageTransformEngine_GT.cc )
casa_add_google_test( MODULES mstransform SOURCES TVI/test/tStatWtTVI.cc  TVI/test/TestUtilsTVI.cc)
