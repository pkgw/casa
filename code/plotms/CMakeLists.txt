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

include( ${QT_USE_FILE} )
set( QT4_LIBRARIES ${QT_LIBRARIES} )

qt4_wrap_cpp( outfiles_qoh
  Gui/PlotMSAveragingWidget.qo.h
  Gui/PlotMSAxisWidget.qo.h
  Gui/PlotMSCalibrationWidget.qo.h
  Gui/PlotMSFlaggingWidget.qo.h
  Gui/PlotMSLoggerWidget.qo.h
  Gui/PlotMSPlotter.qo.h
  Gui/PlotMSSelectionWidget.qo.h
  Gui/PlotMSSummaryDialog.qo.h
  Gui/PlotMSTransformationsWidget.qo.h
  Gui/PlotRangeWidget.qo.h
  Gui/PlotMSPageHeaderDataModel.qo.h
  GuiTabs/PlotMSAnnotatorTab.qo.h
  GuiTabs/PlotMSAxesTab.qo.h
  GuiTabs/PlotMSCacheTab.qo.h
  GuiTabs/PlotMSCalibrationTab.qo.h
  GuiTabs/PlotMSCanvasTab.qo.h
  GuiTabs/PlotMSDataCollapsible.qo.h
  GuiTabs/PlotMSDataSummaryTab.qo.h
  GuiTabs/PlotMSDataSymbolWidget.qo.h
  GuiTabs/PlotMSDataTab.qo.h
  GuiTabs/PlotMSDisplayTab.qo.h
  GuiTabs/PlotMSIterateTab.qo.h
  GuiTabs/PlotMSExportTab.qo.h
  GuiTabs/PlotMSFlaggingTab.qo.h
  GuiTabs/PlotMSOptionsTab.qo.h
  GuiTabs/PlotMSPlotTab.qo.h
  GuiTabs/PlotMSTab.qo.h
  GuiTabs/PlotMSToolsTab.qo.h
  GuiTabs/PlotMSTransformationsTab.qo.h
  Threads/Gui/PlotMSCacheThread.qo.h
  Threads/Gui/PlotMSDrawThread.qo.h
  Threads/Gui/PlotMSExportThread.qo.h
  Threads/Gui/PlotMSThread.qo.h
  )

casa_qt4_wrap_ui( outfiles_ui
  Gui/PlotMSAveragingWidget.ui
  Gui/PlotMSAxisWidget.ui
  Gui/PlotMSCalibrationWidget.ui
  Gui/PlotMSFlaggingWidget.ui
  Gui/PlotMSPlotter.ui
  Gui/PlotMSSummaryDialog.ui
  Gui/PlotMSTransformationsWidget.ui
  Gui/PlotRangeWidget.ui
  GuiTabs/PlotMSAnnotatorTab.ui
  GuiTabs/PlotMSAxesTab.ui
  GuiTabs/PlotMSCacheTab.ui
  GuiTabs/PlotMSCalibrationTab.ui
  GuiTabs/PlotMSCanvasTab.ui
  GuiTabs/PlotMSDataTab.ui
  GuiTabs/PlotMSDataCollapsible.ui
  GuiTabs/PlotMSDataSummaryTab.ui
  GuiTabs/PlotMSDataSymbolWidget.ui
  GuiTabs/PlotMSDisplayTab.ui
  GuiTabs/PlotMSIterateTab.ui
  GuiTabs/PlotMSExportTab.ui
  GuiTabs/PlotMSFlaggingTab.ui
  GuiTabs/PlotMSOptionsTab.ui
  GuiTabs/PlotMSPlotTab.ui
  GuiTabs/PlotMSToolsTab.ui
  GuiTabs/PlotMSTransformationsTab.ui
  )

casa_qt4_add_resources( outfiles_qrc Gui/PlotMSPlotter.qrc OPTIONS -name Gui )
casa_qt4_add_resources( outfiles_qrc GuiTabs/GuiTabs.qrc OPTIONS -name GuiTabs )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/Gui ) 
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/GuiTabs ) 
# because the generated code contain relative paths to the source tree

casa_add_library( plotms
  ${outfiles_qoh} ${outfiles_ui} ${outfiles_qrc}
  Actions/ActionAnnotate.cc
  Actions/ActionAnnotateRectangle.cc
  Actions/ActionAnnotateText.cc
  Actions/ActionCache.cc
  Actions/ActionCacheLoad.cc
  Actions/ActionCacheRelease.cc
  Actions/ActionClearPlotter.cc
  Actions/ActionClearRegions.cc
  Actions/ActionExport.cc
  Actions/ActionExportDialog.cc
  Actions/ActionFactory.cc
  Actions/ActionFlagAll.cc
  Actions/ActionHoldReleaseDrawing.cc
  Actions/ActionInformation.cc
  Actions/ActionIterFirst.cc
  Actions/ActionIterLast.cc
  Actions/ActionIterNext.cc
  Actions/ActionIterPrev.cc
  Actions/ActionLocate.cc
  Actions/ActionMarkRegions.cc
  Actions/ActionPan.cc
  Actions/ActionPlot.cc
  Actions/ActionQuit.cc
  Actions/ActionSelect.cc
  Actions/ActionSelectFlag.cc
  Actions/ActionSelectUnflag.cc
  Actions/ActionStack.cc
  Actions/ActionStackBack.cc
  Actions/ActionStackBase.cc
  Actions/ActionStackForward.cc
  Actions/ActionSubtractRegions.cc
  Actions/ActionSummary.cc
  Actions/ActionSummaryDialog.cc
  Actions/ActionTool.cc
  Actions/ActionTracker.cc
  Actions/ActionTrackerDisplay.cc
  Actions/ActionTrackerHover.cc
  Actions/ActionZoom.cc
  Actions/FlagActionUtil.cc
  Actions/PlotMSAction.cc
  Threads/BackgroundThread.cc
  Threads/CacheThread.cc
  Threads/ExportThread.cc
  Threads/Gui/PlotMSCacheThread.cc
  Threads/Gui/PlotMSDrawThread.cc
  Threads/Gui/PlotMSExportThread.cc
  Threads/Gui/PlotMSThread.cc
  Threads/ThreadController.cc
  Client/ClientFactory.cc
  Client/ClientScript.cc
  Client/Client.cc
  Data/PageHeaderCache.cc
  Data/PlotMSCacheBase.cc
  Data/MSCache.cc
  Data/MSCacheVolMeter.cc
  Data/CalCache.cc
  Data/PlotMSAtm.cc
  Data/PlotMSIndexer.cc
  Data/PlotMSVBAverager.cc
  Gui/PlotMSAnnotator.cc
  Gui/PlotMSAveragingWidget.cc
  Gui/PlotMSAxisWidget.cc
  Gui/PlotMSCalibrationWidget.cc
  Gui/PlotMSFlaggingWidget.cc
  Gui/PlotMSLoggerWidget.cc
  Gui/PlotMSPlotter.cc
  Gui/PlotMSSelectionWidget.cc
  Gui/PlotMSSummaryDialog.cc
  Gui/PlotMSTransformationsWidget.cc
  Gui/PlotRangeWidget.cc
  Gui/PlotMSPageHeaderDataModel.cc
  GuiTabs/PlotMSAnnotatorTab.cc
  GuiTabs/PlotMSAxesTab.cc
  GuiTabs/PlotMSCacheTab.cc
  GuiTabs/PlotMSCalibrationTab.cc
  GuiTabs/PlotMSCanvasTab.cc
  GuiTabs/PlotMSDataCollapsible.cc
  GuiTabs/PlotMSDataSummaryTab.cc
  GuiTabs/PlotMSDataSymbolWidget.cc
  GuiTabs/PlotMSDataTab.cc
  GuiTabs/PlotMSDisplayTab.cc
  GuiTabs/PlotMSIterateTab.cc
  GuiTabs/PlotMSExportTab.cc
  GuiTabs/PlotMSFlaggingTab.cc
  GuiTabs/PlotMSOptionsTab.cc
  GuiTabs/PlotMSPlotTab.cc
  GuiTabs/PlotMSTab.cc
  GuiTabs/PlotMSToolsTab.cc
  GuiTabs/PlotMSTransformationsTab.cc
  PlotMS/PlotMSAveraging.cc
  PlotMS/PlotMS.cc
  PlotMS/PlotMSCalibration.cc
  PlotMS/PlotMSConstants.cc
  PlotMS/PlotMSDBusApp.cc
  PlotMS/PlotMSExportParam.cc
  PlotMS/PlotMSFlagging.cc
  PlotMS/PlotMSIterParam.cc
  PlotMS/PlotMSLabelFormat.cc
  PlotMS/PlotMSPageHeaderParam.cc
  PlotMS/PlotMSParameters.cc
  PlotMS/PlotMSRegions.cc
  PlotMS/PlotMSSelection.cc
  PlotMS/PlotMSTransformations.cc
  PlotMS/PlotMSWatchedParameters.cc
  Plots/PlotMSPage.cc
  Plots/PlotMSPages.cc
  Plots/PlotMSPlot.cc
  Plots/PlotMSPlotManager.cc
  Plots/PlotMSPlotParameterGroups.cc
  Plots/PlotMSPlotParameters.cc
  )

INSTALL (FILES
        Actions/PlotMSAction.h
	DESTINATION include/casacode/plotms/Actions )

INSTALL (FILES
        Threads/Gui/PlotMSCacheThread.qo.h
        Threads/Gui/PlotMSDrawThread.qo.h
        Threads/Gui/PlotMSExportThread.qo.h
        Threads/Gui/PlotMSThread.qo.h
	DESTINATION include/casacode/plotms/Threads/Gui )

INSTALL (FILES
        Data/PlotMSIndexer.h
        Data/PlotMSVBAverager.h
        Data/PageHeaderCache.h
	DESTINATION include/casacode/plotms/Data )

INSTALL (FILES
        Gui/PlotMSAnnotator.h
        Gui/PlotMSAveragingWidget.qo.h
        Gui/PlotMSAxisWidget.qo.h
        Gui/PlotMSCalibrationWidget.qo.h
        Gui/PlotMSFlaggingWidget.qo.h
        Gui/PlotMSLoggerWidget.qo.h
        Gui/PlotMSPlotter.qo.h
        Gui/PlotMSSelectionWidget.qo.h
        Gui/PlotMSSummaryDialog.qo.h
        Gui/PlotMSTransformationsWidget.qo.h
        Gui/PlotRangeWidget.qo.h
        Gui/PlotMSPageHeaderDataModel.qo.h
	DESTINATION include/casacode/plotms/Gui )

INSTALL (FILES
        GuiTabs/PlotMSAnnotatorTab.qo.h
        GuiTabs/PlotMSAxesTab.qo.h
        GuiTabs/PlotMSCacheTab.qo.h
        GuiTabs/PlotMSCalibrationTab.qo.h
        GuiTabs/PlotMSCanvasTab.qo.h
        GuiTabs/PlotMSDataCollapsible.qo.h
        GuiTabs/PlotMSDataSummaryTab.qo.h
        GuiTabs/PlotMSDataSymbolWidget.qo.h
        GuiTabs/PlotMSDataTab.qo.h
        GuiTabs/PlotMSDisplayTab.qo.h
        GuiTabs/PlotMSExportTab.qo.h
        GuiTabs/PlotMSFlaggingTab.qo.h
        GuiTabs/PlotMSIterateTab.qo.h
        GuiTabs/PlotMSOptionsTab.qo.h
        GuiTabs/PlotMSPlotTab.qo.h
        GuiTabs/PlotMSTab.qo.h
        GuiTabs/PlotMSToolsTab.qo.h
        GuiTabs/PlotMSTransformationsTab.qo.h
	DESTINATION include/casacode/plotms/GuiTabs )

INSTALL (FILES
        PlotMS/PlotMS.h
        PlotMS/PlotMSAveraging.h
        PlotMS/PlotMSCalibration.h
        PlotMS/PlotMSConstants.h
        PlotMS/PlotMSDBusApp.h
        PlotMS/PlotMSExportParam.h
        PlotMS/PlotMSFlagging.h
        PlotMS/PlotMSIterParam.h
        PlotMS/PlotMSLabelFormat.h
        PlotMS/PlotMSPageHeaderParam.h
        PlotMS/PlotMSParameters.h
        PlotMS/PlotMSRegions.h
        PlotMS/PlotMSSelection.h
        PlotMS/PlotEngine.h
        PlotMS/PlotMSTransformations.h
        PlotMS/PlotMSWatchedParameters.h
	DESTINATION include/casacode/plotms/PlotMS )

INSTALL (FILES
        Plots/PlotMSPage.h
        Plots/PlotMSPages.h
        Plots/PlotMSPlot.h
        Plots/PlotMSPlotManager.h
        Plots/PlotMSPlotParameterGroups.h
        Plots/PlotMSPlotParameters.h
	DESTINATION include/casacode/plotms/Plots )


casa_add_executable( plotms casaplotms apps/casaplotms/casaplotms.cc )

#TESTING
# The app tests cannot be run on build server -- no display, 
# so plotms tested in python with test_plotms.py (with xvfb).
# Kept one demo dExportPlot to demonstrate running app from C++
# and one test program for just loading cache,
# useful for debugging GTs.
casa_add_demo( plotms test/dExportPlot.cc )
casa_add_demo( plotms test/dCacheData.cc )

casa_add_google_test( MODULES plotms SOURCES test/PlotMSCacheMetadata_GT.cc )
casa_add_google_test( MODULES plotms SOURCES test/PlotMSCacheVis_GT.cc )
casa_add_google_test( MODULES plotms SOURCES test/PlotMSCacheGeom_GT.cc )
casa_add_google_test( MODULES plotms SOURCES test/PlotMSCacheInvalid_GT.cc )
casa_add_google_test( MODULES plotms SOURCES test/PlotMSCacheLtdVis_GT.cc )
casa_add_google_test( MODULES plotms SOURCES test/PlotMSCacheOverlays_GT.cc )
casa_add_google_test( MODULES plotms SOURCES test/PlotMSCacheCallib_GT.cc )
casa_add_google_test( MODULES plotms SOURCES test/PlotMSCacheSelection_GT.cc )
casa_add_google_test( MODULES plotms SOURCES test/PlotMSCacheTransform_GT.cc )
casa_add_google_test( MODULES plotms SOURCES test/PlotMSCacheAveraging_GT.cc )
