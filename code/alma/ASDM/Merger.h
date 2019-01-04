
/*
 * ALMA - Atacama Large Millimeter Array
 * (c) European Southern Observatory, 2002
 * (c) Associated Universities Inc., 2002
 * Copyright by ESO (in the framework of the ALMA collaboration),
 * Copyright by AUI (in the framework of the ALMA collaboration),
 * All rights reserved.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY, without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA 02111-1307  USA
 *
 * Warning!
 *  -------------------------------------------------------------------- 
 * | This is generated code!  Do not modify this file.                  |
 * | If you do, all changes will be lost when the file is re-generated. |
 *  --------------------------------------------------------------------
 *
 * File Merger.h
 */
#ifndef Merger_CLASS
#define Merger_CLASS
#include <map>

#include <alma/ASDM/ASDM.h>

#include <alma/ASDM/SBSummaryTable.h>
#include <alma/ASDM/SBSummaryRow.h>

#include <alma/ASDM/ConfigDescriptionTable.h>
#include <alma/ASDM/ConfigDescriptionRow.h>

#include <alma/ASDM/FieldTable.h>
#include <alma/ASDM/FieldRow.h>

#include <alma/ASDM/StateTable.h>
#include <alma/ASDM/StateRow.h>

#include <alma/ASDM/AntennaTable.h>
#include <alma/ASDM/AntennaRow.h>

#include <alma/ASDM/DataDescriptionTable.h>
#include <alma/ASDM/DataDescriptionRow.h>

#include <alma/ASDM/SwitchCycleTable.h>
#include <alma/ASDM/SwitchCycleRow.h>

#include <alma/ASDM/SourceTable.h>
#include <alma/ASDM/SourceRow.h>

#include <alma/ASDM/FeedTable.h>
#include <alma/ASDM/FeedRow.h>

#include <alma/ASDM/SpectralWindowTable.h>
#include <alma/ASDM/SpectralWindowRow.h>

#include <alma/ASDM/FreqOffsetTable.h>
#include <alma/ASDM/FreqOffsetRow.h>

#include <alma/ASDM/PolarizationTable.h>
#include <alma/ASDM/PolarizationRow.h>

#include <alma/ASDM/ReceiverTable.h>
#include <alma/ASDM/ReceiverRow.h>

#include <alma/ASDM/DopplerTable.h>
#include <alma/ASDM/DopplerRow.h>

#include <alma/ASDM/ProcessorTable.h>
#include <alma/ASDM/ProcessorRow.h>

#include <alma/ASDM/CorrelatorModeTable.h>
#include <alma/ASDM/CorrelatorModeRow.h>

#include <alma/ASDM/CalDeviceTable.h>
#include <alma/ASDM/CalDeviceRow.h>

#include <alma/ASDM/FlagCmdTable.h>
#include <alma/ASDM/FlagCmdRow.h>

#include <alma/ASDM/FocusTable.h>
#include <alma/ASDM/FocusRow.h>

#include <alma/ASDM/HistoryTable.h>
#include <alma/ASDM/HistoryRow.h>

#include <alma/ASDM/ObservationTable.h>
#include <alma/ASDM/ObservationRow.h>

#include <alma/ASDM/PointingTable.h>
#include <alma/ASDM/PointingRow.h>

#include <alma/ASDM/SeeingTable.h>
#include <alma/ASDM/SeeingRow.h>

#include <alma/ASDM/SysCalTable.h>
#include <alma/ASDM/SysCalRow.h>

#include <alma/ASDM/TotalPowerTable.h>
#include <alma/ASDM/TotalPowerRow.h>

#include <alma/ASDM/WeatherTable.h>
#include <alma/ASDM/WeatherRow.h>

#include <alma/ASDM/WVMCalTable.h>
#include <alma/ASDM/WVMCalRow.h>

#include <alma/ASDM/EphemerisTable.h>
#include <alma/ASDM/EphemerisRow.h>

#include <alma/ASDM/ExecBlockTable.h>
#include <alma/ASDM/ExecBlockRow.h>

#include <alma/ASDM/ScanTable.h>
#include <alma/ASDM/ScanRow.h>

#include <alma/ASDM/SubscanTable.h>
#include <alma/ASDM/SubscanRow.h>

#include <alma/ASDM/MainTable.h>
#include <alma/ASDM/MainRow.h>

#include <alma/ASDM/FocusModelTable.h>
#include <alma/ASDM/FocusModelRow.h>

#include <alma/ASDM/GainTrackingTable.h>
#include <alma/ASDM/GainTrackingRow.h>

#include <alma/ASDM/PointingModelTable.h>
#include <alma/ASDM/PointingModelRow.h>

#include <alma/ASDM/CalAmpliTable.h>
#include <alma/ASDM/CalAmpliRow.h>

#include <alma/ASDM/CalDataTable.h>
#include <alma/ASDM/CalDataRow.h>

#include <alma/ASDM/CalReductionTable.h>
#include <alma/ASDM/CalReductionRow.h>

#include <alma/ASDM/CalPhaseTable.h>
#include <alma/ASDM/CalPhaseRow.h>

#include <alma/ASDM/CalSeeingTable.h>
#include <alma/ASDM/CalSeeingRow.h>

#include <alma/ASDM/CalPositionTable.h>
#include <alma/ASDM/CalPositionRow.h>

#include <alma/ASDM/CalPointingTable.h>
#include <alma/ASDM/CalPointingRow.h>

#include <alma/ASDM/CalPointingModelTable.h>
#include <alma/ASDM/CalPointingModelRow.h>

#include <alma/ASDM/CalHolographyTable.h>
#include <alma/ASDM/CalHolographyRow.h>

#include <alma/ASDM/CalAtmosphereTable.h>
#include <alma/ASDM/CalAtmosphereRow.h>

#include <alma/ASDM/CalCurveTable.h>
#include <alma/ASDM/CalCurveRow.h>

#include <alma/ASDM/StationTable.h>
#include <alma/ASDM/StationRow.h>

#include <alma/ASDM/AlmaRadiometerTable.h>
#include <alma/ASDM/AlmaRadiometerRow.h>

#include <alma/ASDM/SquareLawDetectorTable.h>
#include <alma/ASDM/SquareLawDetectorRow.h>

#include <alma/ASDM/CalFocusTable.h>
#include <alma/ASDM/CalFocusRow.h>

#include <alma/ASDM/CalDelayTable.h>
#include <alma/ASDM/CalDelayRow.h>

#include <alma/ASDM/HolographyTable.h>
#include <alma/ASDM/HolographyRow.h>

#include <alma/ASDM/CalBandpassTable.h>
#include <alma/ASDM/CalBandpassRow.h>

#include <alma/ASDM/CalFluxTable.h>
#include <alma/ASDM/CalFluxRow.h>

#include <alma/ASDM/CalFocusModelTable.h>
#include <alma/ASDM/CalFocusModelRow.h>

#include <alma/ASDM/CalGainTable.h>
#include <alma/ASDM/CalGainRow.h>

#include <alma/ASDM/CalPrimaryBeamTable.h>
#include <alma/ASDM/CalPrimaryBeamRow.h>

#include <alma/ASDM/CalWVRTable.h>
#include <alma/ASDM/CalWVRRow.h>

#include <alma/ASDM/AnnotationTable.h>
#include <alma/ASDM/AnnotationRow.h>

#include <alma/ASDM/DelayModelTable.h>
#include <alma/ASDM/DelayModelRow.h>

#include <alma/ASDM/FlagTable.h>
#include <alma/ASDM/FlagRow.h>

#include <alma/ASDM/ScaleTable.h>
#include <alma/ASDM/ScaleRow.h>

#include <alma/ASDM/SysPowerTable.h>
#include <alma/ASDM/SysPowerRow.h>

#include <alma/ASDM/CalAppPhaseTable.h>
#include <alma/ASDM/CalAppPhaseRow.h>

#include <alma/ASDM/DelayModelFixedParametersTable.h>
#include <alma/ASDM/DelayModelFixedParametersRow.h>

#include <alma/ASDM/DelayModelVariableParametersTable.h>
#include <alma/ASDM/DelayModelVariableParametersRow.h>

#include <alma/ASDM/CalAntennaSolutionsTable.h>
#include <alma/ASDM/CalAntennaSolutionsRow.h>

#include <alma/ASDM/PulsarTable.h>
#include <alma/ASDM/PulsarRow.h>



/*\file <alma/ASDM/Merger.h>
    \brief Generated from model's revision "-1", branch ""
*/

namespace asdm {
	class Merger {
		public :
			Merger();
			virtual ~Merger();
			
			void merge(ASDM* ds1, ASDM* ds2);
			
		private:
			ASDM* ds1;
			ASDM* ds2;
			std::map<std::string, Tag> tagTag;
			Tag getTag(const Tag& t, void (Merger::*mergeTables)());
			
			std::map<std::string, int> idId;
			int getId(const std::string& tableName, int id, void (Merger::*mergeTables)()); 
			

			bool hasMergedSBSummary;	

			bool hasMergedConfigDescription;	

			bool hasMergedField;	

			bool hasMergedState;	

			bool hasMergedAntenna;	

			bool hasMergedDataDescription;	

			bool hasMergedSwitchCycle;	

			bool hasMergedSource;	

			bool hasMergedFeed;	

			bool hasMergedSpectralWindow;	

			bool hasMergedFreqOffset;	

			bool hasMergedPolarization;	

			bool hasMergedReceiver;	

			bool hasMergedDoppler;	

			bool hasMergedProcessor;	

			bool hasMergedCorrelatorMode;	

			bool hasMergedCalDevice;	

			bool hasMergedFlagCmd;	

			bool hasMergedFocus;	

			bool hasMergedHistory;	

			bool hasMergedObservation;	

			bool hasMergedPointing;	

			bool hasMergedSeeing;	

			bool hasMergedSysCal;	

			bool hasMergedTotalPower;	

			bool hasMergedWeather;	

			bool hasMergedWVMCal;	

			bool hasMergedEphemeris;	

			bool hasMergedExecBlock;	

			bool hasMergedScan;	

			bool hasMergedSubscan;	

			bool hasMergedMain;	

			bool hasMergedFocusModel;	

			bool hasMergedGainTracking;	

			bool hasMergedPointingModel;	

			bool hasMergedCalAmpli;	

			bool hasMergedCalData;	

			bool hasMergedCalReduction;	

			bool hasMergedCalPhase;	

			bool hasMergedCalSeeing;	

			bool hasMergedCalPosition;	

			bool hasMergedCalPointing;	

			bool hasMergedCalPointingModel;	

			bool hasMergedCalHolography;	

			bool hasMergedCalAtmosphere;	

			bool hasMergedCalCurve;	

			bool hasMergedStation;	

			bool hasMergedAlmaRadiometer;	

			bool hasMergedSquareLawDetector;	

			bool hasMergedCalFocus;	

			bool hasMergedCalDelay;	

			bool hasMergedHolography;	

			bool hasMergedCalBandpass;	

			bool hasMergedCalFlux;	

			bool hasMergedCalFocusModel;	

			bool hasMergedCalGain;	

			bool hasMergedCalPrimaryBeam;	

			bool hasMergedCalWVR;	

			bool hasMergedAnnotation;	

			bool hasMergedDelayModel;	

			bool hasMergedFlag;	

			bool hasMergedScale;	

			bool hasMergedSysPower;	

			bool hasMergedCalAppPhase;	

			bool hasMergedDelayModelFixedParameters;	

			bool hasMergedDelayModelVariableParameters;	

			bool hasMergedCalAntennaSolutions;	

			bool hasMergedPulsar;	
			


			void mergeSBSummary();
			void postMergeSBSummary();			

			void mergeConfigDescription();
			void postMergeConfigDescription();			

			void mergeField();
			void postMergeField();			

			void mergeState();
			void postMergeState();			

			void mergeAntenna();
			void postMergeAntenna();			

			void mergeDataDescription();
			void postMergeDataDescription();			

			void mergeSwitchCycle();
			void postMergeSwitchCycle();			

			void mergeSource();
			void postMergeSource();			

			void mergeFeed();
			void postMergeFeed();			

			void mergeSpectralWindow();
			void postMergeSpectralWindow();			

			void mergeFreqOffset();
			void postMergeFreqOffset();			

			void mergePolarization();
			void postMergePolarization();			

			void mergeReceiver();
			void postMergeReceiver();			

			void mergeDoppler();
			void postMergeDoppler();			

			void mergeProcessor();
			void postMergeProcessor();			

			void mergeCorrelatorMode();
			void postMergeCorrelatorMode();			

			void mergeCalDevice();
			void postMergeCalDevice();			

			void mergeFlagCmd();
			void postMergeFlagCmd();			

			void mergeFocus();
			void postMergeFocus();			

			void mergeHistory();
			void postMergeHistory();			

			void mergeObservation();
			void postMergeObservation();			

			void mergePointing();
			void postMergePointing();			

			void mergeSeeing();
			void postMergeSeeing();			

			void mergeSysCal();
			void postMergeSysCal();			

			void mergeTotalPower();
			void postMergeTotalPower();			

			void mergeWeather();
			void postMergeWeather();			

			void mergeWVMCal();
			void postMergeWVMCal();			

			void mergeEphemeris();
			void postMergeEphemeris();			

			void mergeExecBlock();
			void postMergeExecBlock();			

			void mergeScan();
			void postMergeScan();			

			void mergeSubscan();
			void postMergeSubscan();			

			void mergeMain();
			void postMergeMain();			

			void mergeFocusModel();
			void postMergeFocusModel();			

			void mergeGainTracking();
			void postMergeGainTracking();			

			void mergePointingModel();
			void postMergePointingModel();			

			void mergeCalAmpli();
			void postMergeCalAmpli();			

			void mergeCalData();
			void postMergeCalData();			

			void mergeCalReduction();
			void postMergeCalReduction();			

			void mergeCalPhase();
			void postMergeCalPhase();			

			void mergeCalSeeing();
			void postMergeCalSeeing();			

			void mergeCalPosition();
			void postMergeCalPosition();			

			void mergeCalPointing();
			void postMergeCalPointing();			

			void mergeCalPointingModel();
			void postMergeCalPointingModel();			

			void mergeCalHolography();
			void postMergeCalHolography();			

			void mergeCalAtmosphere();
			void postMergeCalAtmosphere();			

			void mergeCalCurve();
			void postMergeCalCurve();			

			void mergeStation();
			void postMergeStation();			

			void mergeAlmaRadiometer();
			void postMergeAlmaRadiometer();			

			void mergeSquareLawDetector();
			void postMergeSquareLawDetector();			

			void mergeCalFocus();
			void postMergeCalFocus();			

			void mergeCalDelay();
			void postMergeCalDelay();			

			void mergeHolography();
			void postMergeHolography();			

			void mergeCalBandpass();
			void postMergeCalBandpass();			

			void mergeCalFlux();
			void postMergeCalFlux();			

			void mergeCalFocusModel();
			void postMergeCalFocusModel();			

			void mergeCalGain();
			void postMergeCalGain();			

			void mergeCalPrimaryBeam();
			void postMergeCalPrimaryBeam();			

			void mergeCalWVR();
			void postMergeCalWVR();			

			void mergeAnnotation();
			void postMergeAnnotation();			

			void mergeDelayModel();
			void postMergeDelayModel();			

			void mergeFlag();
			void postMergeFlag();			

			void mergeScale();
			void postMergeScale();			

			void mergeSysPower();
			void postMergeSysPower();			

			void mergeCalAppPhase();
			void postMergeCalAppPhase();			

			void mergeDelayModelFixedParameters();
			void postMergeDelayModelFixedParameters();			

			void mergeDelayModelVariableParameters();
			void postMergeDelayModelVariableParameters();			

			void mergeCalAntennaSolutions();
			void postMergeCalAntennaSolutions();			

			void mergePulsar();
			void postMergePulsar();			



		void (Merger::*mergeSBSummaryPtr) () ;

		void (Merger::*mergeConfigDescriptionPtr) () ;

		void (Merger::*mergeFieldPtr) () ;

		void (Merger::*mergeStatePtr) () ;

		void (Merger::*mergeAntennaPtr) () ;

		void (Merger::*mergeDataDescriptionPtr) () ;

		void (Merger::*mergeSwitchCyclePtr) () ;

		void (Merger::*mergeSourcePtr) () ;

		void (Merger::*mergeFeedPtr) () ;

		void (Merger::*mergeSpectralWindowPtr) () ;

		void (Merger::*mergeFreqOffsetPtr) () ;

		void (Merger::*mergePolarizationPtr) () ;

		void (Merger::*mergeReceiverPtr) () ;

		void (Merger::*mergeDopplerPtr) () ;

		void (Merger::*mergeProcessorPtr) () ;

		void (Merger::*mergeCorrelatorModePtr) () ;

		void (Merger::*mergeCalDevicePtr) () ;

		void (Merger::*mergeFlagCmdPtr) () ;

		void (Merger::*mergeFocusPtr) () ;

		void (Merger::*mergeHistoryPtr) () ;

		void (Merger::*mergeObservationPtr) () ;

		void (Merger::*mergePointingPtr) () ;

		void (Merger::*mergeSeeingPtr) () ;

		void (Merger::*mergeSysCalPtr) () ;

		void (Merger::*mergeTotalPowerPtr) () ;

		void (Merger::*mergeWeatherPtr) () ;

		void (Merger::*mergeWVMCalPtr) () ;

		void (Merger::*mergeEphemerisPtr) () ;

		void (Merger::*mergeExecBlockPtr) () ;

		void (Merger::*mergeScanPtr) () ;

		void (Merger::*mergeSubscanPtr) () ;

		void (Merger::*mergeMainPtr) () ;

		void (Merger::*mergeFocusModelPtr) () ;

		void (Merger::*mergeGainTrackingPtr) () ;

		void (Merger::*mergePointingModelPtr) () ;

		void (Merger::*mergeCalAmpliPtr) () ;

		void (Merger::*mergeCalDataPtr) () ;

		void (Merger::*mergeCalReductionPtr) () ;

		void (Merger::*mergeCalPhasePtr) () ;

		void (Merger::*mergeCalSeeingPtr) () ;

		void (Merger::*mergeCalPositionPtr) () ;

		void (Merger::*mergeCalPointingPtr) () ;

		void (Merger::*mergeCalPointingModelPtr) () ;

		void (Merger::*mergeCalHolographyPtr) () ;

		void (Merger::*mergeCalAtmospherePtr) () ;

		void (Merger::*mergeCalCurvePtr) () ;

		void (Merger::*mergeStationPtr) () ;

		void (Merger::*mergeAlmaRadiometerPtr) () ;

		void (Merger::*mergeSquareLawDetectorPtr) () ;

		void (Merger::*mergeCalFocusPtr) () ;

		void (Merger::*mergeCalDelayPtr) () ;

		void (Merger::*mergeHolographyPtr) () ;

		void (Merger::*mergeCalBandpassPtr) () ;

		void (Merger::*mergeCalFluxPtr) () ;

		void (Merger::*mergeCalFocusModelPtr) () ;

		void (Merger::*mergeCalGainPtr) () ;

		void (Merger::*mergeCalPrimaryBeamPtr) () ;

		void (Merger::*mergeCalWVRPtr) () ;

		void (Merger::*mergeAnnotationPtr) () ;

		void (Merger::*mergeDelayModelPtr) () ;

		void (Merger::*mergeFlagPtr) () ;

		void (Merger::*mergeScalePtr) () ;

		void (Merger::*mergeSysPowerPtr) () ;

		void (Merger::*mergeCalAppPhasePtr) () ;

		void (Merger::*mergeDelayModelFixedParametersPtr) () ;

		void (Merger::*mergeDelayModelVariableParametersPtr) () ;

		void (Merger::*mergeCalAntennaSolutionsPtr) () ;

		void (Merger::*mergePulsarPtr) () ;

	};
} // End namespace asdm

#endif  /* Merger_CLASS */
