// ASDM2MSFiller.h: implementation of a casacore::MeasurementSet's filler
// for Francois Viallefond & Frederic Badia ALMA Simulator
//
//  Copyright (C) 2001
//  OBSERVATOIRE DE PARIS - DEMIRM
//  Avenue Denfert Rochereau - 75014 - PARIS
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//
//
//////////////////////////////////////////////////////////////////////
#if !defined(ALMA_ASDM2MSFILLER_H)
#define ALMA_ASDM2MSFILLER_H
//# Includes

#include <casa/aips.h>
#include <casa/Utilities/Assert.h>
#include <tables/Tables.h>
#include <ms/MeasurementSets/MeasurementSet.h>
#include <ms/MeasurementSets/MSAntennaColumns.h>
#include <ms/MeasurementSets/MSDataDescColumns.h>
#include <ms/MeasurementSets/MSFeedColumns.h>
#include <ms/MeasurementSets/MSFieldColumns.h>
#include <ms/MeasurementSets/MSFlagCmdColumns.h>
#include <ms/MeasurementSets/MSHistoryColumns.h>
#include <ms/MeasurementSets/MSMainColumns.h>

#include <ms/MeasurementSets/MSObsColumns.h>
#include <ms/MeasurementSets/MSPointingColumns.h>
#include <ms/MeasurementSets/MSPolColumns.h>
#include <ms/MeasurementSets/MSProcessorColumns.h>
#include <ms/MeasurementSets/MSSourceColumns.h>
#include <ms/MeasurementSets/MSStateColumns.h>
#include <ms/MeasurementSets/MSSpWindowColumns.h>
#include <ms/MeasurementSets/MSSysCalColumns.h>
#include <ms/MeasurementSets/MSWeatherColumns.h>

#include <tables/DataMan/StandardStMan.h>
#include <tables/DataMan/TiledShapeStMan.h>
#include <tables/Tables/SetupNewTab.h>
#include <tables/Tables/TableDesc.h>
#include <tables/Tables/TableRecord.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Cube.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayUtil.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Containers/Block.h>
#include <casa/Containers/OrderedMap.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MBaseline.h>
#include <measures/Measures/Muvw.h>
#include <measures/Measures/MeasTable.h>
#include <measures/Measures/Stokes.h>
#include <measures/Measures/MeasConvert.h>
#include <measures/TableMeasures/TableMeasValueDesc.h>
#include <measures/TableMeasures/TableMeasOffsetDesc.h>
#include <measures/TableMeasures/TableMeasRefDesc.h>
#include <measures/TableMeasures/TableMeasDesc.h>
#include <measures/Measures/MeasConvert.h>
#include <casa/BasicSL/Constants.h>
#include <casa/OS/File.h>
#include <casa/OS/Path.h>
#include <complex>

#include <vector>


namespace casacore{
    class MPosition;
    class MeasFrame;
    class MeasurementSet;
    class MSMainColumns;
}

//# Forward Declarations

namespace casa {
    class TimeRange;
}

//
// A structure to define a range of rows in the Pointing table where the attribute overTheTop is defined and with which value.
//
struct s_overTheTop {
  unsigned int  start;   // The index of the first row of the range. 
  unsigned int  len;     // The number of consecutive rows in the range.
  bool value;   // The value of overTheTop in that range.
};

// Class ddMgr is a utility to help for the management
// of DataDescription, SpectralWindow and Polarization ids.
// Here we provide enough space to store 100 values for 
// each quantity; this is very likeky far beyond the actual
// needs.
class ddMgr {
 private:
  int     numCorr[100];
  int     numChan[100];
  struct  {
    int polId;
    int swId;
  } dd[100];
  
 public:

  ddMgr();

  int setNumCorr(int i, int numChan);
  int setNumChan(int i, int numCorr);

  int getNumCorr(int i);
  int getNumChan(int i);

  int setDD(int i, int polId, int swId); 

  int getPolId(int i);
  int getSwId(int i);
};


 
// Class ASDM2MSFiller
class ASDM2MSFiller {
 private:
  double         itsCreationTime;
  const std::string   itsName;
  int            itsNumAntenna;
  int            itsNumChan;
  int            itsNumCorr;
  casacore::MeasurementSet *itsMS;
  casacore::MSMainColumns  *itsMSCol;
  casacore::String     itsMSPath;
  casacore::Bool     itsWithRadioMeters;     /* Are we building an ALMA casacore::MS ?*/
  casacore::Bool     itsFirstScan;
  casacore::uInt     itsMSMainRow;
  /*casacore::TiledDataStManAccessor itsImWgtAcc;*/
  casacore::Block<casacore::IPosition> itsDataShapes;

  int itsScanNumber;
  int itsNCat;
    
  ddMgr    itsDDMgr;

  int itsCalDeviceNumberOfRows;
  casacore::Table itsMSCalDeviceTable;

  int createMS(const std::string& msName, 
               bool complexData, 
               bool withCompression, 
               const std::string& telName, 
               int maxNumCorr,
               int maxNumChan,
               bool withCorrectedData=false,
	       bool useAsdmStMan4DATA=false);

  const char** getPolCombinations(int numCorr);
    
  static std::map<std::string, casacore::MDirection::Types> string2MDirection;
  static std::map<std::string, casacore::MDirection::Types> string2MDirectionInit();
   
 public:  
  ASDM2MSFiller (const std::string&	name_,
		 double		creation_time_,
		 bool		withRadioMeters,
		 bool		complexData,
		 bool		withCompression,
                 const std::string&  telName, 
                 int            intintmaxNumCorr,
                 int            maxNumChan,
		 bool		withCorrectedData=false,
		 bool           useAsdmStMan4DATA=false);
  
  // Destructor
  ~ASDM2MSFiller();

  const casacore::MeasurementSet* ms();

  int addAntenna(const std::string&	 name_,
		 const std::string&	 station_,
		 double		 lx_,
		 double		 ly_,
		 double		 lz_,
		 double		 offset_x_,
		 double		 offset_y_,
		 double		 offset_z_,
		 float		 dish_diam_);


  void addData (bool                      complexData,
		std::vector<double>            &time_,
		std::vector<int>               &antennaId1_,
		std::vector<int>               &antennaId2_,
		std::vector<int>               &feedId1_,
		std::vector<int>               &feedId2_,
		std::vector<int>               &dataDescId_,
		int                       processorId_,
		int                       fieldId_,
		std::vector<double>            &interval_,
		std::vector<double>            &exposure_,
		std::vector<double>            &timeCentroid_,
		int                       scanNumber_,
		int                       arrayId_,
		int                       observationId_,
		std::vector<int>               &stateId_,
		std::vector<std::pair<int, int> >   &nChanNPol_,
		std::vector<double>            &uvw_,
		std::vector<double>            &weight_,
		std::vector<double>            &sigma_);

  void addData (bool                      complexData,
		std::vector<double>            &time_,
		std::vector<int>               &antennaId1_,
		std::vector<int>               &antennaId2_,
		std::vector<int>               &feedId1_,
		std::vector<int>               &feedId2_,
		std::vector<int>               &dataDescId_,
		int                       processorId_,
		std::vector<int>               &fieldId_,
		std::vector<double>            &interval_,
		std::vector<double>            &exposure_,
		std::vector<double>            &timeCentroid_,
		int                       scanNumber_,
		int                       arrayId_,
		int                       observationId_,
		std::vector<int>               &stateId_,
		std::vector<double>            &uvw_,
		std::vector<std::vector<unsigned int> >      &dataShape_,
		std::vector<float *>           &uncorrectedData_,
		std::vector<float *>           &correctedData_,
		std::vector<unsigned int>      &flag_);

  void addData (bool                      complexData,
		std::vector<double>            &time_,
		std::vector<int>               &antennaId1_,
		std::vector<int>               &antennaId2_,
		std::vector<int>               &feedId1_,
		std::vector<int>               &feedId2_,
		std::vector<int>               &dataDescId_,
		int                       processorId_,
		std::vector<int>               &fieldId_,
		std::vector<double>            &interval_,
		std::vector<double>            &exposure_,
		std::vector<double>            &timeCentroid_,
		int                       scanNumber_,
		int                       arrayId_,
		int                       observationId_,
		std::vector<int>               &stateId_,
		std::vector<double>            &uvw_,
		std::vector<std::vector<unsigned int> >      &dataShape_,
		std::vector<float *>           &data_,
		std::vector<unsigned int>      &flag_,
		std::vector<double>            &weight_,
		std::vector<double>            &sigma_);
  
  int  addDataDescription(int spectral_window_id_,
			  int polarizarion_id_);

  int  addUniqueDataDescription(int spectral_window_id_,
				int polarizarion_id_);

  int  exists(char *path);
  casacore::String msPath();


  void addFeed(int      antenna_id_,
	       int      feed_id_,
	       int      spectral_window_id_,
	       double   time_,
	       double   interval_,
	       int      num_receptors_,
	       int      beam_id_,
	       std::vector<double> &   beam_offset_,
	       std::vector<std::string> & pol_type_,
	       std::vector<std::complex<float> > & polarization_response_,
	       std::vector<double>&   position_,  // Must be a 3 elements vector !!!
	       std::vector<double>&   receptor_angle_);
  
  void addField( const std::string&			name_,
		 const std::string&			code_,
		 double				time_,
		 unsigned int                   num_poly_,
		 std::vector<std::vector<double> >&	delay_dir_,
		 std::vector<std::vector<double> >&	phase_dir_,
		 std::vector<std::vector<double> >&	reference_dir_,
		 const std::string&			direction_code_,
		 int				source_id_);

  void updateEphemerisIdInField(std::vector<std::pair<int, int> >& idxEphemerisId_v);

  void addFlagCmd(double	time_,
		  double	interval_,
		  const std::string& type_,
		  const std::string& reason_,
		  int		level_,
		  int		severity_,
		  int		applied_,
		  std::string&	command_);

  void addHistory( double		time_,
		   int			observation_id_,
		   const std::string&	message_,
		   const std::string&	priority_,
		   const std::string&	origin_,
		   int			object_id_,
		   const std::string&	application_,
		   const std::string&	cli_command_,
		   const std::string&	app_parms_ );

  void addObservation(const std::string&		telescopeName_,
		      double			startTime_,
		      double			endTime_,
		      const std::string&		observer_,
		      const std::vector<std::string>&	log_,
		      const std::string&		schedule_type_,
		      const std::vector<std::string>&	schedule_,
		      const std::string&		project_,
		      double			release_date_);

  void addPointingSlice(unsigned int                  n_row_,
			std::vector<int>&                  antenna_id_,
			std::vector<double>&               time_,
			std::vector<double>&               interval_,
			std::vector<double>&               direction_,
			std::vector<double>&               target_,
			std::vector<double>&               pointing_offset_,
			std::vector<double>&               encoder_,
			std::vector<bool>&                 tracking_,
			bool                          overTheTopExists4All_,
			std::vector<bool>&                 v_overTheTop_,
			std::vector<s_overTheTop>&         v_s_overTheTop_);

  int  addPolarization(int num_corr_,
		       std::vector<int>& corr_type_,
		       std::vector<int>& corr_product_);

  int addUniquePolarization(int num_corr_,
			    //			    const std::vector<casacore::Stokes::StokesTypes>& corr_type_,
			    const std::vector<int>& corr_type_,
			    const std::vector<int>& corr_product_);

  void addProcessor(std::string& type_,
		    std::string& sub_type_,
		    int  type_id_,
		    int  mode_id_);

  void addSource(int             source_id_,
		 double          time_,
		 double          interval_,
		 int             spectral_window_id_,
		 int             num_lines_,
		 std::string&         name_,
		 int             calibration_group_,
		 std::string&         code_,
		 std::vector<double>& direction_,
		 std::string&         direction_code_,
		 std::vector<double>& position_,
		 std::vector<double>& proper_motion_,
		 std::vector<std::string>& transition_,
		 std::vector<double>& rest_frequency_,
		 std::vector<double>& sysvel_);
		 
  int  addSpectralWindow(int			num_chan_,
			 const std::string&          name_,
			 double			ref_frequency_,
			 const std::vector<double>&	chan_freq_,
			 const std::vector<double>&	chan_width_,
			 int			meas_freq_ref_,
			 const std::vector<double>&	effective_bw_,
			 const std::vector<double>&	resolution_,
			 double			total_bandwidth_,
			 int			net_sideband_,
			 int			bbc_no_,
			 int			if_conv_chain_,
			 int			freq_group_,
			 const std::string&		freq_group_name_,
			 int			num_assoc_,
			 const std::vector<int>&	assoc_sp_id_,
			 const std::vector<std::string>&	assoc_nature_);

  int  addUniqueState(bool sig_,
		      bool ref_,
		      double cal_,
		      double load_,
		      unsigned int sub_scan_,
		      std::string& obs_mode_,
		      bool flag_row_);
  
  
  void addState(bool    sig_,
		bool    ref_,
		double  cal_,
		double  load_,
		int     sub_scan_,
		std::string& obs_mode_);
  
  void addSysCal(int    antenna_id,
		 int    feed_id,
		 int    spectral_window_id,
		 double time_,
		 double interval_,
		 int    numReceptor_,
		 int    numChan_,
		 std::pair<bool, std::vector<float> >& tcal_spectrum_pair,
		 std::pair<bool, bool>&           tcal_flag_pair,
		 std::pair<bool, std::vector<float> >& trx_spectrum_pair,
		 std::pair<bool, bool>&           trx_flag_pair,
		 std::pair<bool, std::vector<float> >& tsky_spectrum_pair,
		 std::pair<bool, bool>&           tsky_flag_pair,
		 std::pair<bool, std::vector<float> >& tsys_spectrum_pair,
		 std::pair<bool, bool>&           tsys_flag_pair,
		 std::pair<bool, std::vector<float> >& tant_spectrum_pair,
		 std::pair<bool, bool>&           tant_flag_pair,
		 std::pair<bool, std::vector<float> >& tant_tsys_spectrum_pair,
		 std::pair<bool, bool>&           tant_tsys_flag_pair);

 void addWeather(int				antenna_id_,
		  double			time_,
		  double			interval_,
		  const std::pair<bool, float>&	pressure_opt_,
		  const std::pair<bool, float>&	relHumidity_opt_,
		  const std::pair<bool, float>&	temperature_opt_,
		  const std::pair<bool, float>&	windDirection_opt_,
		  const std::pair<bool, float>&	windSpeed_opt_,
		  const std::pair<bool, float>&	dewPoint_opt_,
		  int				wx_station_id_,
		  std::vector<double>&		wx_station_position_);

  /**
   * Add one row in the casacore::MS CALDEVICE table.
   *
   * @param antennaId the index in the ANTENNA table of the antenna for which this row is defined.
   * @param feedId the index in the FEED table of the feeds for which this row is defined.
   * @param spectralWindowId the index in the SPECTRAL WINDOW table of the spectral window for which this row is defined.
   * @param time the midpoint time of measurement.
   * @param interval the interval of measurement.
   * @param numCalload the number of calibration loads.
   * @param calLoadNames a vector of strings.
   */
  void addCalDevice(int				antennaId,
		    int				feedId,
		    int				spectralWindowId,
		    double			time,
		    double			interval,
		    unsigned int		numCalLoad,
		    std::vector<std::string>		calloadNames,
		    unsigned int		numReceptor,
		    std::vector<std::vector<float> >&	calEff,
		    std::vector<std::vector<float> >&	noiseCal,
		    std::vector<double >&		temperatureLoad);

  /**
   * Adds one row in the casacore::MS SYSPOWER table.
   *
   *
   * @param antennaId the index in the ANTENNA table of the antenna for which this row is defined.
   * @param feedId the index in the FEED table of the feeds for which this row is defined.
   * @param spectralWindowId the index in the SPECTRAL WINDOW table of the spectral window for which this row is defined.
   * @param time the midpoint time of measurement.
   * @param interval the interval of measurement.
   * @param numReceptor a null value will be interpreted as "all the optional attributes" are absent, otherwise it will be considered as
   * as the number of useful values to read in the next three vectors. More precisely, for any of the parameters switchedPowerDifference, 
   * switchedPowerSum and requantizedGain, if its size is null then the parameter is considered as "absent" and ignored otherwise
   * its first numReceptor values will be copied into the casacore::MS CalDevice table corresponding field. If the size of one the parameters is not null
   * but smaller than numReceptor then the code will crash miserably.
   *
   * @param switchedPowerDifference a vector of float numbers containing the switched power differences in its numReceptor first elements on
   * the basis of one value per receptor. If the size of the vector is null then it'll be ignored.
   * @param switchedPowerSum a vector of float numbers containing the switched power sums in its numReceptor first elements on
   * the basis of one value per receptor. If the size of the vector is null then it'll be ignored. 
   * @param temperatureLoad a vector of float number containing the requantizer gains on the basis of one value per receptor. If the size
   * of the vector is null then it'll be ignored. 
   */
  void addSysPower(int			antennaId,
		   int			feedId,
		   int			spectralWindowId,
		   double		time,
		   double		interval,
		   unsigned int         numReceptor,
		   std::vector<float>&	switchedPowerDifference,
		   std::vector<float>&	switchedPowerSum,
		   std::vector<float>&	requantizerGain); 

  void addSysPowerSlice(unsigned int	nRow,
			std::vector<int>&    antennaId,
			std::vector<int>&	spectralWindowId,
			std::vector<int>&	feedId,
			std::vector<double>&	time,
			std::vector<double>&	interval,
			unsigned int    numReceptor,
			std::vector<float>&	switchedPowerDifference,
			std::vector<float>&	switchedPowerSum,
			std::vector<float>&	requantizerGain);

  void end();
};

#endif
