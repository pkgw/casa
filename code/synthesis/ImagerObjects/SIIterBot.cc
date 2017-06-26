//# SIIterBot.cc: This file contains the implementation of the SIIterBot class.
//#
//#  CASA - Common Astronomy Software Applications (http://casa.nrao.edu/)
//#  Copyright (C) Associated Universities, Inc. Washington DC, USA 2011, All rights reserved.
//#  Copyright (C) European Southern Observatory, 2011, All rights reserved.
//#
//#  This library is free software; you can redistribute it and/or
//#  modify it under the terms of the GNU Lesser General Public
//#  License as published by the Free software Foundation; either
//#  version 2.1 of the License, or (at your option) any later version.
//#
//#  This library is distributed in the hope that it will be useful,
//#  but WITHOUT ANY WARRANTY, without even the implied warranty of
//#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//#  Lesser General Public License for more details.
//#
//#  You should have received a copy of the GNU Lesser General Public
//#  License along with this library; if not, write to the Free Software
//#  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//#  MA 02111-1307  USA
//# $Id: $

#include <synthesis/ImagerObjects/SIIterBot.h>
#include <casadbus/session/DBusSession.h>
#include <casadbus/utilities/Conversion.h>

/* Include file for the lock guard */
#include <mutex>

/* Records Interface */
#include <casa/Containers/Record.h>
#include <math.h>						// For FLT_MAX

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

	////////////////////////////////////
	/// SIIterBot_callback implementation ///
	////////////////////////////////////
	void SIIterBot_callback::interactionRequired( bool val ) {
		std::lock_guard<std::recursive_mutex> guard(mutex); 
		if ( adaptor != 0 )
			adaptor->interactionRequired(val);
	}

	void SIIterBot_callback::addHandler( SIIterBot_adaptor *adapt ) {
		std::lock_guard<std::recursive_mutex> guard(mutex); 
		if ( adaptor == 0 )
			adaptor = adapt;
		
	}
	void SIIterBot_callback::removeHandler( SIIterBot_adaptor *adapt ) {
		std::lock_guard<std::recursive_mutex> guard(mutex); 
		if ( adaptor == adapt )
			adaptor = 0;
	}

  
	////////////////////////////////////
	/// SIIterBot_state implementation ///
	////////////////////////////////////
  
	// All SIIterBot_states must have 'type' and 'name' defined.
	SIIterBot_state::SIIterBot_state( SHARED_PTR<SIIterBot_callback> cb ) :
						itsDescription("no description is currently available..."),
						itsMinPsfFraction(0.05),
						itsMaxPsfFraction(0.8),
						itsMaxPsfSidelobe(0.0),
						itsPeakResidual(0.0),
						itsPrevPeakResidual(-1.0),
						itsInitPeakResidual(0.0),
						itsMinPeakResidual(1e+9),
						itsMinorCyclePeakResidual(0.0),
						itsPeakResidualNoMask(0.0),
						itsPrevPeakResidualNoMask(-1.0),
						itsMinPeakResidualNoMask(1e+9),
						itsMadRMS(0.0),
						itsMaskSum(-1.0),
						itsPrevMajorCycleCount(0),
						itsControllerCount(0),
						itsNiter(0),
						itsCycleNiter(0),
						itsInteractiveNiter(0),
						itsThreshold(0),
						itsCycleThreshold(0.0),
						itsInteractiveThreshold(0.0),
						itsIsCycleThresholdAuto(true),
						itsIsThresholdAuto(false),
						itsCycleFactor(1.0),
						itsLoopGain(0.1),
						itsStopFlag(false),
						itsPauseFlag(false),
						itsInteractiveMode(false),
						itsUpdatedModelFlag(false),
						itsIterDone(0),
						itsInteractiveIterDone(0),
						itsMaxCycleIterDone(0),
						itsMajorDone(0),
                                                itsStopCode(0),
						itsNSummaryFields(6),
						itsSummaryMinor(IPosition(2,6,0)),
						itsSummaryMajor(IPosition(1,0)),
						callback(cb)
	{
		LogIO os( LogOrigin("SIIterBot_state",__FUNCTION__,WHERE) );
	}
    
	SIIterBot_state::~SIIterBot_state() {
	  //		fprintf( stderr, ">>>>>>\t\tSIIterBot_state::~SIIterBot_state(0x%p)\n", this );
	  //		fflush( stderr );
	}

	bool SIIterBot_state::interactiveInputRequired( ) {
		std::lock_guard<std::recursive_mutex> guard(recordMutex); 
		return ( itsInteractiveMode &&
				 ( itsMaxCycleIterDone+itsInteractiveIterDone>=itsInteractiveNiter ||
				   itsPeakResidual <= itsInteractiveThreshold ||
				   itsPauseFlag ) );
	}

	void SIIterBot_state::waitForInteractiveInput() {
		/* Check that we have at least one controller */
		if (getNumberOfControllers() == 0) {
			/* Spawn a Viewer set up for interactive */
		}
		cout << "UU : setup interaction" << endl;
		std::unique_lock<std::mutex> lock(interactionMutex);
		if(!interactionPending) {
			interactionPending = true;
			pushDetails( );
			if ( callback )
				callback->interactionRequired(interactionPending);
		}
		cout << "UU : about to wait" << endl;
		/* Wait on Condition variable */
		while (interactionPending) {
			interactionCond.wait(lock);
		}

		cout << "UU : returned from wait" << endl;
		if (updateNeeded) {
			updateNeeded = false;
			if ( callback )
				callback->interactionRequired(false);
			pushDetails();
		}
	}

	int SIIterBot_state::cleanComplete(Bool lastcyclecheck){
		std::lock_guard<std::recursive_mutex> guard(recordMutex);    

		LogIO os( LogOrigin("SIIterBot_state",__FUNCTION__,WHERE) );

		//		printOut("FromcleanComplete ", false);

		int stopCode=0;

		Float usePeakRes;
		if(lastcyclecheck==True){usePeakRes = itsMinorCyclePeakResidual; }
		else{usePeakRes = itsPeakResidual; }
                
		//		cout << "itsMajorDone="<<itsMajorDone<<" itsIterDone="<<itsIterDone<< " itsInitPeakResidual="<<itsInitPeakResidual<<" itsPeakResidual="<<itsPeakResidual <<" itsPrevPeakResidual : " <<  itsPrevPeakResidual << " itsStopFlag="<<itsStopFlag<<endl;

		if( itsPeakResidual>0 && itsPrevPeakResidual>0 && 
		    fabs(itsPeakResidual - itsPrevPeakResidual)/fabs(itsPrevPeakResidual) > 2.0 )
		  {
		    os << "[WARN] Peak residual (within the mask) increased from " << itsPrevPeakResidual << " to " << itsPeakResidual << LogIO::POST;
		  }

		/// This may interfere with some other criterion... check.
		if ( itsMajorDone==0 && itsIterDone==0 ) { stopCode=0; }
		else if ( itsIterDone >= itsNiter || 
		     itsPeakResidual <= itsThreshold ||
		     itsStopFlag )
		  {
		    //		    os << "Reached global stopping criteria : ";

		    if( itsIterDone >= itsNiter ) { stopCode=1; }
		    //os << "Numer of iterations. "; // (" << itsIterDone << ") >= limit (" << itsNiter << ")" ;
		    if( usePeakRes <= itsThreshold ) {stopCode=2; }
		    //os << "Peak residual (" << itsPeakResidual << ") <= threshold(" << itsThreshold << ")";
		    if( itsStopFlag ) {stopCode=3;}
		      //os << "Forced stop. ";
		    //		    os << LogIO::POST;

		    //return true;
		  }
		else // not converged yet... but....if nothing has changed in this round... also stop
		  {
		    
		    if (itsMaskSum==0.0)
		      {
			//cout << "(7) Mask is all zero.Stopping" << endl;
			stopCode = 7;
		      }
		    // Nothing has changed across the last set of minor cycle iterations and major cycle.
		    else if( itsIterDone>0 && 
			     //itsMaskSum>0.0 &&
			     			     (itsMajorDone>itsPrevMajorCycleCount) && 
			     fabs(itsPrevPeakResidual - itsPeakResidual)<1e-10) 
		      {stopCode = 4;}
		    
                    // another non-convergent condition: diverging (relative increase is more than 5 times across one major cycle)
                    else if ( itsIterDone > 0 && 
			      fabs(itsPeakResidualNoMask-itsPrevPeakResidualNoMask)/fabs(itsPrevPeakResidualNoMask)  > 3.0) 
                      {
			//cout << "(5) Peak res (no mask) : " << itsPeakResidualNoMask 
			//     << "  Dev from prev peak res " << itsPrevPeakResidualNoMask << endl; 
			stopCode = 5;}

		    // divergence check, 5 times increase from the minimum peak residual so far (across all previous major cycles).
		    else if ( itsIterDone > 0 && 
			      (fabs(itsPeakResidualNoMask)-itsMinPeakResidualNoMask)/itsMinPeakResidualNoMask  > 3.0 )
                      {
			//cout << "(6) Peak res (no mask): " << itsPeakResidualNoMask 
			//    <<  "    Dev from min peak res " << itsMinPeakResidualNoMask << endl; 
			stopCode = 6;
		      }
		    
		  }

		/*
		if( lastcyclecheck==False)
		  {
		    cout << "*****" << endl;
		    cout << "Peak residual : " << itsPeakResidual << "  No Mask : " << itsPeakResidualNoMask << endl;
		    cout << "Prev Peak residual : " << itsPrevPeakResidual << "  No Mask : " << itsPrevPeakResidualNoMask << endl;
		    cout << "Min Peak residual : " << itsMinPeakResidual << "  No Mask : " << itsMinPeakResidualNoMask << endl;
		  }
		*/

		//		os << "Peak residual : " << itsPeakResidual << " and " << itsIterDone << " iterations."<< LogIO::POST;
		//cout << "cleancomp : stopcode : " << stopCode << endl;

		//cout << "peak res : " << itsPeakResidual << "   itsminPR : " << itsMinPeakResidual << endl;

		if( lastcyclecheck==False)
		  {
		    if( fabs(itsPeakResidual) < itsMinPeakResidual ) 
		      {itsMinPeakResidual = fabs(itsPeakResidual);}
		    
		    itsPrevPeakResidual = itsPeakResidual;


		    if( fabs(itsPeakResidualNoMask) < itsMinPeakResidualNoMask ) 
		      {itsMinPeakResidualNoMask = fabs(itsPeakResidualNoMask);}
		    
		    itsPrevPeakResidualNoMask = itsPeakResidualNoMask;

		    itsPrevMajorCycleCount = itsMajorDone;

		  }
		
		itsStopCode=stopCode;
		return stopCode;
	}


	Record SIIterBot_state::getMinorCycleControls(){
		LogIO os( LogOrigin("SIIterBot_state",__FUNCTION__,WHERE) );
		std::lock_guard<std::recursive_mutex> guard(recordMutex);    

		/* If autocalc, compute cyclethresh from peak res, cyclefactor and psf sidelobe 
		   Otherwise, the user has explicitly set it (interactively) for this minor cycle */
		if( itsIsCycleThresholdAuto == true ) { 
		  updateCycleThreshold(); 
		  //cout << "Updating cyc thresh" << endl; 
		}
		//		else { 
		//		  cout << "NOT updating cyc thresh" << endl; 
		//		}
		itsIsCycleThresholdAuto = true; /* Reset this, for the next round */

		/* Now that we have set the threshold, zero the peak residual 
		   so it can be found again after the minor cycles */
		//	itsInitPeakResidual = itsPeakResidual;
		//	itsPeakResidual = 0;


		/* This returns a record suitable for initializing the minor cycle
		   controls. */
		Record returnRecord;

		/* The minor cycle will stop based on the cycle parameters. */
		Int maxCycleIterations = itsCycleNiter;
		Float cycleThreshold     = itsCycleThreshold;
		maxCycleIterations = min(maxCycleIterations, itsNiter - itsIterDone);
		cycleThreshold = max(cycleThreshold, itsThreshold);
		/*
		if (itsInteractiveMode) {
			maxCycleIterations = min(maxCycleIterations, itsInteractiveNiter);
			cycleThreshold = max(cycleThreshold, itsInteractiveThreshold);
		}
		*/
		returnRecord.define( RecordFieldId("cycleniter"),  maxCycleIterations);
		returnRecord.define( RecordFieldId("cyclethreshold"), cycleThreshold);
		returnRecord.define( RecordFieldId("loopgain"), itsLoopGain);

		return returnRecord;
	}

	void SIIterBot_state::mergeCycleInitializationRecord(Record& initRecord){
		std::lock_guard<std::recursive_mutex> guard(recordMutex);  
    
		itsPeakResidual = max( itsPeakResidual,
							   initRecord.asFloat(RecordFieldId("peakresidual")) );
		itsMaxPsfSidelobe =  max( itsMaxPsfSidelobe, initRecord.asFloat(RecordFieldId("maxpsfsidelobe")) );

		itsPeakResidualNoMask = max( itsPeakResidualNoMask, initRecord.asFloat(RecordFieldId("peakresidualnomask")));
		itsMadRMS = max( itsMadRMS, initRecord.asFloat(RecordFieldId("madrms")) );
		
		///itsMaskSum += initRecord.asFloat(RecordFieldId("masksum"));
		/*
		  It has been reset to -1.0.
		  If no masks have changed, it should remain at -1.0
		  If any mask has changed, the sum will come in, and should be added to this.
		 */
		Float thismasksum = initRecord.asFloat(RecordFieldId("masksum"));
		if( thismasksum != -1.0 )
		  {
		    if ( itsMaskSum == -1.0 ) itsMaskSum = thismasksum;
		    else itsMaskSum += thismasksum;
		  }

		if ( itsPrevPeakResidual == -1.0 ) itsPrevPeakResidual = itsPeakResidual;
		if ( itsPrevPeakResidualNoMask == -1.0 ) itsPrevPeakResidualNoMask = itsPeakResidualNoMask;

	}


	void SIIterBot_state::mergeCycleExecutionRecord( Record& execRecord ){
		std::lock_guard<std::recursive_mutex> guard(recordMutex);  

		LogIO os( LogOrigin("SIIterBot_state",__FUNCTION__,WHERE) );

		mergeMinorCycleSummary( execRecord.asArrayDouble( RecordFieldId("summaryminor")) );

		itsIterDone += execRecord.asInt(RecordFieldId("iterdone"));

		itsMaxCycleIterDone = max( itsMaxCycleIterDone, execRecord.asInt(RecordFieldId("maxcycleiterdone")) );

		itsMinorCyclePeakResidual = max( itsPeakResidual, execRecord.asFloat(RecordFieldId("peakresidual")) );
  
		itsUpdatedModelFlag |=execRecord.asBool( RecordFieldId("updatedmodelflag") );

		os << "Completed " << itsIterDone << " iterations." << LogIO::POST;
		//with peak residual "<< itsPeakResidual << LogIO::POST;
	}

	void SIIterBot_state::mergeMinorCycleSummary( const Array<Double>& summary ){
		std::lock_guard<std::recursive_mutex> guard(recordMutex);  
    
		IPosition cShp = itsSummaryMinor.shape();
		IPosition nShp = summary.shape();

		if( cShp.nelements() != 2 || cShp[0] != itsNSummaryFields ||
			nShp.nelements() != 2 || nShp[0] != itsNSummaryFields ) 
			throw(AipsError("Internal error in shape of global minor-cycle summary record"));

		itsSummaryMinor.resize( IPosition( 2, itsNSummaryFields, cShp[1]+nShp[1] ) ,true );

		for (unsigned int row = 0; row < nShp[1]; row++) {
			// iterations done
			itsSummaryMinor( IPosition(2,0,cShp[1]+row) ) = itsIterDone + summary(IPosition(2,0,row));  
			// peak residual
			itsSummaryMinor( IPosition(2,1,cShp[1]+row) ) = summary(IPosition(2,1,row)); 
			// model flux
			itsSummaryMinor( IPosition(2,2,cShp[1]+row) ) = summary(IPosition(2,2,row));
			// cycle threshold
			itsSummaryMinor( IPosition(2,3,cShp[1]+row) ) = summary(IPosition(2,3,row)); 
			// mapper id
			itsSummaryMinor( IPosition(2,4,cShp[1]+row) ) = summary(IPosition(2,4,row)); 
			// chunk id (channel/stokes)
			itsSummaryMinor( IPosition(2,5,cShp[1]+row) ) = summary(IPosition(2,5,row)); 
		}
	}
  
#ifdef INTERACTIVE_ITERATION
	void SIIterBot_state::controlUpdate( const std::map<std::string,DBus::Variant>& updatedParams ) {
		Record controlRecord=dbus::toRecord(updatedParams);
		setControlsFromRecord(controlRecord);
		{
			std::lock_guard<std::mutex> lock(interactionMutex);
			updateNeeded=true;
		}
	}
#endif

	void SIIterBot_state::interactionComplete() {
	  cout << "UU : Interaction completed" << endl;
		changePauseFlag(false);
		itsInteractiveIterDone = 0;    
		{
			std::lock_guard<std::mutex> lock(interactionMutex);
			interactionPending=false;
			updateNeeded=true;
		}
		interactionCond.notify_all();
	}

	std::string SIIterBot_state::getDescription( ) {
		std::lock_guard<std::recursive_mutex> guard(descriptionMutex);
		return itsDescription;
	}

	void SIIterBot_state::setDescription( const std::string &value ) {
		std::lock_guard<std::recursive_mutex> guard(descriptionMutex);    
		itsDescription = value;
	}

#ifdef INTERACTIVE_ITERATION
	std::map<std::string,DBus::Variant> SIIterBot_state::getDetails( ) {
		return dbus::fromRecord(getDetailsRecord());
	}
#endif

	void SIIterBot_state::pushDetails() {
#ifdef INTERACTIVE_ITERATION
/*FIXME    detailUpdate(fromRecord(getDetailsRecord())); */
#endif
	}

	void SIIterBot_state::pushSummary( ) {
		//std::lock_guard<std::recursive_mutex> guard(countMutex);
		std::cout << __FUNCTION__ << "executing" << std::endl;
		//    summaryUpdate();
	}

	DBus::Variant SIIterBot_state::getSummary( ) {
		std::cout << __FUNCTION__ << " executing" << std::endl;
		return DBus::Variant( );
	}

	int SIIterBot_state::getNumberOfControllers( ) {
		std::lock_guard<std::recursive_mutex> guard(recordMutex);    
		return itsControllerCount;
	}

	bool SIIterBot_state::incrementController( ) {
		std::lock_guard<std::recursive_mutex> guard(recordMutex);    
		itsControllerCount++;
		return true;
	}

	bool SIIterBot_state::decrementController( ) {
		std::lock_guard<std::recursive_mutex> guard(recordMutex);    
		itsControllerCount--;
		return true;
	} 

  /* ------------ End of runtime parameter getters -------- */

	void SIIterBot_state::incrementMajorCycleCount( ) {
		std::lock_guard<std::recursive_mutex> guard(recordMutex);
		itsPrevMajorCycleCount = itsMajorDone;
		itsMajorDone++;

		/* Interactive iteractions update */ 
		itsInteractiveIterDone += itsMaxCycleIterDone;
		
		resetMinorCycleInitInfo();
	}

	void SIIterBot_state::resetMinorCycleInitInfo( ) {
		/* Get ready to do the minor cycle */
		itsPeakResidual = 0;
		itsPeakResidualNoMask = 0;
		itsMaxPsfSidelobe = 0;
		itsMaxCycleIterDone = 0;
		itsMaskSum = -1.0;
	}

	Int SIIterBot_state::getMajorCycleCount( ) {
		std::lock_guard<std::recursive_mutex> guard(recordMutex);
		return itsMajorDone;
	}
  
	Record SIIterBot_state::getSummaryRecord( ) {
		LogIO os( LogOrigin("SIIterBot_state",__FUNCTION__,WHERE) );
		std::lock_guard<std::recursive_mutex> guard(recordMutex);
		Record returnRecord = getDetailsRecord();
		returnRecord.define( RecordFieldId("summaryminor"), itsSummaryMinor);
		returnRecord.define( RecordFieldId("summarymajor"), itsSummaryMajor);
		return returnRecord;
	}

  // TODO : Optimize this storage and resizing ? Or call this only now and then... ?

	void SIIterBot_state::addSummaryMajor( ) {
		std::lock_guard<std::recursive_mutex> guard(recordMutex);

		IPosition shp = itsSummaryMajor.shape();
		if( shp.nelements() != 1 ) 
			throw(AipsError("Internal error in shape of major-cycle summary record"));

		itsSummaryMajor.resize( IPosition( 1, shp[0]+1 ) , true );
		itsSummaryMajor( IPosition(1, shp[0] ) ) = itsIterDone;
	}
  
	void SIIterBot_state::updateCycleThreshold( ) {
		std::lock_guard<std::recursive_mutex> guard(recordMutex);    
    
		Float psffraction = itsMaxPsfSidelobe * itsCycleFactor;
    
		psffraction = max(psffraction, itsMinPsfFraction);
		psffraction = min(psffraction, itsMaxPsfFraction);
    
		itsCycleThreshold = itsPeakResidual * psffraction;
		pushDetails();
	}

//   Float SIIterBot_state::getMaxPsfSidelobe()
//   {
//     std::lock_guard<std::recursive_mutex> guard(recordMutex);    
//     return itsMaxPsfSidelobe;
//   }

//   void SIIterBot_state::setMaxPsfSidelobe(Float maxSidelobe)
//   {
//     std::lock_guard<std::recursive_mutex> guard(recordMutex);    
//     itsMaxPsfSidelobe = maxSidelobe;
//   }

//   Float SIIterBot_state::getMaxPsfFraction()
//   {
//     std::lock_guard<std::recursive_mutex> guard(recordMutex);    
//     return itsMaxPsfFraction;
//   }

//   void SIIterBot_state::setMaxPsfFraction(Float maxPsfFraction)
//   {
//     std::lock_guard<std::recursive_mutex> guard(recordMutex);    
//     itsMaxPsfFraction = maxPsfFraction;
//   }

//   Float SIIterBot_state::getMinPsfFraction()
//   {
//     std::lock_guard<std::recursive_mutex> guard(recordMutex);    
//     return itsMinPsfFraction;
//   }

//   void SIIterBot_state::setMinPsfFraction(Float minPsfFraction)
//   {
//     std::lock_guard<std::recursive_mutex> guard(recordMutex);    
//     itsMinPsfFraction = minPsfFraction;
//   }

	Record SIIterBot_state::getDetailsRecord(){
		LogIO os( LogOrigin("SIIterBot_state",__FUNCTION__,WHERE) );
		std::lock_guard<std::recursive_mutex> guard(recordMutex);
		Record returnRecord;

		/* Control Variables */
		returnRecord.define(RecordFieldId("niter"), itsNiter);
		returnRecord.define(RecordFieldId("cycleniter"), itsCycleNiter);
		returnRecord.define(RecordFieldId("interactiveniter"),itsInteractiveNiter);

		returnRecord.define( RecordFieldId("threshold"),  itsThreshold);    
		if( itsIsCycleThresholdAuto == true )  updateCycleThreshold();
		itsIsCycleThresholdAuto = true; /* Reset this, for the next round */

		returnRecord.define( RecordFieldId("cyclethreshold"),itsCycleThreshold);
		returnRecord.define( RecordFieldId("interactivethreshold"), itsInteractiveThreshold);  

		returnRecord.define( RecordFieldId("loopgain"),  itsLoopGain);
		returnRecord.define( RecordFieldId("cyclefactor"), itsCycleFactor);

		/* Status Reporting Variables */
		returnRecord.define( RecordFieldId("iterdone"),  itsIterDone);
		returnRecord.define( RecordFieldId("cycleiterdone"),  itsMaxCycleIterDone);
		returnRecord.define( RecordFieldId("interactiveiterdone"), 
							 itsInteractiveIterDone + itsMaxCycleIterDone);
    
		returnRecord.define( RecordFieldId("nmajordone"),  itsMajorDone);
		returnRecord.define( RecordFieldId("maxpsfsidelobe"), itsMaxPsfSidelobe);
		returnRecord.define( RecordFieldId("maxpsffraction"), itsMaxPsfFraction);
		returnRecord.define( RecordFieldId("minpsffraction"), itsMinPsfFraction);
		returnRecord.define( RecordFieldId("interactivemode"), itsInteractiveMode);

		returnRecord.define( RecordFieldId("stopcode"), itsStopCode);

		/* report clean's state */
		returnRecord.define( RecordFieldId("cleanstate"), itsStopFlag ? "stopped" :
		                                                  itsPauseFlag ? "paused" : "running" );
		return returnRecord;
	}

	void SIIterBot_state::changeNiter( Int niter ) {
		std::lock_guard<std::recursive_mutex> guard(recordMutex);    
		itsNiter = niter;
		//cout << "UUU : change niter : " << niter << endl;
	}

	void SIIterBot_state::changeCycleNiter( Int cycleniter ) {
		std::lock_guard<std::recursive_mutex> guard(recordMutex);
		if (cycleniter <= 0)  
			itsCycleNiter = itsNiter;
		else
			itsCycleNiter = cycleniter;
	}

	void SIIterBot_state::changeInteractiveNiter( Int interactiveNiter ) {
		std::lock_guard<std::recursive_mutex> guard(recordMutex);    
		itsInteractiveNiter = interactiveNiter;
	}

	void SIIterBot_state::changeThreshold( Float threshold ) {
		std::lock_guard<std::recursive_mutex> guard(recordMutex);    

		if( threshold == -1.0 ) {
		  itsThreshold = 0.0;
		  itsIsThresholdAuto = true;
		}
		else { 
		  itsThreshold = threshold; 
		  itsIsThresholdAuto = false;
		}
	}

	void SIIterBot_state::changeCycleThreshold( Float cyclethreshold ) {
		std::lock_guard<std::recursive_mutex> guard(recordMutex);
		itsCycleThreshold = cyclethreshold;
		itsIsCycleThresholdAuto = false;
	}

	void SIIterBot_state::changeInteractiveThreshold( Float interactivethreshold ) {
		std::lock_guard<std::recursive_mutex> guard(recordMutex);
		itsInteractiveThreshold = interactivethreshold;
	}

	void SIIterBot_state::changeLoopGain( Float loopgain ) {
		std::lock_guard<std::recursive_mutex> guard(recordMutex);
		itsLoopGain = loopgain;
	}

	void SIIterBot_state::changeCycleFactor( Float cyclefactor ) {
		std::lock_guard<std::recursive_mutex> guard(recordMutex);    
		itsCycleFactor = cyclefactor;
	}

	void SIIterBot_state::changeInteractiveMode( const bool& interactiveEnabled ) {
		std::lock_guard<std::recursive_mutex> guard(recordMutex);    
		itsInteractiveMode = interactiveEnabled;
	}

	void SIIterBot_state::changePauseFlag( const bool& pauseEnabled ){
		std::lock_guard<std::recursive_mutex> guard(recordMutex);    
		itsPauseFlag = pauseEnabled;
		cout << "UUU : changed pause flag to " << pauseEnabled << endl;
	}

	void SIIterBot_state::changeStopFlag(const bool& stopEnabled){
		std::lock_guard<std::recursive_mutex> guard(recordMutex);    
		itsStopFlag = stopEnabled;
	}

	void SIIterBot_state::changeMinPsfFraction( Float minpsffraction ) {
		std::lock_guard<std::recursive_mutex> guard(recordMutex);    
		itsMinPsfFraction = minpsffraction;
	}

	void SIIterBot_state::changeMaxPsfFraction( Float maxpsffraction ) {
		std::lock_guard<std::recursive_mutex> guard(recordMutex);    
		itsMaxPsfFraction = maxpsffraction;
	}

	void SIIterBot_state::setControlsFromRecord( Record &recordIn ) {
		LogIO os( LogOrigin("SIIterBot_state",__FUNCTION__,WHERE) );
		std::lock_guard<std::recursive_mutex> guard(recordMutex);

		/* Note it is important that niter get set first as we catch
		   negative values in the cycleniter, and set it equal to niter */
		if (recordIn.isDefined("niter"))
			changeNiter(recordIn.asInt( RecordFieldId("niter")));

		if (recordIn.isDefined("cycleniter"))
			changeCycleNiter(recordIn.asInt( RecordFieldId("cycleniter")));

		if (recordIn.isDefined("interactiveniter"))
			changeInteractiveNiter( recordIn.asInt(RecordFieldId("interactiveniter")) );
    
		if (recordIn.isDefined("threshold")) 
		    changeThreshold( readThreshold( recordIn, "threshold" ) );
		
		if (recordIn.isDefined("cyclethreshold")) 
		    changeCycleThreshold( readThreshold( recordIn, "cyclethreshold" ) );

		if (recordIn.isDefined("interactivethreshold")) 
			changeInteractiveThreshold(recordIn.asFloat(RecordFieldId("interactivethreshold")));

		if (recordIn.isDefined("loopgain")) 
			changeLoopGain(recordIn.asDouble( RecordFieldId("loopgain")));;

		if (recordIn.isDefined("cyclefactor"))
			changeCycleFactor(recordIn.asFloat( RecordFieldId("cyclefactor")));

		if (recordIn.isDefined("interactive"))
			changeInteractiveMode(recordIn.asBool(RecordFieldId("interactive")));

		if (recordIn.isDefined("minpsffraction"))
			changeMinPsfFraction(recordIn.asFloat( RecordFieldId("minpsffraction")));

		if (recordIn.isDefined("maxpsffraction"))
			changeMaxPsfFraction(recordIn.asFloat( RecordFieldId("maxpsffraction")));

		//		printOut("After Setting : ", false);

	}

        Float SIIterBot_state::readThreshold( Record recordIn, String id )  {
		LogIO os( LogOrigin("SIIterBot_state",__FUNCTION__,WHERE) );
		std::lock_guard<std::recursive_mutex> guard(recordMutex);    
		// Threshold can be a variant, either Float or String(with units).
		Float fthresh=0.0;
		// If a number, treat it as a number in units of Jy.
		if( recordIn.dataType(id) == TpFloat || 
		    recordIn.dataType(id) == TpDouble || 
		    recordIn.dataType(id) == TpInt )
		  { fthresh = recordIn.asFloat( RecordFieldId(id)); }
		// If a string, try to convert to a Quantity
		else if( recordIn.dataType(id) == TpString )
		  {
		    Quantity thresh; 
		    // If it cannot be converted to a Quantity.... complain, and use zero.
		    if( ! casacore::Quantity::read( thresh, recordIn.asString( RecordFieldId(id) ) ) )
		      {os << LogIO::WARN << "Cannot parse threshold value. Setting to zero." << LogIO::POST;  
			fthresh=0.0;}
		    // If converted to Quantity, get value in Jy. 
		    // ( Note : This does not check for wrong units, e.g. if the user says '100m' ! )
		    else { fthresh = thresh.getValue(Unit("Jy")); }
		  }
		// If neither valid datatype, print a warning and use zero.
		else {os << LogIO::WARN << id << " is neither a number nor a string Quantity. Setting to zero." << LogIO::POST;
		  fthresh=0.0; }

		return fthresh;
	}
  
	/* Print out contents of the IterBot. For debugging. */
	void SIIterBot_state::printOut( String prefix, Bool verbose ) {
		if( verbose == true ) {
			cout << prefix << " : " 
				 << " ItsNiter=" << itsNiter
				 << " itsCycleNiter=" << itsCycleNiter
				 << " itsInteractiveNiter=" << itsInteractiveNiter
				 << " itsThreshold=" << itsThreshold
				 << " itsCycleThreshold=" << itsCycleThreshold
				 << " itsInteractiveThreshold=" << itsInteractiveThreshold
				 << " itsCycleFactor=" << itsCycleFactor
				 << " itsLoopGain=" << itsLoopGain
				 << " itsStopFlag=" << itsStopFlag
				 << " itsPauseFlag=" << itsPauseFlag
				 << " itsInteractiveMode=" << itsInteractiveMode
				 << " itsUpdatedModelFlag=" << itsUpdatedModelFlag
				 << " itsPeakResidual=" << itsPeakResidual
				 << " itsMaxPsfSidelobe=" << itsMaxPsfSidelobe
				 << " itsMinPsfFraction=" << itsMinPsfFraction
				 << " itsMaxPsfFraction=" << itsMaxPsfFraction
				 << " itsIterdone=" << itsIterDone
				 << " itsInteractiveIterDone=" << itsInteractiveIterDone
				 << " itsMaxCycleIterDone=" << itsMaxCycleIterDone
				 << " itsMajorDone=" << itsMajorDone
			         << " itsStopCode=" << itsStopCode 
				 << endl;
		} else {
			cout << prefix << " : " 
				 << " ItsNiter=" << itsNiter
				 << " itsCycleNiter=" << itsCycleNiter
				 << " itsThreshold=" << itsThreshold
				 << " itsCycleThreshold=" << itsCycleThreshold
				 << " itsStopFlag=" << itsStopFlag
				 << " itsPeakResidual=" << itsPeakResidual
				 << " itsIterdone=" << itsIterDone
				 << endl;
		}

	}

	SIIterBot_adaptor::SIIterBot_adaptor( SHARED_PTR<SIIterBot_state> s, const std::string &bus_name, const std::string &object_path) :
#ifdef INTERACTIVE_ITERATION
				dbus::address(bus_name),
				DBus::ObjectAdaptor( DBusSession::instance().connection( ), object_path ),
#endif
				state(s) {
#ifdef INTERACTIVE_ITERATION
		state->acceptCallbacks(this);
#else
                (void)bus_name;
		(void)object_path;
#endif
	}

	SIIterBot_adaptor::~SIIterBot_adaptor() {
		fprintf( stderr, ">>>>>>\t\tSIIterBot_adaptor::~SIIterBot_adaptor(0x%p)\n", this );
		fflush( stderr );
#ifdef INTERACTIVE_ITERATION
		state->denyCallbacks(this);
		disconnect();
		
#endif
	}

} //# NAMESPACE CASA - END

