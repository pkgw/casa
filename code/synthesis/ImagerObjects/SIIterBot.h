//# SIIterBot.h: This file contains the interface definition SIIterBot class
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

#ifndef SYNTHESIS_SIITERBOT
#define SYNTHESIS_SIITERBOT

#include <casadbus/utilities/BusAccess.h>
// .casarc interface
#include <casa/System/AipsrcValue.h>

// System utilities (for profiling macros)
#include <casa/OS/HostInfo.h>
#include <sys/time.h>
#if defined(DBUS_CPP)
#include <dbus-cpp/dbus.h> /*for DBus::Variant... probably can be removed with *_adaptor class*/
#else
#include <dbus-c++/dbus.h>
#endif


// Boost Libraries for mutex and noncopyable semantics
#include <mutex>
#include <condition_variable>

// Include files for the DBus Service
//#include <casadbus/interfaces/SynthImager.adaptor.h>

#ifdef INTERACTIVE_ITERATION
#include <casadbus/interfaces/SynthImager.adaptor.h>
#endif

namespace casacore{

	class Record;
}

namespace casa { //# NAMESPACE CASA - BEGIN


	class SIIterBot_adaptor;

	class SIIterBot_callback {
		public:
			SIIterBot_callback( ) : adaptor(0) { }
			~SIIterBot_callback( ) { }
			void interactionRequired(bool);

			void addHandler( SIIterBot_adaptor* );
			void removeHandler( SIIterBot_adaptor* );

		private:
			std::recursive_mutex mutex;
			SIIterBot_adaptor *adaptor;

	};

	class SIIterBot_state {
		private:
			// make SIIterBot_state uncopyable...
			SIIterBot_state( const SIIterBot_state & );
			SIIterBot_state &operator=( const SIIterBot_state & );

		public:
			SIIterBot_state( SHARED_PTR<SIIterBot_callback> );
			~SIIterBot_state( );

			/****
			***** allow or deny callbacks which are funneled through dbus,
			***** a shared pointer is explicitly NOT used here to avoid a cycle.
			****/
			void acceptCallbacks( SIIterBot_adaptor *siba ) { callback->addHandler(siba); }
			void denyCallbacks( SIIterBot_adaptor *siba ) { callback->removeHandler(siba); }

			/* Wait for an Interactive Clean Cycle */
			virtual bool interactiveInputRequired();
			virtual void waitForInteractiveInput(); 

			// virtual bool majorCycleRequired(casacore::Float currentPeakResidual);

	                virtual int cleanComplete(casacore::Bool lastcyclecheck=casacore::False);
    
			/* --- Functions for interacting with Minor Cycle Control --- */
			virtual casacore::Record getMinorCycleControls();
			virtual void   mergeCycleInitializationRecord(casacore::Record&);
			virtual void   mergeCycleExecutionRecord(casacore::Record&);
    

			//void mergeSubIterBot(SISubIterBot& subIterBot);

			/* ------ Begin functions for runtime parameter modification ------ */
			/* These are the control variables.  We can set them either through these
			   method or through the casacore::Record interface */
			void changeNiter( casacore::Int niter );
			void changeCycleNiter( casacore::Int cycleniter );
			void changeInteractiveNiter(casacore::Int interactiveniter );

			void changeThreshold( casacore::Float threshold );
			void changeCycleThreshold( casacore::Float cyclethreshold );
			void changeInteractiveThreshold( casacore::Float cyclethreshold );
			
			void changeLoopGain(casacore::Float loopgain );
			void changeCycleFactor( casacore::Float cyclefactor);

	                void changeMinPsfFraction(casacore::Float minpsffraction);
	                void changeMaxPsfFraction(casacore::Float maxpsffraction);

			void changeInteractiveMode(const bool& interactiveEnabled);
			void changePauseFlag(const bool& pauseEnabled);
			void changeStopFlag(const bool& stopEnabled);

			/* As a convience the controls can also be updated from a Record.  
			   The following fields are supported:
			   - niter
			   - cycleniter
			   - interactiveniter
			   - threshold
			   - cyclethreshold
			   - interactivethreshold
			   - loopgain
			   - cyclefactor
			*/
			void setControlsFromRecord(casacore::Record &recordIn);
         	        casacore::Float readThreshold( casacore::Record recordIn, casacore::String id );

			virtual casacore::Record getDetailsRecord();
			//casacore::Record getSubIterBotRecord();

			/* ------ END Functions for runtime parameter modification -------*/

			/* Call these functions to keep track of cycles */

			/* Note:  Incrementing the Major cycle will reset the cycleIterDone */
			void incrementMajorCycleCount();
	                void resetMinorCycleInitInfo();

			casacore::Int getMajorCycleCount();



			/* This will calculate and set a new cycle threshold based
			   on the Peak Residual and the current values of the PSF values.*/
			void updateCycleThreshold();
    
			/* Values for the control of the cycle threshold */
			/*     void setMaxPsfSidelobe( casacore::Float maxpsfsidelobe ); */
			/*     casacore::Float getMaxPsfSidelobe(); */

			/*     void setMaxPsfFraction(casacore::Float maxpsffraction ); */
			/*     casacore::Float getMaxPsfFraction(); */

			/*     void setMinPsfFraction(casacore::Float minpsffraction ); */
			/*     casacore::Float getMinPsfFraction(); */
			//casacore::Int getRemainingNiter();
			//casacore::Int getCompletedNiter();

			void addSummaryMajor();

			/* -------- DBus Interface ------------- */

			/* Publish the current details from the Iterbot to all clients */
			void pushDetails();

			/* Publish the current summary from the Iterbot to all clients */
			void pushSummary();

			/* These are the fuctions we provide through the dbus interface */
			bool incrementController();
			bool decrementController();
			int  getNumberOfControllers();

			std::string getDescription( );
			void setDescription( const std::string &value );

#ifdef INTERACTIVE_ITERATION
			std::map<std::string,DBus::Variant> getDetails();
#endif
			DBus::Variant getSummary();
    
			void interactionComplete();

#ifdef INTERACTIVE_ITERATION
			void controlUpdate(const std::map<std::string, DBus::Variant>& parameters);
#endif
			//// START /// Functions for runtime parameter modification

			/* Methods to get the state of the iterbot as a casacore::Record*/

			/* This record has all of the fields associated with the detail record 
			   but adds
			   - summaryminor
			   - summarymajor
			*/
			virtual casacore::Record getSummaryRecord();

		protected:

			virtual void mergeMinorCycleSummary(const casacore::Array<casacore::Double>& summary);

			/* For testing/debugging. Print out the contents of the iterbot */
			void printOut(casacore::String prefix, casacore::Bool verbose); 

			std::string itsDescription;

			casacore::Float itsMinPsfFraction;
			casacore::Float itsMaxPsfFraction;
			casacore::Float itsMaxPsfSidelobe;
			casacore::Float itsPeakResidual;
	                casacore::Float itsPrevPeakResidual;
	                casacore::Float itsInitPeakResidual;
	                casacore::Float itsMinPeakResidual;
	                casacore::Float itsMinorCyclePeakResidual;

			casacore::Float itsPeakResidualNoMask;
			casacore::Float itsPrevPeakResidualNoMask;
			casacore::Float itsMinPeakResidualNoMask;
	  
	                casacore::Float itsMadRMS;
	                casacore::Float itsMaskSum;
	  
 	                casacore::Int itsPrevMajorCycleCount;

			/* The number of Controllers Currently Connected */
			int    itsControllerCount;


			/* A recursive mutex which provides access control to the
			   underlying state variables
			*/
			std::recursive_mutex recordMutex;
			std::recursive_mutex descriptionMutex; /** protects itsDescription **/

			/* Control Variables */
			casacore::Int    itsNiter;
			casacore::Int    itsCycleNiter;
			casacore::Int    itsInteractiveNiter;

			casacore::Float itsThreshold;
			casacore::Float itsCycleThreshold;
			casacore::Float itsInteractiveThreshold;
	  
	        casacore::Bool itsIsCycleThresholdAuto;
	  casacore::Bool itsIsThresholdAuto;

			casacore::Float itsCycleFactor;
			casacore::Float itsLoopGain;
    
			casacore::Bool  itsStopFlag;
			casacore::Bool  itsPauseFlag;
			casacore::Bool  itsInteractiveMode;
			casacore::Bool  itsUpdatedModelFlag;

			casacore::Int   itsIterDone;
			casacore::Int   itsInteractiveIterDone;
			casacore::Int   itsMaxCycleIterDone;
			casacore::Int   itsMajorDone;
	        casacore::Int   itsStopCode;

			/*
			  A condition variable used when we're waiting for interaction to 
			  complete
			*/
			bool                      interactionPending;
			bool                      updateNeeded;
			std::condition_variable   interactionCond; 
			std::mutex                interactionMutex; 


			/* Summary variables */
			casacore::Int itsNSummaryFields;
			casacore::Array<casacore::Double> itsSummaryMinor;
			casacore::Array<casacore::Int>    itsSummaryMajor;

			SHARED_PTR<SIIterBot_callback> callback;
	};

	class SIIterBot_adaptor
#ifdef INTERACTIVE_ITERATION
		: private dbus::address,
		  private edu::nrao::casa::SynthesisImager_adaptor,
		  public DBus::IntrospectableAdaptor,
		  public DBus::ObjectAdaptor
#endif
		{
			public:
				SIIterBot_adaptor( SHARED_PTR<SIIterBot_state> state, const std::string &bus_name, const std::string &object_path );
				~SIIterBot_adaptor();

				bool incrementController( )	{ return state->incrementController( ); }
				bool decrementController( )	{ return state->decrementController( ); }

		  void interactionRequired( const bool &val ) {
#ifdef INTERACTIVE_ITERATION
					edu::nrao::casa::SynthesisImager_adaptor::interactionRequired( val );
#else
					(void)val;  // To get the compiler to not warn...
#endif
				}
				void controlUpdate(const std::map< std::string, ::DBus::Variant >& newParams)
											{
#ifdef INTERACTIVE_ITERATION
												state->controlUpdate(newParams);
#else 
												(void)newParams;// to avoid compiler warning
											       
#endif
											}
				void interactionComplete( )
											{ state->interactionComplete( ); }
				void changePauseFlag(const bool& pauseState)
											{ state->changePauseFlag(pauseState); }
				void changeStopFlag(const bool& stopState)
											{ state->changeStopFlag(stopState); }
				void changeInteractiveMode(const bool& interactiveMode)
											{ state->changeInteractiveMode(interactiveMode); }
				std::string getDescription( )
											{ return state->getDescription( ); }
				std::map< std::string, ::DBus::Variant > getDetails( )
											{
#ifdef INTERACTIVE_ITERATION
												return state->getDetails( );
#else
												return std::map< std::string, ::DBus::Variant >( );
#endif
											}
				::DBus::Variant getSummary( )
											{ return state->getSummary( ); }

			private:
				SHARED_PTR<SIIterBot_state> state;

		};
    
} //# NAMESPACE CASA - END

#endif /* FLAGDATAHANDLER_H_ */
