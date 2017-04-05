//# SynthesisImagerVi2.cc: Implementation of SynthesisImager.h
//# Copyright (C) 1997-2016
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This program is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by the Free
//# Software Foundation; either version 2 of the License, or (at your option)
//# any later version.
//#
//# This program is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
//# more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with this program; if not, write to the Free Software Foundation, Inc.,
//# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id$
#include <casa/Exceptions/Error.h>
#include <casa/iostream.h>
#include <casa/sstream.h>

#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayLogical.h>


#include <casa/Logging.h>
#include <casa/Logging/LogIO.h>
#include <casa/Logging/LogMessage.h>
#include <casa/Logging/LogSink.h>
#include <casa/Logging/LogMessage.h>
#include <casa/System/ProgressMeter.h>

#include <casa/OS/DirectoryIterator.h>
#include <casa/OS/File.h>
#include <casa/OS/HostInfo.h>
#include <casa/OS/Path.h>
//#include <casa/OS/Memory.h>

#include <lattices/LRegions/LCBox.h>

#include <measures/Measures/MeasTable.h>

#include <ms/MeasurementSets/MSHistoryHandler.h>
#include <ms/MeasurementSets/MeasurementSet.h>
#include <ms/MSSel/MSSelection.h>


#include <synthesis/ImagerObjects/SIIterBot.h>
#include <synthesis/ImagerObjects/SynthesisImagerVi2.h>

#include <synthesis/ImagerObjects/SynthesisUtilMethods.h>
#include <synthesis/ImagerObjects/SIImageStore.h>
#include <synthesis/ImagerObjects/SIImageStoreMultiTerm.h>

#include <synthesis/MeasurementEquations/ImagerMultiMS.h>
#include <synthesis/MeasurementEquations/VPManager.h>
#include <imageanalysis/Utilities/SpectralImageUtil.h>
#include <msvis/MSVis/MSUtil.h>
#include <msvis/MSVis/VisSetUtil.h>
#include <msvis/MSVis/VisImagingWeight.h>

#include <synthesis/TransformMachines2/GridFT.h>
#include <synthesis/TransformMachines2/WPConvFunc.h>
#include <synthesis/TransformMachines2/WProjectFT.h>
#include <synthesis/TransformMachines2/VisModelData.h>
#include <synthesis/TransformMachines2/AWProjectFT.h>
#include <synthesis/TransformMachines2/HetArrayConvFunc.h>
#include <synthesis/TransformMachines2/MosaicFTNew.h>
#include <synthesis/TransformMachines2/MultiTermFTNew.h>
#include <synthesis/TransformMachines2/AWProjectWBFTNew.h>
#include <synthesis/TransformMachines2/AWConvFunc.h>
#include <synthesis/TransformMachines2/AWConvFuncEPJones.h>
#include <synthesis/TransformMachines2/NoOpATerm.h>
#include <synthesis/TransformMachines/WProjectFT.h>
#include <casadbus/viewer/ViewerProxy.h>
#include <casadbus/plotserver/PlotServerProxy.h>
#include <casacore/casa/Utilities/Regex.h>
#include <casacore/casa/OS/Directory.h>
#include <msvis/MSVis/VisibilityIteratorImpl2.h>
//#include <casadbus/utilities/BusAccess.h>
//#include <casadbus/session/DBusSession.h>

#include <sys/types.h>
#include <unistd.h>
#include <iomanip>


using namespace std;

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

  SynthesisImagerVi2::SynthesisImagerVi2() : SynthesisImager(), vi_p(0), fselections_p(nullptr) {

    mss_p.resize(0);
  }

  SynthesisImagerVi2::~SynthesisImagerVi2(){
    for (uInt k=0; k < mss_p.nelements(); ++k){
      if(mss_p[k])
	delete mss_p[k];
    }
  }

  Bool SynthesisImagerVi2::selectData(const SynthesisParamsSelect& selpars){
 LogIO os( LogOrigin("SynthesisImagerVi2","selectData",WHERE) );

    try
      {

    //Respect the readonly flag...necessary for multi-process access
    MeasurementSet thisms(selpars.msname, TableLock(TableLock::AutoNoReadLocking),
				selpars.readonly ? Table::Old : Table::Update);
    thisms.setMemoryResidentSubtables (MrsEligibility::defaultEligible());

    useScratch_p=selpars.usescratch;
    readOnly_p = selpars.readonly;
    //    cout << "**************** usescr : " << useScratch_p << "     readonly : " << readOnly_p << endl;
    //if you want to use scratch col...make sure they are there
    if(selpars.usescratch && !selpars.readonly){
      VisSetUtil::addScrCols(thisms, true, false, true, false);
      refim::VisModelData::clearModel(thisms);
    }

    if(!selpars.incrmodel && !selpars.usescratch && !selpars.readonly)
      refim::VisModelData::clearModel(thisms, selpars.field, selpars.spw);

    unlockMSs();


    os << "MS : " << selpars.msname << " | ";

    //Some MSSelection 
    //If everything is empty (which is valid) it will throw an exception..below
    //So make sure the main defaults are not empy i.e field and spw
    MSSelection thisSelection;
    if(selpars.field != ""){
      thisSelection.setFieldExpr(selpars.field);
      os << "Selecting on fields : " << selpars.field << " | " ;//LogIO::POST;
    }else
      thisSelection.setFieldExpr("*");
    if(selpars.spw != ""){
	thisSelection.setSpwExpr(selpars.spw);
	os << "Selecting on spw :"<< selpars.spw  << " | " ;//LogIO::POST;
    }else
      thisSelection.setSpwExpr("*");
    
    if(selpars.antenna != ""){
      Vector<String> antNames(1, selpars.antenna);
      // thisSelection.setAntennaExpr(MSSelection::nameExprStr( antNames));
      thisSelection.setAntennaExpr(selpars.antenna);
      os << "Selecting on antenna names : " << selpars.antenna << " | " ;//LogIO::POST;
	
    }            
    if(selpars.timestr != ""){
	thisSelection.setTimeExpr(selpars.timestr);
	os << "Selecting on time range : " << selpars.timestr << " | " ;//LogIO::POST;	
      }
    if(selpars.uvdist != ""){
      thisSelection.setUvDistExpr(selpars.uvdist);
      os << "Selecting on uvdist : " << selpars.uvdist << " | " ;//LogIO::POST;	
    }
    if(selpars.scan != ""){
      thisSelection.setScanExpr(selpars.scan);
      os << "Selecting on scan : " << selpars.scan << " | " ;//LogIO::POST;	
    }
    if(selpars.obs != ""){
      thisSelection.setObservationExpr(selpars.obs);
      os << "Selecting on Observation Expr : " << selpars.obs << " | " ;//LogIO::POST;	
    }
    if(selpars.state != ""){
      thisSelection.setStateExpr(selpars.state);
      os << "Selecting on Scan Intent/State : " << selpars.state << " | " ;//LogIO::POST;	
    }
    // if(selpars.taql != ""){
    // 	thisSelection.setTaQLExpr(selpars.taql);
    // 	os << "Selecting via TaQL : " << selpars.taql << " | " ;//LogIO::POST;	
    // }
    os << "[Opened " << (readOnly_p?"in readonly mode":(useScratch_p?"with scratch model column":"with virtual model column"))  << "]" << LogIO::POST;
    TableExprNode exprNode=thisSelection.toTableExprNode(&thisms);
    if(!(exprNode.isNull()))
      {
	mss_p.resize(mss_p.nelements()+1, false, true);
    
	MeasurementSet thisMSSelected0 = MeasurementSet(thisms(exprNode));

	if(selpars.taql != "")
	  {
	    MSSelection mss0;
	    mss0.setTaQLExpr(selpars.taql);
	    os << "Selecting via TaQL : " << selpars.taql << " | " ;//LogIO::POST;	

	    TableExprNode tenWithTaQL=mss0.toTableExprNode(&thisMSSelected0);
	    MeasurementSet thisMSSelected1 = MeasurementSet(thisMSSelected0(tenWithTaQL));
	    //mss4vi_p[mss4vi_p.nelements()-1]=MeasurementSet(thisms(exprNode));
	    mss_p[mss_p.nelements()-1]=new MeasurementSet(thisMSSelected1);
	  }
	else
	  mss_p[mss_p.nelements()-1]=new MeasurementSet(thisMSSelected0);
	  
	os << "  NRows selected : " << (mss_p[mss_p.nelements()-1])->nrow() << LogIO::POST;
	unlockMSs();
      }
    else{
      throw(AipsError("Selection for given MS "+selpars.msname+" is invalid"));
    }
    
    ///Channel selection
    {
      if(!fselections_p) fselections_p=new FrequencySelections();
      Matrix<Int> chanlist = thisSelection.getChanList(mss_p[mss_p.nelements()-1]);
      Matrix<Double> freqList=thisSelection.getChanFreqList(mss_p[mss_p.nelements()-1]);
      //cerr << std::setprecision(12) << "FreqList " << freqList << endl;
      IPosition shape = freqList.shape();
      uInt nSelections = shape[0];
      if(selpars.freqbeg==""){
	   // Going round the problem of CAS-8829
        /*vi::FrequencySelectionUsingChannels channelSelector;

        channelSelector.add(thisSelection, mss_p[mss_p.nelements()-1]);

        fselections_p.add(channelSelector);
        */
        ////////////////////////////
        Double lowfreq;
        Double topfreq;
	//cerr << "chanlist " << chanlist << "\n freqlis " << freqList << endl;
        MFrequency::Types freqFrame=MFrequency::castType(ROMSColumns(*mss_p[mss_p.nelements()-1]).spectralWindow().measFreqRef()(Int(chanlist(0,0))));
        vi::FrequencySelectionUsingFrame channelSelector(freqFrame);
	///temporary variable as we carry that for tunechunk
	selFreqFrame_p=freqFrame;
    	  for(uInt k=0; k < nSelections; ++k){
	    //The getChanfreqList is wrong for beg and end..going round that too.
	    Vector<Double> freqies=ROMSColumns(*mss_p[mss_p.nelements()-1]).spectralWindow().chanFreq()(Int(chanlist(k,0)));
	    Vector<Double> reso=ROMSColumns(*mss_p[mss_p.nelements()-1]).spectralWindow().resolution()(Int(chanlist(k,0)));
            
	    if(freqList(k,3) < 0.0){
	      topfreq=freqies(chanlist(k,1));
	      lowfreq=freqies(chanlist(k,2));
	      //lowfreq=freqList(k,2); //+freqList(k,3)/2.0;
	      //topfreq=freqList(k, 1); //-freqList(k,3)/2.0;
	    }
	    else{
	      lowfreq=freqies(chanlist(k,1));
	      topfreq=freqies(chanlist(k,2));
	      //lowfreq=freqList(k,1); //-freqList(k,3)/2.0;
	      //topfreq=freqList(k, 2); //+freqList(k,3)/2.0;
	    }
	    //cerr << std::setprecision(12) << "Dat lowFreq "<< lowfreq << " topfreq " << topfreq << endl; 
            //channelSelector.add(Int(freqList(k,0)), lowfreq, topfreq);
	    andFreqSelection(mss_p.nelements()-1, Int(freqList(k,0)), lowfreq, topfreq, freqFrame);
          }
    	  //fselections_p->add(channelSelector);
          //////////////////////////////////
      }
      else{

	//////More workaroung CAS-8829
	MFrequency::Types freqFrame=MFrequency::castType(ROMSColumns(*mss_p[mss_p.nelements()-1]).spectralWindow().measFreqRef()(Int(freqList(0,0))));
	
    	  Quantity freq;
    	  Quantity::read(freq, selpars.freqbeg);
    	  Double lowfreq=freq.getValue("Hz");
    	  Quantity::read(freq, selpars.freqend);
    	  Double topfreq=freq.getValue("Hz");
    	 
	  ////Work aroun CAS-8829
	  if(vi_p) 
	    VisBufferUtil::getFreqRangeFromRange(lowfreq, topfreq,  selpars.freqframe, lowfreq,  topfreq, *vi_p, freqFrame);
	  //cerr << "lowFreq "<< lowfreq << " topfreq " << topfreq << endl;
	  vi::FrequencySelectionUsingFrame channelSelector((vi_p ? freqFrame :selpars.freqframe));
    	  for(uInt k=0; k < nSelections; ++k){
	    //cerr << "lowFreq "<< lowfreq << " topfreq " << topfreq << endl;
            //channelSelector.add(Int(freqList(k,0)), lowfreq, topfreq);
	    andFreqSelection((mss_p.nelements()-1), Int(freqList(k,0)), lowfreq, topfreq, vi_p ?freqFrame : selpars.freqframe);
          }
    	  //fselections_p->add(channelSelector);

      }
    }
    writeAccess_p=writeAccess_p && !selpars.readonly;
    createVisSet(writeAccess_p);

    /////// Remove this when the new vi/vb is able to get the full freq range.
    mssFreqSel_p.resize();
    mssFreqSel_p  = thisSelection.getChanFreqList(NULL,true);
   
    //// Set the data column on which to operate
    // TT: added checks for the requested data column existace 
    //    cout << "Using col : " << selpars.datacolumn << endl;
    if( selpars.datacolumn.contains("data") || selpars.datacolumn.contains("obs") ) 
      {    if( thisms.tableDesc().isColumn("DATA") ) { datacol_p = FTMachine::OBSERVED; }
           else { os << LogIO::SEVERE <<"DATA column does not exist" << LogIO::EXCEPTION;}
      }
    else if( selpars.datacolumn.contains("corr") ) { datacol_p = FTMachine::CORRECTED; }
    else { os << LogIO::WARN << "Invalid data column : " << selpars.datacolumn << ". Using corrected (or observed if corrected doesn't exist)" << LogIO::POST;  datacol_p =  FTMachine::CORRECTED;}


    dataSel_p.resize(dataSel_p.nelements()+1, true);

    dataSel_p[dataSel_p.nelements()-1]=selpars;


    unlockMSs();
      }
    catch(AipsError &x)
      {
	unlockMSs();
	throw( AipsError("Error in selectData() : "+x.getMesg()) );
      }

    return true;



  }

  void SynthesisImagerVi2::andFreqSelection(const Int msId, const Int spwId,  const Double freqBeg, const Double freqEnd, const MFrequency::Types frame){
    
   
    Int key=msId;
   
    Bool isDefined=False;
    FrequencySelectionUsingFrame frameSel(frame);
    for (uInt k =0; k<freqBegs_p.size(); ++k){ 
      // cerr <<freqBegs_p[k].first  << " == " << key << " && " << freqSpws_p[k].second<< " == " << spwId << " && " << freqBeg << " < " << freqEnds_p[k].second<< " && " << freqEnd << " > " << freqBegs_p[k].second << endl;
	if((freqBegs_p[k].first == key || key <0 ) && (freqSpws_p[k].second==spwId || spwId <0)  && (freqBeg < freqEnds_p[k].second) && (freqEnd > freqBegs_p[k].second)){
	isDefined=True;
	//cerr << k << " inside freqBegs " << freqBegs_p[k].second << "  " << freqBeg << endl;  
	if(freqBegs_p[k].second < freqBeg)
	  freqBegs_p[k].second=freqBeg;
	if(freqEnds_p[k].second > freqEnd)
	  freqEnds_p[k].second=freqEnd;
	if(msId < 0) key=freqBegs_p[k].first;
	//cerr << "modified " <<  freqBegs_p[k].second << "   "  <<  freqEnds_p[k].second << endl;
      }
	//cerr << "added " << k << " freqBegs " << freqBegs_p[k].second << "  " << freqEnds_p[k].second << endl;  
	frameSel.add(freqSpws_p[k].second ,  freqBegs_p[k].second, freqEnds_p[k].second);
    }
    if(!isDefined && msId >=0){
      //cerr << "undefined " << key << " freqBegs "  << freqBeg << "  " << freqEnd << endl;  
      freqBegs_p.push_back(make_pair(key, freqBeg));
      freqEnds_p.push_back(make_pair(key, freqEnd));
      freqSpws_p.push_back(make_pair(key, spwId));
      frameSel.add(spwId,  freqBeg, freqEnd);
    }
    CountedPtr<vi::FrequencySelections> copyFsels=fselections_p->clone();
    uInt nMSs=copyFsels->size() <=msId ? msId+1 : copyFsels->size();
    //cerr << "nms " << nMSs << endl;
    fselections_p=new FrequencySelections();
    for (uInt k=0;  k < nMSs ; ++k){
      if(k==uInt(key)){
	fselections_p->add(frameSel);
	//cerr <<"framesel " << frameSel.toString() << endl;
      }
      else{
	const FrequencySelectionUsingFrame& thissel= static_cast<const FrequencySelectionUsingFrame &> (copyFsels->get(k));
	//cerr <<"framesel orig " << thissel.toString() << endl;
	fselections_p->add(thissel);

      }
    }
    
 

  }

  void SynthesisImagerVi2::tuneChunk(const Int gmap){
    

    CoordinateSystem cs=itsMappers.imageStore(gmap)->getCSys();
    IPosition imshape=itsMappers.imageStore(gmap)->getShape();
    Double minFreq=SpectralImageUtil::worldFreq(cs, 0.0);
    Double maxFreq=SpectralImageUtil::worldFreq(cs,imshape(3)-1);
    if(maxFreq < minFreq){
      Double tmp=minFreq;
      minFreq=maxFreq;
      maxFreq=tmp;
    }
    Int spectralIndex=cs.findCoordinate(Coordinate::SPECTRAL);
    SpectralCoordinate spectralCoord=cs.spectralCoordinate(spectralIndex);
    MFrequency::Types intype=spectralCoord.frequencySystem(True);
    VisBufferUtil::getFreqRangeFromRange(minFreq, maxFreq,  intype, minFreq,  maxFreq, *vi_p, selFreqFrame_p);
    maxFreq+=fabs(spectralCoord.increment()(0))/2.0;
    minFreq-=fabs(spectralCoord.increment()(0))/2.0;
    if(minFreq < 0.0) minFreq=0.0;
    
    auto copyFreqBegs=freqBegs_p;
    auto copyFreqEnds=freqEnds_p;
    auto copyFreqSpws=  freqSpws_p;
    andFreqSelection(-1, -1, minFreq, maxFreq, selFreqFrame_p);
    vi_p->setFrequencySelection (*fselections_p);

    freqBegs_p=copyFreqBegs;
    freqEnds_p=copyFreqEnds;
    freqSpws_p=copyFreqSpws;
    



  }


Bool SynthesisImagerVi2::defineImage(SynthesisParamsImage& impars, 
			   const SynthesisParamsGrid& gridpars)
  {

    LogIO os( LogOrigin("SynthesisImagerVi2","defineImage",WHERE) );
    if(mss_p.nelements() ==0)
      os << "SelectData has to be run before defineImage" << LogIO::EXCEPTION;

    CoordinateSystem csys;
    CountedPtr<refim::FTMachine> ftm, iftm;


    try
      {

	os << "Define image coordinates for [" << impars.imageName << "] : " << LogIO::POST;

	csys = impars.buildCoordinateSystem( *vi_p );
	IPosition imshape = impars.shp();

	os << "Impars : start " << impars.start << LogIO::POST;
	os << "Shape : " << imshape << "Spectral : " << csys.spectralCoordinate().referenceValue() << " at " << csys.spectralCoordinate().referencePixel() << " with increment " << csys.spectralCoordinate().increment() << LogIO::POST;

	if( (itsMappers.nMappers()==0) || 
	    (impars.imsize[0]*impars.imsize[1] > itsMaxShape[0]*itsMaxShape[1]))
	  {
	    itsMaxShape=imshape;
	    itsMaxCoordSys=csys;
	  }
        itsNchan = imshape[3];
        itsCsysRec = impars.getcsys();
	/*
	os << "Define image  [" << impars.imageName << "] : nchan : " << impars.nchan 
	   //<< ", freqstart:" << impars.freqStart.getValue() << impars.freqStart.getUnit() 
	   << ", start:" << impars.start
	   <<  ", imsize:" << impars.imsize 
	   << ", cellsize: [" << impars.cellsize[0].getValue() << impars.cellsize[0].getUnit() 
	   << " , " << impars.cellsize[1].getValue() << impars.cellsize[1].getUnit() 
	   << LogIO::POST;
	*/
      }
    catch(AipsError &x)
      {
	os << "Error in building Coordinate System and Image Shape : " << x.getMesg() << LogIO::EXCEPTION;
      }

	
    try
      {
	os << "Set Gridding options for [" << impars.imageName << "] with ftmachine : " << gridpars.ftmachine << LogIO::POST;

	itsVpTable=gridpars.vpTable;
	itsMakeVP= ( gridpars.ftmachine.contains("mosaicft") ||
		             gridpars.ftmachine.contains("awprojectft") )?False:True;

	createFTMachine(ftm, iftm, gridpars.ftmachine, impars.nTaylorTerms, gridpars.mType, 
			gridpars.facets, gridpars.wprojplanes,
			gridpars.padding,gridpars.useAutoCorr,gridpars.useDoublePrec,
			gridpars.convFunc,
			gridpars.aTermOn,gridpars.psTermOn, gridpars.mTermOn,
			gridpars.wbAWP,gridpars.cfCache,gridpars.doPointing,
			gridpars.doPBCorr,gridpars.conjBeams,
			gridpars.computePAStep,gridpars.rotatePAStep,
			gridpars.interpolation, impars.freqFrameValid, 1000000000,  16, impars.stokes,
			impars.imageName);

      }
    catch(AipsError &x)
      {
	os << "Error in setting up FTMachine() : " << x.getMesg() << LogIO::EXCEPTION;
      }

    try
      {

	
		appendToMapperList(impars.imageName,  csys,  impars.shp(),
			   ftm, iftm,
			   gridpars.distance, gridpars.facets, gridpars.chanchunks,impars.overwrite,
			   gridpars.mType, gridpars.padding, impars.nTaylorTerms, impars.startModel);
	
	imageDefined_p=true;
      }
    catch(AipsError &x)
      {
	os << "Error in adding Mapper : "+x.getMesg() << LogIO::EXCEPTION;
      }

    return true;
  }


 Bool SynthesisImagerVi2::weight(const String& type, const String& rmode,
			       const Quantity& noise, const Double robust,
			       const Quantity& fieldofview,
			       const Int npixels, const Bool multiField,
			       const String& filtertype, const Quantity& filterbmaj,
			       const Quantity& filterbmin, const Quantity& filterbpa   )
  {
    LogIO os(LogOrigin("SynthesisImagerVi2", "weight()", WHERE));

       try {
    	//Int nx=itsMaxShape[0];
    	//Int ny=itsMaxShape[1];
	 Quantity cellx=Quantity(itsMaxCoordSys.increment()[0], itsMaxCoordSys.worldAxisUnits()[0]);
	 Quantity celly=Quantity(itsMaxCoordSys.increment()[1], itsMaxCoordSys.worldAxisUnits()[1]);
	 os << LogIO::NORMAL // Loglevel INFO
	    << "Set imaging weights : " ; //<< LogIO::POST;
	 
	 if (type=="natural") {
	   os << LogIO::NORMAL // Loglevel INFO
	      << "Natural weighting" << LogIO::POST;
	   imwgt_p=VisImagingWeight("natural");
	 }
      else if (type=="radial") {
	os << "Radial weighting" << LogIO::POST;
    	  imwgt_p=VisImagingWeight("radial");
      }
      else{
    	  if(!imageDefined_p)
    		  throw(AipsError("Need to define image"));
    	  Int nx=itsMaxShape[0];
    	  Int ny=itsMaxShape[1];
    	  Quantity cellx=Quantity(itsMaxCoordSys.increment()[0], itsMaxCoordSys.worldAxisUnits()[0]);
    	  Quantity celly=Quantity(itsMaxCoordSys.increment()[1], itsMaxCoordSys.worldAxisUnits()[1]);
    	  if(type=="superuniform"){
    		  if(!imageDefined_p) throw(AipsError("Please define image first"));
    		  Int actualNpix=npixels;
    		  if(actualNpix <=0)
    			  actualNpix=3;
    		  os << LogIO::NORMAL // Loglevel INFO
    				  << "SuperUniform weighting over a square cell spanning ["
    				  << -actualNpix
    				  << ", " << actualNpix << "] in the uv plane" << LogIO::POST;
    		  imwgt_p=VisImagingWeight(*vi_p, rmode, noise, robust, nx,
    				  ny, cellx, celly, actualNpix,
    				  actualNpix, multiField);
    	  }
    	  else if ((type=="robust")||(type=="uniform")||(type=="briggs")) {
   		  if(!imageDefined_p) throw(AipsError("Please define image first"));
    		  Quantity actualFieldOfView_x(fieldofview), actualFieldOfView_y(fieldofview) ;
    		  Int actualNPixels_x(npixels),actualNPixels_y(npixels) ;
    		  String wtype;
    		  if(type=="briggs") {
    			  wtype = "Briggs";
    		  }
    		  else {
    			  wtype = "Uniform";
    		  }
    		  if(actualFieldOfView_x.get().getValue()==0.0&&actualNPixels_x==0) {
    			  actualNPixels_x=nx;
    			  actualFieldOfView_x=Quantity(actualNPixels_x*cellx.get("rad").getValue(),"rad");
    			  actualNPixels_y=ny;
    			  actualFieldOfView_y=Quantity(actualNPixels_y*celly.get("rad").getValue(),"rad");
    			  os << LogIO::NORMAL // Loglevel INFO
    					  << wtype
    					  << " weighting: sidelobes will be suppressed over full image"
    					  << LogIO::POST;
    		  }
    		  else if(actualFieldOfView_x.get().getValue()>0.0&&actualNPixels_x==0) {
    			  actualNPixels_x=nx;
    			  actualNPixels_y=ny;
    			  os << LogIO::NORMAL // Loglevel INFO
    					  << wtype
    					  << " weighting: sidelobes will be suppressed over specified field of view: "
			                  << actualFieldOfView_x.get("arcsec").getValue() << " arcsec by " 
			                  << actualFieldOfView_y.get("arcsec").getValue()  << " arcsec" << LogIO::POST;
    		  }
    		  else if(actualFieldOfView_x.get().getValue()==0.0&&actualNPixels_x>0) {
    			  actualFieldOfView_x=Quantity(actualNPixels_x*cellx.get("rad").getValue(),"rad");
    			  actualFieldOfView_y=Quantity(actualNPixels_y*celly.get("rad").getValue(),"rad");
    			  os << LogIO::NORMAL // Loglevel INFO
    					  << wtype
    					  << " weighting: sidelobes will be suppressed over full image field of view: "
			                  << actualFieldOfView_x.get("arcsec").getValue() << " arcsec by " 
    					  << actualFieldOfView_y.get("arcsec").getValue() << " arcsec" << LogIO::POST;
    		  }
    		  else {
    			  os << LogIO::NORMAL // Loglevel INFO
    					  << wtype
    					  << " weighting: sidelobes will be suppressed over specified field of view: "
			                  << actualFieldOfView_x.get("arcsec").getValue() << " arcsec by " 
    					  << actualFieldOfView_y.get("arcsec").getValue() << " arcsec" << LogIO::POST;
    		  }
    		  os << LogIO::DEBUG1
		     << "Weighting used " << actualNPixels_x << " by " << actualNPixels_y << " uv pixels."
		     << LogIO::POST;
    		  Quantity actualCellSize_x(actualFieldOfView_x.get("rad").getValue()/actualNPixels_x, "rad");
    		  Quantity actualCellSize_y(actualFieldOfView_y.get("rad").getValue()/actualNPixels_y, "rad");


		  //		  cerr << "rmode " << rmode << " noise " << noise << " robust " << robust << " npixels " << actualNPixels << " cellsize " << actualCellSize << " multifield " << multiField << endl;
		  //		  Timer timer;
		  //timer.mark();
		  //Construct imwgt_p with old vi for now if old vi is in use as constructing with vi2 is slower 


		  imwgt_p=VisImagingWeight(*vi_p, wtype=="Uniform" ? "none" : rmode, noise, robust,
                                 actualNPixels_x, actualNPixels_y, actualCellSize_x,
                                 actualCellSize_y, 0, 0, multiField);

		  /*
		  if(rvi_p !=NULL){
		    imwgt_p=VisImagingWeight(*rvi_p, rmode, noise, robust,
                                 actualNPixels, actualNPixels, actualCellSize,
                                 actualCellSize, 0, 0, multiField);
		  }
		  else{
		    ////This is slower by orders of magnitude as of 2014/06/25
		    imwgt_p=VisImagingWeight(*vi_p, rmode, noise, robust,
                                 actualNPixels, actualNPixels, actualCellSize,
                                 actualCellSize, 0, 0, multiField);
		  }
		  */
		    //timer.show("After making visweight ");

    	  }
    	  else {
    		  //this->unlock();
    		  os << LogIO::SEVERE << "Unknown weighting " << type
    				  << LogIO::EXCEPTION;
    		  return false;
    	  }
      }
	 
	 //// UV-Tapering
	 //cout << "Taper type : " << filtertype << " : " << (filtertype=="gaussian") <<  endl;
	 if( filtertype == "gaussian" ) {
	   //	   os << "Setting uv-taper" << LogIO::POST;
	   imwgt_p.setFilter( filtertype,  filterbmaj, filterbmin, filterbpa );
	 }
	 vi_p->useImagingWeight(imwgt_p);
      ///////////////////////////////
	 
	 
	 ///	 return true;
	 
       }
       catch(AipsError &x)
	 {
	   throw( AipsError("Error in Weighting : "+x.getMesg()) );
	 }
       
       return true;
  }

void SynthesisImagerVi2::appendToMapperList(String imagename,  
					   CoordinateSystem& csys, 
					   IPosition imshape,
					    CountedPtr<refim::FTMachine>& ftm,
					    CountedPtr<refim::FTMachine>& iftm,
					   Quantity distance, 
					   Int facets,
					   Int chanchunks,
					   const Bool overwrite,
					   String mappertype,
					   Float padding,
					   uInt ntaylorterms,
					   Vector<String> startmodel)
    {
      LogIO log_l(LogOrigin("SynthesisImagerVi2", "appendToMapperList(ftm)"));
      //---------------------------------------------
      // Some checks..
      
      if(facets > 1 && itsMappers.nMappers() > 0)
	log_l << "Facetted image has to be the first of multifields" << LogIO::EXCEPTION;

     if(chanchunks<1)
	{
	  log_l << "Automatically calculate chanchunks";
	  log_l << " using imshape : " << imshape << LogIO::POST;

	  // Do calculation here.
	  // This runs once per image field (for multi-field imaging)
	  // This runs once per cube partition, and will see only its own partition's shape
	  chanchunks=1;

          CompositeNumber cn(uInt(imshape[0] * 2));
          // heuristic factors multiplied to imshape based on gridder
          size_t fudge_factor = 15;
          if (ftm->name()=="MosaicFTNew") {
              fudge_factor = 15;
          }
          else if (ftm->name()=="GridFT") {
              fudge_factor = 9;
          }

          size_t required_mem = fudge_factor * sizeof(Float);
          for (size_t i = 0; i < imshape.nelements(); i++) {
              // gridding pads image and increases to composite number
              if (i < 2) {
                  required_mem *= cn.nextLargerEven(Int(padding*Float(imshape[i])-0.5));
              }
              else {
                  required_mem *= imshape[i];
              }
          }

          // get number of tclean processes running on the same machine
          size_t nlocal_procs = 1;
          if (getenv("OMPI_COMM_WORLD_LOCAL_SIZE")) {
              std::stringstream ss(getenv("OMPI_COMM_WORLD_LOCAL_SIZE"));
              ss >> nlocal_procs;
          }
          // assumes all processes need the same amount of memory
          required_mem *= nlocal_procs;

          Double usr_memfrac, usr_mem;
          AipsrcValue<Double>::find(usr_memfrac, "system.resources.memfrac", 80.);
          AipsrcValue<Double>::find(usr_mem, "system.resources.memory", -1.);
          Double memory_avail;
          if (usr_mem > 0.) {
              memory_avail = usr_mem * 1024. * 1024.;
          }
          else {
              memory_avail = HostInfo::memoryTotal(false) * (usr_memfrac / 100.) * 1024.;
          }

          // compute required chanchunks to fit into the available memory
          chanchunks = (int)std::ceil((Double)required_mem / memory_avail);
          if (imshape.nelements() == 4 && imshape[3] < chanchunks) {
              chanchunks = imshape[3];
              // TODO make chanchunks a divisor of nchannels?
          }
          chanchunks = chanchunks < 1 ? 1 : chanchunks;

	  log_l << "Required memory " << required_mem / nlocal_procs / 1024. / 1024. / 1024.
                 << "\nAvailable memory " << memory_avail / 1024. / 1024 / 1024.
                 << " (rc: memory fraction " << usr_memfrac << "% rc memory " << usr_mem / 1024.
                 << ")\n" << nlocal_procs << " other processes on node\n"
                 << "Setting chanchunks to " << chanchunks << LogIO::POST;
	}

      if( imshape.nelements()==4 && imshape[3]<chanchunks )
	{
	  log_l << LogIO::WARN << "An image with " << imshape[3] << " channel(s) cannot be divided into " << chanchunks << " chunks. Please set chanchunks=1 or choose chanchunks<nchan." << LogIO::EXCEPTION;
	}

      if(chanchunks > 1 && itsMappers.nMappers() > 0)
	log_l << "Channel chunking is currently not supported with multi(outlier)-fields. Please submit a feature request if needed." << LogIO::EXCEPTION;

      if(chanchunks > 1) itsDataLoopPerMapper=true;
      
      AlwaysAssert( ( ( ! (ftm->name()=="MosaicFTNew" && mappertype=="imagemosaic") )  && 
      		      ( ! (ftm->name()=="AWProjectWBFTNew" && mappertype=="imagemosaic") )) ,
		    AipsError );
      //---------------------------------------------

      // Create the ImageStore object
      CountedPtr<SIImageStore> imstor;
      ROMSColumns msc(*(mss_p[0]));
      imstor = createIMStore(imagename, csys, imshape, overwrite,msc, mappertype, ntaylorterms, distance,facets, iftm->useWeightImage(), startmodel );

      // Create the Mappers
      if( facets<2 && chanchunks<2) // One facet. Just add the above imagestore to the mapper list.
	{
	  itsMappers.addMapper(  createSIMapper( mappertype, imstor, ftm, iftm, ntaylorterms) );
	}
      else // This field is facetted. Make a list of reference imstores, and add all to the mapper list.
	{

	  if ( facets>1 && chanchunks==1 )
	    {
	      // Make and connect the list.
	      Block<CountedPtr<SIImageStore> > imstorList = createFacetImageStoreList( imstor, facets );
	      for( uInt facet=0; facet<imstorList.nelements(); facet++)
		{
		  CountedPtr<refim::FTMachine> new_ftm, new_iftm;
		  if(facet==0){ new_ftm = ftm;  new_iftm = iftm; }
		  else{ new_ftm=ftm->cloneFTM();  new_iftm=iftm->cloneFTM(); }
		  itsMappers.addMapper(createSIMapper( mappertype, imstorList[facet], new_ftm, new_iftm, ntaylorterms));
		}
	    }// facets
	  else if ( facets==1 && chanchunks>1 )
	    {
	      // Make and connect the list.
	      Block<CountedPtr<SIImageStore> > imstorList = createChanChunkImageStoreList( imstor, chanchunks );
	      for( uInt chunk=0; chunk<imstorList.nelements(); chunk++)
		{
		  
		  CountedPtr<refim::FTMachine> new_ftm, new_iftm;
		  if(chunk==0){ 
		    new_ftm = ftm;  
		    new_iftm = iftm; }
		  else{ 
		    new_ftm=ftm->cloneFTM();  
		    new_iftm=iftm->cloneFTM(); }
		 
		  itsMappers.addMapper(createSIMapper( mappertype, imstorList[chunk], new_ftm, new_iftm, ntaylorterms));
		}
	    }// chanchunks
	  else
	    {
	      throw( AipsError("Error in requesting "+String::toString(facets)+" facets on a side with " + String::toString(chanchunks) + " channel chunks.  Support for faceting along with channel chunking is not yet available. Please submit a feature-request if you need multiple facets as well as chanchunks. ") );
	    }

	}// facets or chunks

    }

  /////////////////////////
 void SynthesisImagerVi2::runMajorCycle(const Bool dopsf, 
				      const Bool savemodel)
  {
    LogIO os( LogOrigin("SynthesisImagerVi2","runMajorCycle",WHERE) );

    //    cout << "Savemodel : " << savemodel << "   readonly : " << readOnly_p << "   usescratch : " << useScratch_p << endl;

    Bool savemodelcolumn = savemodel && !readOnly_p && useScratch_p;
    Bool savevirtualmodel = savemodel && !readOnly_p && !useScratch_p;

    if( savemodelcolumn ) os << "Saving model column" << LogIO::POST;
    if( savevirtualmodel ) os << "Saving virtual model" << LogIO::POST;

    SynthesisUtilMethods::getResource("Start Major Cycle");

    itsMappers.checkOverlappingModels("blank");

    {
      vi::VisBuffer2* vb=vi_p->getVisBuffer();
      vi_p->originChunks();
      vi_p->origin();
      Double numcoh=0;
      for (uInt k=0; k< mss_p.nelements(); ++k)
	numcoh+=Double(mss_p[k]->nrow());
      ProgressMeter pm(1.0, numcoh, 
			 dopsf?"Gridding Weights and PSF":"Major Cycle", "","","",true);
	Int cohDone=0;


    	if(!dopsf)itsMappers.initializeDegrid(*vb);
    	itsMappers.initializeGrid(*vb,dopsf);
	SynthesisUtilMethods::getResource("After initGrid for all mappers");

    	for (vi_p->originChunks(); vi_p->moreChunks();vi_p->nextChunk())
    	{

	  for (vi_p->origin(); vi_p->more(); vi_p->next())
    		{
		  //if (SynthesisUtilMethods::validate(*vb)==SynthesisUtilMethods::NOVALIDROWS) break; // No valid rows in this VB
		  //		  cerr << "nRows "<< vb->nRow() << "   " << max(vb->visCube()) <<  endl;
		  if (SynthesisUtilMethods::validate(*vb)!=SynthesisUtilMethods::NOVALIDROWS)
		    {
    			if(!dopsf) {
			  { Cube<Complex> mod(vb->nCorrelations(), vb->nChannels(), vb->nRows(), Complex(0.0));
			    vb->setVisCubeModel(mod); 
			  }
			  itsMappers.degrid(*vb, savevirtualmodel );
			  if(savemodelcolumn && writeAccess_p ){
			    //Darn not implented
			    vi_p->writeVisModel(vb->visCubeModel());
			    //static_cast<VisibilityIteratorImpl2 *> (vi_p->getImpl())->writeVisModel(vb->visCubeModel());

			    // Cube<Complex> tt=vb->visCubeModel();
			    // tt = 20.0;
			    // cout << "Vis:" << tt << endl;
			    // static_cast<VisibilityIteratorImpl2 *> (vi_p->getImpl())->writeVisModel(tt);
			  }
    			}
    			itsMappers.grid(*vb, dopsf, (refim::FTMachine::Type)datacol_p);
			
			cohDone += vb->nRows();
			pm.update(Double(cohDone));
		    }
    		}
    	}

	// cerr << "VI2 data: " << cohDone << endl;
	// exit(0);
    	//cerr << "IN SYNTHE_IMA" << endl;
    	//VisModelData::listModel(rvi_p->getMeasurementSet());
	SynthesisUtilMethods::getResource("Before finalize for all mappers");
    	if(!dopsf) itsMappers.finalizeDegrid(*vb);
    	itsMappers.finalizeGrid(*vb, dopsf);

    }

    itsMappers.checkOverlappingModels("restore");

    unlockMSs();

    SynthesisUtilMethods::getResource("End Major Cycle");

  }// end runMajorCycle

 
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /// The mapper loop is outside the data iterator loop.
  /// This is for cases where the image size is large compared to the RAM and
  /// where data I/O is the relatively minor cost.
  void SynthesisImagerVi2::runMajorCycle2(const Bool dopsf, 
				      const Bool savemodel)
  {
    LogIO os( LogOrigin("SynthesisImagerVi2","runMajorCycle2",WHERE) );

    //    cout << "Savemodel : " << savemodel << "   readonly : " << readOnly_p << "   usescratch : " << useScratch_p << endl;

    Bool savemodelcolumn = savemodel && !readOnly_p && useScratch_p;
    Bool savevirtualmodel = savemodel && !readOnly_p && !useScratch_p;

    if( savemodelcolumn ) os << "Saving model column" << LogIO::POST;
    if( savevirtualmodel ) os << "Saving virtual model" << LogIO::POST;

    itsMappers.checkOverlappingModels("blank");

    Bool resetModel=False;
    if( savemodelcolumn && writeAccess_p)
      {
	resetModel=True;
	os << "Iterating through the model column to reset it to zero" << LogIO::POST;
	vi::VisBuffer2* vb=vi_p->getVisBuffer();
    	vi_p->originChunks();
    	vi_p->origin();
	Double numcoh=0;
	for (uInt k=0; k< mss_p.nelements(); ++k)
	  numcoh+=Double(mss_p[k]->nrow());
	ProgressMeter pm(1.0, numcoh, 
			 dopsf?"Seting model column to zero":"pre-Major Cycle", "","","",True);
	Int cohDone=0;
    	for (vi_p->originChunks(); vi_p->moreChunks();vi_p->nextChunk())
	  {
	    
	    for (vi_p->origin(); vi_p->more(); vi_p->next())
	      {
		if (SynthesisUtilMethods::validate(*vb)!=SynthesisUtilMethods::NOVALIDROWS)
		  {
		    { Cube<Complex> mod(vb->nCorrelations(), vb->nChannels(), vb->nRows(), Complex(0.0));
			    vb->setVisCubeModel(mod); 
		    }
		    vi_p->writeVisModel(vb->visCubeModel());
		    
		  }
		cohDone += vb->nRows();;
		pm.update(Double(cohDone));
	      }
	  }
      }// setting model to zero

    
    for(Int gmap=0;gmap<itsMappers.nMappers();gmap++)
       {
	 os << "Running major cycle for chunk : " << gmap << LogIO::POST;

	 SynthesisUtilMethods::getResource("Start Major Cycle for mapper"+String::toString(gmap));
	 CountedPtr<vi::FrequencySelections> copyFsels=fselections_p->clone();
	 tuneChunk(gmap);
	 vi::VisBuffer2* vb=vi_p->getVisBuffer();
	 vi_p->originChunks();
	 vi_p->origin();
	 Double numcoh=0;
	 for (uInt k=0; k< mss_p.nelements(); ++k)
	   numcoh+=Double(mss_p[k]->nrow());


	 ProgressMeter pm(1.0, numcoh, 
			  dopsf?"Gridding Weights and PSF":"Major Cycle", "","","",true);
	Int cohDone=0;


	itsMappers.getFTM2(gmap, False)->reset();
	itsMappers.getFTM2(gmap, True)->reset();

    	if(!dopsf){
	  itsMappers.initializeDegrid(*vb, gmap);
		  //itsMappers.getMapper(gmap)->initializeDegrid(*vb);
	}
	itsMappers.initializeGrid(*vb,dopsf, gmap);
		//itsMappers.getMapper(gmap)->initializeGrid(*vb,dopsf);

	SynthesisUtilMethods::getResource("After initialize for mapper"+String::toString(gmap));

    	for (vi_p->originChunks(); vi_p->moreChunks();vi_p->nextChunk())
    	{

	  for (vi_p->origin(); vi_p->more(); vi_p->next())
	    {
	      //if (SynthesisUtilMethods::validate(*vb)==SynthesisUtilMethods::NOVALIDROWS) break; // No valid rows in this VB
	      //		  cerr << "nRows "<< vb->nRow() << "   " << max(vb->visCube()) <<  endl;
	      if (SynthesisUtilMethods::validate(*vb)!=SynthesisUtilMethods::NOVALIDROWS)
		{
		  if(!dopsf) {
		    if(resetModel==False) 
		      { 
			Cube<Complex> mod(vb->nCorrelations(), vb->nChannels(), vb->nRows(), Complex(0.0));
			vb->setVisCubeModel(mod); 
		      }
		    itsMappers.degrid(*vb, savevirtualmodel, gmap );
		    //itsMappers.getMapper(gmap)->degrid(*vb); //, savevirtualmodel );
		    if(savemodelcolumn && writeAccess_p ){
		      vi_p->writeVisModel(vb->visCubeModel());
		      //vi_p->writeBackChanges(vb);
		      // static_cast<VisibilityIteratorImpl2 *> (vi_p->getImpl())->writeVisModel(vb->visCubeModel());
		    }

		  }
		  itsMappers.grid(*vb, dopsf, (refim::FTMachine::Type)(datacol_p), gmap);
		  //itsMappers.getMapper(gmap)->grid(*vb, dopsf, datacol_p);
		  cohDone += vb->nRows();
		  pm.update(Double(cohDone));
		}
	    }
    	}
    	//cerr << "IN SYNTHE_IMA" << endl;
    	//VisModelData::listModel(rvi_p->getMeasurementSet());

	SynthesisUtilMethods::getResource("Before finalize for mapper"+String::toString(gmap));
	
    	if(!dopsf) 
	  {
	    itsMappers.finalizeDegrid(*vb,gmap);
	    //itsMappers.getMapper(gmap)->finalizeDegrid();
	  }
	itsMappers.finalizeGrid(*vb, dopsf,gmap);
    	//itsMappers.getMapper(gmap)->finalizeGrid(*vb, dopsf);
	
	//	itsMappers.getMapper(gmap)->releaseImageLocks();
	itsMappers.getMapper(gmap)->imageStore()->releaseComplexGrids();        
	
	SynthesisUtilMethods::getResource("End Major Cycle for mapper"+String::toString(gmap));
	fselections_p=copyFsels;
       }// end of mapper loop
    vi_p->setFrequencySelection(*fselections_p);

    itsMappers.checkOverlappingModels("restore");

    unlockMSs();

    SynthesisUtilMethods::getResource("End Major Cycle");

  }// end runMajorCycle2


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void SynthesisImagerVi2::predictModel(){
    LogIO os( LogOrigin("SynthesisImagerVi2","predictModel ",WHERE) );

    os << "---------------------------------------------------- Predict Model ---------------------------------------------" << LogIO::POST;
    
    Bool savemodelcolumn = !readOnly_p && useScratch_p;
    Bool savevirtualmodel = !readOnly_p && !useScratch_p;

    if( savemodelcolumn ) os << "Saving model column" << LogIO::POST;
    if( savevirtualmodel ) os << "Saving virtual model" << LogIO::POST;

    itsMappers.checkOverlappingModels("blank");


    {
      vi::VisBuffer2* vb = vi_p->getVisBuffer();;
      vi_p->originChunks();
      vi_p->origin();
      Double numberCoh=0;
      for (uInt k=0; k< mss_p.nelements(); ++k)
	numberCoh+=Double(mss_p[k]->nrow());

      ProgressMeter pm(1.0, numberCoh, "Predict Model", "","","",true);
      Int cohDone=0;

      itsMappers.initializeDegrid(*vb);
      for (vi_p->originChunks(); vi_p->moreChunks();vi_p->nextChunk())
	{
	  
	  for (vi_p->origin(); vi_p->more(); vi_p->next())
	    {
	      //if (SynthesisUtilMethods::validate(*vb)==SynthesisUtilMethods::NOVALIDROWS) break; //No valid rows in this MS
	      //if !usescratch ...just save
	      vb->setVisCubeModel(Complex(0.0, 0.0));
	      itsMappers.degrid(*vb, savevirtualmodel);
	      if(savemodelcolumn && writeAccess_p )
		vb->setVisCubeModel(vb->visCubeModel());

	      //	      cout << "nRows "<< vb->nRow() << "   " << max(vb->modelVisCube()) <<  endl;
	      cohDone += vb->nRows();
	      pm.update(Double(cohDone));

	    }
	}
      itsMappers.finalizeDegrid(*vb);
    }

    itsMappers.checkOverlappingModels("restore");
    unlockMSs();
   
  }// end of predictModel
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


CountedPtr<SIMapper> SynthesisImagerVi2::createSIMapper(String mappertype,  
							   CountedPtr<SIImageStore> imagestore,
							CountedPtr<refim::FTMachine> ftmachine,
							CountedPtr<refim::FTMachine> iftmachine,
						       uInt /*ntaylorterms*/)
  {
    LogIO os( LogOrigin("SynthesisImagerVi2","createSIMapper",WHERE) );
    
    CountedPtr<SIMapper> localMapper;

    try
      {
	
	if( mappertype == "default" || mappertype == "multiterm" )
	  {
	    localMapper = new SIMapper( imagestore, ftmachine, iftmachine );
	  }
	else if( mappertype == "imagemosaic") // || mappertype == "mtimagemosaic" )
	  {
	    localMapper = new SIMapperImageMosaic( imagestore, ftmachine, iftmachine );
	  }
	else
	  {
	    throw(AipsError("Unknown mapper type : " + mappertype));
	  }

      }
    catch(AipsError &x) {
	throw(AipsError("Error in createSIMapper : " + x.getMesg() ) );
      }
    return localMapper;
  }
  

void SynthesisImagerVi2::unlockMSs()
  {
    LogIO os( LogOrigin("SynthesisImagerVi2","unlockMSs",WHERE) );
    for(uInt i=0;i<mss_p.nelements();i++)
      { 
	os << LogIO::NORMAL2 << "Unlocking : " << (mss_p[i])->tableName() << LogIO::POST;
	MeasurementSet *ms_l = 	const_cast<MeasurementSet* >(mss_p[i]);
	ms_l->unlock(); 
	ms_l->antenna().unlock();
	ms_l->dataDescription().unlock();
	ms_l->feed().unlock();
	ms_l->field().unlock();
	ms_l->observation().unlock();
	ms_l->polarization().unlock();
	ms_l->processor().unlock();
	ms_l->spectralWindow().unlock();
	ms_l->state().unlock();
	//
	// Unlock the optional sub-tables as well, if they are present
	//
	if(!(ms_l->source().isNull()))     ms_l->source().unlock();
	if(!(ms_l->doppler().isNull()))    ms_l->doppler().unlock();
	if(!(ms_l->flagCmd().isNull()))    ms_l->flagCmd().unlock();
	if(!(ms_l->freqOffset().isNull())) ms_l->freqOffset().unlock();
	if(!(ms_l->history().isNull()))    ms_l->history().unlock();
	if(!(ms_l->pointing().isNull()))   ms_l->pointing().unlock();
	if(!(ms_l->sysCal().isNull()))     ms_l->sysCal().unlock();
	if(!(ms_l->weather().isNull()))    ms_l->weather().unlock();
      }
  }
  void SynthesisImagerVi2::createFTMachine(CountedPtr<refim::FTMachine>& theFT, 
					   CountedPtr<refim::FTMachine>& theIFT, 
					   const String& ftname,
					   const uInt nTaylorTerms,
					   const String mType,
					   const Int facets,            //=1
					   //------------------------------
					   const Int wprojplane,        //=1,
					   const Float padding,         //=1.0,
					   const Bool useAutocorr,      //=false,
					   const Bool useDoublePrec,    //=true,
					   const String gridFunction,   //=String("SF"),
					//------------------------------
					   const Bool aTermOn,          //= true,
					   const Bool psTermOn,         //= true,
					   const Bool mTermOn,          //= false,
					const Bool wbAWP,            //= true,
					   const String cfCache,        //= "",
					   const Bool doPointing,       //= false,
					   const Bool doPBCorr,         //= true,
					   const Bool conjBeams,        //= true,
					const Float computePAStep,         //=360.0
					   const Float rotatePAStep,          //=5.0
					   const String interpolation,  //="linear"
					   const Bool freqFrameValid, //=true
					   const Int cache,             //=1000000000,
					   const Int tile,               //=16
					   const String stokes, //=I
					   const String imageNamePrefix
					   )

  {
    LogIO os( LogOrigin("SynthesisImagerVi2","createFTMachine",WHERE));

    if(ftname=="gridft"){
      if(facets >1){
	theFT=new refim::GridFT(cache, tile, gridFunction, mLocation_p, phaseCenter_p, padding, useAutocorr, useDoublePrec);
	theIFT=new refim::GridFT(cache, tile, gridFunction, mLocation_p, phaseCenter_p, padding, useAutocorr, useDoublePrec);

      }
      else{
	theFT=new refim::GridFT(cache, tile, gridFunction, mLocation_p, padding, useAutocorr, useDoublePrec);
	theIFT=new refim::GridFT(cache, tile, gridFunction, mLocation_p, padding, useAutocorr, useDoublePrec);
      }
    }
    else if(ftname== "wprojectft"){
     Double maxW=-1.0;
     Double minW=-1.0;
     Double rmsW=-1.0;
     if(wprojplane <1)
       casa::refim::WProjectFT::wStat(*vi_p, minW, maxW, rmsW);
    if(facets >1){
      theFT=new refim::WProjectFT(wprojplane,  phaseCenter_p, mLocation_p,
			   cache/2, tile, useAutocorr, padding, useDoublePrec, minW, maxW, rmsW);
      theIFT=new refim::WProjectFT(wprojplane,  phaseCenter_p, mLocation_p,
			    cache/2, tile, useAutocorr, padding, useDoublePrec, minW, maxW, rmsW);
    }
    else{
      theFT=new refim::WProjectFT(wprojplane,  mLocation_p,
			   cache/2, tile, useAutocorr, padding, useDoublePrec, minW, maxW, rmsW);
      theIFT=new refim::WProjectFT(wprojplane,  mLocation_p,
			    cache/2, tile, useAutocorr, padding, useDoublePrec, minW, maxW, rmsW);
    }
    CountedPtr<refim::WPConvFunc> sharedconvFunc=static_cast<refim::WProjectFT &>(*theFT).getConvFunc();
      //static_cast<WProjectFT &>(*theFT).setConvFunc(sharedconvFunc);
    static_cast<refim::WProjectFT &>(*theIFT).setConvFunc(sharedconvFunc);
    }
    else if ((ftname == "awprojectft") || (ftname== "mawprojectft") || (ftname == "protoft")) {
      createAWPFTMachine(theFT, theIFT, ftname, facets, wprojplane, 
			 padding, useAutocorr, useDoublePrec, gridFunction,
			 aTermOn, psTermOn, mTermOn, wbAWP, cfCache, 
			 doPointing, doPBCorr, conjBeams, computePAStep,
			 rotatePAStep, cache,tile,imageNamePrefix);
    }
    else if ( ftname == "mosaic" || ftname== "mosft" || ftname == "mosaicft" || ftname== "MosaicFT"){

      createMosFTMachine(theFT, theIFT, padding, useAutocorr, useDoublePrec, rotatePAStep, stokes);
    }
    else
      {
	throw( AipsError( "Invalid FTMachine name : " + ftname ) );
      }
    /* else if(ftname== "MosaicFT"){

       }*/



    ///////// Now, clone and pack the chosen FT into a MultiTermFT if needed.
    if( mType=="multiterm" )
      {
	AlwaysAssert( nTaylorTerms>=1 , AipsError );

	CountedPtr<refim::FTMachine> theMTFT = new refim::MultiTermFTNew( theFT , nTaylorTerms, true/*forward*/ );
	CountedPtr<refim::FTMachine> theMTIFT = new refim::MultiTermFTNew( theIFT , nTaylorTerms, false/*forward*/ );

	theFT = theMTFT;
	theIFT = theMTIFT;
      }




    ////// Now, set the SkyJones if needed, and if not internally generated.
    if( mType=="imagemosaic" && 
	(ftname != "awprojectft" && ftname != "mawprojectft" && ftname != "proroft") )
      {
	CountedPtr<refim::SkyJones> vp;
	ROMSColumns msc(*(mss_p[0]));
	Quantity parang(0.0,"deg");
	Quantity skyposthreshold(0.0,"deg");
	vp = new refim::VPSkyJones(msc, true,  parang, BeamSquint::NONE,skyposthreshold);

	Vector<CountedPtr<refim::SkyJones> > skyJonesList(1);
	skyJonesList(0) = vp;
	theFT->setSkyJones(  skyJonesList );
	theIFT->setSkyJones(  skyJonesList );

      }

    //// For mode=cubedata, set the freq frame to invalid..
    // get this info from buildCoordSystem
    //theFT->setSpw( tspws, false );
    //theIFT->setSpw( tspws, false );
    theFT->setFrameValidity( freqFrameValid );
    theIFT->setFrameValidity( freqFrameValid );

    //// Set interpolation mode
    theFT->setFreqInterpolation( interpolation );
    theIFT->setFreqInterpolation( interpolation );
    /* vi_p has chanselection now
    //channel selections from spw param
    theFT->setSpwChanSelection(chanSel_p);
    theIFT->setSpwChanSelection(chanSel_p);
    */
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void SynthesisImagerVi2::createAWPFTMachine(CountedPtr<refim::FTMachine>& theFT, CountedPtr<refim::FTMachine>& theIFT, 
					   const String&,// ftmName,
					   const Int,// facets,            //=1
					   //------------------------------
					   const Int wprojPlane,        //=1,
					   const Float,// padding,         //=1.0,
					   const Bool,// useAutocorr,      //=false,
					   const Bool useDoublePrec,    //=true,
					   const String,// gridFunction,   //=String("SF"),
					   //------------------------------
					   const Bool aTermOn,          //= true,
					   const Bool psTermOn,         //= true,
					   const Bool mTermOn,          //= false,
					   const Bool wbAWP,            //= true,
					   const String cfCache,        //= "",
					   const Bool doPointing,       //= false,
					   const Bool doPBCorr,         //= true,
					   const Bool conjBeams,        //= true,
					   const Float computePAStep,   //=360.0
					   const Float rotatePAStep,    //=5.0
					   const Int cache,             //=1000000000,
					   const Int tile,               //=16
					   const String imageNamePrefix
					)

  {
    LogIO os( LogOrigin("SynthesisImagerVi2","createAWPFTMachine",WHERE));

    if (wprojPlane<=1)
      {
	os << LogIO::NORMAL
	   << "You are using wprojplanes=1. Doing co-planar imaging (no w-projection needed)" 
	   << LogIO::POST;
	os << LogIO::NORMAL << "Performing WBA-Projection" << LogIO::POST; // Loglevel PROGRESS
      }
    // if((wprojPlane>1)&&(wprojPlane<64)) 
    //   {
    // 	os << LogIO::WARN
    // 	   << "No. of w-planes set too low for W projection - recommend at least 128"
    // 	   << LogIO::POST;
    // 	os << LogIO::NORMAL << "Performing WBAW-Projection" << LogIO::POST; // Loglevel PROGRESS
    //   }

    // CountedPtr<ATerm> apertureFunction = createTelescopeATerm(mss4vi_p[0], aTermOn);
    // CountedPtr<PSTerm> psTerm = new PSTerm();
    // CountedPtr<WTerm> wTerm = new WTerm();
    
    // //
    // // Selectively switch off CFTerms.
    // //
    // if (aTermOn == false) {apertureFunction->setOpCode(CFTerms::NOOP);}
    // if (psTermOn == false) psTerm->setOpCode(CFTerms::NOOP);

    // //
    // // Construct the CF object with appropriate CFTerms.
    // //
    // CountedPtr<ConvolutionFunction> tt;
    // tt = AWProjectFT::makeCFObject(aTermOn, psTermOn, true, mTermOn, wbAWP);
    // CountedPtr<ConvolutionFunction> awConvFunc;
    // //    awConvFunc = new AWConvFunc(apertureFunction,psTerm,wTerm, !wbAWP);
    // if ((ftmName=="mawprojectft") || (mTermOn))
    //   awConvFunc = new AWConvFuncEPJones(apertureFunction,psTerm,wTerm,wbAWP);
    // else
    //   awConvFunc = new AWConvFunc(apertureFunction,psTerm,wTerm,wbAWP);

    ROMSObservationColumns msoc((mss_p[0])->observation());
    String telescopeName=msoc.telescopeName()(0);
    CountedPtr<refim::ConvolutionFunction> awConvFunc = refim::AWProjectFT::makeCFObject(telescopeName, 
									   aTermOn,
									   psTermOn, (wprojPlane > 1),
									   mTermOn, wbAWP, conjBeams);
    //
    // Construct the appropriate re-sampler.
    //
    CountedPtr<refim::VisibilityResamplerBase> visResampler;
    //    if (ftmName=="protoft") visResampler = new ProtoVR();
    //elsef
    visResampler = new refim::AWVisResampler();
    //    CountedPtr<VisibilityResamplerBase> visResampler = new VisibilityResampler();

    //
    // Construct and initialize the CF cache object.
    //


    // CountedPtr<CFCache> cfCacheObj = new CFCache();
    // cfCacheObj->setCacheDir(cfCache.data());
    // //    cerr << "Setting wtImagePrefix to " << imageNamePrefix.c_str() << endl;
    // cfCacheObj->setWtImagePrefix(imageNamePrefix.c_str());
    // cfCacheObj->initCache2();

    CountedPtr<refim::CFCache> cfCacheObj;
      

    //
    // Finally construct the FTMachine with the CFCache, ConvFunc and
    // Re-sampler objects.  
    //
    Float pbLimit_l=1e-3;
    theFT = new refim::AWProjectWBFTNew(wprojPlane, cache/2, 
			      cfCacheObj, awConvFunc, 
			      visResampler,
			      /*true */doPointing, doPBCorr, 
			      tile, computePAStep, pbLimit_l, true,conjBeams,
			      useDoublePrec);

    cfCacheObj = new refim::CFCache();
    cfCacheObj->setCacheDir(cfCache.data());
    //    cerr << "Setting wtImagePrefix to " << imageNamePrefix.c_str() << endl;
    cfCacheObj->setWtImagePrefix(imageNamePrefix.c_str());
    cfCacheObj->initCache2();

    theFT->setCFCache(cfCacheObj);
    

    Quantity rotateOTF(rotatePAStep,"deg");
    static_cast<refim::AWProjectWBFTNew &>(*theFT).setObservatoryLocation(mLocation_p);
    static_cast<refim::AWProjectWBFTNew &>(*theFT).setPAIncrement(Quantity(computePAStep,"deg"),rotateOTF);

    // theIFT = new AWProjectWBFT(wprojPlane, cache/2, 
    // 			       cfCacheObj, awConvFunc, 
    // 			       visResampler,
    // 			       /*true */doPointing, doPBCorr, 
    // 			       tile, computePAStep, pbLimit_l, true,conjBeams,
    // 			       useDoublePrec);

    // static_cast<AWProjectWBFT &>(*theIFT).setObservatoryLocation(mLocation_p);
    // static_cast<AWProjectWBFT &>(*theIFT).setPAIncrement(Quantity(computePAStep,"deg"),rotateOTF);

    theIFT = new refim::AWProjectWBFTNew(static_cast<refim::AWProjectWBFTNew &>(*theFT));

    os << "Sending frequency selection information " <<  mssFreqSel_p  <<  " to AWP FTM." << LogIO::POST;
    theFT->setSpwFreqSelection( mssFreqSel_p );
    theIFT->setSpwFreqSelection( mssFreqSel_p );
    

  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void SynthesisImagerVi2:: createMosFTMachine(CountedPtr<refim::FTMachine>& theFT,CountedPtr<refim::FTMachine>&  theIFT, const Float /*padding*/, const Bool useAutoCorr, const Bool useDoublePrec, const Float rotatePAStep, const String stokes){
    
    LogIO os(LogOrigin("SynthesisImagerVi2", "createMosFTMachine",WHERE));
   
    ROMSColumns msc(vi_p->ms());
    String telescop=msc.observation().telescopeName()(0);
    Bool multiTel=False;
    Int msid=0;
     for(vi_p->originChunks(); vi_p->moreChunks(); vi_p->nextChunk()){
       if(((vi_p->getVisBuffer())->msId() != msid) && telescop !=  ROMSColumns(vi_p->ms()).observation().telescopeName()(0)){
	 msid=(vi_p->getVisBuffer())->msId();
	 multiTel=True;
       }
     }
    vi_p->originChunks();
  
  

    PBMath::CommonPB kpb;
    Record rec;
    getVPRecord( rec, kpb, telescop );
   

    if(rec.empty()){os << LogIO::SEVERE << "Cannot proceed with mosaicft gridder without a valid PB model" << LogIO::POST; }
    
    /*
   VPManager *vpman=VPManager::Instance();
    PBMath::CommonPB kpb;
    PBMath::enumerateCommonPB(telescop, kpb);
    Record rec;
    vpman->getvp(rec, telescop);
    */

   refim::VPSkyJones* vps=NULL;
    if(rec.asString("name")=="COMMONPB" && kpb !=PBMath::UNKNOWN ){
      vps= new refim::VPSkyJones(msc, true, Quantity(rotatePAStep, "deg"), BeamSquint::GOFIGURE, Quantity(360.0, "deg"));
      /////Don't know which parameter has pb threshold cutoff that the user want 
      ////leaving at default
      ////vps.setThreshold(minPB);
      
    }
    else{
      PBMath myPB(rec);
      String whichPBMath;
      PBMathInterface::namePBClass(myPB.whichPBClass(), whichPBMath);
      os  << "Using the PB defined by " << whichPBMath << " for beam calculation for telescope " << telescop << LogIO::POST;
      vps= new refim::VPSkyJones(telescop, myPB, Quantity(rotatePAStep, "deg"), BeamSquint::GOFIGURE, Quantity(360.0, "deg"));
      kpb=PBMath::DEFAULT;
    }
   
    
    theFT = new refim::MosaicFTNew(vps, mLocation_p, stokes, 1000000000, 16, useAutoCorr, 
		      useDoublePrec);
    PBMathInterface::PBClass pbtype=((kpb==PBMath::EVLA) || multiTel)? PBMathInterface::COMMONPB: PBMathInterface::AIRY;
    if(rec.asString("name")=="IMAGE")
       pbtype=PBMathInterface::IMAGE;
    ///Use Heterogenous array mode for the following
    if((kpb == PBMath::UNKNOWN) || (kpb==PBMath::OVRO) || (kpb==PBMath::ACA)
       || (kpb==PBMath::ALMA) || (kpb==PBMath::EVLA) || multiTel){
      CountedPtr<refim::SimplePBConvFunc> mospb=new refim::HetArrayConvFunc(pbtype, "");
      static_cast<refim::MosaicFTNew &>(*theFT).setConvFunc(mospb);
    }
    ///////////////////make sure both FTMachine share the same conv functions.
    theIFT= new refim::MosaicFTNew(static_cast<refim::MosaicFTNew &>(*theFT));

    
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //// Get/Set Weight Grid.... write to disk and read

  /// todo : do for full mapper list, and taylor terms.
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// todo : do for full mapper list, and taylor terms.
  
  Bool SynthesisImagerVi2::setWeightDensity( )
  {
    LogIO os(LogOrigin("SynthesisImagerVi2", "setWeightDensity()", WHERE));
    try
      {
	Block<Matrix<Float> > densitymatrices(itsMappers.nMappers());
	for (uInt fid=0;fid<densitymatrices.nelements();fid++)
	  {
	    Array<Float> arr;
	    itsMappers.imageStore(fid)->gridwt(0)->get(arr,true);
	    densitymatrices[fid].reference( arr );
	    //cout << "Density shape (set) for f " << fid << " : " << arr.shape() << " : " << densitymatrices[fid].shape() << endl;
	  }


	imwgt_p.setWeightDensity( densitymatrices );
	vi_p->useImagingWeight(imwgt_p);
	itsMappers.releaseImageLocks();

      }
    catch (AipsError &x)
      {
	throw(AipsError("In setWeightDensity : "+x.getMesg()));
      }
    return true;
  }
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void SynthesisImagerVi2::createVisSet(const Bool /*writeAccess*/)
  {
    LogIO os( LogOrigin("SynthesisImagerVi2","createVisSet",WHERE) );
    //cerr << "mss_p num" << mss_p.nelements() <<  " sel  " << fselections_p->size() << endl;
    if(mss_p.nelements() > uInt(fselections_p->size()) && (fselections_p->size() !=0)){
      throw(AipsError("Discrepancy between Number of MSs and Frequency selections"));
    }
    vi_p=new vi::VisibilityIterator2(mss_p, vi::SortColumns(), true); //writeAccess);

    if(fselections_p->size() !=0){
      CountedPtr<vi::FrequencySelections> tmpfselections=new FrequencySelections();
      //Temporary fix till we get rid of old vi and we can get rid of tuneSelect
      if(uInt(fselections_p->size()) > mss_p.nelements()){
	for(uInt k=0 ; k <  mss_p.nelements(); ++k){
	  tmpfselections->add(fselections_p->get(k));
	}
      }
      else{
	tmpfselections=fselections_p;
      }
      ////end of fix for tuneSelectdata 
      vi_p->setFrequencySelection (*tmpfselections);

    }
    //
    vi_p->originChunks();
    vi_p->origin();
  }// end of createVisSet

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Method to run the AWProjectFT machine to seperate the CFCache
  // construction from imaging.  This is done by splitting the
  // operation in two steps: (1) run the FTM in "dry" mode to create a
  // "blank" CFCache, and (2) re-load the "blank" CFCache and "fill"
  // it.
  //
  // If someone can get me (SB) out of the horrible statc_casts in the
  // code below, I will be most grateful (we are out of it! :-)).
  //
  void SynthesisImagerVi2::dryGridding(const Vector<String>& cfList)
  {
    LogIO os( LogOrigin("SynthesisImagerVi2","dryGridding",WHERE) );

    Int cohDone=0, whichFTM=0;
    (void)cfList;
    // If not an AWProject-class FTM, make this call a NoOp.  Might be
    // useful to extend it to other projection FTMs -- but later.
    String ftmName = ((*(itsMappers.getFTM2(whichFTM)))).name();

    if (!((itsMappers.getFTM2(whichFTM,true))->isUsingCFCache())) return;

    os << "---------------------------------------------------- Dry Gridding ---------------------------------------------" << LogIO::POST;

    //
    // Go through the entire MS in "dry" mode to set up a "blank"
    // CFCache.  This is done by setting the AWPWBFT in dryrun mode
    // and gridding.  The process of gridding emits CFCache, which
    // will be "blank" in a dry run.
    {
      vi::VisBuffer2* vb=vi_p->getVisBuffer();
      vi_p->originChunks();
      vi_p->origin();
      Double numberCoh=0;
      for (uInt k=0; k< mss_p.nelements(); ++k)
	numberCoh+=Double(mss_p[k]->nrow());

      ProgressMeter pm(1.0, numberCoh, "dryGridding", "","","",true);

      itsMappers.initializeGrid(*vb);
    
      // Set the gridder (iFTM) to run in dry-gridding mode
      (itsMappers.getFTM2(whichFTM,true))->setDryRun(true);

      Bool aTermIsOff=False;
      {
	CountedPtr<refim::FTMachine> ftm=itsMappers.getFTM2(whichFTM,True);
	CountedPtr<refim::ConvolutionFunction> cf=ftm->getAWConvFunc();
	aTermIsOff = cf->getTerm("ATerm")->isNoOp();
      }

      os << "Making a \"blank\" CFCache"
	 << (aTermIsOff?" (without the A-Term)":"")
	 << LogIO::WARN << LogIO::POST;

      // Step through the MS.  This triggers the logic in the Gridder
      // to determine all the CFs that will be required.  These empty
      // CFs are written to the CFCache which can then be filled via
      // a call to fillCFCache().
      for (vi_p->originChunks(); vi_p->moreChunks();vi_p->nextChunk())
	{
	  for (vi_p->origin(); vi_p->more(); vi_p->next())
	    {
	      if (SynthesisUtilMethods::validate(*vb)!=SynthesisUtilMethods::NOVALIDROWS) 
		{
		  itsMappers.grid(*vb, true, refim::FTMachine::OBSERVED, whichFTM);
		  cohDone += vb->nRows();
		  pm.update(Double(cohDone));
		  // If there is no term that depends on time, don't iterate over the entire data base
		  if (aTermIsOff) break;
		}
	    }
	}
    }
    if (cohDone == 0) os << "No valid rows found in dryGridding." << LogIO::EXCEPTION << LogIO::POST;
    // Unset the dry-gridding mode.
    (itsMappers.getFTM2(whichFTM,true))->setDryRun(false);

    //itsMappers.checkOverlappingModels("restore");
    unlockMSs();
    //fillCFCache(cfList);
  }
  //
  // Re-load the CFCache from the disk using the supplied list of CFs
  // (as cfList).  Then extract the ConvFunc object (which was setup
  // in the FTM) and call it's makeConvFunction2() to fill the CF.
  // Finally, unset the dry-run mode of the FTM.
  //
  void SynthesisImagerVi2::fillCFCache(const Vector<String>& cfList,
				       const String& ftmName,
				       const String& cfcPath,
				       const Bool& psTermOn,
				       const Bool& aTermOn,
				       const Bool& conjBeams)
    {
      LogIO os( LogOrigin("SynthesisImagerVi2","fillCFCache",WHERE) );
      // If not an AWProject-class FTM, make this call a NoOp.  Might be
      // useful to extend it to other projection FTMs -- but later.
      // String ftmName = ((*(itsMappers.getFTM(whichFTM)))).name();

      if (!ftmName.contains("awproject") and
	  !ftmName.contains("multitermftnew")) return;
      //if (!ftmName.contains("awproject")) return;
      
      os << "---------------------------------------------------- fillCFCache ---------------------------------------------" << LogIO::POST;

      //String cfcPath = itsMappers.getFTM(whichFTM)->getCacheDir();
      //String imageNamePrefix=itsMappers.getFTM(whichFTM)->getCFCache()->getWtImagePrefix();

      //cerr << "Path = " << path << endl;

      // CountedPtr<AWProjectWBFTNew> tmpFT = new AWProjectWBFTNew(static_cast<AWProjectWBFTNew &> (*(itsMappers.getFTM(whichFTM))));


      Float dPA=360.0,selectedPA=2*360.0;
      if (cfList.nelements() > 0)
      {
	CountedPtr<refim::CFCache> cfCacheObj = new refim::CFCache();
	  //Vector<String> wtCFList; wtCFList.resize(cfList.nelements());
	  //for (Int i=0; i<wtCFList.nelements(); i++) wtCFList[i] = "WT"+cfList[i];
	  //Directory dir(path);
	  Vector<String> cfList_p=cfList;//dir.find(Regex(Regex::fromPattern("CFS*")));
	  Vector<String> wtCFList_p;
	  wtCFList_p.resize(cfList_p.nelements());
	  for (Int i=0; i<(Int)wtCFList_p.nelements(); i++) wtCFList_p[i]="WT"+cfList_p[i];

	  //cerr << cfList_p << endl;
      	  cfCacheObj->setCacheDir(cfcPath.data());

	  os << "Re-loading the \"blank\" CFCache for filling" << LogIO::WARN << LogIO::POST;

      	  cfCacheObj->initCacheFromList2(cfcPath, cfList_p, wtCFList_p,
      					 selectedPA, dPA,1);
	  // tmpFT->setCFCache(cfCacheObj);
	  Vector<Double> uvScale, uvOffset;
	  Matrix<Double> vbFreqSelection;
	  CountedPtr<refim::CFStore2> cfs2 = CountedPtr<refim::CFStore2>(&cfCacheObj->memCache2_p[0],false);//new CFStore2;
	  CountedPtr<refim::CFStore2> cfwts2 =  CountedPtr<refim::CFStore2>(&cfCacheObj->memCacheWt2_p[0],false);//new CFStore2;

	  //
	  // Get whichFTM from itsMappers (SIMapperCollection) and
	  // cast it as AWProjectWBFTNew.  Then get the ConvFunc from
	  // the FTM and cast it as AWConvFunc.  Finally call
	  // AWConvFunc::makeConvFunction2().
	  //
	  // (static_cast<AWConvFunc &> 
	  //  (*(static_cast<AWProjectWBFTNew &> (*(itsMappers.getFTM(whichFTM)))).getAWConvFunc())
	  //  ).makeConvFunction2(String(path), uvScale, uvOffset, vbFreqSelection,
	  // 		       *cfs2, *cfwts2);

	  // This is a global methond in AWConvFunc.  Does not require
	  // FTM to be constructed (which can be expensive in terms of
	  // memory footprint).
	  refim::AWConvFunc::makeConvFunction2(String(cfcPath), uvScale, uvOffset, vbFreqSelection,
					       *cfs2, *cfwts2, psTermOn, aTermOn, conjBeams);
      	}
      //cerr << "Mem used = " << itsMappers.getFTM(whichFTM)->getCFCache()->memCache2_p[0].memUsage() << endl;
      //(static_cast<AWProjectWBFTNew &> (*(itsMappers.getFTM(whichFTM)))).getCFCache()->initCache2();
    }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void SynthesisImagerVi2::reloadCFCache()
  {
      LogIO os( LogOrigin("SynthesisImagerVi2","reloadCFCache",WHERE) );
      Int whichFTM=0;
      String ftmName = ((*(itsMappers.getFTM2(whichFTM)))).name();
      if (!ftmName.contains("AWProject")) return;

      os << "-------------------------------------------- reloadCFCache ---------------------------------------------" << LogIO::POST;
      String path = itsMappers.getFTM2(whichFTM)->getCacheDir();
      String imageNamePrefix=itsMappers.getFTM2(whichFTM)->getCFCache()->getWtImagePrefix();

      CountedPtr<refim::CFCache> cfCacheObj = new refim::CFCache();
      cfCacheObj->setCacheDir(path.data());
      cfCacheObj->setWtImagePrefix(imageNamePrefix.c_str());
      cfCacheObj->initCache2();
      
      // This assumes the itsMappers is always SIMapperCollection.
      for (whichFTM = 0; whichFTM < itsMappers.nMappers(); whichFTM++)
	{
	  (static_cast<refim::AWProjectWBFTNew &> (*(itsMappers.getFTM2(whichFTM)))).setCFCache(cfCacheObj,true); // Setup iFTM
	  (static_cast<refim::AWProjectWBFTNew &> (*(itsMappers.getFTM2(whichFTM,false)))).setCFCache(cfCacheObj,true); // Set FTM
	}
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Bool SynthesisImagerVi2::makePB()
  {
      LogIO os( LogOrigin("SynthesisImagerVi2","makePB",WHERE) );

      if( itsMakeVP==False )
	{
	  os << LogIO::NORMAL1 << "Not making .pb by direct evaluation. The gridder will make a .weight and a .pb will be computed from it." << LogIO::POST;
	  // Check that the .weight exists.. ?

	  return False;
	}
      else
	{
	  Bool doDefaultVP = itsVpTable.length()>0 ? False : True;

	  CoordinateSystem coordsys=itsMappers.imageStore(0)->getCSys();
	  String telescope=coordsys.obsInfo().telescope();
	  
	  if (doDefaultVP) {
	    
	    ROMSAntennaColumns ac(mss_p[0]->antenna());
	    Double dishDiam=ac.dishDiameter()(0);
	    if(!allEQ(ac.dishDiameter().getColumn(), dishDiam))
	      os << LogIO::WARN
		 << "The MS has multiple antenna diameters ..PB could be wrong "
		 << LogIO::POST;
	    return makePBImage( telescope, False, dishDiam);
	  }
	  else{
	    return makePBImage(telescope );	
	  }
	  
	}
 
      return False;
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Bool SynthesisImagerVi2::makePrimaryBeam(PBMath& pbMath)
  {
    LogIO os( LogOrigin("SynthesisImagerVi2","makePrimaryBeam",WHERE) );

    os << "vi2 : Evaluating Primary Beam model onto image grid(s)" << LogIO::POST;

    itsMappers.initPB();

    vi::VisBuffer2* vb = vi_p->getVisBuffer();
    vi_p->originChunks();
    vi_p->origin();
    Int fieldCounter=0;
    Vector<Int> fieldsDone;

    for(vi_p->originChunks(); vi_p->moreChunks(); vi_p->nextChunk())
      {
	for (vi_p->origin(); vi_p->more(); vi_p->next())
	  {
	    Bool fieldDone=False;
	    for (uInt k=0;  k < fieldsDone.nelements(); ++k)
	      fieldDone=fieldDone || (vb->fieldId()(0)==fieldsDone(k));
	    if(!fieldDone){
	      ++fieldCounter;
	      fieldsDone.resize(fieldCounter, True);
	      fieldsDone(fieldCounter-1)=vb->fieldId()(0);
	      
	      itsMappers.addPB(*vb,pbMath);
	      
	    }
	  }
      }
    itsMappers.releaseImageLocks();
    unlockMSs();

    return True;
  }// end makePB




} //# NAMESPACE CASA - END

