// -*- C++ -*-
//# InteractiveMasking.cc: Implementation of the InteractiveMasking class
//# Copyright (C) 1997,1998,1999,2000,2001,2002,2003
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id$
//
#include <synthesis/ImagerObjects/InteractiveMasking.h>
#include <images/Images/PagedImage.h>
#include <tables/Tables/Table.h>

#include <images/Images/PagedImage.h>
#include <imageanalysis/ImageAnalysis/ImageStatsCalculator.h>
#include <casa/Logging/LogIO.h>

using namespace casacore;
namespace casa{
  Bool clone(const String& imageName, const String& newImageName)
  {
    //if(!valid()) return false;
    // This is not needed if(!assertDefinedImageParameters()) return false;
    LogIO os(LogOrigin("InteractiveMasking", "clone()", WHERE));
    try {
      PagedImage<Float> oldImage(imageName);
      PagedImage<Float> newImage(TiledShape(oldImage.shape(), 
					    oldImage.niceCursorShape()), oldImage.coordinates(),
				 newImageName);
      newImage.set(0.0);
      newImage.table().flush(true, true);
    } catch (AipsError x) {
      os << LogIO::SEVERE << "Exception: " << x.getMesg() << LogIO::POST;
      return false;
    } 
    return true;
  }
  
  bool new_interactive_clean_callback::callback( const DBus::Message &msg ) {
    if (msg.is_signal("edu.nrao.casa.viewer","interact")) {
      DBus::MessageIter ri = msg.reader( );
      ::operator >>(ri,result_);
      casa::DBusSession::instance( ).dispatcher( ).leave( );
    }
    return true;
  }

  /*  
  Int InteractiveMasking::interactivemask(const String& image, const String& mask, 
					  Int& niter, Int& ncycles, String& thresh, 
					  const Bool forceReload)
  {
    
    LogIO os(LogOrigin("InteractiveMasking", "interactivemask()", WHERE));
    if(Table::isReadable(mask)) {
      if (! Table::isWritable(mask)) {
	os << LogIO::WARN << "Mask image is not modifiable " << LogIO::POST;
	return false;
      }
      //we should regrid here if image and mask do not match
    }
    else{
      clone(image, mask);
    }
    
    os << "Loading image: " << image << " mask: " << mask << LogIO::POST;

    if ( viewer_p == 0 ) {
      std::list<std::string> args;
      args.push_back("--oldregions");
      viewer_p = dbus::launch<ViewerProxy>(args);
      if ( viewer_p == 0 ) {
	os << LogIO::WARN << "failed to launch viewer gui" << LogIO::POST;
	return false;
      }
    }
    if ( clean_panel_p == 0) {
      dbus::variant panel_id = viewer_p->panel( "clean" );
      if ( panel_id.type() != dbus::variant::INT ) {
	os << LogIO::WARN << "failed to create clean panel" << LogIO::POST;
	return false;
      }
      clean_panel_p = panel_id.getInt( );
    }
    
    if ( image_id_p == 0 || mask_id_p == 0 || forceReload ) {
      //Make sure image left after a "no more" is pressed is cleared
      if(forceReload && image_id_p !=0)
	prev_image_id_p=image_id_p;
      if(forceReload && mask_id_p !=0)
	prev_mask_id_p=mask_id_p;
      if(prev_image_id_p){
	viewer_p->unload( prev_image_id_p );
      }
      if(prev_mask_id_p)
	viewer_p->unload( prev_mask_id_p );
      prev_image_id_p=0;
      prev_mask_id_p=0;
      dbus::variant image_id = viewer_p->load(image, "raster", clean_panel_p);
      if ( image_id.type() != dbus::variant::INT ) {
	os << LogIO::WARN << "failed to load image" << LogIO::POST;
	return false;
      }
      image_id_p = image_id.getInt( );
      
      dbus::variant mask_id = viewer_p->load(mask, "contour", clean_panel_p);
      if ( mask_id.type() != dbus::variant::INT ) {
	os << "failed to load mask" << LogIO::WARN << LogIO::POST;
	return false;
      }
      mask_id_p = mask_id.getInt( );
    } else {
      //viewer_p->reload( clean_panel_p );
      viewer_p->reload(image_id_p);
      viewer_p->reload(mask_id_p);
    }
    
    
    casa::dbus::record options;
    options.insert("niter", niter);
    options.insert("ncycle", ncycles);
    options.insert("threshold", thresh);  
    viewer_p->setoptions(options, clean_panel_p);
    
    new_interactive_clean_callback *mycb = new new_interactive_clean_callback( );
    DBus::MessageSlot filter;
    filter = new DBus::Callback<new_interactive_clean_callback,bool,const DBus::Message &>( mycb, &new_interactive_clean_callback::callback );
    casa::DBusSession::instance( ).connection( ).add_filter( filter );
    casa::dbus::variant res = viewer_p->start_interact( dbus::variant(), clean_panel_p);
    
    //casa::DBusSession::instance( ).dispatcher( ).set_responsiveness(10000.0, 10.0);
    casa::DBusSession::instance( ).dispatcher( ).enter( );
    casa::DBusSession::instance( ).connection( ).remove_filter( filter );
    casa::dbus::variant interact_result = mycb->result( );
    delete mycb;
    
    
    int result = 0;
    if ( interact_result.type() == dbus::variant::RECORD ) {
      const dbus::record  &rec = interact_result.getRecord( );
      for ( dbus::record::const_iterator iter = rec.begin(); iter != rec.end(); ++iter ) {
	if ( iter->first == "action" ) {
	  if ( iter->second.type( ) != dbus::variant::STRING ) {
	    os << "ill-formed action result" << LogIO::WARN << LogIO::POST;
	    return false;
	  } else {
	    const std::string &action = iter->second.getString( );
	    if ( action == "stop" )
	      result = 2;
	    else if ( action == "no more" )
	      result = 1;
	    else if ( action == "continue" )
	      result = 0;
	    else {
	      os << "ill-formed action result" << LogIO::WARN << LogIO::POST;
	      return false;
	    }
	  }
	} else if ( iter->first == "ncycle" ) {
	  if ( iter->second.type( ) != dbus::variant::INT ) {
	    os << "ill-formed ncycle result" << LogIO::WARN << LogIO::POST;
	    return false;
	  } else {
	    ncycles = iter->second.getInt( );
	  }
	} else if ( iter->first == "niter" ) {
	  if ( iter->second.type( ) != dbus::variant::INT ) {
	    os << "ill-formed niter result" << LogIO::WARN << LogIO::POST;
	    return false;
	  } else {
	    niter = iter->second.getInt( );
	  }
	} else if ( iter->first == "threshold" ) {
	  if ( iter->second.type( ) != dbus::variant::STRING ) {
	    os << "ill-formed threshold result" << LogIO::WARN << LogIO::POST;
	    return false;
	  } else {
	    thresh = iter->second.getString( );
	  }
	}
      }
    } else {
      os << "failed to get a vaild result for viewer" << LogIO::WARN << LogIO::POST;
      return false;
    }
    prev_image_id_p=image_id_p;
    prev_mask_id_p=mask_id_p;
    
    if(result==1){
      //Keep the image up but clear the next time called
      image_id_p=0;
      mask_id_p=0;
    }
    if(result==2){
      //clean up
      //viewer_p->close(clean_panel_p);
      //viewer_p->done();
      //delete viewer_p;
      //viewer_p=0;
      //viewer_p->unload(image_id_p);
      //viewer_p->unload(mask_id_p);
      //Setting clean_panel_p to 0 seems to do the trick...the above stuff 
      // like done causes a crash after a call again...have to understand that
      viewer_p->unload(image_id_p);
      viewer_p->unload(mask_id_p);
      viewer_p->close(clean_panel_p);
      clean_panel_p=0;
      image_id_p=0;
      mask_id_p=0;
    }
    
    // return 0 if "continue"
    // return 1 if "no more interaction"
    // return 2 if "stop"
    return result;
  }
  */

  Int InteractiveMasking::interactivemask(const String& image, const String& mask, 
					  Int& niter, Int& cycleniter, 
					  String& thresh, String& cyclethresh, 
					  const Bool forceReload)
  {
    
    LogIO os(LogOrigin("InteractiveMasking", "interactivemask()", WHERE));
    if(Table::isReadable(mask)) {
      if (! Table::isWritable(mask)) {
	os << LogIO::WARN << "Mask image is not modifiable " << LogIO::POST;
	return false;
      }
      //we should regrid here if image and mask do not match
    }
    else{
      clone(image, mask);
    }
    
    //os << "Loading image: " << image << " mask: " << mask << LogIO::POST;

    Double startmask = maskSum(mask);

    if ( viewer_p == 0 ) {
      std::list<std::string> args;
      args.push_back("--oldregions");
      viewer_p = dbus::launch<ViewerProxy>(args);
      if ( viewer_p == 0 ) {
	os << LogIO::WARN << "failed to launch viewer gui" << LogIO::POST;
	return false;
      }
    }
    if ( clean_panel_p == 0) {
      dbus::variant panel_id = viewer_p->panel( "clean2" );
      if ( panel_id.type() != dbus::variant::INT ) {
	os << LogIO::WARN << "failed to create clean panel" << LogIO::POST;
	return false;
      }
      clean_panel_p = panel_id.getInt( );
    }
    
    if ( image_id_p == 0 || mask_id_p == 0 || forceReload ) {
      //Make sure image left after a "no more" is pressed is cleared
      if(forceReload && image_id_p !=0)
	prev_image_id_p=image_id_p;
      if(forceReload && mask_id_p !=0)
	prev_mask_id_p=mask_id_p;
      if(prev_image_id_p){
	viewer_p->unload( prev_image_id_p );
      }
      if(prev_mask_id_p)
	viewer_p->unload( prev_mask_id_p );
      prev_image_id_p=0;
      prev_mask_id_p=0;
      dbus::variant image_id = viewer_p->load(image, "raster", clean_panel_p);
      if ( image_id.type() != dbus::variant::INT ) {
	os << LogIO::WARN << "failed to load residual image in viewer : " << image << LogIO::POST;
	return false;
      }
      image_id_p = image_id.getInt( );
      
      dbus::variant mask_id = viewer_p->load(mask, "contour", clean_panel_p);
      if ( mask_id.type() != dbus::variant::INT ) {
	os << "failed to load mask in viewer : " << mask << LogIO::WARN << LogIO::POST;
	return false;
      }
      mask_id_p = mask_id.getInt( );
    } else {
      //viewer_p->reload( clean_panel_p );
      viewer_p->reload(image_id_p);
      viewer_p->reload(mask_id_p);
    }
    
    
    casa::dbus::record options;
    options.insert("niter", niter);
    options.insert("cycleniter", cycleniter);
    options.insert("threshold", thresh);  
    options.insert("cyclethreshold", cyclethresh);  
    viewer_p->setoptions(options, clean_panel_p);
    
    new_interactive_clean_callback *mycb = new new_interactive_clean_callback( );
    DBus::MessageSlot filter;
    filter = new DBus::Callback<new_interactive_clean_callback,bool,const DBus::Message &>( mycb, &new_interactive_clean_callback::callback );
    casa::DBusSession::instance( ).connection( ).add_filter( filter );
    casa::dbus::variant res = viewer_p->start_interact( dbus::variant(), clean_panel_p);
    
    //casa::DBusSession::instance( ).dispatcher( ).set_responsiveness(10000.0, 10.0);
    casa::DBusSession::instance( ).dispatcher( ).enter( );
    casa::DBusSession::instance( ).connection( ).remove_filter( filter );
    casa::dbus::variant interact_result = mycb->result( );
    delete mycb;
    
    
    int result = 1;
    if ( interact_result.type() == dbus::variant::RECORD ) {
      const dbus::record  &rec = interact_result.getRecord( );
      for ( dbus::record::const_iterator iter = rec.begin(); iter != rec.end(); ++iter ) {
	if ( iter->first == "action" ) {
	  if ( iter->second.type( ) != dbus::variant::STRING ) {
	    os << "ill-formed action result" << LogIO::WARN << LogIO::POST;
	    return false;
	  } else {
	    const std::string &action = iter->second.getString( );
	    if ( action == "stop" )
	      result = 3;
	    else if ( action == "no more" )
	      result = 2;
	    else if ( action == "continue" )
	      result = 1;
	    else {
	      os << "ill-formed action result" << LogIO::WARN << LogIO::POST;
	      return false;
	    }
	  }
	} else if ( iter->first == "cycleniter" ) {
	  if ( iter->second.type( ) != dbus::variant::INT ) {
	    os << "ill-formed cycleniter result" << LogIO::WARN << LogIO::POST;
	    return false;
	  } else {
	    cycleniter = iter->second.getInt( );
	  }
	} else if ( iter->first == "niter" ) {
	  if ( iter->second.type( ) != dbus::variant::INT ) {
	    os << "ill-formed niter result" << LogIO::WARN << LogIO::POST;
	    return false;
	  } else {
	    niter = iter->second.getInt( );
	  }
	} else if ( iter->first == "threshold" ) {
	  if ( iter->second.type( ) != dbus::variant::STRING ) {
	    os << "ill-formed threshold result" << LogIO::WARN << LogIO::POST;
	    return false;
	  } else {
	    thresh = iter->second.getString( );
	  }
	} else if ( iter->first == "cyclethreshold" ) {
	  if ( iter->second.type( ) != dbus::variant::STRING ) {
	    os << "ill-formed cyclethreshold result" << LogIO::WARN << LogIO::POST;
	    return false;
	  } else {
	    cyclethresh = iter->second.getString( );
	  }
	}
	//	if( niter<cycleniter ) 
	// { os << "niter must be >= cycleniter" << LogIO::WARN << LogIO::POST; return false; }
      }
    } else {
      os << "failed to get a vaild result for viewer" << LogIO::WARN << LogIO::POST;
      return false;
    }
    prev_image_id_p=image_id_p;
    prev_mask_id_p=mask_id_p;
    
    if(result==2){
      //Keep the image up but clear the next time called
      image_id_p=0;
      mask_id_p=0;
    }

    /*
    if(result==3){
      //clean up
      //viewer_p->close(clean_panel_p);
      //viewer_p->done();
      //delete viewer_p;
      //viewer_p=0;
      //viewer_p->unload(image_id_p);
      //viewer_p->unload(mask_id_p);
      //Setting clean_panel_p to 0 seems to do the trick...the above stuff 
      // like done causes a crash after a call again...have to understand that
      viewer_p->unload(image_id_p);
      viewer_p->unload(mask_id_p);
      viewer_p->close(clean_panel_p);
      clean_panel_p=0;
      image_id_p=0;
      mask_id_p=0;
    }
    */

    //    viewer_p->unload(image_id_p);
    //    viewer_p->unload(mask_id_p);
    image_id_p=0;
    mask_id_p=0;


    Double endmask = maskSum(mask);

    
    if( startmask != endmask)
      {
	result = -1*result;
        LogIO os( LogOrigin("InteractiveMasking","interactivemask",WHERE) );
	os << "[" << mask << "] Mask modified from " << startmask << " pixels to " << endmask << " pixels " << LogIO::POST;
      }
    //    else
      //cout << " Mask did not change. Result :  " << result << endl;
    

    // return 1 if "continue"
    // return 2 if "no more interaction"
    // return 3 if "stop"
    return result;
  }


  Float InteractiveMasking::maskSum(const String & maskname)
  {

    PagedImage<Float> *mask = new PagedImage<Float> (maskname);

    LatticeExprNode msum( sum( *mask ) );
    Float maskSum = msum.getFloat();

    mask->unlock();
    mask->tempClose();

    delete mask;
    mask = NULL;

    //cout << "Mask Sum : " << maskSum << endl;

    return maskSum;

    /*
    Record*  regionPtr=0;
    String LELmask("");

    ImageStatsCalculator imcalc( imask, regionPtr, LELmask, False);
    Record thestats = imcalc.statistics();

    Array<Double> msum;
    thestats.get(RecordFieldId("sum"), msum);

    cout << "mask sum : " << msum << endl;

    AlwaysAssert( msum.nelements()>0, AipsError );
    
    return sum(msum);
    */

  }

};
