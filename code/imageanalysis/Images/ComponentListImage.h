//# ComponentListImage.h
//# Copyright (C) 1994,1995,1996,1997,1998,1999,2000,2001,2003
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

#ifndef IMAGES_COMPONENTlISTIMAGE_H
#define IMAGES_COMPONENTlISTIMAGE_H


//# Includes
#include <casacore/casa/aips.h>
#include <casacore/images/Images/ImageInterface.h>
#include <casacore/images/Images/ImageAttrHandlerCasa.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/TableRecord.h>

#include <components/ComponentModels/ComponentList.h>

namespace casa {

// <summary>
// Read, store, and manipulate an astronomical image based on a component list.
// </summary>

// <prerequisite>
//   <li> <linkto class=CoordinateSystem>CoordinateSystem</linkto>
//   <li> <linkto class=ImageInterface>ImageInterface</linkto>
//   <li> <linkto class=Lattice>Lattice</linkto>
//   <li> <linkto class=LatticeIterator>LatticeIterator</linkto>
//   <li> <linkto class=LatticeNavigator>LatticeNavigator</linkto>
//   <li> <linkto class=ImageRegion>ImageRegion</linkto>
// </prerequisite>

// <etymology>

// </etymology>

// <synopsis> 

// </synopsis> 

// <example>

// </example>

// <motivation>

// </motivation>

class ComponentListImage: public casacore::ImageInterface<casacore::Float> {
public:

    // Create a new ComponentListImage from a component list, coordinate system, and shape.
    // The input component list will be copied (via the ComponentList::copy() method, so
    // that the constructed object has a unique ComponentList that is not referenced by
    // Exceptions thrown if shape.size() != 4 or csys.nPixelAxes != 4 and csys is missing
    // any of Direction, Spectral, or Stokes coordinates. If tableName is not empty,
    // a new table by that name is created on disk to persistently store the image;
    // otherwise, the image is not persistent.
    ComponentListImage(
        const ComponentList& compList, const casacore::CoordinateSystem& csys,
        const casacore::IPosition& shape, const casacore::String& tableName=""
    );

    // Construct an object by reading a persistent ComponentListImage from disk.
    explicit ComponentListImage(const casacore::String& filename);

    ComponentListImage(const ComponentListImage& image);

    ~ComponentListImage();

    casacore::ImageInterface<casacore::Float>* cloneII() const;

    bool doGetSlice(
        casacore::Array<casacore::Float>& buffer, const casacore::Slicer& section
    );

    void doPutSlice (
        const casacore::Array<casacore::Float>& buffer,
        const casacore::IPosition& where, const casacore::IPosition& stride
    );

    const casacore::LatticeRegion* getRegionPtr() const;

    casacore::String imageType() const;

    casacore::String name (bool stripPath=false) const;

    casacore::Bool ok() const;

    void resize (const casacore::TiledShape& newShape);

    // Flushes the new coordinate system to disk if the table exists is writable.
    casacore::Bool setCoordinateInfo (const casacore::CoordinateSystem& coords);

    casacore::Bool setImageInfo(const casacore::ImageInfo& info);

    casacore::Bool setMiscInfo(const casacore::RecordInterface& newInfo);

    casacore::Bool setUnits (const casacore::Unit& newUnits);

    casacore::IPosition shape() const;


    /*
  
  
  // Assignment operator (reference semantics).
  PagedImage<T>& operator= (const PagedImage<T>& other);
  
  // Make a copy of the object (reference semantics).
  virtual ImageInterface<T>* cloneII() const;

  // Return the name of this derived class.
  static String className()
    { return "PagedImage"; }

  // Get the image type (returns name of this derived class).
  virtual String imageType() const;

  // A PagedImage is always persistent.
  virtual Bool isPersistent() const;

  // A PagedImage is always paged to disk.
  virtual Bool isPaged() const;

  // Is the PagedImage writable?
  virtual Bool isWritable() const;

  // Does the image object use a pixelmask?
  virtual Bool hasPixelMask() const;

  // Get access to the pixelmask used.
  // An exception is thrown if the image does not use a pixelmask.
  // <group>
  virtual const Lattice<Bool>& pixelMask() const;
  virtual Lattice<Bool>& pixelMask();
  // </group>

  // Get a pointer the default pixelmask object used with this image.
  // It returns 0 if no default pixelmask is used.
  virtual const LatticeRegion* getRegionPtr() const;

  // Set the default pixelmask to the mask with the given name
  // (which has to exist in the "masks" group).
  // If the image table is writable, the setting is persistent by writing
  // the name as a keyword.
  // If the given regionName is the empty string,
  // the default pixelmask is unset.
  virtual void setDefaultMask (const String& maskName);

  // Use the mask as specified.
  // If a mask was already in use, it is replaced by the new one.
  virtual void useMask (MaskSpecifier = MaskSpecifier());

  // Function to change the name of the Table file on disk.
  // PagedImage is given a file name at construction time.  You may change
  // that name here.
  void rename (const String& newName);

  // Return the current Table name. By default this includes the full path. 
  // the path preceding the file name can be stripped off on request.
  virtual String name (Bool stripPath=False) const;

  // Return the current TableColumn row number.
  uInt rowNumber() const;

  // Return the shape of the image.
  virtual IPosition shape() const;

  // Change the shape of the image (N.B. the data is thrown away).
  virtual void resize (const TiledShape& newShape);

  // Function which extracts an array from the map.
  virtual Bool doGetSlice (Array<T>& buffer, const Slicer& theSlice);
  
  // Function to replace the values in the map with soureBuffer.
  virtual void doPutSlice (const Array<T>& sourceBuffer,
			   const IPosition& where,
			   const IPosition& stride);

  // Replace every element, x, of the lattice with the result of f(x).
  // you must pass in the address of the function -- so the function
  // must be declared and defined in the scope of your program.  
  // Both versions of apply require a function that accepts a single 
  // argument of type T (the Lattice template actual type) and returns
  // a result of the same type.  The first apply expects a function with
  // an argument passed by value; the second expects the argument to
  // be passed by const reference.  The first form ought to run faster
  // for the built-in types, which may be an issue for large Lattices
  // stored in memory, where disk access is not an issue.
  // <group>
  virtual void apply (T (*function)(T));
  virtual void apply (T (*function)(const T& ));
  virtual void apply (const Functional<T,T>& function);
  // </group>

  // Add a lattice to this image.
  PagedImage<T>& operator+= (const Lattice<T>& other);

  // Function which sets the units associated with the image
  // pixels (i.e. the "brightness" unit). <src>setUnits()</src> returns
  // False if it cannot set the unit for some reason (e.g. the underlying
  // file is not writable).
  virtual Bool setUnits (const Unit& newUnits);

  // Return the table holding the data.
  Table& table()
    { return map_p.table(); }

  // Flushes the new coordinate system to disk if the table is writable.
  virtual Bool setCoordinateInfo (const CoordinateSystem& coords);

  // Check for symmetry in data members.
  virtual Bool ok() const;

  // These are the true implementations of the paran operator.
  // <note> Not for public use </note>
  // <group>
  virtual T getAt (const IPosition& where) const;
  virtual void putAt (const T& value, const IPosition& where);
  // </group>

  // Replace the miscinfo in the PagedImage.
  // It can fail if, e.g., the underlying table is not writable.
  virtual Bool setMiscInfo (const RecordInterface& newInfo);

  // The ImageInfo object contains some miscellaneous information about the
  // image, which unlike that stored in MiscInfo, has a standard list of
  // things, such as the restoring beam.
  // Note that setImageInfo REPLACES the information with the new information.
  // It can fail if, e.g., the underlying table is not writable.
  virtual Bool setImageInfo(const ImageInfo& info);

  // Get access to the attribute handler.
  // If a handler keyword does not exist yet, it is created if
  // <src>createHandler</src> is set.
  // Otherwise the handler is empty and no groups can be created for it.
  virtual ImageAttrHandler& attrHandler (Bool createHandler=False);

  // Remove a region/mask belonging to the image from the given group
  // (which can be Any).
  // If a mask removed is the default mask, the image gets unmasked.
  // <br>Optionally an exception is thrown if the region does not exist.
  virtual void removeRegion (const String& name,
			     RegionHandler::GroupType = RegionHandler::Any,
			     Bool throwIfUnknown = True);

  // This is the implementation of the letter for the envelope Iterator
  // class. <note> Not for public use </note>.
  virtual LatticeIterInterface<T>* makeIter
                                 (const LatticeNavigator& navigator,
				  Bool useRef) const;

  // Returns the maximum recommended number of pixels for a cursor. This is
  // the number of pixels in a tile. 
  virtual uInt advisedMaxPixels() const;

  // Help the user pick a cursor for most efficient access.
  virtual IPosition doNiceCursorShape (uInt maxPixels) const;

  // Maximum size - not necessarily all used. In pixels.
  virtual uInt maximumCacheSize() const;

  // Set the maximum (allowed) cache size as indicated.
  virtual void setMaximumCacheSize (uInt howManyPixels);

  // Set the cache size as to "fit" the indicated path.
  virtual void setCacheSizeFromPath (const IPosition& sliceShape,
  			             const IPosition& windowStart,
			             const IPosition& windowLength,
			             const IPosition& axisPath);
    
  // Set the actual cache size for this Array to be be big enough for the
  // indicated number of tiles. This cache is not shared with PagedArrays
  // in other rows and is always clipped to be less than the maximum value
  // set using the setMaximumCacheSize member function.
  // tiles. Tiles are cached using a first in first out algorithm. 
  virtual void setCacheSizeInTiles (uInt howManyTiles);

  // Clears and frees up the caches, but the maximum allowed cache size is 
  // unchanged from when setCacheSize was called
  virtual void clearCache();

  // Report on cache success.
  virtual void showCacheStatistics (ostream& os) const;

  // Handle the (un)locking.
  // Unlocking also unlocks the logtable and a possible mask table.
  // Locking only locks the image itself.
  // <group>
  virtual Bool lock (FileLocker::LockType, uInt nattempts);
  virtual void unlock();
  virtual Bool hasLock (FileLocker::LockType) const;
  // </group>

  // Resynchronize the PagedImage object with the table contents.
  // The logtable and possible mask table are also synchronized if
  // they do not have a readlock.
  // <br>This function is only useful if no read-locking is used, ie.
  // if the table lock option is UserNoReadLocking or AutoNoReadLocking.
  // In that cases the table system does not acquire a read-lock, thus
  // does not synchronize itself automatically.
  virtual void resync();

  // Flush the data.
  virtual void flush();

  // Close the Image and associated files temporarily.
  // It'll be reopened automatically when needed or when
  // <src>reopen</src> is called explicitly.
  virtual void tempClose();

  // If needed, reopen a temporarily closed Image.
  virtual void reopen();
*/

private:

    ComponentList _cl;
    casacore::IPosition _shape = casacore::IPosition(4, 1);
    std::shared_ptr<casacore::Lattice<bool>> _mask = nullptr;

    void _openLogTable();

    void _reopenRW();

    void _restoreAll(const casacore::TableRecord& rec);

    void _restoreImageInfo(const casacore::TableRecord& rec);

    void _restoreMiscInfo(const casacore::TableRecord& rec);

    void _restoreUnits(const casacore::TableRecord& rec);

    /*
  // Function to return the internal Table object to the RegionHandler.
  static Table& getTable (void* imagePtr, Bool writable);

  // This must be called in every constructor and place where the image
  // is attached to a new image.
  void attach_logtable();
  void open_logtable();
  void restoreUnits (const TableRecord& rec);
  void restoreMiscInfo (const TableRecord& rec);
  void restoreImageInfo (const TableRecord& rec);
  void restoreAll (const TableRecord& rec);

  void check_conformance (const Lattice<T>& other);
  void reopenRW();
  void setTableType();
  void applyMaskSpecifier (const MaskSpecifier&);
  void applyMask (const String& maskName);
  void makePagedImage (const TiledShape& mapShape,
		       const CoordinateSystem& coordinateInfo,
		       const String& nameOfNewFile,
		       const TableLock& lockOptions,
		       uInt rowNumber);
  void makePagedImage (const String& filename, const TableLock& lockOptions,
		       const MaskSpecifier&, uInt rowNumber);

  const Table& table() const
    { return const_cast<PagedImage<T>*>(this)->table(); }


  PagedArray<T>  map_p;
  LatticeRegion* regionPtr_p;
  ImageAttrHandlerCasa itsAttrHandler;

  //# Make members of parent class known.
public:
  using ImageInterface<T>::logSink;
  using ImageInterface<T>::logger;
  using ImageInterface<T>::imageInfo;
  using ImageInterface<T>::coordinates;
  using ImageInterface<T>::getDefaultMask;
  using ImageInterface<T>::hasRegion;
  using ImageInterface<T>::getImageRegionPtr;
protected:
  using ImageInterface<T>::setCoordsMember;
  using ImageInterface<T>::setMiscInfoMember;
  using ImageInterface<T>::setLogMember;
  using ImageInterface<T>::setUnitMember;
  using ImageInterface<T>::setImageInfoMember;
  */
};

/*
//# A nasty - the column name is hard-coded into this function, needs to
//# be centralized somewhere.
// Determine the pixel type in the PagedImage contained in
// <src>fileName</src>.  If the file doesn't appear to be a Table or cannot
// be opened, TpOther is returned.
// <group name="pixeltype")
    DataType imagePixelType(const String& fileName);
// </group>


//# Declare extern templates for often used types.
#ifdef AIPS_CXX11
  extern template class PagedImage<Float>;
  extern template class PagedImage<Complex>;
#endif
*/

}
/*
#ifndef CASACORE_NO_AUTO_TEMPLATES
#include <casacore/images/Images/PagedImage.tcc>
#endif //# CASACORE_NO_AUTO_TEMPLATES

*/
#endif
