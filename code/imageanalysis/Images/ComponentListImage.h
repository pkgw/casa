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
#include <casacore/casa/Quanta/MVAngle.h>
#include <casacore/images/Images/ImageInterface.h>
#include <casacore/images/Images/TempImage.h>
#include <casacore/lattices/Lattices/TempLattice.h>
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

    static const casacore::String IMAGE_TYPE;

    // Create a new ComponentListImage from a component list, coordinate system, and shape.
    // The input component list will be copied (via the ComponentList::copy() method, so
    // that the constructed object has a unique ComponentList that is not referenced by
    // Exceptions thrown if shape.size() != 4 or csys.nPixelAxes != 4 or csys is missing
    // any of Direction, Spectral, or Stokes coordinates. If tableName is not empty,
    // a new table by that name is created on disk to persistently store the image;
    // otherwise, the image is not persistent. The brightness units are automatically
    // set to Jy/pixel.
    ComponentListImage(
        const ComponentList& compList, const casacore::CoordinateSystem& csys,
        const casacore::IPosition& shape, const casacore::String& tableName="",
        const casacore::Bool doCache=casacore::True
    );

    // Construct an object by reading a persistent ComponentListImage from disk.
    explicit ComponentListImage(
        const casacore::String& filename,
        const casacore::Bool doCache=casacore::True
    );

    ComponentListImage(const ComponentListImage& image);

    ~ComponentListImage();

    ComponentListImage& operator=(const ComponentListImage& other);

    // <group>
    // These override the methods in Lattice and always throw an exception because
    // there is not a general way to do this for a component list.
    void apply (casacore::Float (*function)(casacore::Float));
    void apply (casacore::Float (*function)(const casacore::Float&));
    void apply (const casacore::Functional<casacore::Float, casacore::Float>& function);
    // </group>

    casacore::ImageInterface<casacore::Float>* cloneII() const;

    // get a reference to the underlying ComponentList
    const ComponentList& componentList() const;

    casacore::Bool doGetMaskSlice(
        casacore::Array<casacore::Bool>& buffer, const casacore::Slicer& section
    );

    // <group>
    // These methods override those in Lattice and always throw exceptions as there is not
    // general way to do this for component lists.
    void copyData (const casacore::Lattice<casacore::Float>&);
    void copyDataTo (casacore::Lattice<casacore::Float>&) const;
    // </group>

    bool doGetSlice(
        casacore::Array<casacore::Float>& buffer, const casacore::Slicer& section
    );

    void doPutSlice (
        const casacore::Array<casacore::Float>& buffer,
        const casacore::IPosition& where, const casacore::IPosition& stride
    );

    const casacore::LatticeRegion* getRegionPtr() const;

    casacore::Bool hasPixelMask() const;

    casacore::String imageType() const;

    casacore::Bool isMasked() const;

    casacore::Bool isPersistent() const;

    // Overrides the LatticeBase method. The method refers to the Lattice
    // characteristic and therefore always returns False.
    casacore::Bool isWritable() const;

    casacore::String name (bool stripPath=false) const;

    casacore::Bool ok() const;

    // Open an existing ComponentListImage. This gets registered with ImageOpener. It simply calls and
    // returns the result of new ComponentListImage(name).
    static casacore::LatticeBase* openCLImage(const casacore::String& name, const casacore::MaskSpecifier&);

    const casacore::Lattice<casacore::Bool>& pixelMask() const;

    static void registerOpenFunction();

    // overrides ImageInterface method
    void removeRegion(
        const casacore::String& name,
        casacore::RegionHandler::GroupType=casacore::RegionHandler::Any,
        casacore::Bool throwIfUnknown=casacore::True
    );

    void resize(const casacore::TiledShape& newShape);

    // Overrides Lattice method. Always throws exception; there is no way to represent a source with constant
    // intensity everywhere.
    void set(const casacore::Float& value);

    // controls if pixel values, directions, and frequencies are cached. If doCache is false,
    // existing cache values are destroyed if they exist.
    void setCache(casacore::Bool doCache);

    // Flushes the new coordinate system to disk if the table exists is writable.
    casacore::Bool setCoordinateInfo (const casacore::CoordinateSystem& coords);

    void setDefaultMask(const casacore::String& regionName);

    casacore::Bool setImageInfo(const casacore::ImageInfo& info);

    casacore::Bool setMiscInfo(const casacore::RecordInterface& newInfo);

    // An exception is thrown if newUnits are anything but Jy/pixel.
    casacore::Bool setUnits (const casacore::Unit& newUnits);

    // return the shape of the image.
    casacore::IPosition shape() const;

    // Overrides ImageInterface method.
    void useMask(casacore::MaskSpecifier=casacore::MaskSpecifier());

protected:
    //Overrides Lattice. Always throws exception; arbitrary math cannot be represented with components.
    void handleMath(const casacore::Lattice<casacore::Float>& from, int oper);

    //Overrides Lattice. Always throws exception; arbitrary math cannot be represented with components.
    void handleMathTo (casacore::Lattice<casacore::Float>& to, int oper) const;

private:

    ComponentList _cl, _modifiedCL;

    casacore::IPosition _shape = casacore::IPosition(4, 1);
    std::shared_ptr<casacore::Lattice<bool>> _mask = nullptr;
    casacore::Int _latAxis, _longAxis, _freqAxis, _stokesAxis;
    casacore::MVAngle _pixelLongSize, _pixelLatSize;
    casacore::MeasRef<casacore::MDirection> _dirRef;
    casacore::Matrix<std::shared_ptr<casacore::MVDirection>> _dirVals;
    casacore::MeasRef<casacore::MFrequency> _freqRef;
    casacore::Vector<std::shared_ptr<casacore::MVFrequency>> _freqVals;
    std::shared_ptr<
        std::map<casacore::IPosition, casacore::Float, casacore::IPositionComparator>
    > _ptSourcePixelVals = nullptr;
    casacore::Vector<casacore::Int> _pixelToIQUV;
    casacore::Bool _cache = casacore::True;
    std::shared_ptr<casacore::TempImage<casacore::Float>> _valueCache = nullptr;

    void _applyMask (const casacore::String& maskName);

    void _applyMaskSpecifier (const casacore::MaskSpecifier& spec);

    void _cacheCoordinateInfo(const casacore::CoordinateSystem& csys);

    // improve performance by caching point source pixel values. This is always done;
    // the value of _cache has no effect on it.
    void _computePointSourcePixelValues();

    void _deleteCache();

    void _fillBuffer(
        casacore::Array<casacore::Float>& buffer, const casacore::IPosition& chunkShape,
        const casacore::IPosition& secStart, casacore::Bool lookForPtSources,
        casacore::uInt nFreqs, const casacore::Cube<casacore::Double>& pixelVals
    ) const;

    casacore::Bool _findPixel(
        casacore::Cube<casacore::Double>& values, const casacore::IPosition& pixelPosition,
        const casacore::DirectionCoordinate& dirCoord, const SkyComponent& point,
        const casacore::Vector<casacore::MVFrequency>& freqValues
    ) const;

    casacore::Bool _findPixelIn3x3Box(
        casacore::IPosition& pixelPosition, casacore::Cube<casacore::Double>& values,
        const casacore::DirectionCoordinate& dirCoord, const SkyComponent& point,
        const casacore::Vector<casacore::MVFrequency>& freqValues
    ) const;

    casacore::Vector<casacore::MVFrequency> _getAllFreqValues(casacore::uInt nFreqs);

    void _getDirValsDoCache(
        casacore::Vector<casacore::MVDirection>& dirVals,
        const casacore::IPosition& secStart, casacore::uInt endLong,
        casacore::uInt endLat, const casacore::DirectionCoordinate& dirCoord
    );

    void _getDirValsNoCache(
        casacore::Vector<casacore::MVDirection>& dirVals,
        const casacore::IPosition& secStart, casacore::uInt endLong,
        casacore::uInt endLat, const casacore::DirectionCoordinate& dirCoord
    ) const;

    void _getFreqValsDoCache(
        casacore::Vector<casacore::MVFrequency>& freqVals,
        const casacore::IPosition& secStart, casacore::uInt nFreqs,
        const casacore::SpectralCoordinate& specCoord
    );

    void _getFreqValsNoCache(
        casacore::Vector<casacore::MVFrequency>& freqVals,
        const casacore::IPosition& secStart, casacore::uInt nFreqs,
        const casacore::SpectralCoordinate& specCoord
    ) const;

    void _initCache();

    void _openLogTable();

    void _reopenRW();

    void _resize(const casacore::TiledShape& shape);

    void _restoreAll(const casacore::TableRecord& rec);

    void _restoreImageInfo(const casacore::TableRecord& rec);

    void _restoreMiscInfo(const casacore::TableRecord& rec);

    void _restoreUnits(const casacore::TableRecord& rec);

};

}

#endif
