// -*- C++ -*-
//# Framework independent implementation file for ms..
//# Copyright (C) 2006-2007-2008
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
//# @author
//# @version
//////////////////////////////////////////////////////////////////////////////

#include <msmetadata_cmpt.h>

#include <tools/ms/msmetadata_forward.h>

#include <casa/Containers/Record.h>
#include <casa/BasicSL/STLIO.h>
#include <casa/Logging/LogIO.h>
#include <casa/Quanta/QuantumHolder.h>
#include <casa/Quanta/QLogical.h>
#include <casa/Quanta/QVector.h>
#include <measures/Measures/MeasureHolder.h>
#include <measures/Measures/MDirection.h>
#include <ms/MeasurementSets/MeasurementSet.h>
#include <ms/MSOper/MSMetaData.h>
#include <ms/MSOper/MSKeys.h>

#include <msvis/MSVis/MSChecker.h>

#include <algorithm>
#include <regex>

#include <casa/namespace.h>

#define _ORIGIN *_log << LogOrigin("msmetadata_cmpt.cc", __func__, __LINE__);
// common method scaffold

#define COMMA ,

#define _FUNC(BODY) \
	try { \
		_ORIGIN; \
		_isAttached(); \
		BODY \
	} \
	catch (const AipsError& x) { \
		_handleException(x); \
	}

#define _FUNC2(BODY) \
	try { \
		_ORIGIN; \
		BODY \
	} \
	catch (const AipsError& x) { \
		_handleException(x); \
	}

using namespace casacore;
using namespace casa;

namespace casac {

template<class T> struct counting_iterator : public std::iterator<std::input_iterator_tag, T> {
    bool operator!=(const counting_iterator &other) { return val != other.val; }
    T operator*( ) { return val; }
    void operator++( ) { ++val; }
    counting_iterator( T v ) : val(v) { }
    T val;
};

msmetadata::msmetadata() : _log(new LogIO()) {}

msmetadata::~msmetadata() {}

vector<int> msmetadata::almaspws(
	bool chavg, bool fdm, bool sqld, bool tdm, bool wvr, bool complement
) {
	_FUNC(
	    std::set<uInt> x;
	    if (chavg) {
	    	std::set<uInt> y = _msmd->getChannelAvgSpw();
	    	x.insert(y.begin(), y.end());
	    }
	    if (fdm) {
	    	std::set<uInt> y = _msmd->getFDMSpw();
	    	x.insert(y.begin(), y.end());
	    }
	    if (sqld) {
	    	std::set<uInt> y = _msmd->getSQLDSpw();
	    	x.insert(y.begin(), y.end());
	    }
	    if (tdm) {
	    	std::set<uInt> y = _msmd->getTDMSpw();
	    	x.insert(y.begin(), y.end());
	    }
	    if (wvr) {
	    	std::set<uInt> y = _msmd->getWVRSpw();
	    	x.insert(y.begin(), y.end());
	    }
		if (complement) {
			uInt nspw = _msmd->nSpw(true);
			set<uInt> allSpws(
				counting_iterator<int>(0),
				counting_iterator<int>(nspw)
			);
			vector<uInt> mycompl(nspw);
			vector<uInt>::iterator begin = mycompl.begin();
			vector<uInt>::iterator end = std::set_difference(
				allSpws.begin(), allSpws.end(), x.begin(),
				x.end(), begin
			);
			mycompl.resize(end - begin);
			return _vectorUIntToVectorInt(mycompl);
		}
		return _setUIntToVectorInt(x);
	)
	return vector<int>();
}

vector<String> msmetadata::_vectorStdStringToVectorString(
	const vector<string>& inset
) {
	vector<String> outset;
	for(string el: inset) {
		outset.push_back(el);
	}
	return outset;
}

record* msmetadata::antennadiameter(const variant& antenna) {
	_FUNC(
		variant::TYPE type = antenna.type();
		int antID;
		if (type == variant::INT) {
			antID = antenna.toInt();
			_checkAntennaId(antID, false);
		}
		else if (type == variant::STRING) {
			antID = _msmd->getAntennaID(antenna.toString());
		}
		else {
			ThrowCc(
				"Unsupported type for input parameter antenna. "
				"Supported types are int and string"
			);
		}
		Record rec;
		Quantum<Vector<Double> > out = _msmd->getAntennaDiameters();
        if (antID >= 0) {
            QuantumHolder qh(casacore::Quantity(out.getValue()[antID], out.getFullUnit()));
		    qh.toRecord(rec);
        }
        else {
            auto i = 0;
            auto unit = out.getFullUnit();
            auto values = out.getValue();
            for (const auto& diam: values) {
                QuantumHolder qh(casacore::Quantity(diam, unit));
                Record rec1;
                qh.toRecord(rec1);
                rec.defineRecord(String::toString(i), rec1);
                ++i;
            }
        }
        return fromRecord(rec);
	)
	return nullptr;
}

vector<int> msmetadata::antennaids(
	const variant& names, const variant& mindiameter,
	const variant& maxdiameter, int obsid
) {
	_FUNC(
		_checkObsId(obsid, False);
		vector<String> myNames;
		// because variants default to boolvecs even when we specify elsewise in the xml.
		casacore::Quantity mind = mindiameter.type() == variant::BOOLVEC ?
			casacore::Quantity(0, "m") : casaQuantity(mindiameter);
		casacore::Quantity maxd = maxdiameter.type() == variant::BOOLVEC ?
			casacore::Quantity(1, "pc") :casaQuantity(maxdiameter);
		ThrowIf(
			! mind.isConform("m"),
			"mindiameter must have units of length"
		);
		ThrowIf(
			! maxd.isConform("m"),
			"maxdiameter must have units of length"
		);
        set<uInt> filteredAnts;
        auto filterObsID = obsid >= 0;
        if (filterObsID) {
            filteredAnts = _antsForObsID(obsid);
        }
		vector<Int> foundIDs;
		variant::TYPE type = names.type();
		if (type == variant::BOOLVEC) {
		    if (filterObsID) {
                foundIDs = _setUIntToVectorInt(filteredAnts);
                if (foundIDs.empty()) {
                    return vector<int>();
                }
            }
            else {
	            Vector<Int> x(_msmd->nAntennas());
			    indgen(x, 0);
			    foundIDs = x.tovector();
            }
		}
		else if (type == variant::STRING) {
			myNames.push_back(names.toString());
		}
		else if (type == variant::STRINGVEC) {
			myNames = _vectorStdStringToVectorString(names.toStringVec());
		}
		else {
			ThrowCc(
				"Unsupported type for parameter names. Must "
				"be either a string or string array"
			);
		}
		vector<casacore::String> allNames COMMA foundNames;
		set<casacore::String> foundSet;
		if (foundIDs.empty()) {
			for(String name: myNames) {
				Bool expand = name.find('*') != std::string::npos;
				if (expand) {
					if (allNames.empty()) {
						std::map<String COMMA uInt> namesToIDsMap;
						allNames = _msmd->getAntennaNames(namesToIDsMap);
					}
					vector<String> matches = _match(allNames, name);
					for(String match: matches) {
						if (foundSet.find(match) == foundSet.end()) {
							foundNames.push_back(match);
							foundSet.insert(match);
						}
					}
				}
				else {
					if (foundSet.find(name) == foundSet.end()) {
						foundNames.push_back(name);
						foundSet.insert(name);
					}
				}
			}
            auto allFoundIDs = _msmd->getAntennaIDs(foundNames);
            uInt i = 0;
            vector<uInt> r;
            for (const auto& idSet: allFoundIDs) {
                if (filterObsID) {
                    _intersection(r, idSet, filteredAnts);
                }
                else {
                    auto nEntries = idSet.size();
                    if (nEntries > 1) {
                        *_log << LogIO::WARN << "Antenna " << foundNames[i]
                            << " is listed in the ANTENNA table " << nEntries
                            << " times and is associated with IDs " << idSet
                            << ". Only the largest ID associated with it "
                            << "will be returned" << LogIO::POST;
                    }
                    foundIDs.push_back(*idSet.rbegin());
                    ++i;
                }
            }
            if (filterObsID) {
                foundIDs = _vectorUIntToVectorInt(r);
            }
		}
		Quantum<Vector<Double> > diams = _msmd->getAntennaDiameters();
		casacore::Quantity maxAntD(max(diams.getValue()), diams.getUnit());
		casacore::Quantity minAntD(min(diams.getValue()), diams.getUnit());
		if (mind > minAntD || maxd < maxAntD) {
			vector<Int> newList;
			String unit = diams.getUnit();
			Vector<Double> v = diams.getValue();
			for(uInt id: foundIDs) {
				casacore::Quantity d(v[id], unit);
				if (d >= mind && d <= maxd) {
					newList.push_back(id);
				}
			}
			foundIDs = newList;
		}
		return foundIDs;
	)
	return vector<int>();
}

std::set<uInt> msmetadata::_antsForObsID(Int obsid) const {
    auto ssprops = _msmd->getSubScanProperties(True);
    std::set<uInt> filteredAnts;
    for (const auto& prop: *ssprops) {
        if (prop.first.obsID == obsid) {
            auto ants = prop.second.antennas;
            filteredAnts.insert(ants.begin(), ants.end());
        }
    }
    return filteredAnts;
}

vector<string> msmetadata::antennanames(const variant& antennaids) {
	_FUNC(
		vector<uInt> myIDs;
		variant::TYPE type = antennaids.type();
		if (type == variant::BOOLVEC) {
			// do nothing, all names will be returned.
		}
		else if (type == variant::INT) {
			Int id = antennaids.toInt();
			ThrowIf(id < 0, "Antenna ID must be nonnegative");
			myIDs.push_back(id);
		}
		else if (type == variant::INTVEC) {
			vector<Int> tmp = antennaids.toIntVec();
			vector<Int>::const_iterator end = tmp.end();
			for (
				vector<Int>::const_iterator iter=tmp.begin();
				iter!=tmp.end(); iter++
			) {
				ThrowIf(*iter < 0, "Antenna ID must be nonnegative");
				myIDs.push_back(*iter);
			}
		}
		else {
			ThrowCc(
				"Unsupported type for parameter antennaids. "
				"Must be either an integer or integer array"
			);
		}
		std::map<String COMMA uInt> namesToIDsMap;
		return _vectorStringToStdVectorString(
			_msmd->getAntennaNames(namesToIDsMap COMMA myIDs)
		);
	)
	return vector<string>();
}

record* msmetadata::antennaoffset(const variant& which) {
	_FUNC(
		variant::TYPE type = which.type();
		Quantum<Vector<Double> > out;
		if (type == variant::INT) {
			out = _msmd->getAntennaOffset(which.toInt());
		}
		else if (type == variant::STRING) {
			out = _msmd->getAntennaOffset(which.toString());
		}
		else {
			ThrowCc(
				"Unsupported type for input parameter which. "
				"Supported types are int and string"
			);
		}
		Vector<Double> v = out.getValue();
		String u = out.getUnit();
		QuantumHolder longitude(casacore::Quantity(v[0], u));
		QuantumHolder latitude(casacore::Quantity(v[1], u));
		QuantumHolder elevation(casacore::Quantity(v[2], u));
		Record x;
		Record outRec;
		longitude.toRecord(x);
		outRec.defineRecord("longitude offset", x);
		latitude.toRecord(x);
		outRec.defineRecord("latitude offset", x);
		elevation.toRecord(x);
		outRec.defineRecord("elevation offset", x);
		return fromRecord(outRec);
	)
	return 0;
}

record* msmetadata::antennaposition(const variant& which) {
	_FUNC(
		variant::TYPE type = which.type();
		ThrowIf(
			type != variant::BOOLVEC && type != variant::INT
			&& type != variant::STRING,
			"Unsupported type for which, must be either int or string"
		);
		MeasureHolder out;
		if (type == variant::BOOLVEC || type == variant::INT) {
			out = MeasureHolder(
				_msmd->getAntennaPositions(
					type == variant::BOOLVEC ? vector<uInt>(1, 0)
					: vector<uInt>(1, which.toInt())
				)[0]
			);
		}
		else {
			out = MeasureHolder(
				*_msmd->getAntennaPositions(
					vector<String>(1, which.toString())
				)[0].rbegin()
			);
		}
		Record outRec;
		out.toRecord(outRec);
		return fromRecord(outRec);
	)
	return 0;
}

vector<string> msmetadata::antennastations(const variant& ants, int obsid) {
    _FUNC(
		_checkObsId(obsid, False);
        std::set<uInt> filteredAnts;
        auto filterObsID = obsid >= 0;
        if (filterObsID) {
            filteredAnts = _antsForObsID(obsid);
            if (filteredAnts.empty()) {
                return vector<string>();
            }
        }
		variant::TYPE type = ants.type();
		if (type == variant::INT) {
			int id = ants.toInt();
            vector<uInt> ids;
			if (id >= 0) {
			    if (filterObsID && filteredAnts.find(id) == filteredAnts.end()) {
                    // ant ID not included in obsid
                    return vector<string>();
                }
	            ids.push_back(id);
			}
            else if (filterObsID) {
                ids = _setUIntToVectorUInt(filteredAnts);
            }
			return _vectorStringToStdVectorString(
				_msmd->getAntennaStations(ids)
			);
		}
		else if (type == variant::INTVEC) {
			vector<Int> ids = ants.toIntVec();
			if (ids.empty() || (ids.size() == 1 && ids[0] < 0)) {
				return _vectorStringToStdVectorString(
					_msmd->getAntennaStations(
                        filterObsID ? _setUIntToVectorUInt(filteredAnts) : vector<uInt>()
                    )
				);
			}
			else if (*min_element(ids.begin(), ids.end()) < 0) {
				throw AipsError("No antenna ID may be less than zero when multiple IDs specified.");
			}
			else if (! filterObsID) {
				return _vectorStringToStdVectorString(
					_msmd->getAntennaStations(
						_vectorIntToVectorUInt(ants.toIntVec())
					)
				);
			}
            else {
                auto goodAnts = _intersection(_vectorIntToVectorUInt(ants.toIntVec()), filteredAnts);
                if (goodAnts.empty()) {
                    return vector<string>();
                }
                else {
                    return _vectorStringToStdVectorString(
				        _msmd->getAntennaStations(goodAnts)
				    );
                }
            }
		}
		else if (type == variant::STRING) {
            if (filterObsID) {
                auto antids = _msmd->getAntennaIDs(ants.toString());
                vector<uInt> ids;
                _intersection(ids, antids, filteredAnts);
			    return _vectorStringToStdVectorString(
				    _msmd->getAntennaStations(ids)
			    );
            }
            else {
                vector<String> names(1, String(ants.toString()));
			    return _vectorStringToStdVectorString(
				    *_msmd->getAntennaStations(names).rbegin()
			    );
            }
		}
		else if (type == variant::STRINGVEC) {
            if (filterObsID) {
                auto antIDs = _vectorIntToVectorUInt(antennaids(ants, "0m", "1pc", obsid));
                return _vectorStringToStdVectorString(_msmd->getAntennaStations(antIDs));
            }
            else {
                auto inputAnts = _vectorStdStringToVectorString(ants.toStringVec());
                auto allStations = _msmd->getAntennaStations(inputAnts);
                vector<string> lastStations;
                uInt i = 0;
                for (const auto& myStations: allStations) {
                    auto n = myStations.size();
                    if (n > 1) {
                        auto antIDs = _setUIntToVectorUInt(_msmd->getAntennaIDs(inputAnts[i]));
                        *_log << LogIO::WARN << "Antenna " << inputAnts[i] << " is listed "
                            << n << " times in the ANTENNA table. Its associated stations are "
                            << myStations << " which are associated with antenna IDs " << antIDs
                            << ",respectively. Only the station associated with the last ID ("
                            << *antIDs.rbegin() << ") will be returned." << LogIO::POST;
                    }
                    // This works in the n > 1 case because myStations is a vector ordered by
                    // antenna ID, not a set
                    lastStations.push_back(*myStations.rbegin());
                    ++i;
                }
                return lastStations;
            }
		}
		else if (type == variant::BOOLVEC) {
            vector<uInt> myAnts;
            if (filterObsID) {
                myAnts = _setUIntToVectorUInt(filteredAnts);
            }
			return _vectorStringToStdVectorString(_msmd->getAntennaStations(myAnts));
		}
		else {
			throw AipsError(
				"Unsupported type (" + ants.typeString() + ") for ants."
			);
		}
    )
	return vector<string>();
}

vector<int> msmetadata::antennasforscan(int scan, int obsid, int arrayid) {
	_FUNC(
		auto scanKeys = _getScanKeys(scan, obsid, arrayid);
		std::set<int> myres;
		for (auto scanKey: scanKeys) {
			auto t = _msmd->getAntennasForScan(scanKey);
			myres.insert(t.begin(), t.end());
		}
		return vector<int>(myres.begin(), myres.end());
	)
	return vector<int>();
}

variant* msmetadata::bandwidths(const variant& spws) {
	_FUNC(
		variant::TYPE type = spws.type();
		ThrowIf(
			type != variant::INT && type != variant::INTVEC
			&& type != variant::BOOLVEC,
			"Unsupported type for spws (" + spws.typeString()
			+ "), must be either an int or an array of ints"
		)

		int mymax = 0;
		if (type == variant::INT) {
			mymax = spws.toInt();
		}
		else if (type == variant::INTVEC) {
			Vector<Int> tmp(spws.toIntVec());
			ThrowIf(
				min(tmp) < 0,
				"When specified as an array, no element of spws may be < 0"
			);
			mymax = max(tmp);
		}
		_checkSpwId(mymax, false);
		vector<Double> bws = _msmd->getBandWidths();
		if (type == variant::BOOLVEC) {
			return new variant(bws);
		}
		if (type == variant::INT) {
			if (spws.toInt() < 0) {
				return new variant(bws);
			}
			else {
				return new variant(bws[spws.toInt()]);
			}
		}
		else {
			vector<Double> ret;
			for(auto id : spws.toIntVec()) {
				ret.push_back(bws[id]);
			}
			return new variant(ret);
		}
	)
	return 0;
}

int msmetadata::baseband(int spw) {
	_FUNC (
		_checkSpwId(spw, true);
		return _msmd->getBBCNos()[spw];
	)
	return 0;
}

variant* msmetadata::baselines() {
	_FUNC (
		Matrix<Bool> baselines = _msmd->getUniqueBaselines();
		vector<bool> values = baselines.tovector();
		vector<int> shape = baselines.shape().asStdVector();
		return new variant(values, shape);
	)
	return 0;
}

vector<int> msmetadata::chanavgspws() {
	_FUNC (
			/*
		*_log << LogIO::WARN << "This method is deprecated and will be removed. "
			<< "Use almaspws(chavg=true) instead." << LogIO::POST;
			*/
		return _setUIntToVectorInt(_msmd->getChannelAvgSpw());
	)
	return vector<int>();
}

vector<double> msmetadata::chaneffbws(int spw, const string& unit, bool asvel) {
    _FUNC(
        _checkSpwId(spw, true);
        string myunit = unit;
        if (myunit.empty()) {
            myunit = asvel ? "km/s" : "Hz";
        }
        vector<QVector<Double> > qvd = _msmd->getChanEffectiveBWs(asvel);
        return qvd[spw].getValue(Unit(myunit), true).tovector();
    )
    return vector<double>();
}

vector<double> msmetadata::chanfreqs(int spw, const string& unit) {
	_FUNC (
		_checkSpwId(spw, true);
		return _msmd->getChanFreqs()[spw].getValue(Unit(unit)).tovector();
	)
	return vector<double>();
}

vector<double> msmetadata::chanres(int spw, const string& unit, bool asvel) {
    _FUNC(
        _checkSpwId(spw, true);
        string myunit = unit;
        if (myunit.empty()) {
            myunit = asvel ? "km/s" : "Hz";
        }
        vector<QVector<Double> > qvd = _msmd->getChanResolutions(asvel);
        return qvd[spw].getValue(Unit(myunit), true).tovector();
    )
    return vector<double>();
}

vector<double> msmetadata::chanwidths(int spw, const string& unit) {
	_FUNC (
		_checkSpwId(spw, true);
		return _msmd->getChanWidths()[spw].getValue(Unit(unit), true).tovector();
	)
	return vector<double>();
}

bool msmetadata::close() {
	_FUNC2(
		_msmd.reset(0);
		_ms.reset(0);
		return true;
	)
	return false;
}

variant* msmetadata::corrprodsforpol(int polid) {
	_FUNC(
		_checkPolId(polid, true);
		Array<Int> prods = _msmd->getCorrProducts()[polid];
		return new variant(
			vector<int>(prods.begin(), prods.end()),
			prods.shape().asStdVector()
		);
	)
	return NULL;
}

vector<int> msmetadata::corrtypesforpol(int polid) {
	_FUNC(
		_checkPolId(polid, true);
		return _msmd->getCorrTypes()[polid];
	)
	return vector<int>();
}

vector<int> msmetadata::datadescids(int spw, int pol) {
	_FUNC(
		_checkSpwId(spw, false);
		_checkPolId(pol, false);
		auto mymap = _msmd->getSpwIDPolIDToDataDescIDMap();
		vector<int> ddids;
		for(auto iter : mymap) {
			uInt myspw = iter.first.first;
			uInt mypol = iter.first.second;
			uInt ddid = iter.second;
			if (
				(spw < 0 || (Int)myspw == spw)
				&& (pol < 0 || (Int)mypol == pol)
			) {
				ddids.push_back(ddid);
			}
		}
		std::sort (ddids.begin(),ddids.end());
		return ddids;
	)
	return vector<int>();
}

bool msmetadata::done() {
	_FUNC2(
		_msmd.reset(0);
		_ms.reset(0);
		return true;
	)
	return false;
}

record* msmetadata::effexposuretime() {
	_FUNC(
		return fromRecord(
			QuantumHolder(_msmd->getEffectiveTotalExposureTime()).toRecord()
		);
	)
	return 0;
}

record* msmetadata::exposuretime(
	int scan, int spwid, int polid, int obsid, int arrayid
) {
	_FUNC(
		_checkSpwId(spwid, true);
		_checkObsId(obsid, true);
		if (polid >= 0) {
			_checkPolId(polid, true);
		}
		else {
			std::set<uInt> polids = _msmd->getPolarizationIDs(
				(uInt)obsid, arrayid, scan, (uInt)spwid
			);
			ThrowIf(
				polids.empty(),
				"This dataset has no records for the specified scan and spwid"
			);
			if (polids.size() == 1) {
				polid = *(polids.begin());
			}
			else {
				*_log << LogIO::WARN	<< "This scan and spwID has multiple "
					<< " polarization IDs which are " << polids
					<< ". You must specify one of those."
					<< LogIO::POST;
				return nullptr;
			}
		}
		std::map<std::pair<uInt COMMA uInt> COMMA uInt> ddidMap
		    = _msmd->getSpwIDPolIDToDataDescIDMap();
		std::pair<uInt COMMA uInt> mykey;
		mykey.first = spwid;
		mykey.second = polid;
		ThrowIf(
			ddidMap.find(mykey) == ddidMap.end(),
			"MS has no data description ID for spectral window ID "
			+ String::toString(spwid) + " and polarization ID "
			+ String::toString(polid)
		);
		uInt dataDescID = ddidMap[mykey];
		map<ScanKey COMMA MSMetaData::FirstExposureTimeMap> map2
		    = _msmd->getScanToFirstExposureTimeMap(false);
		ScanKey scanKey;
		scanKey.arrayID = arrayid;
		scanKey.obsID = obsid;
		scanKey.scan = scan;
		ThrowIf(
			map2.find(scanKey) == map2.end()
			|| map2[scanKey].find(dataDescID) == map2[scanKey].end(),
			"MS has no records for scan number " + String::toString(scan)
			+ ", spectral window ID " + String::toString(spwid)
		    + ", and polarization ID " + String::toString(polid)
		);
		return fromRecord(
		    QuantumHolder(map2[scanKey][dataDescID].second).toRecord()
		);
	)
	return nullptr;
}

vector<int> msmetadata::fdmspws() {
	_FUNC(
			/*
		*_log << LogIO::WARN << __FUNCTION__ << " is deprecated and will be removed. "
			<< "Use almaspws(fdm=true) instead." << LogIO::POST;
			*/
		return _setUIntToVectorInt(_msmd->getFDMSpw());
	)
	return vector<int>();
}


vector<string> msmetadata::fieldnames() {
	_FUNC(
		return _vectorStringToStdVectorString(_msmd->getFieldNames());
	)
	return vector<string>();
}

variant* msmetadata::fieldsforintent(
    const string& intent, const bool asnames
) {
    _FUNC(
        std::set<Int> ids; 
        auto expand = intent.find('*') != std::string::npos;
        if (intent == "*" && _msmd->getIntents().empty()) {
            auto nFields = _msmd->nFields();
            for (uInt i=0; i<nFields; ++i) {
                ids.insert(i);
            }
        }
        else if (expand) {
            auto mymap = _msmd->getIntentToFieldsMap();
            ids = _idsFromExpansion(mymap, intent);
        }
        else {
            ids = _msmd->getFieldsForIntent(intent);
        }
        variant *x;
        if (ids.size() == 0) { 
            *_log << LogIO::WARN << "No intent " << (expand ? "matching '" : "'") 
                << intent << "' exists in this dataset." << LogIO::POST;
            x = asnames
                ? new variant(vector<string>(0))
                : new variant(vector<int>(0));
        }
        else {
            x = asnames
                ? new variant(_fieldNames(ids))
                : new variant(_setIntToVectorInt(ids));
        }
        return x;
    )
    return nullptr;
}

vector<int> msmetadata::fieldsforname(const string& name) {
	_FUNC(
		if (name.empty()) {
			return _setIntToVectorInt(_msmd->getUniqueFiedIDs());
		}
		return _setIntToVectorInt(_msmd->getFieldIDsForField(name));
	)
	return vector<int>();
}

variant* msmetadata::fieldsforscan(int scan, bool asnames, int obsid, int arrayid) {
	_FUNC(
		ThrowIf(
			scan < 0,
			"Scan number must be nonnegative."
		);
		auto scanKeys = _getScanKeys(scan, obsid, arrayid);
		std::set<int> ids;
		for (auto scanKey: scanKeys) {
			auto t = _msmd->getFieldsForScan(scanKey);
			ids.insert(t.begin(), t.end());
		}
		if (asnames) {
			return new variant(_fieldNames(ids));
		}
		else {
			return new variant(
				_setIntToVectorInt(ids)
			);
		}
	)
	return nullptr;
}

variant* msmetadata::fieldsforscans(
	const vector<int>& scans, const bool asnames,
	int obsid, int arrayid, bool asmap
) {
	_FUNC(
		ThrowIf(
			! scans.empty()
            && *std::min_element (scans.begin(), scans.end()) < 0,
			"All scan numbers must be non-negative"
		);
        if (asmap) {
            ThrowIf(
                obsid < 0 || arrayid < 0,
                "When asmap is true, both obsid and arrayid must be nonnegative"
            );
            ArrayKey ak;
            ak.obsID = obsid;
            ak.arrayID = arrayid;
            auto subScanKeys = _msmd->getSubScanKeys(ak);
            std::map<int COMMA set<int>> mymap;
            for (const auto& ss : subScanKeys) {
                mymap[ss.scan].insert(ss.fieldID);
            }
            if (! scans.empty()) {
                auto tmpmap = mymap;
                for (auto& p : mymap) {
                    if (std::find(scans.begin(), scans.end(), p.first) == scans.end()) {
                        tmpmap.erase(p.first);
                    }
                }
                mymap = tmpmap;
            }
            record ret;
            for (const auto& p : mymap) {
                if (asnames) {
                    ret.insert(to_string(p.first), _fieldNames(p.second));
                }
                else {
                    ret.insert(to_string(p.first), _setIntToVectorInt(p.second));
                }
            }
            return new variant(ret);
        }
        else {
		    auto scanKeys = _getScanKeys(scans, obsid, arrayid);
		    std::set<int> ids;
		    for (auto scanKey: scanKeys) {
			    auto t = _msmd->getFieldsForScan(scanKey);
			    ids.insert(t.begin(), t.end());
		    }
		    if (asnames) {
			    return new variant(_fieldNames(ids));
		    }
		    else {
			    return new variant(
				    _setIntToVectorInt(ids)
			    );
            }
		}
	)
	return nullptr;
}

variant* msmetadata::fieldsforsource(int sourceID, bool asnames) {
	_FUNC(
		if (asnames) {
			auto res = _msmd->getFieldNamesForSourceMap();
			if (res.find(sourceID) == res.end()) {
				return new variant(vector<string>());
			}
			else {
				return new variant(
					_setStringToVectorString(res[sourceID])
				);
			}
		}
		else {
			auto res = _msmd->getFieldsForSourceMap();
			if (res.find(sourceID) == res.end()) {
				return new variant(vector<int>());
			}
			else {
				return new variant(
					_setIntToVectorInt(res[sourceID])
				);
			}
		}
	)
	return nullptr;
}

record* msmetadata::fieldsforsources(bool asnames) {
    _FUNC(
        auto_ptr<record> ret(new record());
        if (asnames) {
            auto mymap = _msmd->getFieldNamesForSourceMap();
            for (const auto& p: mymap) {
                ret->insert(
                    String::toString(p.first), _setStringToVectorString(p.second)
                );
            }
        }
        else {
            auto mymap = _msmd->getFieldsForSourceMap();
            for (const auto& p: mymap) {
                ret->insert(
                    String::toString(p.first), _setIntToVectorInt(p.second)
                );
            }
        }
        return ret.release();
    )
    return nullptr;
}

variant* msmetadata::fieldsforspw(const int spw, const bool asnames) {
	_FUNC(
		_checkSpwId(spw, true);
		if (asnames) {
			return new variant(
				_setStringToVectorString(_msmd->getFieldNamesForSpw(spw))
			);
		}
		else {
			return new variant(
				_setIntToVectorInt(_msmd->getFieldIDsForSpw(spw))
			);
		}
	)
	return 0;
}

vector<int> msmetadata::fieldsfortimes(const double center, const double tol) {
	_FUNC(
		return _setIntToVectorInt(_msmd->getFieldsForTimes(center, tol));
	)
	return vector<int>();
}

vector<string> msmetadata::intents() {
	_FUNC(
		return _setStringToVectorString(_msmd->getIntents());
	)
	return vector<string>();
}

vector<string> msmetadata::intentsforfield(const variant& field) {
	_FUNC(
		Int id = -1;
		switch (field.type()) {
		case variant::STRING:
			id = *(_msmd->getFieldIDsForField(field.toString()).begin());
			break;
		case variant::INT:
			id = field.toInt();
			break;
		default:
			*_log << "Unsupported type for field which must be "
				<< "a nonnegative int or string." << LogIO::EXCEPTION;
		}
		if (id < 0) {
			throw AipsError("field must be nonnegative if an int.");
		}
		return _setStringToVectorString(_msmd->getIntentsForField(id));
	)
	return vector<string>();
}

vector<string> msmetadata::intentsforscan(int scan, int obsid, int arrayid) {
	_FUNC(
		if (scan < 0) {
			throw AipsError("Scan number must be nonnegative.");
		}
		auto scanKeys = _getScanKeys(scan, obsid, arrayid);
		std::set<String> intents;
		for (auto scanKey: scanKeys) {
			auto t = _msmd->getIntentsForScan(scanKey);
			intents.insert(t.begin(), t.end());
		}
		return _setStringToVectorString(intents);
	)
	return vector<string>();
}

vector<string> msmetadata::intentsforspw(int spw) {
	_FUNC(
		_checkSpwId(spw, true);
		return _setStringToVectorString(_msmd->getIntentsForSpw(spw));
	)
	return vector<string>();
}

double msmetadata::meanfreq(int spw, const string& unit) {
	_FUNC (
		_checkSpwId(spw, true);
		return _msmd->getMeanFreqs()[spw].getValue(Unit(unit));
	)
	return 0;
}

vector<string> msmetadata::namesforfields(const variant& fieldids) {
	_FUNC(
		variant::TYPE myType = fieldids.type();
		vector<uInt> fieldIDs;
		if (myType == variant::INT) {
			Int id = fieldids.toInt();
			if (id < 0) {
				throw AipsError("Field ID must be nonnegative.");
			}
			fieldIDs.push_back(id);
		}
		else if (myType == variant::INTVEC) {
			vector<Int> kk = fieldids.toIntVec();
			if (min(Vector<Int>(kk)) < 0 ) {
				throw AipsError("All field IDs must be nonnegative.");
			}
			fieldIDs = _vectorIntToVectorUInt(kk);
		}
		else if (fieldids.size() != 0) {
			throw AipsError(
				"Unsupported type for fieldids. It must be a nonnegative integer or nonnegative integer array"
			);
		}
		return _vectorStringToStdVectorString(
			_msmd->getFieldNamesForFieldIDs(fieldIDs)
		);
	)
	return vector<string>();
}

vector<string> msmetadata::namesforspws(const variant& spwids) {
	_FUNC(
		variant::TYPE myType = spwids.type();
		vector<uInt> spwIDs;
		if (myType == variant::INT) {
			Int id = spwids.toInt();
			ThrowIf(id < 0, "Spectral window ID must be nonnegative.");
			ThrowIf(
				id >= (Int)_msmd->nSpw(true),
				"Spectral window ID must be less than total number of spws"
			);
			spwIDs.push_back(id);
		}
		else if (myType == variant::INTVEC) {
			vector<Int> kk = spwids.toIntVec();
			Vector<Int> xx(kk);
			ThrowIf(
				min(xx) < 0,
				"All spectral window IDs must be nonnegative."
			);
			ThrowIf(
				max(xx) >= (Int)_msmd->nSpw(true),
				"All spectral window IDs must be less than "
				"the total number of spws"
			);
			spwIDs = _vectorIntToVectorUInt(kk);
		}
		else if (
			(myType == variant::STRING && spwids.toString().empty())
			|| myType == variant::BOOLVEC
		) {
			return _vectorStringToStdVectorString(
				_msmd->getSpwNames()
			);
		}
		else if (spwids.size() != 0) {
			ThrowCc(
				"Unsupported type for spwids. It must be a "
				"nonnegative integer or nonnegative integer array"
			);
		}
		vector<String> allNames = _msmd->getSpwNames();
		vector<String> names;
		for(auto i : spwIDs) {
			names.push_back(allNames[i]);
		}
		return _vectorStringToStdVectorString(names);
	)
	return vector<string>();
}

string msmetadata::name() {
	_FUNC(
		return _ms->tableName();
	)
	return "";
}

int msmetadata::nantennas() {
	_FUNC(
		return _msmd->nAntennas();
	)
	return 0;
}

int msmetadata::narrays() {
	_FUNC(
		return _msmd->nArrays();
	)
	return 0;
}

int msmetadata::nbaselines(bool doAC) {
	_FUNC(
		return _msmd->nBaselines(doAC);
	)
	return 0;
}

int msmetadata::nchan(int spw) {
	_FUNC (
		_checkSpwId(spw, true);
		return _msmd->nChans()[spw];
	)
	return 0;
}

variant* msmetadata::ncorrforpol(int polid) {
	_FUNC(
		_checkPolId(polid, false);
		vector<Int> ncorr = _msmd->getNumCorrs();
		if (polid < 0) {
			return new variant(ncorr);
		}
		return new variant(ncorr[polid]);
	)
	return NULL;
}

int msmetadata::nfields() {
	_FUNC(
		return _msmd->nFields();
	)
	return 0;

}

int msmetadata::nobservations() {
	_FUNC(
		return _msmd->nObservations();
	)
	return 0;
}

int msmetadata::nscans() {
	_FUNC(
		return _msmd->nScans();
	)
	return 0;
}

int msmetadata::nsources() {
	_FUNC(
		return _msmd->nUniqueSourceIDsFromSourceTable();
	)
	return 0;
}

int msmetadata::nspw(bool includewvr) {
	_FUNC(
		return _msmd->nSpw(includewvr);
	)
	return 0;
}

int msmetadata::nstates() {
	_FUNC(
		return _msmd->nStates();
	)
	return 0;
}

double msmetadata::nrows(const bool ac, const bool flagged) {
	_FUNC(
		if (ac) {
			return flagged ? _msmd->nRows() : _msmd->nUnflaggedRows();
		}
		else {
			return flagged ? _msmd->nRows(MSMetaData::CROSS)
				: _msmd->nUnflaggedRows(MSMetaData::CROSS);
		}
	)
	return 0;
}

vector<string> msmetadata::observers() {
	_FUNC(
		return _vectorStringToStdVectorString(_msmd->getObservers());
	)
	return vector<string>();
}

vector<string> msmetadata::observatorynames() {
	_FUNC(
		return _vectorStringToStdVectorString(_msmd->getObservatoryNames());
	)
	return vector<string>();
}

record* msmetadata::observatoryposition(const int which) {
	_FUNC(
		MeasureHolder out(_msmd->getObservatoryPosition(which));
		Record outRec;
		String error;
		if (!out.toRecord(error, outRec)) {
			error += "Failed to generate position.\n";
			*_log << LogIO::SEVERE << error << LogIO::POST;
			return 0;
		}
		return fromRecord(outRec);
	)
	return 0;
}

::casac::record* msmetadata::phasecenter(const int fieldid, const ::casac::record& epoch){
	::casac::record *rval=0;
    _FUNC(   
    	PtrHolder<Record> ep(toRecord(epoch));
	  	Record outRec;
	  	MeasureHolder mh;
	  	String err;
	  	if (ep->nfields() == 0) {
	  		mh = MeasureHolder (_msmd->phaseDirFromFieldIDAndTime(fieldid));
	  	}
	  	else {
	  		MeasureHolder ephold;
	  		ThrowIf(! ephold.fromRecord(err, *ep), "Epoch cannot be converted \n" + err);
	  		ThrowIf(! ephold.isMEpoch(), "Epoch parameter is not an MEpoch  \n");
	  		mh = MeasureHolder (_msmd->phaseDirFromFieldIDAndTime(fieldid, ephold.asMEpoch()));
	  	}
	  	ThrowIf(! mh.toRecord(err, outRec), "Could not convert phasecenter \n" + err);
	  	rval = fromRecord(outRec);
    )
    return rval;
}

record* msmetadata::pointingdirection(int rowid, bool const interpolate, int const initialrow) {
	_FUNC(
		Int ant1 COMMA ant2;
		Double time;
		
		std::pair<casacore::MDirection COMMA casacore::MDirection> pDirs = _msmd->getPointingDirection(
			ant1, ant2, time, rowid, interpolate, initialrow
		);
		MeasureHolder m1(pDirs.first);
		MeasureHolder m2(pDirs.second);
		Record ret;
		ret.define("time", time);
		Record ant1Rec COMMA ant2Rec COMMA m1rec COMMA m2rec;
		ant1Rec.define("id", ant1);
		m1.toRecord(m1rec);
		ant1Rec.defineRecord("pointingdirection", m1rec);
		ret.defineRecord("antenna1", ant1Rec);
		ant2Rec.define("id", ant2);
		m2.toRecord(m2rec);
		ant2Rec.defineRecord("pointingdirection", m2rec);
		ret.defineRecord("antenna2", ant2Rec);
		return fromRecord(ret);
	);
	return NULL;
}

void msmetadata::_init(const casacore::MeasurementSet *const &ms, const float cachesize) {
    _msmd.reset(new MSMetaData(ms, cachesize));
    _msmd->setForceSubScanPropsToCache(true);
    _msmd->setShowProgress(true);
}

bool msmetadata::open(const string& msname, const float cachesize) {
	_FUNC2(

        _ms.reset(new MeasurementSet(msname));
        *_log << LogIO::NORMAL << "Performing internal consistency checks on "
            << msname << "..." << LogIO::POST;
        MSChecker msChecker(*_ms);
		msChecker.checkReferentialIntegrity();
		_init(_ms.get(), cachesize);
		return true;
	)
	return false;
}

variant* msmetadata::polidfordatadesc(int ddid) {
	_FUNC(
		vector<uInt> pols = _msmd->getDataDescIDToPolIDMap();
		if (ddid < 0) {
			return new variant(pols);
		}
		int nddids = pols.size();
		ThrowIf(ddid >= nddids, "ddid must be less than " + String::toString(nddids));
		return new variant(pols[ddid]);
	)
	return NULL;
}

vector<string> msmetadata::projects() {
	_FUNC(
		return _vectorStringToStdVectorString(_msmd->getProjects());
	)
	return vector<string>();
}

record* msmetadata::propermotions() {
	_FUNC(
		vector<std::pair<casacore::Quantity COMMA casacore::Quantity> > mu = _msmd->getProperMotions();
		Record rec;
		Record subrec;
		uInt n = mu.size();
		Vector<Record> v(2);
		for (uInt i=0; i<n; ++i) {
			QuantumHolder q0(mu[i].first);
			QuantumHolder q1(mu[i].second);
			q0.toRecord(v[0]);
			q1.toRecord(v[1]);
			subrec.defineRecord("longitude", v[0]);
			subrec.defineRecord("latitude", v[1]);
			rec.defineRecord(casacore::String::toString(i), subrec);
		}
		return fromRecord(rec);
	)
	return NULL;
}

record* msmetadata::refdir(
    const variant& field, const record& epoch
) {
    _FUNC(
        Int id;
        switch (field.type()) {
        case variant::STRING:
            id = *(_msmd->getFieldIDsForField(field.toString()).begin());
            break;
        case variant::INT:
            id = field.toInt();
            break;
        default:
            ThrowCc(
                "Unsupported type for field which must be "
                "a nonnegative int or string."
            );
        }
        unique_ptr<Record> ep(toRecord(epoch));
        Record outRec;
        MeasureHolder mh;
        String err; 
        if (ep->nfields() == 0) { 
            mh = MeasureHolder(_msmd->getReferenceDirection(id));
        }
        else {
            MeasureHolder ephold;
            ThrowIf(
                ! ephold.fromRecord(err, *ep),
                "Epoch cannot be converted \n" + err
            );
            ThrowIf(
                ! ephold.isMEpoch(), "Epoch parameter is not an MEpoch  \n"
            );
            mh = MeasureHolder(_msmd->getReferenceDirection(id, ephold.asMEpoch()));
        }
        ThrowIf(
            ! mh.toRecord(err, outRec),
            "Could not convert reference direction to Record \n" + err
        );
        return fromRecord(outRec);
    )    
    return nullptr;
}

record* msmetadata::reffreq(int spw) {
	_FUNC(
		_checkSpwId(spw, true);
		MeasureHolder freq(_msmd->getRefFreqs()[spw]);
		Record ret;
		freq.toRecord(ret);
		return fromRecord(ret);
	)
	return NULL;
}

variant* msmetadata::restfreqs(int sourceid, int spw) {
    _FUNC(
        _checkSpwId(spw, true);
        ThrowIf(
            sourceid < 0, "sourceid cannot be negative"
        );
        map<SourceKey COMMA SHARED_PTR<vector<MFrequency> > > mymap
            = _msmd->getRestFrequencies();
        SourceKey key;
        key.id = sourceid;
        key.spw = spw;
        ThrowIf (
            mymap.find(key) == mymap.end(),
            "SOURCE table does not contain a row with SOURCE_ID="
            + String::toString(sourceid) + " and SPECTRAL_WINDOW_ID="
            + String::toString(spw)
        );
        SHARED_PTR<vector<MFrequency> > ptr = mymap[key];
        if (ptr) {
            Record mr;
            Record r;
            uInt i = 0;
            for (const auto& freq : *ptr) {
                MeasureHolder mh(freq);
                mh.toRecord(mr);
                r.defineRecord(casacore::String::toString(i), mr);
                ++i;
            }
            return new variant(fromRecord(r));
        }
        else {
            return new variant(false);
        }
    )
    return nullptr;
}

vector<int> msmetadata::scannumbers(int obsid, int arrayid) {
	_FUNC(
		_checkObsId(obsid, false);
		_checkArrayId(arrayid, false);
		return _setIntToVectorInt(_msmd->getScanNumbers(obsid, arrayid));
	)
	return vector<int>();
}

vector<int> msmetadata::scansforfield(
	const variant& field, int obsid, int arrayid
) {
	_FUNC(
		_checkObsId(obsid, false);
		_checkArrayId(arrayid, false);
		switch (field.type()) {
		case variant::INT:
			return _setIntToVectorInt(_msmd->getScansForFieldID(field.toInt(), obsid, arrayid));
			break;
		case variant::STRING:
			return _setIntToVectorInt(_msmd->getScansForField(field.toString(), obsid, arrayid));
			break;
		default:
			throw AipsError("Unsupported type for field parameter.");
		}
	)
	return vector<int>();
}

record* msmetadata::scansforfields(int obsid, int arrayid) {
    _FUNC(
        _checkArrayId(arrayid, true);
        _checkObsId(obsid, true);
        auto fieldToScans = _msmd->getFieldToScansMap();
        std::auto_ptr<record> ret(new record());
        uInt n = fieldToScans.size();
        ArrayKey ak;
        ak.obsID = obsid;
        ak.arrayID = arrayid;
        for (uInt i=0; i<n; ++i) {
            auto scans = filter(fieldToScans[i], ak);
            std::set<uInt> scanNums;
            for (const auto& scan: scans) {
                scanNums.insert(scan.scan);
            }
            ret->insert(
                String::toString(i), _setUIntToVectorInt(scanNums)
            );
        }
        return ret.release();
    )
    return nullptr;
}

vector<int> msmetadata::scansforintent(const string& intent, int obsid, int arrayid) {
	_FUNC(
		_checkObsId(obsid, false);
		_checkArrayId(arrayid, false);
		Bool expand = intent.find('*') != std::string::npos;
		if (expand) {
			auto mymap = _msmd->getIntentToScansMap();
			auto scanKeys = _idsFromExpansion(mymap, intent);
			auto doAllObs = obsid < 0;
			auto doAllArrays = arrayid < 0;
			std::set<Int> myScanNumbers;
			for (const auto k : scanKeys) {
				if (
					(doAllObs || k.obsID == obsid)
					&& (doAllArrays || k.arrayID == arrayid)
				) {
					myScanNumbers.insert(k.scan);
				}
			}
			return _setIntToVectorInt(myScanNumbers);
		}
		else {
			return _setIntToVectorInt(_msmd->getScansForIntent(intent, obsid, arrayid));
		}
	)
	return vector<int>();
}

vector<int> msmetadata::scansforspw(const int spw, int obsid, int arrayid) {
	_FUNC(
		_checkSpwId(spw, true);
		_checkObsId(obsid, false);
		_checkArrayId(arrayid, false);
		return _setIntToVectorInt(_msmd->getScansForSpw(spw, obsid, arrayid));
	)
	return vector<int>();
}

record* msmetadata::scansforspws(int obsid, int arrayid) {
	_FUNC(
	    _checkArrayId(arrayid, true);
	    _checkObsId(obsid, true);
        auto spwToScans = _msmd->getSpwToScansMap();
		auto_ptr<record> ret(new record());
        uInt n = spwToScans.size();
   		ArrayKey ak;
	    ak.obsID = obsid;
	    ak.arrayID = arrayid;
        for (uInt i=0; i<n; ++i) {
		    auto scans = filter(spwToScans[i], ak);
            std::set<uInt> scanNums;
            for (const auto& scan: scans) {
                scanNums.insert(scan.scan);
            }
            ret->insert(
                String::toString(i), _setUIntToVectorInt(scanNums)
            );
		}
		return ret.release();
	)
	return nullptr;
}

vector<int> msmetadata::scansforstate(int state, int obsid, int arrayid) {
	_FUNC(
		ThrowIf(
			state < 0, "State must be nonnegative."
		);
		_checkObsId(obsid, false);
		_checkArrayId(arrayid, false);
		return _setIntToVectorInt(_msmd->getScansForState(state, obsid, arrayid));
	)
	return vector<int>();
}

vector<int> msmetadata::scansfortimes(
	double center, double tol, int obsid, int arrayid
) {
	_FUNC(
		_checkObsId(obsid, false);
		_checkArrayId(arrayid, false);
		return _setIntToVectorInt(
			_msmd->getScansForTimes(center, tol, obsid, arrayid)
		);
	)
	return vector<int>();
}

vector<string> msmetadata::schedule(int obsid) {
	_FUNC(
		_checkObsId(obsid, true);
		return _vectorStringToStdVectorString(_msmd->getSchedules()[obsid]);
	)
	return vector<string>();
}

int msmetadata::sideband(int spw) {
	_FUNC (
		_checkSpwId(spw, true);
        Int ret = _msmd->getNetSidebands()[spw] == 2 ? 1 : -1;
        return ret;
	)
	return 0;
}

record* msmetadata::sourcedirs() {
	_FUNC(
		std::vector<casacore::MDirection> mdirs = _msmd->getSourceDirections();
		uInt i = 0;
		vector<casacore::MDirection>::const_iterator iter = mdirs.begin();
		vector<casacore::MDirection>::const_iterator end = mdirs.end();
		Record r;
		Record mr;
		while (iter != end) {
			MeasureHolder mh(*iter);
			mh.toRecord(mr);
			r.defineRecord(casacore::String::toString(i), mr);
			++iter;
			++i;
		}
		return fromRecord(r);
	)
	return nullptr;
}

record* msmetadata::sourcetimes() {
	_FUNC(
		auto times = _msmd->getSourceTimes();
        Record r;
		uInt i = 0;
        Unit u = times->getUnit();
        for (const auto& t: times->getValue()) {
            casacore::Quantity q(t, u);
            QuantumHolder qh(q);
            r.defineRecord(String::toString(i), qh.toRecord());
            ++i;
        }
		return fromRecord(r);
	)
	return nullptr;
}

int msmetadata::sourceidforfield(int field) {
	_FUNC(
		_checkFieldId(field, true);
		return _msmd->getFieldTableSourceIDs()[field];
	)
	return 0;
}

vector<int> msmetadata::sourceidsfromsourcetable() {
	_FUNC(
		return _msmd->getSourceTableSourceIDs();
	)
	return vector<int>();
}

vector<string> msmetadata::sourcenames() {
	_FUNC(
		return _vectorStringToStdVectorString(_msmd->getSourceNames());
	)
	return vector<string>();
}

variant* msmetadata::spwfordatadesc(int ddid) {
	_FUNC(
		vector<uInt> spws = _msmd->getDataDescIDToSpwMap();
		if (ddid < 0) {
			return new variant(spws);
		}
		int nddids = spws.size();
		ThrowIf(ddid >= nddids, "ddid must be less than " + String::toString(nddids));
		return new variant(spws[ddid]);
	)
	return NULL;
}

variant* msmetadata::spwsforbaseband(int bb, const string& sqldmode) {
	_FUNC(
		String mode = sqldmode;
		mode.downcase();
		MSMetaData::SQLDSwitch sqld;
		if (mode.startsWith("i")) {
			sqld = MSMetaData::SQLD_INCLUDE;
		}
		else if (mode.startsWith("e")) {
			sqld = MSMetaData::SQLD_EXCLUDE;
		}
		else if (mode.startsWith("o")) {
			sqld = MSMetaData::SQLD_ONLY;
		}
		else {
			throw AipsError(
				"Unsupported sqldmode " + sqldmode
				+ ". Must be either i(nclude), e(xclude), or o(nly)."
			);
		}
		map<uInt COMMA set<uInt> > x = _msmd->getBBCNosToSpwMap(sqld);
		if (bb >= 0) {
			if (x.find(bb) == x.end()) {
				return new variant(vector<int>(0));
			}
			else {
				return new variant(_setUIntToVectorInt(x[bb]));
			}
		}
		else {
			record out;
			std::map<uInt COMMA std::set<uInt> >::const_iterator end = x.end();
			for (
				std::map<uInt COMMA std::set<uInt> >::const_iterator iter=x.begin();
				iter!=end; iter++
			) {
				vector<int> v = _setUIntToVectorInt(iter->second);
				out.insert(String::toString(iter->first), variant(v));
			}
			return new variant(out);
		}
	)
	return 0;
}

vector<int> msmetadata::spwsforintent(const string& intent) {
	_FUNC(
		Bool expand = intent.find('*') != std::string::npos;
		if (expand) {
			std::map<String COMMA std::set<uInt> > mymap = _msmd->getIntentToSpwsMap();
			std::set<Int> ids = _idsFromExpansion(mymap, intent);
			return _setIntToVectorInt(ids);
		}
		else {
			vector<int> x = _setUIntToVectorInt(_msmd->getSpwsForIntent(intent));
			if (x.size() == 0) {
				*_log << LogIO::WARN << "Intent " << intent
					<< " does not exist in this dataset." << LogIO::POST;
			}
			return x;
		}
	)
	return vector<int>();
}

vector<int> msmetadata::spwsforfield(const variant& field) {
	_FUNC(
		switch (field.type()) {
		case variant::INT:
			return _setUIntToVectorInt(_msmd->getSpwsForField(field.toInt()));
			break;
		case variant::STRING:
			return _setUIntToVectorInt(_msmd->getSpwsForField(field.toString()));
			break;
		default:
			throw AipsError("Unacceptable type for field parameter.");
		}
	)
	return vector<int>();
}

record* msmetadata::spwsforfields() {
	_FUNC(
        auto mymap = _msmd->getFieldsToSpwsMap();
        record* ret = new record();
        for (const auto& elem: mymap) {
            ret->insert(
                String::toString(elem.first),
                vector<uInt>(elem.second.begin(), elem.second.end())
            );
        }
        return ret;
	)
	return nullptr;
}

record* msmetadata::spwsfornames(const variant& names) {
    _FUNC(
        variant::TYPE myType = names.type();
        auto allNames = _msmd->getSpwNames();
        set<string> spwNames;
        if (
            (myType == variant::STRING && names.toString().empty())
            || myType == variant::BOOLVEC
        ) {
            spwNames.insert(allNames.begin(), allNames.end());
        }
        else if (myType == variant::STRING) {
            spwNames.insert(names.toString());
        }
        else if (myType == variant::STRINGVEC) {
            auto x = names.toStringVec();
            spwNames.insert(x.begin(), x.end());
        }
        else if (names.size() != 0) {
            ThrowCc(
                "Unsupported type for spwids. It must be a "
                "nonnegative integer or nonnegative integer array"
            );
        }
        unique_ptr<record> rec(new record());
        uInt id = 0;
        for (const auto& name : allNames) {
            if (
                spwNames.find(name) != spwNames.end()
            ) {
                rec->insert(name, vector<int>(1, id));
                spwNames.erase(name);
            }
            else if (rec->find(name) != rec->end()) {
                auto v = rec->at(name).asIntVec();
                v.push_back(id);
                (*rec)[name] = v;
            }
            ++id;
        }
        for (const auto& name : spwNames) {
            *_log << LogIO::WARN << "Specified spw named "
                << name << " is not present in MS"
                << LogIO::POST;
        }
        return rec.release();
    )
    return nullptr;
}

vector<int> msmetadata::spwsforscan(int scan, int obsid, int arrayid) {
	_FUNC(
		if (scan < 0) {
			throw AipsError("Scan must be nonnegative");
		}
		auto scanKeys = _getScanKeys(scan, obsid, arrayid);
		std::set<uInt> spws;
		for (const auto scanKey : scanKeys) {
			auto t = _msmd->getSpwsForScan(scanKey);
			spws.insert(t.begin(), t.end());
		}
		return _setUIntToVectorInt(spws);
	)
	return vector<int>();
}

record* msmetadata::spwsforscans(int obsid, int arrayid) {
	_FUNC(
	    _checkArrayId(arrayid, true);
	    _checkObsId(obsid, true);
        auto scanToSpws = _msmd->getScanToSpwsMap();
		auto_ptr<record> ret(new record());
		std::set<ScanKey> allScans;
		for (const auto& p : scanToSpws) {
		    allScans.insert(p.first);
		}
		ArrayKey ak;
	    ak.obsID = obsid;
	    ak.arrayID = arrayid;
		auto scans = filter(allScans, ak);
		for (const auto& sk : scans) {
		    ret->insert(
		        String::toString(sk.scan), _setUIntToVectorInt(scanToSpws[sk])
		    );
		}
		return ret.release();
	)
	return nullptr;
}

vector<int> msmetadata::statesforscan(int scan, int obsid, int arrayid) {
	_FUNC(
		ThrowIf(scan < 0, "Scan number must be nonnegative");
		_checkObsId(obsid, false);
		_checkArrayId(arrayid, false);
		return _setIntToVectorInt(_msmd->getStatesForScan(obsid, arrayid, scan));
	)
	return vector<int>();
}

record* msmetadata::statesforscans(int obsid, int arrayid) {
	_FUNC(
		_checkObsId(obsid, false);
		_checkArrayId(arrayid, false);
        std::map<ScanKey COMMA std::set<Int> > mymap = _msmd->getScanToStatesMap();
        ArrayKey arrayKey;
        arrayKey.obsID = obsid;
        arrayKey.arrayID = arrayid;
        auto scanKeys = _msmd->getScanKeys();
        auto filtered = filter(_msmd->getScanKeys(), arrayKey);
        unique_ptr<record> ret(new record());
        for (const auto& f: filtered) {
            ret->insert(to_string(f.scan), _setIntToVectorInt(mymap[f]));
        }
        return ret.release();
	)
	return nullptr;
}

record* msmetadata::timerangeforobs(int obsid) {
	_FUNC(
		_checkObsId(obsid, true);
		auto range = _msmd->getTimeRangesOfObservations()[obsid];
		MeasureHolder begin(range.first);
		MeasureHolder end(range.second);
		Record ret COMMA beginRec COMMA endRec;
		begin.toRecord(beginRec);
		end.toRecord(endRec);
		ret.defineRecord("begin", beginRec);
		ret.defineRecord("end", endRec);
		return fromRecord(ret);
	)
	return NULL;
}

vector<double> msmetadata::timesforfield(int field) {
	_FUNC(
		if (field < 0) {
			throw AipsError("Field ID must be nonnegative");
		}
		return _setDoubleToVectorDouble(_msmd->getTimesForField(field));
	)
	return vector<double>();
}

record* msmetadata::summary() {
	_FUNC(
		return fromRecord(_msmd->getSummary());
	)
	return NULL;
}

vector<double> msmetadata::timesforintent(const string& intent) {
	_FUNC(
		return _setDoubleToVectorDouble(_msmd->getTimesForIntent(intent));
	)
	return vector<double>();
}

variant* msmetadata::timesforscan(int scan, int obsid, int arrayid, bool perspw) {
	_FUNC(
		if (scan < 0) {
			throw AipsError("Scan number must be nonnegative");
		}
		auto scanKeys = _getScanKeys(scan, obsid, arrayid);
		if (perspw) {
			std::map<uInt COMMA std::set<Double> > spwToTimes;
			Record ret;
			for (const auto scanKey : scanKeys) {
				auto mymap = _msmd->getSpwToTimesForScan(scanKey);
				for (const auto& kv : mymap) {
					auto spw = kv.first;
					auto times = kv.second;
					if (spwToTimes.find(spw) == spwToTimes.end()) {
						spwToTimes[spw] = times;
					}
					else {
						spwToTimes[spw].insert(times.begin(), times.end());
					}
				}
			}
			for (const auto& kv : spwToTimes) {
				ret.define(
					casacore::String::toString(kv.first) COMMA
					Vector<Double>(_setDoubleToVectorDouble(kv.second)
					)
				);
			}
			SHARED_PTR<record> rec(fromRecord(ret));
			return new variant(*rec);
		}
		else {
			std::set<Double> times;
			for (const auto& scanKey : scanKeys) {
				auto t = _msmd->getTimesForScan(scanKey);
				times.insert(t.begin(), t.end());
			}
			return new variant(_setDoubleToVectorDouble(times));
		}
	)
	return nullptr;
}

vector<double> msmetadata::timesforscans(const variant& scans, int obsid, int arrayid) {
	_FUNC(
	    vector<int> myscans;
	    if (scans.type() == variant::INT) {
	        myscans.push_back(scans.toInt());
	    }
	    else if (scans.type() == variant::INTVEC) {
	        myscans = scans.toIntVec();
	    }
	    else {
	        ThrowCc("scans must either be an int or an array of ints");
	    }
	    ThrowIf(
	        *std::min_element(myscans.begin(), myscans.end()) < 0,
	        "All scan numbers must be nonnegative"
	    );
		auto scanKeys = _getScanKeys(myscans, obsid, arrayid);
		return _setDoubleToVectorDouble(_msmd->getTimesForScans(scanKeys));
	)
	return vector<double>();
}

vector<int> msmetadata::tdmspws() {
	_FUNC(
			/*
		*_log << LogIO::WARN << __FUNCTION__ << " is deprecated and will be removed. "
			<< "Use almaspws(tdm=true) instead." << LogIO::POST;
			*/
		return _setUIntToVectorInt(_msmd->getTDMSpw());
	)
	return vector<int>();
}

variant* msmetadata::timesforspws(const variant& spws) {
    _FUNC(
        auto type = spws.type();
        std::vector<Int> myspws;
        Int nspw = _msmd->nSpw(True);
        if (type == variant::BOOLVEC) {
            myspws.resize(nspw);
            std::iota(myspws.begin(), myspws.end(), 0);
        }
        else if (type == variant::INT) {
            auto spw = spws.toInt();
            if (spw >= 0) {
                myspws.push_back(spw);
            }
            else {
                std::iota(myspws.begin(), myspws.end(), 0);
            }
        }
        else if(type == variant::INTVEC) {
            myspws = spws.toIntVec();
            Vector<Int> test(myspws);
            ThrowIf(min(test) < 0, "If spws is an array, all values must be non-negative");
            ThrowIf(
                max(test) >= nspw,
                "Values in spws array must all be less than the number "
                "of spectral windows, which is " + String::toString(nspw)
            );
        }
        else {
            ThrowCc("Unsupported type for spws");
        }
        variant ret;
        auto vec = _msmd->getTimesForSpws(True);
        if (myspws.size() == 1) {
            auto spw = myspws[0];
            return new variant(vector<Double>(vec[spw].begin(), vec[spw].end()));
        }
        else {
            record x;
            for (const auto& spw: myspws) {
                vector<Double> times(vec[spw].begin(), vec[spw].end());
                x[String::toString(spw)] = times;
            }
            return new variant(x);
        }
    )
    return nullptr;
}

variant* msmetadata::transitions(int sourceid, int spw) {
    _FUNC(
        _checkSpwId(spw, true);
        ThrowIf(
            sourceid < 0, "sourceid cannot be negative"
        );
        map<SourceKey COMMA SHARED_PTR<vector<String> > > mymap
            = _msmd->getTransitions();
        SourceKey key;
        key.id = sourceid;
        key.spw = spw;
        ThrowIf (
            mymap.find(key) == mymap.end(),
            "SOURCE table does not contain a row with SOURCE_ID="
            + String::toString(sourceid) + " and SPECTRAL_WINDOW_ID="
            + String::toString(spw)
        );
        SHARED_PTR<vector<String> > ptr = mymap[key];
        if (ptr) {
        	vector<string> v = _vectorStringToStdVectorString(*ptr);
            return new variant(v);
        }
        else {
            return new variant(false);
        }
    )
    return nullptr;
}


vector<int> msmetadata::wvrspws(bool complement) {
	_FUNC(
			/*
		*_log << LogIO::WARN << __FUNCTION__ << " is deprecated and will be removed. "
			<< "Use almaspws(tdm=true) instead." << LogIO::POST;
			*/
		vector<int> wvrs = _setUIntToVectorInt(_msmd->getWVRSpw());
		if (complement) {
			vector<int> nonwvrs(
				counting_iterator<int>(0),
				counting_iterator<int>(_msmd->nSpw(true)));
			vector<int>::iterator begin = nonwvrs.begin();
			for_each(wvrs.rbegin(), wvrs.rend(), [&](int spw){ nonwvrs.erase(begin + spw); });
			return nonwvrs;
		}
		else {
			return wvrs;
		}
	)
	return vector<int>();
}

msmetadata::msmetadata(
	const MeasurementSet *const &ms, const float cachesize
) : _msmd(), _ms(), _log(new LogIO()) {
	_init(ms, cachesize);
}

vector<string> msmetadata::_fieldNames(const set<int>& ids) {
	if (*min_element(ids.begin(), ids.end()) < 0) {
		throw AipsError("All provided IDs must be greater than 0");
	}
	return _vectorStringToStdVectorString(
		_msmd->getFieldNamesForFieldIDs(vector<uInt>(ids.begin(), ids.end()))
	);
}

bool msmetadata::_isAttached(const bool throwExceptionIfNotAttached) const {
	if (_msmd.get() != 0) {
		return true;
	}
	else if (throwExceptionIfNotAttached) {
		throw AipsError("Tool is not attached to an MS. Use open()");
	}
	return false;
}

void msmetadata::_handleException(const AipsError& x) const {
	*_log << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
	RETHROW(x);
}

std::vector<double> msmetadata::_setDoubleToVectorDouble(
	const std::set<casacore::Double>& inset
) {
	vector<double> output;
	std::copy(inset.begin(), inset.end(), std::back_inserter(output));
	return output;
}

std::vector<std::string> msmetadata::_setStringToVectorString(
	const std::set<casacore::String>& inset
) {
	vector<string> output;
	std::copy(inset.begin(), inset.end(), std::back_inserter(output));
    return output;
}

std::vector<int> msmetadata::_setUIntToVectorInt(const std::set<casacore::uInt>& inset) {
	vector<int> output;
    std::copy(inset.begin(), inset.end(), std::back_inserter(output));
	return output;
}

std::vector<uInt> msmetadata::_setUIntToVectorUInt(const std::set<casacore::uInt>& inset) {
	vector<uInt> output;
    std::copy(inset.begin(), inset.end(), std::back_inserter(output));
	return output;
}

std::vector<int> msmetadata::_setIntToVectorInt(const std::set<casacore::Int>& inset) {
	vector<int> output;
    std::copy(inset.begin(), inset.end(), std::back_inserter(output));
	return output;
}

std::vector<std::string> msmetadata::_vectorStringToStdVectorString(const std::vector<casacore::String>& inset) {
	vector<string> output;
	std::copy(inset.begin(), inset.end(), std::back_inserter(output));
	return output;
}

std::vector<int> msmetadata::_vectorUIntToVectorInt(const std::vector<uInt>& inset) {
	vector<int> output;
	std::copy(inset.begin(), inset.end(), std::back_inserter(output));
	return output;
}

std::vector<uint> msmetadata::_vectorIntToVectorUInt(const std::vector<Int>& inset) {
	vector<uint> output;
	std::copy(inset.begin(), inset.end(), std::back_inserter(output));
	return output;
}

void msmetadata::_checkAntennaId(int id, bool throwIfNegative) const {
	ThrowIf(
		id >= (int)_msmd->nAntennas() || (throwIfNegative && id < 0),
		"Antenna ID " + String::toString(id)
		+ " out of range, must be less than "
		+ String::toString((int)_msmd->nAntennas())
	);
}

void msmetadata::_checkArrayId(int id, bool throwIfNegative) const {
	ThrowIf(
		id >= (int)_msmd->nArrays() || (throwIfNegative && id < 0),
		"Array ID " + String::toString(id)
		+ " out of range, must be less than "
		+ String::toString((int)_msmd->nArrays())
	);
}

void msmetadata::_checkFieldId(int id, bool throwIfNegative) const {
	ThrowIf(
		id >= (int)_msmd->nFields() || (throwIfNegative && id < 0),
		"Antenna ID " + String::toString(id)
		+ " out of range, must be less than "
		+ String::toString((int)_msmd->nFields())
	);
}

void msmetadata::_checkObsId(int id, bool throwIfNegative) const {
	ThrowIf(
		id >= (int)_msmd->nObservations() || (throwIfNegative && id < 0),
		"Observation ID " + String::toString(id)
		+ " out of range, must be less than "
		+ String::toString((int)_msmd->nObservations())
	);
}

void msmetadata::_checkSpwId(int id, bool throwIfNegative) const {
	ThrowIf(
		id >= (int)_msmd->nSpw(true) || (throwIfNegative && id < 0),
		"Spectral window ID " + String::toString(id)
		+ " out of range, must be "
		+ (throwIfNegative ? "nonnegative and " : "") + "less than "
		+ String::toString((int)_msmd->nSpw(true))
	);
}

void msmetadata::_checkPolId(int id, bool throwIfNegative) const {
	ThrowIf(
		id >= (int)_msmd->nPol(),
		"Polarization ID " + String::toString(id)
		+ " out of range, must be less than "
		+ String::toString((int)_msmd->nPol())
	);
	ThrowIf(throwIfNegative && id < 0, "Polarization ID cannot be negative");
}

std::set<ScanKey> msmetadata::_getScanKeys(int scan, int obsid, int arrayid) const {
	_checkObsId(obsid, false);
	_checkArrayId(arrayid, false);
	if (obsid >= 0 && arrayid >= 0) {
		std::set<ScanKey> scanKey;
		ScanKey x;
		x.scan = scan;
		x.obsID = obsid;
		x.arrayID = arrayid;
		scanKey.insert(x);
		return scanKey;
	}
	else {
		ArrayKey ak;
		ak.obsID = obsid;
		ak.arrayID = arrayid;
		auto scanKeys = _msmd->getScanKeys(ak);
		std::set<ScanKey> myKeys;
		for (const auto k : scanKeys) {
			if (k.scan == scan) {
				myKeys.insert(k);
			}
		}
		ThrowIf(myKeys.empty(), "No matching scans found");
		return myKeys;
	}
}

std::set<ScanKey> msmetadata::_getScanKeys(
	const vector<int>& scans, int obsid, int arrayid
) const {
	_checkObsId(obsid, false);
	_checkArrayId(arrayid, false);
    std::set<ScanKey> myKeys;
    if (scans.empty()) {
        ArrayKey akey;
        akey.obsID = obsid;
        akey.arrayID = arrayid;
        myKeys = _msmd->getScanKeys(akey);
    }
    else if (obsid >= 0 && arrayid >= 0) {
	    ScanKey x;
	    for (auto scan : scans) {
			x.scan = scan;
			myKeys.insert(x);
        }
	}
	else {
		ArrayKey ak;
		ak.obsID = obsid;
		ak.arrayID = arrayid;
		auto scanKeys = _msmd->getScanKeys(ak);
		for (const auto& k : scanKeys) {
			if (
				std::find(scans.begin(), scans.end(), k.scan) != scans.end()
			) {
				myKeys.insert(k);
			}
		}
	}
	ThrowIf(myKeys.empty(), "No matching scans found");
	return myKeys;
}

template <class T>
std::set<T> msmetadata::_idsFromExpansion(
	const std::map<String, std::set<T> >& mymap, const String& matchString
) {
	std::set<T> ids;
	std::regex re;
	re.assign(_escapeExpansion(matchString));
	for(auto kv : mymap) {
		if (std::regex_match(kv.first, re)) {
			ids.insert(kv.second.begin(), kv.second.end());
		}
	}
	return ids;
}

std::set<Int> msmetadata::_idsFromExpansion(
	const std::map<String, std::set<uInt> >& mymap, const String& matchString
) {
	std::set<Int> ids;
	std::regex re;
	re.assign(_escapeExpansion(matchString));
	for(auto kv : mymap) {
		if (std::regex_match(kv.first, re)) {
			ids.insert(kv.second.begin(), kv.second.end());
		}
	}
	return ids;
}

template <class T> vector<T> msmetadata::_intersection(
    const vector<T>& v, const std::set<T>& s
) {
    vector<T> r;
    for (const auto& e: v) {
        if (s.find(e) != s.end()) {
            r.push_back(e);
        }
    }
    return r;
}

template <class T> void msmetadata::_intersection(
    vector<T>& r, const std::set<T>& s1, const std::set<T>& s2
) {
    set<int> intersect;
    set_intersection(
        s1.begin(), s1.end(), s2.begin(), s2.end(),
        std::back_inserter(r)
    );
}

std::vector<casacore::String> msmetadata::_match(
	const vector<casacore::String>& candidates, const casacore::String& matchString
) {
	vector<casacore::String> matches;
	std::regex re;
	re.assign(_escapeExpansion(matchString));
	for(auto candidate : candidates) {
		if (std::regex_match(candidate, re)) {
			matches.push_back(candidate);
		}
	}
	return matches;
}

std::string msmetadata::_escapeExpansion(const casacore::String& stringToEscape) {
	const std::regex esc("[\\^\\.\\$\\|\\(\\)\\[\\]\\+\\?\\/\\\\]");
	const std::string rep("\\\\\\1");
	std::string result = regex_replace(
		stringToEscape, esc, rep, std::regex_constants::match_default | std::regex_constants::format_sed
	);
	const std::regex expand("\\*");
	const std::string rep1(".*");
	return regex_replace(
		result, expand, rep1, std::regex_constants::match_default | std::regex_constants::format_sed
	);
}

} // casac namespace
