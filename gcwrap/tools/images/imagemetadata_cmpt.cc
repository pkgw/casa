#include <imagemetadata_cmpt.h>

#include <casa/Containers/ValueHolder.h>
#include <casa/Utilities/PtrHolder.h>
#include <imageanalysis/ImageAnalysis/ImageFactory.h>
#include <imageanalysis/ImageAnalysis/ImageMetaDataRW.h>
#include <images/Images/ImageOpener.h>

#include <stdcasa/version.h>

#include <casa/namespace.h>

#include <memory>

using namespace std;

using namespace casacore;
using namespace casa;

namespace casac {

const String imagemetadata::_class = "imagemetadata";

imagemetadata::imagemetadata() :
_log(new LogIO()), _mdf(), _mdc() {}

imagemetadata::~imagemetadata() {}

bool imagemetadata::close() {
	_mdf.reset();
	_mdc.reset();
	return true;
}

bool imagemetadata::add(const string& key, const variant& value) {
	try {
		_exceptIfDetached();
		std::unique_ptr<const ValueHolder> vh(casa::toValueHolder(value));
		if (_mdf) {
		    return _mdf->add(key, *vh);
		}
		else if(_mdc) {
		    return _mdc->add(key, *vh);
		}
		else {
		    ThrowCc("Logic error");
		}
	}
	catch (const AipsError& x) {
		*_log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
			<< LogIO::POST;
		RETHROW(x);
	}
}

bool imagemetadata::done() {
	return close();
}

variant* imagemetadata::get(const string& key) {
	try {
		_exceptIfDetached();
		return casa::fromValueHolder(
			_mdf ? _mdf->getFITSValue(key) : _mdc->getFITSValue(key)
		);
	}
	catch (const AipsError& x) {
		*_log << LogIO::SEVERE << "Exception Reported: "
			<< x.getMesg() << LogIO::POST;
		RETHROW(x);
	}
}

record* imagemetadata::list(bool verbose) {
	try {
		_exceptIfDetached();
		return casa::fromRecord(_mdf ? _mdf->toRecord(verbose) : _mdc->toRecord(verbose));
	}
	catch (const AipsError& x) {
		*_log << LogIO::SEVERE << "Exception Reported: "
			<< x.getMesg() << LogIO::POST;
		RETHROW(x);
	}
}

bool imagemetadata::open(const std::string& infile) {
	try {
		if (! _log.get()) {
			_log.reset(new LogIO());
		}
        SPIIF imageF;
        SPIIC imageC;
        std::tie(imageF, imageC, std::ignore, std::ignore)
            = ImageFactory::fromFile(infile);
        if (imageF) {
			_mdf.reset(new ImageMetaDataRW<Float>(imageF));
		}
		else if (imageC) {
			_mdc.reset(new ImageMetaDataRW<Complex>(imageC));
		}
        else {
            ThrowCc("Image type not yet supported");
        }
		return true;
	}
	catch (const AipsError& x) {
		*_log << LogIO::SEVERE << "Exception Reported: "
			<< x.getMesg() << LogIO::POST;
		RETHROW(x);
	}
}

bool imagemetadata::remove(
	const string& key, const variant& value
) {
	try {
		_exceptIfDetached();
		if (String(key) == ImageMetaDataConstants::MASKS) {
			return _mdf ? _mdf->removeMask(value.toString())
			    : _mdc->removeMask(value.toString());
		}
		else {
			return _mdf ? _mdf->remove(key) : _mdc->remove(key);
		}
	}
	catch (const AipsError& x) {
		*_log << LogIO::SEVERE << "Exception Reported: "
			<< x.getMesg() << LogIO::POST;
		RETHROW(x);
	}
}

bool imagemetadata::set(const string& key, const variant& value) {
	try {
		_exceptIfDetached();
		std::unique_ptr<const ValueHolder> vh(toValueHolder(value));
		return _mdf ? _mdf->set(key, *vh) : _mdc->set(key, *vh);
	}
	catch (const AipsError& x) {
		*_log << LogIO::SEVERE << "Exception Reported: "
			<< x.getMesg() << LogIO::POST;
		RETHROW(x);
	}
}

void imagemetadata::_exceptIfDetached() const {
	ThrowIf(
	    ! (_mdf || _mdc),
	    "Tool is not attached to a metadata object. Call open() first."
	);
}

} // casac namespace
