private:

std::unique_ptr<casacore::LogIO> _log;

std::shared_ptr<casa::ImageMetaDataRW<casacore::Float> > _mdf;
std::shared_ptr<casa::ImageMetaDataRW<casacore::Complex> > _mdc;

static const casacore::String _class;

void _exceptIfDetached() const;

