casacore::LogIO *itsLog;
casa::ImagePol *itsImPol;

std::shared_ptr<casacore::Record> _getRegion(
    const casac::variant& region, bool nullIfEmpty
) const;

