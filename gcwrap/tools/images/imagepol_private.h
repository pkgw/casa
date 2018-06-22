casacore::LogIO *itsLog;
casa::ImagePol *itsImPol;

SHARED_PTR<casacore::Record> _getRegion(
    const casac::variant& region, bool nullIfEmpty
) const;

