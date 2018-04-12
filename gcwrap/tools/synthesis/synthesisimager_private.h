casa::SynthesisImager *itsImager;

casacore::String checkStr(std::string instr);
casacore::LogIO *itsLog;
casa::SynthesisImager* makeSI(bool forceNew=false, bool oldvi=false);
