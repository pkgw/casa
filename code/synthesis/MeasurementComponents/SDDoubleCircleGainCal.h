/*
 * SDDoubleCircleGainCal.h
 *
 *  Created on: Jun 3, 2016
 *      Author: nakazato
 */

#ifndef SYNTHESIS_MEASUREMENTCOMPONENTS_SDDOUBLECIRCLEGAINCAL_H_
#define SYNTHESIS_MEASUREMENTCOMPONENTS_SDDOUBLECIRCLEGAINCAL_H_

#include <synthesis/MeasurementComponents/StandardVisCal.h>
#include <casacore/casa/BasicSL/String.h>

namespace casa {

class VisSet;
class VisEquation;

class SDDoubleCircleGainCal final : public GJones {
public:
  SDDoubleCircleGainCal(VisSet& vs);
  SDDoubleCircleGainCal(const MSMetaInfoForCal& msmc);
  virtual ~SDDoubleCircleGainCal();

  // Return type name as string (ditto)
//  virtual casacore::String typeName() {
//    return "SDGAIN_OTFD";
//  }
//  virtual casacore::String longTypeName() {
//    return "SDGAIN_OTFD (Single Dish gain calibration for double circle fast scan";
//  }

  // Return the parameter type
  // so far single dish calibration is real
//  virtual VisCalEnum::VCParType parType() { return VisCalEnum::REAL; }

  // Frequency-dependent Parameters?
  virtual casacore::Bool freqDepPar() { return true; };

  // useGenericGatherForSolve must return true to migrate VI/VB2 based implementation
  virtual casacore::Bool useGenericGatherForSolve() override {
    return true;
  }
  // Do not use generic data gathering mechanism for solve
  virtual casacore::Bool useGenericSolveOne() override {
    return false;
  }

  // Set the solving parameters
  virtual void setSolve() override;
  virtual void setSolve(const casacore::Record& solve) override;

  // Report solve info/params, e.g., for logging
  virtual casacore::String solveinfo() override;

  // Post solve tinkering
  virtual void globalPostSolveTinker() override;

  // Self- gather and/or solve prototypes
  //  (triggered by useGenericGatherForSolve=F or useGenericSolveOne=F)
  virtual void selfGatherAndSolve(VisSet& vs, VisEquation& ve) override;
  virtual void selfSolveOne(SDBList &sdbs) override;

  // specific keepNCT
  virtual void keepNCT() override;

private:
  template<class Accessor>
  void executeDoubleCircleGainCal(casacore::MeasurementSet const &ms);

  casacore::Double central_disk_size_;
  casacore::Bool smooth_;
  casacore::uInt currAnt_;
};

} // namespace casa END

#endif /* SYNTHESIS_MEASUREMENTCOMPONENTS_SDDOUBLECIRCLEGAINCAL_H_ */
