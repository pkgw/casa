#ifndef _ATM_WATERVAPORRADIOMETER_H
#define _ATM_WATERVAPORRADIOMETER_H
/*******************************************************************************
 * ALMA - Atacama Large Millimiter Array
 * (c) Instituto de Estructura de la Materia, 2009
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 *
 * "@(#) $Id: ATMWaterVaporRadiometer.h Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#ifndef __cplusplus
#error This is a C++ include file and cannot be used from plain C
#endif

#include "ATMCommon.h"
#include "ATMPercent.h"
#include "ATMTemperature.h"

#include <vector>

using std::vector;

ATM_NAMESPACE_BEGIN

/*! \brief This is an axuliary class to SkyStatus that allows to identify the spectral windows
 *         that actually correspond to a dedicated water vapor radiometry system, and also to
 *         associate to those spectral windows some instrumental caracteristics such as the signal
 *         side band gain, the sky coupling or the spillover temperature.
 */
class WaterVaporRadiometer
{
public:

  WaterVaporRadiometer() {}

  /** Class constructor with identifiers of radiometric channels. Sky Coupling = 1.0, Sideband Gain: 50% (no spilloverTemperature given) */
  WaterVaporRadiometer(const std::vector<unsigned int> &IdChannels);

  /** Class constructor with identifiers of radiometric channels. Sideband Gain: 50% (no spilloverTemperature given) */
  WaterVaporRadiometer(const std::vector<unsigned int> &IdChannels,
                       const std::vector<double> &skyCoupling);

  inline WaterVaporRadiometer(const std::vector<unsigned int> &IdChannels,
                              double skyCoupling)
  {
    WaterVaporRadiometer(IdChannels, std::vector<double> (IdChannels.size(),
                                                     skyCoupling));
  }

  /** Class constructor with identifiers of radiometric channels. Sky Coupling = 1.0 (no spilloverTemperature given) */
  WaterVaporRadiometer(const std::vector<unsigned int> &IdChannels,
                       const std::vector<Percent> &signalGain);

  /* WaterVaporRadiometer(vector<unsigned int> IdChannels, Percent signalGain); */

  /** Full class constructor (no spilloverTemperature given) */
  WaterVaporRadiometer(const std::vector<unsigned int> &IdChannels,
                       const std::vector<double> &skyCoupling,
                       const std::vector<Percent> &signalGain);

  /** Class constructor with identifiers of radiometric channels. Sky Coupling = 1.0, Sideband Gain: 50% (spilloverTemperature given) */
  WaterVaporRadiometer(const std::vector<unsigned int> &IdChannels,
                       const Temperature &spilloverTemperature);

  /** Class constructor with identifiers of radiometric channels. Sideband Gain: 50% (spilloverTemperature given) */
  WaterVaporRadiometer(const std::vector<unsigned int> &IdChannels,
                       const std::vector<double> &skyCoupling,
                       const Temperature &spilloverTemperature);

  /** Class constructor with identifiers of radiometric channels. Sky Coupling = 1.0 (spilloverTemperature given) */
  WaterVaporRadiometer(const std::vector<unsigned int> &IdChannels,
                       const std::vector<Percent> &signalGain,
                       const Temperature &spilloverTemperature);

  /** Full class constructor (spilloverTemperature given) */
  WaterVaporRadiometer(const std::vector<unsigned int> &IdChannels,
                       const std::vector<double> &skyCoupling,
                       const std::vector<Percent> &signalGain,
                       const Temperature &spilloverTemperature);

  virtual ~WaterVaporRadiometer();

  /** Accessor to identifiers of radiometric channels (vector of int) */
  inline std::vector<unsigned int> getIdChannels() const { return IdChannels_; }

  /** Accessor to Sky coupling of the different radiometric channels (vector of double) */
  inline std::vector<double> getSkyCoupling() const { return skyCoupling_; }

  /** Multiplier of the Sky coupling vector  */
  inline void multiplySkyCoupling(double factor)
  {
    for(unsigned int i = 0; i < skyCoupling_.size(); i++) {
      skyCoupling_[i] = skyCoupling_[i] * factor;
    }
  }
  /** Multiplier of one  Sky coupling channel  */
  inline void multiplySkyCouplingChannel(unsigned int ichan, double factor)
  {
    skyCoupling_[ichan] = skyCoupling_[ichan] * factor;
  }

  /** Setter of the Sky coupling vector to a single Sky Coupling value */
  inline void setSkyCoupling(double factor)
  {
    /* for(unsigned int i=0; i<skyCoupling_.size(); i++){ skyCoupling_[i] = factor;} */
    for(std::vector<double>::iterator it(skyCoupling_.begin()),
        itmax(skyCoupling_.end()); it != itmax; ++it) {
      (*it) = factor;
    }
  }
  inline void setSkyCoupling(unsigned int i, double factor)
  {
    if(i < skyCoupling_.size()) {
      skyCoupling_[i] = factor;
    } else { /* throw(Error()); */
    }
  }

  /** Accessor to signalGain of the different radiometric channels (vector of percent) */
  inline std::vector<Percent> getsignalGain() const { return signalGain_; }
  /** Accessor to SpilloverTemperature */
  inline Temperature getSpilloverTemperature() const { return spilloverTemperature_; }
  /** Setter of SpilloverTemperature */
  inline void setSpilloverTemperature(Temperature spilloverTemperature) { spilloverTemperature_ = spilloverTemperature; }

protected:
  std::vector<unsigned int> IdChannels_;
  std::vector<double> skyCoupling_;
  std::vector<Percent> signalGain_;
  Temperature spilloverTemperature_;

private:
  //     /** Inaccessible Copy Constructor */
  //     WaterVaporRadiometer(const WaterVaporRadiometer&);

  //      /** Inaccessible assignment operator */
  //      WaterVaporRadiometer& operator=(const WaterVaporRadiometer&);


}; // class WaterVaporRadiometer

ATM_NAMESPACE_END

#endif /*!_ATM_WATERVAPORRADIOMETER_H*/
