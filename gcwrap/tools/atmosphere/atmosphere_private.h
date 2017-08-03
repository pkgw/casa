/* Private parts for atmosphere */

/* 
   Assert ATM type is in permissive range of enum, typeAtm_t
 */
void check_atmtype_enum(int atmtype);

/* 
   Assert int value is positive or zero.
   This function is necessary because tool interface does not
   support parameter in unsigned integer.
 */
void assert_unsigned_int(int value);

/* 
   Assert SPW ID is valid number in pSpectralGrid.
 */
void assert_spwid(int spwid);
/* 
   Assert SPW and CHAN IDs are valid number in pSpectralGrid.
 */
void assert_spwid_and_channel(int spwid, int chan);


typedef unsigned int (atm::SpectralGrid::*SpGridSingleIdFuncInt) (unsigned int) const;
typedef atm::Frequency (atm::SpectralGrid::*SpGridSingleIdFuncFreq) (unsigned int) const;
// helper functions to invoke ATM functions in SpectralGrid class which take one integer id as the parameter
// returns int 
int DoSpGridSingleIdFuncInt(SpGridSingleIdFuncInt func, int spwid);
// returns quantity
casac::Quantity DoSpGridSingleIdFuncQuantum(SpGridSingleIdFuncFreq func, int spwid, std::string qunits);

// helper functions to invoke ATM functions in RefractiveIndexProfile and SkyStatus classes
// for atmosphere functions which take two integer ids as paramters
// return a double
template<typename Func, typename ClassType>
double doTwoIdATMFuncDouble(Func func, ClassType obj, int nc, int spwid);
// return a quantity
template<typename Func, typename ClassType>
casac::Quantity doTwoIdATMFuncQuantum(Func func, ClassType pobj, int nc, int spwid, std::string units);

// helper functions to invoke ATM functions in RefractiveIndexProfile class
// for atmosphere functions which take two integer ids as paramters and return a quantity
template<typename Func>
casac::Quantity doRIPThreeIdFuncQuantum(Func func, int nl, int nf, int spwid, std::string units);

// clean-up members allocated by new
inline void cleanUp() {
    if (pSkyStatus != 0) {
      delete pSkyStatus;
      pSkyStatus = 0;
    }
    if (pRefractiveIndexProfile != 0) {
      delete pRefractiveIndexProfile;
      pRefractiveIndexProfile = 0;
    }
    if (pSpectralGrid != 0) {
      delete pSpectralGrid;
      pSpectralGrid = 0;
    }
    if (pAtmProfile != 0) {
      delete pAtmProfile;
      pAtmProfile = 0;
    }
};

atm::AtmProfile *pAtmProfile;
atm::SpectralGrid *pSpectralGrid;
atm::RefractiveIndexProfile *pRefractiveIndexProfile;
atm::SkyStatus *pSkyStatus;
casacore::LogIO *itsLog;

