#ifndef _ASDMVERBATIMFILLER_H_
#define _ASDMVERBATIMFILLER_H_
#include "ASDMTableBase.h"
#include "ASDMTables.h"

#include <set>
#include <alma/ASDM/ASDM.h>

class ASDMVerbatimFiller {
public:
    virtual ~ASDMVerbatimFiller();
    //  ASDMVerbatimFiller(casacore::MS* ms_p, const std::set<const asdm::ASDM_TABLE_BASE*>& table); 
    ASDMVerbatimFiller(casacore::MS* ms_p, const std::set<asdm::ASDM_TABLE_BASE*>& table); 
    void fill(const asdm::ASDM& asdm);  
private:
    std::set<asdm::ASDM_TABLE_BASE*> table_;
    ASDMVerbatimFiller();
};
#endif // _ASDMVERBATIMFILLER_H_

