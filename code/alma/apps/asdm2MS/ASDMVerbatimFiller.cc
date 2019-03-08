#include <alma/apps/asdm2MS/ASDMVerbatimFiller.h>


ASDMVerbatimFiller::ASDMVerbatimFiller() {;}
//ASDMVerbatimFiller::ASDMVerbatimFiller(MS* ms_p, const std::set<const asdm::ASDM_TABLE_BASE*>& table) {
ASDMVerbatimFiller::ASDMVerbatimFiller(casacore::MS* ms_p, const std::set<asdm::ASDM_TABLE_BASE*>& table) {
  table_ = table;
  for(std::set<asdm::ASDM_TABLE_BASE*>::iterator iter = table_.begin();
      iter != table_.end(); ++iter)
    (*iter)->buildAndAttachTable(ms_p);  
}

ASDMVerbatimFiller::~ASDMVerbatimFiller() {;}

void ASDMVerbatimFiller::fill(const asdm::ASDM& asdm) {
    for (std::set<asdm::ASDM_TABLE_BASE*>::const_iterator iter = table_.begin(); iter!=table_.end(); ++iter)
    (*iter)->fill(asdm);
}
