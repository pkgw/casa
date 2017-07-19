#ifndef PAGEHEADERDATAMODEL_H_
#define PAGEHEADERDATAMODEL_H_

#include <casa/Utilities/CountedPtr.h>

//
namespace casa {

class PageHeaderDataModel {
protected:
	PageHeaderDataModel() {}
	virtual ~PageHeaderDataModel() {}
};

typedef casacore::CountedPtr<PageHeaderDataModel> PageHeaderDataModelPtr;

}


#endif /* PAGEHEADERDATAMODEL_H_ */
