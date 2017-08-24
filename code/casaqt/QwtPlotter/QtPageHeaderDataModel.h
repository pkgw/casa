#ifndef QTPAGEHEADERDATAMODEL_H_
#define QTPAGEHEADERDATAMODEL_H_

#include <QtCore/QAbstractItemModel>

#include <graphics/GenericPlotter/PageHeaderDataModel.h>

namespace casa {

// Qt-based implementation of PageHeaderDataModel
class QtPageHeaderDataModel : public PageHeaderDataModel {
public:
	QtPageHeaderDataModel(QAbstractItemModel *qtModel=nullptr) : qtModel_(qtModel) {}
	QAbstractItemModel *model() { return qtModel_; }
private:
	QAbstractItemModel *qtModel_;
};
typedef casacore::CountedPtr<QtPageHeaderDataModel> QtPageHeaderDataModelPtr;


}



#endif /* QTPAGEHEADERDATAMODEL_H_ */
