#ifndef PLOTMSPAGEHEADERDATAMODEL_H_
#define PLOTMSPAGEHEADERDATAMODEL_H_

#include <plotms/PlotMS/PlotMS.h>
#include <QAbstractTableModel>

namespace casa {

// A Qt Data Model for the table(view) displayed in PlotMS's page header
// Implements header items layout rules specifications:
// - 1. Layout header items in a 2-column table, in their command-line order,
//      filling the table left to right first, then top to bottom ( "Z" order )
// - 2. Header items of a successive plot must be laid out starting from the
//      the beginning of the next empty row

class PlotMSPageHeaderDataModel : public QAbstractTableModel
{
	Q_OBJECT
public:
	static constexpr int physicalColumns = 3;
	static constexpr int logicalColumns = 2;
	PlotMSPageHeaderDataModel(PlotMSApp* app, QObject *parent=0);
    int rowCount(const QModelIndex &parent = QModelIndex()) const ;
    int columnCount(const QModelIndex &parent = QModelIndex()) const;
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;
    inline int lastColumnIndex() const { return physicalColumns -1 ; }
    inline bool isLogicalColumn(int column) const;
    inline int logicalColumn(int column) const;

private:
    PlotMSApp* itsParent_;
};

}





#endif /* PLOTMSPAGEHEADERDATAMODEL_H_ */
