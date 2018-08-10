/*
 * QPHeaderTableWidget.qo.h
 *   GUI Widget for displaying header items in a table
 *
 */

#ifndef QPHEADERTABLEWIDGET_QO_H_
#define QPHEADERTABLEWIDGET_QO_H_

// For PlotMSHeaderTable
#include <QHeaderView>
#include <QTableView>

// For HeaderTableDataModel
#include <QtCore/QAbstractTableModel>
#include <QtCore/QVariant>
#include <QtCore/QModelIndex>


namespace casa {

class HeaderTableDataModel; // Forward declaration

// Widget for displaying header items in a table
class QPHeaderTable : public QTableView {
    Q_OBJECT

public:
    // Constructor which takes optional parent widget.
    QPHeaderTable(QWidget* parent = nullptr);

    // Destructor.
    ~QPHeaderTable();

protected:
	void resizeEvent(QResizeEvent *event);

};

}


#endif /* QPHEADERTABLEWIDGET_QO_H_ */
