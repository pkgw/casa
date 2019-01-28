//# PlotMSIterationTab.qo.h: Plot tab for iterated plot settings 
//# Copyright (C) 2009
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id: $
#ifndef PLOTMSITERATETAB_QO_H_
#define PLOTMSITERATETAB_QO_H_

#include <plotms/GuiTabs/PlotMSIterateTab.ui.h>

#include <plotms/GuiTabs/PlotMSPlotTab.qo.h>
#include <plotms/Plots/PlotMSPlotParameterGroups.h>

#include <array>
#include <QAbstractTableModel>
#include <QIcon>

namespace casa {

//# Forward declarations.
class QtLabelWidget;


// Subclass of PlotMSPlotSubtab to manage plot display parameters.
class PlotMSIterateTab : public PlotMSPlotSubtab, Ui::IterateTab {
    Q_OBJECT

    
public:

    // Constructor which takes the parent tab and plotter.
    PlotMSIterateTab(PlotMSPlotTab* plotTab, PlotMSPlotter* parent);
    
    // Destructor.
    ~PlotMSIterateTab();
    
    
    // Implements PlotMSTab::tabName().
    QString tabName() const { return "Page"; }
    
    // Implements PlotMSPlotSubtab::getValue().
    void getValue(PlotMSPlotParameters& params) const;
    
    // Implements PlotMSPlotSubtab::setValue().
    void setValue(const PlotMSPlotParameters& params);
    
    // Implements PlotMSPlotSubtab::update().
    void update(const PlotMSPlot& plot);

    // Uses the index chooser at the top, with the given number of rows and
    // columns, to manage multi-plot display parameters.
    bool setGridSize(unsigned int nRows,unsigned int nCols);

    //Returns true if a reasonable row and column location has been
    //set (nonzero); false otherwise.
    bool isPlottable() const;

    // Clear user-selected header items
	void clearSelectedItems();

signals:
	void plottableChanged();

private slots:
	//Whether to use a single global axis has changed.
	void globalChanged();
	void locationChanged();
	// Header items edition
	// ---- Add items selected in availableItemsTableView
	void addItems();
	// ---- Remove items selected in selectedItemsTableView
	void removeItems();
	void headerItemsChanged();
	// Icons animation
	void setArrowUpNormal();
	void setArrowUpActive();
	void setArrowDownNormal();
	void setArrowDownActive();

private:
	void hideGridLocation( bool hide );
	void setGridIndices( int rowIndex, int colIndex );
	void setHeaderItems(const PMS_PP_PageHeader* pageHeaderParamsGroup);
	void setAvailable(PageHeaderItemsDef::Item item, bool isAvailable);
	void setButtonIcon(QPushButton *button, Bool isActive);

	//Location of the plot
	int gridRow;
	int gridCol;

	// Page header items
	PageHeaderItems headerItems;

	// Animated buttons
	QIcon arrowUpIcon;
	QIcon arrowDownIcon;
};


class SelectedItemsDataModel : public QAbstractTableModel
{
	Q_OBJECT

public:
	using Items = PageHeaderItemsDef;
	using Item = Items::Item;
	SelectedItemsDataModel(PageHeaderItems &items, QObject *parent =0);
	// Read-Only interface
	int rowCount(const QModelIndex &parent = QModelIndex()) const;
	int columnCount(const QModelIndex &parent = QModelIndex()) const;
	QVariant data(const QModelIndex &index, int role) const;
	// Resize interface
    void prepareToAppend(Item item);
    std::vector<Item>& itemsToRemove() { return itemsToRemove_ ; }
    bool insertRows(int row, int count, const QModelIndex & parent = QModelIndex());
    bool removeRows(int row, int count, const QModelIndex & parent = QModelIndex());

    // Edit interface
    // bool setData(const QModelIndex & index, const QVariant & value, int role = Qt::EditRole);
private:
	PageHeaderItems &items_;
	std::vector<Item> newItem_;
	std::vector<Item> itemsToRemove_;
};


class AvailableItemsDataModel : public QAbstractTableModel
{
	Q_OBJECT

public:
	using Items = PageHeaderItemsDef;
	using Item = Items::Item;
	using ItemsArray = std::array<Item,Items::n_items>;
	AvailableItemsDataModel(const ItemsArray &items, QObject *parent =0);
	// Read-Only interface
	int rowCount(const QModelIndex &parent = QModelIndex()) const;
	int columnCount(const QModelIndex &parent = QModelIndex()) const;
	QVariant data(const QModelIndex &index, int role) const;
private:
	const ItemsArray &items_;
};
    
    

    




}

#endif /* PLOTMSDISPLAYTAB_QO_H_ */
