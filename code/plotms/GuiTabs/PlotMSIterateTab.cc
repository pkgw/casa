//# PlotMSIterateTab.cc: Plot tab to manage plot display parameters.
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
#include <plotms/GuiTabs/PlotMSIterateTab.qo.h>


#include <casaqt/QtUtilities/QtPlotWidget.qo.h>
#include <casaqt/QtUtilities/QtUtilities.h>
#include <plotms/Gui/PlotMSPlotter.qo.h>
#include <plotms/Plots/PlotMSPlot.h>
#include <plotms/Plots/PlotMSPlotParameterGroups.h>
#include <plotms/PlotMS/PlotMSPageHeaderParam.h>

#include <QDebug>
using namespace casacore;
namespace casa {

//////////////////////////////////
// DATA MODELS DEFINITIONS      //
//////////////////////////////////

AvailableItemsDataModel::AvailableItemsDataModel(const ItemsArray &items, QObject *parent)
	: QAbstractTableModel(parent),
	  items_(items)
{

}

// Read-Only interface
int AvailableItemsDataModel::rowCount(const QModelIndex & /* parent */) const {
	return static_cast<int>(Items::n_items);
}
int AvailableItemsDataModel::columnCount(const QModelIndex & /* parent */) const {
	return 2;
}
QVariant AvailableItemsDataModel::data(const QModelIndex &index, int role) const {
	if ( ! index.isValid() ) return QVariant();
	if ( index.row() >= static_cast<int>(items_.size()) ) return QVariant();
	if (role == Qt::DisplayRole) {
		auto item = items_[index.row()];
		switch(index.column()){
		case 0:
			return static_cast<int>(item);
			break;
		case 1:
			return QString::fromUtf8(Items::item2Info.at(item).longLabel().c_str());
			break;
		default:
			return QVariant();
		}
	}

	return QVariant();
}

SelectedItemsDataModel::SelectedItemsDataModel(PageHeaderItems &items, QObject *parent)
	: QAbstractTableModel(parent),
	  items_(items)
{

}

// Read-Only interface
int SelectedItemsDataModel::rowCount(const QModelIndex & /* parent */) const {
	return static_cast<int>(items_.items().size());
}

int SelectedItemsDataModel::columnCount(const QModelIndex & /* parent */) const {
	return 2;
}

QVariant SelectedItemsDataModel::data(const QModelIndex &index, int role) const {
	if ( ! index.isValid() ) return QVariant();
	if ( index.row() >= static_cast<int>(items_.items().size()) ) return QVariant();
	if (role == Qt::DisplayRole) {
		auto item = items_.items()[index.row()];
		switch(index.column()){
		case 0:
			return static_cast<int>(item);
			break;
		case 1:
			return QString::fromUtf8(Items::item2Info.at(item).longLabel().c_str());
			break;
		default:
			return QVariant();
		}
	}

	return QVariant();
}

// Resize interface
void SelectedItemsDataModel::prepareToAppend(Item item){
	newItem_.clear();
	newItem_.push_back(item);
}

// Resize interface
bool SelectedItemsDataModel::insertRows(int row, int count, const QModelIndex & parent)
{
	// http://doc.qt.io/qt-4.8/qabstractitemmodel.html#InsertRows
	// Inserts count rows into the model before the given row
	//   . If row is rowCount(),the rows are appended to any existing rows in the parent

	// Implement support only for appending 1 item at a time
	if (  (count != 1)        ||
	      (row != rowCount()) ||
	      (items_.full())     ||
		  (newItem_.empty())  ) { // Value of new item must be known
		newItem_.clear();
		return false;
	}
	auto newItem = newItem_[0];
	if ( items_.selected(newItem) ) {
		newItem_.clear();
		return false;
	}

	int firstNewRowIndex = rowCount();
	int lastNewRowIndex = firstNewRowIndex;
	// http://doc.qt.io/qt-4.8/qabstractitemmodel.html#beginInsertRows
	// When reimplementing insertRows() in a subclass, you MUST call this function
	// BEFORE inserting data into the model's underlying data store.
	beginInsertRows(parent,firstNewRowIndex,lastNewRowIndex);
	items_.append(newItem);
	// http://doc.qt.io/qt-4.8/qabstractitemmodel.html#endInsertRows
	// When reimplementing insertRows() in a subclass, you MUST call this function
	// AFTER inserting data into the model's underlying data store.
	endInsertRows();
	newItem_.clear();

	return true;
}

bool SelectedItemsDataModel::removeRows(int row, int count, const QModelIndex & parent )
{
	// Implement support for removing 1 item at a time
	if (  (count != 1) || ( row >= rowCount() )  ) return false;

	auto headerItem = items_.items()[row];

	auto firstRowToRemove = row;
	auto lastRowToRemove  = firstRowToRemove;
	beginRemoveRows(parent,firstRowToRemove,lastRowToRemove);
	auto itemWasRemoved = items_.remove(headerItem);
	endRemoveRows();

	return itemWasRemoved;
}


//////////////////////////////////
// PLOTMSDISPLAYTAB DEFINITIONS //
//////////////////////////////////

PlotMSIterateTab::PlotMSIterateTab(PlotMSPlotTab* tab, PlotMSPlotter* parent) 
	: PlotMSPlotSubtab(tab, parent) 
{
    setupUi(this);
    gridRow = 1;
    gridCol = 1;
    
    // Fill list of available iteration axis choices
    // For now, is same as colorize list
    iterationAxisChooser->addItem(PMS::axis(PMS::NONE).c_str());
    iterationAxisChooser->addItem(PMS::axis(PMS::SCAN).c_str());
    iterationAxisChooser->addItem(PMS::axis(PMS::FIELD).c_str());
    iterationAxisChooser->addItem(PMS::axis(PMS::SPW).c_str());
    iterationAxisChooser->addItem(PMS::axis(PMS::BASELINE).c_str());
    iterationAxisChooser->addItem(PMS::axis(PMS::ANTENNA).c_str());
    iterationAxisChooser->addItem(PMS::axis(PMS::TIME).c_str());
    iterationAxisChooser->addItem(PMS::axis(PMS::CORR).c_str());
    /* not yet:
    iterationAxisChooser->addItem(PMS::axis(PMS::ANTENNA1).c_str());
    iterationAxisChooser->addItem(PMS::axis(PMS::ANTENNA2).c_str());
    iterationAxisChooser->addItem(PMS::axis(PMS::CHANNEL).c_str());
    */
    hideGridLocation( true );

    // Set up label defaults.
	itsLabelDefaults_.insert(iterationAxisChooserLabel,  iterationAxisChooserLabel->text());
	itsLabelDefaults_.insert(rowIndexLabel, rowIndexLabel->text());
	itsLabelDefaults_.insert(colIndexLabel, colIndexLabel->text());

	// Page header items
	using Items = PageHeaderItemsDef;
	// ---- Available items
	auto availableItemsDataModel = new AvailableItemsDataModel(Items::items());
	availableItemsDataModel->setParent(availableItemsTableView);
	availableItemsTableView->setModel(availableItemsDataModel);
	availableItemsTableView->hideColumn(0);
	availableItemsTableView->verticalHeader()->hide();
	auto availableItemsH_Header = availableItemsTableView->horizontalHeader();
	availableItemsH_Header->setStretchLastSection(true);
	availableItemsH_Header->hide();
	// ---- Selected items
	auto selectedItemsDataModel = new SelectedItemsDataModel(headerItems);
	selectedItemsDataModel->setParent(selectedItemsTableView);
	selectedItemsTableView->setModel(selectedItemsDataModel);
	selectedItemsTableView->hideColumn(0);
	selectedItemsTableView->verticalHeader()->hide();
	auto selectedItemsH_Header = selectedItemsTableView->horizontalHeader();
	selectedItemsH_Header->setStretchLastSection(true);
	selectedItemsH_Header->hide();

	// ---- Buttons animation
	arrowUpIcon = addItemsButton->icon();
	// -------- Create arrow down icons: vertical flip
	auto arrowUpNormalPixmap = arrowUpIcon.pixmap(QSize(128,128),QIcon::Mode::Normal,QIcon::State::Off);
	QImage arrowDownNormalImg(arrowUpNormalPixmap.toImage().mirrored(false,true));
	arrowDownIcon.addPixmap(QPixmap().fromImage(arrowDownNormalImg),QIcon::Mode::Normal,QIcon::State::Off);

	auto arrowUpActivePixmap = arrowUpIcon.pixmap(QSize(128,128),QIcon::Mode::Active,QIcon::State::On);
	QImage arrowDownActiveImg(arrowUpActivePixmap.toImage().mirrored(false,true));
	arrowDownIcon.addPixmap(QPixmap().fromImage(arrowDownActiveImg),QIcon::Mode::Active,QIcon::State::On);

	auto arrowUpActiveErrPixmap = arrowUpIcon.pixmap(QSize(128,128),QIcon::Mode::Active,QIcon::State::Off);
	QImage arrowDownActiveErrImg(arrowUpActiveErrPixmap.toImage().mirrored(false,true));
	arrowDownIcon.addPixmap(QPixmap().fromImage(arrowDownActiveErrImg),QIcon::Mode::Active,QIcon::State::Off);

	removeItemsButton->setIcon(arrowDownIcon);

    // Connect widgets.
    connect(iterationAxisChooser, SIGNAL(currentIndexChanged(int)), SIGNAL(changed()) );
    connect(globalXCheck, SIGNAL(stateChanged(int)), SLOT(globalChanged()));
    connect(globalYCheck, SIGNAL(stateChanged(int)), SLOT(globalChanged()));
    connect(sharedXCheck, SIGNAL(stateChanged(int)), SIGNAL(changed()));
    connect(sharedYCheck, SIGNAL(stateChanged(int)), SIGNAL(changed()));
    connect(gridRowSpin, SIGNAL(valueChanged(int)), SLOT(locationChanged()));
    connect(gridColSpin, SIGNAL(valueChanged(int)), SLOT(locationChanged()));
    connect(addItemsButton, SIGNAL(clicked()), SLOT(addItems()));
    connect(removeItemsButton, SIGNAL(clicked()), SLOT(removeItems()));
    connect(addItemsButton,SIGNAL(pressed()),this,SLOT(setArrowUpActive()));
    connect(addItemsButton,SIGNAL(released()),this,SLOT(setArrowUpNormal()));
    connect(removeItemsButton,SIGNAL(pressed()),this,SLOT(setArrowDownActive()));
    connect(removeItemsButton,SIGNAL(released()),this,SLOT(setArrowDownNormal()));

    //Initialize the check box indicating whether to use a shared scale and/or axis.
    sharedXCheck->setEnabled( false );
    sharedYCheck->setEnabled( false );
    globalXCheck->setEnabled( true );
    globalYCheck->setEnabled( true );
}



bool PlotMSIterateTab::isPlottable() const {
	bool plottable = true;
	if ( gridRowSpin->value() == 0 || gridColSpin->value() == 0 ){
		plottable = false;
	}
	return plottable;
}

void PlotMSIterateTab::clearSelectedItems(){
	selectedItemsTableView->selectAll();
	removeItems();
}

void PlotMSIterateTab::locationChanged(){
	emit changed();
}

void PlotMSIterateTab::headerItemsChanged(){
	emit changed();
}

void PlotMSIterateTab::addItems() {
	auto labelsIndexes = availableItemsTableView->selectionModel()->selectedIndexes();
	auto selectedItemsModel = dynamic_cast<SelectedItemsDataModel *>(selectedItemsTableView->model());
	if ( selectedItemsModel == nullptr ) return;
	auto availableItemsModel = availableItemsTableView->model();
	bool itemAdded = false;
	for ( const auto & labelIndex : labelsIndexes ) {
		auto itemIndex = availableItemsModel->index(labelIndex.row(),0);
		QVariant itemValue = itemIndex.data(Qt::DisplayRole);
		bool haveIntValue = false;
		auto itemIntValue = itemValue.toInt(&haveIntValue);
		if ( ! haveIntValue ) continue;
		using Item = PageHeaderItemsDef::Item;
		auto newItem = static_cast<Item>(itemIntValue);
		selectedItemsModel->prepareToAppend(newItem);
		auto fromRow = selectedItemsModel->rowCount();
		auto itemInserted = selectedItemsModel->insertRows(fromRow,1);
		if ( itemInserted ) {
			availableItemsTableView->setRowHidden(labelIndex.row(),true);
			itemAdded = true;
		}
	}
	availableItemsTableView->selectionModel()->clear();
	if (itemAdded) headerItemsChanged();
}

void PlotMSIterateTab::setAvailable(PageHeaderItemsDef::Item item, bool isAvailable){

	auto availableItemsDataModel = availableItemsTableView->model();
	const int headerItemsColumn = 0;
	for ( int row = 0; row < availableItemsDataModel->rowCount(); ++row ) {
		auto currentModelIndex = availableItemsDataModel->index(row,headerItemsColumn);
		QVariant currentValue = currentModelIndex.data(Qt::DisplayRole);
		bool haveIntValue = false;
		auto currentIntValue = currentValue.toInt(&haveIntValue);
		if ( ! haveIntValue ) continue;
		auto currentHeaderItem = static_cast<PageHeaderItemsDef::Item>(currentIntValue);
		if ( currentHeaderItem == item ){
			availableItemsTableView->setRowHidden(row, ! isAvailable);
			break;
		}
	}
}

void PlotMSIterateTab::removeItems() {
	auto labelsIndexes = selectedItemsTableView->selectionModel()->selectedIndexes();
	auto selectedItemsModel = dynamic_cast<SelectedItemsDataModel *>(selectedItemsTableView->model());
	if ( selectedItemsModel == nullptr ) return;
	// Deleting items will invalidate selection ...
	auto headerItemsToRemove = selectedItemsModel->itemsToRemove();
	headerItemsToRemove.clear();
	for ( const auto & labelIndex : labelsIndexes ) {
		const int headerItemsColumn = 0;
		auto headerItemModelIndex = selectedItemsModel->index(labelIndex.row(),headerItemsColumn);
		QVariant headerItemIndexValue = headerItemModelIndex.data(Qt::DisplayRole);
		bool haveIntValue = false;
		auto itemIntValue = headerItemIndexValue.toInt(&haveIntValue);
		if ( ! haveIntValue ) continue;
		using HeaderItem = PageHeaderItemsDef::Item;
		auto headerItemToRemove = static_cast<HeaderItem>(itemIntValue);
		headerItemsToRemove.push_back(headerItemToRemove);
	}
	bool itemRemoved = false;
	for ( auto headerItem : headerItemsToRemove ){
		auto foundIndexPair = headerItems.index(headerItem);
		if ( ! foundIndexPair.first ) continue;
		int rowToRemove = static_cast<int>(foundIndexPair.second);
		auto rowRemoved = selectedItemsModel->removeRow(rowToRemove);
		if ( rowRemoved ) {
			setAvailable(headerItem,true);
			itemRemoved = true;
		}
	}
	selectedItemsTableView->selectionModel()->clear();
	if (itemRemoved) headerItemsChanged();
}



bool PlotMSIterateTab::setGridSize(unsigned int nRows,unsigned int nCols){
    //Decide whether we let the user choose the location of the plot
	//or not based on the grid size.
	bool hideGrid = true;
    if ( nRows > 1 || nCols > 1 ){
    	hideGrid = false;
    }
    hideGridLocation( hideGrid );

    //Decide whether the current location exceeds the number of rows/columns
    //available now.
    bool validLocation = true;
    int currentRow = gridRowSpin->value() - 1;
    int currentCol = gridColSpin->value() - 1;
    if ( currentRow >= static_cast<int>(nRows) ||
    		currentCol >= static_cast<int>(nCols) ){
    	validLocation = false;
    }

    //Reset the limits on the spins.
    gridRowSpin->setMaximum( nRows );
    gridColSpin->setMaximum( nCols );
    //emit plottableChanged();
    return validLocation;
}

void PlotMSIterateTab::globalChanged(){
	bool globalYChecked = globalYCheck->isChecked();
	bool globalXChecked = globalXCheck->isChecked();
	sharedYCheck->setEnabled( globalYChecked );
	sharedXCheck->setEnabled( globalXChecked );
	if ( !globalYChecked ){
		sharedYCheck->setChecked( false );
	}
	if ( !globalXChecked ){
		sharedXCheck->setChecked( false );
	}
	emit changed();
}

void PlotMSIterateTab::hideGridLocation( bool hide ) {
    chooserFrame->setVisible(!hide);
}


void PlotMSIterateTab::setGridIndices( int rowIndex, int colIndex ){
	//Note that the rows and columns we display to the user begin with 1,
	//whereas what we store in the code begins at 0.
	int maxRows = gridRowSpin->maximum();
	if ( 0<= rowIndex && rowIndex <= maxRows ){
		gridRowSpin->setValue( rowIndex );
	}
	else {
	    qDebug() << "PlotMSIterateTab::setValue maxRows="<<maxRows<<" rowIndex="<<rowIndex;
	}
	int maxCols = gridColSpin->maximum();
	if ( 0 <= colIndex && colIndex <= maxCols ){
		gridColSpin->setValue( colIndex );
	}
	else {
		 qDebug() << "PlotMSIterateTab::setValue maxCols="<<maxCols<<" colIndex="<<colIndex;
	}
	gridRow = rowIndex;
	gridCol = colIndex;
}

PlotMSIterateTab::~PlotMSIterateTab() { }



void PlotMSIterateTab::getValue(PlotMSPlotParameters& params) const   {
    PMS_PP_Iteration* d = params.typedGroup<PMS_PP_Iteration>();
    if(d == NULL) {
        params.setGroup<PMS_PP_Iteration>();
        d = params.typedGroup<PMS_PP_Iteration>();
    }
        
    d->setIterationAxis(
           PMS::axis(iterationAxisChooser->currentText().toStdString()) );

	d->setGlobalScaleX( globalXCheck->isChecked());
	d->setGlobalScaleY( globalYCheck->isChecked());

	// only set if grid
	d->setCommonAxisX( gridRowSpin->maximum()>1 && sharedXCheck->isChecked());
	d->setCommonAxisY( gridColSpin->maximum()>1 && sharedYCheck->isChecked());

	//Note, we are subtracting 1 because UI values start with 1,
	//but internal values are 0 based.
	int rowSpinIndex = gridRowSpin->value() - 1;
	int colSpinIndex = gridColSpin->value() - 1;
	d->setGridRow( rowSpinIndex );
	d->setGridCol( colSpinIndex );

	PMS_PP_PageHeader* pageHeaderParamsGroup = params.typedGroup<PMS_PP_PageHeader>();
	if (pageHeaderParamsGroup == NULL){
        params.setGroup<PMS_PP_PageHeader>();
        pageHeaderParamsGroup = params.typedGroup<PMS_PP_PageHeader>();
	}
	pageHeaderParamsGroup->setPageHeaderItems(headerItems);
}



void PlotMSIterateTab::setValue(const PlotMSPlotParameters& params) {

  const PMS_PP_Iteration* d = params.typedGroup<PMS_PP_Iteration>();
  
  PlotMSTab::setChooser(iterationAxisChooser, PMS::axis(d->iterationAxis()));

  sharedXCheck->setChecked( d->isCommonAxisX() );
  sharedYCheck->setChecked( d->isCommonAxisY() );
  globalXCheck->setChecked( d->isGlobalScaleX() );
  globalYCheck->setChecked( d->isGlobalScaleY() );
  globalChanged();

  int rowIndex = d->getGridRow();
  int colIndex = d->getGridCol();
  setGridIndices( rowIndex+ 1, colIndex+1 );

  setHeaderItems(params.typedGroup<PMS_PP_PageHeader>());

}

void PlotMSIterateTab::setHeaderItems(const PMS_PP_PageHeader* pageHeaderParamsGroup){
	auto selectedItemsDataModel = dynamic_cast<SelectedItemsDataModel *>(selectedItemsTableView->model());
	if (selectedItemsDataModel == nullptr) return;

	// Clear existing items
	clearSelectedItems();

	// Set new ones if any
	if (pageHeaderParamsGroup == nullptr) return;
	for ( const auto & item : pageHeaderParamsGroup->pageHeaderItems().items() ) {
		selectedItemsDataModel->prepareToAppend(item);
		auto fromRow = selectedItemsDataModel->rowCount();
		selectedItemsDataModel->insertRows(fromRow,1);
		setAvailable(item,false);
	}
}

void PlotMSIterateTab::update(const PlotMSPlot& plot) {    

	const PlotMSPlotParameters &params = plot.parameters();
	PlotMSPlotParameters  newParams(params);
	getValue(newParams);
	
	const PMS_PP_Iteration *d = params.typedGroup<PMS_PP_Iteration>();
	const PMS_PP_Iteration *d2 = newParams.typedGroup<PMS_PP_Iteration>();
	
	if (d==0 || d2==0)  return;

	highlightWidgetText(iterationAxisChooserLabel, d->iterationAxis() != d2->iterationAxis() );

	bool commonXChange = d->isCommonAxisX() != d2->isCommonAxisX();
	bool commonYChange = d->isCommonAxisY() != d2->isCommonAxisY();
	bool globalXChange = d->isGlobalScaleX() != d2->isGlobalScaleX();
	bool globalYChange = d->isGlobalScaleY() != d2->isGlobalScaleY();
	highlightWidgetText(commonAxisLabel, commonXChange || commonYChange );
	highlightWidgetText(globalAxisLabel, globalXChange || globalYChange );

	bool gridRowChanged = false;
	int oldRowCount = d->getGridRow();
	int newRowCount = d2->getGridRow();
	if (oldRowCount != newRowCount ){
		 gridRowChanged = true;
	}
	highlightWidgetText( rowIndexLabel, gridRowChanged );
	bool gridColChanged = false;
	int oldColCount = d->getGridCol();
	int newColCount = d2->getGridCol();
	if (oldColCount != newColCount ){
		 gridColChanged = true;
	}
	highlightWidgetText(colIndexLabel, gridColChanged );

	auto oldPageHeaderParamsGroup = params.typedGroup<PMS_PP_PageHeader>();
	auto newPageHeaderParamsGroup = newParams.typedGroup<PMS_PP_PageHeader>();

	if (oldPageHeaderParamsGroup==NULL || newPageHeaderParamsGroup==NULL) return;

	auto oldItems = oldPageHeaderParamsGroup->pageHeaderItems();
	auto newItems = newPageHeaderParamsGroup->pageHeaderItems();
	bool headerItemsChanged = ( oldItems != newItems );
	highlightWidgetText(selectedItemsLabel, headerItemsChanged );

}

void PlotMSIterateTab::setArrowUpNormal(){
	Bool isActive = false;
	setButtonIcon(addItemsButton,isActive);
}

void PlotMSIterateTab::setArrowUpActive(){
	Bool isActive = true;
	setButtonIcon(addItemsButton,isActive);
}

void PlotMSIterateTab::setArrowDownNormal(){
	Bool isActive = false;
	setButtonIcon(removeItemsButton,isActive);
}

void PlotMSIterateTab::setArrowDownActive(){
	Bool isActive = true;
	setButtonIcon(removeItemsButton,isActive);
}

void PlotMSIterateTab::setButtonIcon(QPushButton *button, Bool isActive){
	Bool useArrowUpIcon = false;
	if ( button == addItemsButton ) useArrowUpIcon = true;
	else if  ( button == removeItemsButton ) useArrowUpIcon = false;
	else return;

	auto & icon = useArrowUpIcon ? arrowUpIcon : arrowDownIcon;

	Bool haveSelection;
	if ( button == addItemsButton )
		haveSelection = availableItemsTableView->selectionModel()->hasSelection();
	else
		haveSelection = selectedItemsTableView->selectionModel()->hasSelection();

	auto iconMode = isActive ? QIcon::Mode::Active : QIcon::Mode::Normal;
	auto iconState = isActive && haveSelection ?
			         QIcon::State::On : QIcon::State::Off;

	auto iconSize = icon.actualSize(QSize(32,32),iconMode,iconState);
	auto iconPixmap = icon.pixmap(iconSize,iconMode,iconState);

	button->setIcon(QIcon(iconPixmap));
}

} //namespace

