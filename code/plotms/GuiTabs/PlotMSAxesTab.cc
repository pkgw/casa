//# PlotMSAxesTab.cc: Plot tab for axes parameters.
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
#include <plotms/GuiTabs/PlotMSAxesTab.qo.h>

#include <casaqt/QtUtilities/QtUtilities.h>
#include <plotms/Gui/PlotMSAxisWidget.qo.h>
#include <plotms/Plots/PlotMSPlot.h>
#include <plotms/Plots/PlotMSPlotParameterGroups.h>

#include <QDebug>
#include <QLayout>

using namespace casacore;
namespace casa {

///////////////////////////////
// PLOTMSAXESTAB DEFINITIONS //
///////////////////////////////

PlotMSAxesTab::PlotMSAxesTab(PlotMSPlotTab* plotTab, PlotMSPlotter* parent) :
        PlotMSPlotSubtab(plotTab, parent) {
    setupUi(this);

    // Setup x-axis widget.
    itsXWidget_ = new PlotMSAxisWidget(PMS::DEFAULT_XAXIS, X_BOTTOM | X_TOP);
    itsXWidget_->axisLabel()->setText("X Axis:");
    itsXWidget_->insertLabelDefaults( itsLabelDefaults_ );
    if (xFrame->layout() != NULL) delete xFrame->layout();
    QVBoxLayout* l = new QVBoxLayout();
    l->setContentsMargins(0, 0, 0, 0);
    l->setSpacing(0);
    l->addWidget(itsXWidget_,Qt::AlignTop);
    l->setSizeConstraint(QLayout::SizeConstraint::SetMaximumSize);
    xFrame->setLayout(l);

    // Initially, hide multiple y-axis support.
    setMultipleAxesYEnabled();

    // Connect widgets.
    connect(itsXWidget_, SIGNAL(axisChanged()), SIGNAL(changed()));
    connect(noneRadio, SIGNAL(toggled(bool)), SLOT(overlayChanged()));
    connect(atmRadio, SIGNAL(toggled(bool)), SLOT(overlayChanged()));
    connect(tskyRadio, SIGNAL(toggled(bool)), SLOT(overlayChanged()));    
    connect(imageSbCheckBox, SIGNAL(toggled(bool)), SLOT(overlayChanged()));    
    connect(addYButton, SIGNAL(clicked()), this, SLOT(addYWidget()));
    connect(removeYButton, SIGNAL(clicked()), this, SLOT(removeYWidget()));
    connect(yAxisCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(yAxisSelected(int)));

    //Add an initial y-widget.
    addYWidget();
}

void PlotMSAxesTab::setMultipleAxesYEnabled(){
    bool multipleYAxes = false;
    if ( itsYWidgets_.size() > 1 ){
        multipleYAxes = true;
    }
    yAxisLabel->setVisible( multipleYAxes );
    yAxisCombo->setVisible( multipleYAxes );
    removeYButton->setVisible( multipleYAxes );
}

void PlotMSAxesTab::yAxisSelected( int index ){
    //Remove the current y-widget from the display and
    //put in the newly selected one.
    QLayout* yAxisLayout = yAxisFrame->layout();
    if ( yAxisLayout == NULL ){
        yAxisLayout = new QVBoxLayout();
        yAxisLayout->setContentsMargins(0, 0, 0, 0);
    } else {
        QLayoutItem* layoutItem = yAxisLayout->takeAt( 0 );
        if ( layoutItem != NULL ){
            QWidget* widget = layoutItem->widget();
            if ( widget != NULL ){
                widget->setParent( NULL );
            }
            delete layoutItem;
        }
    }
    yAxisLayout->addWidget( itsYWidgets_[index]);
    yAxisFrame->setLayout( yAxisLayout );
}

void PlotMSAxesTab::overlayChanged() {
    imageSbCheckBox->setEnabled(atmRadio->isChecked() || tskyRadio->isChecked());
    emit changed();
}

void PlotMSAxesTab::removeYWidget(){

    //Remove the item from the combo
    int removeIndex = yAxisCombo->currentIndex();
    yAxisCombo->removeItem( removeIndex );

    //Remove the widget
    PlotMSAxisWidget* removeWidget = itsYWidgets_[removeIndex];
    itsYWidgets_.removeAt( removeIndex );

    //Reset the selected y-widget.
    int currentIndex = yAxisCombo->currentIndex();
    yAxisSelected( currentIndex );

    setMultipleAxesYEnabled();
    emit yAxisIdentifierRemoved( removeIndex );
    delete removeWidget;
    emit changed();
}

void PlotMSAxesTab::setYAxisLabel( PlotMSAxisWidget* yWidget, int /*index*/ ){
    QString identifier = yWidget->getIdentifier();
    yAxisCombo->addItem( identifier );
}

void PlotMSAxesTab::addYWidget(){
    PlotMSAxisWidget* yWidget = new PlotMSAxisWidget(PMS::DEFAULT_YAXIS, Y_LEFT | Y_RIGHT );
    int index = itsYWidgets_.size();
    yWidget->insertLabelDefaults( itsLabelDefaults_ );
    itsYWidgets_.append( yWidget );
    setYAxisLabel( yWidget, index);
    yAxisCombo->setCurrentIndex( index );

    setMultipleAxesYEnabled();

    connect(yWidget, SIGNAL(axisChanged()), SIGNAL(changed()));
    connect(yWidget, SIGNAL(axisIdentifierChanged(PlotMSAxisWidget*)),
            this, SLOT(axisIdentifierChanged(PlotMSAxisWidget*)));
    //Because we may not have been able to connect to listeners when we
    //had just one y-axis, let them know the first one exists, too.
    if ( index <= 1 ){
        emit yAxisIdentifierChanged( 0, itsYWidgets_[0]->getIdentifier());
    }
    emit yAxisIdentifierChanged( index, yWidget->getIdentifier());
    emit changed();
}

PlotMSAxesTab::~PlotMSAxesTab() {
    while( ! itsYWidgets_.isEmpty() ){
        PlotMSAxisWidget* widget = itsYWidgets_.takeAt(0 );
        widget->setParent( NULL );
        delete widget;
    }
}

bool PlotMSAxesTab::isAxesValid() const {
    bool axesValid = true;
    QList<PlotMSAxisWidget*> uniqueYWidgets;
    bool found = false;
    for ( int i = 0; i < itsYWidgets_.size(); i++ ){
        found = false;
        for ( int j = 0; j < uniqueYWidgets.size(); j++ ){
            if ( uniqueYWidgets[j]->matchesData( itsYWidgets_[i])){
                found = true;
                break;
            }
        }
        if ( found ){
            break;
        } else {
            uniqueYWidgets.append( itsYWidgets_[i]);
        }
    }
    if ( found ){
        axesValid = false;
    }
    return axesValid;
}


void PlotMSAxesTab::getValue(PlotMSPlotParameters& params) const {
    PMS_PP_Cache* c = params.typedGroup<PMS_PP_Cache>();
    PMS_PP_Axes* a = params.typedGroup<PMS_PP_Axes>();
    PMS_PP_MSData* d = params.typedGroup<PMS_PP_MSData>();

    if (c == NULL) {
        params.setGroup<PMS_PP_Cache>();
        c = params.typedGroup<PMS_PP_Cache>();
    }
    if (a == NULL) {
        params.setGroup<PMS_PP_Axes>();
        a = params.typedGroup<PMS_PP_Axes>();
    }
    if (d == NULL) {
        params.setGroup<PMS_PP_MSData>();
        d = params.typedGroup<PMS_PP_MSData>();
    }
    
    //The cache must have exactly as many x-axes as y-axes so we duplicate
    //the x-axis properties here.
    int yAxisCount = itsYWidgets_.size();
    PMS::Axis xAxis = itsXWidget_->axis();
    c->resize( yAxisCount );
    a->resize( yAxisCount );
    for ( int i = 0; i < yAxisCount; i++ ){
        c->setXAxis(xAxis, itsXWidget_->data(), i);
        c->setXInterp(itsXWidget_->interpMethod(),i);
        c->setXFrame(itsXWidget_->refFrame(),i);
        a->setXAxis(itsXWidget_->attachAxis(), i);
        a->setXRange(itsXWidget_->rangeCustom(), itsXWidget_->range(), i);
    }

    for ( int i = 0; i < yAxisCount; i++ ){
        c->setYAxis(itsYWidgets_[i]->axis(), itsYWidgets_[i]->data(), i);
        c->setYInterp(itsYWidgets_[i]->interpMethod(),i);
        c->setYFrame(itsYWidgets_[i]->refFrame(),i);
        a->setYAxis(itsYWidgets_[i]->attachAxis(), i);
        a->setYRange(itsYWidgets_[i]->rangeCustom(), itsYWidgets_[i]->range(), i);
    }

    bool showatm(atmRadio->isChecked());
    bool showtsky(tskyRadio->isChecked());
    bool overlay(showatm || showtsky);
    bool showimage(imageSbCheckBox->isChecked() && overlay);
    c->setShowAtm(showatm);
    c->setShowTsky(showtsky);
    c->setShowImage(showimage);
    if (overlay) {
        // add ATM/TSKY/IMAGESB yaxis "under the hood" for GUI client
        if (xAxis==PMS::CHANNEL || 
            xAxis==PMS::FREQUENCY) { 
            PMS::Axis overlayAxis = (showatm ? PMS::ATM : PMS::TSKY);
            bool foundOverlayAxis(false), foundImageAxis(false);
            const vector<PMS::Axis> yAxes = c->yAxes();
            for (uInt i=0; i<yAxes.size(); ++i) {
                if (yAxes[i] == overlayAxis) {
                    foundOverlayAxis = True;
                } else if (yAxes[i] == PMS::IMAGESB) {
                    foundImageAxis = True;
                }
            }
            if (!foundOverlayAxis) {
                // add axis to Cache axes
                int index = c->numXAxes();
                c->setAxes(xAxis, overlayAxis, c->xDataColumn(0), PMS::DEFAULT_DATACOLUMN, index);
                // set Axes positions
                a->resize(index+1, true);  // copy values for index 0
                a->setAxes(a->xAxis(index-1), Y_RIGHT, index);
                // keep same xaxis range
                a->setXRange(itsXWidget_->rangeCustom(), itsXWidget_->range(), index);
                // set Display symbol color
                PMS_PP_Display* disp = params.typedGroup<PMS_PP_Display>();
                PlotSymbolPtr overlaySymbol = disp->unflaggedSymbol(index);
                overlaySymbol->setSymbol("circle");
                overlaySymbol->setColor("#FF00FF"); // magenta
                disp->setUnflaggedSymbol(overlaySymbol, index);
                PlotSymbolPtr flaggedSymbol = disp->flaggedSymbol();
                disp->setFlaggedSymbol(flaggedSymbol, index);
            }
            if (showimage && !foundImageAxis) {
                // add axis to Cache axes
                int index = c->numXAxes();
                c->setAxes(xAxis, PMS::IMAGESB, c->xDataColumn(0), PMS::DEFAULT_DATACOLUMN, index);
                // set Axes positions
                a->resize(index+1, true);  // copy values for index 0
                a->setAxes(a->xAxis(index-1), Y_RIGHT, index);
                // keep same xaxis range
                a->setXRange(itsXWidget_->rangeCustom(), itsXWidget_->range(), index);
                // set Display symbol color
                PMS_PP_Display* disp = params.typedGroup<PMS_PP_Display>();
                PlotSymbolPtr imageSymbol = disp->unflaggedSymbol(index);
                imageSymbol->setSymbol("circle");
                imageSymbol->setColor("#000000"); // black
                disp->setUnflaggedSymbol(imageSymbol, index);
                PlotSymbolPtr flaggedSymbol = disp->flaggedSymbol();
                disp->setFlaggedSymbol(flaggedSymbol, index);
            }
        }
    }
}

void PlotMSAxesTab::setValue(const PlotMSPlotParameters& params) {
    const PMS_PP_Cache* c = params.typedGroup<PMS_PP_Cache>();
    const PMS_PP_Axes* a = params.typedGroup<PMS_PP_Axes>();
    if(c == NULL || a == NULL)
        return; // shouldn't happen

    // X Widget
    PMS::Axis cacheAxis = c->xAxis();
    PMS::DataColumn cacheColumn =  c->xDataColumn();
    PlotAxis axesAxis = a->xAxis();
    bool axesXRangeSet = a->xRangeSet();
    std::pair<double, double> axesXRange = a->xRange();
    itsXWidget_->setValue(cacheAxis, cacheColumn, axesAxis, axesXRangeSet, axesXRange);
    const auto & xFrames = c->xFrames();
    const auto & xInterps = c->xInterps();
    itsXWidget_->setDirParams(xInterps[0],xFrames[0]);

    // Radio Buttons
    bool atm(c->showAtm()), tsky(c->showTsky()), overlay(atm || tsky);
    atmRadio->setChecked(atm);
    tskyRadio->setChecked(tsky);
    noneRadio->setChecked(!overlay);
    imageSbCheckBox->setEnabled(overlay);
    imageSbCheckBox->setChecked(c->showImage());

    // Y Widgets
    int yAxisCount = a->numYAxes();
    int yWidgetSize = itsYWidgets_.size();
    int minValue = qMin( yAxisCount, yWidgetSize );
    const auto & yFrames = c->yFrames();
    int nFramesY = yFrames.size();
    const auto & yInterps = c->yInterps();
    int nInterpsY = yInterps.size();
    for ( int i = 0; i < minValue; i++ ){
        itsYWidgets_[i]->setValue(c->yAxis(i), c->yDataColumn(i), a->yAxis(i),
            a->yRangeSet(i), a->yRange(i));
        if (i < nInterpsY and i < nFramesY ) {
            itsYWidgets_[i]->setDirParams(yInterps[i],yFrames[i]);
        } else {
            throw AipsError("PlotMS internal error.",AipsError::BOUNDARY);
        }
    }
}

void PlotMSAxesTab::axisIdentifierChanged(PlotMSAxisWidget* axisWidget) {

    int yIndex = itsYWidgets_.indexOf( axisWidget );
    if ( yIndex >= 0 ){
        QString newIdentifier = axisWidget->getIdentifier();
        yAxisCombo->setItemText(yIndex, newIdentifier );
        emit yAxisIdentifierChanged( yIndex, newIdentifier );
    }
}

void PlotMSAxesTab::update(const PlotMSPlot& plot) {
    const PlotMSPlotParameters& params = plot.parameters();
    PlotMSPlotParameters newParams(params);

    // Deal somehow with data structures inconsistencies
    // between the Cache and the GUI
    getValue(newParams);

    // Shortcuts
    const PMS_PP_MSData* data = params.typedGroup<PMS_PP_MSData>();
    const PMS_PP_Cache* c = params.typedGroup<PMS_PP_Cache>(),
                    *cNew = newParams.typedGroup<PMS_PP_Cache>();
    const PMS_PP_Axes* a = params.typedGroup<PMS_PP_Axes>(),
                   *aNew = newParams.typedGroup<PMS_PP_Axes>();

    // Shouldn't happen
    if(data == NULL || c == NULL || cNew == NULL || a == NULL || aNew == NULL)
        throw(AipsError("PlotMSAxesTab::update(): internal error"));

    // Update XWidget's "in cache" checkbox
    auto loadedAxes = plot.cache().loadedAxes();
    auto newXAxis = cNew->xAxis(0);
    auto newXAxisIsLoaded = std::find(loadedAxes.begin(),loadedAxes.end(),
        newXAxis) != loadedAxes.end();
    if (PMS::axisIsRaDec(newXAxis)) {
        auto newXDirParams = cNew->xDirectionParams(0);
        newXAxisIsLoaded = plot.cache().areRaDecAxesLoaded(newXDirParams);
    }
    itsXWidget_->setInCache(newXAxisIsLoaded);

    // Update YWidgets "in cache" checkbox
    int yAxisCount = itsYWidgets_.size();
    for ( int j = 0; j < yAxisCount; j++ ){
        auto newYAxis = cNew->yAxis(j);
        auto newYAxisIsLoaded = std::find(loadedAxes.begin(),loadedAxes.end(),
            newYAxis) != loadedAxes.end();
        if (PMS::axisIsRaDec(newYAxis)) {
            auto newYDirParams = cNew->yDirectionParams(j);
            newYAxisIsLoaded = plot.cache().areRaDecAxesLoaded(newYDirParams);
        }
        itsYWidgets_[j]->setInCache(newYAxisIsLoaded);
    }

    // Highlight XWidget changes
    if (data->isSet()){
        auto xAxisChanged = cNew->xAxis(0) != c->xAxis(0);
        auto xDataAxisParamsChanged = cNew->xDataColumn(0) != c->xDataColumn(0);
        auto xDirectionAxisParamsChanged = cNew->xDirectionParams(0) != c->xDirectionParams(0);
        auto xAxisParamsChanged = xDataAxisParamsChanged || xDirectionAxisParamsChanged;
        auto xInterpChanged = cNew->xInterp(0) != c->xInterp(0);
        auto xRefFrameChanged = cNew->xFrame(0) != c->xFrame(0);
        auto xAxisLocationChanged = aNew->xAxis(0) != a->xAxis(0);
        auto xRangeModeChanged = aNew->xRangeSet(0) != a->xRangeSet(0);
        auto xRangeChanged = aNew->xRange(0) != a->xRange(0);
        highlightWidgetText(itsXWidget_->axisLabel(), xAxisChanged || xAxisParamsChanged);
        highlightWidgetText(itsXWidget_->dataLabel(), xDataAxisParamsChanged);
        highlightWidgetText(itsXWidget_->interpLabel(), xInterpChanged);
        highlightWidgetText(itsXWidget_->refFrameLabel(), xRefFrameChanged);
        highlightWidgetText(itsXWidget_->attachLabel(), xAxisLocationChanged);
        highlightWidgetText(itsXWidget_->rangeLabel(), xRangeModeChanged || xRangeChanged );
    }
    // Highlight overlay changes
    bool atmRequested(atmRadio->isChecked()), tskyRequested(tskyRadio->isChecked());
    bool overlayChanged = (c->showAtm() != atmRequested) || (c->showTsky() != tskyRequested);
    bool sidebandChanged = ((atmRequested || tskyRequested) && (c->showImage() != imageSbCheckBox->isChecked()));
    highlightWidgetText(overlayLabel, overlayChanged || sidebandChanged);

    // Highlight YWidgets changes
    if (data->isSet()){
        for ( int i = 0; i < yAxisCount; i++ ){
            QLabel* axisLabel = itsYWidgets_[i]->axisLabel();
            auto yAxisChanged = c->yAxis(i) != cNew->yAxis(i);
            auto yDataAxisParamsChanged = PMS::axisIsData(c->yAxis(i)) &&
                (c->yDataColumn(i) != cNew->yDataColumn(i));
            auto yDirectionAxisParamsChanged = PMS::axisIsRaDec(c->yAxis(i)) &&
                (c->yDirectionParams(i) != cNew->yDirectionParams(i));
            auto yAxisParamsChanged = yDataAxisParamsChanged || yDirectionAxisParamsChanged;
            highlightWidgetText(axisLabel, yAxisChanged || yAxisParamsChanged);
            highlightWidgetText(itsYWidgets_[i]->dataLabel(),
                c->yDataColumn(i) != cNew->yDataColumn(i));
            highlightWidgetText(itsYWidgets_[i]->interpLabel(),
                c->yInterp(i) != cNew->yInterp(i));
            highlightWidgetText(itsYWidgets_[i]->refFrameLabel(),
                c->yFrame(i) != cNew->yFrame(i));
            highlightWidgetText(itsYWidgets_[i]->attachLabel(),
                a->yAxis() != aNew->yAxis());
            highlightWidgetText(itsYWidgets_[i]->rangeLabel(),
                (a->yRangeSet(i) != aNew->yRangeSet(i)) ||
                (a->yRangeSet(i) && (a->yRange(i) != aNew->yRange(i))));
        }
    }

    // If the user hasn't set a custom range, set defaults for time axis.
    // For time axis, get bounds from cache; else use 0 for other axes types
    pair<Double,Double> timebounds;
    bool isDate;

    if (!itsXWidget_->rangeCustom()) {
        isDate = (cNew->xAxis() == PMS::TIME);
        if (isDate) {
            timebounds = plot.cache().getTimeBounds();
            itsXWidget_->setRange(isDate, timebounds.first, timebounds.second); 
        } else {
            itsXWidget_->setRange(isDate, 0.0, 0.0);
        }
    }

    for (int i = 0; i < yAxisCount; i++ ){
        if (!itsYWidgets_[i]->rangeCustom()) {
            PMS::Axis yData = cNew->yAxis(i);
            isDate = (yData == PMS::TIME);
            if (isDate) {
                timebounds = plot.cache().getTimeBounds();
                itsYWidgets_[i]->setRange(isDate, timebounds.first, timebounds.second); 
                break;
            } else {
                itsYWidgets_[i]->setRange(isDate, 0.0, 0.0);
            }
        }
    }
}

}
