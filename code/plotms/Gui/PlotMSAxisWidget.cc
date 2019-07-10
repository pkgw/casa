//# PlotMSAxisWidget.cc: Widget for choosing a single axis.
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
#include <plotms/Gui/PlotMSAxisWidget.qo.h>
#include <plotms/Gui/PlotMSAxisWidget.ui.h>


#include <plotms/Gui/PlotRangeWidget.qo.h>
#include <plotms/GuiTabs/PlotMSTab.qo.h>
#include <QDebug>
#include <QLayout>

using namespace casacore;
namespace casa {

//////////////////////////////////
// PLOTMSAXISWIDGET DEFINITIONS //
//////////////////////////////////

// Constructors/Destructors //

PlotMSAxisWidget::PlotMSAxisWidget(PMS::Axis defaultAxis, int attachAxes,
		QWidget* parent) : QtEditingWidget(parent) {
	setupUi(this);


	// Setup axes choices
	const vector<String>& axes = PMS::axesStrings();
	String def = PMS::axis(defaultAxis);
	// Hide the last axes (overlays) from the user....
	for (int i = 0; i < PMS::ATM; i++) {
		chooser->addItem(axes[i].c_str());
		if(axes[i] == def) chooser->setCurrentIndex(i);
	}

	// Setup attach axes.
	initPlotAxis( attachAxes );

	// Setup data column choices.
	const vector<String>& data = PMS::dataColumnStrings();
	def = PMS::dataColumn(PMS::DEFAULT_DATACOLUMN);

	for(unsigned int i = 0; i < data.size(); i++) {
		dataChooser->addItem(data[i].c_str());
		if(data[i] == def){
			dataChooser->setCurrentIndex(i);
		}
	}
	// Setup direction parameters choices.
	const auto & interpStrings = PMS::interpMethodStrings();
	auto defaultInterp = PMS::interpMethod(PMS::DEFAULT_INTERPMETHOD);
	for(unsigned int i = 0; i < interpStrings.size(); i++) {
		interpChooser->addItem(interpStrings[i].c_str());
		if(interpStrings[i] == defaultInterp){
			interpChooser->setCurrentIndex(i);
		}
	}
	const auto & refFrameStrings = PMS::coordSystemStrings();
	auto defaultRefFrame = PMS::coordSystem(PMS::DEFAULT_COORDSYSTEM);
	for(unsigned int i = 0; i < refFrameStrings.size(); i++) {
		refFrameChooser->addItem(refFrameStrings[i].c_str());
		if(refFrameStrings[i] == defaultRefFrame){
			refFrameChooser->setCurrentIndex(i);
		}
	}

	// Setup range widget.
	itsRangeWidget_ = new PlotRangeWidget(true);
	itsRangeWidget_->setFixedHeight(itsRangeWidget_->height());

	if (rangeWidgetFrame->layout() != nullptr)
		delete rangeWidgetFrame->layout();
	{
		auto* l = new QVBoxLayout();
		l->setContentsMargins(0, 0, 0, 0);
		l->setSpacing(0);
		l->addWidget(itsRangeWidget_,Qt::AlignTop);
		l->setAlignment(Qt::AlignTop);
		rangeWidgetFrame->setLayout(l);
		rangeFrame->layout()->setAlignment(rangeWidgetFrame,Qt::AlignTop);
		rangeFrame->setFixedHeight(itsRangeWidget_->height());
	}

	setAutoFillBackground( true );
	QPalette pal = palette();
	QColor bgColor( "#F0F0F0" );
	pal.setColor( QPalette::Background, bgColor );
	setPalette( pal );

	axisChanged(chooser->currentText());

	// Connect widgets.
	connect(chooser, SIGNAL(currentIndexChanged(const QString&)),
			SLOT(axisChanged(const QString&)));
	connect( dataChooser, SIGNAL(currentIndexChanged(const QString&)),
			SLOT(axisDataChanged()));
	connect( interpChooser, SIGNAL(currentIndexChanged(const QString&)),
			SLOT(axisInterpChanged()));
	connect( refFrameChooser, SIGNAL(currentIndexChanged(const QString&)),
			SLOT(axisRefFrameChanged()));

	connect(chooser, SIGNAL(currentIndexChanged(int)), SIGNAL(axisChanged()));
	connect(dataChooser, SIGNAL(currentIndexChanged(int)), SIGNAL(axisChanged()));
	connect(interpChooser, SIGNAL(currentIndexChanged(int)), SIGNAL(axisChanged()));
	connect(refFrameChooser, SIGNAL(currentIndexChanged(int)), SIGNAL(axisChanged()));

	if ( attachBottom != NULL ){
		connect(attachBottom, SIGNAL(toggled(bool)), SIGNAL(axisChanged()));
	}
	if ( attachTop != NULL ){
		connect(attachTop, SIGNAL(toggled(bool)), SIGNAL(axisChanged()));
	}
	if ( attachLeft != NULL ){
		connect(attachLeft, SIGNAL(toggled(bool)), SIGNAL(changed()));
	}
	if ( attachRight != NULL ){
		connect(attachRight, SIGNAL(toggled(bool)), SIGNAL(changed()));
	}
	connect(itsRangeWidget_, SIGNAL(changed()), SIGNAL(changed()));

}

PlotMSAxisWidget::~PlotMSAxisWidget() { }

void PlotMSAxisWidget::initPlotAxis(int attachAxes){
	attachLeft->setVisible(attachAxes & Y_LEFT);
	attachRight->setVisible(attachAxes & Y_RIGHT);
	attachBottom->setVisible(attachAxes & X_BOTTOM);
	attachTop->setVisible(attachAxes & X_TOP);
	if ( attachAxes & Y_LEFT ){
		setAttachAxis( Y_LEFT );
	}
	else if ( attachAxes & Y_RIGHT ){
		setAttachAxis( Y_RIGHT );
	}
	else if ( attachAxes & X_BOTTOM ){
		setAttachAxis( X_BOTTOM );
	}
	else {
		setAttachAxis( X_TOP );
	}

}

void PlotMSAxisWidget::insertLabelDefaults( QMap<QLabel*,QString>& map ){
	map.insert(AxisWidget::axisLabel, AxisWidget::axisLabel->text());
	map.insert(AxisWidget::dataLabel,AxisWidget::dataLabel->text());
	map.insert(AxisWidget::attachLabel,AxisWidget::attachLabel->text());
	map.insert(AxisWidget::rangeLabel, AxisWidget::rangeLabel->text());
}

// Public Methods //
PMS::Axis PlotMSAxisWidget::axis() const {
    return PMS::axis(chooser->currentText().toStdString());
}

PMS::DataColumn PlotMSAxisWidget::data() const {
	QString dataText = dataChooser->currentText();
	return PMS::dataColumn(dataText.toStdString());
}

PMS::InterpMethod PlotMSAxisWidget::interpMethod() const {
	QString interpText = interpChooser->currentText();
	return PMS::interpMethod(interpText.toStdString());
}
PMS::CoordSystem PlotMSAxisWidget::refFrame() const {
	QString refFrameText = refFrameChooser->currentText();
	return PMS::coordSystem(refFrameText.toStdString());
}

QString PlotMSAxisWidget::getIdentifier() const {
	QString axisText = chooser->currentText();
	QString dataText = dataChooser->currentText();
	if ( dataText != "data"){
		axisText = axisText + ": "+dataText;
	}
	return axisText;
}

bool PlotMSAxisWidget::matchesData(const PlotMSAxisWidget* other ) const {
	bool matchingData = false;
	if ( other != NULL ){
		if ( other->axis() == this->axis() ){
			if ( other->data() == this->data() ){
				matchingData = true;
			}
		}
	}
	return matchingData;
}

PlotAxis PlotMSAxisWidget::attachAxis() const {
	PlotAxis plotAxis = X_BOTTOM;
	if ( attachTop != NULL && attachTop->isChecked()){
		plotAxis = X_TOP;
	}
	else if ( attachLeft != NULL && attachLeft->isChecked()){
		plotAxis = Y_LEFT;
	}
	else if ( attachRight != NULL && attachRight->isChecked()){
		plotAxis = Y_RIGHT;
	}
	return plotAxis;
}



bool PlotMSAxisWidget::rangeCustom() const {
    return itsRangeWidget_->isCustom();
}

prange_t PlotMSAxisWidget::range() const {
    return itsRangeWidget_->getRange();
}

void PlotMSAxisWidget::setRange(bool isDate, double from, double to) {
	itsRangeWidget_->setRange(isDate, false, from, to); 
}

void PlotMSAxisWidget::setAttachAxis(PlotAxis attachAxis ){
	switch(attachAxis) {
	case X_BOTTOM:
	    if ( attachBottom != NULL ){
	    	attachBottom->setChecked(true);
	    }
	    break;
	case X_TOP:
	    if ( attachTop != NULL ){
	    	attachTop->setChecked(true);
	    }
	    break;
	case Y_LEFT:
	    if ( attachLeft != NULL ){
	    	attachLeft->setChecked(true);
	    }
	    break;
	case Y_RIGHT:
	    if ( attachRight != NULL ){
	    	attachRight->setChecked(true);
	    }
	    break;
	default:
	    break;
	}

}

void PlotMSAxisWidget::setValue(PMS::Axis axis, PMS::DataColumn data,
        PlotAxis attachAxis, bool rangeCustom, prange_t range){
    PlotMSTab::setChooser(chooser, PMS::axis(axis));
    PlotMSTab::setChooser(dataChooser, PMS::dataColumn(data));
    setAttachAxis( attachAxis );
    itsRangeWidget_->setRange(PMS::axisType(axis) == PMS::TTIME, rangeCustom,
            range);
}

void PlotMSAxisWidget::setDirParams(PMS::InterpMethod interp, PMS::CoordSystem refFrame){
	PlotMSTab::setChooser(interpChooser, PMS::interpMethod(interp));
	PlotMSTab::setChooser(refFrameChooser, PMS::coordSystem(refFrame));
}

void PlotMSAxisWidget::setInCache(bool isInCache) {
    inCache->setChecked(isInCache);
}

// Private Slots //
void PlotMSAxisWidget::axisChanged(const QString& value) {
	PMS::Axis currAxis=PMS::axis(value.toStdString());

	// Show Data Column chooser, if necessary
	bool dataAxis = PMS::axisIsData(currAxis);
	AxisWidget::dataLabel->setVisible( dataAxis );
	dataChooser->setVisible( dataAxis );

	// Show Direction parameters, if necessary
	bool needDirParams = PMS::axisIsRaDec(currAxis);
	AxisWidget::dirParamsFrame->setVisible(needDirParams);

	emit axisIdentifierChanged(this);
}

void PlotMSAxisWidget::axisDataChanged(){
	emit axisIdentifierChanged( this );
}

void PlotMSAxisWidget::axisInterpChanged(){
	emit axisIdentifierChanged( this );
}

void PlotMSAxisWidget::axisRefFrameChanged(){
	emit axisIdentifierChanged( this );
}

}
