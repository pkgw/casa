//# PlotMSLabelFormat.cc: Class for generating labels based on a format.
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
#include <plotms/PlotMS/PlotMSLabelFormat.h>

#include <plotms/Plots/PlotMSPlotParameterGroups.h>

using namespace casacore;
namespace casa {

///////////////////////////////////
// PLOTMSLABELFORMAT DEFINITIONS //
///////////////////////////////////

// Static //

const String& PlotMSLabelFormat::TAGSEPARATOR() {
    static const String str("%%"); return str; }

const String& PlotMSLabelFormat::TAG_AXIS() {
    static const String str("axis"); return str; }
const String& PlotMSLabelFormat::TAG_XAXIS() {
    static const String str("xaxis"); return str; }
const String& PlotMSLabelFormat::TAG_YAXIS() {
    static const String str("yaxis"); return str; }

const String& PlotMSLabelFormat::TAG_UNIT() {
    static const String str("unit"); return str; }
const String& PlotMSLabelFormat::TAG_XUNIT() {
    static const String str("xunit"); return str; }
const String& PlotMSLabelFormat::TAG_YUNIT() {
    static const String str("yunit"); return str; }

const String& PlotMSLabelFormat::TAG_IF_UNIT() {
    static const String str("ifunit"); return str; }
const String& PlotMSLabelFormat::TAG_IF_XUNIT() {
    static const String str("ifxunit"); return str; }
const String& PlotMSLabelFormat::TAG_IF_YUNIT() {
    static const String str("ifyunit"); return str; }
const String& PlotMSLabelFormat::TAG_ENDIF_UNIT() {
    static const String str("endifunit"); return str; }
const String& PlotMSLabelFormat::TAG_ENDIF_XUNIT() {
    static const String str("endifxunit"); return str; }
const String& PlotMSLabelFormat::TAG_ENDIF_YUNIT() {
    static const String str("endifyunit"); return str; }

const String& PlotMSLabelFormat::TAG_REFVALUE() {
    static const String str("refvalue"); return str; }
const String& PlotMSLabelFormat::TAG_XREFVALUE() {
    static const String str("xrefvalue"); return str; }
const String& PlotMSLabelFormat::TAG_YREFVALUE() {
    static const String str("yrefvalue"); return str; }

const String& PlotMSLabelFormat::TAG_IF_REFVALUE() {
    static const String str("ifrefvalue"); return str; }
const String& PlotMSLabelFormat::TAG_IF_XREFVALUE() {
    static const String str("ifxrefvalue"); return str; }
const String& PlotMSLabelFormat::TAG_IF_YREFVALUE() {
    static const String str("ifyrefvalue"); return str; }
const String& PlotMSLabelFormat::TAG_ENDIF_REFVALUE() {
    static const String str("endifrefvalue"); return str; }
const String& PlotMSLabelFormat::TAG_ENDIF_XREFVALUE() {
    static const String str("endifxrefvalue"); return str; }
const String& PlotMSLabelFormat::TAG_ENDIF_YREFVALUE() {
    static const String str("endifyrefvalue"); return str; }

String PlotMSLabelFormat::TAG(const String& tag) {
    return TAGSEPARATOR() + tag + TAGSEPARATOR(); }

void PlotMSLabelFormat::addDataToTag( String& tag, PMS::Axis axis, PMS::DataColumn column ){
	if ( PMS::axisIsData(axis) ){
	    if ( column != PMS::DATA ){
	    	String axisData = PMS::dataColumn(column);
	        tag = tag + ":"+axisData;
	    }
	}
}

void PlotMSLabelFormat::addPolnRatioToTag( String& tag, PMS::Axis axis ) {
    switch (axis) {
        case PMS::AMP:
        case PMS::REAL:
        case PMS::IMAG:
        case PMS::GAMP:
        case PMS::GREAL:
        case PMS::GIMAG:
        case PMS::TSYS:
        case PMS::SWP:
        case PMS::SNR:
            tag += " POLN Ratio";
            break;
        case PMS::PHASE:
        case PMS::GPHASE:
        case PMS::DELAY:
            tag += " POLN Difference";
            break;
        case PMS::CORR:
            tag += " Ratio";
        default:
            break;
    }
}

String PlotMSLabelFormat::getLabel(const String& format, PMS::Axis axis,
            PMS::Axis xAxis, vector<PMS::Axis> yAxes, bool refValueSet,
            double refValue, bool xRefValueSet, double xRefValue,
            vector<bool> yRefValueSets, vector<double> yRefValues,
            PMS::DataColumn xData, const vector<PMS::DataColumn>& yData,
            bool polnRatio) {
    stringstream ss;
    
    int yAxisCount = yAxes.size();
    PMS::AxisUnit unit = PMS::axisUnit(axis);
    PMS::AxisUnit xUnit = PMS::axisUnit(xAxis);
    vector<PMS::AxisUnit> yUnits;
    for ( int i = 0; i < yAxisCount; i++ ){
    	yUnits.push_back(PMS::axisUnit(yAxes[i]));
    }
    PlotAxisScale scale = PMS::axisScale(axis);
    PlotAxisScale xScale = PMS::axisScale(xAxis);
    vector<PlotAxisScale> yScales;
    for ( int i = 0; i < yAxisCount; i++ ){
    	yScales.push_back( PMS::axisScale(yAxes[i]) );
    }
    bool isDate = scale == DATE_MJ_DAY || scale == DATE_MJ_SEC,
        xIsDate = xScale == DATE_MJ_DAY || xScale == DATE_MJ_SEC;
    vector<bool> yIsDates;
    for ( int i = 0; i < yAxisCount; i++ ){
        yIsDates.push_back( yScales[i] == DATE_MJ_DAY || yScales[i] == DATE_MJ_SEC );
    }
    
    String tempFormat = format, token, tag;
    bool tokenWasTag, ifUnit = false, ifXUnit = false, ifYUnit = false,
         ifRefValue = false, ifXRefValue = false, ifYRefValue = false;
    
    while(nextToken(tempFormat, token, tokenWasTag)) {
        if(tokenWasTag) {
            tag = "";

            if(PMS::strEq(token, TAG_AXIS(), true)){
            	tag = PMS::axis(axis);
            	addDataToTag( tag, axis, xData );
                if (polnRatio) addPolnRatioToTag(tag, axis);
            }
            else if(PMS::strEq(token, TAG_XAXIS(), true)){
            	tag=PMS::axis(xAxis);
            	addDataToTag( tag, xAxis, xData );
                if (polnRatio) addPolnRatioToTag(tag, xAxis);
            }
            else if(PMS::strEq(token, TAG_YAXIS(), true)){
            	for( int i = 0; i < yAxisCount; i++ ){
            		tag= tag + PMS::axis(yAxes[i]);
            		addDataToTag( tag, yAxes[i], yData[i]);
                    if (polnRatio) addPolnRatioToTag(tag, yAxes[i]);
            		if ( i != yAxisCount - 1 ){
            			tag = tag + ", ";
            		}
            	}
            }
            else if(PMS::strEq(token, TAG_UNIT(), true)){
            	tag = PMS::axisUnit(unit);
            }
            else if(PMS::strEq(token, TAG_XUNIT(), true)){
            	tag = PMS::axisUnit(xUnit);
            }
            else if(PMS::strEq(token, TAG_YUNIT(), true)){
                tag = PMS::axisUnit(yUnits[0]);
            }
            else if(PMS::strEq(token, TAG_IF_UNIT(), true)){
                ifUnit = true;
            }
            else if(PMS::strEq(token, TAG_IF_XUNIT(), true)){
                ifXUnit = true;
            }
            else if(PMS::strEq(token, TAG_IF_YUNIT(), true)){
                ifYUnit = true;
            }
            else if(PMS::strEq(token, TAG_ENDIF_UNIT(), true)){
                ifUnit = false;
            }
            else if(PMS::strEq(token, TAG_ENDIF_XUNIT(), true)){
                ifXUnit = false;
            }
            else if(PMS::strEq(token, TAG_ENDIF_YUNIT(), true)){
                ifYUnit = false;
            }
            else if(PMS::strEq(token, TAG_REFVALUE(), true)) {
                tag = refValueSet ? (isDate ?
                      Plotter::formattedDateString(REFERENCE_DATE_FORMAT,
                                                   refValue, scale) :
                      String::toString(refValue)) : "";
            }
            else if(PMS::strEq(token, TAG_XREFVALUE(), true)) {
                tag = xRefValueSet ? (xIsDate ?
                      Plotter::formattedDateString(REFERENCE_DATE_FORMAT,
                                                   xRefValue, xScale) :
                      String::toString(xRefValue)) : "";
            } else if(PMS::strEq(token, TAG_YREFVALUE(), true)) {
                tag = yRefValueSets[0] ? (yIsDates[0] ?
                      Plotter::formattedDateString(REFERENCE_DATE_FORMAT,
                                                   yRefValues[0], yScales[0]) :
                      String::toString(yRefValues[0])) : "";
            }
            else if(PMS::strEq(token, TAG_IF_REFVALUE(), true))
                ifRefValue = true;
            else if(PMS::strEq(token, TAG_IF_XREFVALUE(), true))
                ifXRefValue = true;
            else if(PMS::strEq(token, TAG_IF_YREFVALUE(), true))
                ifYRefValue = true;
            else if(PMS::strEq(token, TAG_ENDIF_REFVALUE(), true))
                ifRefValue = false;
            else if(PMS::strEq(token, TAG_ENDIF_XREFVALUE(), true))
                ifXRefValue = false;
            else if(PMS::strEq(token, TAG_ENDIF_YREFVALUE(), true))
                ifYRefValue = false;
            else tag = TAG(token);
        }
        else {
        	tag = token;
        }
        
        if((!ifUnit || unit != PMS::UNONE) && (!ifXUnit || xUnit != PMS::UNONE)
           && (!ifYUnit || yUnits[0] != PMS::UNONE) && (!ifRefValue || refValueSet)
           && (!ifXRefValue || xRefValueSet) && (!ifYRefValue || yRefValueSets[0])){
            ss << tag;
        }
    }
    return ss.str();
}

bool PlotMSLabelFormat::nextToken(String& format, String& token,
        bool& tokenWasTag) {
    if(format.size() == 0) {
        token = "";
        tokenWasTag = false;
        return false;
    }
    
    unsigned int i = format.find(TAGSEPARATOR()), j;
    if(i >= format.size() ||
       (j = format.find(TAGSEPARATOR(), i + 1)) >= format.size()) {
        // no more tags left
        token = format;
        tokenWasTag = false;
        format = "";
        return true;
    }

    if(i == 0) {
        // tag is next token
        token= format.substr(TAGSEPARATOR().size(), j - TAGSEPARATOR().size());
        tokenWasTag = true;
        format = format.substr(j + TAGSEPARATOR().size());        
    } else {
        // text is next token
        token = format.substr(0, i);
        tokenWasTag = false;
        format = format.substr(i);
    }
    
    return true;
}

const String PlotMSLabelFormat::REFERENCE_DATE_FORMAT = "%y/%m/%d";


// Non-Static //

PlotMSLabelFormat::PlotMSLabelFormat(const String& f) : format(f) { }

PlotMSLabelFormat::PlotMSLabelFormat(const PlotMSLabelFormat& copy) {
    operator=(copy); }

PlotMSLabelFormat::~PlotMSLabelFormat() { }

String PlotMSLabelFormat::getLabel( PMS::Axis axis, bool refValueSet,
		double refValue, PMS::DataColumn dataColumn, bool polnRatio ) const{
	vector<PMS::Axis> axes(1);
	axes[0] = axis;
	vector<bool> refValueSets(1);
	refValueSets[0] = refValueSet;
	vector<double> refValues(1);
	refValues[0] = refValue;
	vector<PMS::DataColumn> datas(1);
	datas[0] = dataColumn;
	return getLabel( axes, refValueSets, refValues, datas, polnRatio );
}


String PlotMSLabelFormat::getLabel(vector<PMS::Axis> axes,
		vector<bool> refValueSets,
        vector<double> refValues,
        vector<PMS::DataColumn> datas,
        bool polnRatio) const {
	String axisLabel;
	int axesCount = axes.size();
	for ( int i = 0; i < axesCount; i++ ){
		vector<PMS::Axis> singleAxis(1, axes[i]);
		vector<bool> singleRefValueSets(1, refValueSets[i]);
		vector<double> singleRefValues(1, refValues[i]);
		vector<PMS::DataColumn> singleDatas(1, datas[i]);
		singleDatas[0] = datas[i];
		String yLabel = getLabel(format, axes[i], axes[i], singleAxis, 
                    refValueSets[i], refValues[i],
                    refValueSets[i], refValues[i],
                    singleRefValueSets, singleRefValues,
                    singleDatas[0], singleDatas, polnRatio);
		if ( i > 0 ){
			axisLabel.append( ", ");
		}
		axisLabel.append( yLabel );
	}
	return axisLabel;
}

String PlotMSLabelFormat::getLabel( PMS::Axis xAxis, vector<PMS::Axis> yAxis,
        bool xRefValueSet, double xRefValue, vector<bool> yRefValueSet,
        vector<double> yRefValue, PMS::DataColumn xData,
        vector<PMS::DataColumn> yDatas, bool polnRatio) const {
    String titleLabel = getLabel(format, xAxis, xAxis, yAxis, xRefValueSet, xRefValue,
        xRefValueSet, xRefValue, yRefValueSet, yRefValue, xData, yDatas, polnRatio);
    return titleLabel;
}

bool PlotMSLabelFormat::operator==(const PlotMSLabelFormat& other) const {
    return format == other.format; }

PlotMSLabelFormat& PlotMSLabelFormat::operator=(const PlotMSLabelFormat& copy){
    format = copy.format;
    return *this;
}

}
