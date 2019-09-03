//# QPScatterPlot.cc: Qwt implementation of generic ScatterPlot class.
//# Copyright (C) 2008
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
#ifdef AIPS_HAS_QWT

#include <casaqt/QwtPlotter/QPScatterPlot.h>
#include <casaqt/QwtPlotter/QPFactory.h>

#if QWT_VERSION < 0x060000
#include <qwt_legend_item.h>
#endif

#include <QPainter>
#include <QDebug>
#include <iomanip>
using namespace casacore;
namespace casa {

///////////////////////////////
// QPSCATTERPLOT DEFINITIONS //
///////////////////////////////

// Static //

const String QPScatterPlot::CLASS_NAME = "QPScatterPlot";


// Constructors/Destructors //

QPScatterPlot::QPScatterPlot(PlotPointDataPtr data, const String& title):
        m_data(data), m_symbol(QPFactory::defaultPlotSymbol()),
        m_line(QPFactory::defaultPlotLine()),
        m_step(false),
        m_maskedSymbol(QPFactory::defaultPlotMaskedSymbol()),
        m_maskedLine(QPFactory::defaultPlotLine()),
        m_errorLine(QPFactory::defaultPlotLine()),
        m_errorCap(QPFactory::DEFAULT_ERROR_CAP),
        m_maskedStep(false) {
    QPPlotItem::setTitle(title);
    
    if(!data.null()) {
        PlotMaskedPointData* p = dynamic_cast<PlotMaskedPointData*>(&*data);
        if(p != NULL) m_maskedData = PlotMaskedPointDataPtr(p, false);
        PlotErrorData* e = dynamic_cast<PlotErrorData*>(&*data);
        if(e != NULL) m_errorData = PlotErrorDataPtr(e, false);
        PlotBinnedData* b = dynamic_cast<PlotBinnedData*>(&*data);
        if(b != NULL) m_coloredData = PlotBinnedDataPtr(b, false);
    }
    
    setItemAttribute(QwtPlotItem::AutoScale);
    setItemAttribute(QwtPlotItem::Legend, true);
}

QPScatterPlot::QPScatterPlot(const ScatterPlot& copy) :
        m_data(copy.pointData()), m_symbol(QPFactory::defaultPlotSymbol()),
        m_line(QPFactory::defaultPlotLine()),
        m_step(false),
        m_maskedSymbol(QPFactory::defaultPlotMaskedSymbol()),
        m_maskedLine(QPFactory::defaultPlotLine()),
        m_errorLine(QPFactory::defaultPlotLine()),
        m_errorCap(QPFactory::DEFAULT_ERROR_CAP), 
        m_maskedStep(false) {
    QPPlotItem::setTitle(copy.title());
    
    if(!m_data.null()) {
        PlotMaskedPointData* p = dynamic_cast<PlotMaskedPointData*>(&*m_data);
        if(p != NULL) m_maskedData = PlotMaskedPointDataPtr(p, false);
        PlotErrorData* e = dynamic_cast<PlotErrorData*>(&*m_data);
        if(e != NULL) m_errorData = PlotErrorDataPtr(e, false);
        PlotBinnedData* b = dynamic_cast<PlotBinnedData*>(&*m_data);
        if(b != NULL) m_coloredData = PlotBinnedDataPtr(b, false);
    }
    
    setLine(copy.line());
    setSymbol(copy.symbol());
    
    // check for masked plot
    const MaskedScatterPlot* m = dynamic_cast<const MaskedScatterPlot*>(&copy);
    if(m != NULL) {
        setMaskedLine(m->maskedLine());
        setMaskedSymbol(m->maskedSymbol());
    }
    
    // check for error plot
    const ErrorPlot* e = dynamic_cast<const ErrorPlot*>(&copy);
    if(e != NULL) {
        setErrorLine(e->errorLine());
        setErrorCapSize(e->errorCapSize());
    }
    
    setItemAttribute(QwtPlotItem::AutoScale);
    setItemAttribute(QwtPlotItem::Legend, true);
}

QPScatterPlot::~QPScatterPlot() {
    foreach(QPColor* c, m_colors) delete c;
    m_colors.clear();
    logDestruction();
}


// Public Methods //

bool QPScatterPlot::isValid() const {
    return !m_data.null() && m_data->isValid(); }


bool QPScatterPlot::shouldDraw() const {
    bool shouldDraw = isValid() && m_data->size() > 0;
    return shouldDraw;
}

QwtDoubleRect QPScatterPlot::boundingRect() const {
    bool ret;
    double xMin, xMax, yMin, yMax;
    
    // Determine which bounding rect to return.
    if(!m_maskedData.null()) {
        bool showMasked = m_maskedLine.style() != PlotLine::NOLINE ||
                          m_maskedSymbol.symbol() != PlotSymbol::NOSYMBOL,
             showNormal = m_line.style() != PlotLine::NOLINE ||
                          m_symbol.symbol() != PlotSymbol::NOSYMBOL;
        PlotMaskedPointDataPtr data = m_maskedData;
        
        if(showMasked && !showNormal)
            ret = data->maskedMinsMaxes(xMin, xMax, yMin, yMax);
        else if(showNormal && !showMasked)
            ret = data->unmaskedMinsMaxes(xMin, xMax, yMin, yMax);
        else
            ret = data->minsMaxes(xMin, xMax, yMin, yMax);
        
    } else if (!m_data.null()) {
        ret = const_cast<PlotPointDataPtr&>(m_data)->minsMaxes(
                xMin, xMax, yMin, yMax);
    } else ret=0;
    
    // have to switch y min and max for some reason..
    if(!ret) return QwtDoubleRect();
    else return QwtDoubleRect(QPointF(xMin, yMin), QPointF(xMax, yMax));
}

#if QWT_VERSION >= 0x060000
QwtGraphic QPScatterPlot::legendIcon(int /*index*/, const QSizeF& size) const {
    const QPSymbol* plotsymbol(nullptr);
    if (symbolsShown()) {
        plotsymbol = &m_symbol;
    } else if (maskedSymbolsShown()) {
        plotsymbol = &m_maskedSymbol;
    }
    QwtSymbol* iconsymbol = new QwtSymbol(QwtSymbol::HLine);
    QwtGraphic icon = iconsymbol->graphic();
    if (plotsymbol != nullptr) {
        iconsymbol->setPen(plotsymbol->drawPen().color());
        icon = iconsymbol->graphic();
        icon.setDefaultSize(QSizeF(size.width()*2.0, size.height()*2.0));
        QPainter painter(&icon);
        const QRectF rect (0.0, 0.0, size.width()*2.0, size.height()*2.0);
        iconsymbol->drawSymbol(&painter, rect);
    }
    return icon;
}
#else
QWidget* QPScatterPlot::legendItem() const {
    QPen legendPen = m_line.asQPen();
    if ( !linesShown() ){
        QColor penColor = m_symbol.drawPen().color();
        legendPen = QPen(penColor );
    }
    QwtLegendItem* i= new QwtLegendItem(m_symbol, legendPen, qwtTitle());
    i->setIdentifierMode(QwtLegendItem::ShowLine | QwtLegendItem::ShowSymbol |
                         QwtLegendItem::ShowText);
    return i;
}
#endif

bool QPScatterPlot::linesShown() const {
    return m_line.style() != PlotLine::NOLINE; }
void QPScatterPlot::setLinesShown(bool l) {
    if(l != linesShown()) {
        m_line.setStyle(l ? PlotLine::SOLID : PlotLine::NOLINE);
        itemChanged();
    }
}

PlotLinePtr QPScatterPlot::line() const { return new QPLine(m_line); }
void QPScatterPlot::setLine(const PlotLine& line) {
    if(m_line != line) {
        m_line = line;
        itemChanged();
    }
}


PlotPointDataPtr QPScatterPlot::pointData() const {
    return m_data;
}

bool QPScatterPlot::symbolsShown() const {
    return m_symbol.symbol() != PlotSymbol::NOSYMBOL; }
void QPScatterPlot::setSymbolsShown(bool show) {
    if(show != symbolsShown()) {
        m_symbol.setSymbol(show ? PlotSymbol::CIRCLE : PlotSymbol::NOSYMBOL);
        itemChanged();
    }
}

PlotSymbolPtr QPScatterPlot::symbol() const {
#if QWT_VERSION >= 0x060000
    psize_t casasize = m_symbol.size();
    QSize qsize = QSize(casasize.first, casasize.second);
    QPSymbol* symbol = new QPSymbol(m_symbol.style(), 
    m_symbol.drawBrush(), m_symbol.drawPen(), qsize);
    return symbol;
#else
    return new QPSymbol(m_symbol);
#endif
}

void QPScatterPlot::setSymbol(const PlotSymbol& sym) {
    if(sym != m_symbol) {
        m_symbol = sym;
        updateBrushes();
        itemChanged();
    }
}


PlotMaskedPointDataPtr QPScatterPlot::maskedData() const{ return m_maskedData;}

void QPScatterPlot::clearData() {
    // null pointer when data is deleted in another thread
    m_maskedData = PlotMaskedPointDataPtr();
    m_data = PlotPointDataPtr();
}

bool QPScatterPlot::maskedLinesShown() const {
    return m_maskedLine.style() != PlotLine::NOLINE; }
void QPScatterPlot::setMaskedLinesShown(bool s) {
    if(s != maskedLinesShown()) {
        m_maskedLine.setStyle(s ? PlotLine::SOLID : PlotLine::NOLINE);
        itemChanged();
    }
}

PlotLinePtr QPScatterPlot::maskedLine() const {
    return new QPLine(m_maskedLine); }
void QPScatterPlot::setMaskedLine(const PlotLine& line) {
    if(m_maskedLine != line) {
        m_maskedLine = line;
        itemChanged();
    }
}

bool QPScatterPlot::maskedSymbolsShown() const {
    return m_maskedSymbol.symbol() != PlotSymbol::NOSYMBOL; }
void QPScatterPlot::setMaskedSymbolsShown(bool s) {
    if(s != maskedSymbolsShown()) {
        m_maskedSymbol.setSymbol(s ? PlotSymbol::CIRCLE :
                                     PlotSymbol::NOSYMBOL);
        itemChanged();
    }
}

PlotSymbolPtr QPScatterPlot::maskedSymbol() const {
#if QWT_VERSION >= 0x060000
    psize_t casasize = m_maskedSymbol.size();
    QSize qsize = QSize(casasize.first, casasize.second);
    QPSymbol* symbol = new QPSymbol(m_maskedSymbol.style(), 
    m_maskedSymbol.drawBrush(), m_maskedSymbol.drawPen(), qsize);
    return symbol;
#else
    return new QPSymbol(m_maskedSymbol);
#endif
}

void QPScatterPlot::setMaskedSymbol(const PlotSymbol& symbol) {
    if(symbol != m_maskedSymbol) {
        m_maskedSymbol = symbol;
        updateBrushes();
        itemChanged();
    }
}


PlotErrorDataPtr QPScatterPlot::errorData() const { return m_errorData; }

bool QPScatterPlot::errorLineShown() const {
    return m_errorLine.style() != PlotLine::NOLINE; }
void QPScatterPlot::setErrorLineShown(bool s) {
    if(s != errorLineShown()) {
        m_errorLine.setStyle(s ? PlotLine::SOLID : PlotLine::NOLINE);
        itemChanged();
    }
}

PlotLinePtr QPScatterPlot::errorLine() const{ return new QPLine(m_errorLine); }
void QPScatterPlot::setErrorLine(const PlotLine& line) {
    if(m_errorLine != line) {
        m_errorLine = line;
        itemChanged();
    }
}

unsigned int QPScatterPlot::errorCapSize() const { return m_errorCap; }

void QPScatterPlot::setErrorCapSize(unsigned int capSize) {
    if(m_errorCap != capSize) {
        m_errorCap = capSize;
        itemChanged();
    }
}


PlotBinnedDataPtr QPScatterPlot::binnedColorData() const {
    return m_coloredData; }

PlotColorPtr QPScatterPlot::colorForBin(unsigned int bin) const {
    if((int)bin >= m_colors.size() || m_colors[bin] == NULL)
        return PlotColorPtr();
    else return new QPColor(*m_colors[bin]);
}

void QPScatterPlot::setColorForBin(unsigned int bin, const PlotColorPtr color){
    // Add any intermediate entries, if needed.
    for(unsigned int i = m_colors.size(); i < bin; i++)
        m_colors.insert(bin, NULL);
    
    if(color.null()) m_colors.insert(bin, NULL);
    else m_colors.insert(bin, new QPColor(*color));
    
    // Update brushes
    updateBrushes();
}


// Protected Methods //

#if QWT_VERSION >= 0x060000
void QPScatterPlot::draw_(QPainter* p, const QwtScaleMap& xMap,
        const QwtScaleMap& yMap, const QRectF& brect,
        unsigned int drawIndex, unsigned int drawCount) const {
#else
void QPScatterPlot::draw_(QPainter* p, const QwtScaleMap& xMap,
        const QwtScaleMap& yMap, const QRect& brect,
        unsigned int drawIndex, unsigned int drawCount) const {
#endif
    if (!shouldDraw()) return;
    unsigned int n = m_data->size();
    if(n == 0 || drawIndex >= n) {
        return;
    }
        
    if(drawIndex + drawCount > n){
        drawCount = n - drawIndex;
    }
    n = drawIndex + drawCount;

    p->save();
    
    // Draw error lines if needed.
    bool drawErrorLine = !m_errorData.null() &&
                         m_errorLine.style() != PlotLine::NOLINE;
    if(drawErrorLine) {
        p->setPen(m_errorLine.asQPen());
        unsigned int cap = m_errorCap / 2;
        double tempx, tempy, txleft, txright, tybottom, tytop;
        int temp, min, max;
        
        bool drawNormally = m_maskedData.null();
        if(!drawNormally) {
            if(m_maskedSymbol.symbol() != PlotSymbol::NOSYMBOL ||
               m_maskedLine.style() != PlotLine::NOLINE) {
                // masked points are shown, so draw all error bars
                drawNormally = true;
                
            } else {
                // only draw unmasked error lines
                for(unsigned int i = drawIndex; i < n; i++) {
                    if(!m_maskedData->maskedAt(i)) {
                        m_errorData->xyAndErrorsAt(i, tempx, tempy, txleft,
                                txright, tybottom, tytop);
                        
                        min = xMap.transform(tempx - txleft);
                        max = xMap.transform(tempx + txright);
                        temp = yMap.transform(tempy);
                        
                        if(brect.contains(min,temp)||brect.contains(max,temp)){
                            p->drawLine(min, temp, max, temp);                        
                            if(cap > 0) {
                                p->drawLine(min, temp - cap, min, temp + cap);
                                p->drawLine(max, temp - cap, max, temp + cap);
                            }
                        }
                        
                        min = yMap.transform(tempy - tybottom);
                        max = yMap.transform(tempy + tytop);
                        temp = xMap.transform(tempx);
                        
                        if(brect.contains(temp,min)||brect.contains(temp,max)){
                            p->drawLine(temp, min, temp, max);                        
                            if(cap > 0) {
                                p->drawLine(temp - cap, min, temp + cap, min);
                                p->drawLine(temp - cap, max, temp + cap, max);
                            }
                        }
                    }
                }
            }            
        }
            
        // draw all error lines
        if(drawNormally) {
            for(unsigned int i = drawIndex; i < n; i++) {
                m_errorData->xyAndErrorsAt(i, tempx, tempy, txleft,
                        txright, tybottom, tytop);
                
                min = xMap.transform(tempx - txleft);
                max = xMap.transform(tempx + txright);
                temp = yMap.transform(tempy);
                
                if(brect.contains(min, temp) || brect.contains(max, temp)) {
                    p->drawLine(min, temp, max, temp);                
                    if(cap > 0) {
                        p->drawLine(min, temp - cap, min, temp + cap);
                        p->drawLine(max, temp - cap, max, temp + cap);
                    }
                }
                
                min = yMap.transform(tempy - tybottom);
                max = yMap.transform(tempy + tytop);
                temp = xMap.transform(tempx);
                
                if(brect.contains(temp, min) || brect.contains(temp, max)) {
                    p->drawLine(temp, min, temp, max);                
                    if(cap > 0) {
                        p->drawLine(temp - cap, min, temp + cap, min);
                        p->drawLine(temp - cap, max, temp + cap, max);
                    }
                }
            }
        }
    }

    unsigned int numBins(1);
    if (!m_coloredData.null()) numBins = m_coloredData->numBins();
    // Draw normal/masked lines
    bool drawLine = m_line.style() != PlotLine::NOLINE,
         drawMaskedLine = !m_maskedData.null() &&
                          m_maskedLine.style() != PlotLine::NOLINE,
         diffColorLine = !m_coloredData.null() && m_coloredData->isBinned();
    if(drawLine || drawMaskedLine) {
        // connect last point to this point
        double lastx, lasty, thisx, thisy;
        int lastix, lastiy, thisix, thisiy;
        if(!m_maskedData.null()) {
            bool lastmask, thismask, sameMask;
            bool sameBin(true), sameLine(true);
            
            // set the painter's pen only once if possible
            // if unmasked==masked color or only one unmasked/masked line used
            bool samePen(m_line.asQPen() == m_maskedLine.asQPen());
            if (!drawMaskedLine || samePen) p->setPen(m_line.asQPen());
            if (!drawLine) p->setPen(m_maskedLine.asQPen());
           
            // get first point
            m_maskedData->xyAndMaskAt(drawIndex, lastx, lasty, lastmask);
            lastix = xMap.transform(lastx);
            lastiy = yMap.transform(lasty);

            unsigned int thisColorBin;  // for colorized data
            unsigned int lastConnectBin = 0, thisConnectBin = 0; // whether to connect
            if (diffColorLine) {
                lastConnectBin = m_coloredData->connectBinAt(drawIndex);
            }
            for(unsigned int i = drawIndex + 1; i < n; i++) {
                // get this point
                m_maskedData->xyAndMaskAt(i, thisx, thisy, thismask);
                thisix = xMap.transform(thisx);
                thisiy = yMap.transform(thisy);
                // don't connect nan and inf points, or points at same x 
                if (!casacore::isNaN(thisx) && !casacore::isNaN(thisy) &&
                    !casacore::isInf(thisx) && !casacore::isInf(thisy) &&
                    (thisx != lastx)) {

                    bool reversed(m_maskedData->reverseConnect(i));
                    if (diffColorLine) {  // set pen for colorized plot
                        QPen linePen;
                        if (!drawMaskedLine || samePen)
                           linePen = m_line.asQPen();
                        if (!drawLine)
                           linePen = m_maskedLine.asQPen();
                        thisColorBin = m_coloredData->binAt(i);
                        thisConnectBin = m_coloredData->connectBinAt(i);
                        sameBin = (thisConnectBin==lastConnectBin);
                        unsigned int colorBin = thisColorBin % numBins;
                        QBrush coloredBrush = m_coloredBrushes[colorBin];
                        linePen.setBrush(coloredBrush);
                        QColor brushColor = coloredBrush.color();
                        linePen.setColor(brushColor);
                        p->setPen(linePen);
                        // compare connect bins then save for next point
                        if (reversed)
                            sameLine = sameBin && (thisx < lastx); // for colorized data
                        else
                            sameLine = sameBin && (thisx > lastx); // for colorized data
                        lastConnectBin = thisConnectBin;
                    }

                    sameMask = (thismask==lastmask); // only connect same-masked points
                    bool secondpoint(i==drawIndex+1), morepoints(i<n-1);

                    double nextx, nexty;
                    bool nextmask;
                    if (morepoints) {
                        m_maskedData->xyAndMaskAt(i+1, nextx, nexty, nextmask);
                    }

                    if (drawLine && !thismask) { // connect unflagged pts
                        if(drawMaskedLine && !samePen && !diffColorLine)
                            p->setPen(m_line.asQPen()); // set pen for unflagged symbols
                        if(brect.contains(lastix, lastiy) || brect.contains(thisix, thisiy)) {
                            if (m_step) {  // need midpoint for step
                                if (sameMask && sameLine) { // connect to last point
                                    int ixdiff(((lastix + thisix)/2) - lastix); // midpoint for step
                                    // 3 lines: last pt to mid; vertical line; mid to this point
                                    p->drawLine(lastix, lastiy, thisix-ixdiff, lastiy);
                                    p->drawLine(thisix-ixdiff, lastiy, thisix-ixdiff, thisiy);
                                    p->drawLine(thisix-ixdiff, thisiy, thisix, thisiy);

                                    if (secondpoint) { // draw leading line for first point in plot
                                        p->drawLine(lastix-ixdiff, lastiy, lastix, lastiy);
                                    }

                                    // normally, next point draws trailing line for last point unless:
                                    // 1. this is last point in plot or in this line (bin changes)
                                    // 2. same bin but mask changes for next point
                                    if (!morepoints || (morepoints && (m_coloredData->connectBinAt(i+1)!=thisConnectBin))) {
                                        p->drawLine(thisix, thisiy, thisix+ixdiff, thisiy); // last point
                                    }
                                    if ((morepoints && (m_coloredData->connectBinAt(i+1) == thisConnectBin)) &&
                                        (nextmask != thismask)) {
                                        p->drawLine(thisix, thisiy, thisix+ixdiff, thisiy); // mask changes
                                    }
                                } else if ((i < n-1) && (m_coloredData->connectBinAt(i+1)==thisConnectBin))  {
                                    // first point in new line, use next point to draw leading line
                                    if (nextmask==thismask) {
                                        int nextix(xMap.transform(nextx));
                                        int ixdiff(((thisix + nextix)/2) - thisix);
                                        // don't want ix<0
                                        if (ixdiff > thisix) ixdiff=thisix;
                                        p->drawLine(thisix-ixdiff, thisiy, thisix, thisiy);
                                    }
                                }
                            } else if (sameMask && sameLine) { // connect last point to this point
                                p->drawLine(lastix, lastiy, thisix, thisiy);
                            }
                        }
                    } else if(drawMaskedLine && thismask) { // connect flagged points
                        if(drawLine && !samePen && !diffColorLine)
                            p->setPen(m_maskedLine.asQPen());  // set pen for flagged symbols
                        if(brect.contains(lastix, lastiy) || brect.contains(thisix, thisiy)) {
                            if (m_maskedStep) {
                                if (sameMask && sameLine) {
                                    int ixdiff(((lastix + thisix)/2) - lastix); // midpoint for step
                                    p->drawLine(lastix, lastiy, thisix-ixdiff, lastiy);
                                    p->drawLine(thisix-ixdiff, lastiy, thisix-ixdiff, thisiy);
                                    p->drawLine(thisix-ixdiff, thisiy, thisix, thisiy);
                                    if (secondpoint) {
                                        p->drawLine(lastix-ixdiff, lastiy, lastix, lastiy);
                                    }
                                    if (!morepoints || (morepoints && (m_coloredData->connectBinAt(i+1)!=thisConnectBin))) {
                                        p->drawLine(thisix, thisiy, thisix+ixdiff, thisiy);
                                    }
                                    if ((morepoints && (m_coloredData->connectBinAt(i+1) == thisConnectBin)) &&
                                        (nextmask != thismask)) {
                                        p->drawLine(thisix, thisiy, thisix+ixdiff, thisiy);
                                    }
                                } else if ((i < n-1) && (m_coloredData->connectBinAt(i+1)==thisConnectBin))  {
                                    double nextx, nexty;
                                    bool nextmask;
                                    m_maskedData->xyAndMaskAt(i+1, nextx, nexty, nextmask);
                                    if (nextmask==thismask) {
                                        int nextix(xMap.transform(nextx));
                                        int ixdiff(((thisix + nextix)/2) - thisix);
                                        if (ixdiff > thisix) ixdiff=thisix;
                                        p->drawLine(thisix-ixdiff, thisiy, thisix, thisiy);
                                    }
                                }
                            } else if (sameMask && sameLine) { // connect last point to this point
                                p->drawLine(lastix, lastiy, thisix, thisiy);
                            }
                        }
                    }
                }
                lastx = thisx; lasty = thisy;
                lastix = thisix; lastiy = thisiy;
                lastmask = thismask;
            }
        } else {
            p->setPen(m_line.asQPen());
            m_data->xAndYAt(drawIndex, lastx, lasty);
            lastix = xMap.transform(lastx);
            lastiy = yMap.transform(lasty);
            for(unsigned int i = drawIndex + 1; i < n; i++) {
                m_data->xAndYAt(i, thisx, thisy);
                thisix = xMap.transform(thisx); thisiy = yMap.transform(thisy);
                if(brect.contains(lastix, lastiy) || brect.contains(thisix, thisiy))
                    p->drawLine(lastix, lastiy, thisix, thisiy);
                lastix = thisix; lastiy = thisiy;
            }
        }
    }
 
    // Draw normal/masked symbols
    bool drawSymbol = m_symbol.symbol() != PlotSymbol::NOSYMBOL,
         drawMaskedSymbol = !m_maskedData.null() &&
                            m_maskedSymbol.symbol() != PlotSymbol::NOSYMBOL,
         diffColor = !m_coloredData.null() && m_coloredData->isBinned();    
    if(drawSymbol || drawMaskedSymbol) {
        double tempx, tempy;
        if(!m_maskedData.null()) {
            unsigned int ptsToDraw = ( m_maskedData->plotConjugates() ? 2 : 1 );
            bool mask;
            
            const QPen& pen = m_symbol.drawPen(),
                      & mpen = m_maskedSymbol.drawPen();
            const QBrush& brush = m_symbol.drawBrush(),
                        & mbrush = m_maskedSymbol.drawBrush();
            bool samePen = pen == mpen, sameBrush = brush == mbrush;
            
            // set the painter's pen/brush only once if possible
            if(!drawMaskedSymbol || samePen)  // use unmasked pen
                p->setPen(pen);
            else if(!drawSymbol)  // only masked symbols; use masked pen
                p->setPen(mpen);
            if(!drawMaskedSymbol || sameBrush)  // use unmasked brush
                p->setBrush(brush);
            else if(!drawSymbol)  // only masked symbols; use masked brush
                p->setBrush(mbrush);
            
            QSize size = ((QwtSymbol&)m_symbol).size();
            QRect rect(0, 0, size.width(), size.height());
            size = ((QwtSymbol&)m_maskedSymbol).size();
            QRect mRect(0, 0, size.width(), size.height());

#if QWT_VERSION >= 0x060000
            bool symIsPixel(m_symbol.symbol()==PlotSymbol::PIXEL),
            msymIsPixel(m_maskedSymbol.symbol()==PlotSymbol::PIXEL);
            std::vector<QPointF> upoints, mpoints;
            QPoint qpt;
            for(unsigned int i = drawIndex; i < n; i++) {
                m_maskedData->xyAndMaskAt(i, tempx, tempy, mask);
                // don't plot nan and inf !
                if (!casacore::isNaN(tempx) && !casacore::isNaN(tempy) &&
                    !casacore::isInf(tempx) && !casacore::isInf(tempy)) {
                  if(drawSymbol && !mask) { // unflagged points
                    for (unsigned int pt=0; pt<ptsToDraw; ++pt) {
                        // set point position
                        if (pt==0) {
                            qpt = QPoint(xMap.transform(tempx),yMap.transform(tempy));
                        } else { // conjugate for uv plot
                            qpt = QPoint(xMap.transform(-tempx),yMap.transform(-tempy));
                        }
                        QPointF qptf = QPointF(qpt);

                        // check if point is in drawing area
                        rect.moveCenter(qpt);
                        if (!brect.intersects(rect)) {
                            continue;
                        }

                        // if drawing mixed symbols, set back to unmasked pen/brush
                        if(drawMaskedSymbol) {
                            if(!samePen) p->setPen(pen);
                            if(!sameBrush) p->setBrush(brush);
                        }

                        if (diffColor) {
                            // set pen/brush for color bin and draw point/symbol
                            unsigned int colorBin = (m_coloredData->binAt(i)) % numBins;
                            QBrush coloredBrush = m_coloredBrushes[colorBin];
                            QColor brushColor = coloredBrush.color();
                            p->setBrush(coloredBrush);
                            p->setPen(brushColor);
                            if (symIsPixel) {
                                p->drawPoint(qptf);
                            } else {
                                QPSymbol* coloredSym = coloredSymbol(brushColor);
                                coloredSym->drawSymbol(p, qptf);
                                delete coloredSym;
                            }
                        } else {
                            // draw point/symbol in batch mode
                            upoints.push_back(qptf);
                            if (upoints.size() == 15000) {
                                if (symIsPixel) {
                                    p->drawPoints(&upoints[0], upoints.size());
                                } else {
                                    m_symbol.drawSymbols(p, &upoints[0], upoints.size());
                                }
                                upoints.clear();
                            }
                        }
                    }
                  } else if (drawMaskedSymbol && mask) { // flagged points
                    for (unsigned int pt=0; pt<ptsToDraw; ++pt) {
                        // set point position
                        if (pt==0) {
                            qpt = QPoint(xMap.transform(tempx), yMap.transform(tempy));
                        } else {
                            qpt = QPoint(xMap.transform(-tempx), yMap.transform(-tempy));
                        }
                        QPointF qptf = QPointF(qpt);

                        // check if point is in drawing area
                        mRect.moveCenter(qpt);
                        if (!brect.intersects(mRect)) {
                            continue;
                        }

                        // if drawing mixed symbols, set back to masked pen/brush
                        if (drawSymbol) {
                            if(!samePen) p->setPen(mpen);
                            if(!sameBrush) p->setBrush(mbrush);
                        }

                        if (diffColor) {
                            // set pen/brush for color bin and draw point/symbol
                            unsigned int colorBin = (m_coloredData->binAt(i)) % numBins;
                            QBrush coloredBrush = m_coloredBrushes[colorBin];
                            QColor brushColor = coloredBrush.color();
                            p->setBrush(coloredBrush);
                            p->setPen(brushColor);
                            if (symIsPixel) {
                                p->drawPoint(qptf);
                            } else {
                                QPSymbol* coloredSym = coloredSymbol(brushColor);
                                coloredSym->drawSymbol(p, qptf);
                                delete coloredSym;
                            }
                        } else {
                            mpoints.push_back(qptf);
                            if (mpoints.size() == 15000) {
                                if (msymIsPixel) {
                                    p->drawPoints(&mpoints[0], mpoints.size());
								} else {
                                    m_maskedSymbol.drawSymbols(p, &mpoints[0], mpoints.size());
								}
                                mpoints.clear();
                            }
                        }
                    }
                  }
                }
            }
            // draw the rest
            if (!diffColor && drawSymbol) {
                if (symIsPixel) p->drawPoints(&upoints[0], upoints.size());
                else m_symbol.drawSymbols(p, &upoints[0], upoints.size());
            }
            if (!diffColor && drawMaskedSymbol) {
                if (msymIsPixel) p->drawPoints(&mpoints[0], mpoints.size());
                else m_maskedSymbol.drawSymbols(p, &mpoints[0], mpoints.size());
            }
#else
            QPoint qpt;
            for(unsigned int i = drawIndex; i < n; i++) {
                m_maskedData->xyAndMaskAt(i, tempx, tempy, mask);
                // don't plot nan and inf !
                if (!casacore::isNaN(tempx) && !casacore::isNaN(tempy) &&
                    !casacore::isInf(tempx) && !casacore::isInf(tempy)) {
                  if(drawSymbol && !mask) {
                    for (unsigned int pt=0; pt<ptsToDraw; ++pt) {
                        if (pt==0) {
                            qpt = QPoint(xMap.transform(tempx),yMap.transform(tempy));
                        } else { // conjugate for uv plot
                            qpt = QPoint(xMap.transform(-tempx),yMap.transform(-tempy));
                        }
                        rect.moveCenter(qpt);
                        if(!brect.intersects(rect)) continue;
                        if(drawMaskedSymbol) {
                            if(!samePen) p->setPen(pen);
                            if(!sameBrush) p->setBrush(brush);
                        }
                        if(diffColor) {
                            unsigned int colorBin = (m_coloredData->binAt(i)) % numBins;
                            QBrush coloredBrush = m_coloredBrushes[colorBin];
                            QColor brushColor = coloredBrush.color();
                            p->setBrush(coloredBrush);
                            p->setPen(brushColor);
                        }
                        m_symbol.draw(p, rect);
                    }
                  } else if(drawMaskedSymbol && mask) {
                    for (unsigned int pt=0; pt<ptsToDraw; ++pt) {
                        if (pt==0) {
                            qpt = QPoint(xMap.transform(tempx), yMap.transform(tempy));
                        } else {
                            qpt = QPoint(xMap.transform(-tempx), yMap.transform(-tempy));
                        }
                        mRect.moveCenter(qpt);
                        if(!brect.intersects(mRect)) continue;
                        if(drawMaskedSymbol) {
                            if(!samePen) p->setPen(mpen);
                            if(!sameBrush) p->setBrush(mbrush);
                        }
                        if(diffColor) {
                            unsigned int colorBin = (m_coloredData->binAt(i)) % numBins;
                            QBrush coloredBrush = m_coloredBrushes[colorBin];
                            QColor brushColor = coloredBrush.color();
                            p->setBrush(coloredBrush);
                            p->setPen(brushColor);
                        }
                        m_maskedSymbol.draw(p, mRect);
                    }
                  }
                }
            }
#endif
        } else {
            // draw all symbols normally
            const QBrush& brush = m_symbol.drawBrush();
            p->setPen(m_symbol.drawPen());
            p->setBrush(brush);
            
            QSize size = ((QwtSymbol&)m_symbol).size();
            QRect rect(0, 0, size.width(), size.height());
                        
            for(unsigned int i = drawIndex; i < n; i++) {
                m_data->xAndYAt(i, tempx, tempy);
                rect.moveCenter(QPoint(xMap.transform(tempx),
                                       yMap.transform(tempy)));
                if(!brect.intersects(rect)) continue;
                if(diffColor) {
                    unsigned int colorBin = m_coloredData->binAt(i) % numBins;
                    p->setBrush(m_coloredBrushes[colorBin]);
                }
                m_symbol.draw(p, rect);
            }
        }
    }

    p->restore();
}


// Private Methods //

#if QWT_VERSION >= 0x060000
QPSymbol* QPScatterPlot::coloredSymbol(const QColor& color) const {
    // Qwt 6 uses pen color from symbol not painter;
    // Need a non-const symbol
    QPSymbol* coloredSym = new QPSymbol();
    coloredSym->setStyle(m_symbol.style());
    coloredSym->setSize(m_symbol.size().first, m_symbol.size().second);
    coloredSym->setBrush(color);
    coloredSym->setPen(color);
    return coloredSym;
}
#endif

void QPScatterPlot::updateBrushes() {
    m_coloredBrushes.clear();
    if(m_coloredData.null()) return;

    const QBrush& brush = m_symbol.drawBrush();

    //bool allSame = true;
    for(unsigned int i = 0; i < m_coloredData->numBins(); i++) {
        m_coloredBrushes << brush;
        if((int)i < m_colors.size() && m_colors[i] != NULL) {
            m_coloredBrushes[i].setColor(m_colors[i]->asQColor());
            //allSame &= brushes[i].color() == brush.color();
        }
    }
    //if(allSame) diffColor = false;
}

}

#endif
