#include <plotms/Gui/PlotMSPageHeaderDataModel.qo.h>

#include <plotms/Plots/PlotMSPlot.h>
#include <plotms/Plots/PlotMSPlotParameterGroups.h>
#include <plotms/Data/MSCache.h>



namespace casa {

PlotMSPageHeaderDataModel::PlotMSPageHeaderDataModel(PlotMSApp* app , QObject *parent)
	: QAbstractTableModel(parent),
	  itsParent_(app)
{

}


int PlotMSPageHeaderDataModel::columnCount(const QModelIndex & /*parent*/) const {
	return physicalColumns;
}

inline bool PlotMSPageHeaderDataModel::isLogicalColumn(int column) const {
	return column == 0 || column == lastColumnIndex();
}

inline int PlotMSPageHeaderDataModel::logicalColumn(int column) const {
	if ( column == lastColumnIndex() ) return 1;
	return column;
}

int PlotMSPageHeaderDataModel::rowCount(const QModelIndex & /*parent*/) const {
	const auto & plots = itsParent_->getPlotManager().plots();
	int row_count = 0;
	for (size_t k=0; k<plots.size(); k++){
		auto page_header_group = plots[k]->parameters().typedGroup<PMS_PP_PageHeader>();
		if (page_header_group == NULL) continue;
		const auto & plot_header_items = page_header_group->pageHeaderItems().items();
		row_count += static_cast<int>( (plot_header_items.size()+1) / logicalColumns );
	}
	return row_count;
}

QVariant PlotMSPageHeaderDataModel::data(const QModelIndex &index, int role) const
{
	if (role == Qt::DisplayRole ) {
		if ( isLogicalColumn(index.column()) ) {
			const auto & plots = itsParent_->getPlotManager().plots();
			if (plots.empty()) return QVariant();
			// Have plots: find owner of item displayed at location pointed to by index
			// Logical Item Index
			int itemIndex = index.row()*logicalColumns + logicalColumn(index.column());
			int plotItemsStartIndex = 0;
			int plotIndex = 0;
			bool found = false;
			for ( const auto & plot : plots) {
				auto pageHeaderGroup = plot->parameters().typedGroup<PMS_PP_PageHeader>();
				if ( pageHeaderGroup == NULL ) {
					++plotIndex;
					continue;
				}
				const auto & plotHeaderItems = pageHeaderGroup->pageHeaderItems().items();
				if ( plotHeaderItems.empty() ) {
					++plotIndex;
					continue;
				}
			    int plotRows = (plotHeaderItems.size() + 1) / logicalColumns;
			    int nextPlotItemsStartIndex = plotItemsStartIndex + plotRows*logicalColumns;
				if ( plotItemsStartIndex <= itemIndex &&  itemIndex < nextPlotItemsStartIndex ) {
					found = true;
					break;
				} else {
					++plotIndex;
					plotItemsStartIndex = nextPlotItemsStartIndex;
				}
			}
			if ( ! found ) return QVariant();
			// Have owner plot: retrieve header item if any
			int headerItemIndex = itemIndex - plotItemsStartIndex;
			const auto & plot = plots[plotIndex];
			auto pageHeaderGroup = plot->parameters().typedGroup<PMS_PP_PageHeader>();
			const auto & plotHeaderItems = pageHeaderGroup->pageHeaderItems().items();
			// Last cell my be empty
			if (static_cast<size_t>(headerItemIndex) >= plotHeaderItems.size()) return QVariant();
			const auto & headerItem = plotHeaderItems[headerItemIndex];
			// Have header item: retrieve its value if any
			QString cellValue;
			// Header Item Label string
			using Items = PageHeaderItemsDef;
			const auto & headerItemLabel = Items::item2Info.at(headerItem).shortLabel();
			cellValue += QString::fromUtf8(headerItemLabel.c_str());
			cellValue += ":";
			auto & plotCache = plot->cache();
			if (! plotCache.cacheReady()) return cellValue;
			// Header Item Value
			auto plotMSCache = dynamic_cast<MSCache *>(&(plot->cache()));
			if (plotMSCache == nullptr) return cellValue;
			const auto & headerItemValue = plotMSCache->pageHeaderCache().value(headerItem);
			cellValue += QString(" ") + QString::fromUtf8(headerItemValue.c_str());
			return cellValue;
		}
		else { // Not a logical column: unused physical column
			return QVariant();
		}
	}
	else { // role != Qt::DisplayRole
		return QVariant();
	}
}


} // Namespace casa



