#include <limits>
#include <plotms/PlotMS/PlotMSPageHeaderParam.h>


namespace casa {

using Items = PageHeaderItemsDef;
using Item = Items::Item;

const map<HeaderItemInfo,Item>& Items::info2Item = Items::info2Item_();
const map<Item,HeaderItemInfo>& Items::item2Info = Items::item2Info_();
const map<String,Item>& Items::name2Item = Items::name2Item_();


const array<Item,Items::n_items>& Items::items(){
	static array<Item,Items::n_items> items_;
	static bool initialized = false;
	if ( ! initialized ) {
		for ( size_t k=0 ; k<Items::size(); k++){
			auto item = static_cast<Item>(k);
			items_[k] = item;
		}
		initialized = true;
	}
	return items_;
}


bool Items::isItemName(const String& s)
{
	return name2Item.find(s) != name2Item.end();
}

PageHeaderItems::PageHeaderItems(const String& items, const char sep)
{
	if ( items.empty() ) return;
	setItems(items,sep);
}

void PageHeaderItems::clear()
{
	items_.clear();
	item2index_.clear();
}

void PageHeaderItems::setItems(const String& items, const char sep)
{
	clear();
	String item_arr[kMaxItems];
	Int n_items = split(items,item_arr,kMaxItems,sep);

	for (Int k=0; k<n_items; k++) append(item_arr[k]);
}

const vector<Item>&  PageHeaderItems::items() const
{
	return items_;
}

bool PageHeaderItems::full() const
{
	return items_.size() >= kMaxItems;
}

bool PageHeaderItems::selected(Item item) const
{
	return item2index_.find(item) != item2index_.end();
}

bool PageHeaderItems::append(const String& item_string)
{
	if ( full() ) return false;
	if ( ! Items::isItemName(item_string) ) return false;
	Item item = Items::name2Item.find(item_string)->second;
	return append(item);
}

bool PageHeaderItems::append(Item item)
{
	if ( full() ) return false;
	if ( selected(item) ) return false;
	items_.push_back(item);
	item2index_[item] = items_.size() - 1;
	return true;
}

bool PageHeaderItems::remove(HeaderItem item)
{
	auto item2IndexIter = item2index_.find(item);
	if ( item2IndexIter == item2index_.end() ) return false;

	auto itemIndex = item2IndexIter->second;
	if ( itemIndex >= items_.size() ) return false;

	// Remove from vector
	items_.erase(items_.begin() + itemIndex);

	// Update map and remove item from map
	for ( auto & itemIndexPair : item2index_ ) {
		if ( itemIndexPair.second > itemIndex ) --(itemIndexPair.second);
	}
	item2index_.erase(item2IndexIter);
	return true;
}

pair<bool,size_t> PageHeaderItems::index(HeaderItem item) const
{
	pair<bool,size_t> result {false,numeric_limits<size_t>::max()};

	auto iter = item2index_.find(item);
	if ( iter != item2index_.end() ){
		result.first = true;
		result.second = iter->second;
	}
	return result;
}

}




