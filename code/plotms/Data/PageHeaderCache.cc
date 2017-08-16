#include <plotms/Data/PageHeaderCache.h>

#include <sstream>

namespace casa {

using ItemsDef = PageHeaderItemsDef;

void PageHeaderCache::store(HeaderItemData data) {
	auto item_pos = itemIndex_.find(data.item());
	if (item_pos == itemIndex_.end()) {
		// New item: append
		values_.push_back(data);
		itemIndex_[data.item()] = values_.size() - 1;
	}
	else {
		// Existing item: update value
		auto itemIndex = item_pos->second;
		values_[itemIndex] = data;
	}
}

void PageHeaderCache::clear() {
	itemIndex_.clear();
	values_.clear();
}

String PageHeaderCache::toString() const {
	const auto & v = values_;
	if ( v.empty() ) return "";
	std::stringstream ss;
	ss << v[0].value();
	for (size_t k = 1 ; k < v.size(); k++) {
		ss << std::endl;
		ss << v[k].value();
	}
	return ss.str();

}

ostream& operator<<(ostream &os, const PageHeaderCache&  cache) {
	const auto & v = cache.values();
	if ( v.empty() ) return os;
	os << "label_value" << ": " << v[0].value();
	for (size_t k = 1 ; k < v.size(); k++) {
		os << endl;
		os << "label_value" << ": " << v[k].value();
	}
	return os;
}

}


