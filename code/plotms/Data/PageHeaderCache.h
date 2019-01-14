#ifndef PAGEHEADERCACHE_H_
#define PAGEHEADERCACHE_H_

#include <map>
#include <vector>

#include <plotms/PlotMS/PlotMSPageHeaderParam.h>
#include <casacore/casa/BasicSL/String.h>

namespace casa {

using std::vector;
using std::map;

class HeaderItemData {
public:
	using HeaderItem = PageHeaderItemsDef::Item;
	HeaderItemData(HeaderItem item,const String & value)
		: item_(item), value_(value) {}
	HeaderItem item() const { return item_ ; }
	const String& value() const { return value_ ; }
	void setValue(const String& value) { value_ = value ; }
private:
	HeaderItem item_;
	String value_;
};

class PageHeaderCache {
private:
	static const String& emptyString() { static const String empty_ {""};  return empty_; }
public:
	using HeaderItem = PageHeaderItemsDef::Item;
	PageHeaderCache() {}
	void store(HeaderItemData data);
	void clear();
	const String& value(HeaderItem header_item) const {
		const auto it = itemIndex_.find(header_item);
		if ( it == itemIndex_.end() ) return emptyString();
		return values_[it->second].value();  }
	const std::vector<HeaderItemData>& values() const { return values_ ; }
	String toString() const ;

private:
	map<HeaderItem,uInt> itemIndex_;
	std::vector<HeaderItemData> values_;

};

ostream& operator<<(ostream &os, const PageHeaderCache&  cache);

}

#endif
