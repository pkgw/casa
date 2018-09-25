#ifndef PLOTMSPAGEHEADERPARAM_H_
#define PLOTMSPAGEHEADERPARAM_H_

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <array>
#include <algorithm>

#include <casa/namespace.h>
#include <casacore/casa/BasicSL/String.h>


namespace casa {

using namespace std;


class HeaderItemInfo {
public:
	HeaderItemInfo(String name=String(),String shortLabel=String(),String longLabel=String())
		: name_(name),
          shortLabel_(shortLabel),
          longLabel_(longLabel) {}

	const String& name() const { return name_ ; }
	const String& shortLabel() const { return shortLabel_ ; }
	const String& longLabel() const { return longLabel_ ; }
private:
	String name_;
	String shortLabel_;
	String longLabel_;
};

inline bool operator<(const HeaderItemInfo& left, const HeaderItemInfo& right){
	return left.name() < right.name();
}


// Definition of supported Page Header Items
class PageHeaderItemsDef {
public:
	enum class Item {
		begin=0,
	    // Plot Data Storage
	    Filename=begin,
		YColumns,
		// Queries on MeasurementSet
		ms_begin,
		// ---- Queries on Observation table
		Obs_Start_Date=ms_begin,
		Obs_Start_Time,
		Obs_Observer,
		Obs_Project,
		Obs_Telescope_Name,
		// ---- "Target" information
		Target_Name,
		Target_Direction,
		ms_end,
		// Queries on CalTable: currently none
		// cal_begin=ms_end, cal_end=ms_end, end=cal_end
		end=ms_end,
		n_items=end
	};
	static constexpr size_t n_items = static_cast<size_t>(Item::n_items);
	static constexpr size_t size() { return n_items; }

	static const map<HeaderItemInfo,Item>& info2Item_() {
		static const map<HeaderItemInfo,Item> info2item = {
	      // Name                   Short label         Long label                   Item
		  { {"filename"           , "File"            , "Filename"               } , Item::Filename           },
		  { {"y_columns"          , "Y Column(s)"     , "Y Column(s)"            } , Item::YColumns           },
		  { {"obs_start_date"     , "Observation date", "Observation Start Date" } , Item::Obs_Start_Date     },
		  { {"obs_start_time"     , "Observation time", "Observation Start Time" } , Item::Obs_Start_Time     },
		  { {"obs_observer"       , "Observer"        , "Observer"               } , Item::Obs_Observer       },
		  { {"obs_project"        , "Project ID"      , "Project ID"             } , Item::Obs_Project        },
		  { {"obs_telescope_name" , "Telescope"       , "Telescope Name"         } , Item::Obs_Telescope_Name },
		  { {"target_name"        , "Target name"     , "Target Name"            } , Item::Target_Name        },
		  { {"target_direction"   , "Target RA,Dec"   , "Target Direction"       } , Item::Target_Direction   }
		};
		return info2item;
	}

	struct Reverse {
		Reverse() {}
		void operator()(pair<HeaderItemInfo,Item> p)
		{
			item2Info[p.second] = p.first;
		}
		const map<Item,HeaderItemInfo> & result() { return item2Info; }
		map<Item,HeaderItemInfo> item2Info;
	};

	static map<Item,HeaderItemInfo> reverseMap(const map<HeaderItemInfo,Item>& map ) {
		return for_each(map.begin(),map.end(),Reverse()).result();
	}

	static const map<Item,HeaderItemInfo>& item2Info_() {
		static const auto item2info = reverseMap(info2Item_());
		return item2info;
	}

	struct GetName {
		GetName() {}
		void operator()(pair<HeaderItemInfo,Item> p)
		{
			toItem[p.first.name()] = p.second;
		}
		map<String,Item> toItem;
	};

	static const map<String,Item>& name2Item_() {
		static const auto & info = info2Item_();
		static const auto name2item = for_each(info.begin(),info.end(),GetName()).toItem ;
		return name2item;
	}

	static const array<Item,n_items>& items();
	static const array<String,n_items>& itemNames();
	static bool isItemName(const String& s);

	static const map<HeaderItemInfo,Item>& info2Item;
	static const map<Item,HeaderItemInfo>& item2Info;
	static const map<String,Item>& name2Item;

private:
	PageHeaderItemsDef();
};

// Ordered selection of distinct (no duplicates) page header items.
// Enforce a maximum of 10 header items per plot.
class PageHeaderItems {
public:
	static constexpr size_t kMaxItems = 10;
	using HeaderItem = PageHeaderItemsDef::Item;

	PageHeaderItems(const String& items=String(),const char sep=',');
	const std::vector<HeaderItem>& items() const;
	void setItems(const String& items=String(),const char sep=',');
	void clear();
	bool append(const String& item);
	bool append(HeaderItem item);
	bool remove(HeaderItem item);
	bool full() const;
	bool selected(HeaderItem item) const;
	pair<bool,size_t> index(HeaderItem item) const;

    bool operator==(const PageHeaderItems& other) const {
        return items_ == other.items() ; }
    bool operator!=(const PageHeaderItems& other) const {
        return !(operator==(other)); }

private:
	std::vector<HeaderItem> items_;
	map<HeaderItem,size_t> item2index_;
};


}


#endif /* PLOTMSPAGEHEADERPARAM_H_ */
