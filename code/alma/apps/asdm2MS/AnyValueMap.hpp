/*
 * AnyValueMap.hpp
 *
 *  Created on: Jun 3, 2013
 *      Author: alvaro
 */

// Only used by parasdm2MS, could possibly be replaced by a casacore::Record?
// boost::any is in c++17

#ifndef ANYVALUEMAP_HPP_
#define ANYVALUEMAP_HPP_

#include <map>
#include <boost/any.hpp>

template <class T>
class AnyValueMap {

public:
    AnyValueMap(){}

    virtual ~AnyValueMap(){}

private:
    std::map<T, boost::any> container_;

    typedef typename std::map<T, boost::any>::iterator map_iterator;
    typedef typename std::map<T, boost::any>::const_iterator map_const_iterator;

public:

    bool containsKey(const T key) const
    {
        return container_.find(key) != container_.end();
    }

    bool remove(const T key)
    {
        std::map_iterator it = container_.find(key);
        if(it != container_.end())
        {
            container_.erase(it);
            return true;
        }
        return false;
    }

    template <class V>
    V getValue(const T key, const V defaultValue) const
    {
        std::map_const_iterator it = container_.find(key);
        if(it != container_.end())
        {
            return boost::any_cast<V>(it->second);
        }
        return defaultValue;
    }

    template <class V>
    V getValue(const T key) const
    {
        return boost::any_cast<V>(container_.at(key));
    }

    template <class V>
    void setValue(const T key, const V value)
    {
        container_[key] = value;
    }
};

#endif /* ANYVALUEMAP_HPP_ */
