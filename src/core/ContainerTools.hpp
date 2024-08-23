#ifndef __OXSX__CONTAINER_TOOLS__
#define __OXSX__CONTAINER_TOOLS__
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <cassert>
#include <Formatter.hpp>

namespace ContainerTools
{

    template <typename T1, typename T2>
    std::map<typename T1::value_type, typename T2::value_type>
    CreateMap(const T1 &v1_, const T2 &v2_)
    {
        assert(v1_.size() == v2_.size());
        typedef typename T1::const_iterator It2;
        typedef typename T2::const_iterator It3;

        std::map<typename T1::value_type, typename T2::value_type> map;

        It2 it1 = v1_.begin();
        It3 it2 = v2_.begin();

        for (; it1 != v1_.end();)
        {
            map[*it1] = *it2;
            ++it1;
            ++it2;
        }
        return map;
    }

    template <typename T1>
    std::set<typename T1::key_type> GetKeys(const T1 &mp_)
    {
        typedef typename T1::const_iterator MapIt;
        typedef typename T1::key_type KeyType;

        std::set<KeyType> set;
        for (MapIt it = mp_.begin(); it != mp_.end(); ++it)
        {
            set.insert(it->first);
        }
        return set;
    }

    template <typename T1, typename T2>
    std::vector<typename T1::mapped_type> GetValues(const T1 &mp_,
                                                    const T2 &keys_)
    {
        // throws std::out_of_range if there is no matching key
        typedef typename T1::mapped_type SecType;
        typedef typename T2::const_iterator It;

        std::vector<SecType> vals;
        for (It it = keys_.begin(); it != keys_.end(); ++it)
        {
            vals.push_back(mp_.at(*it));
        }

        return vals;
    }

    template <typename T1>
    std::vector<typename T1::mapped_type> GetValues(const T1 &mp_)
    {

        // throws std::out_of_range if there is no matching key
        typedef typename T1::const_iterator It;

        std::vector<typename T1::mapped_type> vals;
        for (It it = mp_.begin(); it != mp_.end(); ++it)
        {
            vals.push_back(it->second);
        }

        return vals;
    }

    template <typename T1, typename T2, typename T3>
    void SetValues(T1 &mp_, const T2 &keys_, const T3 &values_)
    {
        assert(values_.size() == keys_.size());

        typedef typename T2::const_iterator It2;
        typedef typename T3::const_iterator It3;

        It2 it2 = keys_.begin();
        It3 it3 = values_.begin();

        for (; it2 != keys_.end();)
        {
            mp_[*it2] = *it3;
            ++it2;
            ++it3;
        }
    }

    template <typename T1, typename T2>
    bool HasSameKeys(const T1 &m1_, const T2 &m2_)
    {
        return GetKeys(m1_) == GetKeys(m2_);
    }

    template <typename T1, typename T2>
    bool HasKey(const std::map<T1, T2> &mp_, const T1 &key_)
    {
        return mp_.find(key_) != mp_.end();
    }

    template <typename T1, typename T2>
    std::string CompareKeys(const T1 &m1_, const T2 &m2_, const std::string name1_ = "p1", const std::string &name2_ = "p2")
    {
        // nothing doing
        if (HasSameKeys(m1_, m2_))
            return name1_ + " and " + name2_ + " have the same keys";

        typename std::vector<typename T1::key_type> missingFrom2;

        typedef typename T1::const_iterator mapIt;

        for (mapIt it = m1_.begin(); it != m1_.end(); ++it)
        {
            if (!HasKey(m2_, it->first))
                missingFrom2.push_back(it->first);
        }

        std::string returnStr = "Keys in " + name1_ + " but missing from " + name2_ + " : \n";
        for (size_t i = 0; i < missingFrom2.size(); i++)
        {
            returnStr += std::string(missingFrom2.at(i));
            returnStr += ",";
        }

        returnStr += "\n";
        return returnStr;
    }

    template <typename T1>
    std::string ToString(const T1 &c_, const std::string &delimit_ = ", ")
    {
        Formatter f;
        for (typename T1::const_iterator it = c_.begin(); it != c_.end(); ++it)
            f << *it << delimit_;
        return f.str();
    }

    // thanks to https://stackoverflow.com/questions/236129/most-elegant-way-to-split-a-string
    template <typename Out>
    void Split(const std::string &s, char delim, Out result)
    {
        std::stringstream ss;
        ss.str(s);
        std::string item;
        while (std::getline(ss, item, delim))
        {
            if (!item.empty())
                *(result++) = item;
        }
    }

    std::vector<std::string>
    Split(const std::string &s, char delim = ' ');

} // namespace
#endif
