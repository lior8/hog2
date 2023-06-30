#include <iostream>
#include <vector>
#include <algorithm>
#include <set>

template<typename Key, typename Value>
class FastDict {
private:
    std::vector<std::pair<Key, Value>> data;

public:
    void insert(const Key& key, const Value& value) {
        auto it = std::lower_bound(data.begin(), data.end(), key,
            [](const std::pair<Key, Value>& a, const Key& b) {
                return a.first < b;
            });

        data.insert(it, std::make_pair(key, value));
    }

    Value& operator[](const Key& key) {
        auto it = std::lower_bound(data.begin(), data.end(), key,
            [](const std::pair<Key, Value>& a, const Key& b) {
                return a.first < b;
            });

        if (it == data.end() || it->first != key) {
            it = data.insert(it, std::make_pair(key, Value{}));
        }

        return it->second;
    }

    void erase(const Key& key) {
        auto it = std::lower_bound(data.begin(), data.end(), key,
            [](const std::pair<Key, Value>& a, const Key& b) {
                return a.first < b;
            });

        if (it != data.end() && it->first == key) {
            data.erase(it);
        }
    }

    bool contains(const Key& key) const {
        auto it = std::lower_bound(data.begin(), data.end(), key,
            [](const std::pair<Key, Value>& a, const Key& b) {
                return a.first < b;
            });

        return (it != data.end() && it->first == key);
    }

    std::size_t size() const {
        return data.size();
    }
};


template<class Value>
class Cache
{
public:
    std::vector<std::set<Value>> data;
    int cap;
    void insert(const uint64_t &key, const Value &value) {
        auto left = key%cap;
        data[left].insert(value);
    }
    void remove(const uint64_t &key, const Value &value)
    {
        auto left = key%cap;
        data[left].erase(value);
    }
    bool contains(const uint64_t &key, const Value &value)
    {
        auto left = key%cap;
        return data[left].find(value) != data[left].end();
    }
};