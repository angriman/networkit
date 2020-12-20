/*
 *  GroupHarmonicClosenessLocalSearch.hpp
 *
 *  Created on: 20.12.2020
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

// networkit-format

#ifndef NETWORKIT_CENTRALITY_GROUP_HARMONIC_CLOSENESS_LOCAL_SEARCH_HPP_
#define NETWORKIT_CENTRALITY_GROUP_HARMONIC_CLOSENESS_LOCAL_SEARCH_HPP_

#include <limits>
#include <memory>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

namespace GroupHarmonicClosenessLocalSearchDetails {
template <class>
class GroupHarmonicClosenessLocalSearchImpl;
using GroupHarmonicClosenessLocalSearchImplUnweighted =
    GroupHarmonicClosenessLocalSearchImpl<count>;
using GroupHarmonicClosenessLocalSearchImplWeighted =
    GroupHarmonicClosenessLocalSearchImpl<edgeweight>;
} // namespace GroupHarmonicClosenessLocalSearchDetails

class GroupHarmonicClosenessLocalSearch final : public Algorithm {
public:
    GroupHarmonicClosenessLocalSearch(const Graph &G, count k);

    template <class InputIt>
    GroupHarmonicClosenessLocalSearch(const Graph &G, InputIt first, InputIt last);

    ~GroupHarmonicClosenessLocalSearch() override;

    void run() override;

    const std::vector<node> &groupMaxHarmonicCloseness() const;

    count maxSwaps = std::numeric_limits<count>::max();

private:
    const bool weighted;
    std::unique_ptr<
        GroupHarmonicClosenessLocalSearchDetails::GroupHarmonicClosenessLocalSearchImplUnweighted>
        implUnweighted;
    std::unique_ptr<
        GroupHarmonicClosenessLocalSearchDetails::GroupHarmonicClosenessLocalSearchImplWeighted>
        implWeighted;
};

template <class InputIt>
GroupHarmonicClosenessLocalSearch::GroupHarmonicClosenessLocalSearch(const Graph &G, InputIt first,
                                                                     InputIt last)
    : weighted(G.isWeighted()) {
    if (weighted)
        implWeighted =
            std::make_unique<GroupHarmonicClosenessLocalSearchDetails::
                                 GroupHarmonicClosenessLocalSearchImplWeighted>(G, first, last);
    else
        implUnweighted =
            std::make_unique<GroupHarmonicClosenessLocalSearchDetails::
                                 GroupHarmonicClosenessLocalSearchImplUnweighted>(G, first, last);
}

} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_GROUP_HARMONIC_CLOSENESS_LOCAL_SEARCH_HPP_
