/*
 *  GroupHarmonicClosenessLocalSearchImpl.hpp
 *
 *  Created on: 20.12.2020
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

// networkit-format

#ifndef NETWORKIT_CENTRALITY_GROUP_HARMONIC_CLOSENESS_LOCAL_SEARCH_IMPL_HPP_
#define NETWORKIT_CENTRALITY_GROUP_HARMONIC_CLOSENESS_LOCAL_SEARCH_IMPL_HPP_

#include <limits>
#include <vector>

#include <networkit/graph/Graph.hpp>

#include <tlx/container/d_ary_addressable_int_heap.hpp>

namespace NetworKit {
namespace GroupHarmonicClosenessLocalSearchDetails {

template <class Weight>
class GroupHarmonicClosenessLocalSearchImpl final {
    static constexpr Weight infDist = std::numeric_limits<Weight>::max();
    GroupHarmonicClosenessLocalSearchImpl(const Graph &G);

public:
    GroupHarmonicClosenessLocalSearchImpl(const Graph &G, count k);

    template <class InputIt>
    GroupHarmonicClosenessLocalSearchImpl(const Graph &G, InputIt first, InputIt last);

    ~GroupHarmonicClosenessLocalSearchImpl();

    const std::vector<node> &groupMaxHarmonicCloseness() const {
        return group;
    }

    void run(count maxSwaps);

private:
    const Graph *G;
    std::vector<node> group;
    const count k;
    std::vector<bool> inGroup;
    std::vector<std::vector<bool>> visitedGlobal;
    std::vector<Weight> distFromGroup, newDistFromGroup;
    std::vector<std::vector<Weight>> distGlobal;

    struct Less {
        Less(const std::vector<Weight> &dist) : dist(dist) {}
        bool operator()(node u, node v) const noexcept {
            return dist[u] < dist[v];
        }
        private:
        const std::vector<Weight> &dist;
    };

    std::vector<tlx::d_ary_addressable_int_heap<node, 2, Less>> dijkstraHeaps;

    std::vector<node> generateRandomGroup();
    void computeDistFromGroup();

    struct PrunedSSSPResult {
        bool pruned;
        double score;
        PrunedSSSPResult(double score, bool pruned) : pruned(pruned), score(score);
    };

    PrunedSSSPResult prunedSSSP(node source, double maxDecrement);
};

template <class Weight>
constexpr Weight GroupHarmonicClosenessLocalSearchImpl<Weight>::infDist;

} // namespace GroupHarmonicClosenessLocalSearchDetails
} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_GROUP_HARMONIC_CLOSENESS_LOCAL_SEARCH_IMPL_HPP_
