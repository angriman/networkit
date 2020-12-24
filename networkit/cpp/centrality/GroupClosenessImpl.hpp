// networkit-format

#ifndef NETWORKIT_CENTRALITY_GROUP_CLOSENESS_IMPL_HPP_
#define NETWORKIT_CENTRALITY_GROUP_CLOSENESS_IMPL_HPP_

#include <limits>
#include <vector>

#include <networkit/graph/BFS.hpp>
#include <networkit/graph/Dijkstra.hpp>

#include <tlx/container/d_ary_addressable_int_heap.hpp>

namespace NetworKit {
namespace GroupClosenessDetails {

template <class Weight>
class GroupClosenessImpl final {
    static constexpr Weight infDist = std::numeric_limits<Weight>::max();

public:
    GroupClosenessImpl(const Graph &G, count k = 1, Weight maxDist = 0);

    void run();

    const std::vector<node> &groupMaxCloseness() const { return group; }

    template <class InputIt>
    static double scoreOfGroup(const Graph &graph, InputIt first, InputIt last);

private:
    const Graph *G;
    const count k;
    const Weight maxDist;
    Weight minEdgeWeight = 1;

    std::vector<node> group, nearerNodes;
    std::vector<std::vector<node>> predGlobal;
    std::vector<std::vector<node>> nearerNodesGlobal;
    std::vector<Weight> distFromGroup, farness, marginalGain, lowestDist;
    std::vector<std::vector<Weight>> distGlobal;
    std::vector<std::vector<bool>> visitedGlobal;

    struct Greater {
        Greater(const std::vector<Weight> &vec) : vec(vec) {}
        bool operator()(node x, node y) const noexcept { return vec[x] > vec[y]; }

    private:
        const std::vector<Weight> &vec;
    };

    struct Less {
        Less(const std::vector<Weight> &vec) : vec(vec) {}
        bool operator()(node x, node y) const noexcept { return vec[x] < vec[y]; }

    private:
        const std::vector<Weight> &vec;
    };

    tlx::d_ary_addressable_int_heap<node, 2, Greater> candidateNodesPQ{Greater(marginalGain)};
    std::vector<tlx::d_ary_addressable_int_heap<node, 2, Less>> dijkstraHeaps;

    struct PrunedSSSPResult {
        Weight farness;
        bool pruned;
        PrunedSSSPResult(bool pruned, Weight farness) : farness(farness), pruned(pruned) {}
    };

    PrunedSSSPResult prunedSSSPEmptyGroup(node source, double highestClosenessScore);
    node topClosenessNode();

    void computeFarnessLowerBound();

    node findNodeWithHighestMarginalGain();
    Weight computeMargGain(node source);

#ifdef NETWORKIT_SANITY_CHECKS
    void checkTopNode(node u, const std::vector<Weight> &dist, Weight computedFarness) const;
    void checkDistFromGroup() const;
#endif // NETWORKIT_SANITY_CHECKS
};

template <class Weight>
constexpr Weight GroupClosenessImpl<Weight>::infDist;

template <>
template <class InputIt>
double GroupClosenessImpl<count>::scoreOfGroup(const Graph &graph, InputIt first, InputIt last) {
    count distSum = 0;
    Traversal::BFSfrom(graph, first, last, [&distSum](node, count dist) { distSum += dist; });
    return distSum == 0 ? 0
                        : static_cast<double>(graph.numberOfNodes()) / static_cast<double>(distSum);
}

template <>
template <class InputIt>
double GroupClosenessImpl<edgeweight>::scoreOfGroup(const Graph &graph, InputIt first,
                                                    InputIt last) {
    edgeweight distSum = 0;
    Traversal::DijkstraFrom(graph, first, last,
                            [&distSum](node, edgeweight dist) { distSum += dist; });
    return distSum == 0 ? 0 : static_cast<edgeweight>(graph.numberOfNodes()) / distSum;
}
} // namespace GroupClosenessDetails
} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_GROUP_CLOSENESS_IMPL_HPP_
