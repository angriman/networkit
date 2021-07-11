#include <limits>
#include <omp.h>
#include <queue>
#include <stdexcept>
#include <unordered_set>

#include <networkit/auxiliary/VectorComparator.hpp>
#include <networkit/centrality/GroupHarmonicCloseness.hpp>
#include <networkit/centrality/GroupHarmonicClosenessLocalSearch.hpp>
#include <networkit/graph/GraphTools.hpp>

#include <tlx/container/d_ary_addressable_int_heap.hpp>
#include <tlx/unused.hpp>

namespace NetworKit {

namespace {

template <class WeightType>
class GroupHarmonicClosenessLocalSearchImpl final
    : public GroupHarmonicClosenessLocalSearch::GroupHarmonicClosenessLocalSearchInterface {

    static constexpr WeightType infDistance = std::numeric_limits<WeightType>::max();

public:
    template <class InputIt>
    GroupHarmonicClosenessLocalSearchImpl(const Graph &G, InputIt first, InputIt last);

    GroupHarmonicClosenessLocalSearchImpl(const Graph &G, size_t k);

    void run() override;

private:
    const Graph *G;

    // Distance from the group to all other nodes, one vector per thread
    std::vector<std::vector<WeightType>> distFromGroupGlobal;
    // Distance from the group to all other nodes: current and after a swap
    std::vector<WeightType> distFromGroup, nextDistFromGroup;

    void computeDistFromGroup(std::vector<WeightType> &dist);
};

template <class WeightType>
constexpr WeightType GroupHarmonicClosenessLocalSearchImpl<WeightType>::infDistance;

template <class WeightType>
template <class InputIt>
GroupHarmonicClosenessLocalSearchImpl<WeightType>::GroupHarmonicClosenessLocalSearchImpl(
    const Graph &G, InputIt first, InputIt last)
    : GroupHarmonicClosenessLocalSearchInterface(first, last), G(&G) {

    if (group.size() < 2)
        throw std::runtime_error("Error, the group needs to have at least two vertices");
}

template <class WeightType>
GroupHarmonicClosenessLocalSearchImpl<WeightType>::GroupHarmonicClosenessLocalSearchImpl(
    const Graph &G, size_t k)
    : GroupHarmonicClosenessLocalSearchInterface(), G(&G) {

    if (k < 2 || k >= G.numberOfNodes())
        throw std::runtime_error("Error, k needs to be in [2, n - 1]");

    // Initialize group with greedy algorithm
    GroupHarmonicCloseness ghc(G, k);
    ghc.run();
    group = ghc.groupMaxHarmonicCloseness();
}

template <class WeightType>
void GroupHarmonicClosenessLocalSearchImpl<WeightType>::run() {
    const auto n = G->upperNodeIdBound();
    distFromGroupGlobal.resize(omp_get_max_threads(), std::vector<WeightType>(n));
    nextDistFromGroup.resize(n);

    distFromGroup.resize(n, infDistance);
    computeDistFromGroup(distFromGroup);

    // Whether or not a node is in the group
    std::vector<bool> inGroup(n);
    for (node u : group)
        inGroup[u] = true;

    // Decrease and increase in closeness after a swap
    double scoreDecrease, scoreIncrease;

    size_t k = group.size(); // Only used for debugging
    tlx::unused(k);

    // Whether or not a swap led to an improvement of the centrality
    bool anyProgreess;
    do {
        anyProgreess = false;
        node lastRemoved = none;

        assert(group.size() == k);

        for (size_t i = 0; i < group.size(); ++i) {
            lastRemoved = group[i];
            std::swap(group[i], group.back());
            group.pop_back();

            computeDistFromGroup(nextDistFromGroup);

            // Compute the decrease in closeness due to the removal
            G->parallelForNodes([&](node u) {
                if (distFromGroup[u] != nextDistFromGroup[u]) {
                    assert(nextDistFromGroup[u] > 0);
                    if (distFromGroup[u] > 0)
#pragma omp atomic update
                        scoreDecrease += 1. / static_cast<double>(distFromGroup[u]);
                    if (nextDistFromGroup[u] != infDistance)
#pragma omp atomic update
                        scoreDecrease -= 1. / static_cast<double>(nextDistFromGroup[u]);
                }
            });


        }
    } while (anyProgreess);

    hasRun = true;
}

template <> // Unweighted graphs, run a BFS
void GroupHarmonicClosenessLocalSearchImpl<count>::computeDistFromGroup(std::vector<count> &dist) {
    std::queue<node> q;
    for (node u : group) {
        q.push(u);
        dist[u] = 0;
    }

    do {
        const node u = q.front();
        q.pop();

        G->forNeighborsOf(u, [&](node v) {
            if (dist[u] == infDistance) {
                dist[v] = dist[u] + static_cast<count>(defaultEdgeWeight);
                q.push(v);
            }
        });
    } while (!q.empty());
}

template <> // Weighted graphs, run Dijkstra
void GroupHarmonicClosenessLocalSearchImpl<edgeweight>::computeDistFromGroup(
    std::vector<edgeweight> &dist) {
    std::fill(dist.begin(), dist.end(), infDistance);

    tlx::d_ary_addressable_int_heap<node, 2, Aux::LessInVector<edgeweight>> prioQ{dist};
    prioQ.reserve(G->upperNodeIdBound());

    for (node u : group) {
        dist[u] = 0;
        prioQ.push(u);
    }

    do {
        const node u = prioQ.extract_top();

        G->forNeighborsOf(u, [&](node v, edgeweight ew) {
            const edgeweight newDistV = dist[u] + ew;
            if (newDistV < dist[v]) {
                dist[v] = newDistV;
                prioQ.update(v);
            }
        });
    } while (!prioQ.empty());
}

} // namespace

GroupHarmonicClosenessLocalSearch::GroupHarmonicClosenessLocalSearch(const Graph &G,
                                                                     const std::vector<node> &group)
    : weighted(G.isWeighted()) {

    if (weighted)
        impl = std::make_unique<GroupHarmonicClosenessLocalSearchImpl<edgeweight>>(G, group.begin(),
                                                                                   group.end());
    else
        impl = std::make_unique<GroupHarmonicClosenessLocalSearchImpl<count>>(G, group.begin(),
                                                                              group.end());
}

GroupHarmonicClosenessLocalSearch::GroupHarmonicClosenessLocalSearch(const Graph &G, size_t k)
    : weighted(G.isWeighted()) {

    if (weighted)
        impl = std::make_unique<GroupHarmonicClosenessLocalSearchImpl<edgeweight>>(G, k);
    else
        impl = std::make_unique<GroupHarmonicClosenessLocalSearchImpl<count>>(G, k);
}

void GroupHarmonicClosenessLocalSearch::run() {
    impl->run();
    hasRun = true;
}

std::vector<node> GroupHarmonicClosenessLocalSearch::groupMaxHarmonicCloseness() const {
    assureFinished();
    return impl->group;
}

count GroupHarmonicClosenessLocalSearch::numberOfIterations() const {
    assureFinished();
    return impl->nIterations;
}

} // namespace NetworKit
