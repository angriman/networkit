#ifndef NETWORKIT_CENTRALITY_GROUP_HARMONIC_CLOSENESS_LOCAL_SEARCH_HPP_
#define NETWORKIT_CENTRALITY_GROUP_HARMONIC_CLOSENESS_LOCAL_SEARCH_HPP_

#include <memory>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
class GroupHarmonicClosenessLocalSearch : public Algorithm {

public:
    GroupHarmonicClosenessLocalSearch(const Graph &G, size_t k = 1);

    template <class InputIt>
    GroupHarmonicClosenessLocalSearch(const Graph &G, InputIt first, InputIt last)
        : GroupHarmonicClosenessLocalSearch(G, std::vector<node>(first, last), false) {}

    GroupHarmonicClosenessLocalSearch(const Graph &G, const std::vector<node> &group);

    void run() override;

    std::vector<node> groupMaxHarmonicCloseness() const;

    count numberOfIterations() const;

    class GroupHarmonicClosenessLocalSearchInterface : public Algorithm {
    public:
        GroupHarmonicClosenessLocalSearchInterface() = default;

        template <class InputIt>
        GroupHarmonicClosenessLocalSearchInterface(InputIt first, InputIt last)
            : group(first, last) {}


        std::vector<node> group;
        count nIterations;
    };

private:
    const bool weighted;
    // Is always one between GroupHarmonicClosenessLocalSearchImpl<count/edgeweight>,
    // see implementation.
    std::unique_ptr<GroupHarmonicClosenessLocalSearchInterface> impl;
};
} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_GROUP_HARMONIC_CLOSENESS_LOCAL_SEARCH_HPP_
