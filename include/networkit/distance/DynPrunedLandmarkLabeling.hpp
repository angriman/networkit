#ifndef NETWORKIT_DISTANCE_DYN_PRUNED_LANDMARK_LABELING_HPP_
#define NETWORKIT_DISTANCE_DYN_PRUNED_LANDMARK_LABELING_HPP_

// networkit-format

#include <algorithm>
#include <iostream>
#include <limits>
#include <queue>
#include <vector>

#ifndef NDEBUG
#include <memory>

#include <networkit/distance/APSP.hpp>
#endif

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class DynPrunedLandmarkLabeling : public Algorithm {
    static constexpr auto infDist = static_cast<count>(std::numeric_limits<count>::max());
    static constexpr auto maxTS = std::numeric_limits<uint8_t>::max();

    struct Label {
        node node_ = none;
        count distance_ = infDist;

        Label() = default;
        Label(node node_, count distance_) : node_(node_), distance_(distance_) {}
    };

public:
    DynPrunedLandmarkLabeling(const Graph &G);
    void run() override;
    count query(node u, node v, bool prefixal = false, node upperBound = none) const;

    void addEdge(node u, node v);
    void removeEdge(node u, node v);

private:
    const Graph *G;
    std::vector<node> nodesSortedByDegree, updatedNodes, rankOfNode, affectedURanks, affectedVRanks;
    std::vector<bool> visited;
    std::vector<uint8_t> affected;
    std::vector<count> prevDistX, prevDistY, newDistX, newDistY;
    // affected timestamp needs to start from an odd number (last used before reset will be 255)
    uint8_t tsAffected = 1;
    std::vector<std::vector<Label>> labelsOut, labelsIn;
    std::vector<Label> labelsUCopy, labelsVCopy;

    void updateAffected();
    void prunedBFS(node root, node rootRank, bool reverse = false);
    void sortUpdatedLabels();
    void resumeBFS(node k, node startNode, node otherNode, count level, bool reverse = false);

    void computePrevAndNewDists(node x, node y, std::vector<count> &prevDist,
                                std::vector<count> &newDist);

    void computeAffectedNodes(node x, node y, std::vector<node> &affectedNodesRanks,
                              std::vector<count> &prevDist, std::vector<count> &newDist,
                              std::vector<count> &prevDist2);

    void computeNewHubs(node rootRank, uint8_t tsAffectedU, uint8_t tsAffectedV,
                        bool rootInAffectedU);

    void removeAffectedHubs(uint8_t tsAffectedU, uint8_t tsAffectedV);

    node smallestHub(node u, node v) const;

#ifndef NDEBUG
    std::unique_ptr<APSP> apsp, apspNew;
    void checkLabels() const;
    void checkAffectedHubs(uint8_t tsAffectedU, uint8_t tsAffectedV) const;
#endif
};

} // namespace NetworKit
#endif // NETWORKIT_DISTANCE_DYN_PRUNED_LANDMARK_LABELING_HPP_
