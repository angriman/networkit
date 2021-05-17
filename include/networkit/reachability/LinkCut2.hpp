#ifndef NETWORKIT_REACHABILITY_LINK_CUT2_HPP
#define NETWORKIT_REACHABILITY_LINK_CUT2_HPP

#include <vector>

#include <networkit/graph/Graph.hpp>

namespace NetworKit {
class LinkCut2 final {
public:
    LinkCut2(const Graph &G);

    std::vector<double> simulation(count reps, count cutsPerRep);
private:
    const Graph *G;
    node root;

    std::vector<node> sequence, parent, seenNodes;
    std::vector<count> pathLength;
    std::vector<EdgeWithId> nonSTEdges;

    std::vector<short> status;
    std::vector<bool> edgeInST;

    void sampleUST();
    void doLinkCut();
    void buildNodeSequence();
};

} // namespace NetworKit

#endif
