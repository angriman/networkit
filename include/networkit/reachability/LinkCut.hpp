#ifndef NETWORKIT_REACHABILITY_LINK_CUT_HPP
#define NETWORKIT_REACHABILITY_LINK_CUT_HPP

#include <vector>

#include <networkit/graph/Graph.hpp>

namespace NetworKit {
class LinkCut final {
public:
    LinkCut(const Graph &G);

    void doCuts(count n);

private:
    const Graph *G;
    node root;

    std::vector<node> sequence, parent, seenNodes;
    std::vector<count> edgeScores, pathLength;
    std::vector<EdgeWithId> nonSTEdges;

    std::vector<short> status, edgeInST;

    void sampleUST();
    void doLinkCut();
    void buildNodeSequence();
};

} // namespace NetworKit
#endif
