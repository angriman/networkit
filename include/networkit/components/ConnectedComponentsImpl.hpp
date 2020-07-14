// networkit-format
#ifndef NETWORKIT_COMPONENTS_CONNECTED_COMPONENTS_IMPL_HPP_
#define NETWORKIT_COMPONENTS_CONNECTED_COMPONENTS_IMPL_HPP_

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {
namespace ConnectedComponentsDetails {

class ConnectedComponentsImpl {
public:
    ConnectedComponentsImpl(const Graph &G, bool weaklyCC = false);

    void run();
    count componentOfNode(node u) const;
    count numberOfComponents() const;
    Partition getPartition() const;
    std::map<index, count> getComponentSizes() const;
    std::vector<std::vector<node>> getComponents() const;
    Graph extractLargestConnectedComponent(bool compactGraph) const;

protected:
    const Graph *G;
    const bool weaklyCC;
    bool hasRun;
    Partition component;
    count numComponents;
    void assureFinished() const;
};

} // namespace ConnectedComponentsDetails
} // namespace NetworKit

#endif // NETWORKIT_COMPONENTS_CONNECTED_COMPONENTS_IMPL_HPP_
