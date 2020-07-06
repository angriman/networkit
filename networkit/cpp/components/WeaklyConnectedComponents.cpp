// networkit-format

#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

WeaklyConnectedComponents::WeaklyConnectedComponents(const Graph &G)
    : ConnectedComponentsGeneral(G, true) {}

count WeaklyConnectedComponents::numberOfComponents() const {
    return impl->numberOfComponents();
}

count WeaklyConnectedComponents::componentOfNode(node u) const {
    return impl->componentOfNode(u);
}

void WeaklyConnectedComponents::run() {
    impl->run();
}

Partition WeaklyConnectedComponents::getPartition() const {
    return impl->getPartition();
}

std::vector<std::vector<node>> WeaklyConnectedComponents::getComponents() const {
    return impl->getComponents();
}

std::map<index, count> WeaklyConnectedComponents::getComponentSizes() const {
    return impl->getComponentSizes();
}

Graph WeaklyConnectedComponents::extractLargestWeaklyConnectedComponent(const Graph &G,
                                                                        bool compactGraph) {
    ConnectedComponentsDetails::ConnectedComponentsImpl ccImpl(G, true);
    return ccImpl.extractLargestConnectedComponent(compactGraph);
}

} // namespace NetworKit
