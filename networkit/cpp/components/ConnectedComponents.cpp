/*
 * ConnectedComponents.cpp
 *
 *  Created on: Dec 16, 2013
 *      Author: cls
 */

// networkit-format

#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

ConnectedComponents::ConnectedComponents(const Graph &G) : ConnectedComponentsGeneral(G) {}

count ConnectedComponents::numberOfComponents() const {
    return impl->numberOfComponents();
}

count ConnectedComponents::componentOfNode(node u) const {
    return impl->componentOfNode(u);
}

void ConnectedComponents::run() {
    impl->run();
}

Partition ConnectedComponents::getPartition() const {
    return impl->getPartition();
}

std::vector<std::vector<node>> ConnectedComponents::getComponents() const {
    return impl->getComponents();
}

std::map<index, count> ConnectedComponents::getComponentSizes() const {
    return impl->getComponentSizes();
}

Graph ConnectedComponents::extractLargestConnectedComponent(const Graph &G, bool compactGraph) {
    ConnectedComponentsDetails::ConnectedComponentsImpl ccImpl(G);
    return ccImpl.extractLargestConnectedComponent(compactGraph);
}

} // namespace NetworKit
