#include <networkit/components/DynConnectedComponents.hpp>

namespace NetworKit {

DynConnectedComponents::DynConnectedComponents(const Graph &G)
    : impl(new DynConnectedComponentsImpl{G}) {}

void DynConnectedComponents::update(GraphEvent e) {
    impl->update(e);
}

void DynConnectedComponents::updateBatch(const std::vector<GraphEvent> &batch) {
    impl->update(e);
}
} // namespace NetworKit
