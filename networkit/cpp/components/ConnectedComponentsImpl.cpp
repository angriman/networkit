// networkit-format

#include <cassert>
#include <unordered_set>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/components/ConnectedComponentsImpl.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {
namespace ConnectedComponentsDetails {

ConnectedComponentsImpl::ConnectedComponentsImpl(const Graph &G, bool weaklyCC)
    : G(&G), weaklyCC(weaklyCC) {
    if (!weaklyCC && G.isDirected())
        throw std::runtime_error(
            "Error, connected components of directed graphs cannot be "
            "computed, use StronglyConnectedComponents or WeaklyConnectedComponents instead.");
    if (weaklyCC && !G.isDirected())
        throw std::runtime_error("Error, weakly connected components of an undirected graph cannot "
                                 "be computed, use ConnectedComponents instead.");
    hasRun = false;
}

void ConnectedComponentsImpl::assureFinished() const {
    if (!hasRun)
        throw std::runtime_error("Error, call run() method first.");
}

count ConnectedComponentsImpl::numberOfComponents() const {
    assureFinished();
    return numComponents;
}

count ConnectedComponentsImpl::componentOfNode(node u) const {
    assert(component[u] != none);
    assureFinished();
    return component[u];
}

void ConnectedComponentsImpl::run() {
    component = Partition(G->upperNodeIdBound(), none);
    numComponents = 0;
    count visitedNodes = 0;

    std::queue<node> q;

    // perform breadth-first searches
    for (const node u : G->nodeRange()) {
        if (component[u] != none)
            continue;

        component.setUpperBound(numComponents + 1);
        const index c = numComponents;

        q.push(u);
        component[u] = c;

        do {
            const node v = q.front();
            q.pop();
            ++visitedNodes;

            const auto visitNeighbor = [&](node w) -> void {
                if (component[w] == none) {
                    q.push(w);
                    component[w] = c;
                }
            };

            // enqueue neighbors, set component
            G->forNeighborsOf(v, visitNeighbor);
            if (weaklyCC)
                G->forInNeighborsOf(v, visitNeighbor);

        } while (!q.empty());

        ++numComponents;

        if (visitedNodes == G->numberOfNodes())
            break;
    }

    hasRun = true;
}

Partition ConnectedComponentsImpl::getPartition() const {
    assureFinished();
    return this->component;
}

std::vector<std::vector<node>> ConnectedComponentsImpl::getComponents() const {
    assureFinished();

    std::vector<std::vector<node>> result(numComponents);
    G->forNodes([&](node u) { result[component[u]].push_back(u); });

    return result;
}

std::map<index, count> ConnectedComponentsImpl::getComponentSizes() const {
    assureFinished();
    return this->component.subsetSizeMap();
}

Graph ConnectedComponentsImpl::extractLargestConnectedComponent(bool compactGraph) const {
    if (G->isEmpty())
        return *G;

    ConnectedComponentsImpl cc(*G, weaklyCC);
    cc.run();

    const auto compSizes = cc.getComponentSizes();
    if (compSizes.size() == 1) {
        if (compactGraph)
            return GraphTools::getCompactedGraph(*G, GraphTools::getContinuousNodeIds(*G));
        return *G;
    }

    const auto largestCC = std::max_element(
        compSizes.begin(), compSizes.end(),
        [](const std::pair<index, count> &x, const std::pair<index, count> &y) -> bool {
            return x.second < y.second;
        });

    if (compactGraph) {
        std::unordered_map<node, node> continuousNodeIds;
        index nextId = 0;
        G->forNodes([&](const node u) {
            if (cc.componentOfNode(u) == largestCC->first)
                continuousNodeIds[u] = nextId++;
        });

        return GraphTools::getRemappedGraph(
            *G, largestCC->second, [&](const node u) { return continuousNodeIds[u]; },
            [&](const node u) { return cc.componentOfNode(u) != largestCC->first; });

    } else {
        Graph S(*G);
        const auto components = cc.getComponents();
        for (size_t i = 0; i < components.size(); ++i) {
            if (i != largestCC->first)
                for (const node u : components[i])
                    S.removeNode(u);
        }

        return S;
    }
}

} // namespace ConnectedComponentsDetails
} // namespace NetworKit
