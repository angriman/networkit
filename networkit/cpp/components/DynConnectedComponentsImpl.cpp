#include <networkit/components/DynConnectedComponentsImpl.hpp>

namespace NetworKit {
namespace DynConnectedComponentsDetails {

DynConnectedComponentsImpl::DynConnectedComponentsImpl(const Graph &G, bool weaklyCC)
    : ConnectedComponentsImpl(G, weaklyCC) {
    components.resize(G.upperNodeIdBound(), none);
    tmpDistances.resize(G.upperNodeIdBound(), none);
    indexEdges();
    isTree.assign(edgesMap.size(), false);
    hasRun = false;
}

void DynConnectedComponentsImpl::update(GraphEvent event) {
    assureFinished();

    if (event.type == GraphEvent::EDGE_ADDITION)
        addEdge(event.u, event.v);
    else if (event.type == GraphEvent::EDGE_REMOVAL)
        removeEdge(event.u, event.v);
    else
        WARN("Unsupported edge update");
}

void DynConnectedComponentsImpl::updateBatch(const std::vector<GraphEvent> &batch) {
    for (const auto &e : batch)
        update(e);
}

void DynConnectedComponentsImpl::addEdge(node u, node v) {

    std::pair<bool, edgeid> updateResult = updateMapAfterInsertion(u, v);

    // If u and v are already in the same component, we
    // don't have to do anything
    index maxComp = std::max(components[u], components[v]);
    index minComp = std::min(components[u], components[v]);

    if (maxComp == minComp) {
        if (!updateResult.first)
            isTree.push_back(false);
        return;
    }

    // In the other case, we can merge the two components in an undirected
    // graph merge components
    G->parallelForNodes([&](node w) {
        if (components[w] == maxComp)
            components[w] = minComp;
    });

    compSize.find(minComp)->second += compSize.find(maxComp)->second;
    compSize.erase(maxComp);
    componentIds.push(maxComp);

    if (updateResult.first)
        isTree[updateResult.second] = true;
    else
        isTree.push_back(true);
}

void DynConnectedComponentsImpl::removeEdge(node u, node v) {
    const auto eid = edgesMap.at(makeEdge(u, v));

    // This edge removal does not split two components. Nothing to do.
    if (!isTree[eid])
        return;

    isTree[eid] = false; // for coherence, we mark this edge as not valid
    std::fill(tmpDistances.begin(), tmpDistances.end(), none);
    index nextId = nextAvailableComponentId(false);

    std::vector<node> newCmp(components);
    newCmp[u] = nextId;
    count newCmpSize = 0;

    std::queue<node> q;
    q.push(u);
    tmpDistances[u] = 0;

    bool connected = false;

    // Berform BFS from v to check if v reaches u
    do {
        const node s = q.front();
        q.pop();
        ++newCmpSize;

        count d = tmpDistances[s] + 1;
        const auto visitNeighbor = [&](node w) -> bool {
            if (tmpDistances[w] == none) {
                tmpDistances[w] = d;
                if (w == v) { // Found another path from u to v
                    // Backtracks the path from v to u and marks all its
                    // nodes as part of the spanning tree
                    reverseBFS(u, v);
                    connected = true;
                    return true; // Exit from the loop
                }

                newCmp[w] = nextId;
                q.push(w);
            }
            return false;
        };

        // Enqueue not visited neighbors
        for (const node w : G->neighborRange(s))
            if (visitNeighbor(w))
                break;

        if (connected)
            break;

        if (weaklyCC) {
            for (const node w : G->inNeighborRange(s))
                if (visitNeighbor(w))
                    break;

            if (connected)
                break;
        }

    } while (!q.empty());

    if (!connected) {
        index nextId = nextAvailableComponentId();
        compSize.find(components[u])->second -= newCmpSize;
        compSize.emplace(nextId, newCmpSize);
        components = newCmp;
    }
}

void DynConnectedComponentsImpl::reverseBFS(node u, node v) {
    std::queue<node> q;
    q.push(v);

    count d = tmpDistances[v];
    count level = 1;

    do {
        node s = q.front();
        q.pop();

        bool nextEdgeFound = false;
        const auto visitNeighbor = [&](node w) -> void {
            if (w == u) {
                isTree[edgesMap.find(makeEdge(w, s))->second] = true;
                nextEdgeFound = true;
                return;
            }

            if ((tmpDistances[w] != none) && (d == tmpDistances[w] + level)) {
                isTree[edgesMap.find(makeEdge(w, s))->second] = true;
                nextEdgeFound = true;
                q.push(w);
            }
        };

        for (const node w : G->neighborRange(s)) {
            visitNeighbor(w);
            if (nextEdgeFound)
                break;
        }

        if (weaklyCC && !nextEdgeFound)
            for (const node w : G->inNeighborRange(s)) {
                visitNeighbor(w);
                if (nextEdgeFound)
                    break;
            }
        ++level;

    } while (!q.empty());
}
index DynConnectedComponentsImpl::nextAvailableComponentId(bool eraseId) {
    if (componentIds.empty())
        return compSize.size();
    index result = componentIds.front();
    if (eraseId)
        componentIds.pop();
    return result;
}
void DynConnectedComponentsImpl::indexEdges() {
    edgeid eid = 0;
    G->forEdges([&](node u, node v) {
        if (edgesMap.find(makeEdge(u, v)) == edgesMap.end())
            edgesMap.emplace(makeEdge(u, v), eid++);
    });
}

std::pair<bool, edgeid> DynConnectedComponentsImpl::updateMapAfterInsertion(node u, node v) {
    const auto it = edgesMap.find(makeEdge(u, v));

    if (it == edgesMap.end()) {
        edgeid newId = edgesMap.size();
        // Adding edge never deleted before
        edgesMap.emplace(makeEdge(u, v), newId);
        return {false, none};
    }

    return {true, it->second};
}

Edge DynConnectedComponentsImpl::makeEdge(node u, node v) const noexcept {
    return {std::min(u, v), std::max(u, v)};
}

} // namespace DynConnectedComponentsDetails
} // namespace NetworKit
