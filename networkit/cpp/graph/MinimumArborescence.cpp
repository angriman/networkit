// networkit-format
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/graph/MinimumArborescence.hpp>

namespace NetworKit {

namespace {

// An edge of the current contracted graph. Ids identify it in the contraction history.
struct Arc {
    index from;
    index to;
    edgeweight weight;
    index id;
};

// An edge created by a contraction: selecting it means selecting the edge added, which enters a
// contracted cycle, instead of the edge replaced, which had been selected in the previous round.
struct Replacement {
    index added;
    index replaced;
};

} // namespace

MinimumArborescence::MinimumArborescence(const Graph &G, node root) : G(&G), root(root) {
    if (!G.isDirected()) {
        throw std::runtime_error("MinimumArborescence: the graph must be directed.");
    }
    if (!G.hasNode(root)) {
        throw std::runtime_error("MinimumArborescence: the root must be in the graph.");
    }
}

void MinimumArborescence::run() {
    Aux::SignalHandler handler;
    const count z = G->upperNodeIdBound();

    exists = false;
    totalWeight = 0;
    arborescence = Graph(z, true, true);
    for (node u = 0; u < z; ++u) {
        if (!G->hasNode(u)) {
            arborescence.removeNode(u);
        }
    }

    // Node ids may have gaps, the nodes of the contracted graphs are consecutive.
    std::vector<index> positionOf(z, none);
    count numberOfNodes = 0;
    G->forNodes([&](node u) -> void { positionOf[u] = numberOfNodes++; });

    std::vector<Arc> arcs, inputArcs;
    G->forEdges([&](node u, node v, edgeweight weight) -> void {
        if (u != v && v != root) {
            const Arc arc{positionOf[u], positionOf[v], weight, arcs.size()};
            arcs.push_back(arc);
            inputArcs.push_back({u, v, weight, arc.id});
        }
    });

    // Number of rounds that selected each edge, minus the ones that replaced it afterwards.
    std::vector<int64_t> selections(arcs.size(), 0);
    std::vector<Replacement> replacements;

    index contractedRoot = positionOf[root];
    std::vector<index> cheapest, component, visitedBy;

    while (true) {
        handler.assureRunning();

        cheapest.assign(numberOfNodes, none);
        for (index i = 0; i < arcs.size(); ++i) {
            const Arc &arc = arcs[i];
            if (arc.from != arc.to && arc.to != contractedRoot
                && (cheapest[arc.to] == none || arc.weight < arcs[cheapest[arc.to]].weight)) {
                cheapest[arc.to] = i;
            }
        }

        for (index u = 0; u < numberOfNodes; ++u) {
            if (u != contractedRoot) {
                if (cheapest[u] == none) {
                    hasRun = true;
                    return;
                }
                ++selections[arcs[cheapest[u]].id];
            }
        }

        // The cheapest incoming edges form a functional graph: find its cycles.
        component.assign(numberOfNodes, none);
        visitedBy.assign(numberOfNodes, none);
        count numberOfComponents = 0;
        for (index start = 0; start < numberOfNodes; ++start) {
            index u = start;
            while (u != contractedRoot && visitedBy[u] == none) {
                visitedBy[u] = start;
                u = arcs[cheapest[u]].from;
            }

            if (u != contractedRoot && visitedBy[u] == start) {
                for (index v = arcs[cheapest[u]].from; v != u; v = arcs[cheapest[v]].from) {
                    component[v] = numberOfComponents;
                }
                component[u] = numberOfComponents;
                ++numberOfComponents;
            }
        }

        if (numberOfComponents == 0) {
            break;
        }

        for (index u = 0; u < numberOfNodes; ++u) {
            if (component[u] == none) {
                component[u] = numberOfComponents++;
            }
        }

        // Contract the cycles. Entering a node costs the difference to its selected edge.
        std::vector<Arc> contracted;
        for (const Arc &arc : arcs) {
            if (component[arc.from] != component[arc.to]) {
                const Arc &selected = arcs[cheapest[arc.to]];
                const index id = inputArcs.size() + replacements.size();
                replacements.push_back({arc.id, selected.id});
                selections.push_back(0);
                contracted.push_back(
                    {component[arc.from], component[arc.to], arc.weight - selected.weight, id});
            }
        }

        arcs = std::move(contracted);
        contractedRoot = component[contractedRoot];
        numberOfNodes = numberOfComponents;
    }

    // Replay the contractions backwards, so that only edges of the input graph remain selected.
    for (index i = replacements.size(); i > 0; --i) {
        const index id = inputArcs.size() + i - 1;
        const Replacement &replacement = replacements[i - 1];
        selections[replacement.added] += selections[id];
        selections[replacement.replaced] -= selections[id];
    }

    for (const Arc &arc : inputArcs) {
        if (selections[arc.id] > 0) {
            arborescence.addEdge(arc.from, arc.to, arc.weight);
            totalWeight += arc.weight;
        }
    }

    exists = true;
    hasRun = true;
}

} /* namespace NetworKit */
