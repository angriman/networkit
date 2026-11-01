// networkit-format
#include <algorithm>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/matching/BlossomMatcher.hpp>

namespace NetworKit {

namespace {

struct Arc {
    index to;
    edgeweight weight;
};

struct VertexPair {
    index first = none;
    index second = none;
};

// A sub-blossom of a blossom, and the edge from it to the next one along the odd cycle.
struct CycleEntry {
    index blossom;
    VertexPair toNext;
};

enum class Label : uint8_t { Free, Outer, Inner };

// Vertices are the trivial blossoms, with ids below n. In each stage, alternating trees grow from
// the blossoms with a free base: outer blossoms are at even distance from a root, inner ones at
// odd distance. A blossom keeps its base as the first entry of its cycle, and the other entries
// are matched in pairs: the second with the third, and so on.
class Solver {
public:
    Solver(std::vector<std::vector<Arc>> adjacency, edgeweight maxWeight)
        : n(adjacency.size()), adjacency(std::move(adjacency)), mate(n, none),
          dual(n, maxWeight / 2), tolerance(maxWeight * 1e-11), top(n), parent(n + n / 2 + 2, none),
          base(parent.size(), none), blossomDual(parent.size(), 0),
          label(parent.size(), Label::Free), labelEdge(parent.size()), cycle(parent.size()),
          visited(parent.size(), 0) {
        std::iota(top.begin(), top.end(), index{0});
        std::iota(base.begin(), base.begin() + static_cast<std::ptrdiff_t>(n), index{0});
        for (index id = parent.size(); id > n; --id) {
            unusedIds.push_back(id - 1);
        }
    }

    const std::vector<index> &getMates() const { return mate; }

    // Returns false once the matching is optimal.
    bool runStage() {
        std::fill(label.begin(), label.end(), Label::Free);
        std::fill(labelEdge.begin(), labelEdge.end(), VertexPair{});
        outerToFree = {};
        outerToOuter = {};
        queue.clear();
        elapsed = 0;

        for (index v = 0; v < n; ++v) {
            if (mate[v] == none && base[top[v]] == v) {
                labelOuter(top[v], {});
            }
        }

        if (queue.empty()) {
            return false;
        }

        // Dual adjustments may leave the queue empty, as expanding a blossom adds no outer vertices
        // when it is entered at its base. Only they can tell that the matching is optimal.
        while (true) {
            while (!queue.empty()) {
                const index u = queue.back();
                queue.pop_back();
                if (scan(u)) {
                    return true;
                }
            }

            const auto action = adjustDuals();
            if (action.kind == ActionKind::Stop) {
                return false;
            }
            if (action.kind == ActionKind::Expand) {
                expandInner(action.blossom);
            } else if (handleTightEdge(action.edge.first, action.edge.second)) {
                return true;
            }
        }
    }

private:
    enum class ActionKind : uint8_t { Stop, TightEdge, Expand };

    struct Action {
        ActionKind kind;
        VertexPair edge;
        index blossom;
    };

    // Slacks are stored as they will be when the time equals the key, as they shrink over time.
    using Candidate = std::tuple<double, index, index>;
    using CandidateQueue =
        std::priority_queue<Candidate, std::vector<Candidate>, std::greater<Candidate>>;

    template <class Callback>
    void forLeaves(index blossom, Callback &&callback) const {
        std::vector<index> stack{blossom};
        while (!stack.empty()) {
            const index current = stack.back();
            stack.pop_back();
            if (current < n) {
                callback(current);
            } else {
                for (const CycleEntry &entry : cycle[current]) {
                    stack.push_back(entry.blossom);
                }
            }
        }
    }

    double slack(index u, index v, edgeweight weight) const { return dual[u] + dual[v] - weight; }

    double slackOf(index u, index v) const {
        const auto &arcs = adjacency[u];
        const auto arc = std::lower_bound(
            arcs.begin(), arcs.end(), v,
            [](const Arc &candidate, index target) -> bool { return candidate.to < target; });
        return slack(u, v, arc->weight);
    }

    void labelOuter(index blossom, VertexPair edge) {
        label[blossom] = Label::Outer;
        labelEdge[blossom] = edge;
        forLeaves(blossom, [&](index v) -> void { queue.push_back(v); });
    }

    // Returns true if the matching was augmented.
    bool scan(index u) {
        for (const Arc &arc : adjacency[u]) {
            const index v = arc.to;
            if (top[u] == top[v] || label[top[v]] == Label::Inner) {
                continue;
            }

            const double edgeSlack = slack(u, v, arc.weight);
            if (edgeSlack <= tolerance) {
                if (handleTightEdge(u, v)) {
                    return true;
                }
            } else if (label[top[v]] == Label::Free) {
                outerToFree.emplace(edgeSlack + elapsed, u, v);
            } else {
                outerToOuter.emplace(edgeSlack / 2 + elapsed, u, v);
            }
        }
        return false;
    }

    // Handles a tight edge from an outer vertex. Returns true if the matching was augmented.
    bool handleTightEdge(index u, index v) {
        const index blossomOfV = top[v];
        if (label[blossomOfV] == Label::Free) {
            label[blossomOfV] = Label::Inner;
            labelEdge[blossomOfV] = {u, v};
            const index partner = mate[base[blossomOfV]];
            labelOuter(top[partner], {base[blossomOfV], partner});
            return false;
        }

        if (label[blossomOfV] == Label::Outer && top[u] != blossomOfV) {
            const index ancestor = commonAncestor(top[u], blossomOfV);
            if (ancestor == none) {
                augmentFrom(u, v);
                augmentFrom(v, u);
                return true;
            }
            addBlossom(ancestor, u, v);
        }
        return false;
    }

    index outerParent(index outer) const {
        const VertexPair edge = labelEdge[outer];
        return edge.first == none ? none : top[labelEdge[top[edge.first]].first];
    }

    index commonAncestor(index a, index b) {
        ++currentVisit;
        while (a != none || b != none) {
            if (a != none) {
                if (visited[a] == currentVisit) {
                    return a;
                }
                visited[a] = currentVisit;
                a = outerParent(a);
            }
            std::swap(a, b);
        }
        return none;
    }

    // Contracts the odd cycle closed by the tight edge (u, v) between two outer blossoms.
    void addBlossom(index ancestor, index u, index v) {
        const index blossom = unusedIds.back();
        unusedIds.pop_back();

        // From the ancestor down the tree to u, across to v, and back up.
        std::vector<CycleEntry> entries;
        for (index current = top[u]; current != ancestor;) {
            const VertexPair edge = labelEdge[current];
            current = top[edge.first];
            entries.push_back({current, edge});
        }
        std::reverse(entries.begin(), entries.end());
        entries.push_back({top[u], {u, v}});
        for (index current = top[v]; current != ancestor;) {
            const VertexPair edge = labelEdge[current];
            entries.push_back({current, {edge.second, edge.first}});
            current = top[edge.first];
        }

        base[blossom] = base[ancestor];
        blossomDual[blossom] = 0;
        label[blossom] = Label::Outer;
        labelEdge[blossom] = labelEdge[ancestor];
        cycle[blossom] = std::move(entries);

        // Inner vertices have not been scanned yet, unlike outer ones.
        for (const CycleEntry &entry : cycle[blossom]) {
            parent[entry.blossom] = blossom;
            const bool scanLeaves = label[entry.blossom] == Label::Inner;
            forLeaves(entry.blossom, [&](index leaf) -> void {
                top[leaf] = blossom;
                if (scanLeaves) {
                    queue.push_back(leaf);
                }
            });
        }
    }

    static index positionInCycle(const std::vector<CycleEntry> &entries, index child) {
        return static_cast<index>(
            std::find_if(entries.begin(), entries.end(),
                         [&](const CycleEntry &entry) -> bool { return entry.blossom == child; })
            - entries.begin());
    }

    index childContaining(index blossom, index leaf) const {
        index child = leaf;
        while (parent[child] != blossom) {
            child = parent[child];
        }
        return child;
    }

    // Makes a vertex the base of a blossom by rematching along the even path around its cycle.
    void makeBase(index outermost, index vertex) {
        std::vector<std::pair<index, index>> stack{{outermost, vertex}};
        while (!stack.empty()) {
            const auto [blossom, newBase] = stack.back();
            stack.pop_back();
            if (blossom < n) {
                continue;
            }

            auto &entries = cycle[blossom];
            const index child = childContaining(blossom, newBase);
            const index position = positionInCycle(entries, child);
            const index size = entries.size();
            stack.emplace_back(child, newBase);

            // The new base leaves its pair, which shifts the pairs between it and the old base.
            const index first = position % 2 == 0 ? 0 : position + 1;
            const index last = position % 2 == 0 ? position : size;
            for (index i = first; i < last; i += 2) {
                const VertexPair edge = entries[i].toNext;
                mate[edge.first] = edge.second;
                mate[edge.second] = edge.first;
                stack.emplace_back(entries[i].blossom, edge.first);
                stack.emplace_back(entries[(i + 1) % size].blossom, edge.second);
            }

            std::rotate(entries.begin(), entries.begin() + static_cast<std::ptrdiff_t>(position),
                        entries.end());
            base[blossom] = newBase;
        }
    }

    // Flips the path from an outer vertex to the root of its tree.
    void augmentFrom(index vertex, index partner) {
        while (true) {
            const index outer = top[vertex];
            makeBase(outer, vertex);
            mate[vertex] = partner;

            const VertexPair toInner = labelEdge[outer];
            if (toInner.first == none) {
                return;
            }

            const index inner = top[toInner.first];
            const VertexPair entry = labelEdge[inner];
            makeBase(inner, entry.second);
            mate[entry.second] = entry.first;
            vertex = entry.first;
            partner = entry.second;
        }
    }

    // Replaces an inner blossom without dual by the even path through its cycle from the vertex
    // that was reached to the base. The other sub-blossoms leave the tree.
    void expandInner(index blossom) {
        const VertexPair entryEdge = labelEdge[blossom];
        std::vector<CycleEntry> entries = std::move(cycle[blossom]);
        cycle[blossom].clear();
        const index position = positionInCycle(entries, childContaining(blossom, entryEdge.second));

        for (const CycleEntry &entry : entries) {
            const index child = entry.blossom;
            parent[child] = none;
            label[child] = Label::Free;
            labelEdge[child] = {};
            forLeaves(child, [&](index leaf) -> void { top[leaf] = child; });
        }

        const index size = entries.size();
        const bool forward = position % 2 == 1;
        const auto next = [&](index i) -> index {
            return forward ? (i + 1) % size : (i + size - 1) % size;
        };
        const auto edgeBetween = [&](index from, index to) -> VertexPair {
            if (forward) {
                return entries[from].toNext;
            }
            return {entries[to].toNext.second, entries[to].toNext.first};
        };

        index current = position;
        label[entries[current].blossom] = Label::Inner;
        labelEdge[entries[current].blossom] = entryEdge;
        while (current != 0) {
            const index outer = next(current);
            const index inner = next(outer);
            labelOuter(entries[outer].blossom, edgeBetween(current, outer));
            label[entries[inner].blossom] = Label::Inner;
            labelEdge[entries[inner].blossom] = edgeBetween(outer, inner);
            current = inner;
        }

        // Slacks towards these leaves did not change while they were inner, unlike from now on.
        for (const CycleEntry &entry : entries) {
            if (label[entry.blossom] == Label::Free) {
                forLeaves(entry.blossom, [&](index leaf) -> void {
                    for (const Arc &arc : adjacency[leaf]) {
                        if (label[top[arc.to]] == Label::Outer) {
                            outerToFree.emplace(slack(arc.to, leaf, arc.weight) + elapsed, arc.to,
                                                leaf);
                        }
                    }
                });
            }
        }

        blossomDual[blossom] = 0;
        label[blossom] = Label::Free;
        labelEdge[blossom] = {};
        unusedIds.push_back(blossom);
    }

    // Drops candidates that are stale. Returns false if none is left.
    bool validateTop(CandidateQueue &candidates, bool towardsOuter) {
        const double divisor = towardsOuter ? 2 : 1;
        while (!candidates.empty()) {
            const auto [key, u, v] = candidates.top();
            const bool valid = label[top[u]] == Label::Outer
                               && (towardsOuter ? label[top[v]] == Label::Outer && top[u] != top[v]
                                                : label[top[v]] == Label::Free);
            if (!valid) {
                candidates.pop();
                continue;
            }

            // A vertex that was inner for a while has a higher slack than when it was queued.
            const double currentKey = slackOf(u, v) / divisor + elapsed;
            if (currentKey <= key + tolerance) {
                return true;
            }
            candidates.pop();
            candidates.emplace(currentKey, u, v);
        }
        return false;
    }

    Action adjustDuals() {
        Action action{ActionKind::Stop, {}, none};
        double delta = std::numeric_limits<double>::infinity();
        for (index v = 0; v < n; ++v) {
            if (label[top[v]] == Label::Outer) {
                delta = std::min(delta, dual[v]);
            }
        }

        if (validateTop(outerToFree, false)) {
            const auto [key, u, v] = outerToFree.top();
            const double edgeDelta = std::max(0., slackOf(u, v));
            if (edgeDelta < delta) {
                delta = edgeDelta;
                action = {ActionKind::TightEdge, {u, v}, none};
            }
        }
        if (validateTop(outerToOuter, true)) {
            const auto [key, u, v] = outerToOuter.top();
            const double edgeDelta = std::max(0., slackOf(u, v) / 2);
            if (edgeDelta < delta) {
                delta = edgeDelta;
                action = {ActionKind::TightEdge, {u, v}, none};
            }
        }
        for (index v = 0; v < n; ++v) {
            const index blossom = top[v];
            if (blossom >= n && base[blossom] == v && label[blossom] == Label::Inner
                && blossomDual[blossom] / 2 < delta) {
                delta = blossomDual[blossom] / 2;
                action = {ActionKind::Expand, {}, blossom};
            }
        }

        for (index v = 0; v < n; ++v) {
            const index blossom = top[v];
            if (label[blossom] == Label::Outer) {
                dual[v] -= delta;
            } else if (label[blossom] == Label::Inner) {
                dual[v] += delta;
            }
            if (blossom >= n && base[blossom] == v) {
                if (label[blossom] == Label::Outer) {
                    blossomDual[blossom] += 2 * delta;
                } else if (label[blossom] == Label::Inner) {
                    blossomDual[blossom] -= 2 * delta;
                }
            }
        }
        elapsed += delta;

        if (action.kind == ActionKind::Expand) {
            blossomDual[action.blossom] = 0;
        }
        return action;
    }

    const count n;
    const std::vector<std::vector<Arc>> adjacency;
    std::vector<index> mate;
    std::vector<double> dual;
    const double tolerance;

    // Outermost blossom of each vertex.
    std::vector<index> top;

    // The following are indexed by blossom.
    std::vector<index> parent;
    std::vector<index> base;
    std::vector<double> blossomDual;
    std::vector<Label> label;
    // The edge that added an outermost blossom to its tree, ending in the blossom.
    std::vector<VertexPair> labelEdge;
    std::vector<std::vector<CycleEntry>> cycle;
    std::vector<index> unusedIds;

    std::vector<index> visited;
    index currentVisit = 0;

    std::vector<index> queue;
    CandidateQueue outerToFree;
    CandidateQueue outerToOuter;
    double elapsed = 0;
};

} // namespace

BlossomMatcher::BlossomMatcher(const Graph &G) : Matcher(G) {
    if (G.isDirected()) {
        throw std::invalid_argument("BlossomMatcher: the graph must be undirected.");
    }
}

void BlossomMatcher::run() {
    Aux::SignalHandler handler;
    const count z = G->upperNodeIdBound();
    M = Matching(z);

    // Node ids may have gaps, and blossom ids follow the vertex ids.
    std::vector<node> nodes;
    nodes.reserve(G->numberOfNodes());
    std::vector<index> positionOf(z, none);
    G->forNodes([&](node u) -> void {
        positionOf[u] = nodes.size();
        nodes.push_back(u);
    });

    std::vector<std::vector<Arc>> adjacency(nodes.size());
    edgeweight maxWeight = 0;
    G->forEdges([&](node u, node v, edgeweight weight) -> void {
        if (u != v && weight > 0) {
            adjacency[positionOf[u]].push_back({positionOf[v], weight});
            adjacency[positionOf[v]].push_back({positionOf[u], weight});
            maxWeight = std::max(maxWeight, weight);
        }
    });

    for (auto &arcs : adjacency) {
        std::sort(arcs.begin(), arcs.end(), [](const Arc &a, const Arc &b) -> bool {
            return std::tie(a.to, b.weight) < std::tie(b.to, a.weight);
        });
        arcs.erase(std::unique(arcs.begin(), arcs.end(),
                               [](const Arc &a, const Arc &b) -> bool { return a.to == b.to; }),
                   arcs.end());
    }

    Solver solver(std::move(adjacency), maxWeight);
    while (solver.runStage()) {
        handler.assureRunning();
    }

    const auto &mates = solver.getMates();
    for (index i = 0; i < mates.size(); ++i) {
        if (mates[i] != none && i < mates[i]) {
            M.match(nodes[i], nodes[mates[i]]);
        }
    }

    hasRun = true;
}

} /* namespace NetworKit */
