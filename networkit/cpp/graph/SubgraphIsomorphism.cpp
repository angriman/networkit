// networkit-format
#include <algorithm>
#include <stdexcept>
#include <utility>
#include <vector>

#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/graph/SubgraphIsomorphism.hpp>

namespace NetworKit {

namespace {

using NeighborLists = std::vector<std::vector<index>>;

void sortAndDeduplicate(std::vector<index> &values) {
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
}

// Returns the pattern positions in search order: each next node has as many neighbors among the
// previous ones as possible, so that wrong assignments are detected early. Ties are broken in
// favor of the node with more neighbors.
std::vector<index> searchOrder(const NeighborLists &outNeighbors,
                               const NeighborLists &inNeighbors) {
    const count n = outNeighbors.size();
    std::vector<index> order;
    std::vector<count> placedNeighbors(n, 0);
    std::vector<unsigned char> placed(n, 0);
    order.reserve(n);

    for (count step = 0; step < n; ++step) {
        index best = none;
        for (index p = 0; p < n; ++p) {
            if (placed[p]) {
                continue;
            }
            const count degree = outNeighbors[p].size() + inNeighbors[p].size();
            const count bestDegree =
                best == none ? 0 : outNeighbors[best].size() + inNeighbors[best].size();
            if (best == none || placedNeighbors[p] > placedNeighbors[best]
                || (placedNeighbors[p] == placedNeighbors[best] && degree > bestDegree)) {
                best = p;
            }
        }

        order.push_back(best);
        placed[best] = 1;
        for (const index q : outNeighbors[best]) {
            ++placedNeighbors[q];
        }
        for (const index q : inNeighbors[best]) {
            ++placedNeighbors[q];
        }
    }

    return order;
}

} // namespace

SubgraphIsomorphism::SubgraphIsomorphism(const Graph &G, const Graph &P, bool induced,
                                         count maxMatches)
    : G(&G), P(&P), induced(induced), maxMatches(maxMatches) {
    if (P.numberOfNodes() == 0) {
        throw std::invalid_argument("The pattern must have at least one node.");
    }
    if (maxMatches == 0) {
        throw std::invalid_argument("maxMatches must be positive.");
    }
    if (G.isDirected() != P.isDirected()) {
        throw std::invalid_argument(
            "The graph and the pattern must be both directed or undirected.");
    }
}

void SubgraphIsomorphism::run() {
    Aux::SignalHandler handler;
    matches.clear();

    // The pattern is stored with consecutive positions, without parallel edges and self-loops.
    const bool directed = G->isDirected();
    const count patternSize = P->numberOfNodes();
    std::vector<node> patternNodes;
    std::vector<index> positionOf(P->upperNodeIdBound(), none);
    patternNodes.reserve(patternSize);
    P->forNodes([&](node u) -> void {
        positionOf[u] = patternNodes.size();
        patternNodes.push_back(u);
    });

    NeighborLists patternOut(patternSize), patternIn(patternSize);
    std::vector<unsigned char> patternSelfLoop(patternSize, 0);
    P->forEdges([&](node u, node v) -> void {
        const index p = positionOf[u], q = positionOf[v];
        if (p == q) {
            patternSelfLoop[p] = 1;
        } else if (directed) {
            patternOut[p].push_back(q);
            patternIn[q].push_back(p);
        } else {
            patternOut[p].push_back(q);
            patternOut[q].push_back(p);
        }
    });
    for (index p = 0; p < patternSize; ++p) {
        sortAndDeduplicate(patternOut[p]);
        sortAndDeduplicate(patternIn[p]);
    }

    const std::vector<index> order = searchOrder(patternOut, patternIn);

    std::vector<node> assignment(patternSize, none);            // pattern position -> target node
    std::vector<index> assignedTo(G->upperNodeIdBound(), none); // target node -> pattern position
    std::vector<count> seenAt(patternSize, 0);
    count stamp = 0;

    // Checks the neighbors of a candidate target node in one direction against the pattern.
    auto neighborsMatch = [&](const auto &targetNeighbors, node candidate,
                              const std::vector<index> &neighbors, bool &targetSelfLoop) -> bool {
        ++stamp;
        count found = 0;
        for (const node t : targetNeighbors) {
            if (t == candidate) {
                targetSelfLoop = true;
                continue;
            }
            const index q = assignedTo[t];
            if (q == none) {
                continue;
            }
            if (std::binary_search(neighbors.begin(), neighbors.end(), q)) {
                if (seenAt[q] != stamp) {
                    seenAt[q] = stamp;
                    ++found;
                }
            } else if (induced) {
                return false;
            }
        }

        const count assignedNeighbors =
            std::count_if(neighbors.begin(), neighbors.end(),
                          [&](index q) -> bool { return assignment[q] != none; });
        return found == assignedNeighbors;
    };

    auto isFeasible = [&](index p, node candidate) -> bool {
        if (assignedTo[candidate] != none || G->degree(candidate) < patternOut[p].size()) {
            return false;
        }

        bool targetSelfLoop = false;
        if (!neighborsMatch(G->neighborRange(candidate), candidate, patternOut[p],
                            targetSelfLoop)) {
            return false;
        }
        if (directed
            && (G->degreeIn(candidate) < patternIn[p].size()
                || !neighborsMatch(G->inNeighborRange(candidate), candidate, patternIn[p],
                                   targetSelfLoop))) {
            return false;
        }

        return patternSelfLoop[p] ? targetSelfLoop : !(induced && targetSelfLoop);
    };

    // Candidates are the neighbors of an assigned pattern neighbor, or all nodes if there is none.
    auto fillCandidates = [&](index p, std::vector<node> &candidates) -> void {
        candidates.clear();

        node anchor = none;
        bool useInNeighbors = false;
        count fewest = none;
        for (const index q : patternOut[p]) {
            if (assignment[q] != none) {
                const bool reversed = directed;
                const count degree =
                    reversed ? G->degreeIn(assignment[q]) : G->degree(assignment[q]);
                if (degree < fewest) {
                    fewest = degree;
                    anchor = assignment[q];
                    useInNeighbors = reversed;
                }
            }
        }
        for (const index q : patternIn[p]) {
            if (assignment[q] != none && G->degree(assignment[q]) < fewest) {
                fewest = G->degree(assignment[q]);
                anchor = assignment[q];
                useInNeighbors = false;
            }
        }

        if (anchor == none) {
            G->forNodes([&](node u) -> void { candidates.push_back(u); });
        } else if (useInNeighbors) {
            G->forInNeighborsOf(anchor, [&](node, node u) -> void { candidates.push_back(u); });
        } else {
            G->forNeighborsOf(anchor, [&](node u) -> void { candidates.push_back(u); });
        }

        // Parallel edges must not lead to the same match twice.
        std::sort(candidates.begin(), candidates.end());
        candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
    };

    auto storeMatch = [&]() -> void {
        std::vector<node> match(P->upperNodeIdBound(), none);
        for (index p = 0; p < patternSize; ++p) {
            match[patternNodes[p]] = assignment[p];
        }
        matches.push_back(std::move(match));
    };

    // Iterative backtracking, the depth is the number of assigned pattern nodes.
    std::vector<std::vector<node>> candidates(patternSize);
    std::vector<index> nextCandidate(patternSize, 0);
    fillCandidates(order[0], candidates[0]);

    constexpr count stepsBetweenChecks = 1 << 16;
    count steps = 0;
    index depth = 0;

    while (matches.size() < maxMatches) {
        if (++steps % stepsBetweenChecks == 0) {
            handler.assureRunning();
        }

        if (nextCandidate[depth] == candidates[depth].size()) {
            if (depth == 0) {
                break;
            }
            --depth;
            const index p = order[depth];
            assignedTo[assignment[p]] = none;
            assignment[p] = none;
            continue;
        }

        const index p = order[depth];
        const node candidate = candidates[depth][nextCandidate[depth]++];
        if (!isFeasible(p, candidate)) {
            continue;
        }

        assignment[p] = candidate;
        assignedTo[candidate] = p;

        if (depth + 1 == patternSize) {
            storeMatch();
            assignedTo[candidate] = none;
            assignment[p] = none;
        } else {
            ++depth;
            fillCandidates(order[depth], candidates[depth]);
            nextCandidate[depth] = 0;
        }
    }

    hasRun = true;
}

} /* namespace NetworKit */
