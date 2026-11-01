// networkit-format
#include <algorithm>
#include <vector>

#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/centrality/VoteRank.hpp>

namespace NetworKit {

VoteRank::VoteRank(const Graph &G, count k) : G(&G), k(k) {}

void VoteRank::run() {
    Aux::SignalHandler handler;
    const count z = G->upperNodeIdBound();

    spreaders.clear();
    if (G->numberOfEdges() == 0) {
        hasRun = true;
        return;
    }

    const double averageDegree = static_cast<double>(G->numberOfEdges() * (G->isDirected() ? 1 : 2))
                                 / static_cast<double>(G->numberOfNodes());
    const double abilityLoss = 1. / averageDegree;

    std::vector<double> votingAbility(z, 1.), votes(z);
    std::vector<unsigned char> elected(z, 0);

    while (spreaders.size() < k) {
        handler.assureRunning();

        G->balancedParallelForNodes([&](node u) -> void {
            votes[u] = 0.;
            if (!elected[u]) {
                G->forNeighborsOf(u, [&](node v) -> void { votes[u] += votingAbility[v]; });
            }
        });

        node winner = none;
        double mostVotes = 0.;
        G->forNodes([&](node u) -> void {
            if (votes[u] > mostVotes) {
                mostVotes = votes[u];
                winner = u;
            }
        });

        if (winner == none) {
            break;
        }

        spreaders.push_back(winner);
        elected[winner] = 1;
        votingAbility[winner] = 0.;
        G->forNeighborsOf(winner, [&](node v) -> void {
            votingAbility[v] = std::max(votingAbility[v] - abilityLoss, 0.);
        });
    }

    hasRun = true;
}

} /* namespace NetworKit */
