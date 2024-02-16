#include <set>
#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/scd/CombinedSCD.hpp>
#include <networkit/scd/SelectiveCommunityDetector.hpp>

namespace NetworKit {

CombinedSCD::CombinedSCD(const Graph &g, SelectiveCommunityDetector &first,
                         SelectiveCommunityDetector &second)
    : SelectiveCommunityDetector(g), first(first), second(second) {}

std::set<node> CombinedSCD::expandOneCommunity(node s) {
    return second.expandOneCommunity(first.expandOneCommunity(s));
}

std::set<node> CombinedSCD::expandOneCommunity(const std::set<node> &s) {
    return second.expandOneCommunity(first.expandOneCommunity(s));
}

} // namespace NetworKit
