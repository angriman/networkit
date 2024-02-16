#include <networkit/community/LocalCoverEvaluation.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Cover.hpp>

NetworKit::LocalCoverEvaluation::LocalCoverEvaluation(const NetworKit::Graph &G,
                                                      const NetworKit::Cover &C)
    : G(&G), C(&C) {}
