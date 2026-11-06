#ifndef NETWORKIT_IO_MTX_GRAPH_READER_IMPL_HPP_
#define NETWORKIT_IO_MTX_GRAPH_READER_IMPL_HPP_

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>

#include <networkit/io/MTXGraphReader.hpp>
#include <networkit/io/MTXParser.hpp>

namespace NetworKit {

template <class GraphT>
GraphT GenericMTXGraphReader<GraphT>::read(std::string_view path) {
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;

    MTXParser parser{std::string(path)};

    auto header = parser.getHeader();
    auto size = parser.getMatrixSize();
    auto weighted = true;
    auto symmetric = true;

    if (header.field == MTXParser::Field::Pattern)
        weighted = false;
    if (header.symmetry == MTXParser::Symmetry::General)
        symmetric = false;

    const count n = std::max(size.rows, size.columns);
    // The node ids 0, ..., n - 1 must be representable by NodeT and differ from NullNodeId<NodeT>,
    // which marks "no node". Since NullNodeId<NodeT> is the largest NodeT, this allows at most
    // max() nodes.
    static_assert(NullNodeId<NodeT> == std::numeric_limits<NodeT>::max());
    constexpr auto maxNumberOfNodes = static_cast<count>(std::numeric_limits<NodeT>::max());
    if (n > maxNumberOfNodes) {
        throw std::runtime_error("MTX file " + std::string(path) + " has " + std::to_string(n)
                                 + " nodes, but the node type of the graph can hold at most "
                                 + std::to_string(maxNumberOfNodes) + " nodes");
    }

    GraphT G(n, weighted, !symmetric);

    for (auto edge = parser.getNext<EdgeWeightT>(weighted); edge.has_value();
         edge = parser.getNext<EdgeWeightT>(weighted)) {
        if (edge->u >= n || edge->v >= n) {
            throw std::runtime_error("MTX file " + std::string(path) + " has an entry in row "
                                     + std::to_string(edge->u + 1) + " and column "
                                     + std::to_string(edge->v + 1) + ", but the matrix has only "
                                     + std::to_string(n) + " rows and columns");
        }
        G.addEdge(static_cast<NodeT>(edge->u), static_cast<NodeT>(edge->v), edge->weight);
    }
    return G;
}

} /* namespace NetworKit */

#endif // NETWORKIT_IO_MTX_GRAPH_READER_IMPL_HPP_
