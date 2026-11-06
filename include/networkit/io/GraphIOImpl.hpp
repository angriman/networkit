/*
 * GraphIOImpl.hpp
 *
 *  Created on: 09.01.2013
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_IO_GRAPH_IO_IMPL_HPP_
#define NETWORKIT_IO_GRAPH_IO_IMPL_HPP_

#include <fstream>

#include <networkit/io/GraphIO.hpp>

namespace NetworKit {

template <class GraphT>
void GraphIO::writeEdgeList(const GraphT &G, std::string_view path) {
    using NodeT = typename GraphT::NodeT;

    std::ofstream file;
    file.open(path.data());

    // The casts print the ids as numbers also for 8-bit node types.
    G.forEdges([&](const NodeT v, const NodeT w) -> void {
        file << static_cast<count>(v) << " " << static_cast<count>(w) << '\n';
    });

    file.close();
}

template <class GraphT>
void GraphIO::writeAdjacencyList(const GraphT &G, std::string_view path) {
    using NodeT = typename GraphT::NodeT;

    std::ofstream file;
    file.open(path.data());

    G.forNodes([&](const NodeT v) -> void {
        file << static_cast<count>(v);
        G.forNeighborsOf(v, [&](const NodeT x) -> void { file << " " << static_cast<count>(x); });
        file << '\n';
    });
}

} /* namespace NetworKit */

#endif // NETWORKIT_IO_GRAPH_IO_IMPL_HPP_
