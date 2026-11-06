/*
 * GMLGraphWriterImpl.hpp
 *
 *  Created on: 05.08.2013
 *      Author: Stefan Bertsch
 */

#ifndef NETWORKIT_IO_GML_GRAPH_WRITER_IMPL_HPP_
#define NETWORKIT_IO_GML_GRAPH_WRITER_IMPL_HPP_

#include <fstream>

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/io/GMLGraphWriter.hpp>

namespace NetworKit {

template <class GraphT>
void GenericGMLGraphWriter<GraphT>::write(const GraphT &G, std::string_view path) {
    using NodeT = typename GraphT::NodeT;

    std::ofstream file(path.data());
    Aux::enforceOpened(file);

    file << "graph [\n";
    if (G.isDirected()) {
        file << "  directed 1\n";
    }

    // The casts print the ids as numbers also for 8-bit node types.
    G.forNodes([&](NodeT u) -> void {
        file << "  node [\n";
        file << "    id " << static_cast<count>(u) << "\n";
        file << "  ]\n";
    });

    G.forEdges([&](NodeT u, NodeT v) -> void {
        file << "  edge [\n";
        file << "    source " << static_cast<count>(u) << "\n";
        file << "    target " << static_cast<count>(v) << "\n";
        file << "  ]\n";
    });
    file << "]\n";
}

} /* namespace NetworKit */

#endif // NETWORKIT_IO_GML_GRAPH_WRITER_IMPL_HPP_
