/*
 * DotGraphWriterImpl.hpp
 *
 *  Created on: Jun 5, 2013
 *      Author: forigem
 */

#ifndef NETWORKIT_IO_DOT_GRAPH_WRITER_IMPL_HPP_
#define NETWORKIT_IO_DOT_GRAPH_WRITER_IMPL_HPP_

#include <fstream>

#include <networkit/io/DotGraphWriter.hpp>

namespace NetworKit {

template <class GraphT>
void GenericDotGraphWriter<GraphT>::write(const GraphT &G, std::string_view path) {
    using NodeT = typename GraphT::NodeT;

    std::ofstream file{path.data()};

    file << "graph {\n";
    // The casts print the ids as numbers also for 8-bit node types.
    G.forEdges([&](NodeT u, NodeT v) -> void {
        file << static_cast<count>(u) << " -- " << static_cast<count>(v) << ";\n";
    });
    file << "}\n";
}

} /* namespace NetworKit */

#endif // NETWORKIT_IO_DOT_GRAPH_WRITER_IMPL_HPP_
