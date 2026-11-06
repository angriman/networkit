/*
 * SNAPGraphWriterImpl.hpp
 *
 *  Created on: 24.09.2013
 *      Author: cls
 */

#ifndef NETWORKIT_IO_SNAP_GRAPH_WRITER_IMPL_HPP_
#define NETWORKIT_IO_SNAP_GRAPH_WRITER_IMPL_HPP_

#include <fstream>

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/io/SNAPGraphWriter.hpp>

namespace NetworKit {

template <class GraphT>
void GenericSNAPGraphWriter<GraphT>::write(const GraphT &G, std::string_view path) {
    using NodeT = typename GraphT::NodeT;

    std::ofstream file(path.data());
    Aux::enforceOpened(file);

    // write "problem line" - n, m, directed/undirected, weighted/weight type
    file << "p " << G.numberOfNodes() << " " << G.numberOfEdges() << " u u 0\n"; // FIXME: makeshift

    // The casts print the ids as numbers also for 8-bit node types.
    G.forEdges([&](NodeT u, NodeT v) -> void {
        file << static_cast<count>(u) << " " << static_cast<count>(v) << '\n';
    });

    file.close();
}

} /* namespace NetworKit */

#endif // NETWORKIT_IO_SNAP_GRAPH_WRITER_IMPL_HPP_
