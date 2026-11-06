/*
 * ThrillGraphBinaryWriter.hpp
 *
 * @author Michael Hamann
 */

#ifndef NETWORKIT_IO_THRILL_GRAPH_BINARY_WRITER_HPP_
#define NETWORKIT_IO_THRILL_GRAPH_BINARY_WRITER_HPP_

#include <networkit/io/GraphWriter.hpp>

namespace NetworKit {

/**
 * Writes a graph in the format read by GenericThrillGraphBinaryReader. Writing fails if a node id
 * does not fit into 32 bits.
 */
template <class GraphT>
class GenericThrillGraphBinaryWriter final : public GraphWriterBase<GraphT> {
public:
    /**
     * Write the given graph into a binary file at the given path.
     *
     * @param[in] G The graph to write.
     * @param[in] path The path where to write the graph.
     */
    void write(const GraphT &G, std::string_view path) override;
};

using ThrillGraphBinaryWriter = GenericThrillGraphBinaryWriter<Graph>;

} /* namespace NetworKit */

#include <networkit/io/ThrillGraphBinaryWriterImpl.hpp>

#endif // NETWORKIT_IO_THRILL_GRAPH_BINARY_WRITER_HPP_
