#ifndef NETWORKIT_IO_MTX_GRAPH_READER_HPP_
#define NETWORKIT_IO_MTX_GRAPH_READER_HPP_

#include <networkit/io/GraphReader.hpp>

namespace NetworKit {

/**
 * @ingroup io
 *
 * Reader for the matrix market file format documented in
 * https://networkrepository.com/mtx-matrix-market-format.html
 *
 * Does not allow complex fields.
 *
 * The file does not store the node and edge weight types: they are given by @a GraphT.
 * Reading fails if the number of nodes does not fit into GraphT::NodeT, or if GraphT::EdgeWeightT
 * is integral and an edge weight in the file is not an integer representable by it.
 */
template <class GraphT>
class GenericMTXGraphReader final : public GraphReaderBase<GraphT> {
public:
    GenericMTXGraphReader() = default;

    /**
     * Takes a file path as parameter and returns a graph.
     *
     * @param path
     * @return Graph
     */
    GraphT read(std::string_view path) override;
};

using MTXGraphReader = GenericMTXGraphReader<Graph>;

} /* namespace NetworKit */

#include <networkit/io/MTXGraphReaderImpl.hpp>

#endif // NETWORKIT_IO_MTX_GRAPH_READER_HPP_
