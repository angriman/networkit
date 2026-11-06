/*
 * RBGraphReader.hpp
 *
 *  Created on: 16.10.2024
 *      Author: bernlu
 */

#ifndef NETWORKIT_IO_RB_GRAPH_READER_HPP_
#define NETWORKIT_IO_RB_GRAPH_READER_HPP_

#include <networkit/io/GraphReader.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Reader for the Rutherford Boeing (RB) matrix file format as described in
 * http://sparse-files.engr.tamu.edu/files/DOC/rb.pdf.
 *
 * @note currently the reader only supports compressed column format for real, integer, or pattern
 * data types.
 *
 * Reading fails if the number of nodes does not fit into GraphT::NodeT, or if GraphT::EdgeWeightT
 * is integral and an edge weight in the file is not an integer representable by it.
 */
template <class GraphT>
class GenericRBGraphReader final : public GraphReaderBase<GraphT> {
public:
    GenericRBGraphReader() = default;

    /**
     * Takes a file path as parameter and returns a graph.
     *
     * @param path
     * @return Graph
     */
    GraphT read(std::string_view path) override;
};

using RBGraphReader = GenericRBGraphReader<Graph>;

} /* namespace NetworKit */

#include <networkit/io/RBGraphReaderImpl.hpp>

#endif // NETWORKIT_IO_RB_GRAPH_READER_HPP_
