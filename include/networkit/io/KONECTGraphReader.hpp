/*
 * KONECTGraphReader.hpp
 *
 *  Created on: 11.05.2018
 *      Author: Roman Bange
 */

#ifndef NETWORKIT_IO_KONECT_GRAPH_READER_HPP_
#define NETWORKIT_IO_KONECT_GRAPH_READER_HPP_

#include <unordered_map>

#include <networkit/graph/Graph.hpp>
#include <networkit/io/GraphReader.hpp>

namespace NetworKit {
/**
 * @ingroup io
 *
 * Reading fails if the number of nodes or, without remapping, a node id does not fit into
 * GraphT::NodeT, or if GraphT::EdgeWeightT is integral and an edge weight in the file is not an
 * integer representable by it.
 */
template <class GraphT>
class GenericKONECTGraphReader final : public GraphReaderBase<GraphT> {
public:
    /*
     * If the input graph has multiple edges, you can specify on how these edges are handled.
     * Keep in mind that NetworKit node id's start with 0 while most KONECT graphs start with 1.
     * See GraphReader.h for a closer description of the parameters.
     *
     * @param[in]  remapNodes  specifies whether node ids should be remapped if non consecutive
     * @param[in]  handlingmethod  specifies how multiple edges should be handled (only relevant if
     * graph with multiple edges is given)
     */
    GenericKONECTGraphReader(bool remapNodes = false,
                             GraphReaderOptions::MultipleEdgesHandling handlingmethod =
                                 GraphReaderOptions::DISCARD_EDGES);

    /**
     * Given the path of an input file, read the graph contained.
     *
     * @param[in]  path  input file path
     */
    GraphT read(std::string_view path) override;

private:
    bool remapNodes;
    GraphReaderOptions::MultipleEdgesHandling multipleEdgesHandlingMethod;
};

using KONECTGraphReader = GenericKONECTGraphReader<Graph>;

} /* namespace NetworKit */

#include <networkit/io/KONECTGraphReaderImpl.hpp>

#endif // NETWORKIT_IO_KONECT_GRAPH_READER_HPP_
