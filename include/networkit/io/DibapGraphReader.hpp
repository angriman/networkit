/*
 * DibapGraphReader.hpp
 *
 *  Created on: Jun 12, 2013
 *      Author: Henning
 */

#ifndef NETWORKIT_IO_DIBAP_GRAPH_READER_HPP_
#define NETWORKIT_IO_DIBAP_GRAPH_READER_HPP_

#ifndef NETWORKIT_WINDOWS

#include <vector>

#include <networkit/graph/Graph.hpp>
#include <networkit/io/GraphReader.hpp>
#include <networkit/viz/Point.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * TODO: class documentation
 *
 * Reading fails if the number of nodes does not fit into GraphT::NodeT.
 */
template <class GraphT>
class GenericDibapGraphReader final : public GraphReaderBase<GraphT> {
public:
    GenericDibapGraphReader() = default;

    GraphT read(std::string_view path) override;

    const std::vector<Point<coordinate>> &getCoordinates() const noexcept { return coordinates; }

    std::vector<Point<coordinate>> moveCoordinates() { return std::move(coordinates); }

private:
    std::vector<Point<coordinate>> coordinates;
};

using DibapGraphReader = GenericDibapGraphReader<Graph>;

} /* namespace NetworKit */

#include <networkit/io/DibapGraphReaderImpl.hpp>

#endif // NETWORKIT_WINDOWS

#endif // NETWORKIT_IO_DIBAP_GRAPH_READER_HPP_
