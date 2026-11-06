/*
 * ThrillGraphBinaryReaderImpl.hpp
 *
 * @author Michael Hamann
 */

#ifndef NETWORKIT_IO_THRILL_GRAPH_BINARY_READER_IMPL_HPP_
#define NETWORKIT_IO_THRILL_GRAPH_BINARY_READER_IMPL_HPP_

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <networkit/io/ThrillGraphBinaryReader.hpp>

namespace NetworKit {

namespace ThrillGraphBinaryReaderDetails {

inline uint32_t get_uint32(std::ifstream &is) {
    uint32_t result = 0;

    for (size_t i = 0; i < sizeof(uint32_t); ++i) {
        uint32_t u = is.get();
        result |= (u << (i * 8));
    }

    return result;
}

inline uint64_t get_variant(std::ifstream &is) {
    size_t v = 0;
    for (size_t shift_width = 0; shift_width < 64; shift_width += 7) {
        uint64_t u = is.get();

        // The last value is just a single bit (the 64th bit) - throw if there is
        // more than this single bit in the input.
        if (shift_width == 63 && (u & 0xFE)) {
            throw std::overflow_error("Overflow during variant64 decoding.");
        }

        // Read exactly 7 bits from the input
        v |= (u & 0x7F) << shift_width;

        // If the 8th bit is not set, we reached the end of the variant
        if (!(u & 0x80))
            break;
    }

    return v;
}

} // namespace ThrillGraphBinaryReaderDetails

template <class GraphT>
GenericThrillGraphBinaryReader<GraphT>::GenericThrillGraphBinaryReader(count n) : n(n) {}

template <class GraphT>
GraphT GenericThrillGraphBinaryReader<GraphT>::read(std::string_view path) {
    return read(std::vector<std::string>(1, std::string(path)));
}

template <class GraphT>
GraphT GenericThrillGraphBinaryReader<GraphT>::read(const std::vector<std::string> &paths) {
    using NodeT = typename GraphT::NodeT;
    using ThrillGraphBinaryReaderDetails::get_uint32;
    using ThrillGraphBinaryReaderDetails::get_variant;

    // NullNodeId<NodeT> is the largest NodeT and marks "no node", so that at most max() nodes fit.
    static_assert(NullNodeId<NodeT> == std::numeric_limits<NodeT>::max());
    constexpr auto maxNumberOfNodes = static_cast<count>(std::numeric_limits<NodeT>::max());
    if (n > maxNumberOfNodes) {
        throw std::runtime_error("The number of nodes " + std::to_string(n)
                                 + " does not fit into the node type of the graph");
    }

    GraphT G(n);
    // The edges are added at the end, since they may point to nodes that have not been read yet.
    std::vector<std::pair<NodeT, uint32_t>> edges;

    if (!paths.empty()) {
        std::ifstream is;
        count file_index = 0;

        auto next_input = [&]() -> void {
            is.close();
            if (file_index < paths.size()) {
                is.open(paths[file_index++], std::ios::binary);
                // Throw an exception when no byte could be read at any point
                is.exceptions(std::ifstream::failbit | std::ifstream::badbit);
            }
        };

        next_input();

        uint32_t max_id = 0;

        for (count u = 0; is.good() && is.is_open(); ++u) {
            // Add node if it does not exist yet, only one may be missing
            if (u >= static_cast<count>(G.upperNodeIdBound())) {
                if (G.upperNodeIdBound() == NullNodeId<NodeT>) {
                    throw std::runtime_error("The file has more nodes than the node type of the "
                                             "graph can hold");
                }
                G.addNode();
            }

            for (count deg = get_variant(is); deg > 0 && is.good(); --deg) {
                uint32_t v = get_uint32(is);

                max_id = std::max(max_id, v);

                edges.emplace_back(static_cast<NodeT>(u), v);
            }

            if (is.is_open() && (is.peek() == std::char_traits<char>::eof() || !is.good())) {
                next_input();
            }
        }

        if (max_id >= static_cast<count>(G.upperNodeIdBound())) {
            throw std::runtime_error("Maximum read node id larger than number of nodes read.");
        }
    }

    // All node ids are smaller than the upper node id bound (checked above), so they fit NodeT.
    for (const auto &[u, v] : edges) {
        G.addEdge(u, static_cast<NodeT>(v));
    }
    return G;
}

} // namespace NetworKit

#endif // NETWORKIT_IO_THRILL_GRAPH_BINARY_READER_IMPL_HPP_
