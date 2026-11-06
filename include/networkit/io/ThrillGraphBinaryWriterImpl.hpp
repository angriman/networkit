/*
 * ThrillGraphBinaryWriterImpl.hpp
 *
 * @author Michael Hamann
 */

#ifndef NETWORKIT_IO_THRILL_GRAPH_BINARY_WRITER_IMPL_HPP_
#define NETWORKIT_IO_THRILL_GRAPH_BINARY_WRITER_IMPL_HPP_

#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <vector>

#include <networkit/io/ThrillGraphBinaryWriter.hpp>

namespace NetworKit {

template <class GraphT>
void GenericThrillGraphBinaryWriter<GraphT>::write(const GraphT &G, std::string_view path) {
    using NodeT = typename GraphT::NodeT;

    const auto upperNodeIdBound = static_cast<count>(G.upperNodeIdBound());
    if (upperNodeIdBound > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error(
            "Thrill binary graphs only support graphs with up to 2^32-1 nodes.");
    }

    std::ofstream out_stream(path.data(), std::ios::trunc | std::ios::binary);

    std::vector<uint32_t> neighbors;

    for (count u = 0; u < upperNodeIdBound; ++u) {
        neighbors.clear();

        // u < upperNodeIdBound, so that it fits into NodeT.
        const auto uNode = static_cast<NodeT>(u);
        if (G.hasNode(uNode)) {
            G.forEdgesOf(uNode, [&](NodeT v) -> void {
                if (uNode <= v) {
                    neighbors.push_back(static_cast<uint32_t>(v));
                }
            });
        }

        size_t deg = neighbors.size();

        // Write variable length int for the degree
        if (!deg) {
            out_stream << uint8_t(0);
        }

        while (deg) {
            size_t u = deg & 0x7F;
            deg >>= 7;
            out_stream << uint8_t(u | (deg ? 0x80 : 0));
        }

        for (uint32_t v : neighbors) {
            // write neighbor as little endian
            for (size_t i = 0; i < sizeof(uint32_t); ++i) {
                out_stream << uint8_t(v);
                v >>= 8;
            }
        }
    }

    out_stream.close();
}

} // namespace NetworKit

#endif // NETWORKIT_IO_THRILL_GRAPH_BINARY_WRITER_IMPL_HPP_
