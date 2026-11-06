/*
 * NetworkitBinaryReaderImpl.hpp
 *
 * Author: Charmaine Ndolo <charmaine.ndolo@b-tu.de>
 *
 */

#ifndef NETWORKIT_IO_NETWORKIT_BINARY_READER_IMPL_HPP_
#define NETWORKIT_IO_NETWORKIT_BINARY_READER_IMPL_HPP_

#include <atomic>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/io/MemoryMappedFile.hpp>
#include <networkit/io/NetworkitBinaryGraph.hpp>
#include <networkit/io/NetworkitBinaryReader.hpp>

namespace NetworKit {

namespace NetworkitBinaryReaderDetails {

inline const char *accessData(const MemoryMappedFile &source) {
    return source.cbegin();
}

inline const char *accessData(const std::vector<uint8_t> &source) {
    return reinterpret_cast<const char *>(source.data());
}

/**
 * Converts the weight @a value of a file into @a weight. Returns false, leaving @a weight
 * unchanged, if EdgeWeightT is integral and @a value is not an integer representable by it.
 */
template <class EdgeWeightT, class T>
bool toEdgeWeight(T value, EdgeWeightT &weight) {
    if constexpr (std::is_integral_v<EdgeWeightT>) {
        using Limits = std::numeric_limits<EdgeWeightT>;
        if constexpr (std::is_floating_point_v<T>) {
            // max() + 1 is a power of two and thus exact, unlike max() for 64-bit types.
            constexpr double lowest = static_cast<double>(Limits::lowest());
            constexpr double maxPlusOne = 2.0 * static_cast<double>(Limits::max() / 2 + 1);
            const auto asDouble = static_cast<double>(value);
            if (std::trunc(asDouble) != asDouble || asDouble < lowest || asDouble >= maxPlusOne) {
                return false;
            }
        } else if constexpr (std::is_signed_v<T> && std::is_signed_v<EdgeWeightT>) {
            if (value < static_cast<T>(Limits::lowest()) || value > static_cast<T>(Limits::max())) {
                return false;
            }
        } else if constexpr (std::is_signed_v<T>) {
            // EdgeWeightT is unsigned and at most as wide as uint64_t, the type of the file.
            if (value < 0 || static_cast<uint64_t>(value) > static_cast<uint64_t>(Limits::max())) {
                return false;
            }
        } else if (static_cast<uint64_t>(value) > static_cast<uint64_t>(Limits::max())) {
            // T is unsigned, and Limits::max() is positive for any EdgeWeightT.
            return false;
        }
    }
    weight = static_cast<EdgeWeightT>(value);
    return true;
}

} // namespace NetworkitBinaryReaderDetails

template <class GraphT>
GraphT GenericNetworkitBinaryReader<GraphT>::read(std::string_view path) {
    MemoryMappedFile mmfile(path);
    GraphT G = readData(mmfile);
    return G;
}

template <class GraphT>
GraphT GenericNetworkitBinaryReader<GraphT>::readFromBuffer(const std::vector<uint8_t> &data) {
    return readData(data);
}

template <class GraphT>
template <class T>
GraphT GenericNetworkitBinaryReader<GraphT>::readData(const T &source) {
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;
    using NetworkitBinaryReaderDetails::accessData;
    using NetworkitBinaryReaderDetails::toEdgeWeight;

    nkbg::Header header;
    nkbg::WeightFormat weightFormat;

    const char *startIt = accessData(source);
    const char *it = startIt;
    auto readHeader = [&]() -> void {
        memcpy(&header.magic, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.checksum, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.features, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        if (!memcmp("nkbg002", header.magic, 8)) {
            version = 2;
        } else if (!memcmp("nkbg003", header.magic, 8)) {
            version = 3;
        } else if (!memcmp("nkbg004", header.magic, 8)) {
            version = 4;
        } else {
            throw std::runtime_error("Reader expected another magic value");
        }
        directed = (header.features & nkbg::DIR_MASK);
        weightFormat = static_cast<nkbg::WeightFormat>((header.features & nkbg::WGHT_MASK)
                                                       >> nkbg::WGHT_SHIFT);
        indexed = false;
        if (version >= 3) {
            indexed = (header.features & nkbg::INDEX_MASK) >> nkbg::INDEX_SHIFT;
        }
        memcpy(&header.nodes, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.chunks, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.offsetBaseData, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.offsetAdjLists, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.offsetAdjTranspose, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.offsetWeightLists, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.offsetWeightTranspose, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        if (version > 2) {
            memcpy(&header.offsetAdjIdLists, it, sizeof(uint64_t));
            it += sizeof(uint64_t);
            memcpy(&header.offsetAdjIdTranspose, it, sizeof(uint64_t));
            it += sizeof(uint64_t);
        }
    };

    readHeader();
    nodes = header.nodes;
    DEBUG("# nodes here = ", nodes);
    chunks = header.chunks;
    DEBUG("# chunks here = ", chunks);
    if (weightFormat == nkbg::WeightFormat::NONE) {
        weighted = false;
    } else {
        weighted = true;
    }
    // The node ids 0, ..., nodes - 1 must be representable by NodeT and differ from
    // NullNodeId<NodeT>, which marks "no node". Since NullNodeId<NodeT> is the largest NodeT, this
    // allows at most max() nodes.
    static_assert(NullNodeId<NodeT> == std::numeric_limits<NodeT>::max());
    constexpr auto maxNumberOfNodes = static_cast<count>(std::numeric_limits<NodeT>::max());
    if (nodes > maxNumberOfNodes) {
        throw std::runtime_error("The file has " + std::to_string(nodes)
                                 + " nodes, but the node type of the graph can hold at most "
                                 + std::to_string(maxNumberOfNodes) + " nodes");
    }
    GraphT G(nodes, weighted, directed);

    if (indexed)
        G.indexEdges();
    if (nodes == 0) {
        return G;
    }
    // Read base data.
    std::vector<uint8_t> nodeFlags;
    const char *baseIt = startIt + header.offsetBaseData;
    for (uint64_t i = 0; i < nodes; i++) {
        uint8_t flag;
        memcpy(&flag, baseIt, sizeof(uint8_t));
        baseIt += sizeof(uint8_t);
        if (flag & nkbg::DELETED_BIT) {
            G.removeNode(static_cast<NodeT>(i));
        }
    }

    std::vector<uint64_t> firstVert;
    firstVert.push_back(0);
    for (uint64_t ch = 1; ch < chunks; ch++) {
        uint64_t first;
        memcpy(&first, baseIt, sizeof(uint64_t));
        baseIt += sizeof(uint64_t);
        firstVert.push_back(first);
    }
    firstVert.push_back(nodes);
    // Read adjacency lists.
    const char *adjIt = startIt + header.offsetAdjLists;
    const char *transpIt = startIt + header.offsetAdjTranspose;
    const char *adjWghtIt = startIt + header.offsetWeightLists;
    const char *transpWghtIt = startIt + header.offsetWeightTranspose;
    const char *adjIdIt = (version > 2) ? startIt + header.offsetAdjIdLists : nullptr;
    const char *transpIdIt = (version > 2) ? startIt + header.offsetAdjIdTranspose : nullptr;
    uint64_t adjListSize;
    memcpy(&adjListSize, adjIt + (chunks - 1) * sizeof(uint64_t), sizeof(uint64_t));
    uint64_t transposeListSize;
    memcpy(&transposeListSize, transpIt + (chunks - 1) * sizeof(uint64_t), sizeof(uint64_t));

    if (!directed) {
        assert(adjListSize == transposeListSize);
    }
    G.setEdgeCount(unsafe, adjListSize);

    std::atomic<count> selfLoops{0};
    // The largest edge id read so far plus one, which is the upper edge id bound of the graph.
    std::atomic<edgeid> omega{0};
    // Exceptions must not leave the parallel region, so a weight that does not fit EdgeWeightT is
    // only reported after it.
    std::atomic<bool> invalidWeight{false};

    // Raises omega to id + 1 if it is smaller.
    const auto updateOmega = [&](edgeid id) -> void {
        edgeid current = omega.load(std::memory_order_relaxed);
        while (current < id + 1
               && !omega.compare_exchange_weak(current, id + 1, std::memory_order_relaxed)) {
        }
    };

    // Converts a weight of the file, or flags it if EdgeWeightT cannot represent it.
    const auto convertWeight = [&](auto value) -> EdgeWeightT {
        EdgeWeightT weight{1};
        if (!toEdgeWeight(value, weight)) {
            invalidWeight.store(true, std::memory_order_relaxed);
        }
        return weight;
    };

    auto constructGraph = [&](uint64_t c) -> void {
        const uint64_t vertex = firstVert[c];
        uint64_t off = 0;
        uint64_t transpOff = 0;
        uint64_t wghtOff = 0;
        uint64_t transWghtOff = 0;
        uint64_t indexOff = 0;
        uint64_t transIndexOff = 0;
        if (vertex) {
            memcpy(&off, adjIt + (c - 1) * sizeof(uint64_t), sizeof(uint64_t));
            memcpy(&transpOff, transpIt + (c - 1) * sizeof(uint64_t), sizeof(uint64_t));
            memcpy(&wghtOff, adjWghtIt + (c - 1) * sizeof(uint64_t), sizeof(uint64_t));
            memcpy(&transWghtOff, transpWghtIt + (c - 1) * sizeof(uint64_t), sizeof(uint64_t));
            if (indexed) {
                memcpy(&indexOff, adjIdIt + (c - 1) * sizeof(uint64_t), sizeof(uint64_t));
                memcpy(&transIndexOff, transpIdIt + (c - 1) * sizeof(uint64_t), sizeof(uint64_t));
            }
        }

        off += (chunks - 1) * sizeof(uint64_t) + sizeof(uint64_t);
        transpOff += (chunks - 1) * sizeof(uint64_t) + sizeof(uint64_t);
        wghtOff += (chunks - 1) * sizeof(uint64_t);
        transWghtOff += (chunks - 1) * sizeof(uint64_t);
        indexOff += (chunks - 1) * sizeof(uint64_t);
        transIndexOff += (chunks - 1) * sizeof(uint64_t);
        uint64_t n = firstVert[c + 1] - firstVert[c];

        // Helpers to consume weights/ids even if we skip graph construction for a deleted node.
        auto consumeAdjWeight = [&]() -> void {
            switch (weightFormat) {
            case nkbg::WeightFormat::VARINT:
            case nkbg::WeightFormat::SIGNED_VARINT: {
                uint64_t tmp;
                wghtOff +=
                    nkbg::varIntDecode(reinterpret_cast<const uint8_t *>(adjWghtIt + wghtOff), tmp);
            } break;
            case nkbg::WeightFormat::DOUBLE:
                wghtOff += sizeof(double);
                break;
            case nkbg::WeightFormat::FLOAT:
                wghtOff += sizeof(float);
                break;
            case nkbg::WeightFormat::NONE:
                break;
            }
        };

        auto consumeTranspWeight = [&]() -> void {
            switch (weightFormat) {
            case nkbg::WeightFormat::VARINT:
            case nkbg::WeightFormat::SIGNED_VARINT: {
                uint64_t tmp;
                transWghtOff += nkbg::varIntDecode(
                    reinterpret_cast<const uint8_t *>(transpWghtIt + transWghtOff), tmp);
            } break;
            case nkbg::WeightFormat::DOUBLE:
                transWghtOff += sizeof(double);
                break;
            case nkbg::WeightFormat::FLOAT:
                transWghtOff += sizeof(float);
                break;
            case nkbg::WeightFormat::NONE:
                break;
            }
        };

        auto consumeAdjId = [&]() -> void {
            if (!indexed)
                return;
            edgeid tmpId;
            indexOff +=
                nkbg::varIntDecode(reinterpret_cast<const uint8_t *>(adjIdIt + indexOff), tmpId);
        };

        auto consumeTranspId = [&]() -> void {
            if (!indexed)
                return;
            edgeid tmpId;
            transIndexOff += nkbg::varIntDecode(
                reinterpret_cast<const uint8_t *>(transpIdIt + transIndexOff), tmpId);
        };
        uint64_t ignored = 0;
        for (uint64_t i = 0; i < n; i++) {
            // curr < nodes, which fits into NodeT (checked above).
            const auto curr = static_cast<NodeT>(vertex + i);
            uint64_t outNbrs;
            off += nkbg::varIntDecode(reinterpret_cast<const uint8_t *>(adjIt + off), outNbrs);
            uint64_t inNbrs;
            transpOff +=
                nkbg::varIntDecode(reinterpret_cast<const uint8_t *>(transpIt + transpOff), inNbrs);

            // If the node doesn't exist (hole), we must NOT preallocate/add edges,
            // but we MUST still consume the encoded neighbor/weight/id payload to keep offsets
            // aligned.
            if (!G.hasNode(curr)) {
                for (uint64_t j = 0; j < outNbrs; ++j) {
                    off +=
                        nkbg::varIntDecode(reinterpret_cast<const uint8_t *>(adjIt + off), ignored);
                    consumeAdjWeight();
                    consumeAdjId();
                }
                for (uint64_t j = 0; j < inNbrs; ++j) {
                    transpOff += nkbg::varIntDecode(
                        reinterpret_cast<const uint8_t *>(transpIt + transpOff), ignored);
                    consumeTranspWeight();
                    consumeTranspId();
                }
                continue;
            }

            if (!directed) {
                G.preallocateUndirected(curr, outNbrs + inNbrs);
            } else {
                G.preallocateDirected(curr, outNbrs, inNbrs);
            }
            // Read adjacency lists.
            for (uint64_t j = 0; j < outNbrs; j++) {
                edgeid id = 0;
                uint64_t neighbor;
                EdgeWeightT weight{1};
                off += nkbg::varIntDecode(reinterpret_cast<const uint8_t *>(adjIt + off), neighbor);
                // neighbor < nodes, which fits into NodeT.
                const auto add = static_cast<NodeT>(neighbor);
                switch (weightFormat) {
                case nkbg::WeightFormat::VARINT: {
                    uint64_t unsignedWeight;
                    wghtOff += nkbg::varIntDecode(
                        reinterpret_cast<const uint8_t *>(adjWghtIt + wghtOff), unsignedWeight);
                    weight = convertWeight(unsignedWeight);
                } break;
                case nkbg::WeightFormat::DOUBLE: {
                    double doubleWeight;
                    memcpy(&doubleWeight, adjWghtIt + wghtOff, sizeof(double));
                    wghtOff += sizeof(double);
                    weight = convertWeight(doubleWeight);
                } break;
                case nkbg::WeightFormat::SIGNED_VARINT: {
                    uint64_t unsignedWeight;
                    wghtOff += nkbg::varIntDecode(
                        reinterpret_cast<const uint8_t *>(adjWghtIt + wghtOff), unsignedWeight);
                    weight = convertWeight(nkbg::zigzagDecode(unsignedWeight));
                } break;
                case nkbg::WeightFormat::FLOAT: {
                    float floatWeight;
                    memcpy(&floatWeight, adjWghtIt + wghtOff, sizeof(float));
                    wghtOff += sizeof(float);
                    weight = convertWeight(floatWeight);
                } break;
                case nkbg::WeightFormat::NONE:
                    break;
                }
                if (indexed) {
                    indexOff += nkbg::varIntDecode(
                        reinterpret_cast<const uint8_t *>(adjIdIt + indexOff), id);
                    updateOmega(id);
                }
                if (!directed) {
                    if (!G.addPartialEdge(unsafe, curr, add, weight, id, true))
                        WARN("Not adding edge ", curr, "-", add, " since it is already present.");
                } else {
                    if (!G.addPartialOutEdge(unsafe, curr, add, weight, id, true))
                        WARN("Not adding edge ", curr, "-", add, " since it is already present.");
                }
                if (curr == add) {
                    selfLoops.fetch_add(1, std::memory_order_relaxed);
                }
            }
            // Read transpose lists.
            for (uint64_t j = 0; j < inNbrs; j++) {
                uint64_t neighbor;
                edgeid id = 0;
                EdgeWeightT weight{1};
                transpOff += nkbg::varIntDecode(
                    reinterpret_cast<const uint8_t *>(transpIt + transpOff), neighbor);
                // neighbor < nodes, which fits into NodeT.
                const auto add = static_cast<NodeT>(neighbor);
                switch (weightFormat) {
                case nkbg::WeightFormat::VARINT: {
                    uint64_t unsignedWeight;
                    transWghtOff += nkbg::varIntDecode(
                        reinterpret_cast<const uint8_t *>(transpWghtIt + transWghtOff),
                        unsignedWeight);
                    weight = convertWeight(unsignedWeight);
                } break;
                case nkbg::WeightFormat::DOUBLE: {
                    double doubleWeight;
                    memcpy(&doubleWeight, transpWghtIt + transWghtOff, sizeof(double));
                    transWghtOff += sizeof(double);
                    weight = convertWeight(doubleWeight);
                } break;
                case nkbg::WeightFormat::SIGNED_VARINT: {
                    uint64_t unsignedWeight;
                    transWghtOff += nkbg::varIntDecode(
                        reinterpret_cast<const uint8_t *>(transpWghtIt + transWghtOff),
                        unsignedWeight);
                    weight = convertWeight(nkbg::zigzagDecode(unsignedWeight));
                } break;
                case nkbg::WeightFormat::FLOAT: {
                    float floatWeight;
                    memcpy(&floatWeight, transpWghtIt + transWghtOff, sizeof(float));
                    transWghtOff += sizeof(float);
                    weight = convertWeight(floatWeight);
                } break;
                case nkbg::WeightFormat::NONE:
                    break;
                }
                if (indexed) {
                    transIndexOff += nkbg::varIntDecode(
                        reinterpret_cast<const uint8_t *>(transpIdIt + transIndexOff), id);
                    updateOmega(id);
                }
                if (!directed) {
                    if (curr != add) {
                        if (G.addPartialEdge(unsafe, curr, add, weight, id, true))
                            WARN("Not adding edge ", curr, "-", add,
                                 " since it is already present.");
                    }
                } else {
                    if (!G.addPartialInEdge(unsafe, curr, add, weight, id, true))
                        WARN("Not adding edge ", curr, "-", add, " since it is already present.");
                }
            }
        }
    };

// create graph
#pragma omp parallel for
    for (omp_index c = 0; c < static_cast<omp_index>(chunks); c++) {
        constructGraph(c);
    }
    if (invalidWeight) {
        throw std::runtime_error("The file has an edge weight that the edge weight type of the "
                                 "graph cannot represent");
    }
    G.setNumberOfSelfLoops(unsafe, selfLoops);
    if (indexed)
        G.setUpperEdgeIdBound(unsafe, omega);
    return G;
}
} // namespace NetworKit

#endif // NETWORKIT_IO_NETWORKIT_BINARY_READER_IMPL_HPP_
