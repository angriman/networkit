#ifndef NETWORKIT_IO_MTX_PARSER_HPP_
#define NETWORKIT_IO_MTX_PARSER_HPP_

#include <charconv>
#include <cstdlib>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

namespace MTXParserDetails {

/**
 * Parses the next whitespace-separated token of @a line into a T. Integral numbers must be
 * representable by T exactly: unlike a detour through double, which rounds integers above 2^53,
 * they are parsed directly into T.
 *
 * @throws std::runtime_error if the line has no further token, or the token is not a valid T.
 */
template <class T>
T parseValue(std::istringstream &line, const char *what) {
    std::string token;
    if (!(line >> token)) {
        throw std::runtime_error(std::string("MTX line ended where a ") + what + " was expected");
    }
    T value{};
    const char *const first = token.data();
    const char *const last = first + token.size();
    const char *numberEnd = first;
    if constexpr (std::is_integral_v<T>) {
        const auto [ptr, ec] = std::from_chars(first, last, value);
        if (ec == std::errc{}) {
            numberEnd = ptr;
        }
    } else {
        // std::from_chars for floating-point types is too recent for all supported compilers.
        char *past;
        value = static_cast<T>(std::strtod(first, &past));
        numberEnd = past;
    }
    // E.g., the fraction of 1.5 or the sign of -1 for unsigned types are not part of the number.
    if (numberEnd != last) {
        throw std::runtime_error("MTX " + std::string(what) + " " + token
                                 + " is not valid for the requested type");
    }
    return value;
}

} // namespace MTXParserDetails

/**
 * @ingroup io
 * Parser for the MTX file format.
 */
class MTXParser final {

    std::ifstream graphFile;
    std::string currentLine;

public:
    enum class Object { Matrix, Vector };
    enum class Format { Coordinate, Array };
    enum class Field { Real, Double, Complex, Integer, Pattern };
    enum class Symmetry { General, Symmetric, SkewSymmetric, Hermitian };

    std::unordered_map<std::string, Object> objectMap = {{"matrix", Object::Matrix},
                                                         {"vector", Object::Vector}};
    std::unordered_map<std::string, Format> formatMap = {{"coordinate", Format::Coordinate},
                                                         {"array", Format::Array}};
    std::unordered_map<std::string, Field> fieldMap = {{"real", Field::Real},
                                                       {"double", Field::Double},
                                                       {"integer", Field::Integer},
                                                       {"pattern", Field::Pattern}};
    std::unordered_map<std::string, Symmetry> symmetryMap = {
        {"general", Symmetry::General},
        {"symmetric", Symmetry::Symmetric},
        {"skew-symmetric", Symmetry::SkewSymmetric}};

    struct MTXHeader {
        Object object;
        Format format;
        Field field;
        Symmetry symmetry;
    };

    struct MatrixSize {
        count rows;
        count columns;
        count nonzeros;

        MatrixSize(count r, count c, count nz) : rows(r), columns(c), nonzeros(nz) {}
    };

    MTXParser(const std::string &path);

    /**
     * Get the MTX graph file header.
     */
    MTXHeader getHeader();

    /**
     * Get the MTX graph file matrix size header.
     */
    MatrixSize getMatrixSize();

    /**
     * Get the (weighted) edge from the next line in the MTX graph file if present. The node ids
     * are 0-based, and edges of unweighted files get weight 1.
     *
     * @param weighted
     * @throws std::runtime_error if the line is malformed, a node id is 0, or the weight is not a
     * valid EdgeWeightT.
     */
    template <class EdgeWeightT>
    std::optional<WeightedEdgeT<node, EdgeWeightT>> getNext(bool weighted);
};

template <class EdgeWeightT>
std::optional<WeightedEdgeT<node, EdgeWeightT>> MTXParser::getNext(bool weighted) {
    std::string line;
    while (std::getline(graphFile, line)) {
        // check for comment line starting with '%'
        if (!line.empty() && line[0] == '%') {
            throw std::runtime_error(
                "Invalid MTX file structure. No comments allowed after size line.");
        }
        // Skip empty lines.
        if (line.find_first_not_of(" \t\r") == std::string::npos) {
            continue;
        }

        std::istringstream lineStream(line);
        const auto i = MTXParserDetails::parseValue<node>(lineStream, "node id");
        const auto j = MTXParserDetails::parseValue<node>(lineStream, "node id");
        if (i == 0 || j == 0) {
            throw std::runtime_error("MTX node ids start at 1, but the file has node id 0");
        }
        EdgeWeightT w{1};
        if (weighted) {
            w = MTXParserDetails::parseValue<EdgeWeightT>(lineStream, "weight");
        }
        return WeightedEdgeT<node, EdgeWeightT>(i - 1, j - 1, w);
    }
    return std::nullopt;
}
} /* namespace NetworKit */
#endif // NETWORKIT_IO_MTX_PARSER_HPP_
