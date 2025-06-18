#ifndef NETWORKIT_GRAPH_EDGE_LAMBDA_HPP_
#define NETWORKIT_GRAPH_EDGE_LAMBDA_HPP_

#include <cstddef>
#include <type_traits>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/FunctionTraits.hpp>

namespace NetworKit {

namespace EdgeLambdaDetails {

// Return void when the callable does not have an Nth parameter.
template <class TraitsT, size_t N, bool InRange = (N < TraitsT::arity)>
struct SafeArg {
    using type = void;
};

// Enable this specialization only when N is in range.
template <class TraitsT, size_t N>
struct SafeArg<TraitsT, N, true> {
    using type = typename TraitsT::template arg<N>::type;
};

template <class F, size_t N>
using ArgType = std::decay_t<typename SafeArg<Aux::FunctionTraits<F>, N>::type>;

} // namespace EdgeLambdaDetails

/* Compile-time callable dispatcher shared by the graph data structures.
 *
 * Aux::FunctionTraits lets us inspect a callable's arity and parameter types
 * so we can invoke it with the best matching edge callback signature, i.e. one of
 * (u, v, weight, edgeId), (u, v, edgeId), (u, v, weight), (v, weight), (u, v), (v).
 * SafeArg prevents out-of-bounds access to Traits::arg<N> by returning void
 * when the callable has fewer than N + 1 parameters.
 * ArgType normalizes the inspected type with std::decay_t so comparisons work
 * for references, const qualifiers, and other call-compatible forms.
 *
 * The overload selection below is ordered from most specific to least specific.
 * This matters because several lambda signatures overlap, especially when node
 * and EdgeWeightT are the same type.
 *
 * Graphs without edge ids set WithEdgeIds to false: the signatures taking an edge id
 * are then never selected, so they fail to compile instead of silently receiving `id`.
 * This also keeps a (u, v, weight) callback from being mistaken for (u, v, edgeId) when
 * EdgeWeightT is the same type as edgeid.
 *
 * If no branch matches, the static_assert produces a clear compile-time error
 * instead of a long substitution failure.
 */
template <class NodeT, class EdgeWeightT, bool WithEdgeIds = true, class F>
decltype(auto) callEdgeLambda(F &f, NodeT u, NodeT v, EdgeWeightT ew, [[maybe_unused]] edgeid id) {
    using EdgeLambdaDetails::ArgType;
    constexpr size_t arity = Aux::FunctionTraits<F>::arity;

    // Most specific signatures first to avoid ambiguity between overlapping
    // callable shapes and implicit conversions.
    if constexpr (WithEdgeIds && arity >= 3 && std::is_same_v<ArgType<F, 2>, EdgeWeightT>
                  && std::is_same_v<ArgType<F, 3>, edgeid>) {
        return f(u, v, ew, id);
    } else if constexpr (WithEdgeIds && arity >= 2 && std::is_same_v<ArgType<F, 2>, edgeid>
                         && std::is_same_v<ArgType<F, 1>, NodeT>) {
        return f(u, v, id);
    } else if constexpr (arity >= 2 && std::is_same_v<ArgType<F, 2>, EdgeWeightT>) {
        return f(u, v, ew);
    } else if constexpr (arity >= 1 && std::is_same_v<ArgType<F, 1>, EdgeWeightT>
                         && !std::is_same_v<NodeT, EdgeWeightT>) {
        // Treat this as a weight-taking callback only when node ids and weights
        // are distinct types.
        return f(v, ew);
    } else if constexpr (arity >= 1 && std::is_same_v<ArgType<F, 1>, NodeT>) {
        return f(u, v);
    } else if constexpr (arity >= 1 && std::is_same_v<ArgType<F, 1>, EdgeWeightT>) {
        return f(v, ew);
    } else if constexpr (arity >= 1 && std::is_same_v<ArgType<F, 0>, NodeT>) {
        return f(v);
    } else {
        static_assert(!std::is_same_v<F, F>,
                      "Your lambda does not support the required parameters or the parameters "
                      "have the wrong type.");
    }
}

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_EDGE_LAMBDA_HPP_
