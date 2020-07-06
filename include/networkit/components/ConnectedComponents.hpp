/*
 * ConnectedComponents.hpp
 *
 *  Created on: Dec 16, 2013
 *      Author: cls
 */

// networkit-format

#ifndef NETWORKIT_COMPONENTS_CONNECTED_COMPONENTS_HPP_
#define NETWORKIT_COMPONENTS_CONNECTED_COMPONENTS_HPP_

#include <map>
#include <memory>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/components/ConnectedComponentsImpl.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup components
 * Determines the (weakly) connected components of a directed/undirected graph.
 */
class ConnectedComponentsGeneral : public Algorithm {
public:
    /**
     * Get the number of connected components.
     *
     * @return The number of connected components.
     */
    virtual count numberOfComponents() const = 0;

    /**
     * Get the the component in which node @a u is situated.
     *
     * @param[in]	u	The node whose component is asked for.
     */
    virtual count componentOfNode(node u) const = 0;

    /**
     * Get a Partition that represents the components.
     *
     * @return A partition representing the found components.
     */
    virtual Partition getPartition() const = 0;

    /**
     *Return the map from component to size
     */
    virtual std::map<index, count> getComponentSizes() const = 0;

    /**
     * @return Vector of components, each stored as (unordered) set of nodes.
     */
    virtual std::vector<std::vector<node>> getComponents() const = 0;

protected:
    ConnectedComponentsGeneral(const Graph &G, bool weaklyCC = false)
        : impl(new ConnectedComponentsDetails::ConnectedComponentsImpl{G, weaklyCC}) {}

    std::unique_ptr<ConnectedComponentsDetails::ConnectedComponentsImpl> impl;
};

class ConnectedComponents final : public ConnectedComponentsGeneral {
public:
    ConnectedComponents(const Graph &G);

    void run() override;
    count numberOfComponents() const override;
    count componentOfNode(node u) const override;
    Partition getPartition() const override;
    std::map<index, count> getComponentSizes() const override;
    std::vector<std::vector<node>> getComponents() const override;

    /**
     * Constructs a new graph that contains only the nodes inside the largest
     * connected component.
     * @param G            The input graph.
     * @param compactGraph If true, the node ids of the output graph will be compacted
     * (i.e. re-numbered from 0 to n-1). If false, the node ids will not be changed.
     */
    static Graph extractLargestConnectedComponent(const Graph &G, bool compactGraph = false);
};

class WeaklyConnectedComponents final : public ConnectedComponentsGeneral {
public:
    WeaklyConnectedComponents(const Graph &G);

    void run() override;
    count numberOfComponents() const override;
    count componentOfNode(node u) const override;
    Partition getPartition() const override;
    std::map<index, count> getComponentSizes() const override;
    std::vector<std::vector<node>> getComponents() const override;

    /**
     * Constructs a new graph that contains only the nodes inside the largest
     * weakly connected component.
     * @param G            The input graph.
     * @param compactGraph If true, the node ids of the output graph will be compacted
     * (i.e. re-numbered from 0 to n-1). If false, the node ids will not be changed.
     */
    static Graph extractLargestWeaklyConnectedComponent(const Graph &G, bool compactGraph = false);
};

} // namespace NetworKit

#endif // NETWORKIT_COMPONENTS_CONNECTED_COMPONENTS_HPP_
