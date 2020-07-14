// networkit-format
#ifndef NETWORKIT_COMPONENTS_DYN_CONNECTED_COMPONENTS_HPP_
#define NETWORKIT_COMPONENTS_DYN_CONNECTED_COMPONENTS_HPP_

#include <map>
#include <memory>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/base/DynAlgorithm.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/components/DynConnectedComponentsImpl.hpp>
#include <networkit/dynamics/GraphEvent.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup components
 * Determines and updates the connected components of an undirected graph.
 */
// class DynConnectedComponentsGeneral : public ConnectedComponentsGeneral {
//
// protected:
//    std::unique_ptr<DynConnectedComponentsImpl> impl;
//};
// /**
//  * Updates the connected components after an edge insertion or deletion.
//  *
//  * @param[in]	event	The event that happened (edge deletion or
//  * insertion).
//  */
// virtual void update(GraphEvent e) = 0;

// /**
//  * Updates the connected components after a batch of edge insertions or
//  * deletions.
//  *
//  * @param[in] batch	A vector that contains a batch of edge insertions or
//  *					deletions.
//  */
// virtual void updateBatch(const std::vector<GraphEvent> &batch) = 0;

/**
 * Returns the the component in which node @a u is.
 *
 * @param[in]	u	The node.
 */
//    virtual count componentOfNode(node u) const = 0;
//
//    /**
//     * Returns the map from component to size.
//     */
//    virtual std::map<index, count> getComponentSizes();

//    /**
//     * @return Vector of components, each stored as (unordered) set of nodes.
//     */
//    std::vector<std::vector<node>> getComponents();

class DynConnectedComponents final : public ConnectedComponentsGeneral, public DynAlgorithm {
public:
    DynConnectedComponents(const Graph &G);

    void run() override;
    count numberOfComponents() const override;
    count componentOfNode(node u) const override;
    Partition getPartition() const override;
    std::map<index, count> getComponentSizes() const override;
    std::vector<std::vector<node>> getComponents() const override;
    void update(GraphEvent e) override;
    void updateBatch(const std::vector<GraphEvent> &batch) override;

private:
    std::unique_ptr<DynConnectedComponentsDetails::DynConnectedComponentsImpl> impl;
};

class DynWeaklyConnectedComponents final : public ConnectedComponentsGeneral, public DynAlgorithm {
public:
    DynWeaklyConnectedComponents(const Graph &G);

    void run() override;
    count numberOfComponents() const override;
    count componentOfNode(node u) const override;
    Partition getPartition() const override;
    std::map<index, count> getComponentSizes() const override;
    std::vector<std::vector<node>> getComponents() const override;
    void update(GraphEvent e) override;
    void updateBatch(const std::vector<GraphEvent> &batch) override;

private:
    std::unique_ptr<DynConnectedComponentsDetails::DynConnectedComponentsImpl> impl;
};

} // namespace NetworKit

#endif // NETWORKIT_COMPONENTS_DYN_CONNECTED_COMPONENTS_HPP_
