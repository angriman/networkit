#!/usr/bin/env python3
import unittest
import os

import networkit as nk

class TestDistance(unittest.TestCase):

	def setUp(self):
		self.L = nk.readGraph("input/looptest1.gml", nk.Format.GML) #without self-loops
		self.LL = nk.readGraph("input/looptest2.gml", nk.Format.GML) #with self-loops sprinkled in
		self.G = nk.readGraph("input/karate.graph", nk.Format.METIS)

	def testDistanceDiameter(self):
		D = nk.distance.Diameter(self.LL, nk.distance.DiameterAlgo.EstimatedRange, error = 0.1)
		D.run()
		D = nk.distance.Diameter(self.LL, nk.distance.DiameterAlgo.EstimatedSamples, nSamples = 5)
		D.run()
		D = nk.distance.Diameter(self.LL, nk.distance.DiameterAlgo.Exact)
		D.run()


	def testDistanceEccentricity(self):
		E = nk.distance.Eccentricity()
		E.getValue(self.LL, 0)


	def testDistanceEffectiveDiameter(self):
		algo = nk.distance.EffectiveDiameter(self.L)
		algo.run()
		algo = nk.distance.EffectiveDiameter(self.LL)
		algo.run()


	def testDistanceApproxEffectiveDiameter(self):
		algo = nk.distance.EffectiveDiameterApproximation(self.L)
		algo.run()
		algo = nk.distance.EffectiveDiameterApproximation(self.LL)
		algo.run()


	def testDistanceApproxHopPlot(self):
		algo = nk.distance.HopPlotApproximation(self.L)
		algo.run()
		algo = nk.distance.HopPlotApproximation(self.LL)
		algo.run()


	def testDistanceNeighborhoodFunction(self):
		algo = nk.distance.NeighborhoodFunction(self.L)
		algo.run()
		algo = nk.distance.NeighborhoodFunction(self.LL)
		algo.run()


	def testDistanceApproxNeighborhoodFunction(self):
		algo = nk.distance.NeighborhoodFunctionApproximation(self.L)
		algo.run()
		algo = nk.distance.NeighborhoodFunctionApproximation(self.LL)
		algo.run()

	def testDistanceAStar(self):
		# Builds a mesh graph with the given number of rows and columns
		def buildMesh(rows, cols):
			G = nk.Graph(rows * cols, False, False)
			for i in range(rows):
				for j in range(cols):
					if j < cols - 1:
						G.addEdge(i * cols + j, i * cols + j + 1)
					if i < rows - 1:
						G.addEdge(i * cols + j, (i + 1) * cols + j)
			return G

		# Test the AStar algorithm on a mesh with the given number of rows and columns
		def testMesh(rows, cols):
			G = buildMesh(rows, cols)

			# Test A* on the given source-target pair
			def testPair(s, t):

				# Some distance heuristics:

				# Always returns 0, A* degenerates to Dijkstra
				def zeroDist(u):
					return 0

				# Returns the exact distance from u to the target
				def exactDist(u):
					rowU = int(u / cols)
					colU = int(u % cols)
					rowT = int(t / cols)
					colT = int(t % cols)
					return abs(rowU - rowT) + abs(colU - colT)

				# Returns the eucledian distance from u to the target
				def eucledianDist(u):
					rowT = int(t / cols)
					colT = int(t % cols)
					rowDiff = abs(int(u / cols) - rowT)
					colDiff = abs(int(u % cols) - rowT)
					return (rowDiff**2 + colDiff**2)**.5

				# Use BFS as ground truth
				bfs = nk.distance.BFS(G, s, True, False, t).run()

				# Test A* on all the heuristics
				for heu in [zeroDist, exactDist, eucledianDist]:
					heuristics = [heu(u) for u in range(G.numberOfNodes())]
					astar = nk.distance.AStar(G, heuristics, s, t, True)
					astar.run()

					# Test distance of target
					self.assertEqual(astar.getDistance(), bfs.distance(t))

					# Test path
					path = astar.getPath()
					self.assertEqual(len(path), len(bfs.getPath(t)) - 2)
					if len(path) == 0:
						continue
					for i in range(len(path) - 1):
						self.assertTrue(G.hasEdge(path[i], path[i + 1]))

			# Iterate over all possible source-target pairs
			G.forNodePairs(testPair)

		# Test some meshes
		testMesh(10, 10)
		testMesh(21, 5)
		testMesh(9, 18)
		testMesh(7, 1)


if __name__ == "__main__":
	unittest.main()
