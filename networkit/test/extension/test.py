import unittest
import os
import tempfile


class TestExtMETISGraphReader(unittest.TestCase):
	def setUp(self):
		from networkit.graph import Graph
		self.g = Graph(5)
		self.g.addEdge(0,1)
		self.g.addEdge(0,2)
		self.g.addEdge(0,3)
		self.g.addEdge(0,4)

	def test_readAndWrite(self):
		from networkit.graphio import METISGraphReader
		from networkit.graphio import METISGraphWriter
		w = METISGraphWriter()
		with tempfile.NamedTemporaryFile(mode='w') as tmpfile:
			w.write(self.g, tmpfile.name)
			self.assertTrue(os.path.isfile(tmpfile.name))
			r = METISGraphReader()
			testg = r.read(tmpfile.name)
			self.assertEqual(self.g.numberOfNodes(), testg.numberOfNodes())
			self.assertEqual(self.g.numberOfEdges(), testg.numberOfEdges())

