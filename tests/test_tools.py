import unittest
import mock
import tempfile
import json
from cStringIO import StringIO

from .. import tools

class Test_tools_creation(unittest.TestCase):
    def setUp(self):
        self.params = {
            "bmtagger": {"bitmask": "path", "srprism": "path2"},
            "bowtie": {"index": "path2"},
            "random_human": {"percent_human": 30},
            }

    def test_empty_tool_names(self):
        self.assertEqual(tools.create_tools([], self.params), [])

    def test_unknown_tool(self):
        tool_names = ["unknown_and_never_existing_tool"]
        self.assertRaises(KeyError, tools.create_tools, tool_names, self.params)

    def test_can_create_tools(self):
        tool_names = [ "bmtagger", "all_human", "none_human", "random_human", "bowtie" ]
        tool_runners = tools.create_tools(tool_names, self.params)
        self.assertEqual(len(tool_runners), len(tool_names))
        self.assertTrue(isinstance(tool_runners[0], tools.Bmtagger))
        self.assertTrue(isinstance(tool_runners[1], tools.All_human))
        self.assertTrue(isinstance(tool_runners[2], tools.None_human))
        self.assertTrue(isinstance(tool_runners[3], tools.Random_human))
        self.assertTrue(isinstance(tool_runners[4], tools.Bowtie))

