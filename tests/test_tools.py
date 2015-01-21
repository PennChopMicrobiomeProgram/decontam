import unittest

from .. import tools
from ..human_filtering_tools import Bmtagger
from ..human_filtering_tools import All_human
from ..human_filtering_tools import None_human

class Test_tools_creation(unittest.TestCase):

    def test_empty_tool_names(self):
        self.assertRaises(ValueError, tools.create_tools, [])

    def test_unknown_tool(self):
        tool_names = ["unknown_and_never_existing_tool"]
        self.assertRaises(ValueError, tools.create_tools, tool_names)

    def test_can_create_tools(self):
        tool_names = [ "bmtagger", "all_human", "none_human" ]
        tool_runners = tools.create_tools(tool_names) 
        self.assertEqual(len(tool_runners), len(tool_names))
        self.assertTrue(isinstance(tool_runners[0], Bmtagger))
        self.assertTrue(isinstance(tool_runners[1], All_human))
        self.assertTrue(isinstance(tool_runners[2], None_human))



