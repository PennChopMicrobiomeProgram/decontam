import unittest

from decontamlib import main


class Test_tools_creation(unittest.TestCase):
    def setUp(self):
        self.params = {
            "bmtagger": {"bitmask": "path", "srprism": "path2"},
            "bowtie": {"index": "path2"},
            "random_human": {"percent_human": 30},
            }

    def test_empty_tool_names(self):
        self.assertEqual(main.create_tools([], self.params), [])

    def test_unknown_tool(self):
        tool_names = ["unknown_and_never_existing_tool"]
        self.assertRaises(KeyError, main.create_tools, tool_names, self.params)

    def test_can_create_tools(self):
        tool_names = [ "bmtagger", "all_human", "none_human", "random_human", "bowtie" ]
        tool_runners = main.create_tools(tool_names, self.params)
        self.assertEqual(len(tool_runners), len(tool_names))
        self.assertTrue(isinstance(tool_runners[0], main.Bmtagger))
        self.assertTrue(isinstance(tool_runners[1], main.All_human))
        self.assertTrue(isinstance(tool_runners[2], main.None_human))
        self.assertTrue(isinstance(tool_runners[3], main.Random_human))
        self.assertTrue(isinstance(tool_runners[4], main.Bowtie))

