import unittest
from cStringIO import StringIO

from .. import parser

class Test_Tools(unittest.TestCase):

    def test_empty_file(self):
        self.assertRaises(IOError, parser.parse_tool_names, StringIO(""))

    def test_valid_file(self):
        tools = parser.parse_tool_names(StringIO("tool1\ntool2\n"))
        self.assertEquals(len(tools), 2)
        self.assertEquals(tools[0], "tool1")
        self.assertEquals(tools[1], "tool2")


class  Test_Samples(unittest.TestCase):
    
    def test_empty_file(self):
        self.assertTrue(False)

    def test_incorrect_file_format(self):
        self.assertTrue(False)

    def test_file_does_not_exist(self):
        self.assertTrue(False)

    def test_valid_file(self):
        self.assertTrue(False)
