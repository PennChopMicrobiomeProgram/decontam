import unittest
import mock
import tempfile
import json
from cStringIO import StringIO

from .. import tools
from ..human_filtering_tools import Bmtagger
from ..human_filtering_tools import All_human
from ..human_filtering_tools import None_human
from ..human_filtering_tools import Random_human
from ..human_filtering_tools import Bowtie

class Test_tools_creation(unittest.TestCase):
    
    def setUp(self):
        self.sample_json_fh = StringIO(
            '{"bmtagger" : {"bitmask" : "path", "srprism": "path2"}, "bowtie" : {"index" : "path2"}, "random_human" : {"percent_human" : 30}}')
        self.params = json.load(self.sample_json_fh)

    def test_empty_tool_names(self):
        self.assertRaises(ValueError, tools.create_tools, [], self.params)

    def test_unknown_tool(self):
        tool_names = ["unknown_and_never_existing_tool"]
        self.assertRaises(KeyError, tools.create_tools, tool_names, self.params)

    def test_can_create_tools(self):
        tool_names = [ "bmtagger", "all_human", "none_human", "random_human", "bowtie" ]
        tool_runners = tools.create_tools(tool_names, self.params)
        self.assertEqual(len(tool_runners), len(tool_names))
        self.assertTrue(isinstance(tool_runners[0], Bmtagger))
        self.assertTrue(isinstance(tool_runners[1], All_human))
        self.assertTrue(isinstance(tool_runners[2], None_human))
        self.assertTrue(isinstance(tool_runners[3], Random_human))
        self.assertTrue(isinstance(tool_runners[4], Bowtie))

class Test_parameter_file_default(unittest.TestCase):

    def setUp(self):
        self.real_function = tools.json.load
        sample_json_fh = StringIO('{"bmtagger_DEF" : {"bitmask_DEF" : "path"}}')        
        params = json.load(sample_json_fh)
        tools.json.load = mock.Mock(return_value = params)

    def tearDown(self):
        tools.json.load = self.real_function

    def test_read_json_default(self):
        parameters = tools.get_parameters_for_tools()
        self.assertTrue(parameters["bmtagger_DEF"] is not None)
        self.assertTrue(parameters["bmtagger_DEF"]["bitmask_DEF"] is not None)

class Test_parameter_file(unittest.TestCase):
    
    def test_read_json(self):
        sample_json = tempfile.NamedTemporaryFile()        
        sample_json.write('{"bmtagger" : {"bitmask" : "path"}}')        
        sample_json.seek(0)
        parameters = tools.get_parameters_for_tools(sample_json.name)
        self.assertTrue(parameters["bmtagger"] is not None)
        self.assertTrue(parameters["bmtagger"]["bitmask"] is not None)
