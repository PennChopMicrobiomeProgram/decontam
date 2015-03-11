import unittest
from cStringIO import StringIO

from ..human_filtering_tools import Bmtagger
from ..human_filtering_tools import All_human
from ..human_filtering_tools import None_human

bmtagger_output = (
"""#read-id	#tag
M03249:9:000000000-ABY6B:1:1101:12869:1660	F
M03249:9:000000000-ABY6B:1:1101:13579:1677	U
M03249:9:000000000-ABY6B:1:1101:16098:1733	F
M03249:9:000000000-ABY6B:1:1101:13597:1838	F
M03249:9:000000000-ABY6B:1:1101:11973:1857	H
"""
)

class Test_all_human(unittest.TestCase):

    def setUp(self):
        parameters = { "bitmask" : "path"}     
        self.bmtagger = Bmtagger(parameters)
        self.human_annotation = self.bmtagger.parse_bmtagger_output((StringIO(bmtagger_output)))

    def test_parse_number_of_reads(self):
        self.assertEqual(len(self.human_annotation), 5) 

    def test_parse_has_two_columns(self):
        self.assertEqual(len(self.human_annotation[0]), 2) 

    def test_parse_column_values(self):
        self.assertEqual(self.human_annotation[0][0], 
            "M03249:9:000000000-ABY6B:1:1101:12869:1660") 
        self.assertEqual(self.human_annotation[0][1], 0) 

        self.assertEqual(self.human_annotation[4][0], 
            "M03249:9:000000000-ABY6B:1:1101:11973:1857") 
        self.assertEqual(self.human_annotation[4][1], 1) 

class Test_none_human(unittest.TestCase):
    pass

class Test_random(unittest.TestCase):
    pass

class Test_Bmtagger(unittest.TestCase):
    pass
