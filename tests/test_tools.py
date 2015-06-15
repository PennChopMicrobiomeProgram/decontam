import unittest
from cStringIO import StringIO

from decontamlib.tools import (
    _FilteringTool, Bmtagger, None_human,
    )


bmtagger_output = """#read-id	#tag
M03249:9:000000000-ABY6B:1:1101:12869:1660	F
M03249:9:000000000-ABY6B:1:1101:13579:1677	U
M03249:9:000000000-ABY6B:1:1101:16098:1733	F
M03249:9:000000000-ABY6B:1:1101:13597:1838	F
M03249:9:000000000-ABY6B:1:1101:11973:1857	H
"""


class Test_all_human(unittest.TestCase):
    def setUp(self):
        self.human_annotation = list(
            Bmtagger._parse_bmtagger_output(StringIO(bmtagger_output)))

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

class FilteringToolTests(unittest.TestCase):
    def test_get_args(self):
        self.assertEqual(_FilteringTool.get_argnames(), ["index"])
        self.assertEqual(None_human.get_argnames(), [])


if __name__ == "__main__":
    unittest.main()
