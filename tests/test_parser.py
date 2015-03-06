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
        self.assertRaises(IOError, parser.parse_samples, StringIO(""))

    def test_incorrect_file_format(self):
        self.assertRaises(IOError, parser.parse_samples, StringIO("name1\tfastq\n"))

    def test_file_does_not_exist(self):
        parser.os.path.isfile = lambda _ : False 
        self.assertRaises(IOError, parser.parse_samples, 
            StringIO("name1\tAbsent_fastq\tAbsent_annotation\n"))

    def test_valid_file(self):
        parser.os.path.isfile = lambda _ : True 
        samples = parser.parse_samples(
            StringIO("name1\tfastq1\tannotation1\nname2\tfastq2\tannotation2\n"))
        self.assertEqual(len(samples), 2)
        self.assertEqual(samples[0], ("name1", "fastq1", "annotation1"))
        self.assertEqual(samples[1], ("name2", "fastq2", "annotation2"))

class  Test_Annotations(unittest.TestCase):
    
    def test_empty_file(self):
        self.assertRaises(IOError, parser.parse_annotation, StringIO(""))

    def test_incorrect_number_of_columns(self):
        self.assertRaises(IOError, parser.parse_annotation, StringIO("read1\n"))

    def test_incorrect_description_for_human(self):
        self.assertRaises(IOError, parser.parse_annotation, StringIO("read1\tShouldBe[YN]\n"))

    def test_valid_file(self):
        read_annotation = parser.parse_annotation(
            StringIO("read1\tY\nread2\tY\nread3\tN\n"))
        self.assertEqual(len(read_annotation), 3)
        self.assertEqual(read_annotation[0], ("read1", "Y"))
        self.assertEqual(read_annotation[1], ("read2", "Y"))
        self.assertEqual(read_annotation[2], ("read3", "N"))

class Test_Reads(unittest.TestCase):

    def test_reads(self):
        self.assertTrue(False)
