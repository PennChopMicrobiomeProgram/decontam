import unittest
import tempfile
try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

from decontamlib import utils


class TestIdExtract(unittest.TestCase):
    def setUp(self):
        self.fastq = tempfile.NamedTemporaryFile(mode="w")
        self.fastq.write("@1\nACG\n+\nEEE\n@2\nTTTT\n+\nEEEE\n@3\nANNN\n+\nEEEE\n") 
        self.fastq.seek(0)

    def test_ids(self):
        ids = utils.parse_read_ids(self.fastq.name)
        self.assertEqual(len(ids), 3)
        self.assertTrue("1" in ids)
        self.assertTrue("2" in ids)
        self.assertTrue("3" in ids)


class TestIdConsistency(unittest.TestCase):
    def test_one_empty(self):
        self.assertFalse(utils.check_all_read_ids_are_consistent(set(), {'a'}))

    def test_different(self):
        self.assertFalse(utils.check_all_read_ids_are_consistent({'a', 'b'}, {'a'}))
        self.assertFalse(utils.check_all_read_ids_are_consistent({'a'}, {'a', 'b'}))

    def test_the_same(self):
        self.assertTrue(utils.check_all_read_ids_are_consistent({'a', 'b', 'c'}, {'c', 'a', 'b'}))


class TestAddToolSample(unittest.TestCase):
    def setUp(self):
        self.human_annotation = [("read1", 1), ("read2", 1), ("read3", 0)]
        self.wide_annotation = utils.add_tool_sample("tool", "sample", self.human_annotation)

    def test_number_or_rows(self):
        self.assertEqual(len(self.human_annotation), len(list(self.wide_annotation)))

    def test_number_or_columns(self):
        self.assertEqual(len(self.human_annotation[0]) + 2, len(list(self.wide_annotation)[0]))


if __name__ == "__main__":
    unittest.main()
