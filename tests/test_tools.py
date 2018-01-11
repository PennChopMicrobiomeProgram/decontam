import unittest
try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

from decontamlib.tools import (
    _FilteringTool, None_human,
    )


class FilteringToolTests(unittest.TestCase):
    def test_get_args(self):
        self.assertEqual(_FilteringTool.get_argnames(), ["index"])
        self.assertEqual(None_human.get_argnames(), [])


if __name__ == "__main__":
    unittest.main()
