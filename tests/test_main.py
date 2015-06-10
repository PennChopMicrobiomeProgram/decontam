import json
import os
import shutil
import tempfile
import unittest

from decontamlib import main

data_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data")

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
        tool_names = [ "bmtagger", "all_human", "no_human", "random_human", "bowtie" ]
        tool_runners = main.create_tools(tool_names, self.params)
        self.assertEqual(len(tool_runners), len(tool_names))
        self.assertTrue(isinstance(tool_runners[0], main.Bmtagger))
        self.assertTrue(isinstance(tool_runners[1], main.All_human))
        self.assertTrue(isinstance(tool_runners[2], main.None_human))
        self.assertTrue(isinstance(tool_runners[3], main.Random_human))
        self.assertTrue(isinstance(tool_runners[4], main.Bowtie))


class HumanFilterMainTests(unittest.TestCase):
    def setUp(self):
        self.fwd_fp = os.path.join(data_dir, "B5_short_R1.fastq")
        self.rev_fp = os.path.join(data_dir, "B5_short_R2.fastq")

        self.temp_dir = tempfile.mkdtemp()
        self.output_dir = os.path.join(self.temp_dir, "output")
        self.summary_fp = os.path.join(self.temp_dir, "summary.json")
        self.args = [
            "--forward-reads", self.fwd_fp,
            "--reverse-reads", self.rev_fp,
            "--output-dir", self.output_dir,
            "--summary-file", self.summary_fp,
            ]

        self.nonhuman_fwd_fp = os.path.join(
            self.output_dir, "B5_short_R1.fastq")
        self.nonhuman_rev_fp = os.path.join(
            self.output_dir, "B5_short_R2.fastq")
        self.human_fwd_fp = os.path.join(
            self.output_dir, "B5_short_R1_human.fastq")
        self.human_rev_fp = os.path.join(
            self.output_dir, "B5_short_R2_human.fastq")

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_all_human(self):
        config_file = tempfile.NamedTemporaryFile(suffix=".json")
        json.dump({"method": "all_human"}, config_file)
        config_file.seek(0)
        self.args.extend(["--config-file", config_file.name])

        main.human_filter_main(self.args)

        with open(self.summary_fp) as f:
            obs_summary = json.load(f)
        self.assertEqual(obs_summary["data"], {"true": 10})

        self.assertFalse(os.path.exists(self.nonhuman_fwd_fp))
        self.assertFalse(os.path.exists(self.nonhuman_rev_fp))

        with open(self.human_fwd_fp) as f:
            fwd_lines = len(list(f))
        self.assertTrue(fwd_lines, 10 * 4)

        with open(self.human_rev_fp) as f:
            rev_lines = len(list(f))
        self.assertTrue(rev_lines, 10 * 4)

    def test_no_human(self):
        config_file = tempfile.NamedTemporaryFile(suffix=".json")
        json.dump({"method": "no_human"}, config_file)
        config_file.seek(0)
        self.args.extend(["--config-file", config_file.name])

        main.human_filter_main(self.args)

        with open(self.summary_fp) as f:
            obs_summary = json.load(f)
        self.assertEqual(obs_summary["data"], {"false": 10})

        self.assertFalse(os.path.exists(self.human_fwd_fp))
        self.assertFalse(os.path.exists(self.human_rev_fp))

        with open(self.nonhuman_fwd_fp) as f:
            fwd_lines = len(list(f))
        self.assertTrue(fwd_lines, 10 * 4)

        with open(self.nonhuman_rev_fp) as f:
            rev_lines = len(list(f))
        self.assertTrue(rev_lines, 10 * 4)


if __name__ == "__main__":
    unittest.main()

