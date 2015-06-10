import json
import os
import shutil
import tempfile
import unittest

from decontamlib import main

data_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data")


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

