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
        self.output_fps = {
            "nonhuman": ("B5_short_R1.fastq", "B5_short_R2.fastq"),
            "human": ("B5_short_R1_human.fastq", "B5_short_R2_human.fastq"),
            }
        for k, (v1, v2) in self.output_fps.items():
            self.output_fps[k] = (
                os.path.join(self.output_dir, v1),
                os.path.join(self.output_dir, v2),
                )

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

        for fp in self.output_fps["nonhuman"]:
            self.assertFalse(os.path.exists(fp))

        for fp in self.output_fps["human"]:
            with open(fp) as f:
                num_lines = len(list(f))
            self.assertEqual(num_lines, 10 * 4)


    def test_no_human(self):
        config_file = tempfile.NamedTemporaryFile(suffix=".json")
        json.dump({"method": "no_human"}, config_file)
        config_file.seek(0)
        self.args.extend(["--config-file", config_file.name])

        main.human_filter_main(self.args)

        with open(self.summary_fp) as f:
            obs_summary = json.load(f)
        self.assertEqual(obs_summary["data"], {"false": 10})

        for fp in self.output_fps["nonhuman"]:
            with open(fp) as f:
                num_lines = len(list(f))
            self.assertEqual(num_lines, 10 * 4)

        for fp in self.output_fps["human"]:
            self.assertFalse(os.path.exists(fp))


    def test_bowtie(self):
        index_fp = os.path.join(data_dir, "fakehuman")
        config_file = tempfile.NamedTemporaryFile(suffix=".json")
        json.dump({"method": "bowtie2", "index": index_fp}, config_file)
        config_file.seek(0)
        self.args.extend(["--config-file", config_file.name])

        main.human_filter_main(self.args)

        for fp in self.output_fps["human"]:
            self.assertTrue(os.path.exists(fp))

        for fp in self.output_fps["nonhuman"]:
            self.assertTrue(os.path.exists(fp))


    def test_bwa(self):
        index_fp = os.path.join(data_dir, "fakehuman")
        config_file = tempfile.NamedTemporaryFile(suffix=".json")
        json.dump({"method": "bwa", "index": index_fp}, config_file)
        config_file.seek(0)
        self.args.extend(["--config-file", config_file.name])

        main.human_filter_main(self.args)

        for fp in self.output_fps["human"]:
            self.assertTrue(os.path.exists(fp))

        for fp in self.output_fps["nonhuman"]:
            self.assertTrue(os.path.exists(fp))


if __name__ == "__main__":
    unittest.main()

