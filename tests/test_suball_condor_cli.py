import subprocess
import sys
import unittest

BASE_COMMAND = [
    sys.executable,
    "scripts/suball.py",
    "--scheme",
    "Validation",
    "--campaign",
    "Summer22",
    "--year",
    "2022",
    "--DAS_campaign",
    "test",
]


class SuballCondorCliTest(unittest.TestCase):
    def test_condor_only_argument_requires_explicit_condor_flag(self):
        completed = subprocess.run(
            BASE_COMMAND + ["--jobqueue", "workday"],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Condor-only option(s) require --condor", completed.stderr)

    def test_condor_requires_output_base(self):
        completed = subprocess.run(
            BASE_COMMAND + ["--condor"],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("--condor requires --condorOutputBase", completed.stderr)

    def test_help_documents_explicit_condor_mode(self):
        completed = subprocess.run(
            [sys.executable, "scripts/suball.py", "--help"],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0)
        normalized_help = " ".join(completed.stdout.split())
        self.assertIn("--condor", normalized_help)
        self.assertIn("Requires --condorOutputBase", normalized_help)


if __name__ == "__main__":
    unittest.main()
