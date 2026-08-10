import subprocess
import tempfile
import unittest
from unittest.mock import patch

from condor_lxplus import dashboard


class DashboardSubmissionRetryTest(unittest.TestCase):
    def test_pre_submission_failure_remains_retryable(self):
        with tempfile.TemporaryDirectory() as tmp:
            job_dir = f"{tmp}/jobs_test"
            dashboard.ensure_status(job_dir)

            class FailedBeforeSubmission:
                calls = 0

                def __init__(self, *args, **kwargs):
                    type(self).calls += 1
                    self.pid = 12345
                    self.returncode = 1

                def communicate(self):
                    return "", "schedd temporarily unavailable"

            with patch.object(
                subprocess, "Popen", FailedBeforeSubmission
            ), patch.object(dashboard.time, "sleep"):
                result = dashboard._run_automation_command(
                    job_dir,
                    ["simulated-submission"],
                    "checkoutputs",
                    expect_condor_submit=True,
                    max_retries=2,
                )

            self.assertEqual(FailedBeforeSubmission.calls, 2)
            self.assertEqual(result["automation"]["step"], "blocked")
            self.assertIsNone(result["automation"]["unreconciled_submission"])

    def test_successful_submission_is_not_retried_when_status_recording_fails(self):
        cases = (
            (
                "checkoutputs",
                "record_resubmitted_jobs",
                9001,
                2,
                ("unused", ["1", "2"]),
            ),
            ("hadd_submit", "record_hadd_submitted", 9002, 3, ("unused", 3)),
            ("hadd_check", "record_hadd_submitted", 9003, 1, ("unused", 1)),
        )
        for action, recorder_name, cluster_id, submitted_jobs, recorder_args in cases:
            with self.subTest(action=action), tempfile.TemporaryDirectory() as tmp:
                job_dir = f"{tmp}/jobs_test"
                dashboard.ensure_status(job_dir)

                class SuccessfulSubmitThenRecordFailure:
                    calls = 0

                    def __init__(self, *args, **kwargs):
                        type(self).calls += 1
                        self.pid = 12345
                        self.returncode = None

                    def communicate(self):
                        try:
                            getattr(dashboard, recorder_name)(
                                *recorder_args, cluster_id=cluster_id
                            )
                        except RuntimeError as exc:
                            self.returncode = 1
                            return (
                                f"{submitted_jobs} job(s) submitted to cluster "
                                f"{cluster_id}.\n",
                                str(exc),
                            )
                        raise AssertionError("recorder unexpectedly succeeded")

                with patch.object(
                    dashboard,
                    recorder_name,
                    side_effect=RuntimeError("simulated status persistence failure"),
                ) as recorder_mock, patch.object(
                    subprocess, "Popen", SuccessfulSubmitThenRecordFailure
                ):
                    result = dashboard._run_automation_command(
                        job_dir,
                        ["simulated-submission", action],
                        action,
                        expect_condor_submit=True,
                        max_retries=3,
                    )

                self.assertEqual(SuccessfulSubmitThenRecordFailure.calls, 1)
                self.assertEqual(recorder_mock.call_count, 1)
                automation = result["automation"]
                self.assertEqual(automation["step"], "blocked")
                self.assertEqual(
                    automation["blocked_state"],
                    "submitted_but_status_update_failed",
                )
                reconciliation = automation["unreconciled_submission"]
                self.assertEqual(reconciliation["action"], action)
                self.assertEqual(reconciliation["cluster_id"], cluster_id)
                self.assertEqual(reconciliation["submitted_jobs"], submitted_jobs)

                unchanged = dashboard.advance_job_automation(
                    job_dir, condor_status=None, hadd_condor_status=None
                )
                self.assertEqual(
                    unchanged["automation"]["blocked_state"],
                    "submitted_but_status_update_failed",
                )
                self.assertEqual(SuccessfulSubmitThenRecordFailure.calls, 1)


if __name__ == "__main__":
    unittest.main()
