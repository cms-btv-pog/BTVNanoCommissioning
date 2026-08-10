import unittest

import awkward as ak

from BTVNanoCommissioning.utils.AK4_parameters import correction_config
from BTVNanoCommissioning.utils.array_writer import (
    MET_COLLECTION_BY_CAMPAIGN,
    add_canonical_met,
    canonical_met_name,
)

EXPECTED_COLLECTIONS = {
    "Rereco17_94X": "MET",
    "2016preVFP-UL": "MET",
    "2016postVFP-UL": "MET",
    "2017-UL": "MET",
    "2018-UL": "MET",
    "Winter22Run3": "PuppiMET",
    "Summer22": "PuppiMET",
    "Summer22EE": "PuppiMET",
    "Summer23": "PuppiMET",
    "Summer23BPix": "PuppiMET",
    "Summer24": "PuppiMET",
    "Prompt25": "PuppiMET",
}
CAMPAIGN_YEARS = {
    "Rereco17_94X": "2017",
    "2016preVFP-UL": "2016",
    "2016postVFP-UL": "2016",
    "2017-UL": "2017",
    "2018-UL": "2018",
    "Winter22Run3": "2022",
    "Summer22": "2022",
    "Summer22EE": "2022",
    "Summer23": "2023",
    "Summer23BPix": "2023",
    "Summer24": "2024",
    "Prompt25": "2025",
}


class CanonicalMetTest(unittest.TestCase):
    def test_every_supported_fixed_campaign_is_classified(self):
        supported_fixed_campaigns = set(correction_config) - {"prompt_dataMC"}
        self.assertEqual(set(MET_COLLECTION_BY_CAMPAIGN), supported_fixed_campaigns)
        self.assertEqual(MET_COLLECTION_BY_CAMPAIGN, EXPECTED_COLLECTIONS)

    def test_exported_met_uses_the_classified_collection(self):
        for campaign, expected_collection in EXPECTED_COLLECTIONS.items():
            with self.subTest(campaign=campaign):
                events = ak.Array(
                    {
                        "PuppiMET": {"pt": [11.0], "phi": [0.11]},
                        "MET": {"pt": [22.0], "phi": [0.22]},
                    }
                )
                add_canonical_met(events, campaign, CAMPAIGN_YEARS[campaign])
                expected_pt = 11.0 if expected_collection == "PuppiMET" else 22.0
                expected_phi = 0.11 if expected_collection == "PuppiMET" else 0.22
                self.assertEqual(ak.to_list(events.MET_pt), [expected_pt])
                self.assertEqual(ak.to_list(events.MET_phi), [expected_phi])

    def test_prompt_campaign_uses_year_to_select_met(self):
        expected_by_year = {
            "2016": "MET",
            "2017": "MET",
            "2018": "MET",
            "2022": "PuppiMET",
            "2023": "PuppiMET",
            "2024": "PuppiMET",
            "2025": "PuppiMET",
        }
        for year, expected_collection in expected_by_year.items():
            with self.subTest(year=year):
                self.assertEqual(
                    canonical_met_name("prompt_dataMC", year),
                    expected_collection,
                )
                events = ak.Array(
                    {
                        "PuppiMET": {"pt": [11.0], "phi": [0.11]},
                        "MET": {"pt": [22.0], "phi": [0.22]},
                    }
                )
                add_canonical_met(events, "prompt_dataMC", year)
                expected_pt = 11.0 if expected_collection == "PuppiMET" else 22.0
                expected_phi = 0.11 if expected_collection == "PuppiMET" else 0.22
                self.assertEqual(ak.to_list(events.MET_pt), [expected_pt])
                self.assertEqual(ak.to_list(events.MET_phi), [expected_phi])

    def test_unknown_or_ambiguous_classification_raises(self):
        with self.assertRaisesRegex(ValueError, "Unknown MET campaign"):
            canonical_met_name("Summer26", "2026")
        with self.assertRaisesRegex(ValueError, "Cannot classify MET"):
            canonical_met_name("prompt_dataMC", None)
        with self.assertRaisesRegex(ValueError, "Cannot classify MET"):
            canonical_met_name("prompt_dataMC", "2021")


if __name__ == "__main__":
    unittest.main()
