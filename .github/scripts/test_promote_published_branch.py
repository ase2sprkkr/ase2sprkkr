"""Local, network-free policy tests for the PyPI branch promotion script."""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import unittest

from packaging.version import Version


spec = spec_from_file_location("release_promotion", Path(__file__).with_name("promote_published_branch.py"))
promotion = module_from_spec(spec)
spec.loader.exec_module(promotion)


class PromotionPolicyTests(unittest.TestCase):
    def test_beta_and_stable_go_to_distinct_branches(self):
        self.assertEqual(promotion.target_branch(Version("4.0.0-beta4")), "master")
        self.assertEqual(promotion.target_branch(Version("4.0.0")), "release")
        self.assertEqual(promotion.target_branch(Version("4.1.0rc1")), "master")

    def test_only_latest_published_in_same_channel(self):
        versions = promotion.published_versions({
            "3.6.0": [{"yanked": False}],
            "4.0.0b3": [{}],
            "4.0.0b4": [{}],
            "4.0.0b5": [{"yanked": True}],
            "4.0.0": [],
        })
        self.assertEqual(promotion.latest_in_channel(versions, Version("4.0.0b3")), Version("4.0.0b4"))
        self.assertEqual(promotion.latest_in_channel(versions, Version("3.6.0")), Version("3.6.0"))
        self.assertNotIn(Version("4.0.0b5"), versions)

    def test_can_advance_only_to_higher_version(self):
        self.assertTrue(promotion.can_advance(Version("3.6.0"), Version("4.0.0b4")))
        self.assertTrue(promotion.can_advance(Version("4.0.0b3"), Version("4.0.0b4")))
        self.assertFalse(promotion.can_advance(Version("4.0.0b5"), Version("4.0.0b4")))
        self.assertFalse(promotion.can_advance(Version("4.0.0b4"), Version("4.0.0b4")))

    def test_tag_must_match_metadata(self):
        metadata = '[project]\nname = "ase2sprkkr"\nversion = "4.0.0-beta4"\n'
        self.assertEqual(
            promotion.source_version_and_branch("v4.0.0-beta4", metadata),
            (Version("4.0.0b4"), "master"),
        )
        with self.assertRaises(ValueError):
            promotion.source_version_and_branch("v4.0.0-beta3", metadata)

    def test_other_project_refused(self):
        metadata = '[project]\nname = "other"\nversion = "4.0.0"\n'
        with self.assertRaises(ValueError):
            promotion.source_version_and_branch("v4.0.0", metadata)

    def test_unpublished_and_yanked_versions_ignored(self):
        versions = promotion.published_versions({
            "4.0.0b2": [],
            "4.0.0b3": [{"yanked": True}],
            "4.0.0b4": [{"yanked": True}, {"yanked": False}],
            "invalid-version-??": [{}],
        })
        self.assertEqual(versions, {Version("4.0.0b4")})


if __name__ == "__main__":
    unittest.main()
