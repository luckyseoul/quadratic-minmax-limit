"""Independent checker and full node-result provenance binding."""
import hashlib
import json
from functools import lru_cache
from pathlib import Path

from scripts.joint_5p_minus_one_mesh_check import MODES


ROOT = Path(__file__).resolve().parents[1]
REPLAY = ROOT / 'evidence/joint_5p_minus_one_mesh_replay.json'
REPLAY_SHA256 = 'bbbdb1d5aab9c50a8e7a3f8e349147e32ba488b1894142d112a539183484105b'
EXPECTED = {'flat-symbolic': ('soulkiller', 'x86_64'),
            'p1-carries': ('NUKA', 'x86_64'),
            'p3-carries': ('jellyfin', 'x86_64'),
            'arm-floors': ('orin', 'aarch64')}


@lru_cache(maxsize=None)
def live_result(mode):
    return MODES[mode]()


def test_flat_row_identity_is_coefficientwise_and_not_a_catalog():
    row = live_result('flat-symbolic')
    assert all(value == {} for value in row['identity_residuals'].values())
    assert row['common_actual_P_range'] == list(range(10))
    assert row['forced_row_count'] == 5
    assert row['equality_catalog_used'] is False


def test_p1_carried_row_counts_and_new_target():
    rows = live_result('p1-carries')['rows']
    for row in rows:
        assert row['k'] == 5*row['p']-1
        assert [branch['forced_rows_at_least'] for branch in row['branches']] == [4, 8, 7, 7, 6, 6, 5]


def test_p3_sharp_carry_has_two_units_and_zero_quotients_are_explicit():
    for row in live_result('p3-carries')['rows']:
        assert len(row['sharp_branches']) == 4
        assert all(branch['quotient_carry'] == 2 and branch['forced_rows_at_least'] == 7 for branch in row['sharp_branches'])
        assert all(branch['forced_rows_at_least'] == 5 for branch in row['zero_quotient_branches'])
        assert row['p_minus_one_branch']['lift_mass'] == row['p']-1


def test_arm_floor_and_actual_P9_boundary():
    row = live_result('arm-floors')
    assert len(row['flat_rows']) == 60
    assert all(item['forced_rows_at_least'] == 5 for item in row['flat_rows'])
    at_nine = [item for item in row['flat_rows'] if item['actual_P'] == 9]
    assert all(item['least_allowed_Q'] == 0 and not item['forbidden_Q_is_legal'] for item in at_nine)
    assert [item['p'] for item in row['independent_floor_rows']] == [29, 31]


def test_mesh_sources_hosts_capabilities_and_complete_results_are_bound():
    raw = REPLAY.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == REPLAY_SHA256
    replay = json.loads(raw)
    script_hash = hashlib.sha256((ROOT / replay['script']).read_bytes()).hexdigest()
    assert script_hash == replay['script_sha256']
    assert len(replay['runs']) == 4
    assert {run['output']['mode'] for run in replay['runs']} == set(EXPECTED)
    for run in replay['runs']:
        output = run['output']
        assert isinstance(output, dict)
        assert run['exit_code'] == 0 and output['passed'] is True
        assert (output['host'], output['architecture']) == EXPECTED[output['mode']]
        assert output['script_sha256'] == script_hash
        assert output['result'] == json.loads(json.dumps(live_result(output['mode'])))
    assert 'ARM64' in replay['routing']['orin']
    assert replay['local_equality_catalog_used'] is False
