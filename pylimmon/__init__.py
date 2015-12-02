from pylimmon import open_sqlite_file, open_tdb_file, get_tdb_limits, get_safety_limits, DBDIR, TDBDIR
from pylimmon import check_limit_msid, check_state_msid, get_limits, get_states, get_mission_safety_limits
from .version import __version__

print('Using G_LIMMON DB Here:{}'.format(DBDIR))
print('Using TDB Here:{}'.format(TDBDIR))