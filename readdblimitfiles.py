"""
Read in each telemetry database used for limit monitoring, store all data in all serializable
datastructure. Save this file to .pkl and .json files.

Each database must have been converted to sets of CSV files, one for each file. Filenames and
file contents should all be in lower case text.

At this point only the tables related to MSID definitions, limit and expected state definitions,
and calibaration definitions are included.

The resulting datastructure will have the following format:

    all_databases[database_version][msid]['msid']
    all_databases[database_version][msid]['technical_name']
    all_databases[database_version][msid]['data_type']
    all_databases[database_version][msid]['calibration_type']
    all_databases[database_version][msid]['eng_unit']
    all_databases[database_version][msid]['low_raw_count']
    all_databases[database_version][msid]['high_raw_count']
    all_databases[database_version][msid]['total_length']
    all_databases[database_version][msid]['prop']
    all_databases[database_version][msid]['counter_msid']
    all_databases[database_version][msid]['range_msid']
    all_databases[database_version][msid]['calibration_switch_msid']
    all_databases[database_version][msid]['calibration_default_set_num']
    all_databases[database_version][msid]['limit_switch_msid']
    all_databases[database_version][msid]['limit_default_set_num']
    all_databases[database_version][msid]['es_switch_msid']
    all_databases[database_version][msid]['es_default_set_num']
    all_databases[database_version][msid]['owner_id']
    all_databases[database_version][msid]['description']
    all_databases[database_version][msid]['ehs_header_flag']
    all_databases[database_version][msid]['limit_lrvt_location']
    all_databases[database_version][msid]['em_error_description']

    all_databases[database_version][msid]['limit'][set_num]['caution_low']
    all_databases[database_version][msid]['limit'][set_num]['caution_high']
    all_databases[database_version][msid]['limit'][set_num]['warning_low']
    all_databases[database_version][msid]['limit'][set_num]['warning_high']
    all_databases[database_version][msid]['limit'][set_num]['delta']
    all_databases[database_version][msid]['limit'][set_num]['toler']
    all_databases[database_version][msid]['limit'][set_num]['em_all_samp_flag']

    all_databases[database_version][msid]['lim_switch'][set_num]['limit_set_num']
    all_databases[database_version][msid]['lim_switch'][set_num]['low_range']
    all_databases[database_version][msid]['lim_switch'][set_num]['high_range']
    all_databases[database_version][msid]['lim_switch'][set_num]['state_code']

    all_databases[database_version][msid]['cal_switch']['calibration_set_num']
    all_databases[database_version][msid]['cal_switch']['low_range']
    all_databases[database_version][msid]['cal_switch']['high_range']
    all_databases[database_version][msid]['cal_switch']['state_code']

    all_databases[database_version][msid]['point_pair'][set_num][sequence_num]['calibration_set_num']
    all_databases[database_version][msid]['point_pair'][set_num][sequence_num]['sequence_num']
    all_databases[database_version][msid]['point_pair'][set_num][sequence_num]['raw_count']
    all_databases[database_version][msid]['point_pair'][set_num][sequence_num]['eng_unit_value']

    all_databases[database_version][msid]['poly_cal'][set_num]['calibration_set_num']
    all_databases[database_version][msid]['poly_cal'][set_num]['end_unit_low']
    all_databases[database_version][msid]['poly_cal'][set_num]['eng_unit_high']
    all_databases[database_version][msid]['poly_cal'][set_num]['deg']
    all_databases[database_version][msid]['poly_cal'][set_num]['coef0']
    all_databases[database_version][msid]['poly_cal'][set_num]['coef1']
    all_databases[database_version][msid]['poly_cal'][set_num]['coef2']
    all_databases[database_version][msid]['poly_cal'][set_num]['coef3']
    all_databases[database_version][msid]['poly_cal'][set_num]['coef4']
    all_databases[database_version][msid]['poly_cal'][set_num]['coef5']
    all_databases[database_version][msid]['poly_cal'][set_num]['coef6']
    all_databases[database_version][msid]['poly_cal'][set_num]['coef7']
    all_databases[database_version][msid]['poly_cal'][set_num]['coef8']
    all_databases[database_version][msid]['poly_cal'][set_num]['coef9']

    all_databases[database_version][msid]['exp_state']['es_set_num']
    all_databases[database_version][msid]['exp_state']['expected_state']
    all_databases[database_version][msid]['exp_state']['toler']
    all_databases[database_version][msid]['exp_state']['em_all_samp_flag']

    all_databases[database_version][msid]['state_code'][set_num][sequence_num]['calibration_set_num']
    all_databases[database_version][msid]['state_code'][set_num][sequence_num]['sequence_num']
    all_databases[database_version][msid]['state_code'][set_num][sequence_num]['low_raw_count']
    all_databases[database_version][msid]['state_code'][set_num][sequence_num]['high_raw_count']
    all_databases[database_version][msid]['state_code'][set_num][sequence_num]['state_code']

    all_databases[database_version][msid]['es_switch']['es_set_num']
    all_databases[database_version][msid]['es_switch']['low_range']
    all_databases[database_version][msid]['es_switch']['high_range']
    all_databases[database_version][msid]['es_switch']['state_code']
"""

import os
import pandas
import cPickle as pickle
import json


def assignsetvals(db, table, field, sequence=False):
    """Convert TDB table to dictionary.

    Modifies a pre-existing dictionary by adding a TDB table in the desired format.

    :param db: Pre-existing dictionary to update
    :param table: TDB table as a Pandas 2D dataframe
    :param field: Desired name for converted table
    :param sequence: Boolean indicating if a sequence of items exists in `table` such as a set of
                     point pair calibration values

    """
    for row in table.values:
        msid = row[0]
        setnum = int(row[1])

        if not db[msid].has_key(field):
            db[msid].update({field:{}})

        if not sequence:
            db[msid][field][setnum] = dict(zip(table.columns[2:], row[2:]))
        else:
            seq = int(row[2])
            if not db[msid][field].has_key(setnum):
                db[msid][field].update({setnum:{}})

            db[msid][field][setnum][seq] = dict(zip(table.columns[3:], row[3:]))


def readdb(rootdir):
    """Read TDB csv file and return Pandas data frame.

    Modifies a pre-existing dictionary by adding a TDB table in the desired format.

    :param rootdir: String containing the location of the set of TDB csv files

    :returns: Dictionary of TDB tables in Pandas dataframe format 

    """
    tdbframes = {}
    tdbframes['tdbmsid'] = pandas.read_csv(os.path.join(rootdir, 'tdb_msid.csv'))
    tdbframes['tdblimit'] = pandas.read_csv(os.path.join(rootdir, 'tdb_limit.csv'))
    tdbframes['tdblimswitch'] = pandas.read_csv(os.path.join(rootdir, 'tdb_lim_switch.csv'))
    tdbframes['tdbpointpair'] = pandas.read_csv(os.path.join(rootdir, 'tdb_point_pair.csv'))
    tdbframes['tdbpolycal'] = pandas.read_csv(os.path.join(rootdir, 'tdb_poly_cal.csv'))
    tdbframes['tdbcalswitch'] = pandas.read_csv(os.path.join(rootdir, 'tdb_cal_switch.csv'))
    tdbframes['tdbexpstate'] = pandas.read_csv(os.path.join(rootdir, 'tdb_exp_state.csv'))
    tdbframes['tdbesswitch'] = pandas.read_csv(os.path.join(rootdir, 'tdb_es_switch.csv'))
    tdbframes['tdbstatecode'] = pandas.read_csv(os.path.join(rootdir, 'tdb_state_code.csv'))
    return tdbframes


def processdb(tdbframes):
    """Convert a TDB from Pandas dataframe format to a dictionary

    :param tdbframes: TDB in Pandas dataframe format

    :returns: TDB in dictionary format (serializable)

    """
    tdb = {}
    for row in tdbframes['tdbmsid'].values:
        tdb.update({row[0]:dict(zip(tdbframes['tdbmsid'].columns[1:], row[1:]))})
    assignsetvals(tdb, tdbframes['tdblimit'], 'limit')
    assignsetvals(tdb, tdbframes['tdblimswitch'], 'lim_switch')
    assignsetvals(tdb, tdbframes['tdbpointpair'], 'point_pair', sequence=True)
    assignsetvals(tdb, tdbframes['tdbpolycal'], 'poly_cal')
    assignsetvals(tdb, tdbframes['tdbcalswitch'], 'cal_switch')
    assignsetvals(tdb, tdbframes['tdbexpstate'], 'exp_state')
    assignsetvals(tdb, tdbframes['tdbesswitch'], 'es_switch')
    assignsetvals(tdb, tdbframes['tdbstatecode'], 'state_code', sequence=True)
    return tdb


def process_files(rootdir):
    """Return dictionary of all TDB's: P007, P009, P010, P011

    :param rootdir: String containing the location of all TDB directories

    :returns: Dictionary of serializable TDB's 

    """
    tdbframes = readdb(os.path.join(rootdir, 'p007'))
    tdb007 = processdb(tdbframes)

    tdbframes = readdb(os.path.join(rootdir, 'p009'))
    tdb009 = processdb(tdbframes)

    tdbframes = readdb(os.path.join(rootdir, 'p010'))
    tdb010 = processdb(tdbframes)

    tdbframes = readdb(os.path.join(rootdir, 'p011'))
    tdb011 = processdb(tdbframes)

    tdbframes = readdb(os.path.join(rootdir, 'p012'))
    tdb012 = processdb(tdbframes)

    tdbframes = readdb(os.path.join(rootdir, 'p013'))
    tdb013 = processdb(tdbframes)

    tdbframes = readdb(os.path.join(rootdir, 'p014'))
    tdb014 = processdb(tdbframes)

    return {'p007':tdb007, 'p009':tdb009, 'p010':tdb010, 'p011':tdb011, 'p012':tdb012, 
            'p013':tdb013, 'p014':tdb014}


if __name__ == '__main__':
    tdb_all = process_files('./')
    pickle.dump(tdb_all, open('tdb_all.pkl','w'), protocol=2)
    json.dump(tdb_all, open('tdb_all.json','w'))

