import pickle
import numpy as np
import re
import sqlite3
import cPickle as pickle
from scipy import interpolate
from os.path import expanduser
home = expanduser("~")

from Chandra.Time import DateTime
from Ska.engarchive import fetch_eng as fetch

import glimmondb as gdb

axafauto_url = 'http://occweb.cfa.harvard.edu/occweb/FOT/engineering/thermal/AXAFAUTO_RSYNC/'


def is_not_nan(arg):
    try:
        np.isnan(arg)
    except:  # Need to use blanket except, NotImplementedError won't catch
        return True
    return False


def open_sqlite_file():
    try:
        db = sqlite3.connect(axafauto_url + 'G_LIMMON_Archive/glimmondb.sqlite3')
    except:
        db = sqlite3.connect(home + '/AXAFAUTO/G_LIMMON_Archive/glimmondb.sqlite3')
    return db


def open_tdb_file():
    try:
        tdbs = pickle.load(open(axafauto_url + 'TDB_Archive/tdb_all.pkl', 'r'))
    except:
        tdbs = pickle.load(open(home + '/AXAFAUTO/TDB_Archive/tdb_all.pkl', 'r'))
    return tdbs


def get_tdb_limits(msid, dbver='p013', tdbs=None):
    """ Retrieve the TDB limits from a json version of the MS Access database.

    :param msid: String containing the mnemonic name, must correspond to a numeric limit set

    :returns safetylimits: Dictionary of numeric limits with keys: 'warning_low', 'caution_low',
        'caution_high', 'warning_high'


    Returns an empty dict object if there are no limits specified in the
    database
    """

    def assign_sets(dbsets):
        """ Copy over only the limit/expst sets, other stuff is not copied.

        This also adds a list of set numbers.
        """
        limits = {'setkeys': []}
        for setnum in dbsets.keys():
            setnumint = int(setnum) - 1
            limits.update({setnumint: dbsets[setnum]})
            limits['setkeys'].append(setnumint)
        return limits

    def get_tdb(dbver, tdbs):
        if tdbs:
            return tdbs[dbver.lower()]
        else:
            tdbs = open_tdb_file()
            return tdbs[dbver.lower()]

    msid = msid.lower().strip()

    try:
        tdb = get_tdb(dbver, tdbs)

        limits = assign_sets(tdb[msid]['limit'])
        limits['type'] = 'limit'

        if is_not_nan(tdb[msid]['limit_default_set_num']):
            limits['default'] = tdb[msid]['limit_default_set_num'] - 1
        else:
            limits['default'] = 0

        # Add limit switch info if present
        if is_not_nan(tdb[msid]['limit_switch_msid']):
            limits['mlimsw'] = tdb[msid]['limit_switch_msid']

        # Fill in switchstate info if present
        for setkey in limits['setkeys']:
            if 'state_code' in limits[setkey].keys():
                limits[setkey]['switchstate'] = limits[setkey]['state_code']
                _ = limits[setkey].pop('state_code')

        # For now, only the default limit set is returned, this will help with backwards compatibility.
        # Future versions, rewritten for web applications will not have this limitation.
        tdblimits = limits[limits['default']]

    except KeyError:
        print('{} does not have limits in the TDB'.format(msid.upper()))
        tdblimits = {}

    return tdblimits


def find_violation_time_spans(times, booldata):
    booldata = list(booldata)
    booldata.insert(0, False) # Prepend to make indices line up
    idata = np.array(booldata, dtype=type(1))
    d = np.diff(idata)

    starts = d == 1
    stops = d == -1

    starts = list(starts)
    stops = list(stops)

    if idata[-1] == 1:
        stops.insert(-1, True)
    else:
        stops.insert
    
    startinds = np.where(starts)[0]
    stopinds = np.where(stops)[0]
    starts = times[startinds]
    stops = times[stopinds]
    timebounds = zip(starts, stops)
    indexbounds = zip(startinds, stopinds)

    # activesets = [np.unique(activesetids[a:b]) for a,b in zip(startinds, stopinds)]


    return timebounds, indexbounds



#-------------------------------------------------------------------------------------------------
# Code for checking numeric limits
#-------------------------------------------------------------------------------------------------

def get_safety_limits(msid):
    """ Update the current database numeric limits

    :param msid: String containing the mnemonic name, must correspond to a numeric limit set

    :returns safetylimits: Dictionary of numeric limits with keys: 'warning_low', 'caution_low',
        'caution_high', 'warning_high'

    The database limits are replaced with the G_LIMMON limits if the G_LIMMON numeric limits are 
    more permissive. Only those individual limits that are more permissive are replaced, not the 
    entire set.

    The updated limit set is returned in a separate dictionary, regardless of whether or not any 
    updates were made.

    This is meant to generate a pseudo-tdb database limit set. This assumes that the GLIMMON limit
    set reflects the new database limits. However often GLIMMON limits will be set within newly 
    adjusted database limits, so the values returned by this routine may be more conservative than 
    necessary.


    """

    msid = msid.lower().strip()

    # Set the safetylimits dict here. An empty dict is returned if there are no
    # limits specified. This is intended and relied upon later.
    safetylimits = get_tdb_limits(msid)

    # Read the GLIMMON data
    try:
        db = open_sqlite_file()
        cursor = db.cursor()
        cursor.execute("""SELECT a.msid, a.setkey, a.default_set, a.warning_low, 
                          a.caution_low, a.caution_high, a.warning_high FROM limits AS a 
                          WHERE a.setkey = a.default_set AND a.msid = ?
                          AND a.modversion = (SELECT MAX(b.modversion) FROM limits AS b
                          WHERE a.msid = b.msid and a.setkey = b.setkey)""", [msid, ])
        lims = cursor.fetchone()
        glimits = {'warning_low': lims[3], 'caution_low': lims[4], 'caution_high': lims[5],
                   'warning_high': lims[6]}
    except:
        print('{} not in G_LIMMON Database, message generated in pylimmon.get_safety_limits()'
              .format(msid.upper()))
        glimits = {}

    # If there are no limits in the TDB but there are in GLIMMON, use the
    # GLIMMON limits
    if not safetylimits and glimits:
        safetylimits = glimits

    # If there are limits in both GLIMMON and the TDB use the set that
    # is most permissive.
    if safetylimits and glimits:

        if glimits['warning_low'] < safetylimits['warning_low']:
            safetylimits['warning_low'] = glimits['warning_low']
            print('Updated warning low safety limit for %s' % msid)

        if glimits['caution_low'] < safetylimits['caution_low']:
            safetylimits['caution_low'] = glimits['caution_low']
            print('Updated caution low safety limit for %s' % msid)

        if glimits['warning_high'] > safetylimits['warning_high']:
            safetylimits['warning_high'] = glimits['warning_high']
            print('Updated warning high safety limit for %s' % msid)

        if glimits['caution_high'] > safetylimits['caution_high']:
            safetylimits['caution_high'] = glimits['caution_high']
            print('Updated caution high safety limit for %s' % msid)

    return safetylimits


def get_mission_safety_limits(msid, tdbs=None):
    """
    this assumes that glimmon limits can indicate when a safety limit has been adjusted
    """
    def liminterp(tsum, times, limits):
        # Nans are filled in for cases where a limit isn't established until some point after launch
        # The last date for safety limits and for trending limits should be near the current time and
        #  be identical so that one doesn't dominate when it shouldn't.
        f = interpolate.interp1d(times, limits, kind='zero', bounds_error=False, fill_value=np.nan)
        return list(f(tsum))

    limdict = get_limits(msid)
    lastdate = np.max(limdict['limsets'][0]['times'])

    trendinglimits = {'msid': msid, 'warning_low': [], 'caution_low': [], 'caution_high': [],
                      'warning_high': [], 'times': []}

    # Assume the default set is always 0 - I know this is a hack, but it will work for now
    trendinglimits['warning_low'] = limdict['limsets'][0]['warning_low']
    trendinglimits['caution_low'] = limdict['limsets'][0]['caution_low']
    trendinglimits['caution_high'] = limdict['limsets'][0]['caution_high']
    trendinglimits['warning_high'] = limdict['limsets'][0]['warning_high']
    trendinglimits['times'] = limdict['limsets'][0]['times']

    if not tdbs:
        tdbs = open_tdb_file()
    tdbversions = gdb.get_tdb(return_dates=True)
    allsafetylimits = {'warning_low': [], 'caution_low': [], 'caution_high': [],
                       'warning_high': [], 'times': []}
    for ver in np.sort(tdbversions.keys()):
        date = tdbversions[ver]
        safetylimits = get_tdb_limits(msid, dbver=ver, tdbs=tdbs)
        if safetylimits:
            allsafetylimits['warning_low'].append(safetylimits['warning_low'])
            allsafetylimits['caution_low'].append(safetylimits['caution_low'])
            allsafetylimits['caution_high'].append(safetylimits['caution_high'])
            allsafetylimits['warning_high'].append(safetylimits['warning_high'])
            allsafetylimits['times'].append(DateTime(date).secs)

    if len(allsafetylimits['warning_low']) == 0:
        return None

    # Repeat the last limit to prevent NaNs from being entered for safety limits when interpolating
    for key in allsafetylimits.keys():
        allsafetylimits[key].append(allsafetylimits[key][-1])
    allsafetylimits['times'][-1] = lastdate

    tsum = np.sort(np.unique(np.concatenate((trendinglimits['times'], allsafetylimits['times']))))
    for kind in ['warning_low', 'caution_low', 'caution_high', 'warning_high']:
        trendinglimits[kind] = liminterp(tsum, trendinglimits['times'], trendinglimits[kind])
        allsafetylimits[kind] = liminterp(tsum, allsafetylimits['times'], allsafetylimits[kind])
        if 'high' in kind:
            allsafetylimits[kind] = [
                np.nanmax((t, a)) for t, a in zip(trendinglimits[kind], allsafetylimits[kind])]
            # allsafetylimits[kind] = list(np.nanmax((trendinglimits[kind], allsafetylimits[kind]), axis=0))
        else:
            allsafetylimits[kind] = [
                np.nanmin((t, a)) for t, a in zip(trendinglimits[kind], allsafetylimits[kind])]

    allsafetylimits['times'] = list(tsum)

    # Repeat the last limit for the current date to make plotting easier
    # for key in allsafetylimits.keys():
    #     allsafetylimits[key].append(allsafetylimits[key][-1])
    # allsafetylimits['times'][-1] = DateTime().secs

    return allsafetylimits


def get_limits(msid):

    db = open_sqlite_file()
    cursor = db.cursor()
    cursor.execute("""SELECT a.msid, a.setkey, a.datesec, a.mlmenable, a.default_set, a.switchstate, a.mlimsw, 
                      a.caution_high, a.caution_low, a.warning_high, a.warning_low 
                      FROM limits AS a WHERE a.msid=? """, [msid.lower(), ])
    current_limits = cursor.fetchall()
    db.close()

    limdict = {'msid': current_limits[0][0], 'limsets': {}}

    for row in current_limits:
        setnum = row[1]
        if setnum not in limdict['limsets'].keys():
            limdict['limsets'][row[1]] = {'switchstate': [], 'mlmenable': [], 'times': [],
                                          'caution_high': [], 'caution_low': [], 'warning_low': [],
                                          'warning_high': [], 'mlimsw': [], 'default_set': []}

        limdict['limsets'][setnum]['switchstate'].append(row[5])
        limdict['limsets'][setnum]['mlmenable'].append(row[3])
        limdict['limsets'][setnum]['times'].append(row[2])
        limdict['limsets'][setnum]['caution_high'].append(row[7])
        limdict['limsets'][setnum]['caution_low'].append(row[8])
        limdict['limsets'][setnum]['warning_high'].append(row[9])
        limdict['limsets'][setnum]['warning_low'].append(row[10])
        limdict['limsets'][setnum]['mlimsw'].append(row[6])
        limdict['limsets'][setnum]['default_set'].append(row[4])

    # Append data for current time + 24 hours to avoid interpolation errors
    #
    # You count on this being done in get_mission_safety_limits()
    for setnum in limdict['limsets'].keys():
        limdict['limsets'][setnum]['times'].append(DateTime().secs + 24 * 3600)
        limdict['limsets'][setnum]['warning_low'].append(
            limdict['limsets'][setnum]['warning_low'][-1])
        limdict['limsets'][setnum]['caution_low'].append(
            limdict['limsets'][setnum]['caution_low'][-1])
        limdict['limsets'][setnum]['caution_high'].append(
            limdict['limsets'][setnum]['caution_high'][-1])
        limdict['limsets'][setnum]['warning_high'].append(
            limdict['limsets'][setnum]['warning_high'][-1])
        limdict['limsets'][setnum]['mlmenable'].append(limdict['limsets'][setnum]['mlmenable'][-1])
        limdict['limsets'][setnum]['mlimsw'].append(limdict['limsets'][setnum]['mlimsw'][-1])
        limdict['limsets'][setnum]['default_set'].append(
            limdict['limsets'][setnum]['default_set'][-1])
        limdict['limsets'][setnum]['switchstate'].append(
            limdict['limsets'][setnum]['switchstate'][-1])

    return limdict


def check_limit_msid(msid, t1, t2, greta_msid=None):
    """ Check to see if temperatures are within expected numeric limits.

    :param msid: String containing the mnemonic name
    :param t1: String containing the start time in HOSC format (e.g. 2015:174:08:59:00.000)
    :param t2: String containing the stop time in HOSC format (e.g. 2015:174:15:59:30.000)

    :returns combined_sets_check: Dictionary of arrays indicating whether the value at a 
        particular time is within the defined limits (False) or outside the defined limits (True)
        for the following limit types: 'warning_low', 'caution_low', 'caution_high', 'warning_high'

    Violations are flagged as True. Time values are returned in the combined_sets_check dictionary.
    """

    def combine_limit_checks(all_sets_check):

        all_sets_check_keys = all_sets_check.keys()
        currentset = all_sets_check[all_sets_check_keys.pop(0)]
        wh = currentset['warning_high']
        whlim = currentset['warning_high_limit']

        ch = currentset['caution_high']
        chlim = currentset['caution_high_limit']

        cl = currentset['caution_low']
        cllim = currentset['caution_low_limit']

        wl = currentset['warning_low']
        wllim = currentset['warning_low_limit']

        # import readline # optional, will allow Up/Down/History in the console
        # import code
        # vars = globals().copy()
        # vars.update(locals())
        # shell = code.InteractiveConsole(vars)
        # shell.interact()


        setid = np.zeros(len(wh), dtype=np.int) - 1
        # Locations without NANs same for all limit types
        ind = ~np.isnan(whlim) 
        setid[ind] = 0

        for setnum in all_sets_check_keys:
            currentset = all_sets_check[setnum]
            ind = ~np.isnan(currentset['warning_high_limit'])
    
            wh = wh | currentset['warning_high']
            whlim[ind] = currentset['warning_high_limit'][ind]
    
            ch = ch | currentset['caution_high']
            chlim[ind] = currentset['caution_high_limit'][ind]

            cl = cl | currentset['caution_low']
            cllim[ind] = currentset['caution_low_limit'][ind]
        
            wl = wl | currentset['warning_low']
            wllim[ind] = currentset['warning_low_limit'][ind]

            setid[ind] = int(setnum)

        return {'warning_high_violation': wh, 'warning_high_limit': whlim,
                'caution_high_violation': ch, 'caution_high_limit': chlim,
                'caution_low_violation': cl, 'caution_low_limit': cllim, 
                'warning_low_violation': wl, 'warning_low_limit': wllim,
                'active_set_ids':setid}

    def check_limit_set(msid, limdict, setnum, data):
        mlimsws = limdict['limsets'][setnum]['mlimsw']
        switchstates = limdict['limsets'][setnum]['switchstate']
        times = limdict['limsets'][setnum]['times']
        defaults = limdict['limsets'][setnum]['default_set']

        # Mask identifies "Good" values, so start off with all as False (i.e Bad) and fill in 
        # True where appropriate.
        mask = np.array([False] * len(data.times))

        # [:-1] because the last limit definition is just a copy of the previous definition
        items = zip(times[:-1], times[1:], mlimsws[:-1], switchstates[:-1], defaults[:-1])
        for t1, t2, mlimsw, switchstate, default in items:
            ind1 = data.times >= t1
            ind2 = data.times < t2
            time_ind = ind1 & ind2

            if 'none' in mlimsw:
                if default == setnum:
                    mask = mask | time_ind
            else:
                mask_switch = data[mlimsw].vals == switchstate.upper()
                mask_switch = mask_switch & time_ind
                mask = mask_switch | mask

        check = {}
        check['warning_high'], check['warning_high_limit'] = check_limit(msid, limdict, setnum, data, mask, 'warning_high')
        check['caution_high'], check['caution_high_limit'] = check_limit(msid, limdict, setnum, data, mask, 'caution_high')
        check['caution_low'], check['caution_low_limit'] = check_limit(msid, limdict, setnum, data, mask, 'caution_low')
        check['warning_low'], check['warning_low_limit'] = check_limit(msid, limdict, setnum, data, mask, 'warning_low')

        return check

    def check_limit(msid, limdict, setnum, data, mask, limtype):
        limtype = limtype.lower()

        tlim = limdict['limsets'][setnum]['times']
        vlim = limdict['limsets'][setnum][limtype]
        enab = limdict['limsets'][setnum]['mlmenable']

        f = interpolate.interp1d(tlim, vlim, kind='zero', bounds_error=False, fill_value=np.nan)
        intlim = f(data.times)
        if 'high' in limtype:
            limcheck = data[msid].vals > intlim
        else:
            limcheck = data[msid].vals < intlim

        # Make sure violations are not reported when this set was disabled
        f = interpolate.interp1d(tlim, enab, kind='zero', bounds_error=False, fill_value=np.nan)
        enabled = f(data.times) == 1
        limcheck = limcheck & enabled

        # Make sure violations are not reported when this set is not active
        limcheck[~mask] = False

        # Flag durations when this set is not enabled or active with NANs.
        intlim[~enabled] = np.NAN
        intlim[~mask] = np.NAN

        return limcheck, intlim

    def get_min_violation_data(msid, spans):
        ''' Retrieve detailed data during violation time spans. '''
        violation_data = {'extrema':[], 'times':[]}
        for start, stop in spans:
            data = fetch.Msid(msid, start, stop, stat=None)
            extrema_ind = np.argmin(data.vals)
            violation_data['extrema'].append(data.vals[extrema_ind])
            violation_data['times'].append(data.times[extrema_ind])
        return violation_data

    def get_max_violation_data(msid, spans):
        ''' Retrieve detailed data during violation time spans. '''
        violation_data = {'extrema':[], 'times':[]}
        for start, stop in spans:
            data = fetch.Msid(msid, start, stop, stat=None)
            extrema_ind = np.argmax(data.vals)
            violation_data['extrema'].append(data.vals[extrema_ind])
            violation_data['times'].append(data.times[extrema_ind])
        return violation_data

    # MSID names should be in lower case
    msid = msid.lower()
    if not greta_msid:
        # If greta_msid is not defined, then they are the same msid
        greta_msid = msid
    else:
        greta_msid = greta_msid.lower()

    # Query limit information
    limdict = get_limits(greta_msid.lower())

    # Add limit switch msids to msid list
    mlimsw = np.unique([s for setnum in limdict['limsets'].keys()
                        for s in limdict['limsets'][setnum]['mlimsw']])
    mlimsw = list(mlimsw)
    if 'none' in mlimsw:
        mlimsw.remove('none')
    msids = [msid, ]
    if mlimsw:
        msids.extend(mlimsw)

    # Query data, interpolate to minimum time sampling or 0.256 seconds, whichever is larger
    data = fetch.Msidset(msids, t1, t2, stat=None)
    d = np.max(np.min([np.min(np.diff(data[m].times)) for m in msids]), 0.25620782)
    data.interpolate(dt=d)
    for mlimsw_msid in mlimsw:
        data[mlimsw_msid].vals = np.array([s.strip() for s in data[mlimsw_msid].vals])

    # Calculate violations for all limit types (caution high, etc.), for all sets.
    # Violations are only indicated where the set is valid as indicated by MLIMSW, if applicable.
    all_sets_check = {}
    for setnum in limdict['limsets'].keys():
        all_sets_check[setnum] = check_limit_set(msid, limdict, setnum, data)

    # Return boolean arrays for each limit type after compiling the results for each limit set.
    combined_sets_check = combine_limit_checks(all_sets_check)
    combined_sets_check['time'] = data.times

    # import readline # optional, will allow Up/Down/History in the console
    # import code
    # vars = globals().copy()
    # vars.update(locals())
    # shell = code.InteractiveConsole(vars)
    # shell.interact()


    # combined_sets_check has this format coming out of combine_limit_checks()
    # {'warning_high_violation': wh, 'warning_high_limit': whlim,
    #  'caution_high_violation': ch, 'caution_high_limit': chlim,
    #  'caution_low_violation': cl, 'caution_low_limit': cllim, 
    #  'warning_low_violation': wl, 'warning_low_limit': wllim}, setid

    # Produce the final dictionary to return only the time spans where violations occur
    violation_dict = {'any':False}

    # Note below that only the first value for active_set_ids and limits are returned. This could
    # be changed at some point in the future, however it is rare that a set id or limit will
    # change during a limit violation.
    if any(combined_sets_check['warning_low_violation']):
        timebounds, indexbounds = find_violation_time_spans(combined_sets_check['time'], 
                                                      combined_sets_check['warning_low_violation'])
        activesets = [combined_sets_check['active_set_ids'][a] for a,b in indexbounds]
        limitvals = [combined_sets_check['warning_low_limit'][a] for a,b in indexbounds]
        violation_dict.update({'warning_low': {'timespans':timebounds, 'activesets':activesets, 'limits':limitvals}})
        violation_data = get_min_violation_data(msid, timebounds)
        violation_dict['warning_low'].update(violation_data)
        violation_dict['any'] = True

    if any(combined_sets_check['caution_low_violation']):
        timebounds, indexbounds = find_violation_time_spans(combined_sets_check['time'], 
                                                      combined_sets_check['caution_low_violation'])
        activesets = [combined_sets_check['active_set_ids'][a] for a,b in indexbounds]
        limitvals = [combined_sets_check['caution_low_limit'][a] for a,b in indexbounds]
        violation_dict.update({'caution_low': {'timespans':timebounds, 'activesets':activesets, 'limits':limitvals}})
        violation_data = get_min_violation_data(msid,timebounds)
        violation_dict['caution_low'].update(violation_data)
        violation_dict['any'] = True

    if any(combined_sets_check['caution_high_violation']):
        timebounds, indexbounds = find_violation_time_spans(combined_sets_check['time'], 
                                                      combined_sets_check['caution_high_violation'])
        activesets = [combined_sets_check['active_set_ids'][a] for a,b in indexbounds]
        limitvals = [combined_sets_check['caution_high_limit'][a] for a,b in indexbounds]
        violation_dict.update({'caution_high': {'timespans':timebounds, 'activesets':activesets, 'limits':limitvals}})
        violation_data = get_max_violation_data(msid, timebounds)
        violation_dict['caution_high'].update(violation_data)
        violation_dict['any'] = True

    if any(combined_sets_check['warning_high_violation']):
        timebounds, indexbounds = find_violation_time_spans(combined_sets_check['time'], 
                                                      combined_sets_check['warning_high_violation'])
        activesets = [combined_sets_check['active_set_ids'][a] for a,b in indexbounds]
        limitvals = [combined_sets_check['warning_high_limit'][a] for a,b in indexbounds]
        violation_dict.update({'warning_high': {'timespans':timebounds, 'activesets':activesets, 'limits':limitvals}})
        violation_data = get_max_violation_data(msid, timebounds)
        violation_dict['warning_high'].update(violation_data)
        violation_dict['any'] = True

    violation_dict.update({'type':'limit'})

    return violation_dict

#-------------------------------------------------------------------------------------------------
# Code for checking expected states
#-------------------------------------------------------------------------------------------------

def get_states(msid):
    try:
        db.close()
    except:
        pass

    db = open_sqlite_file()
    cursor = db.cursor()
    cursor.execute("""SELECT a.msid, a.setkey, a.datesec, a.mlmenable, a.default_set, 
                          a.switchstate, a.mlimsw, a.expst FROM expected_states AS a WHERE a.msid=? """,
                   [msid.lower(), ])
    current_limits = cursor.fetchall()

    limdict = {'msid': current_limits[0][0], 'limsets': {}}

    for row in current_limits:
        setnum = row[1]
        if setnum not in limdict['limsets'].keys():
            limdict['limsets'][row[1]] = {'switchstate': [], 'mlmenable': [], 'times': [],
                                          'expst': [], 'mlimsw': [], 'default_set': []}

        limdict['limsets'][setnum]['switchstate'].append(row[5])
        limdict['limsets'][setnum]['mlmenable'].append(row[3])
        limdict['limsets'][setnum]['times'].append(row[2])
        limdict['limsets'][setnum]['expst'].append(row[7])
        limdict['limsets'][setnum]['mlimsw'].append(row[6])
        limdict['limsets'][setnum]['default_set'].append(row[4])

    # Append data for current time + 24 hours to avoid interpolation errors
    for setnum in limdict['limsets'].keys():
        limdict['limsets'][setnum]['times'].append(DateTime().secs + 24 * 3600)
        limdict['limsets'][setnum]['expst'].append(limdict['limsets'][setnum]['expst'][-1])
        limdict['limsets'][setnum]['mlmenable'].append(limdict['limsets'][setnum]['mlmenable'][-1])
        limdict['limsets'][setnum]['mlimsw'].append(limdict['limsets'][setnum]['mlimsw'][-1])
        limdict['limsets'][setnum]['default_set'].append(
            limdict['limsets'][setnum]['default_set'][-1])
        limdict['limsets'][setnum]['switchstate'].append(
            limdict['limsets'][setnum]['switchstate'][-1])

    return limdict


def check_state_msid(msid, t1, t2, greta_msid=None):
    """ Check to see if states match expected values.

    :param msid: String containing the mnemonic name
    :param t1: String containing the start time in HOSC format (e.g. 2015:174:08:59:00.000)
    :param t2: String containing the stop time in HOSC format (e.g. 2015:174:15:59:30.000)

    :returns combined_sets_check: Dictionary of arrays indicating whether the value at a 
        particular time violates the expected state (True) or does not (False)

    Violations are flagged as True. Time values are returned in the combined_sets_check dictionary.
    """

    def combine_state_checks(all_sets_check):

        all_sets_check_keys = all_sets_check.keys()
        currentset = all_sets_check[all_sets_check_keys.pop(0)]
        es = currentset['expst']

        eslim = currentset['expst_limit']

        setid = np.zeros(len(es)) - 1
        ind = np.array([True if len(s) > 0 else False for s in eslim])
        setid[ind] = 0

        for setnum in all_sets_check_keys:
            currentset = all_sets_check[setnum]
            ind = np.array([True if len(s) > 0 else False for s in eslim])

            es = es | currentset['expst']
            eslim[ind] = currentset['expst_limit'][ind]

            setid[ind] = setnum

        return {'expected_state_violation':es, 'expected_state':eslim, 'active_set_ids':setid}

    def check_state_set(msid, limdict, setnum, data):
        mlimsws = limdict['limsets'][setnum]['mlimsw']
        switchstates = limdict['limsets'][setnum]['switchstate']
        times = limdict['limsets'][setnum]['times']
        defaults = limdict['limsets'][setnum]['default_set']

        mask = np.array([False] * len(data.times))

        # [:-1] because the last limit definition is just a copy of the previous definition
        items = zip(times[:-1], times[1:], mlimsws[:-1], switchstates[:-1], defaults[:-1])
        for t1, t2, mlimsw, switchstate, default in items:
            ind1 = data.times >= t1
            ind2 = data.times < t2
            time_ind = ind1 & ind2

            if 'none' in mlimsw:
                if default == setnum:
                    mask = mask | time_ind
            else:
                mask_switch = data[mlimsw].vals == switchstate.upper()
                mask_switch = mask_switch & time_ind
                mask = mask_switch | mask

        check = {}
        check['expst'], check['expst_limit'] = check_state(msid, limdict, setnum, data, mask)

        return check

    def check_state(msid, limdict, setnum, data, mask):
        """ Check telemetry over time span for expected states.

        Since the history of expected state changes needs to be considered, the expected state
        for each telemetry point in time needs to be interpolated.
        """

        # Get the history of expected states
        tlim = limdict['limsets'][setnum]['times']
        vlim = limdict['limsets'][setnum]['expst']
        enab = limdict['limsets'][setnum]['mlmenable']

        # Determine the list of unique states in current expst list
        unique_states = np.unique(vlim)
        limids = np.arange(len(unique_states))

        # Generate a numeric representation of this expst history
        # This tells us what the expected states are at each time point
        vlim_numeric = np.zeros(len(vlim))
        for state, limid in zip(unique_states, limids):
            vlim_numeric[vlim == state] = limid

        # get history of expected states interpolated onto telemetry times
        f = interpolate.interp1d(
            tlim, vlim_numeric, kind='zero', bounds_error=False, fill_value=np.nan)
        # This is the list of numeric values representing EXPECTED states
        intlim_numeric = f(data.times)

        # Generate a numeric representation of the data, states not present in limdict are set to -1
        # This tells us what the ACTUAL states are at each time point
        vals_numeric = np.array([-1] * len(data.times))
        for state, limid in zip(unique_states, limids):
            vals_numeric[data[msid].vals == state] = limid

        # Generate boolean array where True marks where a violation occurs
        limcheck = vals_numeric != intlim_numeric

        # Make sure violations are not reported when this set was disabled
        f = interpolate.interp1d(tlim, enab, kind='zero', bounds_error=False, fill_value=np.nan)
        enabled = f(data.times) == 1
        limcheck = limcheck & enabled

        # "mask" tells us when this set is valid, make sure times when this set is not valid do not
        # report a violation
        limcheck[~mask] = False

        # Generate a list of the expected states in character form
        intlim_char = np.array([''] * len(data.times), dtype='S8')
        for state, limid in zip(unique_states, limids):
            intlim_char[intlim_numeric == limid] = state

        # Flag durations when this set is not enabled or active with empty strings
        intlim_char[~enabled] = ''
        intlim_char[~mask] = ''

        return limcheck, intlim_char


    def get_observed_violation_data(msid, spans):
        ''' Retrieve detailed data during violation time spans. '''
        violation_data = {'extrema':[]}
        for start, stop in spans:
            data = fetch.Msid(msid, start, stop, stat=None)
            extrema = np.array([s.lower() for s in np.unique(data.vals)])
            violation_data['extrema'].append(extrema)
        return violation_data




    # MSID names should be in lower case
    msid = msid.lower()
    if not greta_msid:
        # If greta_msid is not defined, then they are the same msid
        greta_msid = msid
    else:
        greta_msid = greta_msid.lower()

    # Query limit information
    limdict = get_states(greta_msid.lower())

    # Add limit switch msids to msid list
    mlimsw = np.unique([s for setnum in limdict['limsets'].keys()
                        for s in limdict['limsets'][setnum]['mlimsw']])
    mlimsw = list(mlimsw)
    if 'none' in mlimsw:
        mlimsw.remove('none')
    msids = [msid, ]
    if mlimsw:
        msids.extend(mlimsw)

    # Query data, interpolate to minimum time sampling or 0.256 seconds, whichever is larger
    data = fetch.Msidset(msids, t1, t2, stat=None)
    d = np.max(np.min([np.min(np.diff(data[m].times)) for m in msids]), 0.25620782)
    data.interpolate(dt=d)
    data[msid].vals = np.array([s.strip().lower() for s in data[msid].vals])
    for mlimsw_msid in mlimsw:
        data[mlimsw_msid].vals = np.array([s.strip() for s in data[mlimsw_msid].vals])

    # Calculate violations for all sets.
    # Violations are only indicated where the set is valid as indicated by MLIMSW, if applicable.
    all_sets_check = {}
    for setnum in limdict['limsets'].keys():
        all_sets_check[setnum] = check_state_set(limdict['msid'], limdict, setnum, data)


    # Compile the results for each set into one (time, boolean).
    combined_sets_check = combine_state_checks(all_sets_check)
    combined_sets_check['time'] = data.times

    # Produce the final dictionary to return only the time spans where violations occur
    violation_dict = {'any':False}

    # Note below that only the first value for active_set_ids and limits are returned. This could
    # be changed at some point in the future, however it is rare that a set id or limit will
    # change during a limit violation.
    if any(combined_sets_check['expected_state_violation']):
        timebounds, indexbounds = find_violation_time_spans(combined_sets_check['time'], 
                                                      combined_sets_check['expected_state_violation'])
        activesets = [combined_sets_check['active_set_ids'][a] for a,b in indexbounds]
        limitvals = [combined_sets_check['expected_state'][a] for a,b in indexbounds]
        violation_data = get_observed_violation_data(msid, timebounds)
        violation_dict.update(violation_data)
        violation_dict.update({'timespans':timebounds, 'activesets':activesets, 'limits':limitvals})
        violation_dict['any'] = True
    violation_dict.update({'type':'expst'})

    return violation_dict










