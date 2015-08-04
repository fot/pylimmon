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

axafauto_url = 'http://occweb.cfa.harvard.edu/occweb/FOT/engineering/thermal/AXAFAUTO_RSYNC/'


def isnotnan(arg):
   try:
       np.isnan(arg)
   except: # Need to use blanket except, NotImplementedError won't catch
       return True
   return False


def opensqlitefile():
    try:
        db = sqlite3.connect(axafauto_url + 'G_LIMMON_Archive/glimmondb.sqlite3')
    except:
        db = sqlite3.connect(home + '/AXAFAUTO/G_LIMMON_Archive/glimmondb.sqlite3')
    return db


def opentdbfile():
    try:
        tdbs = pickle.load(open(axafauto_url + 'TDB_Archive/tdb_all.pkl','r'))
    except:
        tdbs = pickle.load(open(home + '/AXAFAUTO/TDB_Archive/tdb_all.pkl','r'))
    return tdbs


def getTDBLimits(msid, dbver='p013'):
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
        limits = {'setkeys':[]}
        for setnum in dbsets.keys():
            setnumint = int(setnum) - 1
            limits.update({setnumint:dbsets[setnum]})
            limits['setkeys'].append(setnumint)
        return limits

    def get_tdb(dbver):
        tdbs = opentdbfile()
        return tdbs[dbver.lower()]

    msid = msid.lower().strip()

    try:
        tdb = get_tdb(dbver)

        limits = assign_sets(tdb[msid]['limit'])
        limits['type'] = 'limit'

        if isnotnan(tdb[msid]['limit_default_set_num']):
            limits['default'] = tdb[msid]['limit_default_set_num'] - 1
        else:
            limits['default'] = 0

        # Add limit switch info if present
        if isnotnan(tdb[msid]['limit_switch_msid']):
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




def getSafetyLimits(msid):
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
    safetylimits = getTDBLimits(msid)

    # Read the GLIMMON data
    try:
        db = opensqlitefile()
        cursor = db.cursor()
        cursor.execute('''SELECT a.msid, a.setkey, a.default_set, a.warning_low, 
                          a.caution_low, a.caution_high, a.warning_high FROM limits AS a 
                          WHERE a.setkey = a.default_set AND a.msid = ?
                          AND a.modversion = (SELECT MAX(b.modversion) FROM limits AS b
                          WHERE a.msid = b.msid and a.setkey = b.setkey)''', [MSID.lower(),])
        lims = cursor.fetchone()
        glimits = {'warning_low':lims[3], 'caution_low':lims[4], 'caution_high':lims[5], 
                   'warning_high':lims[6]}
    except:
        print('{} not in G_LIMMON Database, message generated in gretafun.getSafetyLimits()'
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
            print('Updated warning low safety limit for %s'%msid)

        if glimits['caution_low'] < safetylimits['caution_low']:
            safetylimits['caution_low'] = glimits['caution_low']
            print('Updated caution low safety limit for %s'%msid)

        if glimits['warning_high'] > safetylimits['warning_high']:
            safetylimits['warning_high'] = glimits['warning_high']
            print('Updated warning high safety limit for %s'%msid)

        if glimits['caution_high'] > safetylimits['caution_high']:
            safetylimits['caution_high'] = glimits['caution_high']
            print('Updated caution high safety limit for %s'%msid)

    return safetylimits


#-------------------------------------------------------------------------------------------------
# Code for checking numeric limits
#-------------------------------------------------------------------------------------------------

def check_limit_msid(msid, t1, t2, greta_msid=None):
    ''' Check to see if temperatures are within expected numeric limits.

    :param msid: String containing the mnemonic name
    :param t1: String containing the start time in HOSC format (e.g. 2015:174:08:59:00.000)
    :param t2: String containing the stop time in HOSC format (e.g. 2015:174:15:59:30.000)

    :returns combined_sets_check: Dictionary of arrays indicating whether the value at a 
        particular time is within the defined limits (False) or outside the defined limits (True)
        for the following limit types: 'warning_low', 'caution_low', 'caution_high', 'warning_high'

    Violations are flagged as True. Time values are returned in the combined_sets_check dictionary.
    '''


    def getlimits(msid):

        db = opensqlitefile()
        cursor = db.cursor()
        cursor.execute('''SELECT a.msid, a.setkey, a.datesec, a.mlmenable, a.default_set, a.switchstate, a.mlimsw, 
                          a.caution_high, a.caution_low, a.warning_high, a.warning_low 
                          FROM limits AS a WHERE a.msid=? ''', [msid.lower(),])
        current_limits = cursor.fetchall()
        db.close()

        limdict = {'msid':current_limits[0][0], 'limsets':{}}

        for row in current_limits:
            setnum = row[1]
            if setnum not in limdict['limsets'].keys():
                limdict['limsets'][row[1]] = {'switchstate':[], 'mlmenable':[], 'times':[], 
                                              'caution_high':[], 'caution_low':[], 'warning_low':[], 
                                              'warning_high':[], 'mlimsw':[], 'default_set':[]}

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
        for setnum in limdict['limsets'].keys():
            limdict['limsets'][setnum]['times'].append(DateTime().secs + 24*3600)
            limdict['limsets'][setnum]['warning_low'].append(limdict['limsets'][setnum]['warning_low'][-1])
            limdict['limsets'][setnum]['caution_low'].append(limdict['limsets'][setnum]['caution_low'][-1])
            limdict['limsets'][setnum]['caution_high'].append(limdict['limsets'][setnum]['caution_high'][-1])
            limdict['limsets'][setnum]['warning_high'].append(limdict['limsets'][setnum]['warning_high'][-1])
            limdict['limsets'][setnum]['mlmenable'].append(limdict['limsets'][setnum]['mlmenable'][-1])
            limdict['limsets'][setnum]['mlimsw'].append(limdict['limsets'][setnum]['mlimsw'][-1])
            limdict['limsets'][setnum]['default_set'].append(limdict['limsets'][setnum]['default_set'][-1])
            limdict['limsets'][setnum]['switchstate'].append(limdict['limsets'][setnum]['switchstate'][-1])

        return limdict


    def combine_limit_checks(all_sets_check):
        
        all_sets_check_keys = all_sets_check.keys()
        currentset = all_sets_check[all_sets_check_keys.pop(0)]
        wh = currentset['warning_high']
        ch = currentset['caution_high']
        cl = currentset['caution_low']
        wl = currentset['warning_low']

        for setnum in all_sets_check_keys:
            currentset = all_sets_check[setnum]
            wh = wh | currentset['warning_high']
            ch = ch | currentset['caution_high']
            cl = cl | currentset['caution_low']
            wl = wl | currentset['warning_low']
        
        return {'warning_high_violation':wh, 'caution_high_violation':ch, 
                'caution_low_violation':cl, 'warning_low_violation':wl}

        
    def check_limit_set(msid, limdict, setnum, data):
        mlimsws = limdict['limsets'][setnum]['mlimsw']
        switchstates = limdict['limsets'][setnum]['switchstate']
        times = limdict['limsets'][setnum]['times']
        defaults = limdict['limsets'][setnum]['default_set']

        mask = np.array([False]*len(data.times))

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
        check['warning_high'] = check_limit(msid, limdict, setnum, data, mask, 'warning_high')
        check['caution_high'] = check_limit(msid, limdict, setnum, data, mask, 'caution_high')
        check['caution_low'] = check_limit(msid, limdict, setnum, data, mask, 'caution_low')
        check['warning_low'] = check_limit(msid, limdict, setnum, data, mask, 'warning_low')

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

        limcheck[~mask] = False

        return limcheck

    if not greta_msid:
        # If greta_msid is not defined, then they are the same msid
        greta_msid = msid

    limdict = getlimits(greta_msid.lower())
    mlimsw = np.unique([s for setnum in limdict['limsets'].keys()
                               for s in limdict['limsets'][setnum]['mlimsw']])
    mlimsw = list(mlimsw)
    if 'none' in mlimsw:
        mlimsw.remove('none')
    msids = [msid, ]
    if mlimsw:
        msids.extend(mlimsw)
        
    data = fetch.Msidset(msids, t1, t2, stat=None)
    data.interpolate()
    for mlimsw_msid in mlimsw:
        data[mlimsw_msid].vals = np.array([s.strip() for s in data[mlimsw_msid].vals])

    all_sets_check = {}
    for setnum in limdict['limsets'].keys():
        all_sets_check[setnum] = check_limit_set(msid, limdict, setnum, data)

    combined_sets_check = combine_limit_checks(all_sets_check)
    combined_sets_check['time'] = data.times

    return combined_sets_check
        

#-------------------------------------------------------------------------------------------------
# Code for checking expected states
#-------------------------------------------------------------------------------------------------

def check_state_msid(msid, t1, t2, greta_msid=None):
    ''' Check to see if states match expected values.

    :param msid: String containing the mnemonic name
    :param t1: String containing the start time in HOSC format (e.g. 2015:174:08:59:00.000)
    :param t2: String containing the stop time in HOSC format (e.g. 2015:174:15:59:30.000)

    :returns combined_sets_check: Dictionary of arrays indicating whether the value at a 
        particular time violates the expected state (True) or does not (False)

    Violations are flagged as True. Time values are returned in the combined_sets_check dictionary.
    '''


    def getstates(msid):
        try:
            db.close()
        except:
            pass

        db = opensqlitefile()
        cursor = db.cursor()
        cursor.execute('''SELECT a.msid, a.setkey, a.datesec, a.mlmenable, a.default_set, 
                          a.switchstate, a.mlimsw, a.expst FROM expected_states AS a WHERE a.msid=? ''', 
                       [msid.lower(),])
        current_limits = cursor.fetchall()

        limdict = {'msid':current_limits[0][0], 'limsets':{}}

        for row in current_limits:
            setnum = row[1]
            if setnum not in limdict['limsets'].keys():
                limdict['limsets'][row[1]] = {'switchstate':[], 'mlmenable':[], 'times':[], 
                                              'expst':[], 'mlimsw':[], 'default_set':[]}

            limdict['limsets'][setnum]['switchstate'].append(row[5])
            limdict['limsets'][setnum]['mlmenable'].append(row[3])
            limdict['limsets'][setnum]['times'].append(row[2])
            limdict['limsets'][setnum]['expst'].append(row[7])
            limdict['limsets'][setnum]['mlimsw'].append(row[6])
            limdict['limsets'][setnum]['default_set'].append(row[4])

        # Append data for current time + 24 hours to avoid interpolation errors
        for setnum in limdict['limsets'].keys():
            limdict['limsets'][setnum]['times'].append(DateTime().secs + 24*3600)
            limdict['limsets'][setnum]['expst'].append(limdict['limsets'][setnum]['expst'][-1])
            limdict['limsets'][setnum]['mlmenable'].append(limdict['limsets'][setnum]['mlmenable'][-1])
            limdict['limsets'][setnum]['mlimsw'].append(limdict['limsets'][setnum]['mlimsw'][-1])
            limdict['limsets'][setnum]['default_set'].append(limdict['limsets'][setnum]['default_set'][-1])
            limdict['limsets'][setnum]['switchstate'].append(limdict['limsets'][setnum]['switchstate'][-1])

        return limdict


    def combine_state_checks(all_sets_check):
        
        all_sets_check_keys = all_sets_check.keys()
        es = all_sets_check[all_sets_check_keys.pop(0)]
        
        for setnum in all_sets_check_keys:
            currentset = all_sets_check[setnum]
            es = es | currentset['expst']
        
        return es

        
    def check_state_set(msid, limdict, setnum, data):
        mlimsws = limdict['limsets'][setnum]['mlimsw']
        switchstates = limdict['limsets'][setnum]['switchstate']
        times = limdict['limsets'][setnum]['times']
        defaults = limdict['limsets'][setnum]['default_set']

        mask = np.array([False]*len(data.times))

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
            
        return check_state(msid, limdict, setnum, data, mask)

        
    def check_state(msid, limdict, setnum, data, mask):
        ''' Check telemetry over time span for expected states.
        
        Since the history of expected state changes needs to be considered, the expected state
        for each telemetry point in time needs to be interpolated.
        '''
        
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
        f = interpolate.interp1d(tlim, vlim_numeric, kind='zero', bounds_error=False, fill_value=np.nan)
        intlim_numeric = f(data.times)

        # Generate a numeric representation of the data, states not present in limdict are set to -1
        # This tells us what the actual states are at each time point
        vals_numeric = np.array([-1]*len(data.times))
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

        return limcheck


    if not greta_msid:
        # If greta_msid is not defined, then they are the same msid
        greta_msid = msid

    limdict = getstates(greta_msid.lower())
    mlimsw = np.unique([s for setnum in limdict['limsets'].keys()
                               for s in limdict['limsets'][setnum]['mlimsw']])

    mlimsw = list(mlimsw)
    if 'none' in mlimsw:
        mlimsw.remove('none')
    msids = [msid, ]
    if mlimsw:
        msids.extend(mlimsw)
    
    try:
        data = fetch.Msidset(msids, t1, t2, stat=None)
    except:
        print msids, t1, t2

    data.interpolate()
    data[msid].vals = np.array([s.strip().lower() for s in data[msid].vals])
    for mlimsw_msid in mlimsw:
        data[mlimsw_msid].vals = np.array([s.strip() for s in data[mlimsw_msid].vals])

    all_sets_check = {}
    for setnum in limdict['limsets'].keys():
        all_sets_check[setnum] = check_state_set(limdict['msid'], limdict, setnum, data)

    combined_sets_check = combine_state_checks(all_sets_check)
    
    return {'expected_state_violation':combined_sets_check, 'time':data.times}
        
    

