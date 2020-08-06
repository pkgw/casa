from taskinit import casalog

def is_rflag_report(item):
    """
    Is this an item from a flagdata report dictionary?

    :param item: an object, normally an item from a dictionary

    :returns: whether item looks like a report from Rflag (type and name = rflag).
    """
    return 'type' in item and item['type'] == 'rflag'\
        and 'name' in item and item['name'] == 'Rflag'

def combine_rflag_subreport(sub_dict, agg_dict):
    """ Produces an aggregated RFlag return dictionary by adding in a sub-report.

    You normally call this function on a sequence of per-subMS RFlag return dictionaries
    to aggregate all the (sub-)reports into an overall report. Then call 
    finalize_agg_rflag_thresholds() to calculate overall timedev/freqdev thresholds.
    The output from this function has the threshold vectors in a list-of-list-of-list
    format which needs to be finalized using finalize_agg_rflag_thresholds().

    Example RFlag return dictionary:
    {'freqdev': array([[1, 0, 3.12e-02], [1, 3, 2.19e-02], [1, 4, 2.42e-02]]),
    'type': 'rflag', 'name': 'Rflag', 'timedev':
    array([[1, 0, 7.09e-03], [1, 3, 5.43e-03], [1, 4, 7.83e-03]]) }

    :param sub_dict: (sub-)report/dictionary returned by RFlag (from one subMS)
    :param agg_dict: aggregated report or dictionary to aggregate 'sub_dict' into

    :returns: RFlag dictionary after aggregating sub_dict into agg_dict
    """
    for key, item in list(sub_dict.items()):
        agg_dict[key] = _aggregate_rflag_item(key, item, agg_dict)

    return agg_dict

def _aggregate_rflag_item(key, item, ret_dict):
    """
    Aggregates a key-item pair into ret_dict, both from RFlag return dictionaries.
    """

    def aggregate_rflag_thresholds(item, ret_item):
        """
        RFlag produces threshold vectors (freqdev or timedev vector) as a 2D numpy
        array with rows:
        [spw_id, field_id, value]
        Example:
        array([[1, 0, 3.12e-02], [1, 3, 2.19e-02], [1, 4, 2.42e-02]])
        In general there is a list of vectors like these for multiple spw_id-field_id pairs.

        This function aggregates such list of vectors produced for different subMS.
        In the aggregation stage, the structure used is a list-of-list-of-list:
        a list with one element for every spw-field pair, holding:
        [spw_id, field_id, [val1, val2, val3] where val1, val2, ... are the thresholds
        for different subMSs. A finalize step is needed to average/median the innermost
        values.
        Using this trick (accumulate threshold values into a list) which is far from ideal
        but I didn't find a more simple solution given the data structure used for the rflag
        reports (a list of dictionaries structured as a dictionary).

        :param item: an RFlag list of threshold vectors to aggregate
        :param ret_item: an RFlag threshold list-of-list-of-list to aggregate into

        :returns: The result of aggregating item into ret_item
        """
        import numpy as np

        def eq_id(row_a, row_b):
            return row_a[0] == row_b[0] and row_a[1] == row_b[1]

        if type(ret_item) is np.ndarray:
            ret_item = ret_item.tolist()
            # Init as list to add sub-reports
            for idx in range(len(ret_item)):
                ret_item[idx][2] = [ret_item[idx][2]]

        # Find a place for every row of the sub report to be added
        for idx_in in range(item.shape[0]):
            found_idx = False
            for ret_idx in range(len(ret_item)):
                if eq_id(item[idx_in], ret_item[ret_idx]):
                    found_idx = True
                    ret_item[ret_idx][2].append(item[idx_in, 2])
                    break
            if not found_idx:
                ret_item.append([item[idx_in, 0], item[idx_in, 1], [item[idx_in, 2]]])

        return ret_item

    if key in ret_dict:
        ret_item = ret_dict[key]
        if not isinstance(ret_item, str):
            # must be either 'freqdev' or 'timedev'
            ret_dict[key] = aggregate_rflag_thresholds(item, ret_item)
        else:
            ret_dict[key] = item

    return ret_dict[key]


def finalize_agg_rflag_thresholds(rflag_dict):
    """
    For the thresholds included in an RFlag return dictionary (timedev and freqdev):
    build a 2D numpy array from a list of lists of lists, calculating a median of 
    thresholds throughout sub-MSs

    :param rflag_dict: RFlag dictionary with the un-finalized list-of-list-of-list
    structure produced by combine_rflag_subreport().

    :returns: the dictionary finalized, that is, with the per-subMS thresholds
    combined, currently using the median of the subMS values.
    """

    def spw_field_comp(x, y):
        """
        Comparator function to sort by (spw_id, field_id) pairs from the first and
        second coords of RFlag threshold vectors (example):
        [1, 0, 3.12e-02] < [1, 3, 2.19e-02]
        [1, 2, 3.12e-02] < [2, 0, 2.19e-02]
        """
        if x[0] < y[0] or (x[0] == y[0] and x[1] < y[1]):
            return  -1
        elif x[0] > y[0] or (x[0] == y[0] and x[1] > y[1]):
            return 1
        else:
            return 0

    import numpy as np

    for key, val in list(rflag_dict.items()):
        if not isinstance(val, str):
            # If the list was empty, we need a dummy (0,3)-shaped array
            if 0 == len(val):
                rflag_dict[key] = np.empty(shape=[0,3])
                continue

            # Choosing median for now. This is an open question from CAS-10202.
            for idx in range(len(val)):
                val[idx] = [val[idx][0], val[idx][1], np.median(val[idx][2])]
            # Sort to match better what is produced when not using parallelization
            val = sorted(val, cmp=spw_field_comp)
            rflag_dict[key] = np.array(val)

    return rflag_dict
