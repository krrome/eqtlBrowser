
def lead_table_formatting(ordered_obj, draw, total_records, filtered_records):
    ret = {"draw":draw, "recordsTotal": total_records, "recordsFiltered": filtered_records, 'data':[]}
    n_rows = len(ordered_obj[list(ordered_obj.keys())[0]])
    for r in range(n_rows):
        data_row = []
        for k in ordered_obj:
            data_row.append(ordered_obj[k][r])
        ret['data'].append(data_row)
    return ret

def lead_bootstrap_formatting(ordered_obj, total_records, filtered_records, reply_to_search = False):
    rows = []
    for vals in zip(*ordered_obj.values()):
        rows.append({k: vals[i] for i, k in enumerate(ordered_obj.keys())})

    ret = rows
    if reply_to_search:
        ret = {"rows": rows, "total": total_records}
    #handles = {k:ordered_obj[k]}
    #if reply_to_search
    return ret