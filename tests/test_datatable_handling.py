from ebrowse.datatable_handling import lead_table_formatting
from collections import OrderedDict

def test_lead_table_formatting():
    sample_dict = OrderedDict(X=["asd", "hakjf"], A=["asdjk", "hajsdk"])
    values = {"draw": 12, "recordsTotal":67, "recordsFiltered": 45}
    ret = lead_table_formatting(sample_dict, values['draw'], values['recordsTotal'], values['recordsFiltered'])
    assert all([ret[k] == values[k] for k in values])
    for col, v in enumerate(sample_dict.values()):
        for row, rv in enumerate(v):
            assert ret['data'][row][col] == rv