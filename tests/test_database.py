from ebrowse.database import get_leads_query, session, LeadEqtlResult, string_keys_lead
import numpy as np

def assert_dict_equality(a,b):
    assert set(a.keys()) == set(b.keys())
    for k in a:
        assert a[k] == b[k]

def test_get_leads_query():
    q = get_leads_query(per_page=5)
    assert q.count() == 5
    r = [el.as_dict() for el in q.all()]

    # Test pagination
    q2_1 = get_leads_query(per_page=2)
    assert q2_1.count() == 2
    r2_1 = [el.as_dict() for el in q2_1.all()]
    q2_1_cp = get_leads_query(per_page=2, this_page=1)
    assert q2_1_cp.count() == 2
    r2_1_cp = [el.as_dict() for el in q2_1_cp.all()]
    q2_2 = get_leads_query(per_page=2, this_page=2)
    assert q2_2.count() == 2
    r2_2 = [el.as_dict() for el in q2_2.all()]
    [assert_dict_equality(el1, el2) for el1, el2 in zip(r[0:2], r2_1)]
    [assert_dict_equality(el1, el2) for el1, el2 in zip(r2_1, r2_1_cp)]
    [assert_dict_equality(el1, el2) for el1, el2 in zip(r[2:4], r2_2)]

    # Test filter equals
    query_gene_symbol = q.first().geneSymbol
    q_eq = get_leads_query(filter_equals={"geneSymbol": query_gene_symbol})
    r_eq = q_eq.all()
    assert all([el.geneSymbol == query_gene_symbol for el in r_eq])

    # Test filter_all
    query_strs = ["ABCC", "2230", "6:T:C", "D4"]
    for query_str in query_strs:
        q_abbc = get_leads_query(filter_equals_any=query_str)
        print(q_abbc.count())
        assert all([any([query_str in getattr(el, k) for k in string_keys_lead]) for el in q_abbc])


    # Test ordering
    from sqlalchemy import func, desc
    col_oi = LeadEqtlResult.probeId
    max_probe, mp_count= session.query(col_oi, func.count(col_oi).label('counts')).group_by(col_oi).order_by(desc('counts')).first()
    q_asc = get_leads_query(filter_equals={"probeId": max_probe}, order_by={"pValue": "ASC"})
    q_desc = get_leads_query(filter_equals={"probeId": max_probe}, order_by={"pValue": "DESC"})
    assert mp_count == q_asc.count()
    assert mp_count == q_desc.count()
    pvs_asc = np.array([el.pValue for el in q_asc])
    pvs_desc = np.array([el.pValue for el in q_desc])
    assert (pvs_asc[np.argsort(pvs_asc)] == pvs_asc).all()
    assert (pvs_desc[np.argsort(pvs_desc)[::-1]] == pvs_desc).all()

    # Test get counts
    rv = get_leads_query(per_page=5, add_counts=True)
    assert isinstance(rv, tuple)
    assert len(rv) == 3
    assert isinstance(rv[1], int)
    assert isinstance(rv[2], int)
