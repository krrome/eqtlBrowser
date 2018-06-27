from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Float, Boolean
from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import validates
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.ext.declarative import ConcreteBase
from ebrowse import PATHS

Base = declarative_base()

transfer_keys_lead = ["beta", "chromosome", "geneSymbol", "pValue", "snpId", "start", "variantId", "fwdIter",
                         "probeId", "cellType", "dataset"]
transfer_keys_detailed = ["beta", "chromosome", "geneSymbol", "pValue", "snpId", "start", "variantId", "fwdIter",
                         "probeId", "cellType", "dataset", "ldR2"]

class LeadEqtlResult(Base):
    """
    "beta": [0.1,0.2,0.3],
    "chromosome": ["12", "12", "12"],
    "geneSymbol": ["APAF1", "IKBIP", "XYZ"],
    "pValue": [1.37476e-10, 1.37476e-10, 1.37476e-10], 
    "snpId": ["rs111604160", "rs111604160", "rs111604160"], 
    "start": [98980788, 98980788, 98980788],
    "variantId": ["12:98980788:A:G", "12:98980788:A:G", "12:98980788:A:G"],
    "fwdIter": [0,1,0],
    "ldR2": None / [.8, .9, .2],
    """
    __tablename__ = "lead_eqtl_result"
    id = Column(Integer, primary_key=True)
    beta = Column(Float)
    chromosome = Column(String(32))
    geneSymbol = Column(String(32), index=True)
    pValue = Column(Float)
    snpId = Column(String(32))
    start = Column(Integer)
    variantId = Column(String(32), index=True)
    fwdIter = Column(Integer)
    probeId = Column(String(32), index=True)
    cellType = Column(String(32), index=True)
    dataset = Column(String(32), index=True)
    detailedEqtlResults = relationship("DetailedEqtlResult", backref="lead_result")

    def as_dict(self):
        ret = {k: getattr(self, k) for k in transfer_keys_lead}
        return ret

    def __repr__(self):
        return "<LeadEqtlResult: {gsym} ({probe}): {var}>".format(gsym=self.geneSymbol,
                                                                  probe=self.probeId, var=self.variantId)


class DetailedEqtlResult(Base):
    __tablename__ = "detailed_eqtl_result"
    id = Column(Integer, primary_key=True)
    beta = Column(Float)
    chromosome = Column(String(32))
    geneSymbol = Column(String(32))
    pValue = Column(Float)
    snpId = Column(String(32))
    start = Column(Integer)
    variantId = Column(String(32))
    fwdIter = Column(Integer)
    ldR2 = Column(Float)
    probeId = Column(String(32))
    cellType = Column(String(32))
    dataset = Column(String(32))
    lead_result_id = Column(Integer, ForeignKey('lead_eqtl_result.id'))

    def as_dict(self):
        ret = {k: getattr(self, k) for k in transfer_keys_detailed}
        return ret


def get_db_uri():
    base_dir = PATHS['mysqldb_basepath']
    port = PATHS['mysqldb_port']
    socket = base_dir + "/socket"
    db_user = PATHS['db_user']
    db_pwd = PATHS['db_pwd']
    uri = "mysql://{dbuser}:{dbpwd}@localhost:{port}/eqtldata?unix_socket={socket}".format(dbuser=db_user, dbpwd=db_pwd,
                                                                            port=port, socket=socket)
    return uri

def _get_session():
    engine = create_engine(get_db_uri(), pool_recycle=3600)
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()
    return session

def append_to_dict(source, dest):
    for k in dest:
        if isinstance(source[k], list):
            dest[k].extend(source[k])
        else:
            dest[k].append(source[k])

session = _get_session()


def get_lead_by_region(chromosome, start, end, cell_type=None, dataset=None):
    LER = LeadEqtlResult
    query = session.query(LER).filter(LER.chromosome == chromosome, LER.start >= start, LER.start <= end)
    if cell_type is not None:
        query = query.filter(LER.cellType == cell_type)
    if dataset is not None:
        query = query.filter(LER.dataset == dataset)
    ret = {k:[] for k in transfer_keys_lead}
    for obj in query:
        append_to_dict(obj.as_dict(), ret)
    return ret


def get_detailed(cell_type, variant_id, probe_id, fwd_iter, dataset = None, by_id = False):
    LER = LeadEqtlResult
    DER = DetailedEqtlResult
    lead_query = session.query(LER).filter(LER.variantId == variant_id, LER.probeId == probe_id,
                                      LER.fwdIter == int(fwd_iter), LER.cellType == cell_type)
    if dataset is not None:
        lead_query = lead_query.filter(LER.dataset == dataset)
    ret = {k: [] for k in transfer_keys_detailed}
    for lead in lead_query:
        if by_id:
            detailed_query = session.query(DER).filter(DER.lead_result_id == lead.id)
        else:
            detailed_query = session.query(DER).filter(DER.lead_result == lead)
        for obj in detailed_query:
            append_to_dict(obj.as_dict(), ret)
    return ret

def get_leads(filter_equals=None, filter_condition=None, order_by=None, per_page=None, this_page=None):

    q = get_leads_query(filter_equals=filter_equals, filter_condition=filter_condition, order_by=order_by, per_page=per_page, this_page=this_page)

    ret = {k: [] for k in transfer_keys_lead}
    for obj in q:
        append_to_dict(obj.as_dict(), ret)
    return ret

def get_leads_query(filter_equals=None, filter_condition=None, order_by=None, per_page=None, this_page=None):
    LER = LeadEqtlResult
    q = session.query(LER)

    if filter_equals is not None:
        for k, v in filter_equals.items():
            q = q.filter(getattr(LER, k) == v)

    if filter_condition is not None:
        raise Exception("Not implemented!")

    if order_by is not None:
        for k, v in order_by.items():
            attr = getattr(LER, k)
            if v == "DESC":
                attr = attr.desc()
            q = q.order_by(attr)

    if per_page is not None:
        start_at = 0
        if this_page is not None:
            start_at = max(0, (this_page -1) * per_page)
        q = q.offset(start_at).limit(per_page)

    return q



"""
from ebrowse.database import get_detailed
import time
start_time = time.time()
for i in range(100):
    _=get_detailed("CD4", "1:229670964:G:A", "4810048", 0, by_id=False)

print("--- %s seconds ---" % (time.time() - start_time))

start_time = time.time()
for i in range(100):
    _=get_detailed("CD4", "1:229670964:G:A", "4810048", 0, by_id=True)

print("--- %s seconds ---" % (time.time() - start_time))

start_time = time.time()
for i in range(100):
    _=get_detailed("CD4", "1:229670964:G:A", "4810048", 0, by_id=False)

print("--- %s seconds ---" % (time.time() - start_time))

"""


