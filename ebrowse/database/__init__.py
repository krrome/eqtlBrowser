from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Float, Boolean
from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy import ForeignKey, or_
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session
from sqlalchemy.orm import validates
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.ext.declarative import ConcreteBase
from ebrowse import PATHS
from contextlib import contextmanager

_use_mysql = False

if _use_mysql:
    Base = declarative_base()
else:
    Base = object



transfer_keys_lead = ["beta", "chromosome", "geneSymbol", "pValue", "snpId", "start", "variantId", "fwdIter",
                         "probeId", "cellType", "dataset"]

string_keys_lead = ["chromosome", "geneSymbol", "snpId", "variantId", "probeId", "cellType", "dataset"]

transfer_keys_detailed = ["beta", "chromosome", "geneSymbol", "pValue", "snpId", "start", "variantId", "fwdIter",
                         "probeId", "cellType", "dataset", "ldR2"]

transfer_keys_colocalisation = ["cellType","trait","probeId","chi2Eqtl","chi2Gwas","comparedInFinemap","indexGroupEqtl",
                                "indexGroupGwas","pvaluesEqtl","pvaluesGwas","variantId","snpProbEqtl","snpProbGwas"]

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
    chromosome = Column(String(32), index=True)
    geneSymbol = Column(String(32), index=True)
    pValue = Column(Float)
    snpId = Column(String(32), index=True)
    start = Column(Integer)
    variantId = Column(String(32), index=True)
    fwdIter = Column(Integer)
    ldR2 = Column(Float)
    probeId = Column(String(32), index=True)
    cellType = Column(String(32), index=True)
    dataset = Column(String(32), index=True)
    lead_result_id = Column(Integer, ForeignKey('lead_eqtl_result.id'), index=True)

    def as_dict(self):
        ret = {k: getattr(self, k) for k in transfer_keys_detailed}
        return ret



class ColocalisationResult(Base):
    __tablename__ = "colocalisation_result"
    id = Column(Integer, primary_key=True)
    cellType = Column(String(32), index=True)
    trait = Column(String(32), index=True)
    probeId = Column(String(32), index=True)
    chi2Eqtl = Column(Float)
    chi2Gwas = Column(Float)
    comparedInFinemap = Column(Boolean)
    indexGroupEqtl = Column(Integer)
    indexGroupGwas = Column(Integer)
    pvaluesEqtl = Column(Float)
    pvaluesGwas = Column(Float)
    variantId = Column(String(32), index=True)
    snpProbEqtl = Column(Float)
    snpProbGwas = Column(Float)

    def as_dict(self):
        ret = {k: getattr(self, k) for k in transfer_keys_colocalisation}
        return ret
    """      chi2_eqtl   chi2_gwas  compared_in_finemap  index_eqtl  index_gwas  \
2793  82.001408  676.727125                 True         522         522

           pos  pvalues_eqtl   pvalues_gwas             snp  snp_log10bf_eqtl  \
2793  56849749  1.359897e-19  3.400000e-149  3:56849749_T_C            9.0006

      snp_log10bf_gwas  snp_prob_eqtl  snp_prob_gwas
2793           12.8503            1.0            1.0

"""


def get_db_uri():
    base_dir = PATHS['mysqldb_basepath']
    port = PATHS['mysqldb_port']
    socket = base_dir + "/socket"
    db_user = PATHS['db_user']
    db_pwd = PATHS['db_pwd']
    uri = "mysql://{dbuser}:{dbpwd}@localhost:{port}/eqtldata?unix_socket={socket}".format(dbuser=db_user, dbpwd=db_pwd,
                                                                            port=port, socket=socket)
    return uri

if _use_mysql:
    engine = create_engine(get_db_uri(), pool_recycle=3600)
    Base.metadata.create_all(engine)
    session_factory = sessionmaker(bind=engine)


@contextmanager
def use_session():
    """Temporarily change the directory
    """
    Session = scoped_session(session_factory)
    try:
        yield Session()
    finally:
        Session.remove()


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




def get_lead_by_region(chromosome, start, end, cell_type=None, dataset=None):
    LER = LeadEqtlResult
    with use_session() as session:
        query = session.query(LER)
        if (chromosome is not None) and (start is not None) and (end is not None):
            query = query.filter(LER.chromosome == chromosome, LER.start >= start, LER.start <= end)
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
    with use_session() as session:
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


def get_leads(filter_equals=None, filter_condition=None, filter_equals_any=None, order_by=None,
                    per_page=None, this_page=None, offset = None, add_counts = True):
    LER = LeadEqtlResult
    ret = {}
    with use_session() as session:
        q = session.query(LER)
        all_entries=None
        num_filtered = None

        if add_counts:
            all_entries = q.count()

        if filter_equals is not None:
            for k, v in filter_equals.items():
                q = q.filter(getattr(LER, k) == v)

        if filter_condition is not None:
            raise Exception("Not implemented!")

        if filter_equals_any is not None:
            # copied from: https://stackoverflow.com/questions/4926757/sqlalchemy-query-where-a-column-contains-a-substring
            if '*' in filter_equals_any or '_' in filter_equals_any:
                filter_equals_any = filter_equals_any.replace('_', '__') \
                    .replace('*', '%') \
                    .replace('?', '_')
            else:
                filter_equals_any = '%{0}%'.format(filter_equals_any)

            filts = [getattr(LER, k).ilike(filter_equals_any) for k in string_keys_lead]
            q = q.filter(or_(*filts))

        if order_by is not None:
            for k, v in order_by.items():
                attr = getattr(LER, k)
                if v == "DESC":
                    attr = attr.desc()
                q = q.order_by(attr)

        if add_counts:
            num_filtered = q.count()

        if per_page is not None:
            start_at = 0
            if this_page is not None:
                start_at = max(0, (this_page -1) * per_page)
            elif offset is not None:
                start_at = offset
            q = q.offset(start_at).limit(per_page)

        ret = {k: [] for k in transfer_keys_lead}
        for obj in q:
            append_to_dict(obj.as_dict(), ret)

    if add_counts:
        return ret, all_entries, num_filtered
    else:
        return ret



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


