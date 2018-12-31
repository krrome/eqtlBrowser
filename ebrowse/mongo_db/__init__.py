import mongoengine
from ebrowse import PATHS
from mongoengine import Document, FloatField, IntField, StringField, EmbeddedDocumentListField, EmbeddedDocument, BooleanField, BinaryField




transfer_keys_lead = ["id", "beta", "betaFM", "chromosome", "geneSymbol", "pValue", "pValueFM", "snpId", "start", "variantId", "fwdIter",
                         "probeId", "cellType", "dataset", "fdr", "fdrFM", "maf", "ensembleGId"]

transfer_keys_detailed = ["beta", "chromosome", "geneSymbol", "pValue", "snpId", "start", "variantId", "fwdIter",
                         "probeId", "cellType", "dataset", "ldR2"]

transfer_keys_colocalisation = ["cellType","trait","probeId","chi2Eqtl","chi2Gwas","comparedInFinemap","indexGroupEqtl",
                                "indexGroupGwas","pvaluesEqtl","pvaluesGwas","variantId","snpProbEqtl","snpProbGwas"]


mongoengine.connect('eqtl_data', host=PATHS['mongo_host'], port=27017)

TOTAL_LEAD_RECORDS = None


def ler_primary_key(in_dict):
    return "_".join(["ler", in_dict["dataset"], in_dict["cellType"], in_dict["probeId"], in_dict["variantId"]])

class Colocalisation(EmbeddedDocument):
    """
    cellType                        CD14
    chi2Eqtl                     238.561
    chi2Gwas                     38.0912
    comparedInFinemap              False
    indexGroupEqtl                   737
    indexGroupGwas                   737
    probeId                      7150753
    pvaluesEqtl              8.09895e-54
    pvaluesGwas                  6.8e-10
    snpProbEqtl                   0.9714
    snpProbGwas                   0.9987
    trait                    myeloid_wbc
    variantId            8:142328719_G_T
    """
    cellType = StringField()
    chi2Eqtl = FloatField()
    chi2Gwas = FloatField()
    comparedInFinemap = BooleanField()
    indexGroupEqtl = IntField()
    indexGroupGwas = IntField()
    probeId = StringField()
    pvaluesEqtl = FloatField()
    pvaluesGwas = FloatField()
    snpProbEqtl = FloatField()
    snpProbGwas = FloatField()
    trait = StringField()
    variantId = StringField()



class LeadEqtlResult(Document):
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
    id = StringField(primary_key=True)# generate a unique ID for this entry
    leadEqtlExpressionId = StringField()
    beta = FloatField()
    betaFM = FloatField()
    chromosome = StringField(max_length=32)
    geneSymbol = StringField(max_length=32)# index
    pValue = FloatField()
    pValueFM = FloatField()
    snpId = StringField(max_length=32)
    start = IntField()
    variantId = StringField(max_length=32)# index
    fwdIter = IntField()
    probeId = StringField(max_length=32)# index
    cellType = StringField(max_length=32)# index
    dataset = StringField(max_length=32)# index
    fdr = FloatField()
    fdrFM = FloatField()
    maf = FloatField()
    ensembleGId = StringField(max_length=32)
    colocalisation = EmbeddedDocumentListField(Colocalisation)
    meta = {'indexes': ['geneSymbol', 'variantId', 'probeId', 'cellType', 'dataset']}

    def as_dict(self):
        ret = {k: getattr(self, k) for k in transfer_keys_lead}
        return ret

    def __repr__(self):
        return "<LeadEqtlResult: {gsym} ({probe}): {var}>".format(gsym=self.geneSymbol,
                                                                  probe=self.probeId, var=self.variantId)

class LeadEqtlExpression(Document):
    id = StringField(primary_key=True)  # generate a unique ID for this entry
    leadEqtlResultId = StringField()
    data = BinaryField()

class TestedProbes(Document):
    probeId = StringField(max_length=32)
    geneSymbol = StringField(max_length=200)
    cellType = StringField(max_length=32)

    def as_dict(self):
        ret = {k: getattr(self, k) for k in self._fields.keys() if k != 'id'}
        return ret

class ProbeExpressionRank(Document):
    cellType = StringField(max_length=32)
    data = BinaryField()


class DetailedEqtlResult(Document):
    beta = FloatField()
    chromosome = StringField(max_length=32)
    geneSymbol = StringField(max_length=32)
    pValue = FloatField()
    snpId = StringField(max_length=32)# index
    start = IntField()
    variantId = StringField(max_length=32)# index
    fwdIter = IntField()
    ldR2 = FloatField()
    probeId = StringField(max_length=32)# index
    cellType = StringField(max_length=32)
    dataset = StringField(max_length=32)
    leadEqtlResultId = StringField() #index
    meta = {'indexes': ['leadEqtlResultId']}

    def as_dict(self):
        ret = {k: getattr(self, k) for k in transfer_keys_detailed}
        return ret


def append_to_dict(source, dest):
    for k in dest:
        if isinstance(source[k], list):
            dest[k].extend(source[k])
        else:
            dest[k].append(source[k])

def get_total_read_records():
    global TOTAL_LEAD_RECORDS
    if TOTAL_LEAD_RECORDS is None:
        TOTAL_LEAD_RECORDS = LeadEqtlResult.objects().count()
    return TOTAL_LEAD_RECORDS


def query_document(DocClass, filter_dict, exact_match):
    str_to_mongo_expr = {">=": "gte", "<=": "lte", "<":"lt", ">":"gt"}
    query = {}
    if filter_dict is not None:
        field_dict = DocClass._fields
        for field_name, field in field_dict.items():
            if field_name in filter_dict:
                if isinstance(field, StringField):
                    if exact_match:
                        query[field_name] = filter_dict[field_name]
                    else:
                        query[field_name+ "__icontains"] = filter_dict[field_name]
                elif isinstance(field, IntField) or isinstance(field, FloatField):
                    prefix = None
                    if filter_dict[field_name].startswith(">=") or filter_dict[field_name].startswith("<="):
                        prefix = filter_dict[field_name][:2]
                        filter_dict[field_name] = filter_dict[field_name][2:].lstrip(" ")
                    elif filter_dict[field_name].startswith(">") or filter_dict[field_name].startswith("<"):
                        prefix = filter_dict[field_name][:1]
                        filter_dict[field_name] = filter_dict[field_name][1:].lstrip(" ")
                    field_name_fmt = field_name
                    if prefix is not None:
                        field_name_fmt += "__" + str_to_mongo_expr[prefix]
                    if len(filter_dict[field_name]) == 0:
                        continue
                    if isinstance(field, IntField):
                        query[field_name_fmt] = int(filter_dict[field_name])
                    else:
                        query[field_name_fmt] = float(filter_dict[field_name])
                elif isinstance(field, BooleanField):
                    query[field_name] = bool(filter_dict[field_name])
                else:
                    pass
    return query

def execute_table_query(DocClass, filter_dict, exact_match, limit = None, offset = None, sort_by = None,
                        sort_order = None):
    ret = None
    query = query_document(DocClass, filter_dict, exact_match)
    objs = DocClass.objects(**query)
    obj_count = objs.count()
    if sort_by is not None and sort_by in list(DocClass._fields.keys()):
        sort_prefix = ""
        if sort_order.lower() in ["desc"]:
            sort_prefix = "-"
        objs = objs.order_by(sort_prefix + sort_by)
    if limit is not None:
        limit = int(limit)
        if offset is not None:
            offset = int(offset)
        else:
            offset = 0
        objs = objs[offset:(offset + limit)]
    for i, obj in enumerate(objs):
        if i == 0:
            ret = {k: [] for k in obj.as_dict()}
        append_to_dict(obj.as_dict(), ret)
    return ret, obj_count


def get_leads(exact_match = True, sort = None, order = None, **kwargs):
    limit, offset = None, None
    if "limit" in kwargs:
        limit = kwargs['limit']
    if "offset" in kwargs:
        offset = kwargs['offset']
    ret, obj_count = execute_table_query(LeadEqtlResult, kwargs['filterDict'], exact_match, limit, offset,
                                         sort_by = sort, sort_order = order)
    return ret, obj_count, 0

def get_tested_probes(exact_match = True, sort = None, order = None, **kwargs):
    limit, offset = None, None
    if "limit" in kwargs:
        limit = kwargs['limit']
    if "offset" in kwargs:
        offset = kwargs['offset']
    ret, obj_count = execute_table_query(TestedProbes, kwargs['filterDict'], exact_match, limit, offset,
                                         sort_by=sort, sort_order=order)
    return ret, obj_count, 0

def get_detailed(cell_type, variant_id, probe_id, fwd_iter, leadEqtlResultId, dataset = None):
    LER = LeadEqtlResult
    DER = DetailedEqtlResult
    ret = {k: [] for k in transfer_keys_detailed}
    if leadEqtlResultId == "":
        query = dict(variantId = variant_id, probeId = probe_id, fwdIter = int(fwd_iter), cellType = cell_type)
        if dataset is not None:
            query['dataset'] = dataset
        lead_query = LER.objects(**query)
        lead_ret = None
        for i, lead in enumerate(lead_query):
            if i ==0:
                lead_ret = lead.as_dict()
            detailed_query = DER.objects(leadEqtlResultId = lead.id)
            for obj in detailed_query:
                append_to_dict(obj.as_dict(), ret)
    else:
        lead_query = LER.objects(id = leadEqtlResultId)
        lead_ret = lead_query[0].as_dict()
        detailed_query = DER.objects(leadEqtlResultId = leadEqtlResultId)
        for obj in detailed_query:
            append_to_dict(obj.as_dict(), ret)
    return ret, lead_ret

def get_lead_by_region(chromosome, start, end, cell_type=None, dataset=None):
    LER = LeadEqtlResult
    query = {}
    if (chromosome is not None) and (start is not None) and (end is not None):
        query.update(dict(chromosome = chromosome, start__gte = start, start__lte = end))
    if cell_type is not None:
        query.update(dict(cellType=cell_type))
    if dataset is not None:
        query.update(dict(dataset=dataset))
    ret = {k:[] for k in transfer_keys_lead}
    for obj in LER.objects(**query):
        append_to_dict(obj.as_dict(), ret)
    return ret


def get_lead_expression(**kwargs):
    import zlib
    import json
    query = query_document(LeadEqtlResult, kwargs, exact_match=True)
    objs = LeadEqtlResult.objects(**query)
    if not objs.count() == 1:
        return None
    lead_obj = objs[0]
    res = LeadEqtlExpression.objects(id=lead_obj.leadEqtlExpressionId)
    if not res.count() == 1:
        return None
    res_dict = json.loads(zlib.decompress(res[0].data).decode("utf8"))
    res_dict['assoc_info'] = lead_obj.as_dict()
    return res_dict

def get_expression_hist(cellType):
    import zlib
    import json
    objs = ProbeExpressionRank.objects(cellType = cellType)
    res_dict = json.loads(zlib.decompress(objs[0].data).decode("utf8"))
    res_dict['cellType'] = cellType
    return res_dict