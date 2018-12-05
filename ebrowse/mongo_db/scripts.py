from ebrowse.database.utils import get_all_results_iter, get_colocalised_hits, get_tested_probes
from ebrowse.mongo_db import LeadEqtlResult, LeadEqtlExpression, DetailedEqtlResult, ler_primary_key, Colocalisation, TestedProbes, ProbeExpressionRank
import mongoengine
import pandas as pd
from tqdm import tqdm




class FinemapData:
    @staticmethod
    def _get_ler_pk(coloc_entry):
        in_dict = coloc_entry.to_dict()
        in_dict["dataset"]="merged"
        return ler_primary_key(in_dict)
    #
    def __init__(self):
        self.keys = []
        self.values = []
        self.all_keys = []
        for coloc in get_colocalised_hits():
            coloc_df = pd.DataFrame(coloc)
            self.values.append(coloc_df.apply(lambda x: x.to_dict(), axis=1).tolist())
            self.keys.append(coloc_df.apply(self._get_ler_pk, axis=1).tolist())
            self.all_keys.extend(self.keys[-1])
    #
    def query(self, lead_id):
        id = lead_id
        if id in self.all_keys:
            index = None
            for i, ks in enumerate(self.keys):
                if id in ks:
                    index = i
                    break
            return self.values[index]
        else:
            return None


def generate_full_sets():
    fd = FinemapData()
    for lead, details in get_all_results_iter():
        lead['id'] = ler_primary_key(lead)  # add the unique ID
        # also add the colocalising hits.
        colocs = fd.query(lead['id'])
        lead['colocalisation'] = []
        if colocs is not None:
            lead['colocalisation'] = [Colocalisation(**coloc) for coloc in colocs]
        lo = LeadEqtlResult(**lead)
        n_details = len(details['variantId'])
        detail_objs = []
        for i in range(n_details):
            ith_entry = {k: v[i] for k, v in details.items()}
            ith_entry['leadEqtlResultId'] = lead['id']
            detail_objs.append(DetailedEqtlResult(**ith_entry))
        yield lo, detail_objs

def fill_eqtl_results_database():
    for lo, detail_objs in tqdm(generate_full_sets()):
        lo.save()
        [el.save() for el in detail_objs]



def get_lead_expr_data():
    """returns a zlib compressed and pickled dict of lists describing the """
    import h5py
    import json
    import pandas as pd
    import zlib
    for cell_type in ["CD14", "CD19", "CD4", "CD8", "mac", "PLA", "CD15"]:
        fh = h5py.File("/nfs/research1/stegle/users/rkreuzhu/webapp_data/expression_data/%s_lead_exprplots.hdf5"%cell_type, "r")
        for k in fh:
            ret = {k2: fh[k][k2][:] if "S" not in fh[k][k2][:].dtype.str else fh[k][k2][:].astype(str)
                   for k2 in ['dot_labels', 'genoLabels', 'geno_col', 'pheno_col']}
            ret['dot_labels'] = pd.Series(ret['dot_labels']).str.replace("Belgium", "CEDAR").replace("Oxford", "WTCHG").replace("Inhouse", "BLUEPRINT").replace("mac", "MAC").replace("PLA", "PLT").values
            ret2 = {k2: ret[k2].tolist() for k2 in ret}
            thebytes = zlib.compress(json.dumps(ret2).encode("utf8"))
            yield k, thebytes
        fh.close()


def get_median_expr_data():
    """returns a zlib compressed and pickled dict of lists describing the """
    import h5py
    import json
    import zlib
    for cell_type in ["CD14", "CD19", "CD4", "CD8", "mac", "PLA", "CD15"]:
        fh = h5py.File("/nfs/research1/stegle/users/rkreuzhu/webapp_data/expression_data/%s_median_expr.hdf5"%cell_type, "r")
        expr_dict = {k: fh[k][:] if "S" not in fh[k][:].dtype.str else fh[k][:].astype(str)
                     for k in ["probe", "hgnc", "med_expr"]}
        fh.close()
        expr_dict = {k: expr_dict[k].tolist() for k in expr_dict}
        ct = cell_type.replace("mac", "MAC").replace("PLA", "PLT")
        thebytes = zlib.compress(json.dumps(expr_dict).encode("utf8"))
        yield {"cellType": ct, "data": thebytes}


def fill_lead_eqtl_expression():
    for k, bytes in tqdm(get_lead_expr_data()):
        k = k.replace("mac", "MAC").replace("PLA", "PLT")
        this_id = "lee" + k[3:]
        if len(LeadEqtlExpression.objects(id=this_id)) != 0:
            continue
        #try:
        LeadEqtlResult.objects(id=k)[0].modify(leadEqtlExpressionId=this_id)
        LeadEqtlExpression(**{"id": this_id, "leadEqtlResultId": k, "data": bytes}).save()
        #except IndexError as e:
        #    continue

def fill_median_expr_data():
    for kwargs in tqdm(get_median_expr_data()):
        ProbeExpressionRank(**kwargs).save()


def fill_tested_probes():
    for kwargs in tqdm(get_tested_probes()):
        kwargs['probeId'] = str(kwargs['probeId'])
        TestedProbes(**kwargs).save()

def remove_highP_assocs():
    DetailedEqtlResult.objects(pValue__gt = 0.05).delete()