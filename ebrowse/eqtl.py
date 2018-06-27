import pandas as pd
import h5py
from . import PATHS


class LeadEqtls(object):
    def __init__(self, source_fn):
        self.dataset = pd.read_csv(source_fn, sep="\t")
        self.dataset['chrom'] = self.dataset['chrom'].astype(str)
        self.dataset['gene'] = self.dataset['gene'].astype(str)
    #
    def get_by_region(self, dataset_id, chromosome, start, end):
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
        :param chromosome: 
        :param start: 
        :param end: 
        :return: 
        """
        subset = self.dataset
        if (chromosome is not None) and (start is not None) and (end is not None):
            subset_sel = self.dataset['chrom'] == str(chromosome).lstrip("chr")
            subset_sel &= self.dataset['pos'] >= start
            subset_sel &= self.dataset['pos'] <= end
            subset = subset.loc[subset_sel,:]
        subset["ldR2"] = None
        subset = subset[["chrom", "pos", "betas_added", "pvadded", "rsids", "rs", "hgnc", "gene", "fwd_sel_itr", "ldR2"]]
        subset.columns = ["chromosome", "start", "beta", "pValue", "snpId", "variantId", "geneSymbol", "probeId", "fwdIter", "ldR2"]
        return subset.to_dict("list")

def assert_valid_objs(obj):
    return [el.decode("utf8") if hasattr(el, "decode") else el for el in obj]


class SubEqtls(object):
    def __init__(self, dataset_id):
        self.dataset_id = dataset_id
        self.fh = h5py.File(PATHS['sample_subeqtl_data'], "r")
        self.fh_ld = pd.HDFStore(PATHS['sample_subeqtl_ld_data'], "r")
    #
    def get(self, variantId, probeId, fwdIter):
        ret = None
        if probeId in self.fh:
            #"chromosome", "start", "beta", "pValue", "snpId", "variantId", "geneSymbol", "probeId", "fwdIter", "ldR2"
            ret = {}
            fwdIter = int(fwdIter)
            if int(fwdIter) < 0 or int(fwdIter) >= self.fh[probeId]['betas_added'].shape[0]:
                return None
            ret['beta'] = self.fh[probeId]['betas_added'][:][fwdIter, :].tolist()
            ret['pValue'] = self.fh[probeId]['pvadded'][:][fwdIter, :].tolist()
            ret['snpId'] = self.fh[probeId]['rs'][:].tolist()  # this should be the actual rsXXXX id
            ret['variantId'] = self.fh[probeId]['rs'][:].tolist()
            snp_sel_lead = assert_valid_objs([ret['variantId'][self.fh[probeId]['iadded'][:][fwdIter]]])[0]
            ret['ldR2'] = (self.fh_ld[probeId][snp_sel_lead].values**2).tolist()
            ret['chromosome'] = self.fh[probeId]['chrom'][:].tolist()
            ret['start'] = self.fh[probeId]['pos'][:].tolist()
            ret['geneSymbol'] = self.fh[probeId]['hgnc'][:].tolist()*len(ret['beta'])
            ret['probeId'] = [probeId]*len(ret['beta'])
            ret['fwdIter'] = [str(fwdIter)]*len(ret['beta'])
            ret = {k:assert_valid_objs(v) for k,v in ret.items()}
        """
        chrom = variantId.split(":")[0]
        centr_bp = int(variantId.split(":")[1])
        subset = {"chromosome":[chrom]*2, "start":[centr_bp-10, centr_bp+10], "beta":[0.3, 0.4], "pValue":[1e-7, 1e-8],
                  "snpId":["rsXYZ"]*2, "variantId":["asdhb", "adadsd"], "geneSymbol":["BLABLA"]*2,"probeId":[probeId]*2,
                  "fwdIter":[fwdIter]*2, "ldR2":[0.3, 0.9]}
        """
        return ret

#se = SubEqtls("PLA")
#se.get("asda","1010035", "0")

