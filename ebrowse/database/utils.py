__author__ = 'roman'
import os
import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm

gt_path_dataset = {}
gt_path_dataset["Belgium"] = "/nfs/research1/stegle/projects/eQTLHIC/data/ULG/Genotypes/reimpute/impute_merged_byChr.hdf5"
gt_path_dataset["Oxford"] = "/nfs/research1/stegle/projects/eQTLHIC/data/Oxford/genotypes/reimpute/impute_merged_byChr.hdf5"
gt_path_dataset["Inhouse"] = "/nfs/research1/stegle/projects/eQTLHIC/data/pla_inhouse/impute_merged_byChr.hdf5"
gt_path_dataset["Cardiogenics"] = "/nfs/research1/stegle/users/rkreuzhu/Cardiogenics_data/temp/impute_merged_byChr.hdf5"
gt_path_dataset["CardiogenicsCBR"] = "/nfs/research1/stegle/users/rkreuzhu/Cardiogenics_data/temp/impute_merged_byChr.hdf5"
gt_path_dataset["WP10"] = "/nfs/research1/stegle/projects/eQTLHIC/data/pla_inhouse/impute_merged_byChr.hdf5"


pt_path_dataset = {}
pt_path_dataset["Belgium"] = "/nfs/research1/stegle/projects/eQTLHIC/data/ULG/Preproc_RK/%s_phenotypes_merged_byChr.hdf5"
pt_path_dataset["Oxford"] = "/nfs/research1/stegle/projects/eQTLHIC/data/Oxford/raw_pheno/prepped/%s_phenotypes_merged_byChr.hdf5"
pt_path_dataset["Inhouse"] = "/nfs/research1/stegle/projects/eQTLHIC/data/pla_inhouse/%s_phenotypes_merged_byChr.hdf5"
pt_path_dataset["Cardiogenics"] = "/nfs/research1/stegle/projects/eQTLHIC/data/Cardiogenics/raw_expr/BeadStudioOutput/prepped/%s_phenotypes_merged_byChr.hdf5"
pt_path_dataset["CardiogenicsCBR"] = "/nfs/research1/stegle/projects/eQTLHIC/data/Cardiogenics/raw_expr/BeadStudioOutput/prepped_CBR/%s_phenotypes_merged_byChr.hdf5"
pt_path_dataset["WP10"] = "/nfs/research1/stegle/projects/eQTLHIC/BPWP10/chromwise/%s_phenotypes_merged_byChr.hdf5"

avail_data_wwp10 = {"Belgium":["CD14", "CD15", "CD19", "CD4", "CD8", "PLA"], "Oxford":["CD14", "CD15", "CD19", "CD14-IFN", "CD14-LPS2", "CD14-LPS24"], "Inhouse":["PLA"], "Cardiogenics":["CD14", "mac"], "WP10":["CD4", "CD14", "CD15"]}
avail_data_full = {"Belgium":["CD14", "CD15", "CD19", "CD4", "CD8", "PLA"], "Oxford":["CD14", "CD15", "CD19", "CD14-IFN", "CD14-LPS2", "CD14-LPS24"], "Inhouse":["PLA"], "Cardiogenics":["CD14", "mac"]}
avail_data = {"Belgium":["CD14", "CD15", "CD19", "CD4", "CD8", "PLA"], "Oxford":["CD14", "CD15", "CD19"], "Inhouse":["PLA"], "Cardiogenics":["CD14", "mac"]}
mega_data = {"Belgium-Oxford-Inhouse-Cardiogenics":["CD14", "CD15", "CD19", "PLA"], "Belgium":["CD4", "CD8"], "Cardiogenics":["mac"]}


all_singals = True

output_bpath = "/nfs/research1/stegle/users/rkreuzhu/pleiotropy_fwdFDR"
if all_singals:
    output_bpath += "_allsig"


def get_leads_allsig(dataset, cell_type, fdr_thresh = 0.01):
    if dataset in ["Belgium-Oxford-Inhouse-Cardiogenics", "Belgium-Oxford-Inhouse-CardiogenicsCBR"]:
        ofolder = "/nfs/research1/stegle/users/rkreuzhu/eqtl3rdgen/%s_%s_win1000000_minMAF0.01_covs_FDRfwd0.01/summary_multi_fwdsel.hdf5" % (
            dataset, cell_type)
    else:
        ofolder = "/nfs/research1/stegle/users/rkreuzhu/eqtl3rdgen/%s_%s_win1000000_covs_FDRfwd0.01/summary_multi_fwdsel.hdf5" % (
        dataset, cell_type)
    hdf = pd.HDFStore(ofolder)
    lead_eqtls = hdf["summary"]
    hdf.close()
    lead_eqtls = lead_eqtls.loc[lead_eqtls["pv_bonf_BH"].values < fdr_thresh, :]
    # clean up and remove SNPs that had not been found significant in the last round of fwd-selection
    if "fwd_sel_itr" in lead_eqtls:
        sel_lines = []
        for gene in tqdm(lead_eqtls["gene"].unique()):
            lines_here = np.where((lead_eqtls["gene"]==gene).values)[0].tolist()
            itrs_here = lead_eqtls["fwd_sel_itr"].iloc[lines_here]
            if itrs_here.max() >= itrs_here.shape[0]:
                first_missing = np.where(~np.in1d(np.arange(itrs_here.max()+1), itrs_here))[0].min()
                sel_lines.extend(np.array(lines_here)[itrs_here < first_missing].tolist())
            else:
                sel_lines.extend(lines_here)
        sel_lines = np.sort(sel_lines)
        lead_eqtls = lead_eqtls.iloc[sel_lines,:]
    #
    return lead_eqtls

def assert_valid_objs(obj):
    return [el.decode("utf8") if hasattr(el, "decode") else el for el in obj]

def get_detailed_results_entry(fh, fh_ld, probeId, fwdIter, cell_type, dataset):
    ret = {}
    fwdIter = int(fwdIter)
    if int(fwdIter) < 0 or int(fwdIter) >= fh[probeId]['betas_added'].shape[0]:
        return None
    ret['beta'] = fh[probeId]['betas_added'][:][fwdIter, :].tolist()
    ret['pValue'] = fh[probeId]['pvadded'][:][fwdIter, :].tolist()
    ret['snpId'] = fh[probeId]['rs'][:].tolist() # correct this
    ret['variantId'] = fh[probeId]['rs'][:].tolist()
    snp_sel_lead = assert_valid_objs([ret['variantId'][fh[probeId]['iadded'][:][fwdIter]]])[0]
    ret['ldR2'] = (fh_ld[probeId][snp_sel_lead].values ** 2).tolist()
    ret['chromosome'] = fh[probeId]['chrom'][:].tolist()
    ret['start'] = fh[probeId]['pos'][:].tolist()
    ret['geneSymbol'] = fh[probeId]['hgnc'][:].tolist() * len(ret['beta'])
    ret['probeId'] = [probeId] * len(ret['beta'])
    ret['fwdIter'] = [str(fwdIter)] * len(ret['beta'])
    ret['cellType'] = [cell_type] * len(ret['beta'])
    ret['dataset'] = [dataset] * len(ret['beta'])
    ret = {k: assert_valid_objs(v) for k, v in ret.items()}
    return ret


def get_all_leads():
    eqtls_df = pd.read_csv("/nfs/research1/stegle/users/rkreuzhu/webapp_data/all_leads_w_rs.txt.gz", sep="\t", compression="gzip")
    eqtls_df['chrom'] = eqtls_df['chrom'].astype(str)
    eqtls_df['gene'] = eqtls_df['gene'].astype(str)
    eqtls_df['dataset'] = "merged"
    eqtls_df = eqtls_df[["chrom", "pos", "betas_added", "pvadded", "rsids", "rs", "hgnc", "gene", "fwd_sel_itr",
                         "cell_type", "dataset"]]
    eqtls_df.columns = ["chromosome", "start", "beta", "pValue", "snpId", "variantId", "geneSymbol", "probeId",
                        "fwdIter", "cellType", "dataset"]
    ret = eqtls_df.to_dict("list")
    ret = {k: assert_valid_objs(v) for k, v in ret.items()}
    return ret


def get_all_results_iter():
    eqtls_df = pd.read_csv("/nfs/research1/stegle/users/rkreuzhu/webapp_data/all_leads_w_rs.txt.gz", sep="\t",
                           compression="gzip")
    eqtls_df['chrom'] = eqtls_df['chrom'].astype(str)
    eqtls_df['gene'] = eqtls_df['gene'].astype(str)
    eqtls_df['dataset'] = "merged"
    eqtls_df = eqtls_df[["chrom", "pos", "betas_added", "pvadded", "rsids", "rs", "hgnc", "gene", "fwd_sel_itr",
                         "cell_type", "dataset", "file"]]
    eqtls_df.columns = ["chromosome", "start", "beta", "pValue", "snpId", "variantId", "geneSymbol", "probeId",
                        "fwdIter", "cellType", "dataset","file"]

    ret = eqtls_df[[c for c in eqtls_df.columns if c != "file"]].to_dict("list")
    ret = {k: assert_valid_objs(v) for k, v in ret.items()}

    for i in range(eqtls_df.shape[0]):
        entries = eqtls_df.iloc[i, :]
        new_path = entries['file'].replace("research2", "research1")
        fh = h5py.File(new_path, "r")
        j = int(new_path.split("/")[-1].split("_")[-1].split(".")[0])
        bpath = "/".join(new_path.rstrip("/").split("/")[:-1])
        ld_file = bpath + "/lead_LDs_%d.hdf5" % j
        fh_ld = pd.HDFStore(ld_file, "r")
        details_dict = get_detailed_results_entry(fh, fh_ld, entries["probeId"], entries["fwdIter"],
                                                  entries["cellType"], entries["dataset"])
        fh.close()
        fh_ld.close()
        yield {k:v[i] for k,v in ret.items()}, details_dict



def get_all_detailed_results():
    eqtls_df = pd.read_csv("/nfs/research1/stegle/users/rkreuzhu/webapp_data/all_leads_w_rs.txt.gz", sep="\t",
                           compression="gzip")
    eqtls_df = eqtls_df.sort_values[["file"]]
    fh = None
    fh_ld = None
    fpath = ""
    for i in range(eqtls_df.shape[0]):
        entries = eqtls_df.iloc[i,:]
        new_path = entries['file']
        if new_path != fpath:
            if fh is not None:
                try:
                    fh.close()
                    fh_ld.close()
                except:
                    pass
            fh = h5py.File(new_path, "r")
            j = int(new_path.split("/")[-1].split("_")[-1].split(".")[0])
            bpath = "/".join(new_path.rstrip("/").split("/")[:-1])
            ld_file = bpath + "/lead_LDs_%d.hdf5" % j
            fh_ld = h5py.File(ld_file, "r")
            fpath = new_path
