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

colocalisation_col_renaming = {"chi2_eqtl":"chi2Eqtl","chi2_gwas":"chi2Gwas","compared_in_finemap":"comparedInFinemap","index_eqtl":"indexGroupEqtl",
 "index_gwas":"indexGroupGwas","pvalues_eqtl":"pvaluesEqtl","pvalues_gwas":"pvaluesGwas","snp":"variantId","snp_prob_eqtl":"snpProbEqtl","snp_prob_gwas":"snpProbGwas"}



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

def get_detailed_results_entry(fh, fh_ld, reestim_fh, probeId, fwdIter, cell_type, dataset):
    ret = {}
    fwdIter = int(fwdIter)
    if int(fwdIter) < 0 or int(fwdIter) >= fh[probeId]['betas_added'].shape[0]:
        return None
    if reestim_fh is not None:
        ret['beta'] = reestim_fh[probeId]['betas_added_re'][:][fwdIter, :].tolist()
        ret['pValue'] = reestim_fh[probeId]['pvadded_re'][:][fwdIter, :].tolist()
    else:
        ret['beta'] = fh[probeId]['betas_added'][:][fwdIter, :].tolist()
        ret['pValue'] = fh[probeId]['pvadded'][:][fwdIter, :].tolist()
    ret['pValue'] = np.clip(fh[probeId]['pvadded'][:][fwdIter, :].tolist(), 1e-310, 1).tolist()
    ret['snpId'] = [None]*len(ret['pValue'])
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
    ret = pd.DataFrame(ret).query('pValue < 0.05').to_dict("list")
    #ret['snpId'] = get_rs_ids_batched(ret['variantId'])  # correct this
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
    eqtls_df = eqtls_df[["chrom", "pos", "betas_added", "betas_added_re", "pvadded", "pvadded_re", "rsids", "rs", "hgnc", "gene", "fwd_sel_itr",
                         "cell_type", "dataset", "file", "ensembl_gid", "maf", "pv_bonf_BH", "pv_bonf_BH_re"]]
    eqtls_df.columns = ["chromosome", "start", "beta", "betaFM", "pValue", "pValueFM", "snpId", "variantId", "geneSymbol", "probeId",
                        "fwdIter", "cellType", "dataset","file", "ensembleGId", "AAF", "fdr", "fdrFM"]
    maf= eqtls_df["AAF"].values
    maf[maf > 0.5] = 1-maf[maf > 0.5]
    eqtls_df['maf'] = maf
    eqtls_df = eqtls_df[[c for c in eqtls_df.columns if c != "AAF"]]
    for c in ['pValue', 'pValueFM', 'fdr', 'fdrFM']:
        eqtls_df[c] = np.clip(eqtls_df[c], 1e-310, 1)
    eqtls_df['cellType'] = eqtls_df['cellType'].replace("PLA", "PLT").replace("mac", "MAC")
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
        reestim_file = new_path[:-len(".hdf5")] + "_BReestim.hdf5"
        reestim_fh = h5py.File(reestim_file, "r")
        details_dict = get_detailed_results_entry(fh, fh_ld, reestim_fh, entries["probeId"], entries["fwdIter"],
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


def get_colocalised_hits():
    from glob import glob
    finemaps ={}
    colocs_ffmt = "/nfs/research1/stegle/users/rkreuzhu/colocalisation_pw/finemap_data/finemap_pp_*.txt"
    for f in glob(colocs_ffmt):
        tks = f.split("/")[-1].split(".")[0].split("_")[2:]
        cell_type = tks[0]
        trait = "_".join(tks[1:-1])
        gene = tks[-1]
        if cell_type not in finemaps:
            finemaps[cell_type] = {}
        if gene not in finemaps[cell_type]:
            finemaps[cell_type][gene] = {}
        df = pd.read_csv(f, sep="\t")
        df['cellType'] = cell_type
        df['probeId'] = gene
        df['trait'] = trait
        df.columns = [colocalisation_col_renaming[col] if col in colocalisation_col_renaming else col for col in
                      df.columns]
        df = df[[col for col in df.columns if col not in ["snp_log10bf_gwas", "snp_log10bf_eqtl", "pos"]]]
        coloc_sel = (df["snpProbEqtl"] * df["snpProbGwas"] > 0.5)
        if coloc_sel.sum() >0:
            ret = df.loc[coloc_sel]
            ret = ret.to_dict("list")
            ret = {k: assert_valid_objs(v) for k, v in ret.items()}
            yield ret


def get_rs_ids(var_ids):
    import cyvcf2
    from tqdm import tqdm
    dbsnp = cyvcf2.VCF("/nfs/research1/stegle/users/rkreuzhu/reference/dbsnp_hg19/All_20180423.vcf.gz")
    rss = []
    for var_id in tqdm(var_ids):
        # var_id = eqtls_df["rs"].iloc[0]
        if isinstance(var_id, bytes):
            var_id = var_id.decode()
        chrom, pos, ref, alt = var_id.split(":")
        region_str = "{0}:{1}-{2}".format(chrom, str(int(pos) - 1), pos)
        variants = dbsnp(region_str)
        rsids = []
        sel_vars = []
        for rec in variants:
            if (str(rec.POS) == pos) and (rec.REF == ref) and (alt in rec.ALT):
                rsids.append("rs" + str(dict(rec.INFO)['RS']))
                sel_vars.append(rec)
        rsids = np.unique(rsids).tolist()
        if len(rsids) == 0:
            rsids.append(".")
        rss.append(",".join(rsids))
    dbsnp.close()
    return rss


def get_rs_ids_batched(var_ids):
    import cyvcf2
    from tqdm import tqdm
    dbsnp = cyvcf2.VCF("/nfs/research1/stegle/users/rkreuzhu/reference/dbsnp_hg19/All_20180423.vcf.gz")
    if isinstance(var_ids[0], bytes):
        var_ids = [var_id.decode() for var_id in var_ids]
    chroms = [id.split(":")[0] for id in var_ids]
    pos = [int(id.split(":")[1]) for id in var_ids]
    import pdb
    pdb.set_trace()
    vars = pd.DataFrame({"chrom": chroms, "pos":pos, "ids":var_ids , "rs":[[]]*len(var_ids)})
    unq_chroms = vars["chrom"].unique()
    for chrom in tqdm(unq_chroms):
        df_sel = vars["chrom"] == chrom
        start, end = vars["pos"].loc[df_sel].min()-1, vars["pos"].loc[df_sel].max()
        region_str = "{0}:{1}-{2}".format(chrom, str(start), str(end))
        variants = dbsnp(region_str)
        for rec in tqdm(variants):
            for alt in rec.ALT:
                id_str = "%s:%d:%s:%s"%(rec.CHROM, rec.POS, rec.REF, alt)
                id_sel = vars["ids"] == id_str
                if id_sel.any():
                    vars.loc[id_sel, "rs"].iloc[0].append("rs" + str(dict(rec.INFO)['RS']))
    import pdb
    pdb.set_trace()
    vars["rs_str"] = vars["rs"].apply(lambda x: ",".join(x) if len(x) != 0 else ".")
    dbsnp.close()
    return vars["rs_str"].tolist()



def get_tested_probes():
    all_tested = pd.read_csv("/nfs/research1/stegle/users/rkreuzhu/webapp_data/all_tested_probes.txt.gz")
    for i, row in all_tested.iterrows():
        yield row.to_dict()

def save_all_tested_probes():
    dfs = []
    for ds, cell_types in mega_data.items():
        for cell_type in cell_types:
            all_tested = get_leads_allsig(ds, cell_type, 2.0)[["gene", "hgnc"]]
            all_tested.columns = ["probeId", "geneSymbol"]
            all_tested["cellType"] = cell_type.replace("PLA", "PLT").replace("mac", "MAC")
            all_tested = all_tested.loc[~all_tested["probeId"].duplicated(),:]
            dfs.append(all_tested)
    all_tested = pd.concat(dfs, axis=0)
    all_tested.to_csv("/nfs/research1/stegle/users/rkreuzhu/webapp_data/all_tested_probes.txt.gz", index=None, compression="gzip")


def get_expression_boxplots():
    # load the mana / ana object and then iterate over the lead results, generate the expression boxplots / export
    # their data
    pass
