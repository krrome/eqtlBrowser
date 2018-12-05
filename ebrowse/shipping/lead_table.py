from collections import OrderedDict



column_order = ["cellType", "chromosome", "geneSymbol", "variantId", "snpId", "pValue", "beta", "fwdIter", "start", "probeId", "dataset"]
column_labels = ["Cell type", "Chr", "HGNC", "Variant", "rs-id", "p-Value", "Effect size", "Indep. sig.", "Position", "Probe ID", "Dataset"]


def relabel_and_order(lead_dict, relabel = True):
    ret = OrderedDict()
    if relabel:
        for k, new_k in zip(column_order, column_labels):
            ret[new_k] = lead_dict[k]
    else:
        for k, new_k in zip(column_order, column_order):
            ret[new_k] = lead_dict[k]

    return ret