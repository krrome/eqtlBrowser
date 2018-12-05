from flask import Flask, request
from flask import jsonify, send_from_directory
from flask import abort
from ebrowse.eqtl import LeadEqtls, SubEqtls
from ebrowse import PATHS
import re
from collections import OrderedDict
from ebrowse.shipping.lead_table import column_order as lead_table_column_order
from ebrowse.shipping.lead_table import column_labels as lead_table_column_labels
from ebrowse.shipping.lead_table import relabel_and_order as lead_relabel_and_order
from ebrowse.mongo_db import get_leads, get_detailed, get_lead_by_region, get_lead_expression, get_tested_probes, get_expression_hist
from ebrowse.datatable_handling import lead_table_formatting, lead_bootstrap_formatting
from flask import render_template
from ebrowse.epi_annotation import get_celltypes, get_tracks

app = Flask(__name__)

def get_return_dict():
    return {"format": "0.1"}

global_address = "http://127.0.0.1:%d"%PATHS['flask_port']


def split_dict_by(in_dict, split_keys = ["dataset", "cellType"]):
    ret = {}
    num_entries = max([len(in_dict[k]) if in_dict[k] is not None else 1 for k in in_dict])
    for row in range(num_entries):
        ret_subobj = ret
        for split_key in split_keys:
            key_value = in_dict[split_key][row]
            if key_value not in ret_subobj:
                ret_subobj[key_value] = {}
            ret_subobj = ret_subobj[key_value]
        for k in in_dict:
            if in_dict[k] is None:
                ret_subobj[k] = None
            else:
                if k not in ret_subobj:
                    ret_subobj[k] = []
                ret_subobj[k].append(in_dict[k][row])
    return ret



@app.route('/eqtl/<dataset_id>')
@app.route('/eqtl/<dataset_id>/<high_detail_varid>')
@app.route('/eqtl/<dataset_id>/<chromosome>/<start>/<end>')
@app.route('/eqtl/<dataset_id>/<chromosome>/<start>/<end>/<high_detail_varid>')
def get_lead_eqtl_old(dataset_id, high_detail_varid=None, chromosome=None, start=None, end=None):
    """
    :param dataset_id: dataset identifier. e.g.: "CD14_WTCHG"
    :param high_detail_varid: "_".join([probe_id, variant_id, fwd_iter])
    :param chromosome: string, no "chr"
    :param start: integer
    :param end: integer
    :return: 
    """
    if (chromosome is not None) and (start is not None) and (end is not None):
        try:
            start = int(start)
            end = int(end)
        except Exception:
            return jsonify(None)
    ret = get_lead_by_region(chromosome, start, end)

    #if ret['ldR2'] is not None and isinstance(ret['ldR2'], list):
    #    if all([el is None for el in ret['ldR2']]):
    #        ret['ldR2'] = None
    dict_out = get_return_dict()
    dict_out['lead_eqtl'] = {dataset_id: ret}
    if high_detail_varid is not None:
        probe_id, variant_id, fwd_iter  = high_detail_varid.split("_")
        if dataset_id not in sub_eqtls:
            sub_eqtls[dataset_id] = SubEqtls(dataset_id)
        detail_ret = sub_eqtls[dataset_id].get(variant_id, probe_id, fwd_iter)
        if detail_ret is not None:
            if detail_ret['ldR2'] is not None and isinstance(detail_ret['ldR2'], list):
                if all([el is None for el in detail_ret['ldR2']]):
                    detail_ret['ldR2'] = None
            dict_out['detailed_eqtl'] = {dataset_id: {variant_id: {probe_id: {fwd_iter: detail_ret}}}}
    return jsonify(dict_out)

@app.route('/eqtl')
def get_lead_eqtl():
    """
    :return: 
    """
    query = {}
    query['cell_type'] = request.args.get('cellType', default="")
    query['dataset'] = request.args.get('dataset', default="")
    query['chromosome'] = request.args.get('chromosome', default="")
    query['start'] = request.args.get('start', default=None, type=int)
    query['end'] = request.args.get('end', default=None, type=int)
    query = {k:v if v != "" else None for k,v in query.items() }

    ret = get_lead_by_region(**query)
    dict_out = get_return_dict()
    ret = split_dict_by(ret)
    dict_out['lead_eqtl'] = ret
    #dict_out['lead_eqtl'] = {"bla":ret}
    return jsonify(dict_out)


# the all_eqtl has to be integrated into the lead_eqtl query as a single probe_variant identifier.
# in igv in the pop-up it is necessary to add a "link" that will set a flag to request the probe_variant and reload the track.
# the detailed_eqtl eqtls are then parsed individually and coloured if possible.

@app.route('/all_eqtl/<dataset_id>/<variantId>/<probeId>/<fwdIter>')
def get_sub_eqtl(dataset_id, variantId, probeId, fwdIter):
    #if dataset_id not in sub_eqtls:
    #    sub_eqtls[dataset_id] = SubEqtls(dataset_id)
    #ret = sub_eqtls[dataset_id].get(variantId, probeId, fwdIter)
    ret = get_detailed(cell_type="CD14", variant_id=variantId, probe_id=probeId, fwd_iter=fwdIter)
    if ret is None:
        return jsonify(None)
    if isinstance(ret['ldR2'], list):
        if all([el is None for el in ret['ldR2']]):
            ret['ldR2'] = None
    dict_out = get_return_dict()
    dict_out['detailed_eqtl'] = {dataset_id:{variantId:{probeId:{fwdIter:ret}}}}
    return jsonify(dict_out)

@app.route('/all_eqtl')
def get_sub_eqtl_get():
    query = {}
    query['leadEqtlResultId'] = request.args.get('id', default="")
    query['cell_type'] = request.args.get('cellType', default="")
    query['dataset'] = request.args.get('dataset', default="")
    query['variant_id'] = request.args.get('variantId', default="")
    query['probe_id'] = request.args.get('probeId', default="")
    query['fwd_iter'] = request.args.get('fwdIter', default="")
    query_args = {k: v if v != "" else None for k, v in query.items()}
    #'dataset_id', 'variantId', 'probeId', and 'fwdIter'
    ret, lead_ret = get_detailed(**query_args)
    if ret is None:
        return jsonify(None)
    if isinstance(ret['ldR2'], list):
        if all([el is None for el in ret['ldR2']]):
            ret['ldR2'] = None
    dict_out = get_return_dict()
    if lead_ret is not None:
        dict_out['detailed_eqtl'] = {
            lead_ret['dataset']: {lead_ret['variantId']: {lead_ret['probeId']: {lead_ret['fwdIter']: ret}}}}
    else:
        dict_out['detailed_eqtl'] = {query['dataset']:{query['variant_id']:{query['probe_id']:{query['fwd_iter']:ret}}}}
    return jsonify(dict_out)


@app.route('/tested_probes')
def get_tested_probes_table():
    import json
    query = {}
    query['search'] = request.args.get('search', default="")
    query['sort'] = request.args.get('sort', default="")
    query['order'] = request.args.get('order', default="")
    query['offset'] = request.args.get('offset', default="")
    query['limit'] = request.args.get('limit', default="")
    query['filterDict'] = json.loads(request.args.get('filter', default="null"))
    results, all_entries, num_filtered = get_tested_probes(exact_match = False, **query)
    ret = lead_bootstrap_formatting(results, all_entries, num_filtered, reply_to_search = True)
    return jsonify(ret)


def parse_datatable_args(in_dict):
    ret = {}
    for k, v in in_dict.items():
        key_elms = [tk for tk in re.split(r'\[|\]', k) if tk != ""]
        rec_el = ret
        for sub_k in key_elms[:-1]:
            if sub_k not in rec_el:
                rec_el[sub_k] = {}
            rec_el = rec_el[sub_k]
        rec_el[key_elms[-1]] = v
    return ret

def parse_datatable_order(order_dict, column_labels):
    ret_order_dict = OrderedDict()
    for k in range(len(order_dict.keys())):
        v = order_dict[str(k)]
        ret_order_dict[column_labels[int(v['column'])]] = v['dir'].upper()
    return ret_order_dict


@app.route('/lead_table')
def get_lead_table():
    args_dict = parse_datatable_args(request.args)
    # most of the args_dict can be ignored: only relevant is
    entries_per_page = int(args_dict['length'])
    draw_key = int(args_dict['draw'])
    order = parse_datatable_order(args_dict['order'], lead_table_column_order)
    search_str = args_dict['search']['value']
    search_is_regex = args_dict['search']['regex'] == 'true'
    start = int(args_dict['start'])
    args_to_search = {"order_by": order, 'offset': start, 'per_page': entries_per_page}

    if search_str != "":
        args_to_search["filter_equals_any"] = search_str

    lead_results, all_entries, num_filtered = get_leads(**args_to_search)
    lead_results_ordered = lead_relabel_and_order(lead_results)
    ret = lead_table_formatting(lead_results_ordered, draw_key, all_entries, num_filtered)
    return jsonify(ret)


@app.route('/new_lead_table')
def get_bootstrap_lead_table():
    import json
    query = {}
    query['search'] = request.args.get('search', default="")
    query['sort'] = request.args.get('sort', default="")
    query['order'] = request.args.get('order', default="")
    query['offset'] = request.args.get('offset', default="")
    query['limit'] = request.args.get('limit', default="")
    query['filterDict'] = json.loads(request.args.get('filter', default="null"))
    lead_results, all_entries, num_filtered = get_leads(exact_match = False, **query)
    ret = lead_bootstrap_formatting(lead_results, all_entries, num_filtered, reply_to_search = True)
    return jsonify(ret)

@app.route('/expr_boxplot_data')
def get_expr_boxplot_data():
    query = {}
    id = request.args.get('id', default=None)
    if id is not None and id !="":
        query["id"] = id
    else:
        for field_name in ['cellType', 'dataset', 'variantId', 'probeId']:
            query[field_name] = request.args.get(field_name, default="")
        query = {k: v if v != "" else None for k, v in query.items()}
        if any([query[k] is None for k in query]):
            abort(404)
    ret = get_lead_expression(**query)
    return jsonify(ret)

@app.route('/expr_boxplot')
def get_expr_boxplot():
    query = {}
    id = request.args.get('id', default=None)
    if id is not None and id != "":
        query["id"] = id
    else:
        abort(404)
    return render_template('scatterplot.html', expr_plot_link = "expr_boxplot_data", **query)


@app.route('/expr_hist_data')
def get_expr_hist_data():
    query = {}
    query["cellType"] = request.args.get('cellType', default=None)
    ret = None
    if query["cellType"] is not None and query["cellType"] !="":
        ret = get_expression_hist(**query)
    return jsonify(ret)


@app.route('/')
def index():
    var_id_col = lead_table_column_order.index("variantId")
    probe_id_col = lead_table_column_order.index("probeId")
    fwd_iter_col = lead_table_column_order.index("fwdIter")
    cell_type_col = lead_table_column_order.index("cellType")
    dataset_col = lead_table_column_order.index("dataset")
    hgnc_col = lead_table_column_order.index("geneSymbol")
    return render_template('table_igv_iframe.html', port=PATHS['flask_port'], colnames=lead_table_column_labels,
                           var_id_col = var_id_col, probe_id_col = probe_id_col, fwd_iter_col=fwd_iter_col,
                           cell_type_col = cell_type_col, dataset_col = dataset_col, hgnc_col= hgnc_col,
                           cell_types = get_celltypes())

@app.route('/igv')
def render_igv():
    return render_template('igv_view.html', port=PATHS['flask_port'])


@app.route('/igvdir/<path:path>')
def send_js(path):
    return send_from_directory(PATHS['igv_repo_dir'], path)

@app.route('/bootstrapTableDir/<path:path>')
def send_bootstrap_table_dir(path):
    return send_from_directory(PATHS['bootstrap_tables_dir'], path)

@app.route('/prototable/<path:path>')
def send_prototable(path):
    return send_from_directory("/nfs/research1/stegle/users/rkreuzhu/webapp_data/ebrowser/table_prototyping", path)

@app.route('/get_epi_tracks')
def get_epi_tracks():
    cell_type = request.args.get('cellType', default="")
    return jsonify(get_tracks(cell_type, global_address + "/epi_tracks"))

@app.route('/epi_tracks/<path:path>')
def epi_tracks(path):
    return send_from_directory(PATHS['epi_dir'], path)



lead_eqtls = LeadEqtls(PATHS['lead_eqtl_table'])
sub_eqtls = {}

if __name__ == '__main__':
    app.run(port=PATHS['flask_port'])