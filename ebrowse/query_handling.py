from flask import Flask, request
from flask import jsonify, send_from_directory
from ebrowse.eqtl import LeadEqtls, SubEqtls
from ebrowse import PATHS
import re
from collections import OrderedDict
from ebrowse.shipping.lead_table import column_order as lead_table_column_order

app = Flask(__name__)

def get_return_dict():
    return {"format": "0.1"}


#### TODO: ADD THE DATASET_ID TO THE RETUNRED DICT!

@app.route('/eqtl/<dataset_id>')
@app.route('/eqtl/<dataset_id>/<high_detail_varid>')
@app.route('/eqtl/<dataset_id>/<chromosome>/<start>/<end>')
@app.route('/eqtl/<dataset_id>/<chromosome>/<start>/<end>/<high_detail_varid>')
def get_lead_eqtl(dataset_id, high_detail_varid=None, chromosome=None, start=None, end=None):
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
    ret = lead_eqtls.get_by_region(dataset_id, chromosome, start, end)
    if ret['ldR2'] is not None and isinstance(ret['ldR2'], list):
        if all([el is None for el in ret['ldR2']]):
            ret['ldR2'] = None
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


# the all_eqtl has to be integrated into the lead_eqtl query as a single probe_variant identifier.
# in igv in the pop-up it is necessary to add a "link" that will set a flag to request the probe_variant and reload the track.
# the detailed_eqtl eqtls are then parsed individually and coloured if possible.

@app.route('/all_eqtl/<dataset_id>/<variantId>/<probeId>/<fwdIter>')
def get_sub_eqtl(dataset_id, variantId, probeId, fwdIter):
    if dataset_id not in sub_eqtls:
        sub_eqtls[dataset_id] = SubEqtls(dataset_id)
    ret = sub_eqtls[dataset_id].get(variantId, probeId, fwdIter)
    if ret is None:
        return jsonify(None)
    if isinstance(ret['ldR2'], list):
        if all([el is None for el in ret['ldR2']]):
            ret['ldR2'] = None
    dict_out = get_return_dict()
    dict_out['detailed_eqtl'] = {dataset_id:{variantId:{probeId:{fwdIter:ret}}}}
    return jsonify(dict_out)


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

    return jsonify(ret)

"""
#https://datatables.net/reference/option/columns.render
for the hrefs use:
$('#example').dataTable( {
  "columnDefs": [ {
    "targets": 4,
    "data": "description",
    "render": function ( data, type, row, meta ) {
      return type === 'display' && data.length > 40 ?
        '<span title="'+data+'">'+data.substr( 0, 38 )+'...</span>' :
        data;
    }
  } ]
} );
"""



@app.route('/igvdir/<path:path>')
def send_js(path):
    return send_from_directory('/Users/roman/igv_development/igv.js/', path)



lead_eqtls = LeadEqtls(PATHS['lead_eqtl_table'])
sub_eqtls = {}

if __name__ == '__main__':
    app.run(port='5002')