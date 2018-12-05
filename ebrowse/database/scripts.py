from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from ebrowse.database import Base, LeadEqtlResult, DetailedEqtlResult, ColocalisationResult, use_session, session_factory, scoped_session
from ebrowse.database.utils import get_all_leads, get_all_results_iter, get_colocalised_hits
from tqdm import tqdm

def get_all_lead_objects():
    leads = get_all_leads()
    n_entries = len(leads.values().__next__())
    lead_objs = []
    for i in range(n_entries):
        ith_entry = {k:v[i] for k,v in leads.items()}
        lead_objs.append(LeadEqtlResult(**ith_entry))
    return lead_objs


def generate_full_sets():
    for lead, details in get_all_results_iter():
        lo = LeadEqtlResult(**lead)
        n_details = len(details['variantId'])
        detail_objs = []
        for i in range(n_details):
            ith_entry = {k: v[i] for k, v in details.items()}
            ith_entry['lead_result'] = lo
            detail_objs.append(DetailedEqtlResult(**ith_entry))
        yield lo, detail_objs




Session = scoped_session(session_factory)
session = Session()

def relabel_cell_types():
    relabel = [("PLA", "PLT"), ("mac", "MAC")]
    with use_session() as session:
        for lable_tuple in relabel:
            query = session.query(LeadEqtlResult).filter(LeadEqtlResult.cellType == lable_tuple[0])
            for lead_obj in query:
                lead_obj.cellType = lable_tuple[1]
            session.commit()
            ctr = 0
            while True:
                print("running iter: %d, %s"%(ctr, str(lable_tuple)))
                ctr += 1
                r_objs = []
                objs = []
                query = session.query(DetailedEqtlResult).filter(DetailedEqtlResult.cellType == lable_tuple[0])
                stop_it = False
                try:
                    objs = query[:1000]
                except Exception as e:
                    try:
                        objs = query[:query.count()]
                        stop_it = True
                    except Exception as e:
                        print(e)
                        stop_it = True
                for detailed_obj in objs:
                    detailed_obj.cellType = lable_tuple[1]
                    r_objs.append(detailed_obj)
                if len(objs) == 0:
                    break
                session.commit()
                for el in r_objs:
                    del el
                if stop_it:
                    break
            session.commit()


def import_finemap_data():
    fmps = get_colocalised_hits()
    col_objs = []
    for entry in fmps:
        for i in range(len(entry[list(entry.keys())[0]])):
            ith_entry = {k: v[i] for k, v in entry.items()}
            col_objs.append(ColocalisationResult(**ith_entry))
    return col_objs

# for updating the rs-ids:
def get_unique_variant_ids():
    import gzip
    with use_session() as session:
        with gzip.open("unique_var_ids.txt.gz", "w") as ofh:
            for id in session.query(DetailedEqtlResult.variantId).distinct():
                ofh.write(id + "\n")


if __name__ == "__main__":
    # Now generate all the necessary objects and commit them
    for lo, detail_objs in tqdm(generate_full_sets()):
        session.add_all([lo] + detail_objs)
        session.commit()

