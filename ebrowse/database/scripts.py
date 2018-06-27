from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from ebrowse.database import Base, LeadEqtlResult, DetailedEqtlResult, session
from ebrowse.database.utils import get_all_leads, get_all_results_iter
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



if __name__ == "__main__":
    # Now generate all the necessary objects and commit them
    for lo, detail_objs in tqdm(generate_full_sets()):
        session.add_all([lo] + detail_objs)
        session.commit()

