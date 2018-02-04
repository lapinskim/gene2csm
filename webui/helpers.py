import re
from config import basedir
import os


def split_ids(id_string):
    r = re.compile(r'\W+')
    out = [e for e in r.split(id_string) if e != '']
    return out


def get_databasefn():
    db_file = os.path.join(basedir, 'danRer.GRCz10.90.db')
    return db_file


def get_fastaidx():
    fa_fn = os.path.join(basedir, 'Danio_rerio.GRCz10.dna_sm.toplevel.fa')
    idx_fn = os.path.join(basedir, 'danRer.GRCz10.dna_sm.toplevel.idx')
    return fa_fn, idx_fn


def get_varfn():
    return os.path.join(basedir, 'danio_rerio.gvf.gz')
