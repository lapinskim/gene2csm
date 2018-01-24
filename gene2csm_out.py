# Imports
import pandas as pd
from decimal import Decimal
import logging


def prepare_output(result, n_store=50, write=True, file_prefix=None):
    '''
    Prepare the data frames, score and store the output.
    '''

    final_dict = {}
    for name, output in result:

        labels = ['coords',
                  'seq',
                  'GC',
                  'ent',
                  'dG_AA',
                  'G_A',
                  'bitscore',
                  'nident']
        df = pd.DataFrame.from_records(output, columns=labels)
        df[['chr', 'start', 'end']] = pd.DataFrame([e for e in df['coords']])
        del df['coords']

        # reset index to match the start position (shifted by multiprocessing)
        df = df.sort_values('start').reset_index(drop=True)

        # TODO: Move scoring to a proper function
        score = []
        df = df.sort_values('ent', ascending=False)
        points = {}
        point = 1
        for e in df['ent'].unique():
            points[e] = point
            point += 1
        for e in df['ent']:
            score.append(points[e])
        score
        df['score'] = score

        blast_score = []
        blast_points = {}
        df = df.sort_values('bitscore')
        num_uniq_bits = len(df['bitscore'].unique())
        point_step = Decimal(max(points.values())) / num_uniq_bits
        blast_point = point_step
        for e in df['bitscore'].unique():
            blast_points[e] = blast_point
            blast_point += point_step
        for e in df['bitscore']:
            blast_score.append(blast_points[e])
        df['score'] += blast_score

        df = df.sort_values('score')

        if n_store <= 0:
            final_dict[name] = df
        else:
            final_dict[name] = df.head(n_store)
    if write:
        for key in final_dict:
            if file_prefix:
                xlsx_fn = file_prefix + key + '.csm.xlsx'
            else:
                xlsx_fn = key + '.csm.xlsx'
            log.info('Writing file {}.', xlsx_fn)
            final_dict[key].to_excel(xlsx_fn)
    return final_dict


# set logging

# create logger with the module name
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
# create console handler for all but info levels
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
# now for info
chi = logging.StreamHandler()
chi.setLevel(logging.INFO)
# now for DEBUG
chd = logging.StreamHandler()
chd.setLevel(logging.DEBUG)
# create formatters and add them to the handlers
basic_formatter = logging.Formatter(
    fmt='%(asctime)s:%(name)s:%(levelname)s: %(message)s', style='{')
info_formatter = logging.Formatter(fmt='%(message)s', datefmt=None, style='{')
ch.setFormatter(basic_formatter)
chi.setFormatter(info_formatter)
chd.setFormatter(basic_formatter)
# add the handlers to the logger
log.addHandler(ch)
log.addHandler(chi)
log.addHandler(chd)
