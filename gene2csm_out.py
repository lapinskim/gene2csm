# Imports
import pandas as pd
from decimal import Decimal


def prepare_output(result, n_store=50, write=True):
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
            xlsx_fn = key + '.csm.xlsx'
            print('Writing file {}.'.format(xlsx_fn))
            final_dict[key].to_excel(xlsx_fn)
    return final_dict
