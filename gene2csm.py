# Imports
import gffutils
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import re
import logging
from multiprocessing import Pool
from collections import Counter
from decimal import Decimal
from pybedtools import BedTool

# Define functions


def revcomp(seq):
    '''
    Reverse complement sequence
    '''

    # check if sequence is RNA
    if seq.upper().find('U') == -1:
        # DNA
        comp_b = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'n': 'n'}
    else:
        # RNA
        comp_b = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
                  'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g', 'n': 'n'}
    rc_seq = ''.join([comp_b[b] for b in seq][::-1])
    return rc_seq


def gc(seq):
    '''
    Return GC content percentage
    '''

    seq = seq.upper()
    gc_content = 100 * Decimal(seq.count("G") + seq.count("C")) / len(seq)
    return gc_content


def transcribe(seq):
    '''
    Transcribe sequence- change T to U
    '''

    trans_seq = seq.replace('T', 'U').replace('t', 'u')
    return trans_seq


def revtranscribe(seq):
    '''
    Reverse transcribe sequence- change U to T
    '''

    trans_seq = seq.replace('U', 'T').replace('u', 't')
    return trans_seq


def get_cov(database, gene):
    '''
    Create genomic sequence coverage by CDS (without STOP codon).
    '''

    # get the length of the gene (1 based coords)
    gene_length = gene.end - gene.start + 1
    log.info('{} loci length: {}.'.format(gene.id, gene_length))
    # create a numpy array of zeros with the length of the gene
    gene_cov = np.zeros(gene_length, dtype=int)
    for cds in database.children(gene, featuretype='CDS', order_by='start'):
            # get relative exon coords in 0-based space
            rel_start = cds.start - gene.start
            rel_end = cds.end - gene.start
            gene_cov[rel_start:rel_end + 1] += 1
    log.info('Maximal genomic sequence coverage by CDS:\
{}.'.format(gene_cov.max()))
    return gene_cov


def get_int(gene, genomic_coverage):
    '''
    Get the genomic intervals (for sequence retrieval)
    '''

    gen_int = []
    cov = 0
    start = 0
    for i, c in enumerate(genomic_coverage):
        # adjust the genomic coords so they would be 0 based,
        # non-inclusive (*similar to BED*)
        if cov == 0 and c != 0:
            # start cov equals zero and now it changed-
            # store start position and coverage
            start = i
            cov = c
        elif cov != 0 and c != cov:
            # coverage was not zero but it changed-
            # store the previous interval , new start position and new coverage
            gen_int += [[gene.chrom,
                         start + gene.start - 1,
                         i + gene.start - 1,
                         cov]]
            start = i
            cov = c
        # else if the coverage is the same- do nothing
    # in the end store also the coverage of the last segment,
    # if the coverage was not equal to 0
    if cov != 0:
        gen_int += [[gene.chrom, start + gene.start - 1, i + gene.start, cov]]
    i_cnt = Counter()
    for e in gen_int:
        i_cnt[e[3]] += 1
    log.info('{} CDS coverage intervals (coverage, count): {}.'.format(
        len(gen_int),
        i_cnt.most_common()))
    return gen_int


def sub_usr(intervals, strand, user_list=None):
    '''
    Subtract user submitted *CDS containing exons* from processing.
    '''

    i_list = intervals
    if user_list is not None:
        # check if exon numbers are within the exon range of the target
        for e in user_list:
            try:
                if not 1 <= e <= len(intervals):
                    raise Exception('Wrong exon number')
            except Exception as e:
                log.exception('Number {} not within the target\'s \
CDS containing exon range: 1 - {}.'.format(e, len(intervals)))
                raise
        # convert to set, to delete, numbered exon only (unique)
        # and reverse sort to delete from the end of the list
        sorted_list = sorted(set(user_list), reverse=True)
        if strand == "-":
            for e in sorted_list:
                del i_list[-e]
        elif strand == "+":
            for e in sorted_list:
                del i_list[e - 1]
                log.info('Omitting user submitted \
*CDS containing exons*: {}.'.format(', '.join(str(e) + '.' for e in
                                              sorted_list[::-1])))
    return i_list


def sub_var(intervals, variation_fn):
    '''
    Subtract variable sequences from the genomic intervals

    For future enhancement consider soft masking the whole genome
    beforehand or splitting the GVF file on chromosomes for speed improvements.
    '''

    # make some basic checks
    try:
        if intervals == []:
            raise Exception('Empty list')
    except Exception as e:
        log.exception('No intervals to process.')
        raise
    try:
        if not os.path.exists(variation_fn):
            raise Exception('File does not exist.')
    except Exception as e:
        log.exception('GVF file does not exist.')
        raise
    try:
        if not variation_fn.endswith(('gvf', 'gvf.gz')):
            raise Exception('Bad extension')
    except Exception as e:
        log.exception('Ensure the variation \
containing file is in a GVF format. \'.gvf\' extension missing.')
        raise

    gvf_sorted_fn = 'sorted.gvf'.join(variation_fn.split('gvf'))
    # sort the variation database for faster processing and store
    # for further use
    if not os.path.exists(gvf_sorted_fn):
        log.info('Sorting and saving the new GVF file as {}.'.format(
            gvf_sorted_fn))
        gvf_bt = BedTool(variation_fn).sort().saveas(gvf_sorted_fn)
    else:
        gvf_bt = BedTool(gvf_sorted_fn)

    log.info('Subtracting variation from {}.'.format(
             gvf_bt[1]['Dbxref'].split(':')[0]))

    # make sure the format is right
    for e in intervals:
        try:
            if not len(e) == 4:
                raise Exception('Wrong format')
        except Exception as e:
            log.exception('Wrong input format')
            raise
    # perform genome arithmetic's
    gi_gv = BedTool(('\t'
                     .join([c, str(s), str(e), '.', str(cov)]) for c, s, e, cov
                     in intervals)).saveas().subtract(gvf_bt,
                                                      sorted=True,
                                                      stream=True)

    gi_novar = [[c, int(s), int(e), int(cov)] for c, s, e, n, cov in gi_gv]
    return gi_novar


def get_seq(fasta_index, intervals, coverage, length):
    '''
    Extract the sequences covered by all/most isoforms
    '''

    seg_lst = []
    for seg in intervals:
        # drop segments not covered by all/almost all isoforms
        if seg[3] >= coverage:
            # drop those that are shorter than indicated length
            seg_len = seg[2] - seg[1]
            if seg_len >= length:
                # extract whole segment sequence from the fasta file
                # (*BED-like coordinates*)
                seg_lst.append([seg,
                                str(fasta_index[seg[0]].seq[seg[1]:seg[2]])])
    seg_len = 0
    for e in seg_lst:
        seg_len += len(e[1])
    log.info('Processing {} nucleotides in {} segments.'.format(
             seg_len,
             len(seg_lst)))
    return seg_lst


class SequenceFilter(object):
    '''
    Filter sequence based on selected patterns

    Specify middle string or starting and ending sequence for filtering.

    TODO: re.search with '^' VS match
    '''

    def __init__(self, mid=None, start=None, end=None,
                 exclusive=True,
                 reject=True):

        beginning = r'^'
        middle = ''
        ending = ''
        if start:
            beginning = r'^{}'.format(start.upper())
        if mid:
            if exclusive:
                middle = (r'(?!{})'.format(r'|'.join(
                    [r'.*' + e.upper() + r'.+' for e in mid])))
            else:
                for pattern in mid:
                    middle = r''.join(
                        [r'(?=.*{}.+)'.format(e.upper()) for e in mid])
        if end:
            ending = r'.*{}$'.format(end.upper())

        self.pattern = re.compile(beginning + middle + ending)
        self.exclusive = exclusive
        self.reject = reject

    def filter(self, sequence):
        if self.reject:
            if self.pattern.search(sequence.upper()):
                return True
            else:
                return False
        else:
            if self.pattern.search(sequence.upper()):
                return False
            else:
                return True


def gen_seq(coord, sequence, length, strand, GC_lims, filters=None, index=0):
    '''
    Generate N long nucleotide stretches of CDS,
    not containing any RepeatMasker marked low complexity regions.

    Apply custom filters.
    '''

    while index <= len(sequence) - length:
        seq = sequence[index:index + length]
        # generate reverse complement, transcribed sequences,
        # depending on gene strand
        if strand == '+':
            rna_seq = transcribe(revcomp(seq))
        else:
            rna_seq = transcribe(seq)
        s_coords = [coord[0], coord[1] + index, coord[1] + index + length]
        index += 1
        # drop the sequence when outside of the GC content limits
        gc_content = round(gc(rna_seq), 2)
        if gc_content < GC_lims[0] or gc_content > GC_lims[1]:
            continue
        # Apply the filters
        filtered = False
        if filters:
            for f in filters:
                if f.filter(rna_seq):
                    filtered = True
                    break
        if filtered:
            continue
        # do not generate sequences with soft masked nucleotides
        if rna_seq.isupper():
            yield [s_coords, rna_seq, gc_content]
    return


def count_seq(segments, length, strand, GC_lims, filters=None):
    '''
    Count the number of sequences for processing
    '''

    GC_high = 0
    GC_low = 0
    filtered = 0
    sm_droped = 0
    good = 0
    for coord, sequence in segments:
        index = 0
        while index <= len(sequence) - length:
            seq = sequence[index:index + length]
            if strand == '+':
                rna_seq = transcribe(revcomp(seq))
            else:
                rna_seq = transcribe(seq)
            index += 1
            # count the sequence when outside of the GC content limits
            gc_content = round(gc(rna_seq), 2)
            if gc_content < GC_lims[0]:
                GC_low += 1
                continue
            if gc_content > GC_lims[1]:
                GC_high += 1
                continue
            flt = False
            if filters:
                for f in filters:
                    if f.filter(rna_seq):
                        flt = True
                        break
            if flt:
                filtered += 1
                continue
            # count sequences with soft masked nucleotides
            if rna_seq.isupper():
                good += 1
            else:
                sm_droped += 1
    total = GC_high + GC_low + sm_droped + filtered + good
    log.info('Valid sequences: {} / {} (GC_low = {}, GC_high = {}, \
filtered = {}, sm = {})'.format(good,
                                total,
                                GC_low,
                                GC_high,
                                filtered,
                                sm_droped))
    return good


def run_RNAcofold(seq1, seq2):
    '''
    Run RNAcofold to estimate the duplex folding change in free energy
    '''

    in_str = seq1 + '&' + seq2
    params = ['-a0',
              '-d2',
              '--noPS',
              '--output-format=D',
              '--csv-noheader']

    proc = subprocess.run(['RNAcofold'] + params,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          input=in_str,
                          encoding='ascii')
    if proc.returncode != 0:
        log.error('Subprocess exception: ' + proc.stderr)
        proc.check_returncode()
    else:
        result = proc.stdout.split(',')
        en = {'dG': Decimal(result[7]),
              'AB': Decimal(result[8]),
              'AA': Decimal(result[9]),
              'BB': Decimal(result[10]),
              'A': Decimal(result[11]),
              'B': Decimal(result[12])}
        return en


def get_tmpfs():
    '''
    Helper function to to pick tmpfs folder and check,
    if it exists on different OSes.

    Might need to be developed further.
    '''

    tmp_dir = os.path.join('/dev', 'shm')
    if not os.path.exists(tmp_dir):
        tmp_dir = os.path.join('/run', 'shm')
    try:
        if not os.path.exists(tmp_dir):
            raise Exception('Directory does not exist.')
    except Exception as e:
        log.exception('TMPFS directory does not exist.')
        raise
    return tmp_dir


def run_RNAfold(transcript_id, transcript_seq):
    '''
    Run RNAfold from ViennaRNA Package.

    Discard its output, feed the Dot Plot file to the mountain.pl
    script to recover the positional entropy for each base of the
    input sequence.
    '''

    params = ['-p1',
              '-d2',
              '--noLP',
              '--noPS']

    in_str = '>{}\n'.format(transcript_id) + transcript_seq

    tmp_dir = get_tmpfs()

    out_fn = os.path.join(tmp_dir, transcript_id + '_dp.ps')
    if not os.path.exists(out_fn):
        proc1 = subprocess.run(['RNAfold'] + params,
                               stdout=subprocess.DEVNULL,
                               stderr=subprocess.PIPE,
                               input=in_str,
                               encoding='ascii',
                               cwd=tmp_dir)
        if proc1.returncode != 0:
            log.error('Subprocess exception: ' + proc1.stderr)
            proc1.check_returncode()

    try:
        if not os.path.exists(out_fn):
            raise Exception('File does not exist.')
    except Exception as e:
        log.exception('Output file does not exist.')
        raise

    proc2 = subprocess.run(['mountain.pl', out_fn],
                           stderr=subprocess.PIPE,
                           stdout=subprocess.PIPE,
                           cwd=tmp_dir)
    if proc2.returncode != 0:
        log.error('Subprocess exception: ' + proc2.stderr)
        proc2.check_returncode()
    else:
        data = proc2.stdout.decode('ascii').split('\n')
        elements = [e.split() for e in data]
        mark = []
        for i, position in enumerate(elements):
                if position != [] and position[0] == '&':
                    mark.append(i)
        energy = [Decimal(v) for k, v in elements[mark[1] + 1:-1]]
        return energy


def get_trans(database, fasta_index, gene):
    '''
    Obtain transcript sequences

    TODO: with biopython this is slow, as it accesses the fasta file
    from the beginning of each record for each exon of every transcript
    and needs a future rewrite

    For speed, use:
    blastdbcmd -entry [tr.id, tr['transcript_version']] -db danRer_e91_allrna
    '''

    for transcript in database.children(gene,
                                        featuretype='transcript',
                                        order_by='start'):
        if transcript['transcript_biotype'][0] == 'protein_coding':
            trans_seq = None
            for exon in database.children(transcript,
                                          featuretype='exon',
                                          order_by='start'):
                if not trans_seq:
                    trans_seq = str(fasta_index[transcript.chrom]
                                    .seq[exon.start - 1:exon.end])
                else:
                    trans_seq += str(fasta_index[transcript.chrom]
                                     .seq[exon.start - 1:exon.end])
            if transcript.strand == '-':
                trans_seq = revcomp(trans_seq)
            yield transcript.id, trans_seq
    return


def pick_transcript(database, fasta_index, gene=None, eid=None):
    '''
    Pick only the longest from protein coding, preferably havana/ensembl+havana
    annotated transcripts and return it's id and sequence;
    if Ensembl id is provided by user return sequence of that particular one.
    '''

    if not eid:
        log.info('No transcript Ensemble ID specified.\nPicking the longest, \
preferably HAVANA annotated transcript for free energy calculations.')
        max_len = 0
        max_id = None
        max_exons = None
        max_source = None
        for transcript in database.children(gene,
                                            featuretype='transcript',
                                            order_by='start'):
            if transcript['transcript_biotype'][0] == 'protein_coding':
                t_len = 0
                t_id = transcript.id
                t_exo = []
                t_sou = transcript['transcript_source'][0]
                for exon in database.children(transcript,
                                              featuretype='exon',
                                              order_by='start'):
                    # calculate the transcript length from it's exons' length,
                    # store the exons' coords,
                    # do not extract the sequence yet
                    t_len += exon.end - exon.start + 1
                    t_exo.append((exon.start - 1, exon.end))
                if t_len > max_len and ((max_source == 'ensembl'
                                         or max_source is None)
                                        or (t_sou == 'havana'
                                            or t_sou == 'ensembl_havana')):
                    max_len = t_len
                    max_id = t_id
                    max_exons = t_exo
                    max_source = t_sou
                elif max_source == 'ensembl'\
                        and (t_sou == 'havana'
                             or t_sou == 'ensembl_havana'):
                    # prefer havana annotated transcripts
                    max_len = t_len
                    max_id = t_id
                    max_exons = t_exo
                    max_source = t_sou
                elif t_len == max_len and t_sou != max_source:
                    # if by chance two transcripts are of equal length,
                    # pick the one with merged annotation;
                    # if both have the same annotation source,
                    # keep the first one
                    if t_sou == 'ensembl_havana':
                        max_len = t_len
                        max_id = t_id
                        max_exons = t_exo
                        max_source = t_sou
                elif t_len == max_len:
                    log.info('W: Discarding transcript {}, source: {}, \
of equal length: {} to {}'.format(
                             t_id,
                             t_sou,
                             t_len,
                             max_id))
        log.info('Target: {}, source: {}, length: {}'.format(
                 max_id,
                 max_source,
                 max_len))
        max_seq = None
        for exon_coord in max_exons:
            if not max_seq:
                max_seq = str(fasta_index[gene.chrom]
                              .seq[exon_coord[0]:exon_coord[1]])
            else:
                max_seq += str(fasta_index[gene.chrom]
                               .seq[exon_coord[0]:exon_coord[1]])
        if gene.strand == '-':
            max_seq = revcomp(max_seq)
        return max_id, max_seq
    else:
        transcript = database[eid]
        trans_seq = None
        for exon in database.children(transcript,
                                      featuretype='exon',
                                      order_by='start'):
            if not trans_seq:
                trans_seq = str(fasta_index[transcript.chrom]
                                .seq[exon.start - 1:exon.end])
            else:
                trans_seq += str(fasta_index[transcript.chrom]
                                 .seq[exon.start - 1:exon.end])
        if transcript.strand == '-':
            trans_seq = revcomp(trans_seq)
        log.info('Target: {}, source: {}, length: {}'.format(
                 transcript.id,
                 transcript['transcript_source'][0],
                 len(trans_seq)))
        return transcript.id, trans_seq


def create_negseqidlst(database, gene_id=None, transcript_id=None):
    '''
    Create temporary file with transcript id's to exclude from blast searches.

    When transcript Id is given, list all transcripts of the parent gene.
    '''

    g_id = gene_id
    if transcript_id:
        g_id = database[transcript_id]['gene_id'][0]

    gene = database[g_id]

    ngsi = ''
    for i, f in enumerate(database.children(gene, featuretype='transcript')):
        if i == 0:
            ngsi += '.'.join([f.id, f['transcript_version'][0]])
        else:
            ngsi += '\n' + '.'.join([f.id, f['transcript_version'][0]])

    tmp_fn = g_id + '.balst.ngsi.lst'
    tmp_file = os.path.join(get_tmpfs(), tmp_fn)
    with open(tmp_file, 'w') as handle:
        handle.write(ngsi)
    return tmp_file


def blast_it(sequence, db, tmp_file=None):
    '''
    Run blastn on a sequence (wordsize 7).
    Return Bitscore, nident

    *Might* be faster to run on all sequences beforehand.
    '''

    try:
        if not isinstance(sequence, str):
            raise Exception('Wrong type')
    except Exception as e:
        log.exception('Input not a string')
        raise

    params = ['-task', 'blastn',
              '-evalue', '1000',
              '-word_size', '7',
              '-db', db,
              '-outfmt', '6 bitscore nident',
              '-num_threads', '1',
              '-max_target_seqs', '1']
    if tmp_file:
        params += ['-negative_seqidlist', tmp_file]

    proc = subprocess.run(['blastn',
                           '-query', '-',
                           '-out', '-'] + params,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          input=sequence,
                          encoding='ascii')

    if proc.returncode != 0:
        print(proc.stderr)
        proc.check_returncode()
    else:
        b_out = proc.stdout.split()
        if len(b_out) == 0:
            bitscore = -1
            nident = 0
            return bitscore, nident
        else:
            bitscore = float(b_out[0])
            nident = int(b_out[1])
            return bitscore, nident


def feed_fun(segments,
             length,
             strand,
             GC_lims,
             transcript_seq,
             entropy,
             filters=None,
             blastdb=None,
             tmp_path=None):
    '''
    Returns an iterator for multiprocessing function.
    '''

    for coordinates, sequence in segments:
        for result in gen_seq(coordinates, sequence, length,
                              strand, GC_lims, filters, index=0):
            yield result + [transcript_seq] + [tmp_path] + [entropy] +\
                    [blastdb]
    return


def processing_fun(input_list):
    '''
    Processing function for multiprocessing pool.
    '''

    coords, c_seq, gc_content, t_seq, ngsi_tmp, entropy, blastdb = input_list
    # get the mean entropy of the target sequence
    # TODO: !!! This is far from perfect
    # t_seq.upper() required in case of working with user-input
    # and targeting soft masked sequence
    c_start = t_seq.upper().find(revcomp(revtranscribe(c_seq)))
    try:
        if c_start == -1:
            raise Exception('Pattern not found')
    except Exception as e:
        log.exception('Target sequence not found: {}'.format(
                      revcomp(revtranscribe(c_seq))))
        raise
    c_pos_ent = entropy[c_start:c_start + len(c_seq)]
    c_mean_ent = sum(c_pos_ent) / len(c_pos_ent)
    # estimate also the change in free energy of monomer binding to itself
    energy = run_RNAcofold(c_seq, c_seq)
    dG_AA = energy['dG']
    G_A = energy['A']
    # blast the sequence against the RNA database
    bs, ni = 0, 0
    if blastdb:
        bs, ni = blast_it(c_seq, blastdb, ngsi_tmp)
    return coords, c_seq, gc_content, c_mean_ent, dG_AA, G_A, bs, ni


def estimate_energy(database,
                    fasta_index,
                    gene,
                    intervals,
                    coverage,
                    length,
                    GC_lims,
                    blastdb=None,
                    ensembl_id=None,
                    filters=None,
                    proc=1,
                    verbose=True,
                    cleanup=True):
    '''
    Estimate free energy change for a single transcript
    '''

    pool = Pool(proc)
    strand = gene.strand
    segments = get_seq(fasta_index, intervals, coverage, length)
    # do not run if there are no segments to process
    if len(segments) == 0:
        return None
    total = count_seq(segments, length, strand, GC_lims, filters)
    ngsi_tmp = create_negseqidlst(database, gene_id=gene.id)
    result_list = []
    transcript_id, transcript_seq = pick_transcript(database,
                                                    fasta_index,
                                                    gene,
                                                    ensembl_id)
    if total != 0:
        # get the positional entropy for the whole sequence/transcript
        entropy = run_RNAfold(transcript_id, transcript_seq)
    else:
        log.info('No valid sequences could be generated.')
        return None

    if not blastdb:
        log.info('No BLAST database path specified')

    iterator = feed_fun(segments,
                        length,
                        strand,
                        GC_lims,
                        transcript_seq,
                        entropy,
                        filters,
                        blastdb,
                        ngsi_tmp)

    done = 0
    for result in pool.imap_unordered(processing_fun, iterator,
                                      chunksize=1):
        if verbose:
            done += 1
            progress = str(done) + '/' + str(total)
            percent = str(round(100 * Decimal(done) / Decimal(total), 1)) + '%'
            out = '{} ({})'.format(progress, percent)
            print(out, end='', flush=True)
            print('\r', end='')
        result_list.append(result)

    psfile = os.path.join(get_tmpfs(), transcript_id + '_dp.ps')
    if cleanup:
        # clean up
        # blast negative seqid list file
        os.remove(ngsi_tmp)
        # RNAfold dot plot PS file
        os.remove(psfile)
    else:
        log.info('Leaving the intermediate files:\n\
{}'.format('\n'.join(e for e in [ngsi_tmp, psfile] if e)))
    return result_list


def parse_fasta(file_name):
    '''
    Return a generator object yielding sequences from the FASTA file

    Joins sequences separated by new lines (only yields on '>'),
    deletes whitespaces within the sequences.
    '''

    with open(file_name, 'r') as handle:
        entry = ''
        new_entry = 1
        for l in handle:
            line = l.strip('\n ')
            if line.startswith('>') and new_entry == 1:
                # start new entry
                entry = line + '\n'
                new_entry = 0
            elif len(line) != 0 and not line.startswith('>') \
                    and new_entry == 0:
                # extend entry, delete all the whitespaces
                entry += line.replace(' ', '')
            elif line.startswith('>') and new_entry == 0:
                # yield previous entry, start next one
                yield entry
                entry = line + '\n'
        # yield the last entry or throw an error if format was wrong
        try:
            if not new_entry == 0:
                raise Exception('Wrong input format')
        except Exception as e:
            log.exception('Sequence not in FASTA format.')
            raise
        yield entry
    return


def processing_input(string_in):
    '''
    Process the input sequence
    Return the elements needed for further steps.
    '''

    # quick checks
    try:
        if len(string_in) == 0:
            raise Exception('Input empty')
    except Exception as e:
        log.exception('Input empty')
        raise
    # take care of white spaces
    sequence = string_in.strip()
    try:
        if not sequence.upper().startswith(('>', 'A', 'G', 'C', 'T')):
            raise Exception('Wrong input format')
    except Exception as e:
        log.exception('Input not a DNA sequence in FASTA format or PLAIN.')
        raise

    # if fasta, get the sequence id from the header
    if sequence.startswith('>'):
        id_line, seq = sequence.split('\n', maxsplit=1)
        # strip it of all non-alphanumeric characters (allow '_')
        pattern = re.compile(r'\W+')
        # first element is empty (because of '>'), second is SeqID
        seq_id = pattern.split(id_line, maxsplit=2)[1]
    else:
        seq_id = 'plain'
        seq = sequence

    # join all the lines if present, also taking care of white spaces
    seq = ''.join(l.strip() for l in seq.split('\n'))

    return seq_id, seq


def check_segid(database, input_id):
    '''
    Helper function to check if the sequence ID is valid
    for negative seqid list file creation.

    Will have to be developed further.
    '''

    try:
        feature = database[input_id]
    except gffutils.FeatureNotFoundError:
        log.warning('Feature \'{}\' not found in the database.'.format(
            input_id))
        return False
    else:
        if feature.featuretype == 'transcript':
            return True
        else:
            log.warning('\'{}\' is not a valid Transcript ID.'.format(
                input_id))
            return False


def estimate_energy_input(input_sequence,
                          length,
                          GC_lims,
                          strand='+',
                          database=None,
                          fasta_index=None,
                          blastdb=None,
                          filters=None,
                          proc=1,
                          verbose=True,
                          cleanup=True):
    '''
    Estimate free energy change for a single input sequence.

    Parameters:
    -----------

    strand- use if input is genomic sequence
    database; fasta_index - must be provided if database Transcript ID
                            lookup is required
    '''

    pool = Pool(proc)
    input_id, input_seq = processing_input(input_sequence)
    # assume that for the input the strand is always the Watson one
    segment = [[['.', 0, len(input_seq)], input_seq]]
    total = count_seq(segment, length, strand, GC_lims, filters)
    ngsi_tmp = None
    # check if the input ID is in the database and if it is a valid
    # transcript ID, if yes use it to create negative segid list
    # for BLAST and fold this particular transcript for energy estimations
    if database and fasta_index and check_segid(database, input_id):
        ngsi_tmp = create_negseqidlst(database, transcript_id=input_id)
        transcript_id, transcript_seq = pick_transcript(database,
                                                        fasta_index,
                                                        eid=input_id)
    else:
        log.info('Using the provided sequence for free energy calculation.')
        transcript_id, transcript_seq = input_id, input_seq

    result_list = []
    if total != 0:
        # get the positional entropy for the whole sequence/transcript
        entropy = run_RNAfold(transcript_id, transcript_seq)
    else:
        log.info('No valid sequences could be generated.')
        return input_id, None

    if not blastdb:
        log.info('No BLAST database path specified')

    iterator = feed_fun(segment,
                        length,
                        strand,
                        GC_lims,
                        transcript_seq,
                        entropy,
                        filters,
                        blastdb,
                        ngsi_tmp)

    done = 0
    for result in pool.imap_unordered(processing_fun, iterator,
                                      chunksize=1):
        if verbose:
            done += 1
            progress = str(done) + '/' + str(total)
            percent = str(round(100 * Decimal(done) / Decimal(total), 1)) + '%'
            out = '{} ({})'.format(progress, percent)
            print(out, end='', flush=True)
            print('\r', end='')
        result_list.append(result)

    psfile = os.path.join(get_tmpfs(), transcript_id + '_dp.ps')
    if cleanup:
        # clean up
        # blast negative seqid list file
        if ngsi_tmp:
            os.remove(ngsi_tmp)
        # RNAfold dot plot PS file
        os.remove(psfile)
    else:
        log.info('Leaving the intermediate files:\n\
{}'.format('\n'.join(e for e in [ngsi_tmp, psfile] if e)))
    return input_id, result_list


def gene2csm(database,
             fasta_index,
             var_db,
             target_lst,
             crRNA_lenght,
             GC_limit,
             blastdb=None,
             n_threads=1,
             filters=None,
             coverage_limit='max',
             use_variation=True,
             exclude_dict=None,
             file_prefix=None,
             verbose=True,
             cleanup=True):
    '''
    Main function to run the program.
    '''

    result = []
    for target in target_lst:
        log.info('\n---\n')
        gene = database[target]
        # important assertion- gene has a specified strand
        try:
            if not (gene.strand == '+' or gene.strand == '-'):
                raise Exception('Wrong symbol')
        except Exception as e:
            log.exception('Target gene has unspecified strand symbol: \
{}.'.format(
                          gene.strand))
            raise
        g_name = gene['gene_name'][0]
        log.info(g_name)
        # if present get an exclusion list for this target
        e_list = None
        if exclude_dict:
            if target in exclude_dict:
                e_list = exclude_dict[target]
        gene_cov = get_cov(database, gene)
        cov_lims = coverage_limit
        if coverage_limit == 'max':
            cov_lims = gene_cov.max()
        if gene.strand == '+':
            plot_cov = gene_cov
        else:
            plot_cov = gene_cov[::-1]
        # plotting
        plt.plot(plot_cov, 'blue')
        plt.show()
        if use_variation:
            gen_int = sub_var(sub_usr(get_int(gene, gene_cov),
                                      gene.strand,
                                      e_list),
                              var_db)
        else:
            gen_int = sub_usr(get_int(gene, gene_cov),
                              gene.strand,
                              e_list)
        output = estimate_energy(database,
                                 fasta_index,
                                 gene,
                                 gen_int,
                                 cov_lims,
                                 crRNA_lenght,
                                 GC_limit,
                                 blastdb=blastdb,
                                 filters=filters,
                                 proc=n_threads,
                                 verbose=verbose,
                                 cleanup=cleanup)
        # do not store the empty results
        if not output:
            log.warning('No valid segments for target {}:{}'.format(target,
                                                                    g_name))
            continue

        if file_prefix:
            fn = file_prefix + target + '.result.csv'
        else:
            fn = target + '.result.csv'
        log.info('\nWriting file {}.'.format(fn))
        with open(fn, 'w') as handle:
            for item in output:
                handle.write(','.join([str(e) for e in item]) + '\n')

        result.append((g_name, output))
        log.info('\nDone.')
    return result


def conf_logger(logger):
    '''
    Configure logger
    '''

    class SingleLevelFilter(logging.Filter):
        '''
        Logging filter to single out specific log levels.
        '''

        def __init__(self, passlevel, reject):
            self.passlevel = passlevel
            self.reject = reject

        def filter(self, record):
            if self.reject:
                return (record.levelno != self.passlevel)
            else:
                return (record.levelno == self.passlevel)

    # do not set the logging level - let it be set by the parent module
    # log.setLevel(logging.DEBUG)
    # set propagate to False, so the message wont be propagated to the ancestor
    # loggers if they get initiated
    logger.propagate = False
    # create console handler for all but info levels
    ch = logging.StreamHandler(sys.stdout)
    # ch.setLevel(logging.DEBUG)
    # create filter for rejection of info level
    f = SingleLevelFilter(logging.INFO, True)
    ch.addFilter(f)
    # create formatters for this handler
    basic_formatter = logging.Formatter(
        fmt='{asctime}:{name}:{levelname}: {message}',
        style='{')
    # add it to the handler
    ch.setFormatter(basic_formatter)
    # now for info
    chi = logging.StreamHandler(sys.stdout)
    # chi.setLevel(logging.INFO)
    fi = SingleLevelFilter(logging.INFO, False)
    chi.addFilter(fi)
    info_formatter = logging.Formatter(fmt='{message}',
                                       datefmt=None,
                                       style='{')
    chi.setFormatter(info_formatter)
    # add the handlers to the logger
    logger.addHandler(ch)
    logger.addHandler(chi)
    return logger


log = conf_logger(logging.getLogger(__name__))

# TODO:
# Better handling of finding sequences within the folded transcript
# * pass the exon coordinates
# * compute the exon lengths
# * compute the sequence relative coordinates within the transcript
#
# things to consider when doing a rewrite
# -- minor gain would be to run RNAcofold -p instead of -a
# and not calculate the dG of AA complex formation

# round it to 2 decimals
#    dG_AA = round(cofold_dict['AA'] - cofold_dict['A'] * 2, 2)
