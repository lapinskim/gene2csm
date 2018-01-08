# Imports
import gffutils
import os
import numpy as np
import matplotlib.pyplot as plt
import subprocess
from multiprocessing import Pool
from collections import Counter
from decimal import Decimal
from pybedtools import BedTool

# Define functions


def revcomp(seq):
    # reverse complement sequence

    assert isinstance(seq, str)
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
    # return GC content percentage

    assert isinstance(seq, str)
    seq = seq.upper()
    gc_content = 100 * Decimal(seq.count("G") + seq.count("C")) / len(seq)
    return gc_content


def transcribe(seq):
    # transcribe sequence- change T to U

    assert isinstance(seq, str)
    trans_seq = seq.replace('T', 'U').replace('t', 'u')
    return trans_seq


def revtranscribe(seq):
    # reverse transcribe sequence- change U to T

    assert isinstance(seq, str)
    trans_seq = seq.replace('U', 'T').replace('u', 't')
    return trans_seq


def get_cov(database, gene):
    # Create genomic sequence coverage by CDS (without STOP codon).

    # get the length of the gene (1 based coords)
    gene_length = gene.end - gene.start + 1
    print('{} loci length: {!s}.'.format(gene.id, gene_length))
    # create a numpy array of zeros with the length of the gene
    gene_cov = np.zeros(gene_length, dtype=int)
    for cds in database.children(gene, featuretype='CDS', order_by='start'):
            # get relative exon coords in 0-based space
            rel_start = cds.start - gene.start
            rel_end = cds.end - gene.start
            gene_cov[rel_start:rel_end + 1] += 1
    print('Maximal genomic sequence coverage by CDS: {!s}.'
          .format(gene_cov.max()))
    return gene_cov


def get_int(gene, genomic_coverage):
    # get the genomic intervals (for sequence retrieval)

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
    print('{!s} CDS coverage intervals (coverage, count): {!s}.'
          .format(len(gen_int), i_cnt.most_common()))
    return gen_int


def sub_user(intervals, strand, user_list=None):
    '''
    Subtract user submitted *CDS containing exons* from processing.
    '''

    i_list = intervals
    if user_list is not None:
        # check if exon numbers are within the exon range of the target
        for e in user_list:
            assert 1 <= e <= len(intervals), 'Number {} not within the target\'s \
CDS containing exon range: 1 - {}.'.format(e, len(intervals))
        # convert to set, to delete, numbered exon only (unique)
        # and reverse sort to delete from the end of the list
        sorted_list = sorted(set(user_list), reverse=True)
        if strand == "-":
            for e in sorted_list:
                del i_list[-e]
        elif strand == "+":
            for e in sorted_list:
                del i_list[e - 1]
        print('Omitting user submitted *CDS containing exons*: {}.'
              .format(', '.join(str(e) + '.' for e in sorted_list[::-1])))
    return i_list


def sub_var(intervals, variation_fn):
    '''
    Subtract variable sequences from the genomic intervals

    For future enhancement consider soft masking the whole genome
    beforehand or splitting the GVF file on chromosomes for speed improvements.
    '''

    # make some basic checks
    assert intervals != [], 'No intervals to process.'
    assert os.path.exists(variation_fn), 'GVF file does not exist.'
    assert variation_fn.endswith(('gvf', 'gvf.gz')), 'Ensure the variation \
containing file is in a GVF format. \'.gvf\' extension missing.'

    gvf_sorted_fn = 'sorted.gvf'.join(variation_fn.split('gvf'))
    # sort the variation database for faster processing and store
    # for further use
    if not os.path.exists(gvf_sorted_fn):
        print('Sorting and saving the new GVF file as {}.'
              .format(gvf_sorted_fn))
        gvf_bt = BedTool(variation_fn).sort().saveas(gvf_sorted_fn)
    else:
        gvf_bt = BedTool(gvf_sorted_fn)

    print('Subtracting variation from {}.'
          .format(gvf_bt[1]['Dbxref'].split(':')[0]))

    # make sure the format is right
    for e in intervals:
        assert len(e) == 4
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
    print('Processing {!s} nucleotides in {!s} segments.'
          .format(seg_len, len(seg_lst)))
    return seg_lst


def gen_seq(coord, sequence, length, strand, GC_lims, index=0):
    '''
    Generate N long nucleotide stretches of CDS,
    not containing any RepeatMasker marked low complexity regions
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
        # do not generate sequences with soft masked nucleotides
        if rna_seq.isupper():
            yield [s_coords, rna_seq, gc_content]
    return


def count_seq(segments, length, GC_lims):
    '''
    Count the number of sequences for processing
    '''

    GC_high = 0
    GC_low = 0
    sm_droped = 0
    good = 0
    for coord, sequence in segments:
        index = 0
        while index <= len(sequence) - length:
            seq = sequence[index:index + length]
            index += 1
            # count the sequence when outside of the GC content limits
            gc_content = round(gc(seq), 2)
            if gc_content < GC_lims[0]:
                GC_low += 1
                continue
            if gc_content > GC_lims[1]:
                GC_high += 1
                continue
            # count sequences with soft masked nucleotides
            if seq.isupper():
                good += 1
            else:
                sm_droped += 1
    total = GC_high + GC_low + sm_droped + good
    print('Valid sequences: {} / {} (GC_low = {}, GC_high = {}, sm = {})'
          .format(good,
                  total,
                  GC_low,
                  GC_high,
                  sm_droped))
    return good


def run_RNAcofold(seq):
    '''
    Run RNAcofold to estimate the duplex folding change in free energy
    '''

    in_str = seq + '&' + seq
    params = ['-a0',
              '-d2',
              '--noLP',
              '--noPS',
              '--output-format=D',
              '--csv-noheader']

    proc = subprocess.run(['RNAcofold'] + params,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          input=in_str,
                          encoding='ascii')
    if proc.returncode != 0:
        print(proc.stderr)
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
    assert os.path.exists(tmp_dir)
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

    proc1 = subprocess.run(['RNAfold'] + params,
                           stdout=subprocess.DEVNULL,
                           stderr=subprocess.PIPE,
                           input=in_str,
                           encoding='ascii',
                           cwd=tmp_dir)
    if proc1.returncode != 0:
        print(proc1.stderr)
        proc1.check_returncode()

    out_fn = os.path.join(tmp_dir, transcript_id + '_dp.ps')
    assert os.path.exists(out_fn)

    proc2 = subprocess.run(['mountain.pl', out_fn],
                           stderr=subprocess.PIPE,
                           stdout=subprocess.PIPE,
                           cwd=tmp_dir)
    if proc2.returncode != 0:
        print(proc2.stderr)
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
        print('No transcript Ensemble ID specified.\nPicking the longest,\
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
                    print('W: Discarding transcript {}, source: {}, \
of equal length: {!s} to {}'.format(t_id,
                                    t_sou,
                                    t_len,
                                    max_id))
        print('Target: {}, source: {}, length: {!s}'
              .format(max_id,
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
        print('Target: {}, source: {}, length: {!s}'
              .format(transcript.id,
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


def blast_it(sequence, tmp_file=None):
    '''
    Run blastn on a sequence (wordsize 7).
    Return Bitscore, nident

    *Might* be faster to run on all sequences beforehand.
    '''

    # TODO: take out the blastdb parameter

    assert isinstance(sequence, str)

    params = ['-task', 'blastn',
              '-word_size', '7',
              '-db', './blastdb/danRer_e91_allrna',
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
             tmp_path=None):
    '''
    Returns an iterator for multiprocessing function.
    '''

    for coordinates, sequence in segments:
        for result in gen_seq(coordinates, sequence, length,
                              strand, GC_lims, index=0):
            yield result + [transcript_seq] + [tmp_path] + [entropy]
    return


def processing_fun(input_list):
    '''
    Processing function for multiprocessing pool.
    '''

    coords, c_seq, gc_content, t_seq, ngsi_tmp, entropy = input_list
    # get the mean entropy of the target sequence
    # TODO: !!! This is far from perfect
    # t_seq.upper() required in case of working with user-input
    # and targeting soft masked sequence
    c_start = t_seq.upper().find(revcomp(revtranscribe(c_seq)))
    assert c_start != -1, 'Target sequence not found: {}'.format(
        revcomp(revtranscribe(c_seq)))
    c_pos_ent = entropy[c_start:c_start + len(c_seq)]
    c_mean_ent = sum(c_pos_ent) / len(c_pos_ent)
    # estimate also the change in free energy of monomer binding to itself,
    energy = run_RNAcofold(c_seq)
    dG_AA = energy['dG']
    G_A = energy['A']
    # blast the sequence against the Danio rerio RNA database
    bs, ni = blast_it(c_seq, ngsi_tmp)
    return coords, c_seq, gc_content, c_mean_ent, dG_AA, G_A, bs, ni


def estimate_energy(database,
                    fasta_index,
                    gene,
                    intervals,
                    coverage,
                    length,
                    GC_lims,
                    ensembl_id=None,
                    proc=1):
    '''
    Estimate free energy change for a single transcript
    '''

    pool = Pool(proc)
    strand = gene.strand
    segments = get_seq(fasta_index, intervals, coverage, length)
    # do not run if there are no segments to process
    if len(segments) == 0:
        return -1
    total = count_seq(segments, length, GC_lims)
    ngsi_tmp = create_negseqidlst(database, gene_id=gene.id)
    result_list = []
    transcript_id, transcript_seq = pick_transcript(database,
                                                    fasta_index,
                                                    gene,
                                                    ensembl_id)
    # get the positional entropy for the whole transcript
    entropy = run_RNAfold(transcript_id, transcript_seq)

    iterator = feed_fun(segments,
                        length,
                        strand,
                        GC_lims,
                        transcript_seq,
                        entropy,
                        ngsi_tmp)

    done = 0
    for result in pool.imap_unordered(processing_fun, iterator,
                                      chunksize=1):
        done += 1
        progress = str(done) + '/' + str(total)
        percent = str(round(100 * Decimal(done) / Decimal(total), 1)) + '%'
        out = '{} ({})'.format(progress, percent)
        print(out, end='', flush=True)
        print('\r', end='')
        result_list.append(result)

    # clean up
    # blast negative seqid list file
    os.remove(ngsi_tmp)
    # RNAfold dot plot PS file
    os.remove(os.path.join(get_tmpfs(), transcript_id + '_dp.ps'))
    return result_list


def processing_input(string):
    '''
    Process the input sequence
    Return the elements needed for further steps.
    '''

    # quick checks
    assert isinstance(string, str), 'Input not a string.'
    assert len(string) != 0, 'Input empty.'
    # take care of white spaces
    sequence = string.strip()
    assert sequence.upper().startswith(('>', 'A', 'G', 'C', 'T')),\
        'Input not a DNA sequence in FASTA format or PLAIN.'

    # if fasta, get the sequence id from the header
    if sequence.startswith('>'):
        id_line, seq = sequence.split('\n', maxsplit=1)
        seq_id = id_line.lstrip('>').split()[0]
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
        print('Feature \'{}\' not found in the database.'.format(input_id))
        return False
    else:
        if feature.featuretype == 'transcript':
            return True
        else:
            print('\'{}\' is not a valid Transcript ID.'.format(input_id))
            return False


def estimate_energy_input(input_sequence,
                          length,
                          GC_lims,
                          strand='+',
                          database=None,
                          fasta_index=None,
                          proc=1):
    '''
    Estimate free energy change for a single transcript.

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
    total = count_seq(segment, length, GC_lims)
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
        print('Using the provided sequence for free energy calculation.')
        transcript_id, transcript_seq = input_id, input_seq

    result_list = []
    # get the positional entropy for the whole sequence/transcript
    entropy = run_RNAfold(transcript_id, transcript_seq)

    iterator = feed_fun(segment,
                        length,
                        strand,
                        GC_lims,
                        transcript_seq,
                        entropy,
                        ngsi_tmp)

    done = 0
    for result in pool.imap_unordered(processing_fun, iterator,
                                      chunksize=1):
        done += 1
        progress = str(done) + '/' + str(total)
        percent = str(round(100 * Decimal(done) / Decimal(total), 1)) + '%'
        out = '{} ({})'.format(progress, percent)
        print(out, end='', flush=True)
        print('\r', end='')
        result_list.append(result)

    # clean up
    # blast negative seqid list file
    if ngsi_tmp:
        os.remove(ngsi_tmp)
    # RNAfold dot plot PS file
    os.remove(os.path.join(get_tmpfs(), transcript_id + '_dp.ps'))
    return result_list


def gene2csm(database,
             fasta_index,
             var_db,
             target_lst,
             crRNA_lenght,
             GC_limit,
             n_threads,
             coverage_limit='max',
             e_list=None):
    '''
    Main function to run the program.
    '''

    result = []
    for target in target_lst:
        print('\n---\n')
        gene = database[target]
        # important assertion- has a specified strand
        assert gene.strand == '+' or gene.strand == '-', 'Target gene has \
unspecified strand symbol: {}.'.format(gene.strand)
        g_name = gene['gene_name'][0]
        print(g_name)
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
        gen_int = sub_var(sub_user(get_int(gene, gene_cov), gene.strand,
                                   e_list), var_db)
        output = estimate_energy(database,
                                 fasta_index,
                                 gene,
                                 gen_int,
                                 cov_lims,
                                 crRNA_lenght,
                                 GC_limit,
                                 proc=n_threads)
        # do not store the empty results
        if output == -1:
            print('No valid segments for target {}:{}'.foramt(target, g_name))
            continue

        fn = target + '.result.csv'
        print('\nWriting file {}.'.format(fn))
        with open(fn, 'w') as handle:
            for item in output:
                handle.write(','.join([str(e) for e in item]) + '\n')

        result.append((g_name, output))
        print('\nDone.')
    return result


# TODO: things to consider when doing a rewrite
# -- minor gain would be to run RNAcofold -p instead of -a
# and not calculate the dG of AA complex formation

# round it to 2 decimals
#    dG_AA = round(cofold_dict['AA'] - cofold_dict['A'] * 2, 2)
