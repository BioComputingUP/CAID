import os
from itertools import groupby


def parse_vertical_format(inputfile, kwargs):
    seq = ''
    state = ''
    scores = []
    with open(inputfile) as f:
        for i, line in enumerate(f):
            if i >= kwargs.get('header_lines'):
                if line[:3] != 'END': # For the dispred2 (Hoque) predictor
                    line = line.rstrip().split()
                    if kwargs.get('residue_pos') is not None:
                        seq += line[kwargs.get('residue_pos')]
                    if kwargs.get('score_pos') is not None:
                        scores.append(line[kwargs.get('score_pos')] if line[kwargs.get('score_pos')] != 'NA' else '0.0')
                    if kwargs.get('status_pos') is not None:
                        state += '1' if line[kwargs.get('status_pos')] == kwargs.get('status_char') else '0'
                    if kwargs.get('status_th') is not None:
                        state += '1' if float(line[kwargs.get('score_pos')]) >= kwargs.get('status_th') else '0'
    return seq if seq else None, state if state else None, scores if scores else None


def parse_predisorder(inputfile, kwargs):
    with open(inputfile) as f:
        seq = next(f).strip()
        state = ''.join(['1' if status == 'D' else '0' for status in next(f).strip()])
        scores = next(f).split()
    return seq, state, scores


def parse_anchor(inputfile, kwargs):
    seq = ''
    state = ''
    scores = []
    with open(inputfile) as f:
        for i, line in enumerate(f):
            if line[0] != '#':
                line=line.strip().split()
                seq += line[1]
                scores.append(line[2])
                state += line[3]
    return seq, state, scores


def parse_foldunfold(inputfile, kwargs):
    regions = []
    with open(inputfile) as f:
        for line in f:
            start, end = line.split(':')[1].split('--')
            regions.append((int(start), int(end)))
    state = ['0' for i in kwargs.get('seq')]
    for start, end in regions:
        for i in range(start -1, end):
            state[i] = '1'
    return None, ''.join(state), None


def parse_isunstruct(inputfile, kwargs):
    seq = ''
    state = ''
    scores = []
    with open(inputfile) as f:
        next(f)
        next(f)
        next(f)
        next(f)
        for i, line in enumerate(f):
            line = line.split()
            if len(line) == 5:
                _, l, status, _, score = line
            elif len(line) == 4:
                _, l, status, score = line
            seq += l
            scores.append(score)
            state += '1' if status == 'U' else '0'
    return seq, state, scores


# TODO check a second threshold can be applicable
def parse_globplot(inputfile, kwargs):
    with open(inputfile) as f:
        scores = eval(next(f))[0]['p']
    return None, ''.join(['1' if score >= 0.0 else '0' for score in scores]), scores


# TODO check some proteins are parsed wrongly (score alone in the last line)
def parse_dflpred(inputfile, kwargs):
    # Score threshold is 0.18 according to documenation
    with open(inputfile) as f:
        next(f)
        seq = next(f).strip()
        scores = next(f).strip().split(',')
    state = ['1' if l.isupper() else '0' for l in seq]
    return seq.upper(), ''.join(state), scores


def parse_disordpbind_rna(inputfile, kwargs):
    ele = _parse_disordpbind(inputfile, kwargs)
    return ele[0], ele[2], ele[3]


def parse_disordpbind_dna(inputfile, kwargs):
    ele = _parse_disordpbind(inputfile, kwargs)
    return ele[0], ele[4], ele[5]


def parse_disordpbind_prot(inputfile, kwargs):
    ele = _parse_disordpbind(inputfile, kwargs)
    return ele[0], ele[6], ele[7]


def parse_disordpbind_all(inputfile, kwargs):
    ele = _parse_disordpbind(inputfile, kwargs)
    return ele[0], ele[1], None


def _parse_disordpbind(inputfile, kwargs):
    # Uppercase means interact with RNA/DNA/protein
    with open(inputfile) as f:
        next(f)
        seq = next(f).strip()
        state_rna = next(f).split(':')[1].strip()
        scores_rna = next(f).split(':')[1].strip()[:-1].split(',')
        state_dna = next(f).split(':')[1].strip()
        scores_dna = next(f).split(':')[1].strip()[:-1].split(',')
        state_prot = next(f).split(':')[1].strip()
        scores_prot = next(f).split(':')[1].strip()[:-1].split(',')

    state = ['1' if l.isupper() else '0' for l in seq]
    return seq.upper(), state, state_rna, scores_rna, state_dna, scores_dna, state_prot, scores_prot


def parse_fmorfpred(inputfile, kwargs):
    # The second state row is just majority vote between iupred and espritz, so it is discarder
    with open(inputfile) as f:
        next(f)
        seq = next(f).strip()
        scores = next(f).strip().split(',')
        state = next(f).strip()
    return seq.upper(), state, scores


def parse_vsl2(inputfile, kwargs):
    seq = ''
    state = ''
    scores = []
    start = False
    with open(inputfile) as f:
        for i, line in enumerate(f):
            if line.strip() == '========================================':
                start = False
            if start:
                ele = line.strip().split()
                seq += ele[1]
                scores.append(ele[2].replace(',', '.'))
                state += '1' if ele[3] == 'D' else '0'
            if line.strip() == '----------------------------------------':
                start = True
    return seq, state, scores


def parse_dis465(inputfile, kwargs):
    scores_dis465, scores_disHL = _parse_disemble(inputfile, kwargs)
    state = ['1' if score >= 0.5 else '0' for score in scores_dis465]
    return None, ''.join(state), scores_dis465


def parse_disHL(inputfile, kwargs):
    rscores_dis465, scores_disHL = _parse_disemble(inputfile, kwargs)
    state = ['1' if score >= 0.086 else '0' for score in scores_disHL]
    return None, ''.join(state), scores_disHL


def _parse_disemble(inputfile, kwargs):
    with open(inputfile) as f:
        ele = eval(next(f))
    scores_dis465 = next(filter(lambda k: k['pred'] == 'dis465', ele))['p']
    # before: scores_dis465 = filter(lambda k: k['pred'] == 'dis465', ele)[0]['p']
    scores_disHL = next(filter(lambda k: k['pred'] == 'disHL', ele))['p']
    # before: scores_disHL = filter(lambda k: k['pred'] == 'disHL', ele)[0]['p']
    return scores_dis465, scores_disHL


def parse_mobidblite(inputfile, kwargs):
    with open(inputfile) as f:
        ele = eval(next(f))['predictions']
    ele = filter(lambda k:k['method'] == 'mobidb_lite', ele)

    if ele:
        first_ele = next(ele)
        state = ['0' for i in first_ele['scores']]
        # changed from: state = ['0' for i in ele[0]['scores']]
        for start, end, ann in first_ele['regions']:
        # changed from: for start, end, ann in ele[0]['regions']:
            if ann[0] == 'D':
                for i in range(start - 1, end):
                    state[i] = '1'
    else:
        state = None
        first_ele = None
    return None, ''.join(state) if state else None, first_ele['scores'] if first_ele else None
    # changed from: return None, ''.join(state) if state else None, ele[0]['scores'] if ele else None


def parse_mcl(inputfile, kwargs):
    seq, state, scores = parse_vertical_format(inputfile, kwargs)
    regions = _state2regions(state)
    regions = _filter_regions(regions, 3)
    state = _regions2state(regions, len(state))
    return seq, state, scores


# TODO check since it gives wrong format
def parse_mcw(inputfile, kwargs):
    seq, state, scores = parse_vertical_format(inputfile, kwargs)
    regions = _state2regions(state)
    regions = _filter_regions(regions, 3)
    state = _regions2state(regions, len(state))
    return seq, state, scores


def parse_jronn(inputfile, kwargs):
    block = []
    name = None
    jronn_data = {}
    with open(inputfile) as f:
        for line in f:
            if line.strip():
                line = line.strip().split()
                if line[0][0] == '>':
                    if name:
                        jronn_data[name] = _parse_jronn_block(block, kwargs)
                    name = line[0][1:]
                    block = []
                else:
                    block.append(line)

    jronn_data[name] = _parse_jronn_block(block, kwargs)
    return jronn_data


def _parse_jronn_block(block, kwargs):
    seq = ''
    state = ''
    scores = []
    for l, score in block:
        seq += l
        scores.append(score.replace(',', '.'))
        state += '1' if float(score.replace(',', '.')) >= kwargs.get('status_th') else '0'
    return seq, state, scores


def _state2regions(state):
    """
    Transform disorder/order string into a list of disorder regions (start/end)
    """
    regions = []
    start = 0
    end = 0
    for i in range(len(state) - 1):
        if state[i] != state[i + 1]:
            if state[i] == "1":
                regions.append((start, end))
            start = i + 1
        end += 1
    # Check on last char
    if state[i] == "1":
        regions.append((start, end))
    return regions


def _filter_regions(regions, len_cutoff=3):
    filtered_regions = []
    for start, end in regions:
        if (end - start + 1) >= len_cutoff:
            filtered_regions.append((start, end))
    return filtered_regions


def _regions2state(regions, state_len):
    state = ['0' for i in range(state_len)]
    for start, end in regions:
        for i in range(start, end + 1):
            # print i
            state[i] = '1'
    return ''.join(state)


# Parse reference sequences
def parse_reference_sequences(fastafile):
    r_seqs = {}
    with open(fastafile) as f:
        for line in f:
            if line[0] == '>':
                name = line.rstrip().split()[0][1:]
                r_seqs.setdefault(name, '')
            else:
                r_seqs[name] += line.strip()
    return r_seqs


def strip_split(string, expected=None):
    if string:
        s = string.strip('\n').rsplit('\t')
        if expected and len(s) != expected:
            # logging.error('expecting %i elements from split, got %i', expected, len(s))
            raise ValueError('expecting %i elements from split, got %i', expected, len(s))
        return s


def parse_prediction_file(predfile):

    with open(predfile) as f:
        faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
        for header in faiter:
            header = next(header).strip()[1:]
            try:
                a = next(faiter)
                positions, residues, scores, states = zip(*map(strip_split, a))
            except ValueError:
                # logging.error('error while parsing prediction %s', header)
                raise ValueError('error while parsing prediction %s', header)

            # pred_entry = PredictionEntry(positions, residues, scores, states)
            yield header, positions, residues, scores, states


def write_prsed_results(methods, ref_seqs, overwrite=True):

    # result_dir = "/mnt/projects"
    # parsed_dir = "/home/damiano/results_parsed"
    result_dir = "/home/damiano/Projects/caid/data/results_raw"
    parsed_dir_err = "/home/damiano/Projects/caid/data/results_parsed_err"
    parsed_dir = "/home/damiano/Projects/caid/data/results_parsed"

    # Parse jron separately
    inputfile = "{}/jronn/results_jronn.txt".format(result_dir)
    jronn_data = parse_jronn(inputfile, {'status_th': 0.5})

    for method in methods:
        # inputpath = "{}/CAID/2018_09/{}/results/{}".format(result_dir, method['group'], method['result_dir'])
        inputpath = "{}/{}".format(result_dir, method['result_dir'])
        errfile = "{}/{}_{}_{}{}.err".format(parsed_dir_err, method['id'], method['group'], method['result_dir'], method.get('label', ''))
        outfile = "{}/{}_{}_{}{}.out".format(parsed_dir, method['id'], method['group'], method['result_dir'], method.get('label', ''))
        if not os.path.isfile(outfile) or overwrite:
            with open(errfile, 'w') as ferr:
                if method['func'] is not None:
                    with open(outfile, 'w') as fout:

                        for disprot_id in ref_seqs:

                            ref_seq = ref_seqs[disprot_id]
                            method['seq'] = ref_seq

                            is_good = True
                            seq = None
                            state = None
                            scores = None

                            if method['result_dir'] == "jronn":
                                seq, state, scores = jronn_data.get(disprot_id, (None, None, None))

                            else:
                                inputfile = "{}/{}{}".format(inputpath, disprot_id, method['extension'])

                                if os.path.isfile(inputfile):
                                    if os.stat(inputfile).st_size != 0:
                                        # seq, state, scores = method['func'](inputfile, method)
                                        try:
                                            seq, state, scores = method['func'](inputfile, method)
                                        except Exception as e:
                                            is_good = False
                                            ferr.write("ERROR NOT PARSED {} {}\n{}\n".format(disprot_id, inputfile, e))
                                    else:
                                        is_good = False
                                        ferr.write("ERROR EMPTY FILE {} {}\n".format(disprot_id, inputfile))
                                else:
                                    is_good = False
                                    ferr.write("ERROR MISSING FILE {} {}\n".format(disprot_id, inputfile))

                            if is_good:

                                # Build negatives when the output is missing (useful for region based predictors like foldunfold)
                                if method.get('make_complement'):
                                    state = '0' * len(ref_seq)

                                # print scores
                                # TODO check this test is enough to find errors (e.g. VSL DP01530)
                                if not scores and not state:
                                    ferr.write("ERROR MISSING CORRECT OUTPUT {}\n".format(disprot_id))
                                    is_good = False
                                if state and len(state) != len(ref_seq):
                                    ferr.write(
                                        "ERROR STATE LEN {} {} ({})\n{}\n{}\n{}\n".format(disprot_id, len(state),
                                                                                          len(ref_seq), ref_seq,
                                                                                          state, seq))
                                    is_good = False
                                if scores and len(scores) != len(ref_seq):
                                    ferr.write("ERROR SCORE LEN {} {} ({})\n".format(disprot_id, len(scores),
                                                                                     len(ref_seq)))
                                    is_good = False
                                if is_good:
                                    # Write the output
                                    fout.write(">{}\n{}\n".format(disprot_id, '\n'.join(["{}\t{}\t{}\t{}".format(
                                        i + 1, ref_seq[i], scores[i] if scores else '', state[i] if state else '')
                                                                                         for i in
                                                                                         range(len(ref_seq))])))
                    # Check the format (consume the generator)
                    for _ in parse_prediction_file(outfile):
                        pass

                else:
                    ferr.write("WARNING PARSER NOT IMPLEMENTED\n")
        else:
            print("WARNING NOT WRITING {}".format(outfile))

    return


if __name__ == "__main__":

    # Mount caid folder somewere
    # sudo mount 172.21.2.101:/volume1/Projects projects
    # From IDRA
    # sudo sshfs dampiove@protein.bio.unipd.it:/projects /mnt/projects/ -o IdentityFile='/home/damiano/.ssh/id_rsa',allow_other

    # Generate some data for the method list below
    # ls | while read line; do if[[-d $line / results]]; then ls $line / results | while read d; do echo "#" $line $d; ls $line / results / $d / DP00003 *; done; fi; done
    method_list = [
        {'id': 'D001', 'group': 'callebaut', 'name': 'PyHCA', 'result_dir': 'pyhca', 'extension': '.out', 'func': parse_vertical_format, 'header_lines': 1, 'residue_pos': 1, 'score_pos': 2, 'status_th': 0.5},
        {'id': 'D002', 'group': 'cheng', 'name': 'Predisorder', 'result_dir': 'predisorder', 'extension': '.out', 'func': parse_predisorder},
        {'id': 'B001', 'group': 'dosztanyi', 'name': 'ANCHOR', 'result_dir': 'anchor', 'extension': '.fasta.out', 'func': parse_anchor},
        {'id': 'D003', 'group': 'dosztanyi', 'name': 'IUPred2A-long', 'result_dir': 'iupred2al', 'extension': '.fasta.out', 'func': parse_vertical_format, 'header_lines': 6, 'residue_pos': 1, 'score_pos': 2, 'status_th': 0.5},
        {'id': 'D004', 'group': 'dosztanyi', 'name': 'IUPred2A-short', 'result_dir': 'iupred2as', 'extension': '.fasta.out', 'func': parse_vertical_format, 'header_lines': 6, 'residue_pos': 1, 'score_pos': 2, 'status_th': 0.5},
        {'id': 'B002', 'group': 'dosztanyi', 'name': 'IUPred2A-long', 'result_dir': 'iupred2al', 'extension': '.fasta.out', 'label': '_binding', 'func': parse_vertical_format, 'header_lines': 6, 'residue_pos': 1, 'score_pos': 3, 'status_th': 0.5},
        {'id': 'B003', 'group': 'dosztanyi', 'name': 'IUPred2A-short', 'result_dir': 'iupred2as', 'extension': '.fasta.out', 'label': '_binding', 'func': parse_vertical_format, 'header_lines': 6, 'residue_pos': 1, 'score_pos': 3, 'status_th': 0.5},
        {'id': 'D005', 'group': 'dosztanyi', 'name': 'IUPred-long', 'result_dir': 'iupredl', 'extension': '.flat.out', 'func': parse_vertical_format, 'header_lines': 0, 'residue_pos': 1, 'score_pos': 2, 'status_th': 0.5},
        {'id': 'D006', 'group': 'dosztanyi', 'name': 'IUPred-short', 'result_dir': 'iupreds', 'extension': '.flat.out', 'func': parse_vertical_format, 'header_lines': 2, 'residue_pos': 1, 'score_pos': 2, 'status_th': 0.5},
        {'id': 'D007', 'group': 'galzitskaya', 'name': 'FoldUnfold', 'result_dir': 'foldunfold', 'extension': '.out', 'func': parse_foldunfold, 'make_complement': True},
        {'id': 'D008', 'group': 'galzitskaya', 'name': 'IsUnstruct', 'result_dir': 'isunstruct', 'extension': '.iul', 'func': parse_isunstruct},
        {'id': 'D009', 'group': 'gibson', 'name': 'GlobPlot', 'result_dir': 'globplot', 'extension': '.flat.out', 'func': parse_globplot},
        {'id': 'B004', 'group': 'gsponer', 'name': 'MoRFchibi-light', 'result_dir': 'mcl', 'extension': '.out', 'func': parse_mcl, 'header_lines': 14, 'residue_pos': 1, 'score_pos': 2, 'status_th': 0.725},
        {'id': 'B005', 'group': 'gsponer', 'name': 'MoRFchibi-web', 'result_dir': 'mcw', 'extension': '.out', 'func': parse_mcw, 'header_lines': 17, 'residue_pos': 1, 'score_pos': 2, 'status_th': 0.725, 'use_conservation': True},
        {'id': 'D010', 'group': 'hoque', 'name': 'DisPredict-2', 'result_dir': 'dispred2', 'extension': '.drp', 'func': parse_vertical_format, 'header_lines': 6, 'residue_pos': 0, 'score_pos': 2, 'status_pos': 1, 'status_char': 'D', 'use_conservation': True},
        {'id': 'D011', 'group': 'jones', 'name': 'DISOPRED-3.1', 'result_dir': 'disopred3', 'extension': '.diso', 'label': '_disorder', 'func': parse_vertical_format, 'header_lines': 3, 'residue_pos': 1, 'score_pos': 3, 'status_pos': 2, 'status_char': '*', 'use_conservation': True},
        {'id': 'B013', 'group': 'jones', 'name': 'DISOPRED-3.1', 'result_dir': 'disopred3', 'extension': '.pbdat', 'label': '_binding', 'func': parse_vertical_format, 'header_lines': 5, 'residue_pos': 1, 'score_pos': 3, 'status_pos': 2, 'status_char': '^', 'use_conservation': True},
        {'id': 'D032', 'group': 'kurgan', 'name': 'DFLpred', 'result_dir': 'dflpred', 'extension': '.out', 'func': parse_dflpred},
        {'id': 'B007', 'group': 'kurgan', 'name': 'DisoRDPbind', 'result_dir': 'disordpbind', 'extension': '', 'label': '_all', 'func': parse_disordpbind_all},
        {'id': 'B008', 'group': 'kurgan', 'name': 'DisoRDPbind-RNA', 'result_dir': 'disordpbind', 'extension': '', 'label': '_rna', 'func': parse_disordpbind_rna},
        {'id': 'B009', 'group': 'kurgan', 'name': 'DisoRDPbind-DNA', 'result_dir': 'disordpbind', 'extension': '', 'label': '_dna', 'func': parse_disordpbind_dna},
        {'id': 'B010', 'group': 'kurgan', 'name': 'DisoRDPbind-protein', 'result_dir': 'disordpbind', 'extension': '', 'label': '_protein', 'func': parse_disordpbind_prot},
        {'id': 'D013', 'group': 'kurgan', 'name': 'fIDPln', 'result_dir': 'fidp', 'extension': '.log.pred', 'label': '_log', 'func': parse_vertical_format, 'header_lines': 0, 'residue_pos': 1, 'score_pos': 2, 'status_pos': 3, 'status_char': '1'},
        {'id': 'D014', 'group': 'kurgan', 'name': 'fIDPnn', 'result_dir': 'fidp', 'extension': '.nn.pred', 'label': '_nn', 'func': parse_vertical_format, 'header_lines': 0, 'residue_pos': 1, 'score_pos': 2, 'status_pos': 3, 'status_char': '1'},
        {'id': 'B011', 'group': 'kurgan', 'name': 'fMoRFpred', 'result_dir': 'fmorfpred', 'extension': '.out', 'func': parse_fmorfpred},
        {'id': 'D015', 'group': 'obradovic', 'name': 'VSL2B', 'result_dir': 'vsl2', 'extension': '.flat.out', 'func': parse_vsl2},
        {'id': 'D016', 'group': 'russel', 'name': 'DisEMBL-HL', 'result_dir': 'disembl', 'extension': '.flat.out', 'label': '_HL', 'func': parse_disHL},
        {'id': 'D017', 'group': 'russel', 'name': 'DisEMBL-465', 'result_dir': 'disembl', 'extension': '.flat.out', 'label': '_465', 'func': parse_dis465},
        {'id': 'B012', 'group': 'sharma', 'name': 'OPAL', 'result_dir': 'opal', 'extension': '.txt', 'func': parse_vertical_format, 'header_lines': 1, 'residue_pos': 1, 'score_pos': 2, 'use_conservation': True},
        {'id': 'D018', 'group': 'tosatto', 'name': 'ESpritz-D', 'result_dir': 'espritzd', 'extension': '.disbin.out', 'func': parse_vertical_format, 'header_lines': 8, 'score_pos': 1, 'status_pos': 0, 'status_char': 'D'},
        {'id': 'D019', 'group': 'tosatto', 'name': 'ESpritz-N', 'result_dir': 'espritzn', 'extension': '.disbin.out', 'func': parse_vertical_format, 'header_lines': 8, 'score_pos': 1, 'status_pos': 0, 'status_char': 'D'},
        {'id': 'D020', 'group': 'tosatto', 'name': 'ESpritz-X', 'result_dir': 'espritzx', 'extension': '.disbin.out', 'func': parse_vertical_format, 'header_lines': 8, 'score_pos': 1, 'status_pos': 0, 'status_char': 'D'},
        {'id': 'D021', 'group': 'tosatto', 'name': 'MobiDB-lite', 'result_dir': 'mobidblite', 'extension': '.fasta.out', 'func': parse_mobidblite},
        {'id': 'D022', 'group': 'vendruscolo', 'name': 'S2D-1', 'result_dir': 's2d1', 'extension': '.out', 'func': parse_vertical_format, 'header_lines': 2, 'residue_pos': 1, 'score_pos': 4, 'status_pos': 5, 'status_char': 'C', 'use_conservation': True},
        {'id': 'D023', 'group': 'vendruscolo', 'name': 'S2D-2', 'result_dir': 's2d2', 'extension': '.out', 'func': parse_vertical_format, 'header_lines': 2, 'residue_pos': 1, 'score_pos': 4, 'status_pos': 6, 'status_char': 'C', 'use_conservation': True},
        {'id': 'D024', 'group': 'vranken', 'name': 'DisoMine', 'result_dir': 'disomine', 'extension': '.out', 'func': parse_vertical_format, 'header_lines': 0, 'residue_pos': 2, 'score_pos': 3, 'status_pos': 4, 'status_char': '1'},
        {'id': 'D025', 'group': 'wallner', 'name': 'RawMSA', 'result_dir': 'rawmsa', 'extension': '.pred', 'func': parse_vertical_format, 'header_lines': 0, 'status_pos': 1, 'status_char': '1', 'use_conservation': True},
        {'id': 'D026', 'group': 'xu', 'name': 'AUCpreD', 'result_dir': 'aucpred', 'extension': '.diso', 'func': parse_vertical_format, 'header_lines': 3, 'residue_pos': 1, 'score_pos': 3, 'status_pos': 2, 'status_char': '*', 'use_conservation': True},
        {'id': 'D027', 'group': 'xu', 'name': 'AUCpreD-np', 'result_dir': 'aucpred_no_profile', 'extension': '.diso', 'func': parse_vertical_format, 'header_lines': 3, 'residue_pos': 1, 'score_pos': 3, 'status_pos': 2, 'status_char': '*'},
        {'id': 'D028', 'group': 'zhou', 'name': 'SPOT-Disorder1', 'result_dir': 'spot-disorder1', 'extension': '.spotd', 'func': parse_vertical_format, 'header_lines': 1, 'residue_pos': 1, 'score_pos': 2, 'status_pos': 3, 'status_char': 'D', 'use_conservation': True},
        {'id': 'D029', 'group': 'zhou', 'name': 'SPOT-Disorder2', 'result_dir': 'spot-disorder2', 'extension': '.spotd2', 'func': parse_vertical_format, 'header_lines': 2, 'residue_pos': 1, 'score_pos': 2, 'status_pos': 3, 'status_char': 'D', 'use_conservation': True},
        {'id': 'D030', 'group': 'zhou', 'name': 'SPOT-Disorder-Single', 'result_dir': 'spot-disorder-single', 'extension': '.spotds', 'func': parse_vertical_format, 'header_lines': 2, 'residue_pos': 1, 'score_pos': 2, 'status_pos': 3, 'status_char': 'D'},
        {'id': 'D031', 'group': 'esnouf', 'name': 'JRONN', 'result_dir': 'jronn', 'extension': '.out', 'func': parse_jronn, 'status_th': 0.5},
        {'id': 'D033', 'group': 'vranken', 'name': 'DynaMine', 'result_dir': 'dynamine', 'extension': '', 'func': None}
    ]



    # Parse reference sequences
    # ref_seqs = parse_reference_sequences("{}/CAID/2018_09/disprot/disprot8_all.fasta".format(result_dir))
    reference_sequences = parse_reference_sequences("/home/damiano/Projects/caid/data/disprot8_all.fasta")

    overwrite_output = False
    # write_prsed_results(method_list, reference_sequences, overwrite_output)

    # TODO calculate DynaMine
    # Print the {id: name} dictionary
    for method in method_list:
        print("{}\t{}\t{}".format(method['id'], method['name'], method.get('use_conservation', False)))


