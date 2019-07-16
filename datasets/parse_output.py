import os
from itertools import groupby
import logging.config
import json
import traceback


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


def _regions2state_mdbl(regions, state_len):
    state = ['0' for i in range(state_len)]
    for start, end, status in regions:
        for i in range(start, end + 1):
            # print i
            state[i-1] = '1'
    return ''.join(state)


# Parse reference sequences
def parse_reference_sequences(fastafile, filter_targets=None):

    targets = []
    with open(filter_targets) as f:
        for line in f:
            targets.append(line.strip())
    targets = frozenset(targets)

    r_seqs = {}
    with open(fastafile) as f:
        for line in f:
            if line[0] == '>':
                name = line.rstrip().split()[0][1:]
                r_seqs.setdefault(name, '')
            else:
                r_seqs[name] += line.strip()

    new_r_seqs = {}
    if filter_targets:
        for name in r_seqs:
            if name in targets:
                new_r_seqs[name] = r_seqs[name]
    else:
        new_r_seqs = r_seqs

    return new_r_seqs


def strip_split(string, expected=None):
    if string:
        s = string.strip('\n').rsplit('\t')
        if expected and len(s) != expected:
            # logging.error('expecting %i elements from split, got %i', expected, len(s))
            raise ValueError('expecting %i elements from split, got %i', expected, len(s))
        return s


###########################################################


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

    state = ['1' if l.isupper() else '0' for l in seq]  # Disorder
    return seq.upper(), state, state_rna, scores_rna, state_dna, scores_dna, state_prot, scores_prot


def parse_fmorfpred(inputfile, kwargs):
    # The second state row is just majority vote between iupred and espritz, so it is discarder
    with open(inputfile) as f:
        next(f)
        seq = next(f).strip()
        scores = next(f).strip().split(',')
        state = next(f).strip()
    return seq.upper(), state, scores


def parse_mobidblite(inputfile, kwargs):
    with open(inputfile) as f:
        obj = json.loads(f.read())
    pred = list(filter(lambda k: k['method'] == kwargs['method_name'], obj['predictions']))
    if pred:
        pred = pred[0]
        return None, _regions2state_mdbl(pred['regions'], len(pred['scores'])), pred['scores']
    else:
        return None, None, None


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


###########################################################


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


def write_parsed_results(methods, ref_seqs, **kwargs):

    # Create local folder
    if not os.path.exists(kwargs.get('parsed_dir')):
        os.makedirs(kwargs.get('parsed_dir'))

    for method in methods:

        # Input file
        inputpath = "{}/{}".format(kwargs.get('raw_results_dir'), method['dir'])
        logging.info("Parsing {} {}".format(inputpath, method['name']))

        if os.path.isdir(inputpath):

            # Output file
            outfile = "{}/{}_{}.out".format(kwargs.get('parsed_dir'), method['id'], method['name'])  # ex. D012_MobiDB-lite.out

            if not os.path.isfile(outfile) or kwargs.get('overwrite_output'):
                with open(outfile, 'w') as fout:

                    for disprot_id in ref_seqs:

                        ref_seq = ref_seqs[disprot_id]
                        method['seq'] = ref_seq

                        is_good = True
                        seq = None
                        state = None
                        scores = None

                        if 'extension' in method:
                            inputfile = "{}/{}{}".format(inputpath, disprot_id, method['extension'])
                        else:
                            inputfile = "{}/{}.out".format(inputpath, disprot_id)

                        if os.path.isfile(inputfile):
                            if os.stat(inputfile).st_size != 0:
                                # seq, state, scores = method['func'](inputfile, method)
                                try:
                                    seq, state, scores = method['func'](inputfile, method)
                                except Exception as e:
                                    is_good = False
                                    logging.error("Failed to parse: {} {} {}\n{}".format(disprot_id, method['name'], inputfile, traceback.format_exc()))
                            else:
                                is_good = False
                                logging.error("Empty file: {} {} {}".format(disprot_id, method['name'], inputfile))
                        else:
                            is_good = False
                            logging.error("Missing file: {} {} {}".format(disprot_id, method['name'], inputfile))

                        if is_good:

                            # print scores
                            # TODO check this test is enough to find errors (e.g. VSL DP01530)
                            if not scores and not state:
                                logging.error("Wrong output format {} {}".format(disprot_id, method['name']))
                                is_good = False
                            if state and len(state) != len(ref_seq):
                                logging.error("Wrong state length {} {} {} ({})\n{}\n{}\n{}".format(disprot_id, method['name'], len(state),
                                                                                      len(ref_seq), ref_seq,
                                                                                      state, seq))
                                is_good = False
                            if scores and len(scores) != len(ref_seq):
                                logging.error("Wrong score length {} {} {} ({})".format(disprot_id, method['name'], len(scores),
                                                                                 len(ref_seq)))
                                is_good = False

                            # write the output
                            if is_good:
                                fout.write(">{}\n{}\n".format(disprot_id, '\n'.join(["{}\t{}\t{}\t{}".format(
                                    i + 1, ref_seq[i], scores[i] if scores else '', state[i] if state else '')
                                                                                     for i in
                                                                                     range(len(ref_seq))])))
                # Check the format (consume the generator)
                for _ in parse_prediction_file(outfile):
                    pass

            else:
                logging.warning("Output file exists. Not writing {}".format(outfile))
        else:
            logging.error("Missing input dir {}".format(inputpath))

        if os.stat(outfile).st_size == 0:
            logging.warning("Empty output file, deleting {}".format(outfile))
            os.remove(outfile)
    return


if __name__ == "__main__":


    method_list = [
        {'name': 'PyHCA', 'func': parse_vertical_format, 'header_lines': 1, 'residue_pos': 1, 'score_pos': 2, 'status_th': 0.5},
        {'name': 'IUPred2A-long', 'func': parse_vertical_format, 'header_lines': 7, 'residue_pos': 1, 'score_pos': 2, 'status_th': 0.5},
        {'name': 'IUPred2A-short', 'func': parse_vertical_format, 'header_lines': 7, 'residue_pos': 1, 'score_pos': 2, 'status_th': 0.5},
        {'name': 'ANCHOR-2', 'func': parse_vertical_format, 'header_lines': 7, 'residue_pos': 1, 'score_pos': 3, 'status_th': 0.5},
        {'name': 'DisPredict-2', 'func': parse_vertical_format, 'header_lines': 6, 'residue_pos': 0, 'score_pos': 2, 'status_pos': 1, 'status_char': 'D'},
        {'name': 'OPAL', 'func': parse_vertical_format, 'header_lines': 1, 'residue_pos': 1, 'score_pos': 2},
        {'name': 'S2D-1', 'func': parse_vertical_format, 'header_lines': 2, 'residue_pos': 1, 'score_pos': 4, 'status_pos': 5, 'status_char': 'C'},
        {'name': 'S2D-2', 'func': parse_vertical_format, 'header_lines': 2, 'residue_pos': 1, 'score_pos': 4, 'status_pos': 6, 'status_char': 'C'},
        {'name': 'DisoMine', 'func': parse_vertical_format, 'header_lines': 0, 'residue_pos': 2, 'score_pos': 3, 'status_pos': 4, 'status_char': '1'},
        {'name': 'RawMSA', 'func': parse_vertical_format, 'header_lines': 0, 'score_pos': 1, 'status_th': 0.5},
        {'name': 'AUCpreD', 'func': parse_vertical_format, 'header_lines': 3, 'residue_pos': 1, 'score_pos': 3, 'status_pos': 2, 'status_char': '*'},
        {'name': 'AUCpreD-np', 'func': parse_vertical_format, 'header_lines': 3, 'residue_pos': 1, 'score_pos': 3, 'status_pos': 2, 'status_char': '*'},
        {'name': 'SPOT-Disorder1', 'func': parse_vertical_format, 'header_lines': 1, 'residue_pos': 1, 'score_pos': 2, 'status_pos': 3, 'status_char': 'D'},
        {'name': 'SPOT-Disorder2', 'func': parse_vertical_format, 'header_lines': 2, 'residue_pos': 1, 'score_pos': 2, 'status_pos': 3, 'status_char': 'D'},
        {'name': 'SPOT-Disorder-Single', 'func': parse_vertical_format, 'header_lines': 2, 'residue_pos': 1, 'score_pos': 2, 'status_pos': 3, 'status_char': 'D'},

        {'name': 'DISOPRED-3.1', 'func': parse_vertical_format, 'header_lines': 3, 'residue_pos': 1, 'score_pos': 3, 'status_pos': 2, 'status_char': '*'},
        {'name': 'DISOPRED-3.1-binding', 'func': parse_vertical_format, 'header_lines': 5, 'residue_pos': 1, 'score_pos': 3, 'status_pos': 2, 'status_char': '^'},
        {'name': 'fIDPln', 'extension': '.log_out', 'func': parse_vertical_format, 'header_lines': 0, 'residue_pos': 1, 'score_pos': 2, 'status_pos': 3, 'status_char': '1'},
        {'name': 'fIDPnn', 'extension': '.nn_out', 'func': parse_vertical_format, 'header_lines': 0, 'residue_pos': 1, 'score_pos': 2, 'status_pos': 3, 'status_char': '1'},

        {'name': 'DisoRDPbind', 'func': parse_disordpbind_all},
        {'name': 'DisoRDPbind-RNA', 'func': parse_disordpbind_rna},
        {'name': 'DisoRDPbind-DNA', 'func': parse_disordpbind_dna},
        {'name': 'DisoRDPbind-protein', 'func': parse_disordpbind_prot},

        {'name': 'Predisorder', 'func': parse_predisorder},
        {'name': 'FoldUnfold', 'func': parse_foldunfold},
        {'name': 'IsUnstruct', 'func': parse_isunstruct},
        {'name': 'MoRFchibi-light', 'func': parse_mcl, 'header_lines': 14, 'residue_pos': 1, 'score_pos': 2, 'status_th': 0.725},
        {'name': 'MoRFchibi-web', 'func': parse_mcw, 'header_lines': 17, 'residue_pos': 1, 'score_pos': 2, 'status_th': 0.725},
        {'name': 'DFLpred', 'func': parse_dflpred},
        {'name': 'fMoRFpred', 'func': parse_fmorfpred},

        {'name': 'MobiDB-lite', 'func': parse_mobidblite, 'method_name': 'mobidb_lite'},
        {'name': 'ANCHOR', 'func': parse_mobidblite, 'method_name': 'anchor'},
        {'name': 'IUPred-long', 'func': parse_mobidblite, 'method_name': 'iupl'},
        {'name': 'IUPred-short', 'func': parse_mobidblite, 'method_name': 'iups'},
        {'name': 'GlobPlot', 'func': parse_mobidblite, 'method_name': 'glo', 'rescale': True},
        {'name': 'VSL2B', 'func': parse_mobidblite, 'method_name': 'vsl'},
        {'name': 'DisEMBL-HL', 'func': parse_mobidblite, 'method_name': 'disHL'},
        {'name': 'DisEMBL-465', 'func': parse_mobidblite, 'method_name': 'dis465'},
        {'name': 'ESpritz-D', 'func': parse_mobidblite, 'method_name': 'espD'},
        {'name': 'ESpritz-N', 'func': parse_mobidblite, 'method_name': 'espN'},
        {'name': 'ESpritz-X', 'func': parse_mobidblite, 'method_name': 'espX'},
        {'name': 'DynaMine', 'func': parse_mobidblite, 'method_name': 'dynamine'},
        {'name': 'JRONN', 'func': parse_mobidblite, 'method_name': 'jronn'}
    ]

    out_log_file = "parse_output.log"

    logging.basicConfig(filename=out_log_file,
                        level=logging.getLevelName("INFO"),
                        format='%(asctime)s - %(process)d - %(levelname)s - %(message)s')

    logging.info("CAID parsing started")

    module_dir = os.path.abspath(os.path.dirname(__file__))
    logging.info("Module running from: {}".format(module_dir))

    # Mount caid folder somewhere
    # sudo apt install nfs-common
    # (obsolete) sudo sshfs dampiove@protein.bio.unipd.it:/projects /mnt/projects/ -o IdentityFile='/home/damiano/.ssh/id_rsa',allow_other
    # (obsolete) sudo mount 172.21.2.103:/volume1/Pleiadi/projects/CAID/2018_09 /home/damiano/Projects/caid_data/remote/
    # sudo mount 172.21.2.103:/volume1/Pleiadi/projects/CAID/results /home/damiano/Projects/caid/data/results/

    with open("../data/method_names.json") as f:
        objs = json.load(f)

    for method in method_list:
        for ele in objs:
            if ele['name'] == method['name']:
                method.update(ele)
                break


    # Parse reference sequences
    reference_sequences = parse_reference_sequences("../data/disprot8_all.fasta", filter_targets="../data/new_id_list.txt")

    write_parsed_results(method_list,
                         reference_sequences,
                         raw_results_dir="/home/damiano/Projects/caid/data/results",
                         parsed_dir="/home/damiano/Projects/caid/data/results_parsed",
                         overwrite_output=False)

    # Print the dictionary {id: name}
    # for method in method_list:
    #     print("{}\t{}\t{}".format(method['id'], method['name'], method.get('use_conservation', False)))


    # TODO chek opposite predictions (see rocs)

    # TODO ifdp has two extentions in the same folder for two different predictors

    # TODO check foldunfold empty file means fully structured or crash

