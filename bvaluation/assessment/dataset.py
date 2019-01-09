import logging
import numpy as np
from itertools import groupby, chain


class ImmutableDict(dict):
    def __init__(self):
        dict.__init__(self)

    def __setitem__(self, key, val):
        logging.error('Trying to set %s[%s] = %s. SET operation has been disabled',
                      self.__class__.__name__, key, val)
        raise PermissionError('SET operation on Reference not allowed')

    def __delitem__(self, key):
        logging.error('Trying to delete %s field from %s. DELETE operation has been disabled',
                      key, self.__class__.__name__)
        raise PermissionError('DELETE operation on Reference not allowed')


class ReferencePool(ImmutableDict):
    def __init__(self, fname, strict_negatives):
        dict.__init__(self)
        self.fname = fname
        self.strict_negatives = strict_negatives
        self.pattern = {'D': 1,
                        'S': 0,
                        'B': 1,
                        '1': 1,
                        '0': 0,
                        '-': np.nan if self.strict_negatives is True else 0}
        self._build_from_file()

    def _build_from_file(self):
        with open(self.fname) as f:
            '''parse a specific format to build the Reference instance'''
            faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
            for header in faiter:
                header = next(header)
                if header[0] != '#':
                    disprot_acc, *desc = header[1:].strip().split()
                    seq, states = map(str.strip, next(faiter))

                    nastates = np.fromiter(map(self.pattern.get, states), dtype=float)
                    ref_entry = ReferenceEntry(desc, seq, nastates)
                    dict.__setitem__(self, disprot_acc, ref_entry)

    def make_pure_reference(self):
        ref = Reference()
        for acc, data in self.items():
            ref.states.append(data['states'])
            ref.add_accession(acc)
        ref.set_merged_states()

        return ref


class ReferenceEntry(ImmutableDict):
    def __init__(self, uniprot_acc, seq, nastates):
        dict.__init__(self)
        self._build_entry(uniprot_acc, seq, nastates)

    def _build_entry(self, u, s, a):
        dict.__setitem__(self, 'uniprot_acc', u)
        dict.__setitem__(self, 'seq', s)
        dict.__setitem__(self, 'states', a)


class PredictionEntry(ImmutableDict):
    def __init__(self, p, r, sc, st):
        dict.__init__(self)
        self.validity_keys = list()

        self._build_entry(p, r, sc, st)
        self._check_consistency()
        self.is_consistent = None

    def _build_entry(self, p, r, sc, st):
        if p:
            dict.__setitem__(self, 'positions', np.array(p, dtype=np.int))
            self.validity_keys.append('positions')
        if r:
            dict.__setitem__(self, 'seq', ''.join(r))
            self.validity_keys.append('seq')
        if sc and sc[0]:
            dict.__setitem__(self, 'scores', np.array(sc, dtype=np.float))
            self.validity_keys.append('scores')
        if st and st[0]:
            dict.__setitem__(self, 'states', np.array(st, dtype=np.float))
            self.validity_keys.append('states')

    def _check_consistency(self):
        lengths = set()
        for key in self.validity_keys:
            if key != 'seq':
                lengths.add(np.shape(self[key][0]))
            else:
                lengths.add(len(self[key]))

        if len(lengths) == 1:
            self.is_consistent = True
        else:
            self.is_consistent = False


class Reference(object):
    def __init__(self):
        self.accessions = list()
        self.accessions_set = set()
        self.states = list()
        self.mstates = None

    def add_accession(self, acc):
        self.accessions.append(acc)
        self.accessions_set.add(acc)

    def set_merged_states(self):
        self.mstates = np.fromiter(chain(*self.states), dtype=np.float)

    def __str__(self):
        accs = 'accessions    {:>8} [{} {} ...]'.format(len(self.accessions), *self.accessions)
        stts = 'states        {:>8} [{} {} ...]   [{} {} ...] ...'.format(len(self.states),
                                                                          *self.states[0][:2],
                                                                          *self.states[1][:2])
        mgst = 'merged states {:>8} [{} {} {} {} {} {} ...]'.format(
            len(self.mstates),
            *self.mstates) if self.mstates is not None else 'merged states {:>8}'.format('nan')

        return '{}\n{}\n{}'.format(accs, stts, mgst)


class Prediction(object):
    def __init__(self):
        self.accessions = list()
        self.accessions_set = set()
        self.states = list()
        self.mstates = None
        self.scores = list()
        self.mscores = None
        self.coverage = None

    def add_accession(self, acc):
        self.accessions.append(acc)
        self.accessions_set.add(acc)

    def set_coverage(self, ref_accs):
        self.coverage = len(ref_accs & self.accessions_set) / len(ref_accs)

    def set_merged_states(self):
        self.mstates = np.fromiter(chain(*self.states), dtype=np.float)
        self.mscores = np.fromiter(chain(*self.scores), dtype=np.float)

    def apply_cutoff(self, cutoff):
        self.states = list(map(lambda a: (a >= cutoff).astype(int), self.scores))
        self.mstates = np.fromiter(chain(*self.states), dtype=np.float)
