import logging
import copy
import numpy as np
from itertools import groupby, chain


class ImmutableDict(dict):
    """
    Sublcass of dict that cannot be modified after instantiation
    """
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
    """
    Reference dict containing annotation for all possible targets.

    Items can only be set at instantiation. Further `keys: values` pairs cannot be set
    """
    def __init__(self, fname: str, undefined_replace_value: int):
        dict.__init__(self)
        self.fname = fname
        self.undef_replace_value = undefined_replace_value
        self.pattern = {'0': 0,
                        '1': 1,
                        '-': np.nan if undefined_replace_value is None else undefined_replace_value}

        self._build_from_file()

    def _build_from_file(self):
        """
        Set self keys parsing a reference file
        """
        with open(self.fname) as f:
            '''parse a specific format to build the Reference instance'''
            faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
            for header in faiter:
                header = next(header)
                if header[0] != '#':
                    acc, *desc = header[1:].strip().split()
                    seq, states = map(str.strip, next(faiter))

                    nastates = np.fromiter(map(self.pattern.get, states), dtype=float)
                    ref_entry = ReferenceEntry(desc, nastates)
                    dict.__setitem__(self, acc, ref_entry)
        logging.debug('ref pool: %s', self)

    def make_pure_reference(self):
        """
        Return a reference object containing all targets in reference pool
        :return: reference instance
        """
        ref = Reference()
        for acc, data in self.items():
            states = data['states']

            if self.undef_replace_value is None:
                states = states[~np.isnan(states)]

            ref.states.append(states)
            ref.add_accession(acc)
        ref.set_merged_states()

        return ref


class ReferenceEntry(ImmutableDict):
    """
    Entry of a reference pool

    Items can only be set at instantiation. Further `keys: values` pairs cannot be set

    :param accession: reference entry accession
    :param nastates: states
    """
    def __init__(self, accession: str, nastates: np.array):
        dict.__init__(self)
        self._build_entry(accession, nastates)

    def _build_entry(self, acc: str, s: np.array):
        dict.__setitem__(self, 'acc', acc)
        dict.__setitem__(self, 'states', s)
        

class PredictionEntry(ImmutableDict):
    """
    Entry for a prediction set

    Items can only be set at instantiation. Further `keys: values` pairs cannot be set

    :param p: positions
    :param sc: scores
    :param st: states
    """
    def __init__(self, p, sc, st):

        dict.__init__(self)
        self.validity_keys = list()

        self._build_entry(p, sc, st)
        self._check_consistency()
        self.is_consistent = None

    def _build_entry(self, p, sc, st):
        if p:
            dict.__setitem__(self, 'positions', np.array(p, dtype=np.int))
            self.validity_keys.append('positions')
        if sc and sc[0]:
            dict.__setitem__(self, 'scores', np.array(sc, dtype=np.float))
            self.validity_keys.append('scores')
        if st and st[0]:
            dict.__setitem__(self, 'states', np.array(st, dtype=np.float))
            self.validity_keys.append('states')

    def _check_consistency(self):
        lengths = set()
        for key in self.validity_keys:
            lengths.add(len(self[key]))

        if len(lengths) == 1:
            self.is_consistent = True
        else:
            self.is_consistent = False


class Reference(object):
    """
    A reference set
    """
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
        logging.debug('ref mstates: %i %s', len(self.mstates), self.mstates)

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
    """
    A prediction set
    """
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
        self.coverage = (len(ref_accs & self.accessions_set), len(ref_accs))

    def set_merged_states(self):
        self.mstates = np.fromiter(chain(*self.states), dtype=np.float)
        self.mscores = np.fromiter(chain(*self.scores), dtype=np.float)
        logging.debug('pred mstates: %i %s', len(self.mstates), self.mstates)
        logging.debug('pred mscores: %i %s', len(self.mscores), self.mscores)

    def apply_cutoff(self, cutoff):
        self.states = [np.greater_equal(s, cutoff).astype(int) for s in self.scores]
        self.mstates = np.fromiter(chain(*self.states), dtype=np.float)

    def zip(self):
        return zip(self.accessions, self.states)


    def copy(self):
        return copy.deepcopy(self)