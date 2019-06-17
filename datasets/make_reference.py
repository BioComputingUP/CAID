#!/usr/bin/env python3
import sys
import os
import numpy as np
import json
import copy
import logging
from datetime import datetime, date
from collections import Counter
from tempfile import NamedTemporaryFile
from itertools import groupby

from launchers import interproscan
from launchers import fess


# TODO: make reference from predicted secondary structure
# TODO: NORS by Rost
# TODO: filter polyproteins


def parse_entries(disprot_entries_file):
    # TODO Fix the date in the collection

    regions = []
    proteins = {}

    # Parse the collection
    with open(disprot_entries_file) as f:
        for line in f:
            entry = json.loads(line)

            if entry['record_type'] == 'protein_record':
                # print entry

                # WARNING some entries have this wrong format:
                # 'last_edit_date': {u'$date': u'2016-08-11T18:07:34.900+0200'}
                # Otherwise:
                # 2018-09-08T15:40:08.750Z
                # 2016-08-23T16:37:21.222+0200
                if isinstance(entry.get('last_edit_date'), dict):
                    entry['last_edit_date'] = entry.get('last_edit_date')['$date']

                # Fix date format
                entry['last_edit_date'] = datetime.strptime(entry['last_edit_date'][:19],
                                                            "%Y-%m-%dT%H:%M:%S")
                entry['creation_date'] = datetime.strptime(entry['creation_date'][:19],
                                                           "%Y-%m-%dT%H:%M:%S")

                proteins[entry['disprot_id']] = entry
                proteins[entry['disprot_id']]['regions'] = []

            elif entry['record_type'] == 'region_record':
                # print entry
                regions.append(entry)

    # Add regions to protein objects
    for region in regions:
        proteins[region['disprot_id']]['regions'].append(region)

    return proteins


def filter_bad_entries(proteins):
    """
    Filter out entries without regions and entries with only ambiguous regions
    """
    ids = list(proteins.keys())
    for disprot_id in ids:  # Use keys to delete "proteins" keys while looping
        if proteins[disprot_id].get('obsolete').get('reason_curator') or proteins[disprot_id].get(
                'obsolete').get('reason_system'):
            del proteins[disprot_id]
        else:
            if len(proteins[disprot_id]['regions']) == 0:
                del proteins[disprot_id]
            else:
                # Classify bad/good regions
                short_regions = []
                good_regions = []
                bad_regions = []
                for region in proteins[disprot_id]['regions']:
                    if region['end'] - region['start'] + 1 < 10:
                        short_regions.append(region)
                    else:
                        if region.get('tags') == [None]:
                            good_regions.append(region)
                        else:
                            bad_regions.append(region)

                # Proteins with only ambiguous regions
                if len(good_regions) == 0:
                    del proteins[disprot_id]

                else:
                    proteins[disprot_id]['regions'] = good_regions

    return proteins


def filter_date(proteins, min_date, max_date):
    ids = list(proteins.keys())
    for disprot_id in ids:
        protein = proteins[disprot_id]
        if protein['creation_date'].date() < min_date or protein['creation_date'].date() > max_date:
            del proteins[disprot_id]
    return proteins


def filter_polyproteins(proteins):
    ids = list(proteins.keys())
    for disprot_id in ids:
        if 'polyprotein' in proteins[disprot_id]['protein_name']:
            del proteins[disprot_id]
    return proteins


# TODO check/rewrite this function
def get_identical_regions(proteins, same_region_outfile):
    """
    Test identical regions in multiple proteins (transfer by homology)
    """

    region_ids = {}
    for disprot_id in proteins:
        for region in proteins[disprot_id]['regions']:
            if region.get('tags') == [None]:
                region_id = "{}-{}-{}-{}".format(region['pmid'], region['method']['id'],
                                                 region['start'], region['end'])
                region_ids.setdefault(region_id, set())
                region_ids[region_id].add(proteins[disprot_id]['disprot_id'])
            else:
                print(disprot_id, region.get('tags'))

    disprot_ids = {}
    for region_id in region_ids:
        if len(region_ids[region_id]) > 1:
            for disprot_id in region_ids[region_id]:
                disprot_ids.setdefault(disprot_id, set())
                disprot_ids[disprot_id].add(region_id)

    with open(same_region_outfile, 'w') as fout:
        fout.write("Identical regions (pmid,method,star,end) in different proteins\n"
                   "Proteins {}\nRegions {}\nUniRef90 {}\nUniRef50 {}\nPMID {}\n\n"
                   "Region ID\tDisProt IDs\n{}\n\n"
                   "DisProt ID\tRegions\n{}\n".format(
            len(set([disprot_id for region_id in region_ids if len(region_ids[region_id]) > 1 for
                     disprot_id in
                     region_ids[region_id]])),
            len(set([region_id for region_id in region_ids if len(region_ids[region_id]) > 1])),
            len(set(
                [proteins[disprot_id].get('uniref90') for region_id in region_ids if
                 len(region_ids[region_id]) > 1 for
                 disprot_id in region_ids[region_id]])),
            len(set(
                [proteins[disprot_id].get('uniref50') for region_id in region_ids if
                 len(region_ids[region_id]) > 1 for
                 disprot_id in region_ids[region_id]])),
            len(set([region_id.split('-')[0] for region_id in region_ids if
                     len(region_ids[region_id]) > 1])),
            "\n".join(sorted(
                ["{}\t{}".format(region_id, ",".join(region_ids[region_id])) for region_id in
                 region_ids if
                 len(region_ids[region_id]) > 1])),
            "\n".join(sorted(
                ["{}\t{}".format(disprot_id, ",".join(list(disprot_ids[disprot_id]))) for disprot_id
                 in disprot_ids]))
        ))


def write_fasta(proteins, outfile):
    with open(outfile, 'w') as fout:
        for disprot_id in proteins:
            fout.write('>{} {} {} {} {} {} {} {} {} {}\n{}\n'.format(
                  proteins[disprot_id]['disprot_id'], proteins[disprot_id]['uniprot_accession'],
                  proteins[disprot_id]['uniprot_id'], proteins[disprot_id].get('uniref90', None),
                  proteins[disprot_id].get('uniref50', None),
                  proteins[disprot_id]['creation_date'], proteins[disprot_id]['last_edit_date'],
                  proteins[disprot_id]['curator'], proteins[disprot_id]['organism'], proteins[disprot_id]['protein_name'],
                  proteins[disprot_id]['sequence']))


def write_sentences(proteins, outfile):
    with open(outfile, 'w') as fout:
        fout.write(
            "#disprot\tuniprot\tpmid\tmethod\tregion_start\tregion_end\ttype_or_section\ttext\n")
        for disprot_id in proteins:
            protein = proteins[disprot_id]
            for region in protein['regions']:
                if region.get('generif'):
                    fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(disprot_id, protein[
                        'uniprot_accession'], region['pmid'], region['method']['name'],
                                                                         region['start'],
                                                                         region['end'],
                                                                         region['generif']['loc'],
                                                                         region['generif'][
                                                                             'text'].encode(
                                                                             'utf-8')))


def write_caid_references_pdb(proteins, outdir, label, ref_names):
    data = {}
    for k in ref_names:
        data.setdefault(k, {})
        data[k]['outf'] = open("{}/{}-{}.txt".format(outdir, label, k), 'w')

    for disprot_id in proteins:
        protein = proteins[disprot_id]

        # Region information
        regions_names = []
        methods = []
        for region in protein['regions']:
            methods.append(region['method']['id'])
            func_names = [f['name'] for f in region['functions']]
            regions_names += func_names

        for k in ref_names:
            key = k.split('_')[0]
            data[k]['state'] = ['-' for i in range(len(protein.get('sequence')))]

            # Assign pdb residues
            for ele in protein.get("mobidb", []).get("pdb_missing_residues", []):
                if ele['ann'] == 'S':
                    # Assign structure for all type of references
                    for i in range(ele['start'] - 1, ele['end']):
                        data[k]['state'][i] = '0'
                elif ele['ann'] == 'D':
                    if key in ['pdb-missing']:
                        for i in range(ele['start'] - 1, ele['end']):
                            data[k]['state'][i] = '1'
                elif ele['ann'] == 'C':
                    pass

            # Assign positive labels (overwrite structure)
            for region in protein['regions']:
                method = region['method']['id']
                func_names = [f['name'] for f in region['functions']]

                for i in range(int(region['start']) - 1, int(region['end'])):
                    if key in ['disprot-all']:
                        data[k]['state'][i] = '1'

                    if any('binding' in func_name.lower() for func_name in func_names):
                        if key in ['disprot-binding']:
                            data[k]['state'][i] = '1'

                    if any('rna' in func_name.lower() for func_name in func_names):
                        if key in ['disprot-binding', 'disprot-binding-rna']:
                            data[k]['state'][i] = '1'

                    if any('dna' in func_name.lower() for func_name in func_names):
                        if key in ['disprot-binding', 'disprot-binding-dna']:
                            data[k]['state'][i] = '1'

                    if any('protein-protein' in func_name.lower() for func_name in func_names):
                        if key in ['disprot-binding', 'disprot-binding-prot']:
                            data[k]['state'][i] = '1'

                    if any('linker' in func_name.lower() for func_name in func_names):
                        if key in ['disprot-linker']:
                            data[k]['state'][i] = '1'

                    if method == 'NMR':
                        if key in ['disprot-primary', 'disprot-nmr']:
                            data[k]['state'][i] = '1'

                    if method == 'XRAY':
                        if key in ['disprot-primary', 'disprot-xray']:
                            data[k]['state'][i] = '1'

                    if method not in ['NMR', 'XRAY']:
                        if key in ['disprot-secondary']:
                            data[k]['state'][i] = '1'

            # Write to file
            if '1' in data[k]['state']:
                data[k]['outf'].write(
                    ">{} {} {} {} {}\n{}\n{}\n".format(disprot_id, protein['uniprot_accession'],
                                                       ','.join(methods), ','.join(regions_names),
                                                       ",".join(
                                                           ["{}-{}:{}".format(ele['start'],
                                                                              ele['end'],
                                                                              ele['ann']) for ele in
                                                            protein.get("mobidb", []).get(
                                                                "pdb_missing_residues", [])]),
                                                       protein['sequence'],
                                                       ''.join(data[k]['state'])))


def write_caid_references_simple(proteins, outdir, label, ref_names):
    data = {}
    for k in ref_names:
        data.setdefault(k, {})
        data[k]['outf'] = open("{}/{}-{}.txt".format(outdir, label, k), 'w')

    for disprot_id in proteins:
        protein = proteins[disprot_id]

        # Region information
        regions_names = []
        methods = []
        for region in protein['regions']:
            methods.append(region['method']['id'])
            func_names = [f['name'] for f in region['functions']]
            regions_names += func_names

        for k in ref_names:
            key = k.split('_')[0]
            data[k]['state'] = ['0' for _ in range(len(protein.get('sequence')))]

            # Assign pdb residues
            for ele in protein.get("mobidb", []).get("pdb_missing_residues", []):
                if ele['ann'] == 'D':
                    if key in ['pdb-missing']:
                        for i in range(ele['start'] - 1, ele['end']):
                            data[k]['state'][i] = '1'

            # Assign positive labels (overwrite structure)
            for region in protein['regions']:
                method = region['method']['id']
                func_names = [f['name'] for f in region['functions']]

                for i in range(int(region['start']) - 1, int(region['end'])):
                    if key in ['disprot-all']:
                        data[k]['state'][i] = '1'

                    if any('binding' in func_name.lower() for func_name in func_names):
                        if key in ['disprot-binding']:
                            data[k]['state'][i] = '1'

                    if any('rna' in func_name.lower() for func_name in func_names):
                        if key in ['disprot-binding', 'disprot-binding-rna']:
                            data[k]['state'][i] = '1'

                    if any('dna' in func_name.lower() for func_name in func_names):
                        if key in ['disprot-binding', 'disprot-binding-dna']:
                            data[k]['state'][i] = '1'

                    if any('protein-protein' in func_name.lower() for func_name in func_names):
                        if key in ['disprot-binding', 'disprot-binding-prot']:
                            data[k]['state'][i] = '1'

                    if any('linker' in func_name.lower() for func_name in func_names):
                        if key in ['disprot-linker']:
                            data[k]['state'][i] = '1'

                    if method == 'NMR':
                        if key in ['disprot-primary', 'disprot-nmr']:
                            data[k]['state'][i] = '1'

                    if method == 'XRAY':
                        if key in ['disprot-primary', 'disprot-xray']:
                            data[k]['state'][i] = '1'

                    if method not in ['NMR', 'XRAY']:
                        if key in ['disprot-secondary']:
                            data[k]['state'][i] = '1'

            # print "{:<25}{}".format(k, ''.join(data[k]['state']))

            # Write to file
            if '1' in data[k]['state']:
                data[k]['outf'].write(
                    ">{} {} {} {} {}\n{}\n{}\n".format(disprot_id, protein['uniprot_accession'],
                                                       ','.join(methods), ','.join(regions_names),
                                                       ",".join(
                                                           ["{}-{}:{}".format(ele['start'],
                                                                              ele['end'],
                                                                              ele['ann']) for ele in
                                                            protein.get("mobidb", []).get(
                                                                "pdb_missing_residues", [])]),
                                                       protein['sequence'],
                                                       ''.join(data[k]['state'])))


def caid_reference_stats(inputdir, label, ref_names):
    data = {}
    c = 0
    for k in ref_names:
        data.setdefault(k, [])
        with open("{}/{}_{}.txt".format(inputdir, label, k)) as f:
            for line in f:
                if line[0] != '#':
                    if line[0] == '>':
                        # name = line.strip().split()[0][1:]
                        c = 0
                    elif c == 1:
                        pass
                    elif c == 2:
                        data[k].append(line.strip())
                    c += 1

    for k in ref_names:

        # print data[k][0]
        # print count_regions(data[k][0])
        # break

        counts_str = Counter()
        regions_str = Counter()
        regions_str_long = 0

        counts_len = Counter()
        regions_len = Counter()
        regions_len_long = 0

        for state in data[k]:
            regdist_str, regdist_str_size = count_regions(state)
            regdist_len, regdist_len_size = count_regions(state.replace('-', '0'))

            regions_str_long += sum(
                [1 if ann == '1' and (end - start + 1) >= 20 else 0 for start, end, ann in
                 regdist_str_size])
            regions_len_long += sum(
                [1 if ann == '1' and (end - start + 1) >= 20 else 0 for start, end, ann in
                 regdist_len_size])

            counts_str.update(Counter(state))
            regions_str.update(Counter(regdist_str))
            counts_len.update(Counter(state.replace('-', '0')))
            regions_len.update(Counter(regdist_len))

        tot_str = sum([counts_str[i] for i in counts_str])
        tot_len = sum([counts_len[i] for i in counts_len])

        print("{} {} {} {} {} {}".format(k, len(data[k]),
                                         tot_len,
                                         "{} {} {} {} {}".format(counts_len['1'], counts_len['0'],
                                                                 regions_len['1'],
                                                                 regions_len['0'],
                                                                 regions_len_long), tot_str,
                                         "{} {} {} {} {}".format(counts_str['1'], counts_str['0'],
                                                                 regions_str['1'],
                                                                 regions_str['0'],
                                                                 regions_str_long)))


def count_regions(state):
    # print state
    c = {}
    regions = []
    start = 0
    for i in range(len(state)-1):
        if state[i] != state[i+1]:
            c.setdefault(state[i], 0)
            c[state[i]] += 1
            regions.append((start, i, state[i]))
            start = i + 1
    # Last region
    c.setdefault(state[-1], 0)
    c[state[-1]] += 1
    regions.append((start, i + 1, state[-1]))

    return c, regions


def blast_similarity(blastfile, new_targets, old_targets):
    data = {}
    with open(blastfile) as f:
        for line in f:
            a, b, identity, al_len, mismatch, gapopen, \
            qstart, qend, start, send, evalue, bitscore = line.strip().split()

            l_a = len(new_targets[a]['sequence']) if new_targets.get(a) else len(
                old_targets[a]['sequence'])
            l_b = len(new_targets[b]['sequence']) if new_targets.get(b) else len(
                old_targets[b]['sequence'])
            identity = 2.0 * (int(al_len) - int(mismatch)) * 100.0 / float(l_a + l_b)

            if a != b:
                data.setdefault(a, {'new': (0.0, None), 'old': (0.0, None)})
                if b in new_targets:
                    if identity > data[a]['new'][0]:
                        data[a]['new'] = (identity, b)
                else:
                    if identity > data[a]['old'][0]:
                        data[a]['old'] = (identity, b)
    return data


def write_fess_reference(proteins, outdir, label, ref_names):
    data = {}
    for k in ref_names:
        data.setdefault(k, {})
        data[k]['outf'] = open("{}/{}-{}.txt".format(outdir, label, k), 'w')

    for disprot_id in proteins:
        protein = proteins[disprot_id]

        f_fasta = NamedTemporaryFile(mode='w', delete=False)
        f_fasta.write(">{} {}\n{}\n".format(disprot_id, protein.get('uniprot_accession'),
                                            protein.get('sequence')))
        f_fasta.close()
        # launch fess here
        out, err = fess.FessLauncher([f_fasta.name, '-b', '/home/marnec/usr/fess/bin64/']).run()
        states = (np.fromiter(
            list(zip(*list(map(lambda s: s.split(), out.decode('utf8').split('\n')[5:]))[:-1]))[-1],
            dtype=float) >= 0.5).astype(int)
        states = list(states.astype(str))

        os.remove(f_fasta.name)

        # Region information
        regions_names = []
        methods = []
        for region in protein['regions']:
            methods.append(region['method']['id'])
            func_names = [f['name'] for f in region['functions']]
            regions_names += func_names

        for k in ref_names:
            # Write to file
            data[k]['state'] = states
            if '1' in data[k]['state']:
                data[k]['outf'].write(
                    ">{} {} {} {} {}\n{}\n{}\n".format(disprot_id, protein['uniprot_accession'],
                                                       ','.join(methods), ','.join(regions_names),
                                                       ",".join(
                                                           ["{}-{}:{}".format(ele['start'],
                                                                              ele['end'],
                                                                              ele['ann']) for ele in
                                                            protein.get("mobidb", []).get(
                                                                "pdb_missing_residues", [])]),
                                                       protein['sequence'],
                                                       ''.join(data[k]['state'])))

def launch_iprscan(proteins):
    with NamedTemporaryFile(mode='w', delete=False) as f_fasta:
        for disprot_id in proteins:
            protein = proteins[disprot_id]

            f_fasta.write(">{} {}\n{}\n".format(disprot_id, protein.get('uniprot_accession'),
                                                protein.get('sequence')))

    ipr = interproscan.InterProScanLauncher([f_fasta.name, '-o', '.'])
    ipr.run(bindir='/home/marnec/lib/interproscan-5.9-50.0')

    return ipr.outfile


def get_ipr_regions(ipr_outfile):
    ipr_regions = {}
    for group in groupby(open(ipr_outfile), key=lambda l: l.split()[0]):
        acc = group[0]
        regs = list(map(lambda s: list(map(int, s.strip('\n').split('\t')[6:8])), group[1]))
        ipr_regions[acc] = regs
    return ipr_regions


def write_caid_references_gene3d(proteins, outdir, label, ref_names):
    data = {}

    ipr_annotation = launch_iprscan(proteins)
    ipr_regs = get_ipr_regions(ipr_annotation)

    for k in ref_names:
        data.setdefault(k, {})
        data[k]['outf'] = open("{}/{}-{}.txt".format(outdir, label, k), 'w')

    for disprot_id in proteins:
        protein = proteins[disprot_id]

        # Region information
        regions_names = []
        methods = []
        for region in protein['regions']:
            methods.append(region['method']['id'])
            func_names = [f['name'] for f in region['functions']]
            regions_names += func_names

        for k in ref_names:
            key = k.split('_')[0]
            data[k]['state'] = ['-' for i in range(len(protein.get('sequence')))]

            # Assign pdb residues
            for ele in protein.get("mobidb", []).get("pdb_missing_residues", []):
                if ele['ann'] == 'S':
                    # Assign structure for all type of references
                    for i in range(ele['start'] - 1, ele['end']):
                        data[k]['state'][i] = '0'
                elif ele['ann'] == 'D':
                    if key in ['pdb-missing']:
                        for i in range(ele['start'] - 1, ele['end']):
                            data[k]['state'][i] = '1'
                elif ele['ann'] == 'C':
                    pass
            # assign negatives from Gene3D
            if protein['disprot_id'] in ipr_regs:
                for r in ipr_regs[protein['disprot_id']]:
                    data[k]['state'][r[0] - 1: r[1]] = ['0'] * (r[1] - r[0] + 1)
            # Assign positive labels (overwrite structure)
            for region in protein['regions']:
                method = region['method']['id']
                func_names = [f['name'] for f in region['functions']]

                for i in range(int(region['start']) - 1, int(region['end'])):
                    if key in ['disprot-all']:
                        data[k]['state'][i] = '1'

                    if any('binding' in func_name.lower() for func_name in func_names):
                        if key in ['disprot-binding']:
                            data[k]['state'][i] = '1'

                    if any('rna' in func_name.lower() for func_name in func_names):
                        if key in ['disprot-binding', 'disprot-binding-rna']:
                            data[k]['state'][i] = '1'

                    if any('dna' in func_name.lower() for func_name in func_names):
                        if key in ['disprot-binding', 'disprot-binding-dna']:
                            data[k]['state'][i] = '1'

                    if any('protein-protein' in func_name.lower() for func_name in func_names):
                        if key in ['disprot-binding', 'disprot-binding-prot']:
                            data[k]['state'][i] = '1'

                    if any('linker' in func_name.lower() for func_name in func_names):
                        if key in ['disprot-linker']:
                            data[k]['state'][i] = '1'

                    if method == 'NMR':
                        if key in ['disprot-primary', 'disprot-nmr']:
                            data[k]['state'][i] = '1'

                    if method == 'XRAY':
                        if key in ['disprot-primary', 'disprot-xray']:
                            data[k]['state'][i] = '1'

                    if method not in ['NMR', 'XRAY']:
                        if key in ['disprot-secondary']:
                            data[k]['state'][i] = '1'

            # Write to file
            if '1' in data[k]['state']:
                data[k]['outf'].write(
                    ">{} {} {} {} {}\n{}\n{}\n".format(disprot_id, protein['uniprot_accession'],
                                                       ','.join(methods),
                                                       ','.join(regions_names), ",".join(
                            ["{}-{}:{}".format(ele['start'], ele['end'], ele['ann']) for ele in
                             protein.get("mobidb", []).get("pdb_missing_residues", [])]),
                                                       protein['sequence'],
                                                       ''.join(data[k]['state'])))
    os.remove(ipr_annotation)


def write_caid_references_gene3dreverse(proteins, outdir, label, ref_names):
    data = {}

    ipr_annotation = launch_iprscan(proteins)
    ipr_regs = get_ipr_regions(ipr_annotation)

    for k in ref_names:
        data.setdefault(k, {})
        data[k]['outf'] = open("{}/{}-{}.txt".format(outdir, label, k), 'w')

    for disprot_id in proteins:
        protein = proteins[disprot_id]

        # Region information
        regions_names = []
        methods = []
        for region in protein['regions']:
            methods.append(region['method']['id'])
            func_names = [f['name'] for f in region['functions']]
            regions_names += func_names

        for k in ref_names:
            data[k]['state'] = ['1' for _ in range(len(protein.get('sequence')))]

            # assign negatives from Gene3D
            if protein['disprot_id'] in ipr_regs:
                for r in ipr_regs[protein['disprot_id']]:
                    data[k]['state'][r[0] - 1: r[1]] = ['0'] * (r[1] - r[0] + 1)

            # Write to file
            if '1' in data[k]['state']:
                data[k]['outf'].write(
                    ">{} {} {} {} {}\n{}\n{}\n".format(disprot_id, protein['uniprot_accession'],
                                                       ','.join(methods),
                                                       ','.join(regions_names), ",".join(
                            ["{}-{}:{}".format(ele['start'], ele['end'], ele['ann']) for ele in
                             protein.get("mobidb", []).get("pdb_missing_residues", [])]),
                                                       protein['sequence'],
                                                       ''.join(data[k]['state'])))
    os.remove(ipr_annotation)


def write_caid_references_simple_nopdb(proteins, outdir, label, ref_names):
    data = {}
    for k in ref_names:
        data.setdefault(k, {})
        data[k]['outf'] = open("{}/{}-{}.txt".format(outdir, label, k), 'w')

    for disprot_id in proteins:
        protein = proteins[disprot_id]

        # Region information
        regions_names = []
        methods = []
        for region in protein['regions']:
            methods.append(region['method']['id'])
            func_names = [f['name'] for f in region['functions']]
            regions_names += func_names

        for k in ref_names:
            # key = k.split('_')[0]
            data[k]['state'] = ['0' for _ in range(len(protein.get('sequence')))]

            # Assign pdb residues
            for ele in protein.get("mobidb", []).get("pdb_missing_residues", []):
                if ele['ann'] == 'S':
                    for i in range(ele['start'] - 1, ele['end']):
                        data[k]['state'][i] = '-'

            # Assign positive labels (overwrite structure)
            for region in protein['regions']:
                # method = region['method']['id']
                # func_names = [f['name'] for f in region['functions']]

                for i in range(int(region['start']) - 1, int(region['end'])):
                    data[k]['state'][i] = '1'

            # Write to file
            if '1' in data[k]['state']:
                data[k]['outf'].write(
                    ">{} {} {} {} {}\n{}\n{}\n".format(disprot_id, protein['uniprot_accession'],
                                                       ','.join(methods), ','.join(regions_names),
                                                       ",".join(
                                                           ["{}-{}:{}".format(ele['start'],
                                                                              ele['end'],
                                                                              ele['ann']) for ele in
                                                            protein.get("mobidb", []).get(
                                                                "pdb_missing_residues", [])]),
                                                       protein['sequence'],
                                                       ''.join(data[k]['state'])))


def write_caid_references_pdbreverse(proteins, outdir, label, ref_names):
    data = {}
    for k in ref_names:
        data.setdefault(k, {})
        data[k]['outf'] = open("{}/{}-{}.txt".format(outdir, label, k), 'w')

    for disprot_id in proteins:
        protein = proteins[disprot_id]

        # Region information
        regions_names = []
        methods = []
        for region in protein['regions']:
            methods.append(region['method']['id'])
            func_names = [f['name'] for f in region['functions']]
            regions_names += func_names

        for k in ref_names:
            # initialize all positives
            data[k]['state'] = ['1' for _ in range(len(protein.get('sequence')))]

            # assign negatives from pdb
            for ele in protein.get("mobidb", []).get("pdb_missing_residues", []):
                if ele['ann'] == 'S':
                    for i in range(ele['start'] - 1, ele['end']):
                        data[k]['state'][i] = '0'

            # Write to file
            if '1' in data[k]['state']:
                data[k]['outf'].write(
                    ">{} {} {} {} {}\n{}\n{}\n".format(disprot_id, protein['uniprot_accession'],
                                                       ','.join(methods), ','.join(regions_names),
                                                       ",".join(
                                                           ["{}-{}:{}".format(ele['start'],
                                                                              ele['end'],
                                                                              ele['ann']) for ele in
                                                            protein.get("mobidb", []).get(
                                                                "pdb_missing_residues", [])]),
                                                       protein['sequence'],
                                                       ''.join(data[k]['state'])))


########################################################################################################################
# UniRef statistic
# egrep -o "UniRef90_[A-Z,0-9]*\b" | sort | uniq -c | sort -n

# Entries as originally annotated by curators
# prot = parse_entries("../data/curated_entries/entries_curators_disprot8.json")

# Entries fixed by Andras
prot = parse_entries("/home/marnec/Projects/CAID/caid/data/entries_curators_20181122.json")
logging.info("parsed entries %i", len(prot))

good_prot = filter_bad_entries(prot)
logging.info("good proteins %i", len(good_prot))

good_prot = filter_polyproteins(good_prot)
logging.info("removed polyproteins %i", len(good_prot))

new_prot = filter_date(copy.copy(good_prot), date(2017, 1, 1), date(2019, 1, 1))
logging.info("new proteins %i", len(new_prot))

old_prot = filter_date(copy.copy(good_prot), date(2015, 1, 1), date(2017, 1, 1))
logging.info("old proteins %i", len(old_prot))


# TODO check/rewrite this function
# get_identical_regions(prot, "../data/curated_entries/same_regions_problems.txt")

# Seqeunces for CAID predictors
# write_fasta(prot, "../data/curated_entries/sequences_all.fasta")

# Dataset fro Patrick
# write_sentences(prot, '../data/sentences.txt')


# CAID ######################
# Types of reference
# reference_names = ['disprot-all', 'pdb_missing',
#                    'disprot-xray', 'disprot-nmr', 'disprot-primary', 'disprot-secondary',
#                    'disprot-linker',
#                    'disprot-binding', 'disprot-binding-dna', 'disprot-binding-rna', 'disprot-binding-prot']

reference_names = ['disprot-all',
                   'disprot-xray', 'disprot-nmr', 'disprot-primary', 'disprot-secondary',
                   'disprot-linker',
                   'disprot-binding', 'pdb-missing']

reference_names_pdb = list(map(lambda s: s + '_pdb', reference_names))
reference_names_gene3d = list(map(lambda s: s + '_gene3d', reference_names))
reference_names_gene3d_reverse = list(map(lambda s: s + '_gene3d', reference_names))
reference_names_simple = list(map(lambda s: s + '_simple', reference_names))
reference_names_counted = list(map(lambda s: s + '_evidence-num', reference_names))

# write_caid_references_simple(new_prot, '../data/disorder', 'new', reference_names_simple)
# write_caid_references_simple(old_prot, '../data/disorder', 'old', reference_names_simple)
# write_caid_references_pdb(new_prot, '../data/disorder', 'new', reference_names_pdb)
# write_caid_references_pdb(old_prot, '../data/disorder', 'old', reference_names_pdb)
# write_caid_references_gene3d(new_prot, '../data/disorder', 'new', reference_names_gene3d)
# write_caid_references_gene3d(old_prot, '../data/disorder', 'old', reference_names_gene3d)

write_caid_references_gene3dreverse(new_prot, '../data/disorder', 'new', ['gene3d-r'])
write_caid_references_pdbreverse(new_prot, '/home/marnec/Projects/CAID/caid/data/disorder', 'new', ['disprot-all_pdb-r'])
# write_caid_references_simple_nopdb(new_prot, '/home/marnec/Projects/CAID/data/disorder', 'new', ['disprot-all_simple-nopdb'])
# write_caid_references_simpleevidencenum(new_prot, '/home/marnec/Projects/CAID/data/disorder', 'new', ['disprot-all_simple-evidence-num'])

# write_fess_reference(new_prot, '../data/disorder', 'new', ['fess'])


# caid_reference_stats('../data/reference', 'new', reference_names)
# caid_reference_stats('../data/reference', 'old', reference_names)


# similarity =  blast_similarity("../data/blast/blast.out", new_prot, old_prot)
#
# for disprot_id in similarity:
#
#     # All
#     # print disprot_id, "new" if new_prot.get(disprot_id) else "old", similarity[disprot_id]['new'][0] if similarity[disprot_id]['new'][0] > similarity[disprot_id]['old'][0] else similarity[disprot_id]['old'][0]
#
#     # New vs new
#     if new_prot.get(disprot_id):
#         print(disprot_id, similarity[disprot_id]['new'][0])


    # # New vs old
    # if new_prot.get(disprot_id):
    #     print disprot_id, similarity[disprot_id]['old'][0]


