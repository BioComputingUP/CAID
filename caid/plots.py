import os
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy.stats import ttest_ind

from caid import parse_args_plots, parse_config, set_logger, load_names

# Get dir where this piece of code is
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))


def barplots(args):
    csv = pd.read_csv(args.csv, sep='\t', index_col=[0, 1])

    baseline = csv[args.metricName].filter(regex='DB[0-9]').sort_values(ascending=False)
    results = csv[args.metricName].filter(regex='D|B[0-9]').sort_values(ascending=False)

    index = 'pred' if args.codes is False else 'code'

    baseline.index = baseline.index.get_level_values('pred')
    results.index = results.index.get_level_values(index)

    if args.secondaryCsv is not None:
        scsv = pd.read_csv(args.secondaryCsv, sep='\t', index_col=[0, 1])
        sbaseline = scsv[args.metricName].filter(regex='DB[0-9]')
        sresults = scsv[args.metricName].filter(regex='D|B[0-9]')
        sbaseline.index = sbaseline.index.get_level_values('pred')
        sresults.index = sresults.index.get_level_values(index)
        baseline = pd.concat([baseline.to_frame(), sbaseline.to_frame()],
                             axis=1, sort=False, join_axes=[baseline.index])
        baseline.columns = [args.metricName, smetric_name]
        results = pd.concat([results.to_frame(), sresults.to_frame()],
                            axis=1, sort=False, join_axes=[results.index])
        results.columns = [args.metricName, smetric_name]

    max_baseline = baseline[args.metricName].max()
    data = pd.concat([results, baseline])

    results_color = ['grey'] * len(results)
    baseline_color = ['black'] * len(baseline)

    legend_handles = [Patch(facecolor='grey', edgecolor='grey',
                            label='Predictors {}'.format(args.metricName)),
                      Patch(facecolor='black', edgecolor='black',
                            label='Baseline {}'.format(args.metricName))]

    fig, ax = plt.subplots()
    colors = results_color + baseline_color
    ylim = (-0.25, 1) if 'mcc' in args.metricName else (0, 1)
    bar = data[args.metricName].plot(kind='bar', color=colors, ylim=ylim, ax=ax,
                                     linewidth=1, edgecolor=colors, legend=True)

    if args.codes is True:
        if len(baseline) != 0:
            lbls = [c if names[c]['cons'] is False else '*' + c for c in
                    data.index[:-len(baseline)]] \
                   + list(data.index[-len(baseline):])
            if args.hideAfter > 0:
                lbls[:-len(baseline)] = [c if i > args.hideAfter
                                         else names[c if c[0] != '*' else c[1:]]['name']
                                         for i, c in enumerate(lbls, 1)]
        else:
            lbls = [c if names[c]['cons'] is False else '*' + c for c in data.index]
            if args.hideAfter > 0:
                lbls = [c if i > args.hideAfter else names[c if c[0] != '*' else c[1:]]['name']
                        for i, c in enumerate(lbls, 1)]

        ax.set_xticklabels(lbls)

    if args.secondaryCsv is not None and smetric_name:
        ax.axhline(max_baseline, linestyle='--', color='black', linewidth=0.5)
        if ylim[0] < 0:
            ax.axhline(0, color='k', linewidth=0.5)

        w = bar.patches[0].get_width()
        for x, y in zip(map(lambda patch: patch.get_x(), bar.patches), data[smetric_name]):
            ax.add_line(Line2D((x, x + w), (y, y), linewidth=1.5, color='yellowgreen'))

        legend_handles.append(Line2D([0], [0], color='yellowgreen', lw=2, label='Default cutoff'))

    ax.legend(handles=legend_handles, prop={'size': 7})
    plt.tight_layout()
    fig.savefig('{}_{}.png'.format(outbname, args.metricName), dpi=300)

    # ls ../results/codes/disorder/*.csv | sed -e 's/_redefined_results.csv//' -e 's/_results.csv//' | sort -u | while read file; do python3 plots.py ${file}_redefined_results.csv -s ${file}_results.csv -d -m mcc; done


def megamatrix(csvfile, codes=True):
    csv = pd.read_csv(csvfile, sep='\t', index_col=[0, 1])
    positive_rankers = csv[['acc', 'mcc', 'precision', 'recall', 'f1',  'precision_n', 'recall_n', 'f1_n', 'avg_acc', 'avg_mcc', 'avg_precision', 'avg_recall', 'avg_f1', 'avg_precision_n', 'avg_recall_n', 'avg_f1_n']]
    negative_rankers = csv[['std_acc', 'std_mcc', 'std_precision', 'std_recall', 'std_f1', 'std_precision_n', 'std_recall_n', 'std_f1_n']]
    p = positive_rankers.rank(axis=0, method='max', ascending=False)
    n = negative_rankers.rank(axis=0, method='max')

    # df = pd.concat([p, n], axis=1)
    df = p
    c = 'code' if codes is True else 'pred'
    indcol = df.apply(np.mean, axis=1).sort_values().index.get_level_values(c)
    # indcol = df.index.get_level_values(c)
    dff = pd.DataFrame(index=indcol, columns=indcol, dtype=float)
    del dff.index.name

    ii = 1 if codes is True else 0
    for indx, rowx in sorted(df.iterrows(), key=lambda r: r[1].mean(), reverse=True):
        for indy, rowy in sorted(df.iterrows(), key=lambda r: r[1].mean(), reverse=True):
            dff.loc[indx[ii]][indy[ii]] = ttest_ind(rowx, rowy).pvalue

    fig, ax = plt.subplots(figsize=(10, 10))
    sns.heatmap(dff, cmap='Blues', annot=True, linecolor='w', xticklabels=True, yticklabels=True, ax=ax, annot_kws={'size': 5}, fmt='.2f',)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig('{}_ranking.png'.format(os.path.splitext(os.path.basename(csvfile))[0]), dpi=300)


if __name__ == '__main__':
    args = parse_args_plots(SCRIPT_DIR)
    conf = parse_config(args.conf)
    set_logger(args.log, args.logLevel)
    smetric_name = '{}_s'.format(args.metricName) if args.secondaryCsv else None
    names = load_names(os.path.join(SCRIPT_DIR, 'caid_names.txt'))

    args.csv = os.path.abspath(args.csv)
    args.secondaryCsv = os.path.abspath(args.secondaryCsv) if args.secondaryCsv else None

    outdir = args.outputDir if args.outputDir else os.path.dirname(args.csv)
    outbname = os.path.join(outdir, os.path.splitext(os.path.basename(args.csv))[0])

    # barplots(args)
    megamatrix(args.csv, args.codes)