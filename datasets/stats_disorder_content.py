import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.ticker import FormatStrFormatter


# Obsolete. Replaced by stats_pairwise_identity.py
# def plot_comparison():
#     df = pd.read_csv('/home/marnec/Projects/CAID/data/blast_distribution.txt', sep=' ')
#     fig, ax = plt.subplots()
#     cons_new = df[df['age'] == 'new']['cons']
#     cons_old = df[df['age'] == 'old']['cons']
#     ax.hist([cons_new, cons_old],
#             color=['grey', 'black'], label=['DisProt 8', 'DisProt 7'], bins=20)
#     median_new = cons_new.median()
#     median_old = cons_old.median()
#     ax.axvline(median_new, linestyle='--', label='Median DisProt 8 ({:.2f}%)'.format(median_new), color='yellowgreen')
#     ax.axvline(median_old, linestyle='--', label='Median DisProt 7 ({:.2f}%)'.format(median_old), color='lightblue')
#     ax.set_xlabel('Max identity')
#     ax.set_ylabel('Frequency')
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig('max_identity.png', dpi=300)
#
#
# def plot_single(fname, oname):
#     df = pd.read_csv(fname, sep=' ')
#     fig, ax = plt.subplots()
#     ax.hist(df['cons'], bins=20, color='grey', linewidth=1, edgecolor='w', label=None)
#     median = df['cons'].median()
#     ax.axvline(median, linestyle='--', label='Median ({:.2f}%)'.format(median), color='yellowgreen')
#     ax.set_xlabel('Max identity')
#     ax.set_ylabel('Frequency')
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig(oname, dpi=300)


def reference(column):
    # df = pd.read_excel('/home/marnec/Projects/CAID/data/reference.xlsx', header=[0, 1, 2], index_col=[0, 1])
    df = pd.read_csv('../data/reference.csv', header=[0, 1, 2], index_col=[0, 1])
    # print(df.loc['DisProt 8'])
    fig, ax = plt.subplots()
    w = 0.2

    if column != 'Proteins':
        # p = df.loc['DisProt 8']['Lenient'][column].plot.bar(stacked=True, ax=ax, width=w, rot=80,
        #                                                        edgecolor='k', color=['lightgrey', 'w'])
        lenient = df.loc['DisProt 8']['Lenient'][column]
        p = range(len(lenient['Positive']))

        ax.bar(np.array(p), lenient['Positive long'],
               width=w, label='Lenient Positive long', edgecolor='k', color='lightgrey', hatch='//')

        ax.bar(np.array(p), lenient['Positive'] - lenient['Positive long'],
               width=w, label='Lenient Positive short', edgecolor='k', color='lightgrey' ,bottom=lenient['Positive long'])

        ax.bar(np.array(p), lenient['Negative'],
               width=w, label='Lenient Negative', color='w', edgecolor='k', bottom=lenient['Positive'])

        ax.bar(np.array(p) + w * 1.5, df.loc['DisProt 8']['Strict'][column]['Positive long'],
               width=w, label='Strict Positive long', color='coral', edgecolor='k', hatch='//')

        ax.bar(np.array(p) + w * 1.5, df.loc['DisProt 8']['Strict'][column]['Positive'] - df.loc['DisProt 8']['Strict'][column]['Positive long'],
               width=w, label='Strict Positive short', color='coral', edgecolor='k', bottom=df.loc['DisProt 8']['Strict'][column]['Positive long'])

        ax.bar(np.array(p) + w * 1.5, df.loc['DisProt 8']['Strict'][column]['Negative'],
               width=w, bottom=df.loc['DisProt 8']['Strict'][column]['Positive'], label='Strict Negative', color='grey', edgecolor='k')

        ax.set_xticklabels(['',] + list(lenient.index), rotation=80)
        l = plt.legend()
        # l.get_texts()[0].set_text('Lenient Positive')
        # l.get_texts()[1].set_text('Lenient Negative')
        ax.set_xlim(-0.5, len(p) - 0.2)

    else:
        df.loc['DisProt 8']['X']['Y'][column].plot.bar(color='lightgrey', edgecolor='k',
                                                       rot=80, width=w)
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=True)
    ax.set_ylabel(column)
    ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
    plt.tight_layout()
    plt.savefig('../data/disorder_content_{}.png'.format(column), dpi=300)


def reference_residues(column):
    # df = pd.read_excel('/home/marnec/Projects/CAID/data/reference.xlsx', header=[0, 1, 2], index_col=[0, 1])
    df = pd.read_csv('../data/reference.csv', header=[0, 1, 2], index_col=[0, 1])
    # print(df.loc['DisProt 8'])
    fig, ax = plt.subplots()
    w = 0.2

    if column != 'Proteins':
        # p = df.loc['DisProt 8']['Lenient'][column].plot.bar(stacked=True, ax=ax, width=w, rot=80,
        #                                                        edgecolor='k', color=['lightgrey', 'w'])
        lenient = df.loc['DisProt 8']['Lenient'][column]
        strict = df.loc['DisProt 8']['Strict'][column]
        p = range(len(lenient['Positive']))

        ax.bar(np.array(p), lenient['Positive'],
               width=w, label='Lenient Positive', edgecolor='k', color='lightgrey')

        ax.bar(np.array(p), lenient['Negative'],
               width=w, label='Lenient Negative', color='w', edgecolor='k', bottom=lenient['Positive'])

        ax.bar(np.array(p) + w * 1.5, df.loc['DisProt 8']['Strict'][column]['Positive'],
               width=w, label='Strict Positive', color='coral', edgecolor='k')

        ax.bar(np.array(p) + w * 1.5, df.loc['DisProt 8']['Strict'][column]['Negative'],
               width=w, bottom=df.loc['DisProt 8']['Strict'][column]['Positive'], label='Strict Negative', color='grey', edgecolor='k')

        ax.set_xticklabels(['',] + list(lenient.index), rotation=80)
        l = plt.legend()
        # l.get_texts()[0].set_text('Lenient Positive')
        # l.get_texts()[1].set_text('Lenient Negative')
        ax.set_xlim(-0.5, len(p) - 0.2)

    else:
        df.loc['DisProt 8']['X']['Y'][column].plot.bar(color='lightgrey', edgecolor='k',
                                                       rot=80, width=w)
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=True)
    ax.set_ylabel(column)
    ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
    plt.tight_layout()
    plt.savefig('../data/disorder_content_{}.png'.format(column), dpi=300)




if __name__ == '__main__':
    # plot_comparison()
    # plot_single('/home/marnec/Projects/CAID/data/blast_distribution_new_vs_new.txt', 'max_identity_nvn.png')
    # plot_single('/home/marnec/Projects/CAID/data/blast_distribution_new_vs_old.txt', 'max_identity_nvo.png')
    # reference('Regions')
    # reference_residues('Residues')
    # reference('Proteins')

    df = pd.read_csv('../data/disorder_content.csv', header=[0], index_col=[0])

