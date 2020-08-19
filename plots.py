import json
import math
import logging
import matplotlib.ticker as ticker
import warnings
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
from itertools import product
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, pearsonr, hmean

warnings.filterwarnings("ignore")


def plot_metric_to_threshold(metric, dataset_perthr_metric, default_thr, opt_threshold, outdir, outbase, names=None):
    logging.info("Plotting metrics progress with threshold: {}; {}".format(metric, outbase))

    if names is not None:
        default_thr.rename(names, inplace=True)
        opt_threshold.rename(names, inplace=True)


    s = math.ceil(len(dataset_perthr_metric)**.5)
    fig, axes = plt.subplots(s, s, sharex=True, sharey=True, figsize=(s*2.5, s*2.5))

    def df_sorter(df):
        if names is None:
            sorter = df.index.get_level_values(0)[0]
        else:
            sorter = names(df.index.get_level_values(0)[0])
        return sorter

    for i, (ax, p) in enumerate(zip(axes.ravel(), sorted(dataset_perthr_metric, key=df_sorter))):
        p.columns = p.columns.values.astype(np.float)
        p = p.sort_index(axis=1)
        method = p.index.droplevel(1).unique()[0]
        method = names(method) if names is not None else method
        p = p.reset_index(level=0, drop=True)
        # x = np.array(dataset_perthr_metric.columns.values, dtype=np.float)
        p.loc[metric].plot(ax=ax, label=None)
        # dataset_perthr_metric.loc[p].sort_index().plot(ax=ax, legend=default_thr[p])
        if i == 0 or i % s == 0:
            ax.set_ylabel(metric.upper())
        if i > s**2-s:
            ax.set_xlabel("Thresholds")

        ax.set_title(method)
        ax.axvline(default_thr[method], linestyle="--", label="Default: {}".format(default_thr[method]), color="orange")
        ax.axvline(opt_threshold[method], linestyle="--", label="Optimal: {}".format(opt_threshold[method]), color="deeppink")
        ax.legend(*zip(*list(zip(*ax.get_legend_handles_labels()))[1:]))

    for x in range(-1, -(s**2 - len(dataset_perthr_metric))-1, -1):
        fig.delaxes(axes.ravel()[x])

    plt.gcf().suptitle("{} progress with threshold".format(metric.upper()), y=1.01)
    fig.tight_layout()
    plt.savefig(outdir / "{}{}ToThr.png".format(outbase + "_" if outbase else outbase, metric), dpi=dpi, bbox_inches="tight")
    plt.close()



def plot_pertarget_permethod_heatmap(metric, target_metrics_preds, target_metrics_bases, outdir, outbase, names=None):
    logging.info("Plotting target heatmap: {}; {}".format(metric, outbase))

    fig, ax = plt.subplots(figsize=(20, 8))

    tgt_metric = target_metrics_preds[metric].unstack().append(target_metrics_bases[metric].unstack())
    if names is not None:
        tgt_metric.rename(names, inplace=True)

    tgt_pred_metrics_avg_col = tgt_metric.mean().sort_values()

    tgt_pred_metrics_avg_row = tgt_metric.mean(axis=1).sort_values(ascending=False)
    ax.set_facecolor('royalblue')
    ax = sns.heatmap(tgt_metric.reindex(tgt_pred_metrics_avg_col.index, axis=1).reindex(tgt_pred_metrics_avg_row.index), ax=ax)
    ax2 = ax.twiny()
    ax2.set_xticks(ax.xaxis.get_ticklocs())
    ax2.set_xticklabels(tgt_pred_metrics_avg_col.values[[list(tgt_pred_metrics_avg_col.index).index(l.get_text()) for l in
                                                     ax.xaxis.get_ticklabels()]].round(3))
    ax2.tick_params(axis='x', rotation=90)

    plt.savefig(outdir / "{}tgtheatmap_{}.png".format(outbase + "_" if outbase else outbase, metric), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_icontent_correlation(metric, predictions, cons, pdbr, gene3dr, outdir, outbase, names=None):
    logging.info("Plotting idcontent correlation: {}; {}".format(metric, outbase))

    predictions = pd.DataFrame({**predictions.to_dict(), **cons.to_dict(), **pdbr.to_dict(), **gene3dr.to_dict()}).dropna()
    n = len(predictions.columns.get_level_values(0).unique())
    n = int(n ** .5) if (n ** .5) % 2 == 0 else int(math.ceil(n ** .5))

    if names is not None:
        predictions.rename(names, axis=1, level=0, inplace=True)

    idcontent_ref = predictions.iloc[:, 0].groupby(level=0).mean()

    fig, axes = plt.subplots(n, n, figsize=(int(n)*2.5, int(n)*2.5), sharey="row")
    axes = axes.flatten()
    for ax, p in zip(axes, predictions.columns.get_level_values(0).unique()):
        x = predictions[(p, "states")].groupby(level=0).mean()
        y = idcontent_ref
        sns.scatterplot(x, y, ax=ax, label="Pearson R = {:.3f}".format(pearsonr(x, y)[0]))
        ax.legend(loc="upper left")
        # ax.set_title(p)
        ax.set_ylabel("Reference")
        ax.set_xlabel(p)
    fig.tight_layout()

    plt.savefig(outdir / "{}icontentcorr_{}.png".format(outbase + "_" if outbase else outbase, metric), dpi=dpi, bbox_inches="tight")
    plt.close()


def plot_methdod_correlation(metric, target_metrics_preds, outdir, outbase):
    logging.info("Plotting methods correlation: {}; {}".format(metric, outbase))

    t = target_metrics_preds[metric].unstack().reindex(
        target_metrics_preds[metric].groupby(level=0).mean().sort_values(ascending=False).index)
    ax = sns.pairplot(t.T)
    plt.savefig(outdir / "{}methodcorr_{}.png".format(outbase + "_" if outbase else outbase, metric), dpi=dpi, bbox_inches="tight")
    plt.close()


def plot_metrics_correlation(resdir, outdir):
    logging.info("Plotting metrics correlation")
    anymetric = resdir / "new-disprot-all_simple.analysis.all.dataset.default.metrics.csv"
    anymetric = pd.read_csv(anymetric, index_col=0)

    fig, ax = plt.subplots(figsize=(8, 8))
    sns.heatmap(anymetric.reindex(anymetric.corr().mean().sort_values().drop('thr').index, axis=1).corr(),
                     cmap='Blues', cbar=False, annot=True, ax=ax)

    plt.savefig(outdir / "metrics_corr.png", dpi=dpi, bbox_inches="tight")
    plt.close()


def plot_metrics_clustermap(resdir, outdir):
    logging.info("Plotting metrics clustermap")
    anymetric = resdir / "new-disprot-all_simple.analysis.all.dataset.default.metrics.csv"
    anymetric = pd.read_csv(anymetric, index_col=0)
    sns.clustermap(anymetric.drop('thr', axis=1).corr(), metric="correlation")
    plt.savefig(outdir / "metrics_cluster.png", dpi=dpi, bbox_inches="tight")
    plt.close()


def plot_average_overall_ranking(metric, metrics_preds, metrics_bases, outdir, outbase, plotfirst=None, names=None, level='target'):
    logging.info("Plotting ranking: {}".format(metric))
    if level == 'target':
        metrics_preds = metrics_preds.groupby(level=0).mean()
        metrics_bases = metrics_bases.groupby(level=0).mean()



    metrics_selection = ['bac', 'f1s', 'fpr', 'mcc', 'ppv', 'tpr', 'tnr']

    pred_ranking = metrics_preds.append(metrics_bases)[metrics_selection]
    pred_ranking = pred_ranking.rank(axis=0, method='max', ascending=False, na_option='bottom')

    pred_ranking = pred_ranking.reindex(pred_ranking.mean(axis=1).sort_values().index)

    if names is not None:
        pred_ranking.rename(names, inplace=True)

    if plotfirst is not None:
        pred_ranking = pred_ranking.head(10)

    cartesian_product = product(pred_ranking.index, pred_ranking.index)
    dat_pred_ranking_pvals = [[*couple, ttest_ind(*pred_ranking.loc[[*couple]].values, equal_var=True)[1]] for
                              couple in cartesian_product]
    dat_pred_ranking_pvals = pd.DataFrame(dat_pred_ranking_pvals).set_index([0, 1]).unstack()
    dat_pred_ranking_pvals.columns = dat_pred_ranking_pvals.columns.droplevel(0)
    dat_pred_ranking_pvals = dat_pred_ranking_pvals.reindex(index=pred_ranking.index,
                                                            columns=pred_ranking.index).round(2)

    fig, axes = plt.subplots(figsize=(len(dat_pred_ranking_pvals) / 2, len(dat_pred_ranking_pvals) / 2))
    ax = sns.heatmap(dat_pred_ranking_pvals, annot=True, mask=np.triu(dat_pred_ranking_pvals), cmap="Reds", cbar=False,
                     center=0.1, vmin=0, vmax=0.2, edgecolor='w')
    ax = sns.heatmap(dat_pred_ranking_pvals, annot=True, mask=np.tril(dat_pred_ranking_pvals), cmap="Greens",
                     cbar=False, center=0.1, vmin=0, vmax=0.2, edgecolor='w')

    plt.savefig(outdir / "{}ranking.{}_{}{}.png".format(outbase + "_" if outbase else outbase, level,
                                                   metric, "_best{}".format(plotfirst) if plotfirst else ""),
                dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_cput_stacked_bars(ax, predictors, names=None):
    cputime = pd.read_csv('../data/dataset_stats/cpu_time.csv', index_col=0, header=[0, 1])
    order = ['prediction', 'hhblits', 'psiblast']
    ax = cputime.mean().unstack()[order].rename(names).reindex(predictors).plot.bar(ax=ax, stacked=True, log=True)

    ax.errorbar(range(len(predictors)),
                cputime.mean().unstack()[order].rename(names).reindex(predictors).sum(axis=1),
                cputime.std().unstack()[order].rename(names).reindex(predictors).apply(
                    lambda r: np.sqrt(np.sum(np.log(r) ** 2)) / len(r), axis=1),
                linewidth=0, elinewidth=2, c='k')
    ax.set_ylabel('$\log_{10}(Seconds)$')
    #     ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: '$10^{{{0}}}$'.format(int(x))))
    ax.set_xlabel(None)


def plot_cput_boxplot(ax, predictors, names):
    cputime = pd.read_csv('../data/dataset_stats/cpu_time.csv', index_col=0, header=[0, 1])
    order = ['prediction', 'hhblits', 'psiblast']

    cputime = cputime.groupby(level=0, axis=1).sum().rename(names, axis=1).reindex(predictors, axis=1)
    cputime.boxplot(rot=90,
                    ax=ax,
                    grid=False,
                    flierprops=dict(marker=',', markerfacecolor='steelblue', markeredgecolor='none', alpha=.1),
                    boxprops=dict(alpha=.7))

    cputimelessthanone = pd.Series(np.nan, index=cputime.median().index)
    cputimelessthanone[cputime.quantile(.75) < 1] = 1
    ax.plot(np.arange(len(cputimelessthanone)) + 1, cputimelessthanone.values, marker='o', markersize=5,
            color='magenta', linestyle='None')
    ax.set_yscale('log')


    ax.set_ylabel('$\log_{10}(Seconds)$')
    # ax.text(-0.15, 1.05, plotlbl, transform=ax.transAxes, size=20, weight='bold')
    ax.set_xlabel(None)
    ax.set_yticks(10 ** np.linspace(0, 4, 5))
    ax.grid(which='major', axis='both', alpha=.1)
    ax.set_ylim(0.5, 10 ** 4)


def plot_dat_tgt_metric_cpuspeed(metric, dat_metr_preds, dat_metr_bases, tgt_metr_preds, tgt_metr_bases, bts_ci_preds, outdir, outbase="", names=None):
    logging.info("Plotting metric barplot: {}, {}".format(metric, outbase))
    fig = plt.figure(figsize=(15, 6))
    gs = fig.add_gridspec(2, 2, height_ratios=[2, 1])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[0, 1], sharey=ax1)
    ax4 = fig.add_subplot(gs[1, 1], sharey=ax2)

    negatives, optm = outbase.split("_")
    optm = optm[3:]

    # plot dataset metrics on left subplot
    dat_m = dat_metr_preds[metric].sort_values(ascending=False).append(
        dat_metr_bases[metric].sort_values(ascending=False))
    yerr = bts_ci_preds.xs(metric, level=1)[["lo", 'hi']].reindex(dat_m.index)
    colors = ['silver'] * len(dat_metr_preds) + ['grey'] * len(dat_metr_bases)

    if names is not None:
        dat_m = dat_m.rename(names)
        yerr = yerr.rename(names)

    dat_m.plot.bar(ax=ax1, color=colors, yerr=yerr)
    # plot_cput_stacked_bars(ax2, dat_m.index, names)

    plot_cput_boxplot(ax2, dat_m.index, names)

    ax1.axhline(dat_metr_bases[metric].max())
    ax1.set_title("{}; Dataset {}; Optimized on: {}".format(dataset_names[negatives], metric.upper(), optm.upper()))
    ax1.set_ylabel(metric.upper())
    ax1.set_xticks([])

    # plot target metrics on right subplot
    n = tgt_metr_preds.groupby(level=0).count().append(tgt_metr_bases.groupby(level=0).count())
    tgt_m = tgt_metr_preds.groupby(level=0).mean()[metric].sort_values(ascending=False).append(
        tgt_metr_bases.groupby(level=0).mean()[metric].sort_values(ascending=False))
    yerr = tgt_metr_preds.groupby(level=0).std()[metric].append(
        tgt_metr_bases.groupby(level=0).std()[metric])[tgt_m.index] / n[metric]**0.5 / 2
    colors = ['silver'] * len(tgt_metr_preds.groupby(level=0).mean()) + ['grey'] * len(tgt_metr_bases.groupby(level=0).mean())

    if names is not None:
        tgt_m = tgt_m.rename(names)
        yerr = yerr.rename(names)

    tgt_m.plot.bar(ax=ax3, color=colors, yerr=yerr)
    plot_cput_boxplot(ax4, tgt_m.index, names)

    ax3.axhline(tgt_metr_bases.groupby(level=0).mean()[metric].max())
    ax3.set_title("{}; Target {}; Optimized on: {}".format(dataset_names[negatives], metric.upper(), optm.upper()))
    ax3.set_ylabel(metric.upper())
    ax3.set_xticks([])

    plt.savefig(outdir / "{}bar_{}.png".format(outbase+"_" if outbase else outbase, metric), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_dataset_target_metric(metric, dat_metr_preds, dat_metr_bases, tgt_metr_preds, tgt_metr_bases, bts_ci_preds, outdir, outbase="", names=None):
    logging.info("Plotting metric barplot: {}, {}".format(metric, outbase))
    fig, axes = plt.subplots(1, 2, figsize=(15, 6), sharey=True)

    negatives, optm = outbase.split("_")
    optm = optm[3:]

    # plot dataset metrics on left subplot
    dat_m = dat_metr_preds[metric].sort_values(ascending=False).append(
        dat_metr_bases[metric].sort_values(ascending=False))
    yerr = bts_ci_preds.xs(metric, level=1)[["lo", 'hi']].reindex(dat_m.index)
    colors = ['silver'] * len(dat_metr_preds) + ['grey'] * len(dat_metr_bases)

    if names is not None:
        dat_m = dat_m.rename(names)
        yerr = yerr.rename(names)

    ax = dat_m.plot.bar(ax=axes[0], color=colors, yerr=yerr)

    ax.axhline(dat_metr_bases[metric].max())
    ax.set_title("Dataset {}; Negatives: {}; Optimized on: {}".format(metric.upper(), negatives, optm.upper()))
    ax.set_ylabel(metric.upper())

    # plot target metrics on right subplot
    n = tgt_metr_preds.groupby(level=0).count().append(tgt_metr_bases.groupby(level=0).count())
    tgt_m = tgt_metr_preds.groupby(level=0).mean()[metric].sort_values(ascending=False).append(
        tgt_metr_bases.groupby(level=0).mean()[metric].sort_values(ascending=False))
    yerr = tgt_metr_preds.groupby(level=0).std()[metric].append(
        tgt_metr_bases.groupby(level=0).std()[metric])[tgt_m.index] / n[metric]**0.5 / 2
    colors = ['silver'] * len(tgt_metr_preds.groupby(level=0).mean()) + ['grey'] * len(tgt_metr_bases.groupby(level=0).mean())

    if names is not None:
        tgt_m = tgt_m.rename(names)
        yerr = yerr.rename(names)

    ax = tgt_m.plot.bar(ax=axes[1], color=colors, yerr=yerr)

    ax.axhline(tgt_metr_bases.groupby(level=0).mean()[metric].max())
    ax.set_title("Average target {}; Negatives: {}; Optimized on: {}".format(metric.upper(), negatives, optm.upper()))
    ax.set_ylabel(metric.upper())

    plt.savefig(outdir / "{}bar_{}.png".format(outbase+"_" if outbase else outbase, metric), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_roc(preds_rocs, cons_roc, pdb_roc, gene3d_roc, random_rocs, outdir, outbase, names=None):
    # procs, proc, groc, dataset, title, names = names, croc = None
    logging.info("Plotting roc")
    fig, ax = plt.subplots(figsize=(7.5, 7.5))
    preds_rocs.rename(names, axis=1, level=0)
    auc_rocs = sorted(preds_rocs.columns.droplevel(2).unique(), key=lambda t: t[1], reverse=True)[:10]
    rocs = preds_rocs.reindex(list(zip(*auc_rocs))[0], axis=1, level=0)
    ax.plot([0, 1], [0, 1], color='k', linestyle='--')

    for p in rocs.columns.get_level_values(0).unique()[:10]:
        ax.plot(*rocs[p].dropna().T.values, label=p)

    ax.plot(*cons_roc.dropna().T.values, label="Naive Conservation", color='silver', linewidth=2)

    for n, m in zip([pdb_roc, gene3d_roc], ['o', 's']):
        idx = n.index.values - 0.5
        ax.plot(*n.iloc[idx[idx > 0].argmin()], marker=m, markeredgecolor='silver', markeredgewidth=2,
                markerfacecolor='w', markersize=10,
                c='w', label=names(n.columns.get_level_values(0)[0]))

    for rr, m in zip(random_rocs, ['*', 'P', 'd']):
        ax.plot(*rr[['fpr', 'tpr']].mean(), marker=m, markeredgecolor='silver', markeredgewidth=2,
                markerfacecolor='w', markersize=10, c='w', label=names(rr.index.get_level_values(0)[0]))
    ax.legend()
    lhandles, llabels = ax.get_legend_handles_labels()
    if cons_roc is not None:
        auc_rocs.extend(cons_roc.rename(names, axis=1, level=0).columns.droplevel(2).unique().tolist())
    pwauc = next(zip(*auc_rocs))
    ax.legend(lhandles,
              ['{}'.format('{} - AUC: {}'.format(names(l), auc_rocs[pwauc.index(l)][1]) if l in pwauc else l) for l in
               llabels],
              bbox_to_anchor=(0., -.37, 1., .102), loc='lower left', ncol=2, mode="expand", borderaxespad=0)


    ax.set_xlabel("FPR")
    ax.set_ylabel("TPR")
    plt.savefig(outdir /"{}roc.png".format(outbase+"_" if outbase else outbase), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_pr(pred_prs, cons_prc, pdb_prc, gene3d_prc, random_prcs, cov, outdir, outbase, sortby="auc", names=None):
    logging.info("Plotting precision-recall curve")
    fig, ax = plt.subplots(figsize=(10.5, 7))

    auc_pr = sorted(pred_prs.columns.droplevel([2, 3]).unique(), key=lambda t: t[1], reverse=True)
    aps_pr = sorted(pred_prs.columns.droplevel([1, 3]).unique(), key=lambda t: t[1], reverse=True)
    fmax_pr = sorted(((p, hmean(pred_prs[p].dropna().values, axis=1).max().round(2)) for p in pred_prs.columns.get_level_values(0).unique()),
                     key=lambda t: t[1], reverse=True)

    sorter = None
    if sortby == 'auc':
        sorter = auc_pr
    elif sortby == 'aps':
        sorter = aps_pr
    elif sortby == 'fmax':
        sorter = fmax_pr

    # select first 10 predictors (based on AUC)
    prcs = pred_prs.reindex(list(zip(*sorter))[0], axis=1, level=0)

    # plot f-score level lines
    r = np.linspace(0, 1, 1000)
    fs = hmean(np.array(np.meshgrid(r, r)).T.reshape(-1, 2), axis=1).reshape(1000, 1000)
    cs = plt.contour(r, r, fs, levels=np.linspace(0.1, 1, 10), colors='silver', alpha=0.7, linewidths=1, linestyles='--')
    ax.clabel(cs, inline=True, fmt='%.1f', fontsize=10, manual=[(l, l) for l in cs.levels[:-1]])

    # plot predictor lines and markers
    for p, z in zip(prcs.columns.get_level_values(0).unique()[:10], range(5, 55, 5)):
        fmax_idx = hmean(prcs[p].dropna().T).argmax()
        lines = ax.plot(*prcs[p].dropna().T.values, label=p, zorder=55-z)
        ax.plot(*prcs[p].dropna().T.values[:, fmax_idx], color='w', marker='o', markerfacecolor=lines[0].get_color(), markersize=10, zorder=55-z)
        ax.plot(*prcs[p].dropna().T.values[:, fmax_idx], color='w', marker='o', markerfacecolor=lines[0].get_color(), zorder=55-z)

    # plot naive conservation
    ax.plot(*cons_prc.dropna().T.values, label="Naive Conservation", color='k', linewidth=1, zorder=5)
    cov['Naive Conservation'] = 1

    # plot naives
    for n, m in zip([pdb_prc, gene3d_prc], ['o', 's']):
        nname = n.columns.get_level_values(0)[0]
        cov[nname] = 1
        ax.plot(*n.loc[1.0], marker=m, markeredgecolor='k', markeredgewidth=1,
                markerfacecolor='w', markersize=8, zorder=60,
                c='w', label=names(nname))

    # plot randoms
    for rprc, m in zip(random_prcs, ['*', 'P', 'd']):
        rname = rprc.index.get_level_values(0)[0]
        cov[rname] = 1
        ax.plot(*rprc[['tpr', 'ppv']].mean(), marker=m, markeredgecolor='k', markeredgewidth=1, zorder=60,
                markerfacecolor='w', markersize=8, c='w', label=names(rname))

    cov = cov.to_dict()
    ax.legend()
    lhandles, llabels = ax.get_legend_handles_labels()
    if cons_roc is not None:
        sorter.extend(cons_roc.rename(names, axis=1, level=0).columns.droplevel(2).unique().tolist())
    pwauc = next(zip(*sorter))
    ax.legend(lhandles,
              ['{}'.format('{} - {}: {}, C: {:.2f}'.format(names(l), sortby.upper(), sorter[pwauc.index(l)][1], cov[l]) if l in pwauc else l) for l in llabels],
              loc='center left', bbox_to_anchor=(1, 0.5))#, mode="expand", borderaxespad=0)

    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    fig.tight_layout()
    plt.savefig(outdir / "{}pr.{}.png".format(outbase + "_" if outbase else outbase, sortby), dpi=dpi, bbox_inches="tight")
    plt.close(fig)

def plot_cputime_to_performance(metric, tgt_pred_metrics, outdir, outbase, names=None):
    logging.info('plotting cputime to {}'.format(metric))
    fig, ax = plt.subplots(figsize=(8, 6))

    cputime = pd.read_csv('../data/dataset_stats/cpu_time.csv', header=[0, 1], index_col=[0]).groupby(level=0, axis=1).sum()
    tgt_pred_metrics = tgt_pred_metrics[metric]
    y = tgt_pred_metrics.groupby(level=0).mean().sort_values(ascending=False)
    x = np.log10(cputime).mean().reindex(y.index).replace([np.inf, -np.inf], -2)

    ax = sns.scatterplot(x=x, y=y, hue=y.index, zorder=50)

    ax.errorbar(x=x,
                y=y,
                yerr=tgt_pred_metrics.groupby(level=0).std() / (len(tgt_pred_metrics) / len(y.index)) ** 0.5,
                xerr=np.log10(cputime).std().reindex(y.index),
                linewidth=0, elinewidth=0.5, c='k', capsize=2, capthick=0.5)

    ax.axvline(0, linestyle='--', linewidth=1, label='1 Second')
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '$10^{{{}}}$'.format(int(x))))

    _ = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)

    pointstoname = dict(zip(*ax.get_legend_handles_labels()[::-1]))
    sorter = y.sort_values(ascending=False)

    ax.legend([pointstoname[k] for k in sorter.index],
              ["{} (F={:.2f})".format(names(k), v) for k, v in sorter.iteritems()],
              bbox_to_anchor=(1, .5),
              loc='center left',
              ncol=2)

    ax.set_ylabel('$F_{max}$')
    ax.set_xlabel('Seconds')
    plt.savefig(outdir / "{}cputime_to_{}.png".format(outbase + "_" if outbase else outbase, metric), dpi=dpi,
                bbox_inches="tight")
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(
        prog='caid-plots', description="automate plots for CAID analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('resultDir', help='directory where CAID predictors results are saved')

    parser.add_argument('baselineDir', help="directory where CAID baselines results are saved")

    parser.add_argument('referenceDir', help="directory where refernces are stored")

    parser.add_argument('datasetStatsDir', help="directory where data is stored")

    parser.add_argument('-o', '--outputDir', default='.',
                        help='directory where the output will be written (default: cwd)')

    parser.add_argument('-d', '--dpi', default=75, help='figures dpi')
    parser.add_argument('-g', '--glob', default='*.txt')
    parser.add_argument('-n', '--names', default=None, help='json file with predictors names')

    parser.add_argument('-l', '--log', type=str, default=None, help='log file')
    parser.add_argument("-ll", "--logLevel", default="ERROR",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help='log level filter. All levels <= choice will be displayed')

    args = parser.parse_args()
    return args


def set_logger(logfile, level):
    handlers = list()
    log_formatter = logging.Formatter('%(asctime)s | %(module)-13s | %(levelname)-8s | %(message)s')

    if logfile:
        file_handler = logging.FileHandler(logfile, 'a')
        file_handler.setFormatter(log_formatter)
        handlers.append(file_handler)
    else:
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(log_formatter)
        handlers.append(console_handler)

    logging.basicConfig(level=level, format=log_formatter, handlers=handlers)


def plot_subset_redundancy(newvsnew, newvsold, outdir):
    logging.info("Plotting subset redundancy")

    axes = pd.concat([newvsnew, newvsold], axis=1).hist(edgecolor='k', figsize=(15, 5), bins=20)
    axes = axes.ravel()
    line = axes[0].axvline(float(newvsnew.median()), linestyle="--", color='darkorange', linewidth=2)
    l = axes[0].legend([line], ["Median: {:.1f} %".format(float(newvsnew.median()))])
    line = axes[1].axvline(float(newvsold.median()), linestyle="--", color='darkorange', linewidth=2)
    l = axes[1].legend([line], ["Median: {:.1f} %".format(float(newvsold.median()))])
    plt.savefig(outdir / "subset_redundancy.png", dpi=dpi, bbox_inches="tight")
    plt.close(plt.gcf())


def plot_dataset_redundancy(totredund, outdir):
    logging.info("Plotting dataset redundancy")
    axes = totredund.hist(by="Age", edgecolor="k", figsize=(15, 5))
    med_new = float(totredund[totredund["Age"] == "CAID"].mean())
    med_old = float(totredund[totredund["Age"] == "DisProt7"].mean())
    line = axes[0].axvline(med_new, color="darkorange", linestyle="--")
    axes[0].legend([line], ["Mean: {:.1f} %".format(med_new)])
    line = axes[1].axvline(med_old, color="darkorange", linestyle="--")
    axes[1].legend([line], ["Mean: {:.1f} %".format(med_old)])
    plt.savefig(outdir / "dataset_redundancy.png", dpi=dpi, bbox_inches="tight")
    plt.close(plt.gcf())


def plot_dataset_counts(counts, outdir):
    logging.info("Plotting dataset counts")
    fig, axes = plt.subplots(2, 3, figsize=(15, 7), sharey=True, sharex="col")
    axes = axes.ravel()
    counts.loc["CAID"].drop("PDB missing").xs("Proteins", level=2, axis=1).plot.barh(ax=axes[0], legend=False)
    axes[0].set_title("CAID - Number of Proteins")
    counts.xs("Regions", level=1, axis=1).loc["CAID"].drop("PDB missing").drop(("PDB", "Positive"), axis=1).plot.barh(
        ax=axes[1])
    axes[1].legend(axes[1].get_legend().get_patches(), ["Postives", "Negatives (Simples)", "Negatives (PDB)"])
    axes[1].set_title("CAID - Number of Regions")

    counts.xs("Residues", level=1, axis=1).loc["CAID"].drop("PDB missing").drop(("PDB", "Positive"), axis=1).plot.barh(
        ax=axes[2])
    axes[2].legend(axes[2].get_legend().get_patches(), ["Postives", "Negatives (Simples)", "Negatives (PDB)"])
    axes[2].set_title("CAID - Number of Residues")

    counts.loc["DisProt 7"].drop("PDB missing").xs("Proteins", level=2, axis=1).plot.barh(ax=axes[3], legend=False)
    axes[3].set_title("DisProt 7 - Number of Proteins")
    counts.xs("Regions", level=1, axis=1).loc["DisProt 7"].drop("PDB missing").drop(("PDB", "Positive"),
                                                                                    axis=1).plot.barh(ax=axes[4])
    axes[4].legend(axes[4].get_legend().get_patches(), ["Postives", "Negatives (Simples)", "Negatives (PDB)"])
    axes[4].set_title("DisProt 7 - Number of Regions")

    counts.xs("Residues", level=1, axis=1).loc["DisProt 7"].drop("PDB missing").drop(("PDB", "Positive"),
                                                                                     axis=1).plot.barh(ax=axes[5])
    axes[5].legend(axes[5].get_legend().get_patches(), ["Postives", "Negatives (Simples)", "Negatives (PDB)"])
    axes[5].set_title("DisProt 7 - Number of Residues")
    plt.savefig(outdir / "dataset_counts.png", dpi=dpi, bbox_inches="tight")
    plt.close(plt.gcf())


def plot_species_count(dst, outdir):
    dst = dst.loc['disprot-disorder']
    ax = dst.data.groupby('species').count().taxon.sort_values().plot.barh(figsize=(4, 18), logx=True)
    plt.tick_params(axis='x', which='both', labeltop='on', labelbottom='on', top='on')
    ax.vlines([1, 10, 100], 0, len(dst), linestyle='--', linewidth=1, color='silver', zorder=0)
    plt.savefig(outdir / "species_counts.png", dpi=dpi, bbox_inches='tight')
    plt.close(plt.gcf())

def get_names(fname):
    names = json.load(open(args.names))

    name = names.get(fname)
    fname = fname.lower()
    if name is None:
        if "cons" in fname:
            name = "Naive Conservation"
        elif "pdb" in fname and 'reverse' in fname:
            name = "Naive PDB"
        elif "gene3d" in fname and 'reverse' in fname:
            name = "Naive Gene3D"
        elif "random" in fname:
            name = "Random"
        elif "dataset" in fname:
            name = "Shuffled dataset"
        elif "target" in fname:
            name = "Shuffled targets"
        elif "fix" in fname:
            name = "Fixed ID content"
        elif "ref" in fname:
            name = "Reference"
    return name
    
dataset_names = {
    'disprot-disorder': 'DisProt',
    'disprot-disorder-pdb-atleast': 'DisProt-PDB',
    'disprot-binding': 'DisProt-Binding',
    'disprot-binding-all': 'DisProt-Binding-All',
    'disprot-binding-disorder': 'DisProt-Binding-Disorder'
}


if __name__ == "__main__":
    args = parse_args()
    set_logger(args.log, args.logLevel)
    logging.getLogger('matplotlib').setLevel(logging.WARNING)

    dpi = args.dpi
    resultdir = Path(args.resultDir)
    baselndir = Path(args.baselineDir)
    outputdir = Path(args.outputDir)
    refdir = Path(args.referenceDir)
    datadir = Path(args.datasetStatsDir)

    get_names = get_names if args.names is not None else None

    # DON'T CHANGE THE ORDER
    basetypes = ["cons", "naive-pdb-atleast-reverse", "naive-gene3d-reverse",   # naive
                 # "fixedposfrc",
                 "random", "shuffledataset", "shuffletargets"]   # random

    plot_metrics_correlation(resultdir, outputdir)
    plot_metrics_clustermap(resultdir, outputdir)

    cons_nvn = pd.read_csv(datadir / "blast_distribution_new_vs_new.txt", index_col=0)
    cons_nvo = pd.read_csv(datadir / "blast_distribution_new_vs_old.txt", index_col=0)
    cons_tot = pd.read_csv(datadir / "blast_distribution.txt", index_col=0)
    counts = pd.read_csv(datadir / "reference.csv", index_col=[0, 1], header=[0, 1, 2])
    refstats_target = pd.read_csv('../data/dataset_stats/references-stats.target.csv', index_col=[0, 1], header=[0, 1])

    plot_subset_redundancy(cons_nvn, cons_nvo, outputdir)
    plot_dataset_redundancy(cons_tot, outputdir)
    plot_dataset_counts(counts, outputdir)
    plot_species_count(refstats_target, outputdir)

    # iterate over file in dir (foreach reference)
    for reference in refdir.glob(args.glob):
        logging.info(reference)

        reference = Path(reference)
        refname = reference.stem
        logging.debug(refname)

        roc_preds_f = resultdir / "{}.analysis.all.dataset._.roc.csv".format(refname)
        roc_preds = pd.read_csv(roc_preds_f, index_col=[0], header=[0, 1, 2])
        roc_bases = [pd.read_csv(baselndir / "{}.{}.all.dataset._.roc.csv".format(refname, b), index_col=[0],
                                 header=[0, 1, 2]) for b in basetypes[1:3]]
        cons_roc = pd.read_csv('../baseline/{}.cons.all.dataset._.roc.csv'.format(refname), index_col=[0], header=[0, 1, 2])
        roc_random_bases = [pd.read_csv('../baseline/{}.{}.all.target.mcc.metrics.csv'.format(refname, b), index_col=0) for b in basetypes[3:]]

        pr_preds_f = resultdir / "{}.analysis.all.dataset._.pr.csv".format(refname)
        pr_preds = pd.read_csv(pr_preds_f, index_col=[0], header=[0, 1, 2, 3])
        pr_bases = [pd.read_csv(baselndir / "{}.{}.all.dataset._.pr.csv".format(refname, b), index_col=[0],
                                header=[0, 1, 2, 3]) for b in basetypes[1:3]]
        cons_pr = pd.read_csv('../baseline/{}.cons.all.dataset._.pr.csv'.format(refname), index_col=[0],
                               header=[0, 1, 2, 3])
        pr_random_bases = [pd.read_csv('../baseline/{}.{}.all.target.mcc.metrics.csv'.format(refname, b), index_col=0)
                            for b in basetypes[3:]]

        plot_roc(roc_preds, cons_roc, *roc_bases, roc_random_bases, outputdir, refname, names=get_names)
        coverage = pd.read_csv(resultdir / '{}.analysis.all.target.default.metrics.csv'.format(refname), index_col=[0,1], header=[0,1])
        coverage = (coverage.groupby(level=0).count().max(axis=1) / np.max(coverage.groupby(level=0).count().values))
        plot_pr(pr_preds, cons_pr, *pr_bases, pr_random_bases, coverage, outputdir, refname, sortby="auc", names=get_names)
        plot_pr(pr_preds, cons_pr, *pr_bases, pr_random_bases, coverage, outputdir, refname, sortby="aps", names=get_names)
        plot_pr(pr_preds, cons_pr, *pr_bases, pr_random_bases, coverage, outputdir, refname, sortby="fmax", names=get_names)

        dataset_metrics_default_f = resultdir / "{}.analysis.all.dataset.default.metrics.csv".format(refname)
        dataset_metrics_default = pd.read_csv(dataset_metrics_default_f, index_col=0)

        codestart = 'D' if 'disorder' in str(refdir) else 'B'
        dataset_metrics_single_preds = [pd.read_csv(f, index_col=[0, 1]) for f in resultdir.glob("{}.analysis.{}*dataset*".format(Path(refname).stem, codestart))]

        for optimized_metric in ["default"] + list(dataset_metrics_default.columns):
            if optimized_metric not in {"aucroc", "aucpr", "thr", "aps"}:

                predictions = pd.read_csv(resultdir / "{}.analysis.all.dataset._.predictions.csv".format(refname), index_col=[0, 1], header=[0, 1])
                naive_preds = [pd.read_csv(baselndir / "{}.{}.all.dataset._.predictions.csv".format(refname, n), index_col=[0, 1], header=[0, 1]) for n in basetypes[:3]]

                dataset_metrics_preds_f = resultdir / "{}.analysis.all.dataset.{}.metrics.csv".format(refname, optimized_metric)
                dataset_metrics_preds = pd.read_csv(dataset_metrics_preds_f, index_col=0)

                bootstr_ci_preds_f = resultdir / "{}.analysis.all.ci.{}.metrics.csv".format(refname, optimized_metric)
                bootstr_ci_preds = pd.read_csv(bootstr_ci_preds_f, index_col=[0, 1])
                target_metrics_preds_f = resultdir / "{}.analysis.all.target.{}.metrics.csv".format(refname, optimized_metric)
                target_metrics_preds = pd.read_csv(target_metrics_preds_f, index_col=[0, 1])

                dataset_metrics_bases = []
                target_metrics_bases = []
                for bt in basetypes:
                    if bt in basetypes[3:]:
                        target_metrics_base_f = baselndir / "{}.{}.avg.target.{}.metrics.csv".format(refname, bt, optimized_metric)
                        target_metrics_bases.append(pd.read_csv(target_metrics_base_f, index_col=[0, 1], header=[0, 1]).xs('mean', level=1, axis=1))
                        dataset_metrics_base_f = baselndir / "{}.{}.avg.dataset.{}.metrics.csv".format(refname, bt, optimized_metric)
                        dataset_metrics_bases.append(pd.read_csv(dataset_metrics_base_f, index_col=0).loc["mean"].to_frame(bt).T)
                    else:
                        target_metrics_base_f = baselndir / "{}.{}.all.target.{}.metrics.csv".format(refname, bt, optimized_metric)
                        target_metrics_bases.append(pd.read_csv(target_metrics_base_f, index_col=[0, 1]))
                        dataset_metrics_base_f = baselndir / "{}.{}.all.dataset.{}.metrics.csv".format(refname, bt, optimized_metric)
                        dataset_metrics_bases.append(pd.read_csv(dataset_metrics_base_f, index_col=[0]))


                target_metrics_bases = pd.concat(target_metrics_bases)
                dataset_metrics_bases = pd.concat(dataset_metrics_bases)

                for metric_to_plot in list(dataset_metrics_default.columns):
                    if metric_to_plot not in {"aucroc", "aucpr", "thr", "aps"}:
                        plot_dataset_target_metric(metric_to_plot, dataset_metrics_preds, dataset_metrics_bases,
                                                   target_metrics_preds, target_metrics_bases, bootstr_ci_preds, outputdir,
                                                   "{}_opt{}".format(refname, optimized_metric), names=get_names)

                        plot_dat_tgt_metric_cpuspeed(metric_to_plot, dataset_metrics_preds, dataset_metrics_bases,
                                                   target_metrics_preds, target_metrics_bases, bootstr_ci_preds, outputdir,
                                                   "{}_opt{}".format(refname, optimized_metric), names=get_names)

                        plot_pertarget_permethod_heatmap(metric_to_plot, target_metrics_preds, target_metrics_bases,
                                                         outputdir, "{}_opt{}".format(refname, optimized_metric), names=get_names)
                        plot_methdod_correlation(metric_to_plot, target_metrics_preds, outputdir, "{}_opt{}".format(refname, optimized_metric))

                plot_average_overall_ranking(optimized_metric, target_metrics_preds, target_metrics_bases, outputdir, refname, names=get_names, level='target')
                plot_average_overall_ranking(optimized_metric, dataset_metrics_preds, dataset_metrics_bases,
                                                 outputdir, refname, names=get_names, level='dataset')
                plot_average_overall_ranking(optimized_metric, target_metrics_preds, target_metrics_bases, outputdir, refname, plotfirst=10, names=get_names, level='target')
                plot_average_overall_ranking(optimized_metric, dataset_metrics_preds, dataset_metrics_bases, outputdir, refname, plotfirst=10, names=get_names, level='dataset')
                plot_icontent_correlation(optimized_metric, predictions, *naive_preds, outputdir, refname, names=get_names)

                if optimized_metric == 'f1s':
                    plot_cputime_to_performance('f1s', target_metrics_preds, outputdir, refname, names=get_names)

                if optimized_metric != "default":
                    plot_metric_to_threshold(optimized_metric,
                                             dataset_metrics_single_preds,
                                             dataset_metrics_default["thr"].round(3),
                                             dataset_metrics_preds["thr"].round(3),
                                             outputdir, refname, names=get_names)


