from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import product
from scipy.stats import ttest_ind, pearsonr
import numpy as np
import math
import warnings
warnings.filterwarnings("ignore")


def plot_metric_to_threshold(metric, dataset_perthr_metric, outdir, outbase):
    s = math.ceil(len(dataset_perthr_metric.index)**.5)
    dataset_perthr_metric.T.plot(subplots=True, layout=(s, s), figsize=(s*2.5, s*2.5), sharex=True, sharey=True)
    plt.gcf().suptitle(metric)
    plt.savefig(outdir / "{}{}ToThr.png".format(outbase + "_" if outbase else outbase, metric), dpi=DPI, bbox_inches="tight")
    plt.close()


def plot_pertarget_permethod_heatmap(metric, target_metrics_preds, outdir, outbase):
    fig, ax = plt.subplots(figsize=(20, 8))
    tgt_pred_metrics_avg = target_metrics_preds[metric].unstack().mean().sort_values()
    ax = sns.heatmap(target_metrics_preds[metric].unstack().reindex(tgt_pred_metrics_avg.index, axis=1), ax=ax)
    ax2 = ax.twiny()
    ax2.set_xticks(ax.xaxis.get_ticklocs())
    ax2.set_xticklabels(tgt_pred_metrics_avg.values[[list(tgt_pred_metrics_avg.index).index(l.get_text()) for l in
                                                     ax.xaxis.get_ticklabels()]].round(3))
    ax2.tick_params(axis='x', rotation=90)

    plt.savefig(outdir / "{}tgtheatmap_{}.png".format(outbase + "_" if outbase else outbase, metric), dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def plot_icontent_correlation(metric, predictions, cons, pdbr, gene3dr, outdir, outbase):

    predictions = pd.DataFrame({**predictions.to_dict(), **cons.to_dict(), **pdbr.to_dict(), **gene3dr.to_dict()}).dropna()
    n = len(predictions.columns.get_level_values(0).unique())
    n = int(n ** .5) if (n ** .5) % 2 == 0 else int(math.ceil(n ** .5))

    idcontent_ref = predictions[("ref", "states")].groupby(level=0).mean()
    fig, axes = plt.subplots(n, n, figsize=(int(n), int(n)))
    axes = axes.flatten()
    for ax, p in zip(axes, predictions.columns.get_level_values(0).unique()):
        x = predictions[(p, "states")].groupby(level=0).mean()
        y = idcontent_ref
        sns.scatterplot(x, y, ax=ax, label=round(pearsonr(x, y)[0], 3))
        ax.legend(loc="upper left")
        ax.set_ylabel("Reference")
        ax.set_xlabel(p)

    plt.savefig(outdir / "{}icontentcorr_{}.png".format(outbase + "_" if outbase else outbase, metric), dpi=DPI, bbox_inches="tight")
    plt.close()


def plot_methdod_correlation(metric, target_metrics_preds, outdir, outbase):
    t = target_metrics_preds[metric].unstack().reindex(
        target_metrics_preds[metric].groupby(level=0).mean().sort_values(ascending=False).index)
    ax = sns.pairplot(t.T)
    plt.savefig(outdir / "{}methodcorr_{}.png".format(outbase + "_" if outbase else outbase, metric), dpi=DPI,
                bbox_inches="tight")
    plt.close()


def plot_metrics_correlation(resdir, outdir):
    anymetric = resdir / "new-disprot-all_simple.analysis.all.dataset.default.metrics.csv"
    anymetric = pd.read_csv(anymetric, index_col=0)

    fig, ax = plt.subplots(figsize=(8, 8))
    sns.heatmap(anymetric.reindex(anymetric.corr().mean().sort_values().drop('thr').index, axis=1).corr(),
                     cmap='Blues', cbar=False, annot=True, ax=ax)

    plt.savefig(outdir / "metrics_corr.png", dpi=DPI, bbox_inches="tight")
    plt.close()


def plot_metrics_clustermap(resdir, outdir):
    anymetric = resdir / "new-disprot-all_simple.analysis.all.dataset.default.metrics.csv"
    anymetric = pd.read_csv(anymetric, index_col=0)
    sns.clustermap(anymetric.drop('thr', axis=1).corr(), metric="correlation")
    plt.savefig(outdir / "metrics_cluster.png", dpi=DPI, bbox_inches="tight")
    plt.close()


def plot_average_overall_ranking(metric, metrics_preds, metrics_bases, outdir, outbase):
    dat_pred_ranking = metrics_preds.append(metrics_bases).rank(axis=0, method='max', ascending=False,
                                                                      na_option='bottom').drop("thr", axis=1)
    dat_pred_ranking = dat_pred_ranking.reindex(dat_pred_ranking.mean(axis=1).sort_values().index)

    cartesian_product = product(dat_pred_ranking.index, dat_pred_ranking.index)
    dat_pred_ranking_pvals = [[*couple, ttest_ind(*dat_pred_ranking.loc[[*couple]].values, equal_var=False)[1]] for
                              couple in cartesian_product]
    dat_pred_ranking_pvals = pd.DataFrame(dat_pred_ranking_pvals).set_index([0, 1]).unstack()
    dat_pred_ranking_pvals.columns = dat_pred_ranking_pvals.columns.droplevel(0)
    dat_pred_ranking_pvals = dat_pred_ranking_pvals.reindex(index=dat_pred_ranking.index,
                                                            columns=dat_pred_ranking.index).round(3)

    fig, axes = plt.subplots(figsize=(len(dat_pred_ranking_pvals) / 2, len(dat_pred_ranking_pvals) / 2))
    ax = sns.heatmap(dat_pred_ranking_pvals, annot=True, mask=np.triu(dat_pred_ranking_pvals), cmap="Reds", cbar=False,
                     center=0.1, vmax=0.2)
    ax = sns.heatmap(dat_pred_ranking_pvals, annot=True, mask=np.tril(dat_pred_ranking_pvals), cmap="Greens",
                     cbar=False, center=0.1, vmax=0.2)

    plt.savefig(outdir / "{}ranking_{}.png".format(outbase + "_" if outbase else outbase, metric), dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def plot_dataset_target_metric(metric, dat_metr_preds, dat_metr_bases, tgt_metr_preds, tgt_metr_bases, bts_ci_preds, outdir, outbase=""):
    fig, axes = plt.subplots(1, 2, figsize=(15, 6), sharey=True)

    # plot dataset metrics on left subplot
    dat_m = dat_metr_preds[metric].sort_values(ascending=False).append(
        dat_metr_bases[metric].sort_values(ascending=False))

    ax = dat_m.plot.bar(
        ax=axes[0],
        color=['silver'] * len(dat_metr_preds) + ['grey'] * len(dat_metr_bases),
        yerr=bts_ci_preds.xs('bac', level=1)[["lo", 'hi']].reindex(dat_m.index)
    )
    ax.axhline(dat_metr_bases[metric].max())

    # plot target metrics on right subplot
    n = target_metrics_preds.groupby(level=0).count().append(target_metrics_bases.groupby(level=0).count())
    tgt_m = tgt_metr_preds.groupby(level=0).mean()[metric].sort_values(ascending=False).append(
        tgt_metr_bases.groupby(level=0).mean()[metric].sort_values(ascending=False))

    ax = tgt_m.plot.bar(
        ax=axes[1],
        color=['silver'] * len(tgt_metr_preds.groupby(level=0).mean()) + ['grey'] * len(
            tgt_metr_bases.groupby(level=0).mean()),
        yerr=tgt_metr_preds.groupby(level=0).std()[metric].append(
            tgt_metr_bases.groupby(level=0).std()[metric])[tgt_m.index] / n[metric]**0.5 / 2)
    ax.axhline(tgt_metr_bases.groupby(level=0).mean()[metric].max())

    plt.savefig(outdir/"{}bar_{}.png".format(outbase+"_" if outbase else outbase, metric), dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def plot_roc(pred_rocs, naive_rocs, outdir, outbase):
    fig, ax = plt.subplots(figsize=(7.5, 7.5))
    # create auc ranking of all predictor + baselines
    auc_rocs = pred_rocs
    for naive_roc in naive_rocs:
        auc_rocs = auc_rocs.join(naive_roc)
    auc_rocs = sorted(auc_rocs.columns.droplevel(2).unique(), key=lambda t: t[1], reverse=True)


    # select first 10 predictors
    rocs = pred_rocs.reindex(list(zip(*auc_rocs))[0], axis=1, level=0).T.head(20)
    # add baselines if they are not among the best 10
    for naive_roc in naive_rocs:
        rocs = rocs.append(naive_roc.T) if naive_roc.index.get_level_values(0)[0] not in rocs.index.get_level_values('predictor').unique() else rocs

    # recreate auc ranking with first 10 predictors and baselines
    auc_rocs = sorted(rocs.index.droplevel(2).unique(), key=lambda t: t[1], reverse=True)
    lines = plt.plot(*rocs.reindex(list(zip(*auc_rocs))[0], level=0).values)
    ax.plot([0, 1], [0, 1], color='k', linestyle='--')
    ax.legend(lines, ["{} - AUC: {}".format(*t) for t in auc_rocs])
    plt.savefig(outdir/"{}roc.png".format(outbase+"_" if outbase else outbase), dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def plot_pr(pred_prs, naive_prs, outdir, outbase, sortby="auc"):
    fig, ax = plt.subplots(figsize=(7.5, 7.5))

    all_prs = pred_prs
    for naive_pr in naive_prs:
        all_prs = all_prs.join(naive_pr)

    auc_pr = sorted(all_prs.columns.droplevel([2, 3]).unique(), key=lambda t: t[1], reverse=True)
    aps_pr = sorted(all_prs.columns.droplevel([1, 3]).unique(), key=lambda t: t[1], reverse=True)

    sorter = auc_pr if sortby.lower() == "auc" else aps_pr
    # select first 10 predictors (based on AUC)
    prcs = pred_prs.reindex(list(zip(*sorter))[0], axis=1, level=0).T.head(20)
    # add baselines if they are not among the best 10
    for naive_pr in naive_prs:
        prcs = prcs.append(naive_pr.T) if naive_pr.index.get_level_values(0)[0] not in prcs.index.get_level_values('predictor').unique() else prcs

    # recreate auc ranking with first 10 predictors and baselines
    sorter_pr = sorted(prcs.index.droplevel([2, 3]).unique(), key=lambda t: t[1], reverse=True)
    lines = ax.plot(*prcs.reindex(list(zip(*sorter_pr))[0], level=0).values)
    ax.legend(lines, ["{} - {}: {}".format(t[0], sortby.upper(), t[1]) for t in sorter_pr])
    plt.savefig(outdir / "{}pr.{}.png".format(outbase + "_" if outbase else outbase, sortby), dpi=DPI, bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__":
    # obtain from cli arg
    DPI = 75
    resultdir = Path("/home/marnec/Projects/CAID/caid/results")
    baselndir = Path("/home/marnec/Projects/CAID/caid/baseline")
    outputdir = Path("/home/marnec/Projects/CAID/caid/plots")
    refdir = Path("/home/marnec/Projects/CAID/caid/data/disorder")

    # DON'T CHANGE THE ORDER
    basetypes = ["cons", "naive-new-pdb-r_simple", "naive-new-gene3d-r_simple",  # naive
                 # "random", "fixedposfrc", "shuffledataset", "shuffletargets"]   # random
                 "fixedposfrc"]

    # plot_metrics_correlation(resultdir, outputdir)
    # plot_metrics_clustermap(resultdir, outputdir)

    # iterate over file in dir (foreach reference)
    for reference in refdir.glob("*.txt"):
        print(reference)
    # reference = "/home/marnec/Projects/CAID/caid/data/disorder/new-disprot-all_simple.txt"

        reference = Path(reference)
        refname = reference.stem
        
        if refname in ["new-disprot-linker_pdb", "new-disprot-linker_gene3d", "new-gene3d-r_simple", "new-pdb-r_simple"]:
            continue

        roc_preds_f = resultdir / "{}.analysis.all.dataset._.roc.csv".format(refname)
        roc_preds = pd.read_csv(roc_preds_f, index_col=[0], header=[0, 1, 2])
        roc_bases = [pd.read_csv(baselndir / "{}.{}.all.dataset._.roc.csv".format(refname, b), index_col=[0], header=[0, 1, 2]) for b in basetypes[:3]]

        pr_preds_f = resultdir / "{}.analysis.all.dataset._.pr.csv".format(refname)
        pr_preds = pd.read_csv(pr_preds_f, index_col=[0], header=[0, 1, 2, 3])
        pr_bases = [pd.read_csv(baselndir / "{}.{}.all.dataset._.pr.csv".format(refname, b), index_col=[0], header=[0, 1, 2, 3]) for b in basetypes[:3]]

        # plot_roc(roc_preds, roc_bases, outputdir, refname)
        # plot_pr(pr_preds, pr_bases, outputdir, refname, sortby="auc")
        # plot_pr(pr_preds, pr_bases, outputdir, refname, sortby="aps")

        dataset_metrics_default_f = resultdir / "{}.analysis.all.dataset.default.metrics.csv".format(refname)
        dataset_metrics_default = pd.read_csv(dataset_metrics_default_f, index_col=0)

        dataset_metrics_single_preds = pd.concat([pd.read_csv(f, index_col=[0,1]) for f in resultdir.glob(Path(refname).stem+"*D*dataset*")], sort=False)

        for optimized_metric in list(dataset_metrics_default.columns) + ["default"]:
            if optimized_metric not in {"aucroc", "aucpr", "thr", "aps"}:

                predictions = pd.read_csv(resultdir / "{}.analysis.all.dataset._.predictions.csv".format(refname), index_col=[0, 1], header=[0, 1])
                # cons_preds = pd.read_csv(resultdir / "{}.cons.all.dataset._.predictions.csv".format(refname), index_col=[0, 1], header=[0, 1])
                naive_preds = [pd.read_csv(baselndir / "{}.{}.all.dataset._.predictions.csv".format(refname, n), index_col=[0, 1], header=[0, 1]) for n in basetypes[:3]]

                dataset_metrics_preds_f = resultdir / "{}.analysis.all.dataset.{}.metrics.csv".format(refname, optimized_metric)
                dataset_metrics_preds = pd.read_csv(dataset_metrics_preds_f, index_col=0)

                bootstr_ci_preds_f = resultdir / "{}.analysis.all.bootstrap.{}.metrics.csv".format(refname, optimized_metric)
                bootstr_ci_preds = pd.read_csv(bootstr_ci_preds_f, index_col=[0, 1])
                target_metrics_preds_f = resultdir / "{}.analysis.all.target.{}.metrics.csv".format(refname, optimized_metric)
                target_metrics_preds = pd.read_csv(target_metrics_preds_f, index_col=[0, 1])

                dataset_metrics_bases = []
                target_metrics_bases = []
                for bt in basetypes:

                    if bt in basetypes[3:]:
                        target_metrics_base_f = baselndir / "{}.{}.avg.target.{}.metrics.csv".format(refname, bt, optimized_metric)
                        target_metrics_bases.append(pd.read_csv(target_metrics_base_f, index_col=0).loc["mean"].to_frame(bt).T)
                        dataset_metrics_base_f = baselndir / "{}.{}.avg.dataset.{}.metrics.csv".format(refname, bt, optimized_metric)
                        dataset_metrics_bases.append(pd.read_csv(dataset_metrics_base_f, index_col=0).loc["mean"].to_frame(bt).T)
                    else:
                        target_metrics_base_f = baselndir / "{}.{}.all.target.{}.metrics.csv".format(refname, bt, optimized_metric)
                        target_metrics_bases.append(pd.read_csv(target_metrics_base_f, index_col=0))
                        dataset_metrics_base_f = baselndir / "{}.{}.all.dataset.{}.metrics.csv".format(refname, bt, optimized_metric)
                        dataset_metrics_bases.append(pd.read_csv(dataset_metrics_base_f, index_col=0))

                target_metrics_bases = pd.concat(target_metrics_bases)
                dataset_metrics_bases = pd.concat(dataset_metrics_bases)

                for metric_to_plot in list(dataset_metrics_default.columns):
                    if metric_to_plot not in {"aucroc", "aucpr", "thr", "aps"}:
                        pass
                        # plot_dataset_target_metric(metric_to_plot, dataset_metrics_preds, dataset_metrics_bases,
                        #                            target_metrics_preds, target_metrics_bases, bootstr_ci_preds, outputdir,
                        #                            "{}_opt{}".format(refname, optimized_metric))
                        #
                        # plot_pertarget_permethod_heatmap(metric_to_plot, target_metrics_preds,
                        #                                  outputdir, "{}_opt{}".format(refname, optimized_metric))
                        # plot_methdod_correlation(metric_to_plot, target_metrics_preds, outputdir, "{}_opt{}".format(refname, optimized_metric))

                # plot_average_overall_ranking(optimized_metric, dataset_metrics_preds, dataset_metrics_bases, outputdir, refname)
                # plot_icontent_correlation(optimized_metric, predictions, *naive_preds, outputdir, refname)
                plot_metric_to_threshold(optimized_metric, dataset_metrics_single_preds.xs(optimized_metric, level=1), outputdir, refname)

    # for body end

