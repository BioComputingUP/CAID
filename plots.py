import math
import logging
import warnings
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
from itertools import product
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, pearsonr

warnings.filterwarnings("ignore")


def plot_metric_to_threshold(metric, dataset_perthr_metric, default_thr, outdir, outbase):
    logging.info("Plotting metrics progress with threshold: {}; {}".format(metric, outbase))

    s = math.ceil(len(dataset_perthr_metric.index)**.5)

    axes = dataset_perthr_metric.T.sort_index().plot(subplots=True, layout=(s, s), figsize=(s*2.5, s*2.5), sharex=True, sharey=True)
    for ax in axes.ravel():
        ax.axvline(float(default_thr[ax.get_legend_handles_labels()[1]]), linestyle="--")
    plt.gcf().suptitle("{} progress with threshold".format(metric.upper()), y=0.9)
    plt.savefig(outdir / "{}{}ToThr.png".format(outbase + "_" if outbase else outbase, metric), dpi=dpi, bbox_inches="tight")
    plt.close()


def plot_pertarget_permethod_heatmap(metric, target_metrics_preds, outdir, outbase):
    logging.info("Plotting target heatmap: {}; {}".format(metric, outbase))

    fig, ax = plt.subplots(figsize=(20, 8))
    tgt_pred_metrics_avg = target_metrics_preds[metric].unstack().mean().sort_values()
    ax = sns.heatmap(target_metrics_preds[metric].unstack().reindex(tgt_pred_metrics_avg.index, axis=1), ax=ax)
    ax2 = ax.twiny()
    ax2.set_xticks(ax.xaxis.get_ticklocs())
    ax2.set_xticklabels(tgt_pred_metrics_avg.values[[list(tgt_pred_metrics_avg.index).index(l.get_text()) for l in
                                                     ax.xaxis.get_ticklabels()]].round(3))
    ax2.tick_params(axis='x', rotation=90)

    plt.savefig(outdir / "{}tgtheatmap_{}.png".format(outbase + "_" if outbase else outbase, metric), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_icontent_correlation(metric, predictions, cons, pdbr, gene3dr, outdir, outbase):
    logging.info("Plotting idcontent correlation: {}; {}".format(metric, outbase))

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
    anymetric = resdir / "new-disprot-all_simple.analysis.all.dataset.default.metrics.csv"
    anymetric = pd.read_csv(anymetric, index_col=0)

    fig, ax = plt.subplots(figsize=(8, 8))
    sns.heatmap(anymetric.reindex(anymetric.corr().mean().sort_values().drop('thr').index, axis=1).corr(),
                     cmap='Blues', cbar=False, annot=True, ax=ax)

    plt.savefig(outdir / "metrics_corr.png", dpi=dpi, bbox_inches="tight")
    plt.close()


def plot_metrics_clustermap(resdir, outdir):
    anymetric = resdir / "new-disprot-all_simple.analysis.all.dataset.default.metrics.csv"
    anymetric = pd.read_csv(anymetric, index_col=0)
    sns.clustermap(anymetric.drop('thr', axis=1).corr(), metric="correlation")
    plt.savefig(outdir / "metrics_cluster.png", dpi=dpi, bbox_inches="tight")
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

    plt.savefig(outdir / "{}ranking_{}.png".format(outbase + "_" if outbase else outbase, metric), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_dataset_target_metric(metric, dat_metr_preds, dat_metr_bases, tgt_metr_preds, tgt_metr_bases, bts_ci_preds, outdir, outbase=""):
    logging.info("Plotting metric barplot: {}, {}".format(metric, outbase))
    fig, axes = plt.subplots(1, 2, figsize=(15, 6), sharey=True)

    negatives, optm = outbase.split("_")[1:]
    optm = optm[3:]


    # plot dataset metrics on left subplot
    dat_m = dat_metr_preds[metric].sort_values(ascending=False).append(
        dat_metr_bases[metric].sort_values(ascending=False))

    ax = dat_m.plot.bar(
        ax=axes[0],
        color=['silver'] * len(dat_metr_preds) + ['grey'] * len(dat_metr_bases),
        yerr=bts_ci_preds.xs(metric, level=1)[["lo", 'hi']].reindex(dat_m.index)
    )
    ax.axhline(dat_metr_bases[metric].max())
    ax.set_title("Dataset {}; Negatives: {}; Optimized on: {}".format(metric.upper(), negatives, optm.upper()))
    ax.set_ylabel(metric.upper())

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
    ax.set_title("Average target {}; Negatives: {}; Optimized on: {}".format(metric.upper(), negatives, optm.upper()))
    ax.set_ylabel(metric.upper())

    plt.savefig(outdir /"{}bar_{}.png".format(outbase+"_" if outbase else outbase, metric), dpi=dpi, bbox_inches="tight")
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
    ax.set_xlabel("FPR")
    ax.set_ylabel("TPR")
    ax.legend(lines, ["{} - AUC: {}".format(*t) for t in auc_rocs])
    plt.savefig(outdir /"{}roc.png".format(outbase+"_" if outbase else outbase), dpi=dpi, bbox_inches="tight")
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
    auc_pr = sorted(prcs.index.droplevel([2, 3]).unique(), key=lambda t: t[1], reverse=True)
    aps_pr = sorted(prcs.index.droplevel([1, 3]).unique(), key=lambda t: t[1], reverse=True)
    sorter = auc_pr if sortby.lower() == "auc" else aps_pr
    lines = ax.plot(*prcs.reindex(list(zip(*sorter))[0], level=0).values)
    ax.legend(lines, ["{} - {}: {}".format(t[0], sortby.upper(), t[1]) for t in sorter])
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    plt.savefig(outdir / "{}pr.{}.png".format(outbase + "_" if outbase else outbase, sortby), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(
        prog='caid-plots', description="automate plots for CAID analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('resultDir', help='directory where CAID predictors results are saved')

    parser.add_argument('baselineDir', help="directory where CAID baselines results are saved")

    parser.add_argument('referenceDir', help="directory where refernces are stored")

    parser.add_argument('dataDir', help="directory where data is stored")

    parser.add_argument('-o', '--outputDir', default='.',
                        help='directory where the output will be written (default: cwd)')

    parser.add_argument('-d', '--dpi', default=75, help='figures dpi')

    parser.add_argument('-l', '--log', type=str, default=None, help='log file')
    parser.add_argument("-ll", "--logLevel", default="ERROR",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help='log level filter. All levels <= choice will be displayed')

    args = parser.parse_args()
    return args


def set_logger(logfile, level):
    handlers = list()
    #log_formatter = logging.Formatter('%(asctime)s | %(module)-13s | %(levelname)-8s | %(message)s')
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

    axes = pd.concat([newvsnew, newvsold], axis=1).hist(edgecolor='k', figsize=(15, 5))
    axes = axes.ravel()
    line = axes[0].axvline(float(newvsnew.median()), linestyle="--", color='darkorange')
    l = axes[0].legend([line], ["Median: {:.1f} %".format(float(newvsnew.median()))])
    line = axes[1].axvline(float(newvsold.median()), linestyle="--", color='darkorange')
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

    

if __name__ == "__main__":
    args = parse_args()
    set_logger(args.log, args.logLevel)
    # obtain from cli arg
    dpi = args.dpi
    # resultdir = Path("/home/marnec/Projects/CAID/caid/results")
    # baselndir = Path("/home/marnec/Projects/CAID/caid/baseline")
    # outputdir = Path("/home/marnec/Projects/CAID/caid/plots")
    # refdir = Path("/home/marnec/Projects/CAID/caid/data/disorder")

    resultdir = Path(args.resultDir)
    baselndir = Path(args.baselineDir)
    outputdir = Path(args.outputDir)
    refdir = Path(args.referenceDir)
    datadir = Path(args.dataDir)

    # DON'T CHANGE THE ORDER
    basetypes = ["cons", "naive-new-pdb-r_simple", "naive-new-gene3d-r_simple",  # naive
                 # "raagendom", "fixedposfrc", "shuffledataset", "shuffletargets"]   # random
                 "fixedposfrc"]

    # plot_metrics_correlation(resultdir, outputdir)
    # plot_metrics_clustermap(resultdir, outputdir)

    cons_nvn = pd.read_csv(datadir / "blast_distribution_new_vs_new.txt", index_col=0)
    cons_nvo = pd.read_csv(datadir / "blast_distribution_new_vs_old.txt", index_col=0)
    cons_tot = pd.read_csv(datadir / "blast_distribution.txt", index_col=0)
    counts = pd.read_csv(datadir / "reference.csv", index_col=[0,1], header=[0,1,2])

    # plot_subset_redundancy(cons_nvn, cons_nvo, outputdir)
    # plot_dataset_redundancy(cons_tot, outputdir)
    # plot_dataset_counts(counts, outputdir)

    # iterate over file in dir (foreach reference)
    for reference in refdir.glob("*.txt"):
        logging.info(reference)
    # reference = "/home/marnec/Projects/CAID/caid/data/disorder/new-disprot-all_simple.txt"

        reference = Path(reference)
        refname = reference.stem
        
        if refname in ["new-disprot-linker_pdb", "new-disprot-linker_gene3d", "new-gene3d-r_simple", "new-pdb-r_simple"]:
            continue

        if refname != "new-disprot-all_pdb":
            continue

        roc_preds_f = resultdir / "{}.analysis.all.dataset._.roc.csv".format(refname)
        roc_preds = pd.read_csv(roc_preds_f, index_col=[0], header=[0, 1, 2])
        roc_bases = [pd.read_csv(baselndir / "{}.{}.all.dataset._.roc.csv".format(refname, b), index_col=[0], header=[0, 1, 2]) for b in basetypes[:3]]

        pr_preds_f = resultdir / "{}.analysis.all.dataset._.pr.csv".format(refname)
        pr_preds = pd.read_csv(pr_preds_f, index_col=[0], header=[0, 1, 2, 3])
        pr_bases = [pd.read_csv(baselndir / "{}.{}.all.dataset._.pr.csv".format(refname, b), index_col=[0], header=[0, 1, 2, 3]) for b in basetypes[:3]]

        # plot_roc(roc_preds, roc_bases, outputdir, refname)
        plot_pr(pr_preds, pr_bases, outputdir, refname, sortby="auc")
        plot_pr(pr_preds, pr_bases, outputdir, refname, sortby="aps")

        dataset_metrics_default_f = resultdir / "{}.analysis.all.dataset.default.metrics.csv".format(refname)
        dataset_metrics_default = pd.read_csv(dataset_metrics_default_f, index_col=0)

        dataset_metrics_single_preds = pd.concat([pd.read_csv(f, index_col=[0, 1]) for f in resultdir.glob(Path(refname).stem + ".analysis.D*dataset*")], sort=False)

        for optimized_metric in ["default"] + list(dataset_metrics_default.columns):
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
                        plot_pertarget_permethod_heatmap(metric_to_plot, target_metrics_preds,
                                                         outputdir, "{}_opt{}".format(refname, optimized_metric))
                        # plot_methdod_correlation(metric_to_plot, target_metrics_preds, outputdir, "{}_opt{}".format(refname, optimized_metric))

                plot_average_overall_ranking(optimized_metric, dataset_metrics_preds, dataset_metrics_bases, outputdir, refname)
                # plot_icontent_correlation(optimized_metric, predictions, *naive_preds, outputdir, refname)
                if optimized_metric != "default":
                    pass
                    plot_metric_to_threshold(optimized_metric, dataset_metrics_single_preds.xs(optimized_metric, level=1), dataset_metrics_default["thr"], outputdir, refname)

    # for body end

