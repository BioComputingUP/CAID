import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def draw_per_instance_heatmap(per_instance_scores, outbasename:str = "heatmap"):
    # parse isntance scores file, set multi index
    instance_table = pd.read_csv(per_instance_scores, header=[0], index_col=[0, 1])

    # select specified metric, sort by column average
    instance_table = instance_table.xs('bal_acc', level=1)\
        .reindex(instance_table.mean().sort_values().index, axis=1)


    # automatically draw heatmap. nan values are rendered transparent
    ax = sns.heatmap(instance_table, cmap='Blues')

    # set background color to stick out over colormap (hihglight nan values)
    ax.set_facecolor('xkcd:salmon')

    # remove automatically generated x axis label
    ax.axes.get_yaxis().get_label().set_visible(False)

    # generate twin axis (ticks will go on the other side of original axis)
    twin_ax = ax.twiny()

    # get the tick locations in data coordinates as a numpy array
    ax2tick_location = ax.xaxis.get_ticklocs()

    # draw xticks of twin axis (on top) at the same location of original
    twin_ax.set_xticks(ax2tick_location)

    # set twin xticks to show mean value of column
    twin_ax.set_xticklabels(instance_table.mean().round(2))

    # rotate twin xtick vertically
    twin_ax.tick_params(axis='x', rotation=90)

    plt.savefig("{}.png".format(outbasename), bbox_inches='tight', dpi=300)
