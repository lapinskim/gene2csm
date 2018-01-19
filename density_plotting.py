import seaborn as sns
import statistics as stat
from decimal import Decimal
from collections import defaultdict
import matplotlib.pyplot as plt
import gene2csm
# for hiding print
import os
import sys

sns.set()
sns.set_style('ticks')
sns.set_context("talk")
# instead of f.tight_layout()
# plt.rcParams.update({'figure.autolayout': True})

# change the default False, if boarders are needed
# plt.rcParams.update({"patch.force_edgecolor": False})


# temporary solution before the logging is implemented
class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self._original_stdout


def create_plots(entropy_list, suptitle=None, style=None):
    '''
    Create figure from density plots
    '''

    # create a grid with two columns if more than one plot is plotted
    n_items = len(entropy_list)
    n_cols = 2
    if n_items % 2 == 0:
        n_rows = n_items / 2
    else:
        n_rows = (n_items // 2) + 1
    # adjust the figure size accordingly to the number of plots in the figure
    plt.figure(figsize=(14, n_rows * 5)).set_tight_layout(True)
    # create all axes and figures
    index = 1
    for item_id, item_ent in entropy_list:
        plt.subplot(n_rows, n_cols, index)
        draw_plot(item_ent, title=item_id)
        index += 1
    # Add a super title if present
    if suptitle:
        if style:
            plt.suptitle(suptitle, style=style, y=1.06, fontsize=18)
        else:
            plt.suptitle(suptitle)
    plt.show()
    return


def draw_plot(entropy, title=None):
    '''
    Draw density plot and calculate the statistics
    for the provided entropy values
    '''

    from matplotlib.offsetbox import AnchoredText

    # calculate statistics
    totaln = len(entropy)
    erange = round(Decimal(max(entropy)) - Decimal(min(entropy)), 2)
    mean = round(Decimal(stat.mean(entropy)), 2)
    pstdev = round(Decimal(stat.pstdev(entropy)), 2)

    # draw density plot,
    # add kernel density line later for further customisation
    ax = sns.distplot(entropy,
                      rug=False,
                      bins=50,
                      kde=False,
                      hist_kws={"color": "black",
                                "alpha": 1,
                                "histtype": "bar",
                                "edgecolor": "white",
                                "linewidth": 0})
    # annotate x axis
    ax.set_xlabel('entropy')
    # annotate the first y axis
    ax.set_ylabel('hist (counts)')
    # create the second y axis
    ax2 = ax.twinx()
    # draw the kernel density line with a white outline
    sns.kdeplot(entropy,
                linewidth=5,
                color="white",
                alpha=1,
                bw=0.2,
                ax=ax2)
    sns.kdeplot(entropy,
                linewidth=2,
                color="black",
                alpha=1,
                bw=0.2,
                ax=ax2)
    # annotate the second y axis
    ax2.set_ylabel('kde (density)')
    # add some styling
    sns.despine(right=False, offset=10, trim=True)

    # create annotation frame
    at = AnchoredText("total number = {}\nrange = {}\nμ = {}\nσ = {}"
                      .format(totaln,
                              erange,
                              mean,
                              pstdev),
                      loc=2, prop=dict(size=12), frameon=True)
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)

    # add title
    if title:
        plt.title(title)
    return ax


def files2plot(fasta_fn, paths, crRNA_len=36, GC_limit=(0, 100), proc=1,
               supstyle=None):
    '''
    Create density plots from the provided files.
    '''

    file_dict = dict(zip([e.rsplit('.')[0] for e in fasta_fn], paths))

    for key in file_dict:
        energy_dict = defaultdict(list)
        # calculate for each sequence in the fasta file

        # hide the prints for plotting- temporary solution
        with HiddenPrints():
            for seq in gene2csm.parse_fasta(file_dict[key]):
                seq_id, result = gene2csm.estimate_energy_input(seq,
                                                                crRNA_len,
                                                                GC_limit,
                                                                proc=7)
                entropy = [float(round(e[3], 3)) for e in result]
                energy_dict[key].append((seq_id, entropy))
        create_plots(energy_dict[key], suptitle=key, style=supstyle)
    return
