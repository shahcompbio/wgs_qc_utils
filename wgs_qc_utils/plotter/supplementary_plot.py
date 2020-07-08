import os
import wgs_analysis.plots.snv as snv_plots
import wgs_analysis.annotation.position as position
from deconstructSigs import DeconstructSigs


def plot_rainfall(chrom, pos, axis):
    """
    rainfall plot
    :param chrom:
    :param pos:
    :param axis:
    :return:
    """
    abs_checker.check_input_is_valid([chrom, pos],
                                     [abs_checker.CheckerTypes.STRING,
                                      abs_checker.CheckerTypes.INT])

    snvs = pd.DataFrame({"chrom":chrom.tolist(), "coord":pos.tolist()})
    snvs = position.annotate_adjacent_distance(snvs)

    snv_plots.snv_adjacent_distance_plot(axis, snvs)

    return axis


def calculate_mutation_class(row):
    subs = DeconstructSigs._standardize_subs(row['ref'], row['alt'])
    trinuc = DeconstructSigs._standardize_trinuc(row['TC'])
    return (trinuc[0] + '[' + subs + ']' + trinuc[2])


def plot_trinucleotide(snv_cn, somatic, fasta_path, sample, tmp="tmp", cn_frac=0):
    data = somatic.merge(snv_cn[['chr', 'pos', 'frac_cn']])
    data['mutation_class'] = data.apply(lambda row: calculate_mutation_class(row), axis=1)

    # context_counts = data.groupby('mutation_class').size().to_dict()
    context_counts = data[data.frac_cn > cn_frac].groupby('mutation_class').size().to_dict()

    ds = DeconstructSigs(
        maf=sample,
        context_counts=context_counts,
        hg19_fasta_path=fasta_path
    )

    dslogger = logging.getLogger("DeconstructSigs")
    dslogger.setLevel(logging.ERROR)

    weights = ds.which_signatures()
    # return ds, weights
    return ds, weights
