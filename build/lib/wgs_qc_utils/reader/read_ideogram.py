import pandas as pd


def read(ideo_file):
    # read things in
    ideogram = pd.read_csv(ideo_file, sep="\t", names=["chrom", "start", "end", "name", "stain"])
    # Colors for different chromosome stains
    ideogram = ideogram.astype({"chrom":str})
#    ideogram["chrom"] = ideogram.chrom.str.split("chr", expand=True)[1]
    ideogram["start"] = ideogram.start/1000000
    ideogram["end"] = ideogram.end/1000000

    color_lookup = {
        'gneg': [1., 1., 1.],
        'gpos25': [.6, .6, .6],
        'gpos50': [.4, .4, .4],
        'gpos75': [.2, .2, .2],
        'gpos100': [0., 0., 0.],
        'acen': [.8, .4, .4],
        'gvar': [.8, .8, .8],
        'stalk': [.9, .9, .9]
    }

    ideogram["color"] = ideogram.stain.apply(lambda x: color_lookup[x])
    ideogram['width'] = ideogram.end - ideogram.start

    return ideogram


def prepare_at_chrom(ideogram, chromosome):
    return ideogram[ideogram.chrom == chromosome]
