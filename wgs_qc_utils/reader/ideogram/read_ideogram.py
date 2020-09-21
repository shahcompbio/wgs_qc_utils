import pandas as pd
from wgs_qc_utils.reader.ideogram import read_ideogram
import os
import pkg_resources


def read(ideogram_file=False):

    if ideogram_file == False:
        stream = pkg_resources.resource_stream(__name__, 'ideogram.txt')
        ideogram = pd.read_csv(stream, sep="\t", encoding='latin-1', names=["chrom", "start", "end", "name", "stain"])

    else:
        ideogram = pd.read_csv(ideogram_file, sep="\t", names=["chrom", "start", "end", "name", "stain"])

    ideogram = ideogram.astype({"chrom":str})
    ideogram["start"] = ideogram.start/1000000
    ideogram["end"] = ideogram.end/1000000
    ideogram["chrom"] = ideogram.chrom.str.lower()

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
    print(ideogram)
    ideogram["color"] = ideogram.stain.apply(lambda x: color_lookup[x])
    ideogram['width'] = ideogram.end - ideogram.start

    return ideogram


def prepare_at_chrom(ideogram, chromosome):
    return ideogram[ideogram.chrom == chromosome]
