import pandas as pd
from matplotlib.lines import Line2D
import numpy as np

# ~~~~~~~~~~~~~~~~~~~~ code taken from Andrew McPherson ~~~~~~~~~~~~~~~~~~~~~~ #

def read(remixt_file, sample_label):
    """
    read in remixt data into pandas dataframe
    :param remixt_file: remixt input file (i.e. results_files_*.h5)
    :param sample_label: label of sample inside h5
    :return: parsed pandas dataframe
    """
    cnv_data = list()
    with pd.HDFStore(remixt_file, 'r') as store:
        prefix = '/copy_number/{pred_tool}/sample_{sample_id}'.format(sample_id=sample_label,
                                                                      pred_tool=sample_label)
        mixture = store['/mix']
        cn = store['/cn']
        # Ensure that major and minor are correct, otherwise swap
        swap = (cn['major_1'] < cn['minor_1']) & (cn['major_2'] < cn['minor_2'])
        cn.loc[swap, ['major_1', 'minor_1']] = cn.loc[swap, ['minor_1', 'major_1']].values
        cn.loc[swap, ['major_2', 'minor_2']] = cn.loc[swap, ['minor_2', 'major_2']].values
        cn.loc[swap, 'major_is_allele_a'] = 1 - cn.loc[swap, 'major_is_allele_a']
        cn['total_1'] = cn['minor_1'] + cn['major_1']
        cn['total_2'] = cn['minor_2'] + cn['major_2']
        if mixture[1] > mixture[2]:
            cn['major'] = cn['major_1']
            cn['minor'] = cn['minor_1']
            cn['total'] = cn['total_1']
        else:
            cn['major'] = cn['major_2']
            cn['minor'] = cn['minor_2']
            cn['total'] = cn['total_2']
        swap = (cn['major'] < cn['minor'])
        cn.loc[swap, ['major', 'minor']] = cn.loc[swap, ['minor', 'major']].values
        assert (cn['major'] >= cn['minor']).all()
        cn['is_subclonal'] = (
                                     (cn['minor_1'] != cn['minor_2']) |
                                     (cn['major_1'] != cn['major_2'])) * 1
        cn['tumour_content'] = mixture[1:].sum()
        cnv_data.append(cn)
    cnv_data = pd.concat(cnv_data, ignore_index=True)
    cnv_data["total_raw_e"] = cnv_data.major_raw_e + cnv_data.minor_raw_e
    return cnv_data