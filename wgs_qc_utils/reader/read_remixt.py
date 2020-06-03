import pandas as pd

# ~~~~~~~~~~~~~~~~~~~~ code from Andrew McPherson ~~~~~~~~~~~~~~~~~~~~~~ #


def read(remixt_file, sample_label):
    """
    read in remixt data into pandas dataframe
    :param remixt_file: remixt input file (i.e. results_files_*.h5)
    :param sample_label: label of sample inside h5
    :return: parsed pandas dataframe
    """
    if not remixt_file:
        return None
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


def prepare_at_chrom(parsed_remixt, chrom):
    """
    prep copy number rdata for plotting at a chrom
    :param remixt: parsed remixt pandas dataframe
    :param chrom: chromosome
    :return: remixt at chromosome
    """
    if not parsed_remixt:
        return None
    return parsed_remixt[parsed_remixt["chromosome"] == chrom]


def make_for_circos(remixt, sample_id, prepped_remixt):
    '''
    parse in the remixt h5 file and then write it out to a csv for use in the circos plot
    :param remixt: remixt h5 filename
    :param sample_id: sample id
    :param prepped_remixt: output filename
    :return:
    '''
    if not remixt:
        return None
    remixt = read(remixt, sample_id)
    remixt.to_csv(prepped_remixt, sep="\t", index=False, header=True)