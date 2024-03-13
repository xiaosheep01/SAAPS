# -*- coding: utf-8 -*-
# Development Time: 2024-02-27 10:54:00
# Developer: XiaoYang

import linecache
import os
import platform
import re
import sys
from functools import reduce
from itertools import cycle
from typing import Any
import constant as ct
from matplotlib import pyplot as plt, patches
import matplotlib.colors as mcolors
import seaborn as sb
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import OneHotEncoder
import numpy as np
import pandas as pd
from colorama import Fore


def files_path(dirname, sys_os):
    """
    Gets all file paths of the folder
    :param sys_os: system operation type
    :param dirname: folder path
    :return: files path list
    """
    file_path_list = []
    dir_path_list = []
    for root, dirs, files in os.walk(dirname):
        # root represents the path to the folder being accessed
        # dirs represents the names of all subfolders in this folder, and the list is stored
        # files said evil all child filename, under the folder list is stored
        for item in files:
            file_path_list.append(os.path.join(root, item))

        for item in dirs:
            dir_path_list.append(os.path.join(root, item))

    # Judgment system type
    # MacOS
    if sys_os == "Darwin":
        print("{}Note:Current operating system: {}, remove hidden file path by "
              "default.{}".format(Fore.YELLOW, platform.uname()[0], Fore.RESET))
        for item in file_path_list:
            item_len = len(item.split(os.sep))
            if item.split(os.sep)[item_len - 1].startswith("."):
                file_path_list.remove(item)

        # Force out.ds_store files
        file_path_list.remove(os.path.join(dirname, ".DS_Store"))

    # Windows
    elif sys_os == "Windows":
        pass

    # Linux
    elif sys_os == "Linux":
        pass

    return file_path_list


def file_type_judge(file_path):
    file_extension = file_path.split('.')[-1].lower()
    # determine file type
    if file_extension == 'csv':
        return 'CSV'
    elif file_extension == 'xlsx':
        return 'XLSX'
    elif file_extension == "txt":
        return "TXT"
    else:
        return 'Unknown'


def files_combine(file_path_list):
    content = []

    for item in file_path_list:
        for line in linecache.getlines(item):
            content.append(line)

    return content


def read_file(file_path):
    content = []
    for line in linecache.getlines(file_path):
        content.append(line)
    return content


def fasta_doc_std(file_list):
    """
    Standardization of fasta documents
    NOTE: the result (list) format is that seq name and seq locate in same line!!!
    :param file_list:file list by function:files_combine
    :return:Standard fasta content
    """
    seq_name_list = []
    seq_dic = {}
    seq_list = []
    seq_name = ""

    for line in file_list:
        line = line.strip()
        if line.startswith(">"):
            seq_name = line
            seq_name_list.append(seq_name)
            seq_dic[seq_name] = []
        else:
            seq_dic[seq_name].append(line)

    for seq_name in seq_name_list:
        seq_list.append(seq_name + "\n" + "".join(seq_dic[seq_name]))           # element: seq_name + "\n" + seq

    return seq_list  # return list


def seq_cut(start_info, end_info, seq_file_list, cut_tag):
    """
    according to the input seq, including start seq and end seq, this is function will locate them in the reference seq,
    and obtain the index (location), then cut the sequences to target fragment.
    :param end_info: the end fragment of reference sequence
    :param start_info: the start fragment of reference sequence
    :param cut_tag: the input data type: str or int
    :param seq_file_list: the aligned sequence file
    :return: target sequence fragment
    """
    result_list = []
    ref_seq = seq_file_list[0].split("\n")[1]

    if cut_tag == "str":
        ref_seq = ref_seq.upper()  # reference seq and convert to upper
        start_res = re.search(start_info, ref_seq)
        end_res = re.search(end_info, ref_seq)
        for line in seq_file_list:
            seq_name, seq = line.split("\n")
            match_res = seq[start_res.start(): end_res.end()].upper()
            result_list.append(seq_name + "\n" + match_res)

    elif cut_tag == "int":
        start_info = start_info - 1
        for line in seq_file_list:
            seq_name, seq = line.split("\n")
            match_res = seq[start_info: end_info]
            result_list.append(seq_name + "\n" + match_res)

    return result_list  # return a list


def list_to_save(file_in_list, file_path):
    with open(file_path, "w") as F:
        for item in file_in_list:
            F.write(item + "\n")


def nucleotide_judge(seq):
    """
    judge the nucleotide type (DNA or RNA)
    :param seq: the seq format is 'string'
    :return: nucleotide type
    """
    seq = seq.upper()
    if seq.find("T") != -1 and seq.find("U") == -1:
        return "DNA"
    elif seq.find("U") != -1 and seq.find("T") == -1:
        return "RNA"
    elif seq.find("T") == -1 and seq.find("U") == -1:
        return "error"


def dna_to_rna(seq):
    seq = seq.upper()
    seq = re.sub("T", "U", seq)
    return seq


def rna_to_dna(seq):
    seq = seq.upper()
    seq = re.sub("U", "T", seq)
    return seq.upper()


def codon_check(seq_list):
    seq_len_list = []
    # judge the type of nucleic acid
    ref_seq = seq_list[0].split("\n")[1]            # the first seq
    nt_tag = nucleotide_judge(ref_seq)
    print("NOTE: Input sequence is %s" % nt_tag)

    for line in seq_list:
        seq_name, seq = line.split("\n")
        seq_length = len(seq)
        seq_len_list.append(seq_length)

    if len(set(seq_len_list)) != 1:
        print(Fore.RED + "ERROR: Sequence length is inconsistent, please check the sequences!")
        sys.exit()

    if seq_len_list[0] % 3 != 0:
        print(Fore.RED + "ERROR: The sequence length does not conform to the law of codon translation, "
                         "please check the sequences!")
        sys.exit()

    return nt_tag


def translate_to_protein(seq_list, nt_tag):  # DNA in default
    result_list = []
    if nt_tag == "error":
        print(Fore.RED + "ERROR: Sequence information not recognized, please check sequence!")
        sys.exit()

    for line in seq_list:
        seq_name, seq = line.split("\n")

        if nt_tag == "RNA":
            seq = rna_to_dna(seq)

        result_list.append(seq_name)
        seq_length = len(seq)
        protein = ""
        for n in range(0, seq_length, 3):  # three nucleotides as a codon
            codon = seq[n: n + 3]
            try:
                amino_acid = ct.translate_table[codon]
            except:
                amino_acid = "?"
            finally:
                pass
            protein = protein + amino_acid
        result_list.append(protein)
    return result_list


def delete_gap(file_list):
    seq_len_list = []
    seq_dict = {}

    for item in file_list:
        seq_name, seq = item.split("\n")
        seq_dict[seq_name] = list(seq)
        seq_len_list.append(len(seq))

    if len(set(seq_len_list)) != 1:
        print(Fore.RED + "ERROR: Sequence length is inconsistent, please check the sequences!")
        sys.exit()

    seq_df = pd.DataFrame(seq_dict).transpose()

    gaps_list = []
    print("NOTE: Site order starts at 0, not 1!")
    for column in seq_df.columns:
        if any("-" in value for value in seq_df[column]):
            print(f"NOTE: This column contains gap:{column}")
            gaps_list.append(column)

    print("\n" + "NOTE: The columns containing gaps have been deleted, and the processed data is following (partial):")
    seq_df.drop(columns=gaps_list, inplace=True)
    print(seq_df.head())

    def combine_columns(row):
        return "".join(row.values)

    melt_seq_df = seq_df.apply(lambda row: combine_columns(row), axis=1)

    result_list = []
    for index, value in melt_seq_df.iteritems():
        result_list.append(index + "\n" + value)

    return result_list


def compute_polymorphism(seq_list):
    """
    this function will compute the polymorphism of the amino acids, and generate a table on each polymorphism site.
    It contains several steps as following.
    Step 1: convert fasta file to a dataframe.
    Step 2: compute the frequency and density of amino acid in every site.
    Step 3: record each variable site in every sequence.
    :param seq_list:target sequence fragment
    :return:a list
    """
    seq_dict = {}
    # seq_list.remove(seq_list[0])  # delete the reference sequence
    for item in seq_list:
        seq_name, seq = item.split("\n")
        seq = list(seq)  # convert string to list
        seq_dict[seq_name] = seq  # convert to a dictionary

    # in this dataframe, row represents amino acids in each site, column represents sequence names
    seq_df = pd.DataFrame(seq_dict)  # convert to a dataframe

    # count the frequency of each amino acid in each row
    def count_row_freq(row):
        value_list = row.to_list()  # value_list: [G, A, C, T, A...]
        count_dict = {}

        for value in value_list:
            if value not in count_dict:
                count_dict[value] = 1
            else:
                count_dict[value] += 1

        freq_dict = {key: "{:.2f}%".format(value / len(value_list) * 100) for key, value in count_dict.items()}

        # the format of concise result looks like this: [(A, 4, 20%)(B, 5, 25%)]
        count_freq_list = []
        for my_string in set(value_list):
            count_freq_list.append("({0}, {1}, {2})".format(
                my_string, str(count_dict[my_string]), str(freq_dict[my_string])
            ))

        # print(inner_detail_list)
        return count_freq_list

    def count_aa_name(row):
        value_list = row.to_list()
        count_dict = dict()
        for value in value_list:
            if value not in count_dict:
                count_dict[value] = 1
            else:
                count_dict[value] += 1

        freq_dict = {key: "{:.2f}%".format(value / len(value_list) * 100) for key, value in count_dict.items()}

        # detail result
        # the format of detail result looks like this: [A:>seq 1,>seq2,>seq3...,B: >seq1, >seq2..]
        # the detail should be represented according to the threshold
        inner_detail_list = []  # detail_list: [(A: >seq_1, >seq_2...)...]
        name_dict = row.to_dict()  # name_dict: {">seq_1": "A",...}
        threshold_value = 70  # if less than threshold, represent sequence names

        for name_key, name_value in name_dict.items():
            seq_name_info = ""  # sequence name information
            aa_freq = float(freq_dict[name_value].split("%")[0])
            aa_freq = int(aa_freq)
            if aa_freq == 100:
                seq_name_info = "All Sequences"

            elif threshold_value < aa_freq < 100:
                seq_name_info = "Most Sequences"

            elif aa_freq <= threshold_value:
                seq_name_info = name_key

            inner_detail_list.append("{}:{}".format(name_value, seq_name_info))  # contain duplicated data
        inner_detail_list = list(set(inner_detail_list))  # remove duplicated data

        # merge the seq name of the same amino acid
        aa_name_dict: dict[Any, Any] = dict()  # aa_name_dict:{A: >seq_1,>seq_2..., G: >seq_3..}
        for inner_item in inner_detail_list:
            aa, name = inner_item.split(":")

            if aa not in aa_name_dict:
                aa_name_dict[aa] = name
            else:
                aa_name_dict[aa] = aa_name_dict[aa] + "," + name

        final_result = list()  # convert the dictionary to a list
        for key, value in aa_name_dict.items():
            final_result.append("({}: {})".format(key, value))  # modify the format of the list

        return final_result

    # a list contain amino acid, count and frequency in each row
    count_freq_list = seq_df.apply(count_row_freq, axis=1)
    detail_list = seq_df.apply(count_aa_name, axis=1)

    # 'concise result' is a list, separated by '\t'
    concise_result = []
    for n in range(len(count_freq_list)):
        site_info = ""
        for item in count_freq_list[n]:
            site_info = site_info + item + "\t"
        concise_result.append("Site {}{}{}".format(str(n + 1), "\t", site_info))

    # 'detailed_result' is a list, separated by ','
    detailed_result = []
    for n in range(len(detail_list)):
        line_info = ""
        for item in detail_list[n]:
            line_info = line_info + item + "\t"
        detailed_result.append("Site {}{}{}".format(str(n + 1), "\t", line_info))

    return seq_df, concise_result, detailed_result


def change_name(seq_list, name_file_path):
    """
    change the name in sequence files
    :param seq_list: sequence file
    :param name_file_path: old and new name list, separated by tab, it could be a csv, xlsx or txt file
    :return: new seq_list with new name
    """
    # judge the name_file_path type
    name_dict = dict()
    if name_file_path.endswith(".txt"):
        contents = linecache.getlines(name_file_path)
        contents = [x.strip() for x in contents if x.strip() != ""]
        for item in contents:
            old_name, new_name = item.split("\t")
            name_dict[old_name] = new_name

    elif name_file_path.endswith(".csv"):
        name_df = pd.read_csv(name_file_path, sep=",")
        for n in range(len(name_df)):
            name_dict[name_df.iloc[n, 0]] = name_df.iloc[n, 1]

    elif name_file_path.endswith(".xlsx"):
        name_df = pd.read_excel(name_file_path)
        for n in range(len(name_df)):
            name_dict[name_df.iloc[n, 0]] = name_df.iloc[n, 1]

    else:
        print(Fore.YELLOW + "Warning: the file format is not supported right now, proceeding in ‘txt’ format!")
        contents = linecache.getlines(name_file_path)
        contents = [x.strip() for x in contents if x.strip() != ""]
        for item in contents:
            old_name, new_name = item.split("\t")
            name_dict[old_name] = new_name

    # replace the old name with new name
    result_list = list()
    for item in seq_list:
        seq_name = item.split("\n")[0]
        new_name = name_dict[seq_name]

        if not new_name.startswith(">"):
            new_name = ">" + new_name

        new_item = item.replace(seq_name, new_name)
        result_list.append(new_item)

    return result_list


def output_name(seq_list):
    result_list = []
    for line in seq_list:
        seq_name = line.split("\n")[0]
        result_list.append(seq_name)

    result_df = pd.DataFrame({"Old_Name": result_list})
    print("NOTE: Seq Name is following (Partial):")
    print(result_df.head())
    return result_df


def shannon_h_value(list_tuple, tag=20):
    """
    calculate shannon IC value
    :param list_tuple: format like this: [(A,80),(B,90)...]
    :param tag: tag has three conditions: 20, 21 and 22
    :return: IC value in list
    """
    h_value = 0
    aa_list = []
    for tuple_item in list_tuple:
        aa_list.append(tuple_item[0])  # restore amino acid into a list
        prob = tuple_item[1] / 100
        h_value = h_value + (-prob * np.log2(prob))

    IC = np.log2(tag) - h_value
    IC = round(IC, 5)
    return [IC]


def shannon_entropy(file_path, maximum=0, minimum=0):
    """
    compute the shannon entropy and filter the result according to the threshold
    :param minimum: the minimum of the threshold
    :param maximum: the maximum of the threshold
    :param file_path: the concise result path, it should be a txt file
    :return: a result restricted in a certain range
    """
    file_content = linecache.getlines(file_path)
    file_content = [x.strip() for x in file_content if x.strip() != ""]
    file_string = "".join(file_content)
    # determine whether there is illegal icons such as '?' or '-' in sequences
    tag = 20
    if "?" in file_string and "-" in file_string:
        print("Sequences contain illegal icon: '?' and '-', and IC tag is 22.")
        tag = 22
    elif "?" in file_string:
        print("Sequences contain illegal icon: '?' or '-', and IC tag is 21.")
        tag = 21
    elif "-" in file_string:
        print("Sequences contain illegal icon: '?' or '-', and IC tag is 21.")
        tag = 21
    elif "?" not in file_string and "-" not in file_string:
        print("Sequences will be processed in default way, and IC tag is 20.")
        tag = 20

    peptide_dict = dict()
    for item in file_content:
        site, residue = item.split("\t", 1)  # 'site' is the residue location
        residue = residue.split("\t")  # 'residue' is the amino acid
        res_list = []

        for values in residue:
            values = values.replace("(", "").replace(")", "").replace(" ", "")
            values = values.split(",")

            first_letter = values[0]
            prob = float(values[2].split("%")[0])
            my_tuple = (first_letter, prob)

            res_list.append(my_tuple)

        infor_content = shannon_h_value(res_list, tag)
        peptide_dict[site] = infor_content

    peptide_df = pd.DataFrame(peptide_dict)
    peptide_df = peptide_df.transpose()
    peptide_df.reset_index(inplace=True)
    peptide_df.columns = ["SITE", "IC"]
    statistics_info = peptide_df["IC"].describe()
    shannon_log_path = os.path.dirname(file_path) + os.sep + "Shannon_log.txt"

    with open(shannon_log_path, "w") as F:
        F.write("Basic Statistic Information of the Shannon Entropy (IC)" + "\n")

    statistics_info.to_csv(shannon_log_path, sep="\t", mode="a")

    # obtain the result in assigned range
    if maximum != 0 and minimum != 0:
        peptide_df = peptide_df[(peptide_df["IC"] >= minimum) & (peptide_df["IC"] <= maximum)]
        filtered_df = peptide_df.copy()
        filtered_df["Pure_Site"] = filtered_df["SITE"].str.extract(r'(\d+)')
        pure_site = filtered_df["Pure_Site"].tolist()

        with open(shannon_log_path, "a") as F:
            F.write("\n" + "\n")
            F.write("The following is the key sites filtered by assigned minimum (%s) and the maximum (%s) :"
                    % (str(minimum), str(maximum)))
            F.write("\n")

            for i, value in enumerate(pure_site, 1):
                F.write(str(value))

                if i % 10 == 0:
                    F.write("\n")
                else:
                    F.write(", ")

            F.write("\n" + "\n" + "The following is the detailed information of the sites:" + "\n")
            filtered_list = filtered_df.values.tolist()
            for row in filtered_list:
                F.write(", ".join(map(str, row)) + "\n")

    return peptide_df


def shannon_filter(file_df, maximum=99999, minimum=-99999):
    filtered_df = file_df[(file_df["IC"] >= minimum) & (file_df["IC"] <= maximum)]
    filtered_df["Pure_Site"] = filtered_df["SITE"].str.extract(r'(\d+)')
    pure_site = filtered_df["Pure_Site"].tolist()

    return pure_site  # This is a list containing key sites


def one_hot_encoding(seq_matrix, site_list):
    """
    extract key sites of the sequence matrix
    :param seq_matrix: the sequence matrix which row and column stand for sequence and residue site, respectively
    :param site_list: a list contains all putative sites
    :return: a dataframe
    """
    filtered_matrix = seq_matrix.iloc[site_list]
    # print(filtered_matrix.head())

    tag = "type1"                                                                               # type1 is default
    if "?" in filtered_matrix.values and "-" in filtered_matrix.values:                         # both "?" and "-"
        tag = "type4"
        print(f"NOTE: Sequences contain illegal icon, IC tag is '{tag}'.")
    elif "?" not in filtered_matrix.values and "-" in filtered_matrix.values:                   # not "?" but "-"
        tag = "type3"
        print(f"NOTE: Sequences contain illegal icon, IC tag is '{tag}'.")
    elif "?" in filtered_matrix.values and "-" not in filtered_matrix.values:                   # not "-" but "?"
        tag = "type2"
        print(f"NOTE: Sequences contain illegal icon, IC tag is '{tag}'.")
    elif "?" not in filtered_matrix.values and "-" not in filtered_matrix.values:               # not '?' or '-'
        tag = "type1"
        print(f"NOTE: Sequences contain no illegal icon, IC tag is '{tag}'.")

    raw_aa_mat = filtered_matrix.transpose()
    raw_aa_mat.columns = [int(col) + 1 for col in raw_aa_mat.columns.tolist()]

    # amino acids array: 1*n  (according to the tag)
    amino_acids_array = np.array(ct.amino_acids_list[tag]).reshape(-1, 1)
    # create a one-hot object
    encoder = OneHotEncoder(sparse=False, dtype=int)
    # one-hot encoding
    onehot_encoded = encoder.fit_transform(amino_acids_array)
    amino_acid_onehot_dict = dict(zip(ct.amino_acids_list[tag], onehot_encoded.tolist()))

    # replace the values according to the amino acid one-hot dictionary
    def transform_element(value):
        return amino_acid_onehot_dict[value]

    filtered_matrix = filtered_matrix.applymap(transform_element)
    filtered_matrix = filtered_matrix.transpose()

    old_columns = filtered_matrix.columns.tolist()
    new_columns = [int(col) + 1 for col in old_columns]
    filtered_matrix.columns = new_columns

    # expand matrix
    # different amino acids according to the 'tag'
    # expand the dimension of each column
    expand_matrix = pd.DataFrame()
    for col in filtered_matrix.columns.tolist():
        col_list = filtered_matrix[col].tolist()
        # The range will automatically change according to the 'tag'
        temp_df = pd.DataFrame(col_list,
                               columns=[str(col) + "_%s" % str(i + 1) for i in range(ct.tag_amino_list[tag])])
        expand_matrix = pd.concat([expand_matrix, temp_df], axis=1)

    expand_matrix.index = filtered_matrix.index

    # site matrix
    site_matrix = pd.DataFrame()
    for column in raw_aa_mat.columns:
        site_matrix[column] = raw_aa_mat[column].apply(lambda x: f'{column}{x}')

    return raw_aa_mat, filtered_matrix, expand_matrix, site_matrix


def pca(onehot_file_path):
    """
    Dimension reduction by principal component analysis
    :param onehot_file_path: the 'OneHot-Transform.csv'
    :return: dimension reduction result in csv file
    """
    filetype = file_type_judge(onehot_file_path)
    if filetype == "CSV":
        onehot_mat = pd.read_csv(onehot_file_path, header=0, index_col=0)
    elif filetype == "XLSX":
        onehot_mat = pd.read_excel(onehot_file_path, header=0, index_col=0)
    else:
        print(Fore.RED + "ERROR! The file type:%s is not supported right now!" % filetype)
        sys.exit()

    # pca function
    dimension = 1  # the minimal dimension
    cumulative_variance_ratio = 0       # the cumulative variance ratio
    while True:
        print(f"Trying to reduct dimensions to {dimension}")
        pca_value = PCA(n_components=dimension)
        new_data = pca_value.fit_transform(onehot_mat)
        # 'explained variance ratio' contains variance ratio of each principal
        variance_ratios = pca_value.explained_variance_ratio_
        cumulative_variance_ratio = np.cumsum(variance_ratios)[-1]

        # if the cumulative variance ratio equals 1, break circle
        if cumulative_variance_ratio >= 1:
            break

        # if the dimensions more than the amount of rows, break circle
        dimension += 1
        if dimension - 1 == onehot_mat.shape[0]:
            break

    print(f"The cumulative variance is {Fore.BLUE}{cumulative_variance_ratio}{Fore.RESET}")

    new_data = pd.DataFrame(new_data)
    new_data.index = onehot_mat.index
    print()
    print("The PCA scores: ", new_data.head(), sep="\n")
    return new_data


def tsne(onehot_file_path, perp, learn_rate, niter, rand_state):
    filetype = file_type_judge(onehot_file_path)
    if filetype == "CSV":
        onehot_mat = pd.read_csv(onehot_file_path, header=0, index_col=0)
    elif filetype == "XLSX":
        onehot_mat = pd.read_excel(onehot_file_path, header=0, index_col=0)
    else:
        print(Fore.RED + "ERROR! The file type:%s is not supported right now!" % filetype)
        sys.exit()

    t_sne = TSNE(n_components=2,
                 perplexity=perp,
                 learning_rate=learn_rate,
                 n_iter=niter,
                 random_state=rand_state)

    tsne_data = t_sne.fit_transform(onehot_mat)
    tsne_data = pd.DataFrame(tsne_data)
    tsne_data.index = onehot_mat.index
    print("NOTE: The t_SNE scores (Partial):", tsne_data.head(), sep="\n")

    return tsne_data


def set_differences(exp_mat, set_info):

    # remove all the chars started with '>' in the index column
    exp_mat.index = exp_mat.index.str.lstrip('>')

    print("Input Polymorphism Matrix (Partial): ")
    print(exp_mat.head())
    print("\n" + "Input Sets Annotation Table (Partial): ")
    print(set_info.head())
    groups = set_info["Group"].unique()             # get all groups information

    group_dict = dict()
    commom_res_dict = dict()
    for item in groups:
        match_df = set_info[set_info["Group"] == item]
        match_species = match_df["Species"].tolist()
        group_dict[item] = match_species

        group_df = exp_mat.loc[group_dict[item]]
        sets_list = [set(row) for index, row in group_df.iterrows()]
        common_elements = reduce(set.intersection, sets_list)
        commom_res_dict[item] = common_elements

    value_list = list()
    print("The intersection result of each group is following:")
    for key, value in commom_res_dict.items():
        value_list.append(value)
        print(key, value)

    def multiple_sets_difference(sets_list):
        result_list = list()
        sets_number = len(sets_list)

        for i in range(sets_number):
            target_set = sets_list[i].copy()
            temp_sets = sets_list.copy()
            temp_sets.pop(i)

            for my_set in temp_sets:
                target_difference = target_set.difference(my_set)
                target_difference = sorted(list(target_difference))

            result_list.append(target_difference)

        return result_list

    def multiple_sets_intersection(sets_list):
        result = reduce(set.intersection, sets_list)
        return result

    difference_list = multiple_sets_difference(value_list)
    intersection_result = multiple_sets_intersection(value_list)
    print("\n" + "The Intersection Result of all Groups is Following:")
    print(intersection_result)

    # differ_dict = dict(zip(groups, difference_list))
    differ_list = list()
    for n in range(len(groups)):
        differ_value = [groups[n], group_dict[groups[n]], difference_list[n]]
        differ_list.append(differ_value)

    print("\n" + "The Difference Result of all Groups is Following:")

    print(differ_list)
    differ_df = pd.DataFrame(differ_list, columns=["Group", "Species", "Sites"])
    return differ_df


def generate_random_colors(num_colors):
    colors = []

    for i in range(num_colors*2):
        # create random colors
        rgb = np.random.rand(1, 3)
        # convert format
        color = mcolors.to_hex(rgb)
        colors.append(color)

    colors = list(set(colors))
    np.random.shuffle(colors)
    random_colors = colors[:num_colors]
    return random_colors


def polymorphism_figure(compute_polymorphism_df,
                        compute_polymorphism_concise_list,
                        figure_name,
                        figure_dpi,
                        transparent_tag,
                        fig_width=14,
                        fig_height=9):

    original_df = compute_polymorphism_df
    original_list = compute_polymorphism_concise_list

    polymorphism_site_list = list()                 # all polymorphism sites
    conserve_site_list = list()                     # all conserve sites

    for item in original_list:
        new_item = item.split("\t")
        new_item = [x.strip() for x in new_item if x.strip() != ""]
        new_item_len = len(new_item)

        var_site = new_item[0].split(" ")[1]
        # conserve sites
        if new_item_len == 2:
            conserve_site_list.append(var_site)
        # polymorphism sites
        elif new_item_len > 2:
            polymorphism_site_list.append(var_site)

    seq_names_list = original_df.columns.tolist()

    # PLOTTING PART
    raw_figure_width = len(polymorphism_site_list)                      # polymorphism sites number
    raw_figure_height = len(seq_names_list)                             # sequences number

    figure_width = fig_width
    figure_height = fig_height
    print("NOTE: The figure width is:%s" % str(figure_width))
    print("NOTE: The figure height is:%s" % str(figure_height))

    fig, ax = plt.subplots(nrows=1,
                           ncols=1,
                           figsize=(figure_width, figure_height),
                           dpi=100)

    x_axis_len = range(1, len(polymorphism_site_list) + 1)
    y_axis_len = range(1, len(seq_names_list) + 1)
    x_axis_label = polymorphism_site_list                               # x-axis label = polymorphism sites
    y_axis_label = seq_names_list                                       # y-axis label = sequence name
    y_axis_label = [x.replace(">", "") for x in y_axis_label]           # delete '>' symbol in the sequence name

    # plot x-axis limit: 0 ~ number of polymorphism sites + 0.5
    plt.xlim(0, raw_figure_width + 1)
    # plot y-axis limit: 0.3 ~ number of sequences + 0.5
    plt.ylim(0.5, raw_figure_height + 0.5)
    plt.xticks(x_axis_len, x_axis_label, rotation=60, size=8)           # plot x-axis label
    plt.yticks(y_axis_len, y_axis_label, size=8)                        # plot y-axis label

    # conceal axis line
    ax.tick_params(bottom=True, left=False)
    # conceal axis border
    for location in ["top", "left", "bottom", "right"]:
        ax.spines[location].set_visible(False)

    # record polymorphism sites
    polymorphism_df = original_df.transpose().copy()
    var_sites_list = [int(x) - 1 for x in polymorphism_site_list]
    plot_df = polymorphism_df.filter(items=var_sites_list)

    rect_width, rect_height = 0.8, 1                                # tiny rectangle width and height
    rect_x_loc, rect_y_loc = 0.6, 0.5                               # tiny rectangle initial location
    text_width, text_height = 1, 1                                  # amino acid text gap
    text_x_loc, text_y_loc = 1, 1                                   # amino acid initial location
    bg_rect_x, bg_rect_y = 0.2, 0.5                                 # background rectangle x and y  initial location
    bg_rect_width, bg_rect_height = raw_figure_width + 0.8, 1       # background rectangle width and height

    # plot background rectangle
    bg_color_list = cycle(["#E8E8E8", "#FFFFFF"])

    def bg_color():
        return next(bg_color_list)

    for n in seq_names_list:

        ax.add_patch(patches.Rectangle((bg_rect_x, bg_rect_y),
                                       bg_rect_width,
                                       bg_rect_height,
                                       facecolor=bg_color(),
                                       fill=True,
                                       alpha=0.6,
                                       edgecolor="none"))

        bg_rect_y = bg_rect_y + bg_rect_height

    # cycle each row
    for row_index, row_value in plot_df.iterrows():
        plot_text_list = row_value.tolist()                         # format:["A", "T", "M", "P"...]

        for aa in plot_text_list:

            ax.add_patch(patches.Rectangle((rect_x_loc, rect_y_loc),
                                           rect_width,
                                           rect_height,
                                           facecolor=ct.aa_paltte[aa],
                                           fill=True,
                                           alpha=0.6))
            ax.text(text_x_loc, text_y_loc, aa, fontsize=9, ha="center", va="center")
            # x-location changing, y location keeps same
            # '0.2' stands for the gap between two tiny rectangles
            rect_x_loc = rect_x_loc + rect_width + 0.2
            text_x_loc = text_x_loc + text_width

        rect_x_loc, text_x_loc = 0.6, 1
        rect_y_loc = rect_y_loc + rect_height
        text_y_loc = text_y_loc + text_height

    plt.xlabel("Residue", fontsize=12)
    plt.tight_layout()
    # save figure
    plt.savefig(fname=figure_name, dpi=figure_dpi, transparent=transparent_tag)


def ic_scatter_figure(shannon_ic_df,
                      fig_name,
                      fig_dpi,
                      tp_tag,
                      color_tag,
                      fig_width=10,
                      fig_height=8,
                      ):

    # x-axis & y-axis information
    x_info = shannon_ic_df.iloc[:, 0]
    y_info = shannon_ic_df.iloc[:, 1]
    x_info = x_info.str.replace("Site ", "")
    x_info = x_info.astype(int)
    y_info = y_info.round(2)

    # create new x-ticks
    x_max = max(x_info)
    x_min = min(x_info)
    new_x_ticks = np.linspace(x_min, x_max, 7)
    new_x_ticks = new_x_ticks.astype(int)

    # figure size
    plt.figure(figsize=(fig_width, fig_height))

    # line plot
    plt.plot(x_info, y_info, linewidth=0.9, zorder=1, color="grey")

    # palette
    colors = color_tag
    if color_tag.upper() == "RANDOM":
        # colors number
        print(Fore.YELLOW + "WARNING: IC Scatter plots currently do not support 'random' palette!")
        print(Fore.YELLOW + "WARNING: Using default option: 'RdYlBu' palette!")
        colors = "RdYlBu"

    elif color_tag.upper() == "RANDOM_2":
        random_palette = np.random.choice(ct.palettes)

        if re.search("/", random_palette):
            colors_list = random_palette.split("/")
            random_palette = np.random.choice(colors_list)

        colors = "".join(random_palette)

    print(f"NOTE: the palette is {colors}")
    # scatter plot
    plt.scatter(x=x_info, y=y_info, c=y_info, cmap=colors, zorder=2)

    # plot xy-ticks and xy-labels
    plt.colorbar(label="IC Value", shrink=0.5, aspect=18)
    plt.xticks(ticks=new_x_ticks, rotation=60, size=8)
    plt.xlabel("Residue", size=12)
    plt.ylabel("IC Value", size=12)

    plt.tight_layout()
    # save figure
    plt.savefig(fname=fig_name, dpi=fig_dpi, transparent=tp_tag)


def cluster_figure(dr_result,
                   point_label,
                   color_tag,
                   fig_name,
                   fig_dpi,
                   tp_tag,
                   fig_width,
                   fig_height):

    pc_1 = dr_result.iloc[:, 0]
    pc_2 = dr_result.iloc[:, 1]
    virus_info = dr_result.iloc[:, -1]

    # palette
    colors = color_tag
    if color_tag.upper() == "RANDOM":
        # colors number
        num_colors = len(virus_info)
        # create colors
        colors = generate_random_colors(num_colors)
        # print colors
        print("\n" + "Generated Random Colors:")
        for n in range(len(colors)):
            print(f"{virus_info[n]}:{colors[n].upper()}")

    elif color_tag.upper() == "RANDOM_2":
        random_palette = np.random.choice(ct.palettes)

        if re.search("/", random_palette):
            colors_list = random_palette.split("/")
            random_palette = np.random.choice(colors_list)

        colors = "".join(random_palette)

    print(f"NOTE: the palette is {colors}")
    # seaborn plot
    # figure size
    plt.figure(figsize=(fig_width, fig_height))
    sb.scatterplot(x=pc_1, y=pc_2, hue=virus_info, palette=colors, size=virus_info, sizes=(50, 50))
    legend_ncol = int(len(virus_info) / 12)

    # plot point label
    if point_label:
        for index in range(len(virus_info)):
            plt.text(pc_1[index], pc_2[index], virus_info[index], fontsize=8, ha='right', va="bottom")

    # plot legend
    plt.legend(loc="center left",
               title="Group",
               title_fontsize=12,
               fontsize="small",
               bbox_to_anchor=(1, 0.5),
               ncol=legend_ncol,
               frameon=False)

    # plot axis label
    plt.xlabel("PC1", size=12)
    plt.ylabel("PC2", size=12)
    plt.tight_layout()

    plt.savefig(fname=fig_name, dpi=fig_dpi, transparent=tp_tag)
    # plt.show()


def output_palettes():
    print("\n" + "NOTE: All supported palettes in the pipeline:")
    print("Below is a complete list of all palette options. "
          "Most palettes can have the suffix '_r' to indicate the same "
          "palette but reversed order. A few palettes can have '_d' appended at the end which indicates a darker "
          "version of the original palette.")
    print()

    print("Random: 'Generates an equal number of random colors based on the number of items in the data'")
    print("Random_2: 'Randomly select one of all the palettes supported by the system'")
    for palette in ct.palettes:
        print(palette)
