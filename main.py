# -*- coding: utf-8 -*-
# Development Time: 2024-02-27 10:53:33
# Developer: XiaoYang

import argparse
import os.path
import platform
import sys
import time
import pandas as pd
import functions as func
from colorama import Fore, init

example_use = r'''
{}********************** Example of use ***********************{}

    {}THE PROJECT IS IN BETA, DO NOT RELEASE!

{}***************************  End  ***************************{}

'''.format(Fore.GREEN, Fore.RESET, Fore.RED, Fore.GREEN, Fore.RESET)


def starts():
    print(Fore.GREEN + "\n" + "===================================================================")

    print("{}>>> {}Name: Single Amino Acid Polymorphism Statistic (SAAPS)".format(Fore.GREEN, Fore.RESET))

    print(
        "{}>>> {}Description: Count SAP in viral protein sequences.".format(Fore.GREEN, Fore.RESET))

    print("{}>>> {}Version: 0.1 (2024-03-06)".format(Fore.GREEN, Fore.RESET))

    print("{}>>> {}Author: Yang Xiao".format(Fore.GREEN, Fore.RESET))

    print("{}>>> {}Email: Fredrik1999@163.com".format(Fore.GREEN, Fore.RESET))

    print(Fore.GREEN + "===================================================================" + "\n" + Fore.RESET)

    def parameters():

        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            prog="SAPS",            # this parameter can change the name of pipeline in help information
            description="",
            epilog=example_use)

        input_group = parser.add_argument_group("Input Options")
        input_group.add_argument("-i", "--input", type=str, metavar="", dest="input",
                                 help="the path of the input file")
        input_group.add_argument("-sr", "--shannon_result", type=str, metavar="", dest="shannon_result",
                                 help="the path of the shannon result")

        basic_group = parser.add_argument_group("Basic Options")

        basic_group.add_argument("-dg", "--delete_gap", action="store_true", dest="dg",
                                 help="delete the columns contained gaps in sequences")
        basic_group.add_argument("-t", "--translate", action="store_true", dest="translate",
                                 help="translate nucleotides to amino acids")
        basic_group.add_argument("-cs", "--cut_start", type=str, metavar="", dest="cut_start",
                                 help="the start fragment of target sequence, no less than 7 amino acids")
        basic_group.add_argument("-ce", "--cut_end", type=str, metavar="", dest="cut_end",
                                 help="the end fragment of target sequence, no less than 7 amino acids")
        basic_group.add_argument("-osn", "--output_seq_name", action="store_true", dest="name_output",
                                 help="output sequence name, in order to change the sequence name more conveniently")
        basic_group.add_argument("-csn", "--change_seq_name", type=str, metavar="", dest="change_name",
                                 help="change the sequence names, please input the replacement name file path")

        # polymorphism options
        poly_group = parser.add_argument_group("Polymorphism Options")
        poly_group.add_argument("-cp", "--compute_ppm", action="store_true", dest="compute",
                                help="compute polymorphism")
        poly_group.add_argument("-pp", "--plot_ppm", action="store_true", dest="plot_polymorphism",
                                help="plot polymorphism figure")

        # plotting options
        plot_group = parser.add_argument_group("Plotting Options")
        plot_group.add_argument("-fw", "--figure_width", type=int, metavar="", dest="fw",
                                help="The width of the figure")
        plot_group.add_argument("-fh", "--figure_height", type=int, metavar="", dest="fh",
                                help="The height of the figure")
        plot_group.add_argument("-fname", "--figure_name", type=str, metavar="", dest="fname",
                                help="The name of figure")
        plot_group.add_argument("-c", "--color", type=str, metavar="", dest="color",
                                help="The palette of figure, Supports all color schemes in 'Seaborn' package such as "
                                     "'Random', 'Random_2', 'Accent', 'Blues', 'BrBG', 'CMRmap', 'GnBu', "
                                     "'RdYlBu', etc. The default is 'RdYlBu'")
        plot_group.add_argument("-pac", "--print_all_colors", action="store_true", dest="pac",
                                help="Print all supported palettes ")
        plot_group.add_argument("-dpi", "--DPI", type=int, metavar="", dest="dpi",
                                help="The dpi of the figure")
        plot_group.add_argument("-fm", "--format", type=str, metavar="", dest="format",
                                help="The format of the figure")
        plot_group.add_argument("-tp", "--transparent", action="store_true", dest="tp",
                                help="Make the background of the figure transparent")
        plot_group.add_argument("-sl", "--scatter_label", action="store_true", dest="sl",
                                help="Displays scatter plot labels, default is False")

        # shannon entropy option
        se_group = parser.add_argument_group("Shannon Entropy Option")
        se_group.add_argument("-sha", "--shannon", action="store_true", dest="shannon",
                              help="compute shannon entropy")
        se_group.add_argument("-smax", "--shannon_max", type=str, metavar="", dest="shannon_max",
                              help="the maximum of the shannon value (IC)")
        se_group.add_argument("-smin", "--shannon_min", type=str, metavar="", dest="shannon_min",
                              help="the minimum of the shannon value (IC)")
        se_group.add_argument("-pic", "--plot_IC", action="store_true", dest="pic",
                              help="plot IC scatter figure")

        # Dimension reduction option
        dr_group = parser.add_argument_group("Dimension Reduction Options")
        dr_group.add_argument("-oh", "--onehot", action="store_true", dest="one_hot",
                              help="One-Hot Encoding")
        dr_group.add_argument("-pca", "--PCA", action="store_true", dest="pca",
                              help="Dimensionality reduction by PCA, "
                                   "please input the OneHot matrix from the '-oh' parameter")
        dr_group.add_argument("-tsne", "--TSNE", action="store_true", dest="tsne",
                              help="Dimensionality reduction by t-SNE, "
                                   "please input the OneHot matrix from the '-oh' parameter")
        dr_group.add_argument("-tpp", "--tSNE_perplexity", type=int, metavar="", dest="tpp",
                              help="The perplexity in t_SNE, default is 3")
        dr_group.add_argument("-tlr", "--tSNE_learn_rate", type=int, metavar="", dest="tlr",
                              help="The learning rate in t_SNE, default is 200")
        dr_group.add_argument("-tni", "--tSNE_n_iter", type=int, metavar="", dest="tni",
                              help="The n_iter in t_SNE, default is 5000")
        dr_group.add_argument("-trs", "--tSNE_random_state", type=int, metavar="", dest="trs",
                              help="The random state in t_SNE, default is 1")
        dr_group.add_argument("-pc", "--plot_cluster", action="store_true", dest="pc",
                              help="Plot clustering scatter plot according to the result of dimensionality reduction")

        # output options
        output_group = parser.add_argument_group("Output Options")
        output_group.add_argument("-o", "--output_dir", type=str, metavar="", dest="output", default="",
                                  help="the directory path of output file")
        output_group.add_argument("-pre", "--prefix", type=str, metavar="", dest="prefix",
                                  help="The prefix of the output file name")

        my_args = parser.parse_args(sys.argv[1:])
        return my_args

    my_args = parameters()
    result_content = ""
    result_path = "~"                           # workpath in default
    output_prefix = ""                          # the prefix of output file
    user_sys_os = platform.uname()[0]           # system operation type
    start_time = time.time()                    # pipeline start time

    # Judge system operation type
    # macOS
    if user_sys_os == "Darwin":
        print("{}Note: The current operating system is: {}, remove hidden file path in "
              "default.{}".format(Fore.CYAN, user_sys_os, Fore.RESET))
    # Windows
    elif user_sys_os == "Windows":
        print("{}Note: The current operating system is: {}, if there are any hidden files please manually delete.{}"
              .format(Fore.CYAN, user_sys_os, Fore.RESET))
    # Linux
    elif user_sys_os == "Linux":
        print("{}Note: The current operating system is: {}, if there are any hidden files please manually delete.{}"
              .format(Fore.CYAN, user_sys_os, Fore.RESET))

    # parameter instructions
    # '-o' parameter
    # output directory
    if my_args.output:
        if my_args.output == "":                                                # blank path
            print(Fore.RED + "Error! The output path is not specified!")
            sys.exit()

        result_path = os.path.abspath(my_args.output)
        if os.path.exists(result_path):
            pass

        elif not os.path.exists(result_path):
            print("%sWarning: The current path does not exist, automatically create the path: %s %s" %
                  (Fore.YELLOW, Fore.BLUE, os.path.realpath(my_args.output)))
            os.makedirs(result_path)

    print("NOTE: The output path is %s'%s'" % (Fore.BLUE, result_path))

    # '-pre' parameter
    # the prefix of output
    if my_args.prefix:
        output_prefix = my_args.prefix

    # '-i' parameter
    if my_args.input:
        if my_args.shannon:
            pass
        elif my_args.one_hot:
            pass
        elif my_args.pca:
            pass
        elif my_args.tsne:
            pass
        else:
            if os.path.isdir(my_args.input):                                # input is a folder
                print(Fore.RED + "Error!The folder is currently not supported!")
                sys.exit()

            elif os.path.isfile(my_args.input):                             # input is a file
                my_file = func.read_file(my_args.input)
                standard_fasta = func.fasta_doc_std(my_file)
                result_content = standard_fasta
            else:
                print(Fore.RED + "Error! Please make sure that the input is a file, "
                                 "instead of links and so on!")
                sys.exit()

    # '-dg' parameter
    # delete the gaps
    if my_args.dg:
        dg_result_path = result_path + os.sep + output_prefix + "NoGapSeq.txt"
        dg_result = func.delete_gap(result_content)
        func.list_to_save(dg_result, dg_result_path)

        print("NOTE: Sequence processing completed!")

    # '-cs' and '-ce' parameters
    # the principle of obtaining target sequence fragment: the first sequence in the sequence file will be regard as the
    # reference sequence, and according to the location of the start and end fragment of the reference sequence, the
    # remaining sequences will be cut.
    # '-cs' and '-ce' represent the start and end fragment of the reference sequence, note the length of fragment should
    # not be less than 8 amino acid.
    if my_args.cut_start and my_args.cut_end:
        # judge the type of input
        try:
            start_info, end_info = int(my_args.cut_start), int(my_args.cut_end)
            # start_info = start_info - 1
            # end_info = end_info - 1
            input_tag = "int"
        except ValueError:
            start_info, end_info = my_args.cut_start.upper(), my_args.cut_end.upper()       # upper the seq
            input_tag = "str"

        print("Input information of cutting function:")
        print("'-cs':%s" % start_info)                           # print the input seq or site
        print("'-ce':%s" % end_info)

        if input_tag == "str":

            # judge the length of the input fragments
            # To get an accurate, a certain length is required
            if any([len(start_info) < 5, len(end_info) < 5]):                 # fragment length is less than 5
                print("%sERROR: the seq length ('-cs' or '-ce') you input is less than 5, "
                      "the length about 5~10 is recommended" % Fore.YELLOW)
                sys.exit()
            elif any([5 <= len(start_info) <= 10, 5 <= len(end_info) <= 10]):     # fragment length is in 5~10
                pass
            elif any([len(start_info) > 10, len(end_info) > 10]):                 # fragment length is more than 10
                print("%sWarning: the seq length ('-cs' or '-ce') you input is more than 10, "
                      "the length about 5~10 is recommended" % Fore.YELLOW)

        target_seq = func.seq_cut(start_info, end_info, result_content, input_tag)
        result_content = target_seq                                     # receive a list
        cut_result_path = result_path + os.sep + output_prefix + "SeqCut.txt"           # txt in default
        func.list_to_save(result_content, cut_result_path)

    # '-ose' parameter
    # output sequence name
    if my_args.name_output:
        name_result = func.output_name(result_content)
        name_result_path = result_path + os.sep + output_prefix + "SeqName.csv"
        name_result.to_csv(name_result_path, index=False)                                         # csv in default

    # '-t' parameter
    # translate the nucleotides to amino acids
    if my_args.translate:
        nt_tag = func.codon_check(result_content)
        translated_result = func.translate_to_protein(result_content, nt_tag)
        trans_result_path = result_path + os.sep + output_prefix + "Translated.txt"     # txt in default
        func.list_to_save(translated_result, trans_result_path)

    # '-csn' parameter
    # change the sequence names
    if my_args.change_name:
        change_result = func.change_name(result_content, my_args.change_name)
        csn_result_path = result_path + os.sep + output_prefix + "NewName.txt"
        func.list_to_save(change_result, csn_result_path)

    # '-cp' parameter
    # compute the polymorphism
    if my_args.compute:
        computation_result = func.compute_polymorphism(result_content)                  # it should be a dataframe
        seq_matrix, concise_result, detail_result = computation_result                 # it is a tuple
        # output concise result
        concise_result_path = result_path + os.sep + output_prefix + "Concise_Result.txt"
        func.list_to_save(concise_result, concise_result_path)
        # output detail result
        detail_result_path = result_path + os.sep + output_prefix + "Detail_Result.txt"
        func.list_to_save(detail_result, detail_result_path)
        # output sequence matrix
        seq_matrix_path = result_path + os.sep + output_prefix + "Seq_Matrix.csv"
        seq_matrix.to_csv(seq_matrix_path, index=False)

        # '-pp' parameter
        # plot polymorphism
        if my_args.plot_polymorphism:
            print(f"NOTE: Activating Polymorphism Plotting function!")
            figure_format = "pdf"
            figure_name = result_path + os.sep + output_prefix + "SAPs" + "." + figure_format
            dpi = 300
            transparent = False
            figure_width = 14
            figure_height = 9

            if my_args.format:
                figure_format = my_args.format
            if my_args.fname:
                figure_name = result_path + os.sep + output_prefix + my_args.fname + "." + figure_format
            if my_args.dpi:
                dpi = my_args.dpi
            if my_args.tp:
                transparent = True
            if my_args.fw:
                figure_width = my_args.fw
            if my_args.fh:
                figure_height = my_args.fh

            func.polymorphism_figure(seq_matrix, concise_result, figure_name, dpi, transparent,
                                     fig_width=figure_width,
                                     fig_height=figure_height)
            print("\n" + "NOTE: Plotting Completed!")

    # '-sha' parameter
    # compute shannon entropy
    # input should be an aligned fasta file
    if my_args.shannon:
        # parameter '-smax' and '-smin'
        if my_args.shannon_max and my_args.shannon_min:
            shannon_maximum = float(my_args.shannon_max)
            shannon_minimum = float(my_args.shannon_min)
            result_content = func.shannon_entropy(my_args.input, shannon_maximum, shannon_minimum)

        else:
            result_content = func.shannon_entropy(my_args.input)

        # result in 'csv' format in default
        shannon_result_path = result_path + os.sep + output_prefix + "Shannon_IC_Result.csv"
        result_content.to_csv(shannon_result_path, index=False, header=True)

        # plot IC scatter
        if my_args.pic:
            print(f"NOTE: Activating IC Scatter Plotting function!")
            figure_format = "pdf"
            figure_name = result_path + os.sep + output_prefix + "IC-Plot" + "." + figure_format
            dpi = 300
            transparent = False
            figure_width = 10
            figure_height = 8
            palette = "RdYlBu"

            if my_args.format:
                figure_format = my_args.format
            if my_args.fname:
                figure_name = result_path + os.sep + output_prefix + my_args.fname + "." + figure_format
            if my_args.dpi:
                dpi = my_args.dpi
            if my_args.tp:
                transparent = True
            if my_args.fw:
                figure_width = my_args.fw
            if my_args.fh:
                figure_height = my_args.fh
            if my_args.color:
                palette = my_args.color

            func.ic_scatter_figure(result_content,
                                   fig_name=figure_name,
                                   fig_dpi=dpi,
                                   tp_tag=transparent,
                                   fig_width=figure_width,
                                   fig_height=figure_height,
                                   color_tag=palette)                  # input a dataframe
            print("\n" + "NOTE: Plotting Completed!")

    # '-sr' parameter
    # create global variable: key_sites
    key_sites = list()
    if my_args.shannon_result:
        file_type = func.file_type_judge(my_args.shannon_result)
        # determine the file type
        if file_type == "CSV":
            shannon_df = pd.read_csv(my_args.shannon_result, header=0)
        elif file_type == "XLSX":
            shannon_df = pd.read_excel(my_args.shannon_result, header=0)
        else:
            print(Fore.RED + "ERROR! The file type:%s is not supported right now!" % file_type)
            sys.exit()

        # '-smax' and '-smin' parameters
        if all([my_args.shannon_max, my_args.shannon_min]) is True:
            maximum = float(my_args.shannon_max)
            minimum = float(my_args.shannon_min)
            key_sites = func.shannon_filter(shannon_df, maximum, minimum)
        else:
            key_sites = func.shannon_filter(shannon_df)

    # '-oh' parameter
    # one-hot coding
    if my_args.one_hot:
        file_type = func.file_type_judge(my_args.input)

        if file_type == "CSV":
            seq_mat = pd.read_csv(my_args.input, header=0)
        elif file_type == "XLSX":
            seq_mat = pd.read_excel(my_args.input, header=0)
        else:
            print(Fore.RED + "ERROR! The file type:%s is not supported right now!" % file_type)
            sys.exit()

        # index starts from 0, so the value in key_sites need to minus 1
        key_sites = [int(x)-1 for x in key_sites]
        # filter the seq matrix
        # One-Hot Encoding
        raw_aa_matrix, onehot_matrix_concise, onehot_matrix, site_matrix = func.one_hot_encoding(seq_mat, key_sites)
        onehot_result_path = result_path + os.sep + output_prefix

        raw_aa_matrix.to_csv(onehot_result_path + "OneHot-Original.csv", index=True, sep=",", header=True)
        onehot_matrix_concise.to_csv(onehot_result_path + "OneHot-ConciseMatrix.csv", index=True, sep=",", header=True)
        onehot_matrix.to_csv(onehot_result_path + "OneHot-Transform.csv", index=True, sep=",", header=True)
        site_matrix.to_csv(onehot_result_path + "OneHot-ForIntersect.csv", index=True, sep=",", header=True)

    # '-pca' parameter
    # dimension reduction
    if my_args.pca:
        print(f"NOTE: Activating PCA function!")
        pca_result = func.pca(my_args.input)
        pca_result_path = result_path + os.sep + output_prefix + "PCA.csv"
        pca_result.to_csv(pca_result_path)
        print("NOTE: PCA Completed!")

        # Plot Scatter
        if my_args.pc:
            print(f"NOTE: Activating PCA Plotting function!")
            pca_result["Virus"] = pca_result.index
            scatter_label = False
            figure_format = "pdf"
            figure_name = result_path + os.sep + output_prefix + "PCA-Plot" + "." + figure_format
            dpi = 300
            transparent = False
            figure_width = 8
            figure_height = 6
            palette = "RdYlBu"

            if my_args.format:
                figure_format = my_args.format
            if my_args.fname:
                figure_name = result_path + os.sep + output_prefix + my_args.fname + "." + figure_format
            if my_args.dpi:
                dpi = my_args.dpi
            if my_args.tp:
                transparent = True
            if my_args.fw:
                figure_width = my_args.fw
            if my_args.fh:
                figure_height = my_args.fh
            if my_args.sl:
                scatter_label = True
            if my_args.color:
                palette = my_args.color

            func.cluster_figure(pca_result,
                                point_label=scatter_label,
                                fig_name=figure_name,
                                fig_dpi=dpi,
                                tp_tag=transparent,
                                fig_width=figure_width,
                                fig_height=figure_height,
                                color_tag=palette
                                )

            print("\n" + "NOTE: Plotting Completed!")

    # '-tsne' parameter
    # dimension reduction
    if my_args.tsne:
        print("NOTE: Activating t-SNE function!")
        # tSNE default parameters
        perplexity = 3
        learning_rate = 200
        n_iter = 5000
        random_state = 1

        # tSNE parameters
        if my_args.tpp:
            perplexity = my_args.tpp
        if my_args.tlr:
            learning_rate = my_args.tlr
        if my_args.tni:
            n_iter = my_args.tni
        if my_args.trs:
            random_state = my_args.trs

        print("\n" + "t-SNE Parameters:")
        print(f"Perplexity: {perplexity}, "
              f"Learning Rate: {learning_rate}, "
              f"n_iter: {n_iter}, "
              f"Random State: {random_state}")

        tsne_result = func.tsne(my_args.input,
                                perp=perplexity,
                                learn_rate=learning_rate,
                                niter=n_iter,
                                rand_state=random_state)

        tsne_result_path = result_path + os.sep + output_prefix + "tSNE.csv"
        tsne_result.to_csv(tsne_result_path)
        print("NOTE: t-SNE Completed!")

        # PLOTTING
        if my_args.pc:
            # Plot Scatter
            print(f"NOTE: Activating PCA Plotting function!")
            tsne_result["Virus"] = tsne_result.index
            scatter_label = False
            figure_format = "pdf"
            figure_name = result_path + os.sep + output_prefix + "tSNE-Plot" + "." + figure_format
            dpi = 300
            transparent = False
            figure_width = 8
            figure_height = 6
            palette = "RdYlBu"

            if my_args.format:
                figure_format = my_args.format
            if my_args.fname:
                figure_name = result_path + os.sep + output_prefix + my_args.fname + "." + figure_format
            if my_args.dpi:
                dpi = my_args.dpi
            if my_args.tp:
                transparent = True
            if my_args.fw:
                figure_width = my_args.fw
            if my_args.fh:
                figure_height = my_args.fh
            if my_args.sl:
                scatter_label = True
            if my_args.color:
                palette = my_args.color

            func.cluster_figure(tsne_result,
                                point_label=scatter_label,
                                fig_name=figure_name,
                                fig_dpi=dpi,
                                tp_tag=transparent,
                                fig_width=figure_width,
                                fig_height=figure_height,
                                color_tag=palette)

            print("\n" + "NOTE: Plotting Completed!")

    if my_args.pac:
        func.output_palettes()
        sys.exit()

    end_time = time.time()
    execution_time = end_time - start_time
    print("NOTE: This operation takes {}{:.3f}{} seconds!".format(Fore.BLUE, execution_time, Fore.RESET))
    sys.exit()


if __name__ == "__main__":
    init(wrap=True, autoreset=True)
    starts()
