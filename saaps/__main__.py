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

  (1) Calculate single amino acid polymorphisms in a sequence file:
      SAAPS -i your/sequence/path -c -o your/output/path

  (2) Calculate and plot single amino acid polymorphisms in a sequence file:
      SAAPS -i your/sequence/path -c -cp -o your/output/path
      
  Tip: the input(file) and output(directory) is recommended absolute path.

  Above are just simple examples, detailed instruction in website: {}https://github.com/xiaosheep01/SAAPS
{}***************************  End  ***************************{}

'''.format(Fore.GREEN, Fore.RESET, Fore.CYAN, Fore.GREEN, Fore.RESET)


def starts():
    print(Fore.GREEN + "\n" + "===================================================================")

    print("{}>>> {}Name: Single Amino Acid Polymorphism Statistic (SAAPS)".format(Fore.GREEN, Fore.RESET))

    print(
        "{}>>> {}Description: Count SAP in protein alignment sequences".format(Fore.GREEN, Fore.RESET))

    print("{}>>> {}Version: 1.4.2 (2024-03-26)".format(Fore.GREEN, Fore.RESET))

    print("{}>>> {}Author: Yang Xiao".format(Fore.GREEN, Fore.RESET))

    print("{}>>> {}Email: Fredrik1999@163.com".format(Fore.GREEN, Fore.RESET))

    print(Fore.GREEN + "===================================================================" + "\n" + Fore.RESET)

    def parameters():

        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            prog="SAAPS",            # this parameter can change the name of pipeline in help information
            description="",
            epilog=example_use)

        input_group = parser.add_argument_group("Input Options")
        input_group.add_argument("-i", "--input", type=str, metavar="", dest="input",
                                 help="the path of the input file")
        input_group.add_argument("-sr", "--shannon_result", type=str, metavar="", dest="shannon_result",
                                 help="the path of the shannon result")

        basic_group = parser.add_argument_group("Basic Options")
        basic_group.add_argument("-dg", "--delete_gap", action="store_true", dest="dg",
                                 help="delete the columns containing gaps in sequences")
        basic_group.add_argument("-t", "--translate", action="store_true", dest="translate",
                                 help="translate nucleotides to amino acids")
        basic_group.add_argument("-cs", "--cut_start", type=str, metavar="", dest="cut_start",
                                 help="the start fragment/site of target sequence, no less than 7 amino acids")
        basic_group.add_argument("-ce", "--cut_end", type=str, metavar="", dest="cut_end",
                                 help="the end fragment/site of target sequence, no less than 7 amino acids")
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
        dr_group.add_argument("-aaindex", "--AAindex", action="store_true", dest="aaindex",
                              help="AAindex Encoding")
        dr_group.add_argument("-pca", "--PCA", action="store_true", dest="pca",
                              help="Dimensionality reduction by PCA, "
                                   "please input the OneHot matrix from the '-oh' parameter")
        dr_group.add_argument("-tsne", "--TSNE", action="store_true", dest="tsne",
                              help="Dimensionality reduction by t-SNE, "
                                   "please input the OneHot matrix from the '-oh' parameter")
        dr_group.add_argument("-ica", "--ICA", action="store_true", dest="ica",
                              help="Dimensionality reduction by ICA, "
                                   "please input the OneHot matrix from the '-oh' parameter")
        dr_group.add_argument("-fa", "--FactorAnalysis", action="store_true", dest="fa",
                              help="Dimensionality reduction by FactorAnalysis, "
                                   "please input the OneHot matrix from the '-oh' parameter")
        dr_group.add_argument("-mds", "--MDS", action="store_true", dest="mds",
                              help="Dimensionality reduction by MDS, "
                                   "please input the OneHot matrix from the '-oh' parameter")
        dr_group.add_argument("-isomap", "--ISOMAP", action="store_true", dest="isomap",
                              help="Dimensionality reduction by Isomap, "
                                   "please input the OneHot matrix from the '-oh' parameter")
        dr_group.add_argument("-spectremb", "--SPECTREMB", action="store_true", dest="spectremb",
                              help="Dimensionality reduction by Spectremb, "
                                   "please input the OneHot matrix from the '-oh' parameter")
        dr_group.add_argument("-lle", "--LLE", action="store_true", dest="lle",
                              help="Dimensionality reduction by LLE, "
                                   "please input the OneHot matrix from the '-oh' parameter")
        dr_group.add_argument("-kpca", "--kernel_PCA", action="store_true", dest="kpca",
                              help="Dimensionality reduction by Kernel PCA, "
                                   "please input the OneHot matrix from the '-oh' parameter")
        dr_group.add_argument("-rfp", "--RandomForest_PCA", action="store_true", dest="rfp",
                              help="Dimensionality reduction by Kernel RandomForest and PCA, "
                                   "please input the OneHot matrix from the '-oh' parameter")
        dr_group.add_argument("-td", "--tSNE_dimension", type=int, metavar="", dest="td",
                              help="The dimension in t_SNE, default is 2")
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
        dr_group.add_argument("-pc2", "--plot_cluster2", action="store_true", dest="pc2",
                              help="test plotting")

        # Clustering Options
        cluster_group = parser.add_argument_group("Clustering Options")
        cluster_group.add_argument("-dbs", "--DBSCAN", action="store_true", dest="dbscan",
                                   help="Density-Based Spatial Clustering of Applications with Noise (DBSCAN)")
        cluster_group.add_argument("-eps", "--epsilon", type=float, metavar="", dest="eps",
                                   help="The epsilon in DBSCAN, default is 0.5")
        cluster_group.add_argument("-minsap", "--min_samples", type=int, metavar="", dest="minsap",
                                   help="The min_samples state in DBSCAN, default is 6")
        cluster_group.add_argument("-ops", "--OPTICS", action="store_true", dest="optics",
                                   help="Ordering Points To Identify the Clustering Structure (OPTICS)")
        cluster_group.add_argument("-km", "--K_Mean", action="store_true", dest="kmean",
                                   help="K-Mean")
        cluster_group.add_argument("-agg", "--Agglomerative", action="store_true", dest="agglo",
                                   help="Agglomerative")
        cluster_group.add_argument("-nc", "--n_cluster", type=int, metavar="", dest="nc",
                                   help="the number of cluster in K_Mean and Agglomerative")
        cluster_group.add_argument("-pts", "--plot_3d_scatter", action="store_true", dest="pts",
                                   help="plot three-dimension scatter of input data")

        # output options
        output_group = parser.add_argument_group("Output Options")
        output_group.add_argument("-o", "--output_dir", type=str, metavar="", dest="output", default="",
                                  help="The directory path of output file, "
                                       "the default is the current user document folder")
        output_group.add_argument("-pre", "--prefix", type=str, metavar="", dest="prefix",
                                  help="The prefix of the output file name")

        my_args = parser.parse_args(sys.argv[1:])
        return my_args

    my_args = parameters()
    result_content = ""
    result_path = os.path.join(os.path.expanduser("~"), "Documents")        # workpath in default
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

    print("NOTE: The output path is %s'%s'%s" % (Fore.BLUE, result_path, Fore.RESET))

    # '-pre' parameter
    # the prefix of output
    if my_args.prefix:
        output_prefix = my_args.prefix

    # '-i' parameter
    if my_args.input:
        if any([my_args.shannon, my_args.one_hot, my_args.pca, my_args.tsne, my_args.dbscan, my_args.optics,
                my_args.kmean, my_args.agglo, my_args.pts, my_args.aaindex, my_args.isomap, my_args.spectremb,
                my_args.lle, my_args.kpca, my_args.rfp, my_args.ica, my_args.fa, my_args.mds]):
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
    # '-cs' and '-ce' represent the start and end fragment/site of the reference sequence, note the length of
    # fragment should not be less than 8 amino acid.
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

        print("NOTE: Input information of cutting function:")
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
                      "the length about 5~10 is recommended%s" % (Fore.YELLOW, Fore.RESET))

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
        csn_result_path = result_path + os.sep + output_prefix + "NewNameSeq.txt"
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
            if my_args.format:
                figure_format = my_args.format

            figure_name = result_path + os.sep + output_prefix + "SAPs" + "." + figure_format
            dpi = 300
            transparent = False
            figure_width = 14
            figure_height = 9

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
            if my_args.format:
                figure_format = my_args.format

            figure_name = result_path + os.sep + output_prefix + "IC-Plot" + "." + figure_format
            dpi = 300
            transparent = False
            figure_width = 10
            figure_height = 8
            palette = "RdYlBu"

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
        print("NOTE: Activating OneHot Encoding!")
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
        print("NOTE: OneHot Encoding Completed!")

    if my_args.aaindex:
        print("NOTE: Activating AAindex Encoding!")
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
        raw_aa_matrix, AAindex_matrix_concise, AAindex_matrix, site_matrix = func.aaindex_encoding(seq_mat, key_sites)
        AAindex_result_path = result_path + os.sep + output_prefix

        raw_aa_matrix.to_csv(AAindex_result_path + "AAindex-Original.csv", index=True, sep=",", header=True)
        AAindex_matrix_concise.to_csv(AAindex_result_path + "AAindex-ConciseMatrix.csv", index=True, sep=",", header=True)
        AAindex_matrix.to_csv(AAindex_result_path + "AAindex-Transform.csv", index=True, sep=",", header=True)
        site_matrix.to_csv(AAindex_result_path + "AAindex-ForIntersect.csv", index=True, sep=",", header=True)
        print("NOTE: AAindex Encoding Completed!")

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
            if my_args.format:
                figure_format = my_args.format

            figure_name = result_path + os.sep + output_prefix + "PCA-2d-Plot" + "." + figure_format
            dpi = 300
            transparent = False
            figure_width = 8
            figure_height = 6
            palette = "RdYlBu"

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

        if my_args.pc2:
            print(f"{Fore.LIGHTYELLOW_EX}NOTE: Testing Part Start{Fore.RESET}")
            figure_name = result_path + os.sep + output_prefix + "PCA-3d-Plot" + ".pdf"
            func.clustering_plot2(pca_result, fig_name=figure_name)

            print(f"{Fore.LIGHTYELLOW_EX}NOTE: Testing Part Ends{Fore.RESET}")

    # '-tsne' parameter
    # dimension reduction
    if my_args.tsne:
        print("NOTE: Activating t-SNE function!")
        # tSNE default parameters
        dimension = 2
        perplexity = 3
        learning_rate = 200
        n_iter = 5000
        # random_state = 1

        # tSNE parameters
        if my_args.td:
            dimension = my_args.td
        if my_args.tpp:
            perplexity = my_args.tpp
        if my_args.tlr:
            learning_rate = my_args.tlr
        if my_args.tni:
            n_iter = my_args.tni
        if my_args.trs:
            random_state = my_args.trs

        print("\n" + "t-SNE Parameters:")
        print(f"Dimension:{dimension}, "
              f"Perplexity: {perplexity}, "
              f"Learning Rate: {learning_rate}, "
              f"n_iter: {n_iter}, ")

        tsne_result = func.tsne(my_args.input,
                                n_dimension=dimension,
                                perp=perplexity,
                                learn_rate=learning_rate,
                                niter=n_iter)

        tsne_result_path = result_path + os.sep + output_prefix + "tSNE.csv"
        tsne_result.to_csv(tsne_result_path)
        print("NOTE: t-SNE Completed!")

        # PLOTTING
        if my_args.pc:
            # Plot Scatter
            print(f"NOTE: Activating tSNE Plotting function!")
            tsne_result["Virus"] = tsne_result.index
            scatter_label = False
            figure_format = "pdf"
            if my_args.format:
                figure_format = my_args.format

            figure_name = result_path + os.sep + output_prefix + "tSNE-Plot" + "." + figure_format
            dpi = 300
            transparent = False
            figure_width = 8
            figure_height = 6
            palette = "RdYlBu"

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

        if my_args.pc2:
            print(f"{Fore.LIGHTYELLOW_EX}NOTE: Testing Part Start{Fore.RESET}")
            figure_name = result_path + os.sep + output_prefix + "tSNE-Plot" + ".pdf"
            func.clustering_plot2(tsne_result, fig_name=figure_name)

            print(f"{Fore.LIGHTYELLOW_EX}NOTE: Testing Part Ends{Fore.RESET}")

    # '-isomap' parameter
    # dimension reduction
    if my_args.isomap:
        print("NOTE: Activating Isomap function!")
        # isomap default parameters
        dimension = 2

        # isomap parameters
        if my_args.td:
            dimension = my_args.td

        print("\n" + "Isomap Parameters:")
        print(f"Dimension:{dimension}")

        iso_result = func.isomapping(my_args.input, n_dimension=dimension)

        iso_result_path = result_path + os.sep + output_prefix + "isomap.csv"
        iso_result.to_csv(iso_result_path)
        print("NOTE: Isomap Completed!")

        # PLOTTING
        if my_args.pc:
            # Plot Scatter
            print(f"NOTE: Activating PCA Plotting function!")
            iso_result["Virus"] = iso_result.index
            scatter_label = False
            figure_format = "pdf"
            if my_args.format:
                figure_format = my_args.format

            figure_name = result_path + os.sep + output_prefix + "Isomap-Plot" + "." + figure_format
            dpi = 300
            transparent = False
            figure_width = 8
            figure_height = 6
            palette = "RdYlBu"

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

            func.cluster_figure(iso_result,
                                point_label=scatter_label,
                                fig_name=figure_name,
                                fig_dpi=dpi,
                                tp_tag=transparent,
                                fig_width=figure_width,
                                fig_height=figure_height,
                                color_tag=palette)

            print("\n" + "NOTE: Plotting Completed!")

        if my_args.pc2:
            print(f"{Fore.LIGHTYELLOW_EX}NOTE: Testing Part Start{Fore.RESET}")
            figure_name = result_path + os.sep + output_prefix + "Isomap-Plot" + ".pdf"
            func.clustering_plot2(iso_result, fig_name=figure_name)

            print(f"{Fore.LIGHTYELLOW_EX}NOTE: Testing Part Ends{Fore.RESET}")

    # '-spectremb'
    if my_args.spectremb:
        print("NOTE: Activating Spectremb function!")
        # spectremb default parameters
        dimension = 2

        # spectremb parameters
        if my_args.td:
            dimension = my_args.td

        print("\n" + "Spectremb Parameters:")
        print(f"Dimension:{dimension}")

        spectremb_result = func.spectremb(my_args.input, n_dimension=dimension)

        spectremb_result_path = result_path + os.sep + output_prefix + "spectremb.csv"
        spectremb_result.to_csv(spectremb_result_path)
        print("NOTE: Spectremb Completed!")

        # PLOTTING
        if my_args.pc:
            # Plot Scatter
            print(f"NOTE: Activating Spectremb Plotting function!")
            spectremb_result["Virus"] = spectremb_result.index
            scatter_label = False
            figure_format = "pdf"
            if my_args.format:
                figure_format = my_args.format

            figure_name = result_path + os.sep + output_prefix + "Spectremb-Plot" + "." + figure_format
            dpi = 300
            transparent = False
            figure_width = 8
            figure_height = 6
            palette = "RdYlBu"

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

            func.cluster_figure(spectremb_result,
                                point_label=scatter_label,
                                fig_name=figure_name,
                                fig_dpi=dpi,
                                tp_tag=transparent,
                                fig_width=figure_width,
                                fig_height=figure_height,
                                color_tag=palette)

            print("\n" + "NOTE: Plotting Completed!")

        if my_args.pc2:
            print(f"{Fore.LIGHTYELLOW_EX}NOTE: Testing Part Start{Fore.RESET}")
            figure_name = result_path + os.sep + output_prefix + "Spectremb-Plot" + ".pdf"
            func.clustering_plot2(spectremb_result, fig_name=figure_name)

            print(f"{Fore.LIGHTYELLOW_EX}NOTE: Testing Part Ends{Fore.RESET}")

    if my_args.lle:
        print("NOTE: Activating LLE function!")
        # LLE default parameters
        dimension = 2

        # LLE parameters
        if my_args.td:
            dimension = my_args.td

        print("\n" + "LLE Parameters:")
        print(f"Dimension:{dimension}")

        LLE_result = func.my_lle(my_args.input, n_dimension=dimension)

        LLE_result_path = result_path + os.sep + output_prefix + "LLE.csv"
        LLE_result.to_csv(LLE_result_path)
        print("NOTE: LLE Completed!")

        # PLOTTING
        if my_args.pc:
            # Plot Scatter
            print(f"NOTE: Activating LLE Plotting function!")
            LLE_result["Virus"] = LLE_result.index
            scatter_label = False
            figure_format = "pdf"
            if my_args.format:
                figure_format = my_args.format

            figure_name = result_path + os.sep + output_prefix + "LLE-Plot" + "." + figure_format
            dpi = 300
            transparent = False
            figure_width = 8
            figure_height = 6
            palette = "RdYlBu"

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

            func.cluster_figure(LLE_result,
                                point_label=scatter_label,
                                fig_name=figure_name,
                                fig_dpi=dpi,
                                tp_tag=transparent,
                                fig_width=figure_width,
                                fig_height=figure_height,
                                color_tag=palette)

            print("\n" + "NOTE: Plotting Completed!")

        if my_args.pc2:
            print(f"{Fore.LIGHTYELLOW_EX}NOTE: Testing Part Start{Fore.RESET}")
            figure_name = result_path + os.sep + output_prefix + "LLE-Plot" + ".pdf"
            func.clustering_plot2(LLE_result, fig_name=figure_name)

            print(f"{Fore.LIGHTYELLOW_EX}NOTE: Testing Part Ends{Fore.RESET}")

    if my_args.kpca:
        print("NOTE: Activating Kernel PCA function!")
        # kPCA default parameters
        dimension = 2

        # kPCA parameters
        if my_args.td:
            dimension = my_args.td

        print("\n" + "kPCA Parameters:")
        print(f"Dimension:{dimension}")

        kPCA_result = func.kernelPCA(my_args.input, n_dimension=dimension)

        kPCA_result_path = result_path + os.sep + output_prefix + "kPCA.csv"
        kPCA_result.to_csv(kPCA_result_path)
        print("NOTE: kPCA Completed!")

        # PLOTTING
        if my_args.pc:
            # Plot Scatter
            print(f"NOTE: Activating kPCA Plotting function!")
            kPCA_result["Virus"] = kPCA_result.index
            scatter_label = False
            figure_format = "pdf"
            if my_args.format:
                figure_format = my_args.format

            figure_name = result_path + os.sep + output_prefix + "kPCA-Plot" + "." + figure_format
            dpi = 300
            transparent = False
            figure_width = 8
            figure_height = 6
            palette = "RdYlBu"

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

            func.cluster_figure(kPCA_result,
                                point_label=scatter_label,
                                fig_name=figure_name,
                                fig_dpi=dpi,
                                tp_tag=transparent,
                                fig_width=figure_width,
                                fig_height=figure_height,
                                color_tag=palette)

            print("\n" + "NOTE: Plotting Completed!")

        if my_args.pc2:
            print(f"{Fore.LIGHTYELLOW_EX}NOTE: Testing Part Start{Fore.RESET}")
            figure_name = result_path + os.sep + output_prefix + "kPCA-Plot" + ".pdf"
            func.clustering_plot2(kPCA_result, fig_name=figure_name)

            print(f"{Fore.LIGHTYELLOW_EX}NOTE: Testing Part Ends{Fore.RESET}")

    if my_args.rfp:
        print("NOTE: Activating Random Forest + PCA function!")
        # rfp default parameters
        dimension = 2

        # rfp parameters
        if my_args.td:
            dimension = my_args.td

        print("\n" + "rfp Parameters:")
        print(f"Dimension:{dimension}")

        rfp_result = func.randomForest_PCA(my_args.input, n_dimension=dimension)

        rfp_result_path = result_path + os.sep + output_prefix + "rfp.csv"
        rfp_result.to_csv(rfp_result_path)
        print("NOTE: rfp Completed!")

        # PLOTTING
        if my_args.pc:
            # Plot Scatter
            print(f"NOTE: Activating rfp Plotting function!")
            rfp_result["Virus"] = rfp_result.index
            scatter_label = False
            figure_format = "pdf"
            if my_args.format:
                figure_format = my_args.format

            figure_name = result_path + os.sep + output_prefix + "rfp-Plot" + "." + figure_format
            dpi = 300
            transparent = False
            figure_width = 8
            figure_height = 6
            palette = "RdYlBu"

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

            func.cluster_figure(rfp_result,
                                point_label=scatter_label,
                                fig_name=figure_name,
                                fig_dpi=dpi,
                                tp_tag=transparent,
                                fig_width=figure_width,
                                fig_height=figure_height,
                                color_tag=palette)

            print("\n" + "NOTE: Plotting Completed!")

        if my_args.pc2:
            print(f"{Fore.LIGHTYELLOW_EX}NOTE: Testing Part Start{Fore.RESET}")
            figure_name = result_path + os.sep + output_prefix + "rfp-Plot" + ".pdf"
            func.clustering_plot2(rfp_result, fig_name=figure_name)

            print(f"{Fore.LIGHTYELLOW_EX}NOTE: Testing Part Ends{Fore.RESET}")

    # ICA
    if my_args.ica:
        print("NOTE: Activating ICA function!")
        # ica default parameters
        dimension = 2

        # ica parameters
        if my_args.td:
            dimension = my_args.td

        print("\n" + "ICA Parameters:")
        print(f"Dimension:{dimension}")

        ica_result = func.ica(my_args.input, n_dimension=dimension)

        ica_result_path = result_path + os.sep + output_prefix + "ica.csv"
        ica_result.to_csv(ica_result_path)
        print("NOTE: ICA Completed!")

        # PLOTTING
        if my_args.pc:
            # Plot Scatter
            print(f"NOTE: Activating PCA Plotting function!")
            ica_result["Virus"] = ica_result.index
            scatter_label = False
            figure_format = "pdf"
            if my_args.format:
                figure_format = my_args.format

            figure_name = result_path + os.sep + output_prefix + "ICA-Plot" + "." + figure_format
            dpi = 300
            transparent = False
            figure_width = 8
            figure_height = 6
            palette = "RdYlBu"

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

            func.cluster_figure(ica_result,
                                point_label=scatter_label,
                                fig_name=figure_name,
                                fig_dpi=dpi,
                                tp_tag=transparent,
                                fig_width=figure_width,
                                fig_height=figure_height,
                                color_tag=palette)

            print("\n" + "NOTE: Plotting Completed!")

        if my_args.pc2:
            print(f"{Fore.LIGHTYELLOW_EX}NOTE: Testing Part Start{Fore.RESET}")
            figure_name = result_path + os.sep + output_prefix + "ICA-Plot" + ".pdf"
            func.clustering_plot2(ica_result, fig_name=figure_name)

            print(f"{Fore.LIGHTYELLOW_EX}NOTE: Testing Part Ends{Fore.RESET}")

    if my_args.fa:
        print("NOTE: Activating Factor Analysis function!")
        # Factor Analysis default parameters
        dimension = 2

        # Factor Analysis parameters
        if my_args.td:
            dimension = my_args.td

        print("\n" + "Factor Analysis Parameters:")
        print(f"Dimension:{dimension}")

        fa_result = func.factor_analysis(my_args.input, n_dimension=dimension)

        fa_result_path = result_path + os.sep + output_prefix + "fa.csv"
        fa_result.to_csv(fa_result_path)
        print("NOTE: Factor Analysis Completed!")

        # PLOTTING
        if my_args.pc:
            # Plot Scatter
            print(f"NOTE: Activating PCA Plotting function!")
            fa_result["Virus"] = fa_result.index
            scatter_label = False
            figure_format = "pdf"
            if my_args.format:
                figure_format = my_args.format

            figure_name = result_path + os.sep + output_prefix + "FA-Plot" + "." + figure_format
            dpi = 300
            transparent = False
            figure_width = 8
            figure_height = 6
            palette = "RdYlBu"

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

            func.cluster_figure(fa_result,
                                point_label=scatter_label,
                                fig_name=figure_name,
                                fig_dpi=dpi,
                                tp_tag=transparent,
                                fig_width=figure_width,
                                fig_height=figure_height,
                                color_tag=palette)

            print("\n" + "NOTE: Plotting Completed!")

        if my_args.pc2:
            print(f"{Fore.LIGHTYELLOW_EX}NOTE: Testing Part Start{Fore.RESET}")
            figure_name = result_path + os.sep + output_prefix + "FA-Plot" + ".pdf"
            func.clustering_plot2(fa_result, fig_name=figure_name)

            print(f"{Fore.LIGHTYELLOW_EX}NOTE: Testing Part Ends{Fore.RESET}")

    if my_args.mds:
        print("NOTE: Activating MDS function!")
        # Factor Analysis default parameters
        dimension = 2

        # Factor Analysis parameters
        if my_args.td:
            dimension = my_args.td

        print("\n" + "MDS Parameters:")
        print(f"Dimension:{dimension}")

        mds_result = func.mds(my_args.input, n_dimension=dimension)

        mds_result_path = result_path + os.sep + output_prefix + "mds.csv"
        mds_result.to_csv(mds_result_path)
        print("NOTE: MDS Completed!")

        # PLOTTING
        if my_args.pc:
            # Plot Scatter
            print(f"NOTE: Activating PCA Plotting function!")
            mds_result["Virus"] = mds_result.index
            scatter_label = False
            figure_format = "pdf"
            if my_args.format:
                figure_format = my_args.format

            figure_name = result_path + os.sep + output_prefix + "MDS-Plot" + "." + figure_format
            dpi = 300
            transparent = False
            figure_width = 8
            figure_height = 6
            palette = "RdYlBu"

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

            func.cluster_figure(mds_result,
                                point_label=scatter_label,
                                fig_name=figure_name,
                                fig_dpi=dpi,
                                tp_tag=transparent,
                                fig_width=figure_width,
                                fig_height=figure_height,
                                color_tag=palette)

            print("\n" + "NOTE: Plotting Completed!")

        if my_args.pc2:
            print(f"{Fore.LIGHTYELLOW_EX}NOTE: Testing Part Start{Fore.RESET}")
            figure_name = result_path + os.sep + output_prefix + "mds-Plot" + ".pdf"
            func.clustering_plot2(mds_result, fig_name=figure_name)

            print(f"{Fore.LIGHTYELLOW_EX}NOTE: Testing Part Ends{Fore.RESET}")

    if my_args.pts:
        result_dir = os.path.dirname(my_args.input)
        figure_name = os.path.join(result_dir, "3d-scatter.pdf")

        file_type = func.file_type_judge(my_args.input)
        # determine the file type
        if file_type == "CSV":
            raw_df = pd.read_csv(my_args.input, header=0, index_col=0)
        elif file_type == "XLSX":
            raw_df = pd.read_excel(my_args.input, header=0)
        else:
            print(Fore.RED + "ERROR! The file type:%s is not supported right now!" % file_type)
            sys.exit()

        func.plot_3d_scatter(dataframe=raw_df, fig_name=figure_name)

    if my_args.dbscan:
        print("NOTE: Activating DBSCAN(Density-Based Spatial Clustering of Applications with Noise)")
        dbscan_result_path = result_path + os.sep + output_prefix + "DBSCAN.csv"

        file_type = func.file_type_judge(my_args.input)
        # determine the file type
        if file_type == "CSV":
            raw_df = pd.read_csv(my_args.input, header=0, index_col=0)
        elif file_type == "XLSX":
            raw_df = pd.read_excel(my_args.input, header=0)
        else:
            print(Fore.RED + "ERROR! The file type:%s is not supported right now!" % file_type)
            sys.exit()

        eps = 0.5
        minsap = 6
        if my_args.eps:
            eps = my_args.eps
        if my_args.minsap:
            minsap = my_args.minsap

        dbscan_res = func.dbscan(dimension_mat=raw_df, result_path=dbscan_result_path, epsilon=eps, minPts=minsap)
        dbscan_res.to_csv(dbscan_result_path)
        print("NOTE: DBSCAN is done!")

    if my_args.optics:
        print("NOTE: Activating Ordering Points To Identify the Clustering Structure (OPTICS)")
        optics_result_path = result_path + os.sep + output_prefix + "OPTICS.csv"

        file_type = func.file_type_judge(my_args.input)
        # determine the file type
        if file_type == "CSV":
            raw_df = pd.read_csv(my_args.input, header=0, index_col=0)
        elif file_type == "XLSX":
            raw_df = pd.read_excel(my_args.input, header=0)
        else:
            print(Fore.RED + "ERROR! The file type:%s is not supported right now!" % file_type)
            sys.exit()

        func.optics(dimension_mat=raw_df, result_path=optics_result_path)
        print("NOTE: OPTICS is done!")

    if my_args.kmean:
        print("NOTE: Activating K-Mean")
        kmean_result_path = result_path + os.sep + output_prefix + "KMean.csv"
        n_cluster = 3
        file_type = func.file_type_judge(my_args.input)
        # determine the file type
        if file_type == "CSV":
            raw_df = pd.read_csv(my_args.input, header=0, index_col=0)
        elif file_type == "XLSX":
            raw_df = pd.read_excel(my_args.input, header=0)
        else:
            print(Fore.RED + "ERROR! The file type:%s is not supported right now!" % file_type)
            sys.exit()

        if my_args.nc:
            n_cluster = my_args.nc

        func.kmean(dimension_mat=raw_df, result_path=kmean_result_path, cluster=n_cluster)
        print("NOTE: K-Mean is done!")

    if my_args.agglo:
        print("NOTE: Activating Agglomerative!")
        agglo_result_path = result_path + os.sep + output_prefix + "Agglomerative.csv"
        n_cluster = 3
        file_type = func.file_type_judge(my_args.input)
        # determine the file type
        if file_type == "CSV":
            raw_df = pd.read_csv(my_args.input, header=0, index_col=0)
        elif file_type == "XLSX":
            raw_df = pd.read_excel(my_args.input, header=0)
        else:
            print(Fore.RED + "ERROR! The file type:%s is not supported right now!" % file_type)
            sys.exit()

        if my_args.nc:
            n_cluster = my_args.nc

        func.agglomerative(dimension_mat=raw_df, result_path=agglo_result_path, cluster=n_cluster)
        print("NOTE: Agglomerative is done!")

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
