# -*- coding: utf-8 -*-
# Development Time: 2024-02-27 14:06:38
# Developer: XiaoYang

# DNA codon translation table
translate_table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    '---': '-',
}

aa_paltte = {"A": "#F08080", "I": "#FF7F50", "L": "#E9967A",
             "M": "#F0E68C", "V": "#9ACD32", "F": "#90EE90",
             "W": "#66CDAA", "Y": "#3CB371", "N": "#20B2AA",
             "C": "#40E0D0", "S": "#B0E0E6", "T": "#4682B4",
             "Q": "#87CEEB", "D": "#9370DB", "E": "#7B68EE",
             "G": "#DDA0DD", "P": "#DB7093", "R": "#FF69B4",
             "H": "#FAFAD2", "K": "#FFDAB9", "-": "#A9A9A9",
             "*": "#A52A2A", "?": "#808080"}

# 20 aa, no illegal icon
# 21 aa, including '?'
# 21 aa, including '-'
# 22 aa, including '?' and '-'
amino_acids_list = {"type1": ['A', 'R', 'N', 'D', 'C',
                              'Q', 'E', 'G', 'H', 'I',
                              'L', 'K', 'M', 'F', 'P',
                              'S', 'T', 'W', 'Y', 'V'],
                    "type2": ['A', 'R', 'N', 'D', 'C',
                              'Q', 'E', 'G', 'H', 'I',
                              'L', 'K', 'M', 'F', 'P',
                              'S', 'T', 'W', 'Y', 'V', "?"],
                    "type3": ['A', 'R', 'N', 'D', 'C',
                              'Q', 'E', 'G', 'H', 'I',
                              'L', 'K', 'M', 'F', 'P',
                              'S', 'T', 'W', 'Y', 'V', "-"],
                    "type4": ['A', 'R', 'N', 'D', 'C',
                              'Q', 'E', 'G', 'H', 'I',
                              'L', 'K', 'M', 'F', 'P',
                              'S', 'T', 'W', 'Y', 'V', "?", "-"]}

# Tag dictionary
tag_amino_list = {"type1": 20,
                  "type2": 21,
                  "type3": 21,
                  "type4": 22}

# palettes
palettes = ["Accent", "Blues/Blues_r/Blues_d", "BrBG/BrBg_r", "BuGn/BuBn_r", "BuPu/BuPu_r", "CMRmap",
            "Dark2/Dark2_r", "GnBu/GnBu_r", "Greens/Green_r", "RdYlBu", "Greys/Greys_r", "OrRd/OrRd_r",
            "Oranges/Oranges_r/Oranges_d", "PRGn/PRGn_r", "Paired/Paired_r", "Pastel1/Pastel1_r", "Pastel2/Pastel2_r",
            "PiYG/PiYG_r", "PuBu/PuBu_r", "PuOr/PuOr_r", "PuRd/PuRd_r", "Purples/Purplrs_r/Purples_d", "RdBu/RdBu_r",
            "RdGy/RdGy_r", "RdPu/RdPu_r", "RdYlBu/RdYlBu_r", "RdYlGn/RdYlGn_r", "Reds/Reds_r", "Set1/Set1_r",
            "Set2/Set2_r", "Set3/Set3_r", "Spectral", "Wistia", "YlGn/YlGn_r", "YlGnBu/YlGnBu_r", "YlOrBr/YlOrBr_r",
            "afmhot/afmhot_r", "autumn/autumn_r", "binary/binary_r", "bone/bone_r", "brg", "bwr", "bwr_r",
            "cividis/cividis_r", "cool/cool_r", "coolwarm", "copper/copper_r", "crest/crest_r", "cubehelix/cubehelix_r",
            "flag/flare_r", "gist_earth/gist_earth_r", "gist_gray/gist_gray_r", "gist_heat/gist_heat_r",
            "gist_ncar/gist_ncar_r", "gist_rainbow/gist_rainbow_r", "gist_stern/gist_stern_r", "gist_yarg/gist_yarg_r",
            "gnuplot/gnuplot_r", "gnuplot2/gnuplot2_r", "hot/hot_r", "hsv/hsv_r", "icefire/icefire_r",
            "inferno/inferno_r", "magma/magma_r", "mako/mako_r", "nipy_spectral/nipy_spectral_r", "ocean/ocean_r",
            "pink/pink_r", "plasma/plasma_r", "prism", "rainbow", "rocket", "seismic/seismic_r", "spring/spring_r",
            "summer/summer_r", "tab10/tab10_r", "tab20/tab20_r", "tab20b/tab20b_r", "tab20c/tab20c_r",
            "terrain/terrain_r", "twilight/twilight_r", "twilight_shifted/twilight_shifted_r", "viridis/viridis_r",
            "vlag/vlag_r", "winter/winter_r"]
