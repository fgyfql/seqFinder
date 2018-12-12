#!/usr/bin/env python
import os, sys
def main():
    options = parse_command_line()
    ranked_seq_score_list = do_align(options)

    f = open(options.output_file, 'w')
    for each in ranked_seq_score_list:
        if each[1][0] == 0:
            continue
        f.write("****************\n")
        f.write("aligned protein name: %s\n" %each[0].strip('>'))
        f.write("fitting score: %f\n" % each[1][0])
        f.write("aligned sequence:\n")
        f.write('\t#column 1: residue number of the C-alpha model\n')
        f.write('\t#column 2: residue number of the aligned protein sequence\n')
        f.write('\t#column 3: residue name of the aligned protein sequence\n')
        for eac in each[1][1]:
            f.write('\t')
            for ea in eac:
                f.write(str(ea)+'\t')
            f.write('\n')
        f.write('\n\n')
    f.close()

def do_align(options):
    annotation_file = options.annotation_file
    annotation_data = annotation_read(annotation_file)
    annotation_pattern = annotation_pattern_extract(annotation_data)
    sequences_dict = sequences_read(options.fasta_file)
    windows = window_maker(annotation_pattern)
    seq_score_dict = score_calculate(windows, sequences_dict)
    ranked_seq_score_list = score_rank(seq_score_dict)
    return ranked_seq_score_list
        
def score_rank(seq_score_dict):
    seq_score_list = []
    for each in seq_score_dict:
        seq_score_list.append((each, seq_score_dict[each]))
    ranked_seq_score_list = [seq_score_list[0]]
    for each in seq_score_list:
        if each in ranked_seq_score_list:
            continue
        for i in range(len(ranked_seq_score_list)):
            if i == len(ranked_seq_score_list) - 1 and each[1][0] > ranked_seq_score_list[i][1][0]:
                ranked_seq_score_list.append(each)
                break
            if each[1][0] <= ranked_seq_score_list[i][1][0]:
                ranked_seq_score_list.insert(i, each)
                break
    return ranked_seq_score_list

def score_calculate(windows, sequences_dict):
    score_dict = {}
    for each in sequences_dict:
        score = find_max_score(windows, sequences_dict[each])
        score_dict[each] = score
    return score_dict

def find_max_score(windows, sequence):
    max_score = 0
    aligned_seq_list = []
    for window in windows:
        if len(window[0]) > len(sequence)-1:
            continue
        for i in range(1, len(sequence)-len(window[0])+1):
            count_list = []
            aligned_list = []
            for j in range(len(window[0])):
                if window[0][j] != 'X':
                    aligned_list.append([window[1][j], i+j+1, sequence[i+j]])
                    dif = window[0][j] - amino_acid_size_num(sequence[i+j])
                    if dif >=2:
                        count_list = []
                        break
                    elif dif == 1:
                        count_list.append(0)
                    elif dif <1:
                        count_list.append(6+dif)
            if count_list == []:
                continue
            score = sum(count_list)/float(len(count_list))/6.0
            if score > max_score:
                max_score = score
                aligned_seq_list = aligned_list
    return [max_score, aligned_seq_list]

def amino_acid_size_num(aa):
    amino_acid_size_grouped_dict = {'G':0,'A':1,'S':1,'C':1,'P':1,'V':2,'T':2,'I':3,'L':3,'D':3,'N':3,'E':4,'Q':4,'M':4,'H':4,'K':4,'R':5,'F':5,'Y':5,'W':6}
    return amino_acid_size_grouped_dict[aa]

def annotation_pattern_extract(annotation_data):
    pattern = [[], []]
    for i in range(len(annotation_data)):
        if i == 0:
            if annotation_data[i][1] == 'X':
                pattern[0].append('X')
            else:
                pattern[0].append(amino_acid_size_num(annotation_data[i][1]))
            pattern[1].append(annotation_data[i][0])
            continue
        if i > 0:
            gap = int(annotation_data[i][0])-int(annotation_data[i-1][0])-1
            if gap>0:
                pattern[0].append('g-%d' %gap)
                pattern[1].append('NA')
        if annotation_data[i][1] == 'X':
            pattern[0].append('X')
        else:
            pattern[0].append(amino_acid_size_num(annotation_data[i][1]))
        pattern[1].append(annotation_data[i][0])
    return pattern

def window_maker(annotation_pattern):
    windows = [[[],[]]]
    for k in range(len(annotation_pattern[0])):
        if annotation_pattern[0][k] == 'X' or type(annotation_pattern[0][k]) == int:
            for i in range(len(windows)):
                windows[i][0].append(annotation_pattern[0][k])
                windows[i][1].append(annotation_pattern[1][k])

        elif 'g-' in annotation_pattern[0][k]:
            gap = int(annotation_pattern[0][k].split('-')[-1])
            delta_gap = delta(gap)
            if gap > 6:
                half_delta_gap = int(round(delta_gap/2.0))
            else:
                half_delta_gap = delta_gap/2
            tmp = []
            for i in range(-half_delta_gap, delta_gap+1):
                tmp2 = ['X']*(gap + i)
                tmp3 = ['NA']*(gap + i)
                tmp.append([tmp2,tmp3])
            windows_tmp = []
            for i in range(len(windows)):
                for j in range(len(tmp)):
                    e0 = windows[i][0] + tmp[j][0]
                    e1 = windows[i][1] + tmp[j][1]
                    windows_tmp.append([e0, e1])
            windows = windows_tmp
    return windows

def sequences_read(fn):
    d = {}
    f = open(fn, 'r')
    line = f.readlines()
    for i in range(len(line)):
        if '>' in line[i]:
            key = line[i].strip()
            d[key] = line[i+1].strip()
    return d

def delta(gap):
    if gap == 0:
        deltaDif = 0
    if gap <= 2 and gap > 0:
        deltaDif = 1
    if gap > 2 and gap <= 6:
        deltaDif = 2
    if gap >6 and gap <= 15:
        deltaDif = 3
    if gap > 15:
        deltaDif = int(round(0.2*dif, 0))
    return deltaDif

def annotation_read(fn):
    f = open(fn, 'r')
    line = f.readlines()
    data = []
    for each in line:
        data.append(each.split())
    data = data_sort(data)
    return data

def data_sort(data):
    data_sorted = []
    data_tmp_dict = {}
    for each in data:
        data_tmp_dict[int(each[0])] = each
    keylist = data_tmp_dict.keys()
    keylist.sort()
    for key in keylist: 
        data_sorted.append(data_tmp_dict[key])
    return data_sorted

def parse_command_line():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--annotation_file', metavar="<file>", type=str, help='annotation file.')
    parser.add_argument('--fasta_file', metavar="<file>", type=str, help='a fasta file containg all the sequences.')
    parser.add_argument('--output_file', metavar="<file>", type=str, help='output file.')
    
    if len(sys.argv)<4:
        parser.print_help()
        sys.exit(-1)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()
