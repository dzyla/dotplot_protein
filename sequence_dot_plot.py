import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm


def check_arg():
    conc = sys.argv[1:]
    if '-f' and '-s1' and '-s2' and '-w' and '-t' not in conc:
        text = 'DotPlot plotting script use it like this: \n\nsequence_dotplot.py -f fasta.file -s1 1 -s2 0 -w 10 -t 50\n\n-f fasta file with fasta sequences\n-s1 / -s2 sequence number in the file (if the same sequence use the same number)\n-w windo size of calulated aligement\n-t treshold of measured values in % of 100 corresponding aminoacids in sequence'
        print(text)
        exit(0)
    if len(conc) == 0:
        print(text)
        exit(0)
    if '-f' in conc:
        f = conc.index('-f')
        try:
            file = conc[f + 1]
        except IndexError:
            print('\nFasta file not privided.')
            exit(0)
    if '-s1' in conc:
        s1 = conc.index('-s1')
        try:
            seq_num1 = int(conc[s1 + 1])
        except IndexError:
            print('argument missing!')
            exit(0)
    if '-s2' in conc:
        s2 = conc.index('-s2')
        try:
            seq_num2 = int(conc[s2 + 1])
        except IndexError:
            print('argument missing!')
            exit(0)
    if '-w' in conc:
        w = conc.index('-w')
        try:
            ww = float(conc[w + 1])
        except IndexError:
            print('argument missing!')
            exit(0)
    if '-t' in conc:
        t = conc.index('-t')
        try:
            tres = float(conc[t + 1])
        except IndexError:
            print('argument missing!')
            exit(0)
        return file, seq_num1, seq_num2, ww, tres


def blosum62():
    blosum62_matrix = pd.read_csv('blosum62.txt')
    blosum62_matrix = pd.DataFrame.from_records(blosum62_matrix)
    blosum62_matrix.set_index('Name', inplace=True)
    return blosum62_matrix


def align_value(blos_matrix, aa1, aa2):
    aa1 = str(aa1)
    aa2 = str(aa2)
    return blos_matrix.loc[aa1, aa2]


def get_sequence(file, seq_num1, seq_num2):
    count = []
    for seq, rec in enumerate(SeqIO.parse(file, 'fasta')):
        if seq == seq_num1:
            sequence1 = rec.seq
            name1 = rec.description
        if seq == seq_num2:
            sequence2 = rec.seq
            name2 = rec.description
        count.append(seq)
    if seq_num2 not in count or seq_num1 not in count:
        print('number of sequence not in the file. Exiting...')
        exit(1)
    return sequence1, sequence2, name1, name2

def dot_plot_range(seq1, seq2, blos_matrix, ww, threshold):
    time1 = time.time()
    matrix = np.zeros((len(seq1), len(seq2)), np.float)

    # ww - the windows size, treshold in %
    ww1 = -ww / 2 + 1
    ww2 = ww / 2 + 1
    start_seq1 = -ww1
    end_seq1 = len(seq1) - ww2 + 1
    start_seq2 = -ww1
    end_seq2 = len(seq2) - ww2 + 1

    # count = 0
    # old code without tqdm
    # for i,a in enumerate(np.arange(start_seq1,end_seq1,1)):
    #     if end_seq1 > 100:
    #         if int((i/(end_seq1-start_seq1)))%5 == 0 and count > 0:
    #             print(str(int(count/end_seq1))+' / 100%')

    #        print(int(i/(end_seq1-start_seq1)*100))

    for a in tqdm(np.arange(start_seq1, end_seq1, 1)):
        for b in np.arange(start_seq2, end_seq2):
            seq = 0
            total = 0
            for k in np.arange(ww1, ww2, 1):
                seq += align_value(blos_matrix, (seq1[int(a + k)]), (seq2[int(b + k)]))
                total += align_value(blos_matrix, (seq1[int(a + k)]), (seq1[int(a + k)]))
            if threshold > seq / total * 100:
                seq = 0
            matrix[int(a), int(b)] = seq/total*100
    time2 = time.time()
    print('done in: ' + str(round(time2 - time1, 2)) + ' s')
    return matrix


def plot(matrix, name1, name2,ww,tres):
    plt.figure(figsize=(10, 10))
    plt.imshow(matrix, cmap='hot', interpolation='nearest')
    plt.xlabel(name1)
    plt.ylabel(name2)
    plt.title('Dotplot of:\n\n' + name1 + '\n' + name2+'\n\n with window size: '+str(ww)+' and threshold '+str(tres))
    legend = plt.colorbar(fraction=0.046, pad=0.04)
    legend.set_label('% of identity')
    plt.savefig('dotplot.png')
    plt.show()


# load arguments from command line
file, seq1_num, seq2_num, ww, tres = check_arg()

# do the magic
seq1, seq2, name1, name2 = get_sequence(file, seq1_num, seq2_num)
blosum_matrix = blosum62()
plot(dot_plot_range(seq1, seq2, blosum_matrix, ww, tres), name1, name2,ww,tres)
