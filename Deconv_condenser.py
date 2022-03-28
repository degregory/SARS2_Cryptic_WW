#!/bin/env python3

# Writen by Devon Gregory and Rose Kantor
# University of Missouri and University of California - Berkley
# Distributed under GNU GENERAL PUBLIC LICENSE v3
# last edit on 20220327

import os
import sys
import pandas as pd
import numpy as np
from plotnine import *
import re

omicron_mutations = [ ## Original From Rose
    'T19I',
    'LPPA24-27S---',
    'A67V',
    'IHV68-70I--',
    'T95I',
    'GVYY142-145D---',
    'G142D',
    'NL211-212I-',
    'V213G',
    '215EPE',
    'G339D',
    'S371L',
    'S373P',
    'S375F',
    'K417N',
    'N440K',
    'G446S',
    'S477N',
    'T478K',
    'E484A',
    'Q493R',
    'G496S',
    'Q498R',
    'N501Y',
    'Y505H',
    'T547K'
]
omicron_positions = [
    '19',
    '24',
    '25',
    '26',
    '27',
    '67',
    '68',
    '69',
    '70',
    '95',
    '142',
    '143',
    '144',
    '145',
    '211',
    '212',
    '213',
    '215',
    '339',
    '371',
    '373',
    '375',
    '417',
    '440',
    '446',
    '477',
    '478',
    '484',
    '493',
    '496',
    '498',
    '501',
    '505',
    '547'
    ]

def mut_pos_stringer(muts):
    mut_pos_string = ''
    for mut in muts.split(' '):
        mut_pos_string += mut.split('(')[0].strip('ATCGN') + ' '
    return(mut_pos_string)

def load_deconv(file_name): ## From Rose
    '''ingest file and change column names, add seq IDs,
    filter out sequences with fewer than 5 mutations, melt to long format
    returns the filtered raw file with frequencies and the long format file'''
    # load data, rename columns, drop the rows called "Covariant" and "Reference"
    df_raw = pd.read_csv(file_name, sep='\t')
    df_raw = df_raw.rename(columns={'Unnamed: 0':'mutations'})
    # df_raw = df_raw[~df_raw.mutations.isin(['Covariant', 'Reference'])]

    # make column with count of mutations and then filter out seqs with < 5 mutations
    # df_raw['mutations_count'] = df_raw.mutations.apply(lambda x: len(x.split(' ')))
    # df_raw = df_raw[df_raw.mutations_count >= 5]
    # df_raw = df_raw.drop(columns='mutations_count') # now drop this col

    # drop columns that are all NA
    df_raw = df_raw.dropna('columns', how='all')


    df_raw['mut_poses'] = df_raw.mutations.apply(lambda x: mut_pos_stringer(x))
    df_raw = df_raw.sort_values(by='mut_poses').reset_index().drop(columns=['index'])
    df_raw = df_raw.drop(columns=['mut_poses'])

    # pat_E484 = re.compile(r'\(E484(\D)\)')
    # pat_N440 = re.compile(r'\(N440(\D)\)')
    # df_raw['E484'] = df_raw.mutations.apply(lambda x: bool(pat_E484.search(x)))
    # df_raw['N440'] = df_raw.mutations.apply(lambda x: bool(pat_N440.search(x)))

    # df_raw = df_raw.sort_values(by='E484').reset_index().drop(columns=['index'])
    # df_raw = df_raw.drop(columns=['E484'])
    # df_raw = df_raw.sort_values(by='N440').reset_index().drop(columns=['index'])
    # df_raw = df_raw.drop(columns=['N440'])
    # print(df_raw)

    df_raw.fillna(0)
    first_date_dict = {}

    for row in df_raw.itertuples():
        for i in range(2, len(row)):
            if row[i] > 0:
                first_date_dict[row[1]] = i
                break

    df_raw['first_date'] = df_raw.mutations.apply(lambda x: first_date_dict[x])

    df_raw = df_raw.sort_values(by='first_date').reset_index().drop(columns=['index'])
    df_raw = df_raw.drop(columns=['first_date'])
    df_raw = df_raw[::-1]

    df_raw.insert(0,'seq_id', range(0,len(df_raw)))
    df = df_raw.melt(id_vars=['mutations', 'seq_id'], var_name='sample_id', value_name='abundance').dropna()
    df['sample'] = df.sample_id.apply(lambda x: '-'.join(x.split('-')[1:]))
    # df['Date'] = df.sample_id.apply(lambda x: x.split('_')[1])
    # df = df.drop(columns=['sample_id'])


    return df_raw, df

def condensePMs(seq_line_dict, min_samp):
    PMs = {}
    samps = seq_line_dict['samps']
    passed_lines = []
    for line in seq_line_dict['lines']:
        split_line = line.strip("\n\r").split("\t")
        seq_len = len(split_line[0].split(" "))
        if seq_len > 3:
            discard = 0
            if seq_len < 6:
                if ("N501Y" in split_line[0] and "A570D" in split_line[0]) or ("L452R" in split_line[0] and "T478K" in split_line[0]):
                    discard = 1
                if ("203-208Del" in split_line[0] and "429-431Del" in split_line[0]) or ("425A(G142D)" in split_line[0] and "467-472Del" in split_line[0]):
                    discard = 1
                if ("K417" in split_line[0] and "E484K" in split_line[0] and "N501Y" in split_line[0]): # or ("425A(G142D)" in split_line[0] and "467-472Del" in split_line[0]):
                    discard = 1

            Omi_matches = 0
            for mut in omicron_mutations:
                if mut in split_line[0]:
                    Omi_matches += 1
            if discard == 0 and not Omi_matches > (seq_len * .7):
                passed_lines.append(split_line)
                samp_hits = []
                total = 0
                maxim = 0
                for i in range(1, len(samps)):
                    if not split_line[i] == "":
                        samp_hits.append(samps[i])
                        total += float(split_line[i])
                        if float(split_line[i]) > maxim:
                            maxim = float(split_line[i])

                for PM in split_line[0].split(" "):
                    try:
                        PMs[PM]['total'] += total
                        if maxim > PMs[PM]["max"]:
                            PMs[PM]["max"] = maxim
                    except:
                        PMs[PM] = {'total' : total,
                                    'max' : maxim,
                                    'samps' : []}
                    for hit in samp_hits:
                        if not hit in PMs[PM]['samps']:
                            PMs[PM]['samps'].append(hit)

    condensed_seq = {}
    new_samps = ['']
    for line in passed_lines:
        new_seq = []
        for PM in line[0].split(" "):   
            if (len(PMs[PM]['samps']) >= min_samp and PMs[PM]['max'] > (.005) and PMs[PM]['total'] > (.005)) or 'del' in PM.lower() or '-)' in PM:
                if '(fs)' in PM or '*' in PM:
                    continue
                if '(' in PM and not '-' in PM:
                    if PM.split('(')[1][0] == PM.split('(')[1][-2]:
                        # print(PM)
                        continue
                else:
                    if PM[0] == PM[-1]:
                        print(PM)
                        continue
                    
                position = ''
                for c in PM.split('(')[1].split('-')[0]:
                    if c.isdigit():
                        position += c
                if int(position) <= 410 and int(position) >= 403:
                    continue
                if int(position) <= 587 and int(position) >= 580:
                    continue
                # try:
                    # if (not (PM.endswith("*)") and PM.split('(')[1][0] == PM[-2])) and not '-' in PM.split('(')[1]:
                        # new_seq.append(PM)
                # except:
                new_seq.append(PM)

        new_seq = " ".join(new_seq)
        # if not new_seq:
            # new_seq = 'Reference'
        try:
            condensed_seq[new_seq]['total'] += 1
        except:
            condensed_seq[new_seq] = {'total' : 1,
                                        'samps' : {}}

        for i in range(1, len(samps)):
            if not line[i] == "":
                try:
                    condensed_seq[new_seq]['samps'][samps[i]] += float(line[i])
                except:
                    condensed_seq[new_seq]['samps'][samps[i]] = float(line[i])
    new_seq_lines = []
    for seq in sorted(list(condensed_seq.keys())):
        if (len(condensed_seq[seq]['samps']) > 1 or sum(condensed_seq[seq]['samps'].values()) >= .02) and seq:
            new_seq_line = seq
            for i in range(1, len(samps)):
                try:
                    new_seq_line += f"\t{condensed_seq[seq]['samps'][samps[i]]}"
                except:
                    new_seq_line += "\t"
            new_seq_line += "\n"
            new_seq_lines.append(new_seq_line)
    return(new_seq_lines)

for file in os.listdir(os.getcwd()):
    ## break
    if file.endswith('Covar_Deconv.tsv') and not 'Condensed_' in file:
        samps = []
        lines = []

        collect_fh = open(file,"r")
        for line in collect_fh:
            split_line = line.strip("\n\r").split("\t")
            if split_line[0] == "":
                samps = split_line
            elif not split_line[0]== "Covariant":
                lines.append(line)
        collect_fh.close()



        for i in range(1, 4):
            precount = len(lines)
            postcount = 0
            min_samp = i
            inf_loop_sheild = 10
            precond_dict = {'lines' : lines,
                        'samps' : samps}
            while precount > postcount:
                precount = len(precond_dict['lines'])
                cond_lines = condensePMs(precond_dict, min_samp)
                postcount = len(cond_lines)
                precond_dict['lines'] = cond_lines
                inf_loop_sheild -= 1
                if inf_loop_sheild < 1:
                    break

            if cond_lines:
                condensed_fh = open("Condensed_min"+str(min_samp)+"_"+file,"w")
                condensed_fh.write("\t".join(samps))
                condensed_fh.write("\n")
                for line in cond_lines:
                    condensed_fh.write(line)
                condensed_fh.close()

for file in os.listdir(os.getcwd()):
    ## break
    if file.endswith('Covar_Deconv.tsv'):
        if file.startswith('Condensed_'):

            ####
            ####  Original plotting from Rose
            ####
            df_raw, df = load_deconv(file)


            df['PMs'] = df.mutations.apply(lambda x: x.split(' '))
            df = df.drop(columns=['mutations'])



            pat_SCP = re.compile(r'\(((\D)(\d+)(\w))\)')
            pat_ins = re.compile(r'insert\w+\(((\d+)(\w+))\)')
            pat_ins1 = re.compile(r'insert\w+\(((\w)(\d+)(\w+))\)')
            pat_del = re.compile(r'\(((\w)(\d+))-\)')
            pat_MCP = re.compile(r'\(((\D+)(\d+)-(\d+)(\D+))\)')
            new_entries = []
            for row in df.itertuples():
                for pm in row.PMs:
                    Omi = 'Non-Omicron'
                    if 'reference' in pm.lower():
                        print('ref fround in ' + pm)

                    # if ':' in pm:
                        # print(pm)
                        # continue
                    if pat_ins.search(pm):
                        exp = pat_ins.search(pm).groups()
                        position = exp[1]
                        wt = '-'
                        wt_pos = '-'+exp[1]
                        mut = exp[2]
                        if exp[0] in omicron_mutations:
                            Omi = 'Omicron Mutation'
                        elif exp[1] in omicron_positions:
                            Omi = 'Omicron Position'
                        if len(mut) > 1:
                            i = 0
                            for AA in mut:
                                new_entries.append([row.seq_id,  wt_pos+'.'+str(i), int(position), wt, mut[i], Omi])
                                i += 1
                            continue
                    elif pat_ins1.search(pm):
                        exp = pat_ins1.search(pm).groups()
                        position = exp[2]
                        wt = exp[1]
                        mut = exp[3]
                        if not wt == mut[0]:
                            if wt+position+mut[0] in omicron_mutations:
                                Omi = 'Omicron Mutation'
                            elif exp[2] in omicron_positions:
                                Omi = 'Omicron Position'
                            new_entries.append([row.seq_id, exp[1]+exp[2], int(position), wt, mut[0], Omi])
                        Omi = 'Non-Omicron'
                        wt = '-'
                        mut = mut[1:]

                        if position+mut in omicron_mutations:
                            Omi = 'Omicron Mutation'
                        elif position in omicron_positions:
                            Omi = 'Omicron Position'
                        i = 1
                        for AA in mut:
                            newposition = str(int(position)+1) +'.'+str(i)
                            new_entries.append([row.seq_id, '-'+newposition, float(newposition), wt, AA, Omi])
                            i += 1
                        continue

                    elif pat_del.search(pm):
                        exp = pat_del.search(pm).groups()
                        position = exp[2]
                        wt = exp[1]
                        wt_pos = exp[1]+exp[2]
                        mut = 'Δ'
                        if exp[0] in omicron_mutations:
                            Omi = 'Omicron Mutation'
                        elif exp[2] in omicron_positions:
                            Omi = 'Omicron Position'
                    elif pat_SCP.search(pm):
                        exp = pat_SCP.search(pm).groups()
                        position = exp[2]
                        if int(position) <= 410 and int(position) >= 403:
                            continue
                        if int(position) <= 587 and int(position) >= 580:
                            continue
                        wt = exp[1]
                        wt_pos = exp[1]+exp[2]
                        mut = exp[3]
                        if wt == mut:
                            continue
                        if exp[0] in omicron_mutations:
                            Omi = 'Omicron Mutation'
                        elif exp[2] in omicron_positions:
                            Omi = 'Omicron Position'
                    elif pat_MCP.search(pm):
                        exp = pat_MCP.search(pm).groups()
                        position = exp[2]
                        wt = exp[1]
                        wt_pos = 'x'
                        mut = exp[4]
                        # if len(mut) == len(wt):
                        for i in range(0, len(wt)):
                            try:
                                if wt[i] == mut[i]:
                                    continue
                                Omi = 'Non-Omicron'
                                if (wt[i]+str(int(position)+i)+mut[i]) in omicron_mutations:
                                    Omi = 'Omicron Mutation'
                                elif str(int(position)+i) in omicron_positions:
                                    Omi = 'Omicron Position'
                                new_entries.append([row.seq_id,  wt[i]+str(int(position)+i), int(position)+i, wt[i], mut[i].replace('-', 'Δ'), Omi])
                            except:
                                if Omi == 'Non-Omicron':
                                    if (wt[i]+str(int(position)+i)+'-') in omicron_mutations:
                                        Omi = 'Omicron Mutation'
                                    elif str(int(position)+i) in omicron_positions:
                                        Omi = 'Omicron Position'
                                new_entries.append([row.seq_id,  wt[i]+str(int(position)+i), int(position)+i, wt[i], 'Δ', Omi])
                                #print(pm)
                            # print(new_entries[-1])
                            Omi = 'Non-Omicron'
                        if len(mut) > len(wt):
                            mut = mut[i+1:]
                            wt = '-'
                            position = int(position) + i + 1
                            if str(position)+mut in omicron_mutations:
                                Omi = 'Omicron Mutation'
                            elif str(position) in omicron_positions:
                                Omi = 'Omicron Position'
                            i = 1
                            for AA in mut:
                                newposition = str(int(position)+1) +'.'+str(i)
                                new_entries.append([row.seq_id,  '-'+newposition, float(newposition), wt, AA, Omi])
                                i += 1

                        continue
                    else:
                        print(pm)
                        continue
                    new_entries.append([row.seq_id,  wt_pos, int(position), wt, mut, Omi])
            #try:
            df_long = pd.DataFrame.from_records(new_entries, columns = ['seq_id', 'wt_pos', 'position', 'wildtype_aa', 'mutation_aa', 'Omicron Residues'])

            # order the data by wildtype position for plotting nicely, From Rose
            sorting = df_long[['position', 'wt_pos']].copy()
            sorting.position = pd.to_numeric(sorting.position)
            ordered_list = list(sorting.sort_values('position')['wt_pos'].unique())
            df_long.wt_pos = pd.Categorical(df_long.wt_pos, ordered=True, categories=ordered_list)

            seq_num = df_raw.shape[0]
            samp_num = df_raw.shape[1] - 2

            pm_num = df_long.wt_pos.nunique()
            colors = ['#111111', '#d55e00', '#009e73']
            xsize = pm_num /2.75
            ysize = seq_num /2.75

            # while (xsize > 25 or ysize > 25 ):
                # xsize = xsize / 2
                # ysize = ysize / 2

            fig = (ggplot(df_long, aes(x='wt_pos', y='seq_id', color='Omicron Residues', limitsize=False))+
             geom_tile(aes(width=.9, height=.9), size=1, linetype='solid', fill='white')+ #
             geom_text(aes(label='mutation_aa'), color='black')+
             scale_color_manual(values=colors, drop=False, breaks=['Non-Omicron', 'Omicron Mutation', 'Omicron Position'])+
             xlab("Position")+
             ylab('')+
             theme_classic()+
             theme(legend_position='left', figure_size=(xsize,ysize), axis_text_x=element_text(angle=45, hjust=0.5), axis_text_y=element_text(color='white'), axis_line_y=element_line(color='white'), axis_ticks_major_y=element_line(color='white'),axis_ticks_minor_y=element_line(color='white'))) # , panel_grid=element_line(color='#111111', linetype='solid', size=0.1)

            fig.save(filename = file[:-4]+'1.png')



            df['log10_Percent'] = np.log10(100 * df.abundance)


            xsize = samp_num /2.75
            ysize = seq_num /2.75

            # while (xsize > 25 or ysize > 25 ):
                # xsize = xsize / 2
                # ysize = ysize / 2

            fig2 = (ggplot(df, aes(x='sample', y='seq_id', fill='log10_Percent', limitsize=False))+
             geom_tile(color='white')+ # "white" adds gridlines here
             xlab('Sample')+
             scale_fill_cmap(limits=[-1, 2])+ #, breaks=[0, .5, 1])+
             ylab('')+
             theme_classic()+
             theme(figure_size=(xsize,ysize), axis_text_x=element_text(angle=80, hjust=0.5), axis_text_y=element_text(color='white'))) # 

            fig2.save(filename = file[:-4]+'2.png')


