#!/bin/env python3

import os
import sys


WWTP_dict = {   2 : [1, "BB", "Queens"],
                4 : [2, "HP", "Bronx"],
                13 : [3, "TI", "Queens"],
                14 : [4, "WI", "Man. & Bronx"],
                6 : [5, "NC", "Man. Qu. & BL"],
                7 : [6, "NR", "Manhattan"],
                10 : [7, "OB", "Staten Is."],
                12 : [8, "PR", "Staten Is."],
                8 : [9, "RH", "Brooklyn"],
                1 : [10, "26W", "Brooklyn"],
                3 : [11, "CI", "Brooklyn"],
                5 : [12, "JA", "Queens"],
                11 : [13, "OH", "Brooklyn"],
                9 : [14, "RK", "Queens"]
}


WWTP_dict2 = {   "BB" : 2,
                "HP" : 4,
                "TI" : 13,
                "WI" : 14,
                "NC" : 6,
                "NR" : 7,
                "OB" : 10,
                "PR" : 12,
                "RH" : 8,
                "26W" : 1,
                "CI" : 3,
                "Ci" : 3,
                "JA" : 5,
                "OH" : 11,
                "RK" : 9
}

cwdpath = os.getcwd()

deconv_dict_dict = {
'NTD'  : {},
'RBD'  : {},
'S1S2'  : {}
}

wwtp_deconv_dict_dict = {
'NTD'  : {},
'RBD'  : {},
'S1S2'  : {}
}

all_deconv = {
'NTD'  : {},
'RBD'  : {},
'S1S2'  : {}
}

wwtp_deconv = {
'NTD'  : {},
'RBD'  : {},
'S1S2'  : {}
}

all_deconv_sampnames = []


for subdir, dirs, files in os.walk(os.getcwd()):
    for file in files:
        sampname = ""
        if file.endswith('_deconv.tsv'):
            if "RBD" in file:
                amp = "RBD"
            if "NTD" in file:
                amp = "NTD"
            if "S1S2" in file:
                amp = "S1S2"
            splitfile = file.split('_')
            wwtp = splitfile[0]
            sampname = '-'.join(splitfile[0:3])
            if sampname in all_deconv_sampnames:
                count = 2
                while (sampname+'-'+str(count)) in all_deconv_sampnames:
                    count += 1
                sampname = sampname+'-'+str(count)
                # print(f"duplicate sample name changed: {sampname}")
            all_deconv_sampnames.append(sampname)
            deconv_dict_dict[amp][sampname] = {}
            try:
                wwtp_deconv_dict_dict[amp][wwtp][sampname] = {}
            except:
                wwtp_deconv_dict_dict[amp][wwtp] = {sampname : {}}

            try:
                samp=open(os.path.join(subdir, file), "r")
            except:
                print("can't open "+file)
            else:
                for line in samp:
                    splitline = line.strip("\n\r").split("\t")
                    try:
                        if not splitline[1] == 'Count':
                            if float(splitline[2]) >= .001:
                                deconv_dict_dict[amp][sampname][splitline[0]] = [splitline[1], splitline[2]]
                                try:
                                    all_deconv[amp][splitline[0]] += 1
                                except:
                                    all_deconv[amp][splitline[0]] = 1
                                wwtp_deconv_dict_dict[amp][wwtp][sampname][splitline[0]] = [splitline[1], splitline[2]]
                                try:
                                    wwtp_deconv[amp][wwtp][splitline[0]] += 1
                                except:
                                    try:
                                        wwtp_deconv[amp][wwtp][splitline[0]] = 1
                                    except:
                                        wwtp_deconv[amp][wwtp] = {splitline[0] : 1}
                    except:
                        pass
            samp.close()

        # if file.endswith('_chim_rm.tsv'):

for amp in deconv_dict_dict:
    # if len(deconv_dict_dict[amp]) > 0:

        # Col_deconv_fh = open('Collected_'+amp+'_Covar_Deconv.tsv',"w")
        # sorted_deconvs = sorted(all_deconv[amp])
        # Col_deconv_fh.write("\t")
        # for sampline in deconv_dict_dict[amp]:
            # Col_deconv_fh.write(sampline+"\t")
        # Col_deconv_fh.write("\nCovariant\t")
        # for sampline in deconv_dict_dict[amp]:
            # Col_deconv_fh.write("Abundance\t")
        # Col_deconv_fh.write("\n")
        # for deconv in sorted_deconvs:
            # if all_deconv[amp][deconv] > 0:
                # Col_deconv_fh.write(deconv+"\t")
                # for sample in deconv_dict_dict[amp]:
                    # try:
                        # Col_deconv_fh.write(deconv_dict_dict[amp][sample][deconv][1]+"\t")
                    # except:
                        # Col_deconv_fh.write("\t")
                # Col_deconv_fh.write("\n")
        # Col_deconv_fh.close()

    for wwtp in wwtp_deconv_dict_dict[amp]:
        if len(wwtp_deconv_dict_dict[amp][wwtp]) > 0:
            Col_deconv_fh = open(wwtp+'_Collected_'+amp+'_Covar_Deconv.tsv',"w")
            sorted_deconvs = sorted(wwtp_deconv[amp][wwtp])
            Col_deconv_fh.write("\t")
            sorted_samps = sorted(wwtp_deconv_dict_dict[amp][wwtp])
            for sampline in sorted_samps:
                Col_deconv_fh.write(sampline+"\t")
            Col_deconv_fh.write("\nCovariant\t")
            for sampline in sorted_samps:
                Col_deconv_fh.write("Abundance\t")
            Col_deconv_fh.write("\n")
            for deconv in sorted_deconvs:
                if wwtp_deconv[amp][wwtp][deconv] > 0:
                    Col_deconv_fh.write(deconv+"\t")
                    for sample in sorted_samps:
                        try:
                            Col_deconv_fh.write(wwtp_deconv_dict_dict[amp][wwtp][sample][deconv][1]+"\t")
                        except:
                            Col_deconv_fh.write("\t")
                    Col_deconv_fh.write("\n")
            Col_deconv_fh.close()
