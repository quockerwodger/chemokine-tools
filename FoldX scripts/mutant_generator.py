# -*- coding: utf-8 -*-

#Input values
raw_seq = "APLLESQRSNSEEKANFCSTHNDEVYARFRLQMRVGVRHSPLYTPSNMCMLDIEDSVEDIEESTEKEYASTATGEAAGVNVSVALVGEGVSIPFSYIGLGFNPSLEDSYLYVNVSSRAPWVKQTSDLSANGGWGIKQVLEKELLAIQIGCDNQKFPEEPTTTPPSPVTTTLSSTTPDLNEENTENTPTTTGASVDRKRNPADIDFSLLVDPRCVTSVDLHVELRDACIDYKQESPLSLKGKYGDGELVKKEIKDVGKNHNMCSLNLNPGN"
test_mutants = {67:['E', 'G', 'A', 'M'],
                86:['V', 'G', 'A', 'S'],
                143:['L', 'G', 'A', 'S'],
                147:['Q', 'I', 'Y', 'F'],
                204:['D', 'R', 'P', 'K'],
                207:['L', 'G', 'A', 'C'],
                242:['Y', 'A', 'G']
                }

seq_len = len(raw_seq)
res_no = 0
seq = []
mutant_list = []

#Creating mutants. Cutoff must be updated before running
while res_no < seq_len:
    if (res_no+1) < 67:
        seq.append(raw_seq[(res_no)])
    if (res_no+1) == 67:
        for item in test_mutants[(res_no+1)]:
            mut_seq = []
            mut_seq.extend(seq)
            mut_seq.extend(item)
            mutant_list.append(mut_seq)
    if (res_no+1) > 67:
        if (res_no+1) not in test_mutants:
            for item in mutant_list:
                item.extend(raw_seq[res_no])
        else:
            temp_list = []
            for item in test_mutants[(res_no+1)]:
                for mutant in mutant_list:
                        mut_seq = []
                        mut_seq.extend(mutant)
                        mut_seq.extend(item)
                        temp_list.append(mut_seq)
            mutant_list = temp_list
    print(res_no)
    res_no += 1
final_list = []
for item in mutant_list:
    new_item = ''.join(item)
    final_list.append(new_item)

with open('mutant_file.txt', 'w') as f:
    for line in final_list:
        f.write("%s\n" % line)