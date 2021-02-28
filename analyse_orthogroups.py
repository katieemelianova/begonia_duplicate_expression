import os

def fasta_dict(i):
	f = open(i)
	L = {}

	for line in f:
		if line.startswith(">"):
			C_seq = ''
			C_split_line = line.split(' ')
			C_name = C_split_line[0]
			C_name = C_name.rstrip()
			C_name = C_name.lstrip('>')
		else:
			C_seq = C_seq + str(line)
			C_seq = C_seq.rstrip()
			C_seq = C_seq.upper()

		L[C_name] = C_seq
	return(L)

def write_orthogroup_fastas(orthogroup, orthogroup_name):
    for i in orthogroup:
        o_group = orthogroups_dict[i.rstrip("\n")]
        con1 = o_group[0].replace(",", "").lstrip("CON_")
        con2 = o_group[1].replace(",", "").lstrip("CON_")
        ple1 = o_group[2].replace(",", "").lstrip("PLE_")
        ple2 = o_group[3].replace(",", "").lstrip("PLE_")

        directory = os.path.join(orthogroup_name, i.rstrip('\n'))
        if not os.path.exists(directory):
            os.makedirs(directory)

        con1_of = os.path.join(directory, con1 + ".fasta")
        con2_of = os.path.join(directory, con2 + ".fasta")
        ple1_of = os.path.join(directory, ple1 + ".fasta")
        ple2_of = os.path.join(directory, ple2 + ".fasta")
        with open(con1_of, 'w') as of:
            of.write(">CON_" + con1 + '\n' + con_trinity[con1])
        with open(con2_of, 'w') as of:
            of.write(">CON_" + con2 + '\n' + con_trinity[con2])
        with open(ple1_of, 'w') as of:
            of.write(">PLE_" + ple1 + '\n' + ple_trinity[ple1])
        with open(ple2_of, 'w') as of:
            of.write(">PLE_" + ple2 + '\n' + ple_trinity[ple2])

# read in the fasta files used to make orthogroups
con_trinity = fasta_dict("con_Trinity.longest.fasta")
ple_trinity = fasta_dict("ple_Trinity.longest.fasta")

# make a dict of all orthogroups
with open("Orthogroups/Orthogroups.tsv") as openfile:
    all_orthogroups = openfile.readlines()
orthogroups_dict = {}
for l in all_orthogroups:
    o = l.split()[0]
    s = l.split()[1:]
    orthogroups_dict[o] = s

# read in each category of orthogroups
with open("all_similar.orthogroups") as openfile:
    all_similar = openfile.readlines()
with open("one_con_one_div.orthogroups") as openfile:
    one_con_one_div = openfile.readlines()
with open("similar_orthologs.orthogroups") as openfile:
    similar_orthologs = openfile.readlines()

write_orthogroup_fastas(all_similar, "all_similar")
write_orthogroup_fastas(one_con_one_div, "one_con_one_div")
write_orthogroup_fastas(similar_orthologs, "similar_orthologs")






#
