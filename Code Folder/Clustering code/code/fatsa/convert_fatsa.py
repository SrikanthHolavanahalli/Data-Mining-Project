mupro_file = "mupro.data"
protein_fatsa_file = "protein_fatsa.data"
protein_chain_fatsa_file = "protein_chain_fatsa.data"
protein_fatsa_format = {}
protein_chain_fatsa_format = {}

def read_file():
	no_lines = 0
	protein = ""
	protein_chain = ""
	data = open(mupro_file)
	for line in data:
		no_lines = no_lines+1
		if(no_lines%10==1):
			protein = (line.strip())[:-1]
			protein_chain = line
		elif(no_lines%10==3):
			if(not protein_fatsa_format.has_key(protein)):
				protein_fatsa_format[protein] = line
			if(not protein_chain_fatsa_format.has_key(protein_chain)):
				protein_chain_fatsa_format[protein_chain] = line

def write_file():
	output = open(protein_fatsa_file,"w")
	for key in protein_fatsa_format:
		output.write(">"+key+"\n")
		output.write(protein_fatsa_format[key])
	output.close()

	output = open(protein_chain_fatsa_file,"w")
	for key in protein_chain_fatsa_format:
		output.write(">"+key)
		output.write(protein_chain_fatsa_format[key])
	output.close()

read_file()
write_file()
