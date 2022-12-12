# Code to parse through STRONG bin results and determine taxa via remote BLAST

#performs BLAST queries for each STRONG bin
def qblast(files):
    # Import Biopython tools for running remote BLAST searches
    from Bio.Blast import NCBIWWW

    # # Import Biopython SeqIO module to handle reading sequence data
    # from Bio import SeqIO

    # BLAST generates output in XML format and we need to somehow parse it!
    from Bio.Blast import NCBIXML

    for f in files[:1]:
        bin_name=f[f.find("Bin_"):f.find("/haplo")]
        fasta_string = open(f).read()
        
        try:
            print("Blasting", bin_name, flush=True)
            result_handle = NCBIWWW.qblast("blastn", "ref_prok_rep_genomes", fasta_string, megablast = True, hitlist_size = 2, alignments = 2, descriptions = 2)
        except:
            print(bin_name, "timed out")
            
        else:
            #save output as XML file
            xml_name = bin_name + ".xml"
            xml_path = output_folder + "/" + xml_name
            
            with open(xml_path, "w") as out_handle:
                out_handle.write(result_handle.read())
                result_handle.close()
                
            open_result_handle = open(xml_path)
            blast_record = NCBIXML.parse(open_result_handle)
            
            #txt file of BLAST output
            txt_file = output_folder + "/" + bin_name + "_read_me.txt"
            f = open(txt_file, "w")
            
            for item in blast_record:
                f.write(str(item.query) + '\n')
                for alignment in item.alignments:
                    for hsp in alignment.hsps:
                        if hsp.expect < 1e-50:
                            f.write(str(alignment.title))
                            f.write('\n')
            f.close() 
            
            percent_calc(txt_file, bin_name)


# coallates results into a summary file
def percent_calc(txt_file, bin_num): 
    import re
    my_dict = {}

    #parsing each line to get just the species and genus
    #my_dict: keys are taxonomy and the value is the count
    with open(txt_file) as f:
        for line in f:
            if 'COG' and 'nb0' not in line:
                x = re.split('[, /|]', line)
                target = x[5:7]
                identity = " ".join(target)
                if identity not in my_dict:
                    count = 1
                    my_dict[identity] = count
                else:
                    my_dict[identity] += 1
    sums = 0

    #finding the total number of entries
    for k in my_dict.values():
        sums += k
        
    #creating a dictionary with the same keys, but the value is a percentage, not a count. 
    percent_dict = {}
    for key,value in my_dict.items():
        percentage = round(value/sums, 3)*100
        percent_dict[key] = str(percentage)
        
    #write to a text file
    a = open(output_folder + '/bin_compositions.txt', "a")
    a.write(bin_num + '\n')

    for key, value in percent_dict.items():
        a.write(str(key) + '  ' + str(value) + '%' + '\n')
    a.write('\n')
    a.close()

# run code
import os
import glob
import argparse
argParser=argparse.ArgumentParser()
argParser.add_argument('-in','--input_dir', required=True, help='Directory of results from STRONG')
argParser.add_argument('-out','--output_dir', required=True, help='Directory to write BLAST results')
args=argParser.parse_args()

input_folder=args.input_dir
output_folder=args.output_dir


files=glob.glob(input_folder+'Bin*[0-9]/haplotypes_cogs.fna', recursive=True)
if len(files)==0:
    print('Data cannot be found in the specified directory. Please check and run again.')
else:
    if os.path.exists(output_folder)==True:
        output_folder=output_folder[:-1]+'_1/'
    os.mkdir(output_folder)
    qblast(files)
