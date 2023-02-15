import os
import sys
import pandas as pd
import gff_analyser.gffClasses as gff

def strain_exists(name: str, data: list):
    return next((True for organism in data if organism.strain == name), False)


def find_strain(name: str, data: list):
    return next((organism for organism in data if organism.strain == name), None)


        

def add_sequence(data: list, dna_seq: str, fasta_counter: int, printable_seq: str):

    if dna_seq != '':
        tmp_organism: Organism = find_strain(name=data[fasta_counter].strain,
                                                        data=data)
        tmp_organism.fasta += dna_seq
        tmp_organism.printable_fasta += printable_seq
        dna_seq = ''
    return dna_seq

def header_check(file):
    with open(file, 'r') as file:
        return (True if not file.readline().startswith('#') else False)
    

def build_gff3_class(file):

    organism_class_objects = []
    dna_seq = ''
    printable_seq = ''
    fasta = False

    #with open(file) as gff3:
    fasta_extract = False
    fasta_counter = -2
    headerless_file = header_check(file)

    gff3_gen = (row for row in open(file).readlines())


    print('Generating GTF-Objects')

    # Check if the input is a headerless file
    if headerless_file:
        #headerless_file = True
        tmp_organism = gff.Organism()
        tmp_organism.strain = file.split('/')[-1]


        organism_class_objects.append(tmp_organism)

    for line in gff3_gen:


        if line.startswith("##sequence-region"):
            strain = line.split(" ")
            tmp_organism = gff.Organism()
            tmp_organism.strain = strain[1]
            organism_class_objects.append(tmp_organism)

        # Possible Bug, if headerless file contains a sequence
        elif strain_exists(line.split("\t")[0], organism_class_objects) or headerless_file:
            gffrow = line.split("\t")
            if not headerless_file:
                tmp_organism = find_strain(name=gffrow[0], data=organism_class_objects)

            tmp_organism.gff_data.append(gff.GffData(gffrow=gffrow))

        elif line.startswith('##FASTA') and not headerless_file:
            fasta_extract = True

        elif line.startswith('>') and not headerless_file:
            fasta_counter += 1
            dna_seq = add_sequence(data=organism_class_objects, dna_seq=dna_seq, fasta_counter=fasta_counter,
                                   printable_seq=printable_seq)

        elif fasta_extract and not line.startswith('>'):
            dna_seq += line.strip("\n")
            if fasta:
                printable_seq += line

    

    add_sequence(data=organism_class_objects, dna_seq=dna_seq, fasta_counter=fasta_counter + 1,
                 printable_seq=printable_seq)


    for element in organism_class_objects:

        element.set_annotated_dna_seq(fasta_extract=fasta, filename=file.split('/')[-1])
        
    print('Objects done')

    return organism_class_objects

def generate_feature_files(gtf_file, fragments, enhancer_bed, blacklisted_bed, threshold, promotor_distance, tss_distance, out='out'):
    
    object_list = build_gff3_class(file=gtf_file)
    
    
    for element in object_list:
        features = element.count_features()

        element.generate_feature_gtf(gffdata_list=object_list, feature_keys=features, out=out)
        element.generate_promotor_gtf(gffdata_list=object_list, promotor_distance=promotor_distance, out=out)
        element.generate_tss_gtf(gffdata_list=object_list, tss_distance=tss_distance, out=out)
        element.generate_gene_body_gtf(gffdata_list=object_list, out=out)
        
        strain = element.strain
        
        # Getting bedtools Path
        bedtools = os.path.join('/'.join(sys.executable.split('/')[:-1]),'bedtools')

        
        
        # Load the bed file into a pandas dataframe
        bed_df = pd.read_table(fragments, header=None, names=["chrom", "start", "end", "feature", "count", "strand"])

        # Filter the dataframe based on a count threshold
        
        filtered_df = bed_df[bed_df["count"] > threshold]

        filename = f"{fragments}.peak{threshold}.bed"
        # Print the filtered dataframe
        filtered_df.to_csv(filename ,sep='\t', index=None) 
        
        peakfile = f"{fragments}.peak{threshold}.bed"

        # Getting the peak-GTF
        intersect_cmd = f'{bedtools} intersect -a {gtf_file} -b {peakfile} > {out}{strain}.peak{threshold}.gtf'
        os.system(intersect_cmd)
        print("Peaks done")

        # Getting the enhancer-GTF 
        intersect_cmd = f'{bedtools} intersect -a {gtf_file} -b {enhancer_bed} > {out}{strain}.enhancer.gtf'
        os.system(intersect_cmd)
        print("Enhancers done")

        # Getting the exclusion_list-GTF
        intersect_cmd = f'{bedtools} intersect -a {gtf_file}  -b {blacklisted_bed} > {out}{strain}.blacklisted.gtf'
        os.system(intersect_cmd)
        print("Blacklist done")


    