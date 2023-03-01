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

def header_check(file_to_check):
    with open(file_to_check, 'r') as file_to_check:
        return (True if not file_to_check.readline().startswith('#') else False)
    

def build_gff3_class(file: list):
    organism_class_objects = []
    for path in file:
        dna_seq = ''
        printable_seq = ''
        fasta = False


        fasta_extract = False
        fasta_counter = -2
        headerless_file = header_check(path)

        gff3_gen = (row for row in open(path).readlines())

        print('Generating GTF-Objects for' + path)

        # Check if the input is a headerless file
        if headerless_file:
            # headerless_file = True
            tmp_organism = gff.Organism()
            tmp_organism.strain = path.split('/')[-1]

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
            element.set_annotated_dna_seq(fasta_extract=fasta, filename=path.split('/')[-1])

        print('Objects done')
        
    return organism_class_objects

def generate_feature_files(gtf_file, fragments, enhancer_bed, blacklisted_bed, threshold, promoter_distance, tss_distance, out):
    
    object_list = build_gff3_class(file=gtf_file)
    
    
    for i, element in enumerate(object_list):
        features = element.count_features()
        
        
        peak_gtf = element.generate_peak_gtf(fragments=fragments, threshold=threshold, gtf_file=gtf_file[i], out=out)
        peak_object_list = build_gff3_class(file=[peak_gtf])
        
        
        # Calculating feature's only for peak filtered lines
        for ele in peak_object_list:
            ele.generate_feature_gtf(feature_keys=features, out=out)
            ele.generate_promoter_gtf(promoter_distance=promoter_distance, out=out)
            ele.generate_tss_gtf(tss_distance=tss_distance, out=out)
            ele.generate_gene_body_gtf(out=out)
            ele.generate_enhancer_gtf(gtf_file=peak_gtf, enhancer_bed=enhancer_bed, out=out)
            ele.generate_blacklisted_region_gtf(gtf_file=peak_gtf, blacklisted_bed=blacklisted_bed, out=out)        
        
        element.generate_feature_gtf(feature_keys=features, out=out)        
        element.generate_promoter_gtf(promoter_distance=promoter_distance, out=out)
        element.generate_tss_gtf(tss_distance=tss_distance, out=out)
        element.generate_gene_body_gtf(out=out)        
        element.generate_enhancer_gtf(gtf_file=gtf_file[i], enhancer_bed=enhancer_bed, out=out)
        element.generate_blacklisted_region_gtf(gtf_file=gtf_file[i], blacklisted_bed=blacklisted_bed, out=out)
        


    