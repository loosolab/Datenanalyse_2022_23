import gc
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

def header_check():
    return True

def build_gff3_class(file):

    organism_class_objects = []
    dna_seq = ''
    printable_seq = ''
    fasta = False

    #with open(file) as gff3:
    fasta_extract = False
    fasta_counter = -2
    headerless_file = False

    gff3_gen = (row for row in open(file).readlines())


    #change

    # Check if the input is a headerless file
    if header_check():
        headerless_file = True
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

        elif line.startswith('##FASTA'):
            fasta_extract = True

        elif line.startswith('>'):
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
        
    gc.collect()

    return organism_class_objects




    