import os
import sys

def get_complementary_string(sequence: str):
    complementary_string = ''
    for base in sequence[::-1]:
        if base == 'A':
            complementary_string += 'T'
        elif base == 'T':
            complementary_string += 'A'
        elif base == 'G':
            complementary_string += 'C'
        elif base == 'C':
            complementary_string += 'G'
        elif base == 'N':
            complementary_string += 'N'
 

    return complementary_string


class GffData:
    def __init__(self, gffrow):
        self._seq_id = gffrow[0]
        self._source = gffrow[1]
        self._feature_type = gffrow[2]
        self._feature_start = int(gffrow[3])
        self._feature_end = int(gffrow[4])
        self._score = gffrow[5]
        self._strand = gffrow[6]
        self._phase = gffrow[7]
        self._attributes = gffrow[8].split(";")
        self._dnaseq = ''

    def get_whole_line(self, start, end):

        tmp_attribute = ''

        # Warum sind in den gff3 files die Attribute ohne ; am Ende aber in den gtf's mit??
        # Bug hier noch in der Darstellung für gff3 files

        # Adding list object to tmp_string to get a printable attribute line
        for entry in self.attributes:
            tmp_attribute += entry

            # Adding at last entry ; for formating purposes
            if entry != entry[-1]:
                tmp_attribute += ';'



        whole_line = "{seq_id}\t{source}\t{feature_type}\t{feature_start}\t{feature_end}\t{score}\t{strand}\t{phase}\t{tmp_attribute}".format(
            seq_id=self.seq_id, source=self.source, feature_type=self.feature_type,
            feature_start=start, feature_end=end, score=self.score,
            strand=self.strand, phase=self.phase, tmp_attribute=tmp_attribute)

        return whole_line

    @property
    def seq_id(self):
        return self._seq_id

    @seq_id.setter
    def seq_id(self, value: str):
        self._seq_id = value

    @property
    def source(self):
        return self._source

    @source.setter
    def source(self, value: str):
        self._source = value

    @property
    def feature_type(self):
        return self._feature_type

    @feature_type.setter
    def feature_type(self, value: str):
        self._feature_type = value

    @property
    def feature_start(self):
        return self._feature_start

    @feature_start.setter
    def feature_start(self, value: int):
        self._feature_start = value

    @property
    def feature_end(self):
        return self._feature_end

    @feature_end.setter
    def feature_end(self, value: int):
        self._feature_end = value

    @property
    def score(self):
        return self._score

    @score.setter
    def score(self, value: str):
        self._score = value

    @property
    def strand(self):
        return self._strand

    @strand.setter
    def strand(self, value: str):
        self._strand = value

    @property
    def phase(self):
        return self._phase

    @phase.setter
    def phase(self, value: str):
        self._phase = value

    @property
    def attributes(self):
        return self._attributes

    @attributes.setter
    def attributes(self, value: list):
        self._attributes = value

    @property
    def dnaseq(self):
        return self._dnaseq

    @dnaseq.setter
    def dnaseq(self, value: str):
        self._dnaseq = value





class Organism:
    def __init__(self):

        self._strain = ""
        self._fasta = ""
        self._printable_fasta = ""
        self._gff_data = []

    def generate_feature_gtf(self, gffdata_list, feature_keys, out='out/'):
    #feature_list from .count_features

        try:
            os.mkdir("{out}".format(out=out))
        except:
            print('out already exist, skipping creation')

        feature_count = -1
        feature_list = list(feature_keys.keys())

        for feature in feature_list:
            feature_count += 1

            print('Generating ' + feature + '-File')
            for element in gffdata_list:

                if not element.strain.endswith('.gtf'):
                    element.strain += '.gtf'


                filename = "{out}{strain}.{feature}.gtf".format(out=out,
                                                                 strain=element.strain.strip('.gtf'),
                                                                 feature=feature_list[feature_count])

                with open(filename, 'a') as gtf_file:
                    for row in element.gff_data:
                        if row.feature_type == feature_list[feature_count]:
                            gtf_file.write(row.get_whole_line(start=row.feature_start, end=row.feature_end))
                            
    def generate_gene_body_gtf(self, gffdata_list, out='out/'):
                     
        try:
            os.mkdir("{out}".format(out=out))
        except:
            print('out already exist, skipping creation')
        
        for element in gffdata_list:
        
            # Getting the gene_body-GTF
            bedtools = os.path.join('/'.join(sys.executable.split('/')[:-1]),'bedtools')
            intersect_cmd = "{bed} subtract -a {out}{strain}.exon.gtf -b {out}{strain}.five_prime_utr.gtf {out}{strain}.three_prime_utr.gtf > {out}{strain}.exon.gene_bodies.gtf".format(bed=bedtools,
                                                                            out=out,
                                                                            strain=element.strain.strip('.gtf'))
            os.system(intersect_cmd)
            print("Gene_bodys done")
                            
    def generate_promotor_gtf(self, gffdata_list, promotor_distance, out='out/'):
        try:
            os.mkdir("{out}".format(out=out))
        except:
            print('out already exist, skipping creation')
            
        

        for element in gffdata_list:

            feature = 'gene'
            filename_promotor = "{out}{strain}.{feature}.promotor{distance}.gtf".format(out=out,
                                                                                         strain=element.strain.strip('.gtf'),
                                                            feature=feature, distance=promotor_distance)
            
            with open(filename_promotor, 'w') as promotor_file:
                print('Generating Promotor-File')
                row_counter = 0
                for row in element.gff_data:
                    row_counter += 1
                    if row.feature_type == feature and row.strand == '+' and row.feature_start - promotor_distance > 0:
                        # for positive strands
                        promotor_file.write(row.get_whole_line(start=row.feature_start - promotor_distance, end=row.feature_start))
                        #promotor_file.write('\n')
                    elif row.feature_type == feature and row.strand == '-':
                        promotor_file.write(row.get_whole_line(start=row.feature_end, end=row.feature_end + promotor_distance))
                        #promotor_file.write('\n')
                        
        

    def generate_tss_gtf(self, gffdata_list, tss_distance=100, out='out/'):
        
        try:
            os.mkdir("{out}".format(out=out))
        except:
            print('out already exist, skipping creation')


            
        
        
        feature = 'gene'
        for element in gffdata_list:
        # Getting the TSS file
            filename_tss = "{out}{strain}.{feature}.TSS{tss_distance}.gtf".format(out=out, strain=element.strain.strip('.gtf'),
                                                                  feature=feature, tss_distance=tss_distance)

         
            with open(filename_tss, 'w') as tss_file:
                print('Generating TSS-File')
                for row in element.gff_data:
                    if row.feature_type == feature and row.strand == '+':
                        # for positive strands
                        tss_file.write(row.get_whole_line(start=row.feature_start, end=row.feature_start + tss_distance))
                        #tss_file.write('\n')
                    elif row.feature_type == feature and row.strand == '-':
                        tss_file.write(row.get_whole_line(start=row.feature_end - tss_distance, end=row.feature_end ))
                        #tss_file.write('\n')


       

    def count_features(self):

        feature_dict = {}
        for entry in self.gff_data:
            if entry.feature_type not in feature_dict.keys():
                feature_dict.update({entry.feature_type: 0})
            feature_dict[entry.feature_type] += 1

        return feature_dict


    def set_annotated_dna_seq(self, fasta_extract: bool, filename):

        if fasta_extract:
            with open(filename.replace('.gff3', '.fa'), 'w') as fasta_file:
                fasta_file.write("##sequence-region {strain}\n".format(strain=self.strain))

        for element in self.gff_data:

            if element.feature_type == 'region':
                element.dnaseq = self.fasta
                if fasta_extract:
                    with open(filename.replace('.gff3', '.fa'), 'a') as fasta_file:
                        fasta_file.write(">{idname}\n{sequence}\n".format(idname=element.attributes[0].split('=')[-1],
                                                                          sequence=self._printable_fasta))

            elif element.strand == '+':
                element.dnaseq = self.fasta[int(element.feature_start) - 1: int(element.feature_end)]

            elif element.strand == '-':
                unreverted_sequence = self.fasta[int(element.feature_start) - 1: int(element.feature_end)]
                element.dnaseq = get_complementary_string(sequence=unreverted_sequence)

    def print_multiple_fasta(self, data_list: list, filename: str):
        with open(filename.replace('.gff3', '.mpfa'), 'w') as file:
            for element in data_list:

                file.write("##sequence-region {strain}\n".format(strain=element.strain))

                for entry in element.gff_data[1:]:
                    file.write(">{idname}\n{sequence}\n".format(idname=entry.attributes[0].split('=')[-1],
                                                                sequence=entry.dnaseq))



    @property
    def strain(self):
        return self._strain

    @strain.setter
    def strain(self, value: str):
        self._strain = value

    @property
    def fasta(self):
        return self._fasta

    @fasta.setter
    def fasta(self, value: str):
        self._fasta = value

    @property
    def printable_fasta(self):
        return self._printable_fasta

    @printable_fasta.setter
    def printable_fasta(self, value: str):
        self._printable_fasta = value

    @property
    def gff_data(self):
        return self._gff_data

    @gff_data.setter
    def gff_data(self, value: GffData):
        self._gff_data = value