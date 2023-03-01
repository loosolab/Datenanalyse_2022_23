import os
import sys
import pandas as pd
import numpy as np
import episcanpy.api as epi

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


def check_for_output_dir(out):
    try:
        os.mkdir("{out}".format(out=out))
    except:
        print('out already exist, skipping creation')


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

    def generate_feature_gtf(self, feature_keys: list, out: str):
 

        check_for_output_dir(out)

        feature_count = -1
        feature_list = list(feature_keys.keys())

        for feature in feature_list:
            feature_count += 1

            print('Generating ' + feature + '-File')

            if not self.strain.endswith('.gtf'):
                self.strain += '.gtf'

            if self.strain.split('.')[-2].startswith('peak'):
                filename = f"{out}{self.strain.strip('.gtf')}_{feature_list[feature_count]}.gtf"
            else:
                filename = f"{out}{self.strain.strip('.gtf')}.{feature_list[feature_count]}.gtf"

            with open(filename, 'a') as gtf_file:
                for row in self.gff_data:
                    if row.feature_type == feature_list[feature_count]:
                        gtf_file.write(row.get_whole_line(start=row.feature_start, end=row.feature_end))

    def generate_gene_body_gtf(self, out: str):

        check_for_output_dir(out)
        if self.strain.split('.')[-2].startswith('peak'):
            filename_body = f"{out}{self.strain.strip('.gtf')}_gene_bodies.gtf"
            infileA = f"{out}{self.strain.strip('.gtf')}_gene.gtf"
            infileB1 = f"{out}{self.strain.strip('.gtf')}_five_prime_utr.gtf"
            infileB2 = f"{out}{self.strain.strip('.gtf')}_three_prime_utr.gtf"
        else:
            filename_body = f"{out}{self.strain.strip('.gtf')}.gene_bodies.gtf"
            infileA = f"{out}{self.strain.strip('.gtf')}.gene.gtf"
            infileB1 = f"{out}{self.strain.strip('.gtf')}.five_prime_utr.gtf"
            infileB2 = f"{out}{self.strain.strip('.gtf')}.three_prime_utr.gtf"
        

        bedtools = os.path.join('/'.join(sys.executable.split('/')[:-1]), 'bedtools')
        intersect_cmd = f"{bedtools} subtract -a {infileA} -b {infileB1} {infileB2} > {filename_body}"
        os.system(intersect_cmd)
        print("Gene_bodys done")

    def generate_promoter_gtf(self, promoter_distance: int, out: str):
        check_for_output_dir(out)

        feature = 'gene'
        
        if self.strain.split('.')[-2].startswith('peak'):
            filename_promotor = f"{out}{self.strain.strip('.gtf')}_promotor{promoter_distance}.gtf"
        else:
            filename_promotor = f"{out}{self.strain.strip('.gtf')}.promotor{promoter_distance}.gtf"
                                                                                     
                                                                                                                                                                                                                                                
        with open(filename_promotor, 'w') as promotor_file:
            print('Generating Promotor-File')
            row_counter = 0
            for row in self.gff_data:
                row_counter += 1
                if row.feature_type == feature and row.strand == '+' and row.feature_start - promoter_distance > 0:
                    # for positive strands
                    promotor_file.write(
                        row.get_whole_line(start=row.feature_start - promoter_distance, end=row.feature_start))
                    # promotor_file.write('\n')
                elif row.feature_type == feature and row.strand == '-':
                    promotor_file.write(
                        row.get_whole_line(start=row.feature_end, end=row.feature_end + promoter_distance))
                    # promotor_file.write('\n')

    def generate_tss_gtf(self, tss_distance: int, out: str):

        check_for_output_dir(out)

        feature = 'gene'
        
        if self.strain.split('.')[-2].startswith('peak'):
            filename_tss = f"{out}{self.strain.strip('.gtf')}_TSS{tss_distance}.gtf"
        else:
            filename_tss = f"{out}{self.strain.strip('.gtf')}.TSS{tss_distance}.gtf"
                                                                               
                                                                               

        with open(filename_tss, 'w') as tss_file:
            print('Generating TSS-File')
            for row in self.gff_data:
                if row.feature_type == feature and row.strand == '+':
                    # for positive strands
                    tss_file.write(
                        row.get_whole_line(start=row.feature_start, end=row.feature_start + tss_distance))
                    # tss_file.write('\n')
                elif row.feature_type == feature and row.strand == '-':
                    tss_file.write(row.get_whole_line(start=row.feature_end - tss_distance, end=row.feature_end))
                    # tss_file.write('\n')

    def generate_noise_filtered_gtf(self, fragments: str, threshold: int, gtf_file: str, out: str):

        check_for_output_dir(out)

        bedtools = os.path.join('/'.join(sys.executable.split('/')[:-1]), 'bedtools')

        # Load the bed file into a pandas dataframe
        bed_df = pd.read_table(fragments, header=None, names=["chrom", "start", "end", "feature", "count", "strand"])

        # Filter the dataframe based on a count threshold
        filtered_df = bed_df[bed_df["count"] > threshold]

        # Print the filtered dataframe
        peakfile = f"{fragments}_{self.strain.strip('.gtf')}.peak{threshold}.bed"
        filtered_df.to_csv(peakfile, sep='\t', index=None)

        

        # Getting the peak-GTF
        intersect_cmd = f"{bedtools} intersect -a {gtf_file} -b {peakfile} > {out}{self.strain.strip('.gtf')}.peak{threshold}.gtf"
        os.system(intersect_cmd)
        print("Peaks done")
        
        
        
        return f"{out}{self.strain.strip('.gtf')}.peak{threshold}.gtf"

    
    
    def generate_peak_gtf(self, adata, gtf_file, out):
        
        bedtools = os.path.join('/'.join(sys.executable.split('/')[:-1]), 'bedtools')

        filename_tmp = adata.obs['tissue'].values[0]
        output_file_path = f"{out}{filename_tmp}_peaks.bed"

        

        with open(output_file_path, 'a') as output_file:

            for entry in adata.var.index:
                entry_list = entry.split('_')
                entry_formated = f"{entry_list[0]}\t{entry_list[1]}\t{entry_list[2]}\t.\t.\n"
                
                output_file.write(entry_formated)

        intersect_cmd = f"{bedtools} intersect -a {gtf_file} -b {output_file_path} > {out}{self.strain.strip('.gtf')}.peak.gtf"
        os.system(intersect_cmd)
        return f"{out}{self.strain.strip('.gtf')}.peak.gtf"
    

    def generate_enhancer_gtf(self, gtf_file: str, enhancer_bed: str, out: str):
        
        if self.strain.split('.')[-2].startswith('peak'):
            filename_enhancer = f"{out}{self.strain.strip('.gtf')}_enhancer.gtf"
        else:
            filename_enhancer = f"{out}{self.strain.strip('.gtf')}.enhancer.gtf"
        
        bedtools = os.path.join('/'.join(sys.executable.split('/')[:-1]), 'bedtools')
        # Getting the enhancer-GTF
        intersect_cmd = f"{bedtools} intersect -a {gtf_file} -b {enhancer_bed} > {filename_enhancer}"
        os.system(intersect_cmd)
        print("Enhancers done")

    def generate_blacklisted_region_gtf(self, gtf_file: str, blacklisted_bed: str, out: str):
        
        if self.strain.split('.')[-2].startswith('peak'):
            filename_blacklisted = f"{out}{self.strain.strip('.gtf')}_blacklisted.gtf"
        else:
            filename_blacklisted = f"{out}{self.strain.strip('.gtf')}.blacklisted.gtf"
            
        bedtools = os.path.join('/'.join(sys.executable.split('/')[:-1]), 'bedtools')
        # Getting the exclusion_list-GTF
        intersect_cmd = f"{bedtools} intersect -a {gtf_file}  -b {blacklisted_bed} > {filename_blacklisted}"
        os.system(intersect_cmd)
        print("Blacklist done")


    def count_features(self):

        feature_dict = {}
        for entry in self.gff_data:
            if entry.feature_type not in feature_dict.keys():
                feature_dict.update({entry.feature_type: 0})
            feature_dict[entry.feature_type] += 1

        return feature_dict

    def set_annotated_dna_seq(self, fasta_extract: bool, filename: str):

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