#Created Nov 2022 Dan Maskell VI6 APHA - daniel.maskell@APHA.gov.uk
#text parsing script tailor built to convert Influenza sequence data used at APHA into the required format for upload to GISAID. No usage created due to the incredibly specific and intended personal use of this script.
#generates a .csv and a .fasta file to be copied directly into the GISAID upload document
#updated 09/01/23 for readability


#imports
from Bio import Seq, SeqIO
import argparse
import os
import sys
from os.path import isfile, isdir, join

#generates csv to mirror GISAID  metadata upload format, can be pasted directly into GISAID upload document
def main(folder):
    folder += "/"
    files = file_spew(folder)
    writing = []
    fastas = []
    uploaded = []
    for file in files:
        data, name = object_creator(file)
        objected = Sample(*data)
        writing.append(line_generator(objected))
        fastas.append(fasta_generator(objected))
        uploaded.append(name)
    big_block_of_text = "\n".join(writing)
    big_block_of_fasta = "\n".join(fastas)
    #replaces missing base "-" with "N"
    big_block_of_fasta = big_block_of_fasta.replace("-", "N")
    big_block_of_names = "\n".join(uploaded)
    with open(folder+"GISAID_data.csv", "w+") as output_file_lines:
        output_file_lines.write("Isolate_Id,Segment_Ids,Isolate_Name,Subtype,Lineage,Passage_History,Location,province,sub_province,Location_Additional_info,Host,Host_Additional_info,Seq_Id (HA),Seq_Id (NA),Seq_Id (PB1),Seq_Id (PB2),Seq_Id (PA),Seq_Id (MP),Seq_Id (NS),Seq_Id (NP),Seq_Id (HE),Seq_Id (P3),Submitting_Sample_Id,Authors,Originating_Lab_Id,Originating_Sample_Id,Collection_Month,Collection_Year,Collection_Date,Antigen_Character,Adamantanes_Resistance_geno,Oseltamivir_Resistance_geno,Zanamivir_Resistance_geno,Peramivir_Resistance_geno,Other_Resistance_geno,Adamantanes_Resistance_pheno,Oseltamivir_Resistance_pheno,Zanamivir_Resistance_pheno,Peramivir_Resistance_pheno,Other_Resistance_pheno,Host_Age,Host_Age_Unit,Host_Gender,Health_Status,Note,PMID")
        output_file_lines.write(big_block_of_text)
    with open(folder+"fastas.txt", "w+") as output_file_fasta:
        output_file_fasta.write(big_block_of_fasta)
    with open(folder+"uploaded.txt", "w+") as output_file_list:
        output_file_list.write(big_block_of_names)
        
#this function finds renamed consensuses
def file_spew(argfolder_input):
	files = [
        os.path.join(dirpath, f)
        for (dirpath, dirnames, filenames) in os.walk(argfolder_input)
        for f in filenames
        ]
	
	consensus = [file for file in files if isfile(file) and file.endswith("rename.fasta")]
	return consensus

#extracts sample metadata from header, and removes |{gene} formatting from header
def object_creator(file_input):
    with open(file_input) as sample:
        records = SeqIO.parse(sample, "fasta")
        for record in records:
            record.description = ""
            #trims AIV feature off, preserves clone for later
            pipe_location = record.id.find("|") + 1
            last_pipe_location = record.id.rfind("|")
            uncut_name = record.id[:last_pipe_location]
            if pipe_location < (len(record.id) / 2):
                record.id = record.id[pipe_location:]
            record_id_list = record.id.split("/")
            country = record_id_list[2]
            species = "Other avian"
            if ("duck" in record.id
                or "Duck" in record.id
                ):
                species = "Duck"
                
            if ("chicken" in record.id
                or "Chicken" in record.id
                ):
                species = "Chicken"
                
            if ("turkey" in record.id
                or "Turkey" in record.id
                ):
                species = "Turkey"
            if ("goose" in record.id
                or "Goose" in record.id
                ):
                species = "Goose"
            until_gene = record.id.rfind("|")
            sliced = record.id[:until_gene]
            date = sliced[-10:]
            record.seq = record.seq.replace("\n", "")
            if "|PB2" in record.id:
                record.id = id_trimmer(record.id)
                strain_name = record.id
                record.id += "|PB2"
                PB2 = record
            if "|PB1" in record.id:
                record.id = id_trimmer(record.id)+"|PB1"
                PB1 = record
            if "|PA" in record.id:
                record.id = id_trimmer(record.id)+"|PA"
                PA = record
            if "|HA" in record.id:
                record.id = id_trimmer(record.id)+"|HA"
                HA = record
            if "|NP" in record.id:
                record.id = id_trimmer(record.id)+"|NP"
                NP = record
            if "|NA" in record.id:
                record.id = id_trimmer(record.id)+"|NA"
                NA = record
            if "|MP" in record.id:
                record.id = id_trimmer(record.id)+"|MP"
                MP = record
            if "|NS" in record.id:
                record.id = id_trimmer(record.id)+"|NS"
                NS = record
    return [PB2, PB1, PA, HA, NP, NA, MP, NS, strain_name, date, species, country], uncut_name

#removes VI6 standard format date - not used in GISAID uploads
def id_trimmer(record_id):
	date_end = record_id.find("|") - 1
	new_id = record_id[:date_end]
	return new_id

#generate Sample class to easily extract data for insert into final csv
class Sample:
    def __init__(self, PB2, PB1, PA, HA, NP, NA, MP, NS, strain, date, species, country):
        self.PB2 = PB2
        self.PB1 = PB1
        self.PA = PA
        self.HA = HA
        self.NP = NP
        self.NA = NA
        self.MP = MP
        self.NS = NS
        self.strain = strain
        self.date = date
        self.host = species
        self.country = country.capitalize()

    def __repr__(self):
        return repr(strain)

#formats data EXACTLY as GISAID upload sheet demands in csv
def line_generator(sample):
    output_line = ",,"
    output_line += f"{sample.strain},"
    output_line += "H5N1,"
    output_line += ","
    output_line += "Original,"
    output_line += f"Europe / United Kingdom / {sample.country},"
    output_line += ","
    output_line += ",,"
    output_line += f"{sample.host},"
    output_line += ","
    output_line += f"{sample.HA.id},{sample.NA.id},{sample.PB1.id},{sample.PB2.id},{sample.PA.id},{sample.MP.id},{sample.NS.id},{sample.NP.id},"
    output_line += ",,,,304,,,,"
    output_line += f"{sample.date}"
    return output_line

#generate formatted fasta string for upload
def fasta_generator(sample):
    PB2 = f"\">{sample.PB2.id}\n{sample.PB2.seq}\"\n"
    PB1 = f"\">{sample.PB1.id}\n{sample.PB1.seq}\"\n"
    PA = f"\">{sample.PA.id}\n{sample.PA.seq}\"\n"
    HA = f"\">{sample.HA.id}\n{sample.HA.seq}\"\n"
    NP = f"\">{sample.NP.id}\n{sample.NP.seq}\"\n"
    NA = f"\">{sample.NA.id}\n{sample.NA.seq}\"\n"
    MP = f"\">{sample.MP.id}\n{sample.MP.seq}\"\n"
    NS = f"\">{sample.NS.id}\n{sample.NS.seq}\""
    return str(PB2+PB1+PA+HA+NP+NA+MP+NS)

if __name__ == "__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("folder", help="folder containing renamed consensuses to send to GISAID")
        args = parser.parse_args()
        main(args.folder)
        print("Please check GISAID_data.csv for datasheet input, and fastas.csv for sequence sheet input")
