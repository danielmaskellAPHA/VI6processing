#update pending to improve readability
#updated 19/12/2022 to remove "A_/_" from headers - DM
#created 06/10/2022 splitting function written by Dan Maskell, GISAID metadata cleaning written by Alex Byrne

from Bio import Seq, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import argparse
import os
import sys
from os import listdir
from os.path import isfile, isdir, join
#arguments
parser = argparse.ArgumentParser()
parser.add_argument("folder", help="directory containing consensuses")
args = parser.parse_args()
argfolder = args.folder+"/"
thepath, the_submission = os.path.split(args.folder)
print("Renaming and splitting "+the_submission+"...")

if not isdir(argfolder):
	print("ERROR: Please provide a directory! \n Terminating...")
	sys.exit()
files = []
pre_files = []
samples = []
metadata_check = []
gisaid = []
def file_spew(dir):
	for file in dir:
		if isfile(file) and file.endswith(".fasta"):
			files.append(file)
		elif isfile(file) and file.endswith(".xls"):
			metadata_check.append(file)
		elif isfile(file) and file.endswith("GISAID.txt"):
			print("Found GISAID")
			gisaid.append(file)

found_files = os.listdir(argfolder)

for file in found_files:
	file2 = argfolder+file
	pre_files.append(file2)
file_spew(pre_files)
for submission in files:
	path, sample = os.path.split(submission)
	samples.append(sample)
if len(samples) == 0:
	print("No samples located in " + argfolder)
	print("Terminating.")
	sys.exit()
if len(gisaid) == 0:
	print("GISAID.txt not found. You should not be using this script.")
	sys.exit()
print(str(len(samples)) + " sample(s) located in " + argfolder)
for sample in samples:
	print(sample)
segments = ["|PB2", "|PB1", "|PA", "|HA", "|NP", "|NA", "|MP", "|NS"]
gene_folder = argfolder+"Genes/"
if os.path.exists(gene_folder) == False:
	os.makedirs(gene_folder)
for seg in segments:
	open(gene_folder+f"{seg[1:]}.fasta", "w+")
def renamer(file):
	#iterates through records, moves up the list by one when finished
	with open(file) as original:
		records = SeqIO.parse(original, "fasta")
		for record in records:
			for seg in segments:
				current_gene = seg
				if seg in record.id:
					gene_path = open(gene_folder+f"{current_gene[1:]}.fasta", "a+")
					record.id = record.id.replace(f"{current_gene}", "")
					record.id = record.id.replace("A_/_", "")
					to_write = record.format("fasta")
					gene_path.write(to_write)
for file in files:
	renamer(file)
for seg in segments:
	print("Created "+f"{seg[1:]}.fasta")
print("Please double check output files to ensure everything is as expected.")
#GISAID metadata section
if len(metadata_check) == 0:
	print("No metadata file found. GISAID Metadata cleaning has not been performed.")
	sys.exit()
metadata = metadata_check[0]
df=pd.read_excel(metadata)
# Delete unecessary columna
df.drop(['HE Segment_Id', 'P3 Segment_Id', 'Lineage', 'Passage_History', 'Submitting_Sample_Id', 'Antigen_Character', 'Originating_Sample_Id', 'Note', 'Update_Date', 'Submission_Date', 'Animal_Vaccin_Product', 'Adamantanes_Resistance_geno', 'Oseltamivir_Resistance_geno', 'Zanamivir_Resistance_geno', 'Peramivir_Resistance_geno', 'Other_Resistance_geno', 'Adamantanes_Resistance_pheno', 'Oseltamivir_Resistance_pheno', 'Zanamivir_Resistance_pheno', 'Peramivir_Resistance_pheno', 'Other_Resistance_pheno', 'Host_Age', 'Host_Age_Unit', 'Host_Gender', 'Patient_Status', 'Zip_Code', 'Outbreak', 'Pathogen_Test_Info', 'Is_Vaccinated', 'Human_Specimen_Source', 'Animal_Specimen_Source', 'Animal_Health_Status', 'Domestic_Status', 'PMID', 'PB2 INSDC_Upload', 'PB1 INSDC_Upload', 'PA INSDC_Upload', 'HA INSDC_Upload', 'NP INSDC_Upload', 'NA INSDC_Upload', 'MP INSDC_Upload', 'NS INSDC_Upload', 'HE INSDC_Upload', 'P3 INSDC_Upload'], axis=1, inplace=True)

# Reorder columns
df_columns = [col for col in df.columns if col != 'Isolate_Name']
df_columns.insert(0, 'Isolate_Name')
df = df[df_columns]
df_columns = [col for col in df.columns if col != 'Subtype']
df_columns.insert(1, 'Subtype')
df = df[df_columns]
df_columns = [col for col in df.columns if col != 'Collection_Date']
df_columns.insert(2, 'Collection_Date')
df = df[df_columns]
df_columns = [col for col in df.columns if col != 'Location']
df_columns.insert(3, 'Location')
df = df[df_columns]
df_columns = [col for col in df.columns if col != 'Host']
df_columns.insert(4, 'Host')
df = df[df_columns]

# Remove duplicated entries based on 'Isolate_Name'
df = df.drop_duplicates(subset=['Isolate_Name'], keep=False)

# Renamed Segment_Id columns to remove spaces
df.rename(columns={'PB2 Segment_Id': 'PB2_Segment_ID'}, inplace=True)
df.rename(columns={'PB1 Segment_Id': 'PB1_Segment_ID'}, inplace=True)
df.rename(columns={'PA Segment_Id': 'PA_Segment_ID'}, inplace=True)
df.rename(columns={'HA Segment_Id': 'HA_Segment_ID'}, inplace=True)
df.rename(columns={'NP Segment_Id': 'NP_Segment_ID'}, inplace=True)
df.rename(columns={'NA Segment_Id': 'NA_Segment_ID'}, inplace=True)
df.rename(columns={'MP Segment_Id': 'MP_Segment_ID'}, inplace=True)
df.rename(columns={'NS Segment_Id': 'NS_Segment_ID'}, inplace=True)

# Rearrange the subtype to remove "A / "
df['Subtype'] = df['Subtype'].str.split(' / ').str[1]

# Adjust the Segment_ID columns to only retain the actual accession ID
df['PB2_Segment_ID'] = df['PB2_Segment_ID'].str.split('|').str[0]
df['PB1_Segment_ID'] = df['PB1_Segment_ID'].str.split('|').str[0]
df['PA_Segment_ID'] = df['PA_Segment_ID'].str.split('|').str[0]
df['HA_Segment_ID'] = df['HA_Segment_ID'].str.split('|').str[0]
df['NP_Segment_ID'] = df['NP_Segment_ID'].str.split('|').str[0]
df['NA_Segment_ID'] = df['NA_Segment_ID'].str.split('|').str[0]
df['MP_Segment_ID'] = df['MP_Segment_ID'].str.split('|').str[0]
df['NS_Segment_ID'] = df['NS_Segment_ID'].str.split('|').str[0]

# Split the Location column into 'Region', 'Country', 'Province_City' and then reorder and drop 'Location'
df['Region'] = df['Location'].str.split(' / ').str[0]
df['Country'] = df['Location'].str.split(' / ').str[1]
df['Province_City'] = df['Location'].str.split(' / ').str[2]
df_columns = [col for col in df.columns if col != 'Region']
df_columns.insert(5, 'Region')
df = df[df_columns]
df_columns = [col for col in df.columns if col != 'Country']
df_columns.insert(6, 'Country')
df = df[df_columns]
df_columns = [col for col in df.columns if col != 'Province_City']
df_columns.insert(7, 'Province_City')
df = df[df_columns]
df.drop(['Location'], axis=1, inplace=True)

# Rename some of the countries to align with augur
df['Country'] = df['Country'].str.replace("Korea, Republic of", "South Korea", regex=False)
df['Country'] = df['Country'].str.replace("Russian Federation", "Russia", regex=False)
df['Country'] = df['Country'].str.replace("Lao, People's Democratic Republic", "Laos", regex=False)
df['Country'] = df['Country'].str.replace("Hong Kong (SAR)", "Hong Kong", regex=False)
df['Country'] = df['Country'].str.replace("Iran, Islamic Republic of", "Iran", regex=False)
df['Country'] = df['Country'].str.replace("Congo, the Democatic Republic of", "Congo", regex=False)
df['Country'] = df['Country'].str.replace("Macedonia, the former Yogoslav Republic of", "Macedonia", regex=False)

# Create a column called location, for UK report cases this will be the AIV number, for GISAID data, I will copy the country over as the Location
df['Location'] = df['Country'].copy()
df_columns = [col for col in df.columns if col != 'Location']
df_columns.insert(7, 'Location')
df = df[df_columns]

# Pull the year from the 'Collection_Date' an add to a separate column then sort by 'Collection_Date' in descending order
df['Year'] = df['Collection_Date'].str.split('-').str[0]
df_columns = [col for col in df.columns if col != 'Year']
df_columns.insert(3, 'Year')
df = df[df_columns]
df = df.sort_values(['Collection_Date'], ascending=[False])

# Replace any spaces and any subtypes in bracket '(subtype)' in the 'Isolate_Name' otherwise this will be a problem later on
df['Isolate_Name'] = df['Isolate_Name'].str.replace(' ', '_', regex=False)
df['Isolate_Name'] = df['Isolate_Name'].str.replace(r'\(.*\)','', regex=True)

# Remove rows with incomplete subtypes or H0N0 and missing collection date
df1 = df[(df['Subtype'] != "H1") & (df['Subtype'] != "H2") & (df['Subtype'] != "H3") & (df['Subtype'] != "H4") & (df['Subtype'] != "H5") & (df['Subtype'] != "H6") & (df['Subtype'] != "H7") & (df['Subtype'] != "H8") & (df['Subtype'] != "H9") & (df['Subtype'] != "H10") & (df['Subtype'] != "H11") & (df['Subtype'] != "H12") & (df['Subtype'] != "H13") & (df['Subtype'] != "H14") & (df['Subtype'] != "H15") & (df['Subtype'] != "H16") & (df['Subtype'] != "H17") & (df['Subtype'] != "H18")]
df1 = df1[(df1['Subtype'] != "N1") & (df1['Subtype'] != "N2") & (df1['Subtype'] != "N3") & (df1['Subtype'] != "N4") & (df1['Subtype'] != "N5") & (df1['Subtype'] != "N6") & (df1['Subtype'] != "N7") & (df1['Subtype'] != "N8") & (df1['Subtype'] != "N9") & (df1['Subtype'] != "N10") & (df1['Subtype'] != "N11")]
df2 = df1.dropna(subset=['Subtype'])
df2 = df2[(df2['Subtype'] != 'H0N0')]
df2 = df2.dropna(subset=['Collection_Date'])
Metadataoutput = os.path.splitext(metadata)[0]
df2.to_excel(Metadataoutput+".cleaned.xlsx", index=False)

# Make a dataframe that will contain all the metadata fields for augur
selected_columns = df2[["Isolate_Name","Subtype","Collection_Date","Year","Host","Region","Country","Province_City","Location"]]
df3 = selected_columns.copy()
df3['strain'] = df3['Isolate_Name'].astype(str)+'_|' + df3['Subtype'].astype(str)+'|_' + df3['Collection_Date'].astype(str)
df3_columns = [col for col in df3.columns if col != 'strain']
df3_columns.insert(0, 'strain')
df3 = df3[df3_columns]
df3.insert(1, 'virus', 'Influenza')
df3.rename(columns={'Subtype': 'subtype'}, inplace=True)
df3.rename(columns={'Collection_Date': 'date'}, inplace=True)
df3.rename(columns={'Year': 'year'}, inplace=True)
df3.rename(columns={'Host': 'host'}, inplace=True)
df3.rename(columns={'Region': 'region'}, inplace=True)
df3.rename(columns={'Country': 'country'}, inplace=True)
df3.rename(columns={'Province_City': 'province_city'}, inplace=True)
df3.rename(columns={'location': 'location'}, inplace=True)
df3.drop('Isolate_Name', axis=1, inplace=True)
df3.to_csv(Metadataoutput+".cleaned.nextstrain.tsv", sep="\t", index=False)

# Find the 'Isolate_Name' references that were removed from the dataframe
df4 = pd.merge(df, df2, how='outer', indicator='Exist')
df4 = df4.loc[df4['Exist'] != 'both']
removedseq = df4['Isolate_Name'].tolist()
Removedseqfile = os.path.splitext(metadata)[0]+'.removedsequences.txt'
with open(Removedseqfile, "w") as output:
    for element in removedseq:
        output.write(element + "\n")

print("GISAID Metadata cleaned.")

