#updated 16/11/22 Dan Maskell daniel.maskell@APHA.gov.uk
#this script renames (.fasta headers feature reference sequence names as an artifact of assembly pipeline) consensus sequence files by parsing a pre-made "strains.csv" file in the format "sampleid,strain-name\n"
#this script is unlikely to work outside of the intended use within VI6
#upon renaming, the script generates a concatenated file for each segment containing all samples processed, for alignment, phylogenetics, etc. - as influenza is a segmented virus
#also generates similar files in a format appropriate for use with the internal VI6 "Mutation Checker"
#no usage created due to the existence of internal bench notes - errors will raise to guide user upon incorrect use
#updated 09/01/23 readability

#imports
from Bio import Seq, SeqIO
import argparse
import os
import sys
from os.path import isfile, isdir, join

#main function
def main(folder):
        #organising path, check for errors, print to terminal
        _thepath, the_submission = os.path.split(folder)
        folder += "/"
        error_check(folder, the_submission)
        print(f"Renaming and splitting {the_submission}...\n")
        
        #now generating dictionary and list of samples	
        dictionary = Dictionary(folder)
        files = file_spew(folder)
        samples = sample_list(files)
        
        #output folder name
        new_location = folder+the_submission+"_Final_Consensus_Sequences/"
        
        #mega-file name
        catfile = new_location+the_submission+".fasta"
        
        #creates output folder
        if os.path.exists(new_location) == False:
                os.makedirs(new_location)
        if os.path.exists(new_location+"Mutation_Checker_Genes/") == False:
                os.makedirs(new_location+"Mutation_Checker_Genes/")

        #how many samples were found in folder - not targeted by csv - allow user to check expected output
        print(f"{len(samples)} sample(s) located in {folder}:\n")
        for sample in samples:
                print(sample)
        
        #opening the megafile, as this is w+ whatever existed previously will be overwritten - however shouldn't run if previous exists
        catfile_write = open(catfile, "w+")
        
        print(f"\nCreated {the_submission}.fasta")

        #segment list
        segments = [
                "|PB2",
                "|PB1",
                "|PA",
                "|HA",
                "|NP",
                "|NA",
                "|MP",
                "|NS"
                ]

        #generates split output files
        for seg in segments:
                to_open = new_location+f"{seg[1:]}.fasta"
                opened = open(to_open, "w+")
                opened.close()
                print(f"Created {seg[1:]}.fasta")
                
        #sorts alphanumerically to hopefully put samples in chronological order
        files_sort = sorted(files)
        
        #now puts them in reverse so that the newest sample appears at the top
        files_r = reversed(files_sort)
        
        #where the magic happens
        with open(catfile, "w+") as catfile_write:
                for file_input in files_r:
                        splat(file_input, dictionary, new_location, segments, catfile_write)

        #done!
        print("Please double check output files to ensure everything is as expected.")
                
        
#errors                
def error_check(argfolder_input, submission_input):
        #directory error
        if not isdir(argfolder_input):
                raise NotDirError()
        #potential wrong folder error
        if not ("NGS" in argfolder_input):
                raise BadFolderError()
        #already run error
        if isfile(argfolder_input+submission_input+"_Final_Consensus_Sequences/"+submission_input+".fasta"):
                raise AlreadyRunError()
        #no strains.csv error	
        if not isfile(argfolder_input+"strains.csv"):
                raise CSVError()
                
	
#this function will find consensus files within a directory
def file_spew(argfolder_input):
	files = [
        os.path.join(dirpath, f)
        for (dirpath, dirnames, filenames) in os.walk(argfolder_input)
        for f in filenames
        ]
	
	consensus = [file for file in files if isfile(file) and file.endswith("iter4_consensus.fasta")]
	return consensus


#ensures csv will be parsed correctly permitting common user mistakes
def line_cleaner(line_input):
        (sample, strain) = line_input.split(",")
        clean_sample = sample.strip('"')
        clean_strain = strain.strip("\n").strip('"')
        return clean_sample, clean_strain

#generates a sample:strain dictionary
def Dictionary(argfolder_input):
	dictionary_output={}
	with open(argfolder_input+"strains.csv") as f:
		for line in f:
			sample, strain = line_cleaner(line)
			dictionary_output[sample] = strain
	return dictionary_output


#primary function	
def splat(file_input, dictionary_input, location_input, segments_input,
	catfile_input):
	#acquires file name from path _path
	_path, sample = os.path.split(file_input)
	#determines output file name
	new_name = sample.replace(".fasta", ".rename.fasta")
	#determines the new absolute path
	new_path = location_input+new_name
	#finds the end of the sample (up to first _)
	name_end = sample.find("_")
	sample_id = sample[:name_end]
	strain = dictionary_input.get(sample_id, "")
	if strain == "":
		return
	#list not yet traversed
	#creates output files
	with open(file_input) as original, open(new_path, "w+") as output:
		current_seg = 0
		records = SeqIO.parse(original, "fasta")
		#iterates through them
		for record in records:
		#determines the current gene for later use
			current_gene = segments_input[current_seg]
			#will append to the output file
			with open(location_input+f"{current_gene[1:]}.fasta", "a+") as gene_path, open(location_input+f"Mutation_Checker_Genes/{current_gene[1:]}_mc.fasta", "a+") as gene_path_mc:
				#this SHOULD halt the function if the strain name is already in the file
				check_id = record.id.find(strain)
				if check_id != -1:
					print("It appears multisplat.py has already been run on this file.")
				#otherwise it runs
				else:
				#deletes everything from the id/description
					record.id = ""
					record.description = ""
					#sets the id to be the strain name - no gene
					record.id = strain
					strain_mc = strain.replace("|", "-")
					#inserts it into the output file
					SeqIO.write(record, gene_path, "fasta")
					#now generates id WITH gene
					record.id = strain+segments_input[current_seg]
					#writes to the rename file
					SeqIO.write(record, output, "fasta")
					#writes to the concatenated mega-file
					SeqIO.write(record, catfile_input, "fasta")
					record.id = strain_mc+segments_input[current_seg]
					SeqIO.write(record, gene_path_mc, "fasta")
					#move to the next segment
			current_seg += 1
			

def sample_list(files_input):
        samples_output = []
        for submission in files_input:
                path, sample = os.path.split(submission)
                name_end = sample.find("_")
                samples_output.append(sample[:name_end])

        if len(samples_output) == 0:
                raise NoSamplesError()
                
        return samples_output

#custom exceptions
class Error(Exception):
        #Base class for other exceptions
        pass
class NotDirError(Error):
        #Raised if provided argument is not a directory
        def __init__(self, message = "Provided argument is not a directory."):
                self.message = message
                super().__init__(self.message)               
class BadFolderError(Error):
        #Raised if provided argument is probably not a submission folder
        def __init__(self, message = "Provided directory does not appear to be a submission folder. (Missing 'NGS')"):
                self.message = message
                super().__init__(self.message)           
class AlreadyRunError(Error):
        #Raised if output files already exist
        def __init__(self, message = "Output files already present."
                     "Please delete output folder before running again."):
                self.message = message
                super().__init__(self.message)
class CSVError(Error):
        #Raised if strains.csv cannot be located
        def __init__(self, message = ("strains.csv could not be located in the submission folder.\n"
                                      "Please create a comma separated (.csv) file containing the list of samples in the first column, and a list of strains in the second."
                                      )):
                self.message = message
                super().__init__(self.message)
class NoSamplesError(Error):
        #Raised if no samples are located
        def __init__(self, message = "No samples were found in the submission folder."):
                self.message = message
                super().__init__(self.message)


if __name__ == "__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("folder", help="directory containing consensuses")
        args = parser.parse_args()
        main(args.folder)
			


