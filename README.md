# VI6processing
A collection of python scripts created to streamline processing of assembled Influenza sequences at VI6 within APHA.
Due to the specific nature of their usage, they are unlikely to work unless used on VI6 data within the VI6 storage architecture.
Could potentially be repurposed.

# MultiSplat (multisplat.py)
Requires "strains.csv" in target folder in the format "sampleid,strain\n"
When pointed to a folder containing consensus sequence files, will alter header to remove assembly artifacts and replace with appropriate strain name.
Infers name by reading to the first "_" in the filename.
Splits by segment upon completion and concatenates.
Also generates data appropriate for use in internal "Mutation Checker" tool

# GISAID Upload Format tool (gisaid_upload_format.py)
Tool for preparing data for upload to GISAID sequence repository.
GISAID requires a VERY specific format for uploads in bulk through an excel file containing metadata and sequence data.
Unlikely to work outside of VI6 due to incredibly specific usage context.

# Coverage Wizard (coverage.sh, coverage_wizardv4.py)
Ambitious tool to compare assembled sequence quality and determine best representative sample.
coverage.sh portion is a simple bash script to run samtools coverage function in bulk
python script parses coverage.sh output and compares quality using a simple weighted scoring system to determine the best representative sample, while also excluding incomplete assemblies. 
Requires "strains.csv" in the format "sampleid,samplegroup\n" to run scoring function, this will be passed over otherwise.
Works for both segmented and non-segmented sequences, though this must be input via argument -- Will default to Influenza structure --
Outputs:
- Reformatted coverage data for internal monitoring
- Contest Logic for human inspection (req. strains.csv)
- Suggested representatives by group (req. strains.csv)
