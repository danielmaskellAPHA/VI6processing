#created Nov 2022 - Dan Maskell APHA
#This script is designed to generate coverage stats for all assembled consensuses within a submission. This submission must be passed as an argument when initiating the script.
#The script will then call a python script to pick the desired data (length, coverage, depth) for all consensuses, and place it into a file that can be easily copied and pasted into the Sequence Data spreadsheet found on teams. 
#Segments that were not assembled will show as a 0 in this sheet.

segments="PB2 PB1 PA HA NP NA MP NS"

usage()
{ echo "
coverage.sh Version $VERSION
Usage: bash $0 -s 'seg1 seg2 seg3...' <folder containing assembled sequences>
	help 	Displays this usage
	
	-s 	A list of segments in the order they appear in the assembly.
		Should be contained in \"\" marks and separated by spaces.
		'full' can be used if the genome is not segmented. [DEFAULT: Influenza A segments]
	
	If you intend to compare the quality of assembled sequences, you must create a 'strains.csv' file within this folder.
	You must list a sample and its comparison group in the first 2 columns of each row.
	Make sure this is a comma separated file.
";
}

while getopts "s:f:" opt ; do
	case $opt in
		s) segments=$OPTARG;;
	esac
done

shift $((OPTIND-1))

if [[ $1 = help ]]; 
then
usage
  exit
else
submission=$1
fi

if [[ $# = 0 ]]; then
usage
  exit
fi

for bam in $submission/*_FinalResults_*/*.bam 

do 

outputfile=$(basename $bam | cut -d . -f 1) 

dir=$(dirname $bam) 

samtools coverage $bam > $dir/$outputfile.coverage.csv 

echo "created $dir/$outputfile.coverage.csv"

done 

echo "Coverage files generated."
echo "Formatting."

#change this to python script install location
python /home/$USER/Desktop/coverage_wizard/coverage_wizardv4.py $submission -s $segments

