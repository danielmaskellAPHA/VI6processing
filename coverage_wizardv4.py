# Working build created 25/11/2022 by Dan Maskell VI6 (daniel.maskell@APHA.gov.uk). 
# This script requires a submission folder as an argument. 
# If no coverage summaries exist in the Submission Folder, this script will not work.
# It is recommended to run coverage.sh, 
# as that script calls this script at its end, generating, formatting, and comparing all coverage.

import argparse
import os


def main(args):
    #folder
    argfolder = args.folder + "/"
    thepath, the_submission = os.path.split(args.folder)
    coverage_path = argfolder + the_submission + "_Coverage"
    if not os.path.isdir(coverage_path):
        os.mkdir(coverage_path)
    coverage_path += "/"
    print("Formatting coverage files for " + the_submission + "...")

    #segments
    segments = args.segments
    if "full" in segments:
        genome = True
    else:
        genome = False
    print(f"Analysing using: {segments}")
    segment_num = len(segments)

    # intermediate variables i.e., set-up function outputs
    dictionary, strains_csv_boolean = strains(input_folder=argfolder)
    files = csvs(argfolder)
    cov_write, sample_pool = parse_and_write(files, dictionary, segments, segment_num, genome)
    if strains_csv_boolean == True:
        DPRs = value_list(dictionary)
    else:
        print("strains.csv missing/blank. Only coverage summary will be produced.")
    #generate header
    # output files
    with open(
        coverage_path + the_submission + "_Coverage_Summary.csv", "w+"
    ) as new_file:
        header = "sample_name,"
        for segment in segments:
            header += f"{segment}reads,{segment}coverage,{segment}depth,"
        new_file.write(header+"\n")
        new_file.writelines(cov_write)
    print(coverage_path + the_submission + "_Coverage_Summary.csv created")

    #if using a strains file, runs contests
    if strains_csv_boolean == True:

        print("Comparing samples for each DPR")
        format_return, output_list, failed_list = stage(DPRs, dictionary, sample_pool, segments)
        with open(coverage_path + the_submission + "_Decisions.csv", "w+") as final:
                final.write(format_return)
        print(coverage_path + the_submission + "_Decisions.csv created")

        with open(coverage_path + the_submission + "_Contest_Logic.txt", "w+") as outf:
            outf.write("SUGGESTED SAMPLES ARE ONLY SUGGESTIONS. PLEASE USE YOUR OWN JUDGEMENT IF YOU DO NOT TRUST THESE CONCLUSIONS.\nIncomplete samples can be found at the bottom of this file.\n\n")
            outf.writelines("\n".join(output_list))
            if len(failed_list) != 0:
                outf.write(f"\nIncomplete samples:\n{failed_list}")
        print("Logic recorded in " + coverage_path + the_submission + "_Contest_Logic.txt")


# function to find all coverage files
def csvs(folder):
    files = [
        os.path.join(dirpath, f)
        for (dirpath, dirnames, filenames) in os.walk(folder)
        for f in filenames
    ]
    coverages = [file for file in files if file.endswith(".coverage.csv")]
    return coverages


# determine the DPR/AIVs in submission
def strains(input_folder):
    strains_boolean = os.path.isfile(input_folder + "strains.csv")
    strain_dict = {}
    if strains_boolean == True:
        with open(input_folder + "strains.csv") as strains_f:
            if strains_f.readline() == "\n":
                return strain_dict, False
            strains_f.seek(0)
            for line in strains_f:
                if "\t" in line:
                    raise CommaSeparationError()
                name, DPR = line.split(",")
                strain_dict[name] = DPR.replace("\n", "")
        return strain_dict, True
    return strain_dict, False

# generate a list of unique DPR/AIVs
def value_list(d):
    lst = []
    for DPR in d.values():
        if (
            DPR not in lst
            and DPR != "Wild"
            ):
            lst.append(DPR.replace("\n", ""))
    return lst


# stores coverage stats inside a Sample object for later use
class Sample:
    def __init__(self, name, segments, d):
        self.name = name
        self.DPR = d.get(name, "Wild")
        self.stats = segments
        gene_string = []
        for gene in self.stats:
            blank_gene = []
            for data in gene:
                blank_gene.append(str(data))
            gene_string.append(blank_gene)
        self.check = [",".join(gene) for gene in gene_string]
        if "0,0,0" in self.check:
            self.verify = False
        else:
            self.verify = True

    def __str__(self):
        return str(self.name)


def split_line_4_7(input_string):
    splitted = input_string.split("\t")
    reads = float(splitted[4])
    coverage = float(splitted[5])
    depth = float(splitted[6])
    conc = [reads, coverage, depth]
    return conc

def determine_segments():
    max_segments = 0
    

# finds coverage stats, generates our Sample objects, and returns a data structure to write to our coverage overview
def parse_and_write(files, input_dict, segments, segment_num, genome=False):
    cov_write = []
    sample_pool = []
    for thing in files:
        with open(thing) as stats:
            _path, name = os.path.split(thing)
            name_length = name.find("_")
            name2 = name[:name_length]
            _header = (stats.readline().split("\t"))[4:7]
            all_values = []
            while len(all_values) < segment_num:
                all_values.append(["0", "0", "0"])

            position = 0
            if genome == True:
                all_values[position] = split_line_4_7(stats.readline())
            else:
                for segment in segments:
                    for line in stats.readlines():
                        if segment in line:
                            all_values[position] = split_line_4_7(line)
                    stats.seek(0)
                    position += 1
            a_v_string = []
            for gene in all_values:
                blank_gene = []
                for data in gene:
                    blank_gene.append(str(data))
                a_v_string.append(blank_gene)
            a_v_printable = ",".join([",".join(gene) for gene in a_v_string])
            strings = [name2, a_v_printable]
            sample_info = [name2, all_values]
            towrite = ",".join(strings)
            cov_write.append(towrite + "\n")
            # unpack the strings list to use in the sample constructor
            new_sample = Sample(*sample_info, input_dict)
            sample_pool.append(new_sample)
    return cov_write, sample_pool

#contest of individual gene
def mini_contest(genes, point_length):
    coverage_values = []
    depth_values = []
    scores = []
    while len(scores) < point_length:
        scores.append(0)
    for item in genes:
        if len(item) == 0:
            raise TabSeparationError()
        coverage_values.append(float(item[1]))
        depth_values.append(float(item[2]))
    coverage_values = [str(i) for i in coverage_values]
    depth_values = [str(i) for i in depth_values]
    #100% coverage scores a 1
    if "100.0" in coverage_values:
        position = 0
        for item in coverage_values:
            if item == '100.0':
                scores[position] += 1
            position += 1
    #if 100% isn't in coverages, the highest scores 1
    else:
        c_highest = str(max([float(i) for i in coverage_values]))
        c_winner = coverage_values.index(str(c_highest))
        scores[c_winner] += 1
    #depth has no upper limit so highest scores 1
    d_highest = str(max([float(i) for i in depth_values]))
    d_winner = depth_values.index(str(d_highest))
    scores[d_winner] += 1
    #highest score at the end will be 2
    return scores


# using stats stored in Sample objects, this groups them by DPR and tests the stats for each gene, 
# determining which Sample object wins overall by the number of winning genes it contains.
# any samples with missing segments will be discarded
# if all samples are discarded it will be suggested that an egg isolate is grown.

def contest(DPR, strain_dict, samples, sample_names, segments):
    output_items = []
    point_array_length_goal = len(samples)
    point_array = []
    while len(point_array) < point_array_length_goal:
        point_array.append(0)
    position = 0
    for segment in segments:
        contest = []
        for sample in samples:
            contest.append(sample.stats[position])
        
        value_string = ""
        
        for participant in contest:
            value_string += f"{participant[1:]}  vs  "
        value_string = value_string[:-5]
        value_string += "\n"
        output_items.append(value_string)
        gene_points = mini_contest(contest, point_array_length_goal)
        if "2" in gene_points:
            big_win = gene_points.index(2)
            point_array[big_win] += 1
        else:
            counter = 0
            for point in gene_points:
                point_array[counter] += point
                counter += 1
        position += 1
                
    #generate a nice list of samples and their points           
    point_array_str = [str(i) for i in point_array]
    point_matrix = [": ".join(list(pair)) for pair in zip(sample_names, point_array_str)]
    point_matrix_formatted = "\n".join(point_matrix)
    output_items.append(f"Points:\n{point_matrix_formatted}\n")

    #determine actual winner
    highest_points = max(point_array)
    #check for draw(s)
    if point_array.count(highest_points) > 1:
        draws = []
        current_sample = 0
        for point in point_array:
            if point == highest_points:
                draws.append((samples[current_sample]).name)
                current_sample += 1
        format_draws = " or ".join(draws)
        DPR_winner = format_draws
    #otherwise this one wins
    else:   
        DPR_winner_position = point_array.index(highest_points)
        DPR_winner = (samples[DPR_winner_position]).name
    #clean output before returning
    output_items.append(f"Suggestion: {DPR_winner}\n")
    output = "\n".join(output_items)
    return [DPR_winner, DPR], output

#stages contests, determines which samples are valid before inserting into DPR test-pool
#returns winner(s) and final outputs at the end
def stage(DPR_input_list, strain_dict, sample_pool, segments):
    output_list = []
    selected = {}
    failed = {}
    for DPR in DPR_input_list:
        output_list.append(f"\n{DPR}\n")
        samples = []
        failed_count = 0
        for sample in sample_pool:
            if DPR == strain_dict.get(sample.name):
                if sample.verify is True:
                    samples.append(sample)
                else:
                    failed[sample.name] = DPR
                    failed_count += 1
        print(f"{DPR}: {len(samples)} complete samples. ({failed_count} incomplete.)")         
        sample_names = []
        for sample in samples:
            sample_names.append(sample.name)
        output_list.append("  vs  ".join(sample_names))

        #testing for solitary sample in DPR, will avoid testing PROVIDED it is valid.
        if len(samples) == 1:
            output_list.append(f"\n{samples[0]} was the only valid sample\n")
            genes = [[str(i) for i in gene[1:]] for gene in samples[0].stats]
            for gene in genes:
                output_list.append(str(gene))
            selected[DPR] = samples[0].name
            output_list.append(f"\nSuggestion: {samples[0]}\n")
        
        #no valid samples
        if len(samples) == 0:
            output_list.append(f"\n{DPR} FAILURE: Grow ISOLATE\n")
            selected[DPR] = "Grow Isolate"

        if len(samples) > 1:
            winner, contest_output = contest(DPR, strain_dict, samples, sample_names, segments)
            output_list.append(contest_output)
            selected[winner[1]] = winner[0]

    #output formatting    
    to_return = []
    for k, v in selected.items():
        thing = str(v + "," + k)
        to_return.append(thing)
    format_return = "\n".join(to_return)
    
    to_return_fails = []
    for k, v in failed.items():
        thing = f"{k} ({v})"
        to_return_fails.append(thing)
    failed_format = "\n".join(to_return_fails)
    
    return format_return, output_list, failed_format

class Error(Exception):
        #Base class for other exceptions
        pass
    
class TabSeparationError(Error):
        #Raised if provided argument is not a directory
        def __init__(self, message = "This script is written to expect tab separated coverage files (csv files can, confusingly, be tab separated) as this is the typical output of samtools.\nPlease ensure all coverage files are tab separated."):
                self.message = message
                super().__init__(self.message)
                
class CommaSeparationError(Error):
        #Raised if provided argument is not a directory
        def __init__(self, message = "This script is written to expect a comma separated strains file (csv files can, confusingly, be tab separated).\nPlease ensure strains.csv is comma separated - or delete it if you do not need comparisons."):
                self.message = message
                super().__init__(self.message) 

if __name__ == "__main__":
    # parse argument (will grab from coverage.sh if used)
    parser = argparse.ArgumentParser()
    parser.add_argument("folder", help="directory containing coverage")
    parser.add_argument("-s", "--segments", help="input segment names separated by a space", action = "store", nargs = "*", default = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"])
    input_args = parser.parse_args()
    main(input_args)
