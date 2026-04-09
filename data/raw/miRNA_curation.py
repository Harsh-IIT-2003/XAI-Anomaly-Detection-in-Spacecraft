#!/usr/bin/env python
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import regex
import sys
from collections import OrderedDict
sys.path.append("/usr/local/lib/python3.7/site-packages/RNA")
import RNA
import os
import matplotlib.pyplot as plt
import numpy as np

def parse_fastq(in_file, collect, library, change_TU):
    """Collect all reads from fasta file in a dictionary together with read count."""
    with open(in_file) as fin:
        print("Library ID = "+str(library))
        for line in fin:
            if line[0] == "@":
                read = next(fin, '').strip("\n")
                if change_TU == 1:
                    read = read.replace("T","U")
                if read not in collect:
                    collect[read] = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, "total": 0}
                    collect[read][library] += 1
                    collect[read]["total"] += 1
                else:
                    collect[read][library] += 1
                    collect[read]["total"] += 1
    return collect


def map_reads(parsed_reads, ref_seq, rc, map_mismatch):
    """ Map reads that match perfectly or with one mismatch """
    mapped = {}
    for read in parsed_reads:
        if read in ref_seq:
            mapped[read] = parsed_reads[read]
            mapped[read]["pos"] = ref_seq.find(read)
            mapped[read]["mis"] = 0
        elif map_mismatch == 1:
            reg = "(" + read + "){s<=1}"
            mismatch = regex.findall(reg, ref_seq)
            if mismatch:
                #Add mismatch location to dictionary
                mindex = (regex.search(reg, ref_seq).fuzzy_changes[0])
                mapped[read] = parsed_reads[read]
                mapped[read]["pos"] = ref_seq.find(mismatch[0])
                if mindex != []:
                    mapped[read]["mis"] = mindex[0]
                else:
                    mapped[read]["mis"] = 0
        elif rc == 1:
            rread = Seq(read.replace("U","T"))
            rcread = str(rread.reverse_complement()).replace("T","U")
            if rcread in ref_seq:
                mapped[read] = parsed_reads[read]
                mapped[read]["pos"] = ref_seq.find(rcread)
                mapped[read]["mis"] = -1
    return mapped


def explore_map(ss, pri_mir, mapped):
    """ Identify the 5' and 3' positions of the highest expressed mapped sRNA on the hairpin and find the opposing positions in the hairpin with 2nt overhang """
    positions = [0,0,0,0]
    highest_mapped = list(OrderedDict(sorted(mapped.items(), key=lambda x: x[1]["total"], reverse=True)).keys())[0]
    positions[0] = mapped[highest_mapped]["pos"]
    positions[1] = mapped[highest_mapped]["pos"]+len(highest_mapped)-1
    counter = 0
    complement = []
    avg = 0
    for i in ss:
        if i == "(":
            counter += 1
            complement.append(counter)
        elif i == ")":
            complement.append(counter)
            counter -= 1
        else:
            complement.append(counter)

    for nucleotide in ["fivep","threep"]:  ##Find the opposing nucleotides with 2nt overhang
        if nucleotide == "fivep":
            index = positions[0]
            comp_value = max(complement[positions[0]],1) #set value of 5'nucleotide to match (minimum 1 because otherwise it will never drop below)
        elif nucleotide == "threep":
            index = positions[1]-2
            comp_value = max(complement[positions[1]-2],1) #set value of 3'nucleotide to match (minimum 1 because otherwise it will never drop below)
        mm_start = 0
        while ss[index] == ".": #if nucleotide is in unmatched region, slide left
            if index == 0: #if it slides to the end of the pri-miRNA break
                break
            mm_start += 1
            index -= 1
        if ss[index] == "(": #if the 3'nucleotide is (, look forward into the hairpin until the value would drop below match
            while complement[index+1] >= comp_value:
                index += 1
                if index == len(ss)-1:
                    break
            while mm_start > 0 and complement[index-1] == comp_value:
                index -= 1
                mm_start -= 1
        elif ss[index] == ")": #if the 3'nucleotide is ), look backward into the hairpin until the value would drop below match
            while complement[index-1] >= comp_value:
                index -= 1
                if index == 0:
                    break
            if mm_start > 0:
                while complement[index+1] == comp_value:
                    index +=1
                while mm_start > 1 and complement[index-1] == comp_value:
                    index -= 1
                    mm_start -= 1
        if nucleotide == "fivep":
            positions[3] = min(index+2,len(ss)-1)
        elif nucleotide == "threep":
            positions[2] = index
    positions = sorted(positions)
    print(positions)
    cased_pri_mir = pri_mir[:positions[0]]+pri_mir[positions[0]].lower()+pri_mir[positions[0]+1:positions[1]]+pri_mir[positions[1]].lower()+pri_mir[positions[1]+1:positions[2]]+pri_mir[positions[2]].lower()+pri_mir[positions[2]+1:positions[3]]+pri_mir[positions[3]].lower()+pri_mir[positions[3]+1:]
    return cased_pri_mir, positions



def precision(ord_mapped, positions, library):
    """ Calculate the miRNA precision in the hairpin in different methods according to Axtell 2018, Fromm 2022 and Kozomara 2014 & 2019 """
    axtell5p = 0
    axtell3p = 0
    axtell_total = 0
    fivep = 0
    threep = 0
    fivearm = 0
    threearm = 0
    axtell_percentage = 0
    fivepercentage = 0
    threepercentage = 0
    total_percentage = 0
    for read in ord_mapped:
        if ord_mapped[read]["mis"] >= 0:           
            if positions[0]-1 <= ord_mapped[read]["pos"] <= positions[0]+1 and positions[1]-1 <= ord_mapped[read]["pos"]+len(read)-1 <= positions[1]+1:
                axtell5p += ord_mapped[read][i]
            elif positions[-2]-1 <= ord_mapped[read]["pos"] <= positions[-2]+1 and positions[-1]-1 <= ord_mapped[read]["pos"]+len(read)-1 <= positions[-1]+1:
                axtell3p += ord_mapped[read][i]
            elif len(positions) == 8:
                if positions[2]-1 <= ord_mapped[read]["pos"] <= positions[2]+1 and positions[3]-1 <= ord_mapped[read]["pos"]+len(read)-1 <= positions[3]+1:
                    axtell5p += ord_mapped[read][i]
                elif positions[4]-1 <= ord_mapped[read]["pos"] <= positions[4]+1 and positions[5]-1 <= ord_mapped[read]["pos"]+len(read)-1 <= positions[5]+1:
                    axtell3p += ord_mapped[read][i]
            if positions[0] == ord_mapped[read]["pos"] and positions[1]-5 <= ord_mapped[read]["pos"]+len(read)-1 <= positions[1]+2:
                fivep += ord_mapped[read][i]
            elif positions[-2] == ord_mapped[read]["pos"] and positions[-1]-5 <= ord_mapped[read]["pos"]+len(read)-1 <= positions[-1]+2:
                threep += ord_mapped[read][i]
            elif len(positions) == 8:
                if positions[2] == ord_mapped[read]["pos"] and positions[3]-5 <= ord_mapped[read]["pos"]+len(read)-1 <= positions[3]+2:
                    fivep += ord_mapped[read][i]
                elif positions[4] == ord_mapped[read]["pos"] and positions[5]-5 <= ord_mapped[read]["pos"]+len(read)-1 <= positions[5]+2:
                    threep += ord_mapped[read][i]

        if (positions[int(len(positions)/2)-1]+positions[int(len(positions)/2)])/2 >= (ord_mapped[read]["pos"]+len(read)/2):
            fivearm += ord_mapped[read][i]
        else:
            threearm += ord_mapped[read][i]
        axtell_total += ord_mapped[read][i]
    if (axtell5p + axtell3p) > 0:
        axtell_percentage = (axtell5p + axtell3p) * 100 / axtell_total
    if fivep > 0:
        fivepercentage = fivep * 100 / fivearm
    if threep > 0:
        threepercentage = threep * 100 / threearm
    if (threep + fivep) > 0:
        total_percentage = (threep + fivep) * 100 / (threearm + fivearm)
    return axtell_percentage, fivepercentage, threepercentage, total_percentage, axtell5p, axtell3p, fivep, threep, axtell_total, fivearm, threearm

def pairing(ss, positions):
    """ Analyze the hairpin to identify the number of paired, unpaired, and bulged nucleotides """
    hairpin = 0
    minloop = 5
    paired = 0
    unpaired = 0
    overhang_unpaired = 0
    bulge = 0
    count = 0
    count_star = 0
    if positions[-1]>positions[-2]+17 and positions[-2]>positions[1]+minloop and positions[1]>positions[0]+17:
        if ss[positions[0]:positions[1]].count(")") == 0 and ss[positions[-2]:positions[-1]].count("(") == 0:
            hairpin = 1

    if hairpin == 1:
        while count <= (positions[int(len(positions)/2-1)] - positions[0] + 2):
            if count <= (positions[int(len(positions)/2-1)] - positions[0]) and count > 1: #These are the bases pairing in the duplex
                if ss[positions[0]-2+count] == '(' and ss[positions[-1]-count_star] == ')':
                    paired += 1
                    count += 1
                    count_star += 1
                elif ss[positions[0]-2+count] == '.' and ss[positions[-1]-count_star] == '.':
                    unpaired += 1
                    count += 1
                    count_star += 1
                elif ss[positions[0]-2+count] == '(' and ss[positions[-1]-count_star] == '.':
                    bulge += 1
                    count_star += 1
                elif ss[positions[0]-2+count] == '.' and ss[positions[-1]-count_star] == ')':
                    bulge += 1
                    count += 1
                elif ss[positions[0]-2+count] == ')' or ss[positions[-1]-count_star] == '(':
                    break
            else: #These are the 3' bases overhanging, which are only counted for paired in mirgenedb
                if ss[positions[0]-2+count] == '(' and ss[positions[-1]-count_star] == ')':
                    paired += 1
                    count += 1
                    count_star += 1
                else:
                    overhang_unpaired += 1
                    count += 1
                    count_star += 1
            #print(paired, unpaired, bulge, count, count_star)
        print(ss[positions[0]-2:positions[int(len(positions)/2-1)]+1])
        print(ss[positions[-1]:positions[int(len(positions)/2)]-3:-1])
        print(paired, unpaired, bulge)
    return(hairpin, paired, unpaired, bulge, overhang_unpaired)

def rescale(mapped, positions, cased_pri_mir, extra_bases):
    """ Trim the hairpin to contain only the number of bases defined in extra_bases, leading up to the miRNA duplex """
    frame = [max(min(positions) - extra_bases,0),min(max(positions) + extra_bases-1,len(cased_pri_mir))]
    print(positions, frame)
    scaled_mapped = {}
    for read in mapped:
        if  mapped[read]["pos"] >= frame[0] and mapped[read]["pos"]+len(read) <= frame[1]:
            scaled_mapped[read] = mapped[read]
            scaled_mapped[read]["pos"] -= frame[0]
            if scaled_mapped[read]["mis"] > 0:
                scaled_mapped[read]["mis"] -= frame[0]
    positions = [x - frame[0] for x in positions]
    cased_pri_mir = cased_pri_mir[frame[0]:frame[1]]
    return scaled_mapped, positions, cased_pri_mir, frame


def map_density(cased_pri_mir, ord_mapped):
    """ Quantify the number of reads mapping to each base in the hairpin for plotting of the density histogram """
    den_lib1 = [0] * len(cased_pri_mir)
    den_lib2 = [0] * len(cased_pri_mir)
    den_lib3 = [0] * len(cased_pri_mir)
    for read in ord_mapped:
        if ord_mapped[read]["mis"] >= 0:
            pos = ord_mapped[read]["pos"]
            den_lib1[pos:(pos+len(read))] = [den_lib1[i]+ord_mapped[read][1] for i in range(pos, (pos+len(read)))]
            den_lib2[pos:(pos+len(read))] = [den_lib2[i]+ord_mapped[read][2] for i in range(pos, (pos+len(read)))]
            den_lib3[pos:(pos+len(read))] = [den_lib3[i]+ord_mapped[read][3] for i in range(pos, (pos+len(read)))]
    return den_lib1, den_lib2 ,den_lib3

def create_output_line(read, pos, len_ref, mis):
    """ Helper function for outputting the reads in the correct position below the pri-miRNA with bracket notation """
    start = "-" * pos
    end = "-" * (len_ref - (pos + len(read)))
    if mis == 0:
        return start + read + end
    elif mis == -1:
        start = "<" * pos
        end = "<" * (len_ref - (pos + len(read)))
        return start + read[::-1] + end
    else:
        #Make the mismatching nucleotide lowercase
        m = str(start + read + end)
        if mis == len(m):
            m = m[:mis]+m[mis].lower()
        elif mis < len(m):
            m = m[:mis]+m[mis].lower()+m[(mis+1):]
        return m

def fold_miRNA(ss, cased_pri_mir):
    """ Create a rough 2D miRNA structure to more easily draw the 2D structures """
    folded_miRNA = ""
    half_length = int(len(ss)/2)
    for c in [".","("]:
        indexes = ([pos for pos, char in enumerate(ss[0:half_length]) if char == c])
        for i in range(half_length+1):
            if i in indexes:
                folded_miRNA += cased_pri_mir[i]
            else:
                folded_miRNA += " "
        folded_miRNA += "\n"
    for c in [")","."]:
        indexes = ([pos for pos, char in enumerate(ss[len(ss):half_length-1:-1]) if char == c])
        for i in range(len(ss)+1-half_length):
            if i in indexes:
                folded_miRNA += cased_pri_mir[len(ss):half_length-1:-1][i]
            else:
                folded_miRNA += " "
        folded_miRNA += "\n"
    return folded_miRNA



parser = argparse.ArgumentParser()
parser.add_argument("-f","--fasta",  required=True, help="A multiline fasta file containing the miRNA candidate hairpin sequences")
parser.add_argument("-i","--in_files", nargs='+',  required=True, help="Input fastq sRNA files to be mapped")
parser.add_argument("-o","--out_folder", required=True, help="Destination folder for the output files")
parser.add_argument("-m","--mode", required=True, help="miRNA curation running modes (default : refine); \n \
    \texplore: Fold and analyze the hairpin, identify where the miRNA duplex is located, and calculate the miRNA precision \n \
    \trefine: After running explore, refine the miRNA analysis and add density plot and 2D structure of the miRNA hairpin to output")
parser.add_argument("-e","--extra", required=True, help="Number of nucleotides added to the 5' end of the miRNA-5p and the 3' end of the miRNA-3p")
parser.add_argument("-id","--identifier", nargs='?', const="", default="", help="Identifier of the run added to the output table for simplifying combined analysis afterwards")
parser.add_argument("-s","--substitutions", default="0", help="Number of substitutions allowed when mapping the sRNA reads to the miRNA-hairpin (0 or 1, default : 0)")
parser.add_argument("-rc","--reversecomplement", default="0", help="Also map reads that are mapping on the opposite DNA strand on the miRNA-hairpin. \n \
    Mainly to avoid wrongly annotating siRNAs as miRNAs which often map sRNAs from longer, precisely paired dsRNA")
args = parser.parse_args()

in_files = args.in_files # single fastq file or multiple files separated with commas
fasta_file = args.fasta # fasta file containing one or multiple 
out_folder = args.out_folder # folder to output the files into
map_mismatch = int(args.substitutions) #use fuzzy_regex to map reads with given mismatch
extra_bases = int(args.extra)
mode = args.mode
if mode != "explore":
    mode = "refine"
rc = int(args.reversecomplement) #also map reverse complement reads (1) or not (0)
change_TU = 1 #Output with Uracil instead of Thymine


if args.identifier:
    identifier = str(args.identifier)
else:
    identifier = ""

##Generate output folder
if not os.path.exists(out_folder):
    os.makedirs(out_folder)

##Read fastq file into dict parsed_reads
parsed_reads = {}
folded_struc = {}
summary = {}
lib = 0
reads = 0
print(in_files)
for file in in_files:
    print("Parsing {}...".format(file))
    lib += 1
    parsed_reads = parse_fastq(file, parsed_reads, lib, change_TU)
    print("done")


##Cycle through fasta sequences to map reads onto
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
for fasta in fasta_sequences:

    name, pri_mir = fasta.id, str(fasta.seq).upper()
    len_pri = len(pri_mir)  
    if change_TU == 1:
        pri_mir = pri_mir.replace("T","U")


    print("Mapping parsed reads to "+name)
    mapped = map_reads(parsed_reads, pri_mir, rc, map_mismatch)
    md = RNA.md()
    md.temperature = 22.0 #Dictyostelium discoideum and other amoebae were grown at 22°C
    (ss, mfe) = RNA.fold_compound(pri_mir, md).mfe()

    if mode == "explore":
        if len(mapped) != 0:
            cased_pri_mir, positions = explore_map(ss, pri_mir, mapped) #cased_pri_mir has lowercase nucleotides to identify the miR-5p 5' and 3' ends and miR-3p 5' and 3' ends
        else:
            cased_pri_mir = pri_mir
            positions = [0,0,0,0]
    elif mode == "refine":
        cased_pri_mir = str(fasta.seq)
        positions = [i for i in range(len_pri) if fasta.seq[i].islower()]

    hairpin, paired, unpaired, bulge, overhang_unpaired = pairing(ss, positions) 
    mapped, positions, cased_pri_mir, frame = rescale(mapped, positions, cased_pri_mir, extra_bases)
    (ss, mfe) = RNA.fold_compound(cased_pri_mir, md).mfe()
    
  
    print("Creating output..")
    axtell5p_reads = ""
    axtell3p_reads = ""
    fivep_reads = ""
    threep_reads = ""
    axtell_total_reads = ""
    fivep_total_reads = ""
    threep_total_reads = ""
    reads_list = [0,0,0,0] #index 0 = miR-5p reads, index 1 = miR-3p reads, index 2 = all 5p-arm reads, index 3 = all 3p-arm reads
    axtell_percentage_list = [0,0,0,0]
    fivepercentage_list = [0,0,0,0]
    threepercentage_list = [0,0,0,0]
    total_percentage_list = [0,0,0,0]
    ord_mapped = OrderedDict(sorted(mapped.items(), key=lambda x: x[1]["pos"], reverse=False))
    ord_mapped = OrderedDict(sorted(ord_mapped.items(), key=lambda x: x[1][1], reverse=True))
    mir_seq = next(iter(ord_mapped),"")
    with open((out_folder+name+".txt"), "w") as fout:
        fout.write(ss+"\t"+str(round(mfe,2))+" kcal/mol\n")
        fout.write("{}\tcount\tlength\n".format(cased_pri_mir.replace("T","U").replace("t","u")))
        for i in range(1,lib+1):
            #Loop through the libraries to get the precision in each library, saved as strings with tabs for output table, and ints for calculation in python
            axtell_percentage_list[i], fivepercentage_list[i], threepercentage_list[i], total_percentage_list[i], axtell5p, axtell3p, fivep, threep, axtell_total, fivearm, threearm = precision(ord_mapped, positions, i)
            axtell5p_reads = axtell5p_reads + str(axtell5p) + "\t"
            axtell3p_reads = axtell3p_reads + str(axtell3p) + "\t"
            fivep_reads = fivep_reads + str(fivep) + "\t"
            reads_list[0] += fivep
            threep_reads = threep_reads + str(threep) + "\t"
            reads_list[1] += threep
            axtell_total_reads = axtell_total_reads + str(axtell_total) + "\t"
            fivep_total_reads = fivep_total_reads + str(fivearm) + "\t"
            reads_list[2] += fivearm
            threep_total_reads = threep_total_reads + str(threearm) + "\t"
            reads_list[3] += threearm
        if reads_list[2] > 0 and reads_list[3] > 0:
            percentage_precision = (reads_list[0]+reads_list[1])*100/(reads_list[2]+reads_list[3])
            fout.write("Libraries combined "+"\t Precision: [ total = "+str(round(percentage_precision,1))+"% | 5p-arm = "+str(round(reads_list[0]*100/reads_list[2],1))+"% | 3p-arm = "+str(round(reads_list[1]*100/reads_list[3],1))+"% ]\n")
        else:
            percentage_precision = 0
        for read in ord_mapped:
            fout.write("{}\t{}\t{}\n".format(create_output_line(read, ord_mapped[read]["pos"], frame[1]-frame[0], ord_mapped[read]["mis"]),
                                                         (ord_mapped[read][1]+ord_mapped[read][2]+ord_mapped[read][3]),
                                                         len(read)))
        for i in range(1,lib+1):
            ord_mapped = OrderedDict(sorted(mapped.items(), key=lambda x: x[1]["pos"], reverse=False))
            ord_mapped = OrderedDict(sorted(ord_mapped.items(), key=lambda x: x[1][i], reverse=True))
            fout.write("Library = "+str(i)+"\t Precision: [ total = "+str(round(total_percentage_list[i],1))+"% | 5p-arm = "+str(round(fivepercentage_list[i],1))+"% | 3p-arm = "+str(round(threepercentage_list[i],1))+"% ]\n")
            for read in ord_mapped:
                if ord_mapped[read][i]>0:
                    fout.write("{}\t{}\t{}\n".format(create_output_line(read, ord_mapped[read]["pos"], frame[1]-frame[0], ord_mapped[read]["mis"]),
                                                         ord_mapped[read][i],
                                                         len(read)))

        if mode == "refine":
            #Add 2D structure
            folded_miRNA = fold_miRNA(ss, cased_pri_mir)
            fout.write("\n"+folded_miRNA)

    if mode == "refine":
        #Plot density of miRNA mapping on miRNA hairpin
        den_lib1, den_lib2, den_lib3 = map_density(cased_pri_mir, ord_mapped)
        x = np.arange(len(den_lib1))
        plot_lib2 = np.add(den_lib2, den_lib3)
        plot_lib1 = np.add(den_lib1, plot_lib2)
        plt.figure(figsize=(15,1.5))
        plt.fill_between(x,plot_lib1, color="#AAAAAA", step="mid", label="library 1")
        plt.fill_between(x,plot_lib2, color="#888888", step="mid", label="library 2")
        plt.fill_between(x,den_lib3, color="#666666", step="mid", label="library 3")
        ticks = [ss[i]+"\n"+cased_pri_mir.upper()[i] for i in range(len(den_lib1))]
        plt.xticks(np.arange(len(den_lib1)), ticks)
        if len(positions) == 8:
            plt.vlines(x=np.add(positions,[-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5]), ymin=0, ymax=max(plot_lib1), colors='black', ls=':', lw=1)
        else:
            plt.vlines(x=np.add(positions,[-0.5,0.5,-0.5,0.5]), ymin=0, ymax=max(plot_lib1), colors='black', ls=':', lw=1)
        colors = ["black"] * positions[0] + ["#ED1C24"] * (positions[int(len(positions)/2-1)]+1-positions[0]) + ["black"] * (positions[int(len(positions)/2)]-1-positions[int(len(positions)/2-1)]) + ["#00AEEF"] * (positions[-1]+1-positions[int(len(positions)/2)]) + ["black"] * (len(cased_pri_mir)-1-positions[-1])
        for ticklabel, tickcolor in zip(plt.gca().get_xticklabels(), colors):
            ticklabel.set_color(tickcolor)
        plt.tight_layout()
        plt.box(False)
        #plt.legend()
        plt.savefig(out_folder+name+".png", dpi = 600)
        plt.clf()
        plt.close()

    folded_struc[name] = {"pri_mir": cased_pri_mir.upper(),"ss": ss, "positions": positions}
    if sum(positions) != 0 and (((positions[int(len(positions)/2-1)]-positions[0])+(positions[-1]-positions[int(len(positions)/2)]))/2) > 0:
        percentage_paired = paired/(((positions[int(len(positions)/2-1)]-positions[0]+2)+(positions[-1]-positions[int(len(positions)/2)])+2)/2)
    else:
        percentage_paired = 0

    summary[name] = {"pri_mir": cased_pri_mir,"ss": ss, "hairpin": hairpin, "positions": positions, "bulges": bulge, "paired_bases": paired, "unpaired_bases": unpaired,  "overhang_unpaired_bases": overhang_unpaired, 
    "percentage_paired": percentage_paired, "energy": str(round(mfe,2)), "mir_seq": mir_seq,
    "axtell5p_reads":axtell5p_reads, "axtell3p_reads":axtell3p_reads, "fivep_reads": fivep_reads, "threep_reads":threep_reads, "axtell_total_reads":axtell_total_reads, "fivep_total_reads":fivep_total_reads, "threep_total_reads":threep_total_reads, "precision": percentage_precision}


with open(out_folder+"miRNA_struc_density.txt", "w") as sout:
    for name in folded_struc:
        sout.write("{}\n{}\n{}\n".format((">"+name),folded_struc[name]["ss"],folded_struc[name]["pri_mir"]))

    for name in folded_struc:
        sout.write(">"+name+"\n"+
        str(folded_struc[name]["positions"][0]+1)+"-"+str(folded_struc[name]["positions"][int(len(positions)/2-1)]+1)+":#D9938D "+
        str(folded_struc[name]["positions"][int(len(positions)/2)]+1)+"-"+str(folded_struc[name]["positions"][-1]+1)+":#92AFC7\n")

with open(out_folder+"combined_analysis.tsv", "w") as cout:
    cout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{}{}{}{}{}{}\n".format("species","miRNA-ID","pri-miRNA","ss","hairpin","positions","bulges","paired_bases","unpaired_bases","overhang_unpaired_bases","percentage_paired","energy","mir_seq",
        "lib1_axtell5p_reads\tlib2_axtell5p_reads\tlib3_axtell5p_reads\t","lib1_axtell3p_reads\tlib2_axtell3p_reads\tlib3_axtell3p_reads\t","lib1_fivep_reads\tlib2_fivep_reads\tlib3_fivep_reads\t","lib1_threep_reads\tlib2_threep_reads\tlib3_threep_reads\t",
        "lib1_axtell_total_reads\tlib2_axtell_total_reads\tlib3_axtell_total_reads\t","lib1_fivep_total_reads\tlib2_fivep_total_reads\tlib3_fivep_total_reads\t","lib1_threep_total_reads\tlib2_threep_total_reads\tlib3_threep_total_reads\t",))
    for name in summary:
        cout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{}{}{}{}{}{}\n".format(identifier,identifier+"_"+name,summary[name]["pri_mir"],
            summary[name]["ss"],summary[name]["hairpin"],summary[name]["positions"],summary[name]["bulges"],summary[name]["paired_bases"],summary[name]["unpaired_bases"],summary[name]["overhang_unpaired_bases"],summary[name]["percentage_paired"],summary[name]["energy"],summary[name]["mir_seq"],
            summary[name]["axtell5p_reads"],summary[name]["axtell3p_reads"],summary[name]["fivep_reads"],summary[name]["threep_reads"], summary[name]["axtell_total_reads"], summary[name]["fivep_total_reads"], summary[name]["threep_total_reads"]))

with open(out_folder+"pri_miRNA.fa", "w") as fasta_out:
    for name in summary:
        if summary[name]["hairpin"] == 1 and (float(summary[name]["energy"])/len(summary[name]["pri_mir"])) < -0.2 and summary[name]["precision"] > 30:
            fasta_out.write(">"+name+"\n"+
                summary[name]["pri_mir"]+"\n")



print("done")
