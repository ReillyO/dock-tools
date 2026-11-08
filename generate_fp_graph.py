import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from collections import defaultdict

def parse_args():
    program_desc="This script is used to plot a set of footprints with varying types out of output."
    parser=argparse.ArgumentParser(program_desc)
    parser.add_argument("-p", "--plot_type", choices=['individual','overlaid'], \
                        help="Type of plotting you would like to do. Defaults to 'individual'", \
                        default='individual')
    parser.add_argument("-n", "--molnames", \
                        help="Filename for a list of molecule names to plot, one per line. Molecules not present in this file will be ignored. Defaults to 'no'", \
                        default='no',type=str)
    parser.add_argument("-fpf", "--footprint_file", \
                        help="Plaintext DOCK6-format multi-footprint output file. Multiple footprints should be concatenated into a single file. Defaults to out_footprint_scored.txt", \
                        default='out_footprint_scored.txt')
    parser.add_argument("-ir", "--input_residues",\
                        help="Filename to input residues to plot on the x axis, with the rest being added to Remainder. Defaults to no; otherwise, use a filename. \
                        Ensure that the input file is just residue indices (ex 1), one per line",\
                        default="no")
    parser.add_argument("-r", "--ref_file", \
                        help="Filename for a reference footprint file. If 'no' provided, the first footprint in -fpf will be used.", \
                        default='no',type=str)
    parser.add_argument("-or", "--output_residues", \
                        help="Option to output important residues for the given molecules into a text file. Defaults to no; otherwise, use a filename.", \
                        default= "no")
    parser.add_argument("-mr","--max_res",\
                        help="Maximum number of residues to plot/write out. Defaults to 50.",\
                        default=50,type=int)
    parser.add_argument("-cut","--cutoff",\
                        help="The absolute cutoff value of energies for consideration of importance. Defaults to 0.2.",
                        default=0.2,type=float)
    parser.add_argument("-ro","--reference_override",\
                        help="[yes, no] Whether the reference peaks should take precedence when selecting residues to plot. Defaults to no.",
                        default="no",type=str)
    parser.add_argument("-fo", "--fileout", help="Specifies the name of the output graph for overlaid, or suffix for individual.\
                        Defaults to 'out'", default="out", type=str)
    parser.add_argument("-V", "--verbose", help="[yes, no] Enables verbose program output.", default="no", type=str)
    
    args = parser.parse_args()
    return args

###########################################
#           Initialize Variables          #
###########################################

args = parse_args()

header_keyword = "###  Molecule:"
header_prefix = "########## "

data_keyword = "resname"

# key: mol_name, val: [vdw, es]
all_perres_data = defaultdict(list) #each entry is tuple [vdw, es]
ref_perres_data = [] #each entry is tuple [vdw, es]
all_mol_data = defaultdict(dict)

# key: mol_name, val: [vdw, es]
remainders = defaultdict(list)

residue_names = []
residue_nums = []
resnums_to_plot = []
resnums_remainder = []

fp_infile = args.footprint_file
fp_reffile = args.ref_file
desired_num_to_plot = args.max_res
cutoff = args.cutoff
out_suff = args.fileout
prioritize_ref = args.reference_override

verbose = True if args.verbose == "yes" else False
plot_all_mols = False
mols_to_include = []
mol_names = []
linenum = 0

# plot parameters
individual = True if args.plot_type == "individual" else False
suptitLS=20
titLS=15
axisLS=12

lineweight=3

#######################################
#           Read in data              #
#######################################


# if user has requested specific mols, read those
if verbose: print("Checking if user requested specific mols")
if args.molnames != 'no':
    with open(args.molnames, 'r') as f:
        for line in f.readlines():
            if len(line.split()) == 1:
                mols_to_include.append(line.strip())
    if len(mols_to_include) == 0:
        print("ERROR: molnames file was specified but is empty! Exiting")
        quit()
else:
    plot_all_mols = True
        

# collect data from file
in_header = True
keep_mol = False
in_first_mol = True
first_mol_ref = True

if args.ref_file == "no": first_mol_ref = True
else: first_mol_ref = False

# reads through and looks for (in order):
#  - empty/all-comment lines (skips)
#  - the molecule name entry line (stores name)
#  - other header entries (currently throws away)
#  - the column label line (throws away)
# if none are met, the line is assumed to be data (stored)
if verbose: print("Reading in footprint data")
with open(fp_infile, 'r') as f:
    mol_data = defaultdict(float)
    perres_data = []
    for line in f.readlines():
        linelist = [x.strip() for x in line.split()]
        if line.strip() == "" or len(linelist) == 1: continue
        if header_keyword in line: 
            mol_data = defaultdict(float)
            in_header = True
            mol_name = linelist[2]
            if plot_all_mols or mol_name in mols_to_include:
                keep_mol = True
                if in_first_mol and len(mol_names) > 0:
                    in_first_mol = False
                mol_names.append(mol_name)
            else: 
                keep_mol = False
        elif in_header and header_prefix in line:
            dtype = linelist[1].strip(":")
            dentry = float(linelist[2])
            mol_data[dtype] = dentry
        elif data_keyword in line:
            in_header = False
        elif not in_header and keep_mol:
            perres_data.append(linelist)
            # if this is true, we are in the first molecule AND
            # the user has not provided a separate reference so we save it as the ref
            if in_first_mol:
                residue_names.append(linelist[0])
                residue_nums.append(int(linelist[1]))
            if first_mol_ref:
                ref_perres_data.append([float(linelist[2]), float(linelist[3])])
            all_perres_data[mol_name].append([float(linelist[5]), float(linelist[6])])
        linenum += 1


# sanity check # residues now that we know the number of residues in the footprints
if verbose: print("Sanity checking residue numbers")
if len(residue_names) < desired_num_to_plot:
    print(f"More residues were requested for plotting ({desired_num_to_plot}) than exist in the footprints ({len(residue_names)})!")
    print("Setting number to plot to the number of residues ({len(residue_names)}).")
    desired_num_to_plot = len(residue_names)


# if the user has provided a reference, we read that in and store it as well
if not first_mol_ref:
    if verbose: print("Reading in reference molecule file")
    linenum = 0
    in_header = True
    keep_mol = True
    in_first_mol = True
    
    with open(fp_reffile, 'r') as f:
        mol_data = defaultdict(float)
        perres_data = []
        for line in f.readlines():
            linelist = [x.strip() for x in line.split()]
            if line.strip() == "" or len(linelist) == 1: continue
            if header_keyword in line: 
                if in_first_mol and linenum > 2:
                    break
                mol_data = defaultdict(float)
                in_header = True
                mol_name = linelist[2]
            elif in_header and header_prefix in line:
                dtype = linelist[1].strip(":")
                dentry = float(linelist[2])
                mol_data[dtype] = dentry
            elif data_keyword in line:
                in_header = False
            elif not in_header:
                ref_perres_data.append([float(linelist[5]), float(linelist[6])])
            linenum += 1

# now that all data is accounted for, we need to figure out what residues
# should be plotted across all the footprints for consistency
if args.input_residues == "no":
    if verbose: print("Determining which residues to plot from cutoff compliance")

    #######################################
    #      Figure out plot priority       #
    #######################################
    
    # key resnum, value count
    count_over_cutoff = defaultdict(int)
    # key resnum, value float sum 
    abs_total_res_energies = defaultdict(float)
    total_res_energies = defaultdict(float)
    
    for i in range(0, len(residue_names)): count_over_cutoff[residue_nums[i]] = 0

    # count number of times each residue is over the designated cutoff
    # so we can determine which ones are most important to plot
    # (total energies stored as a backup in case not enough are over cutoff)
    # (if user only wants to use reference for cutoffs, ignore the rest)
    if not prioritize_ref:
        for mol in all_perres_data:
            for i in range(0, len(residue_names)):
                for val in all_perres_data[mol][i]:
                    if abs(val) > cutoff:
                        count_over_cutoff[residue_nums[i]] += 1
                    abs_total_res_energies[residue_nums[i]] += abs(val)
   
    # also count for reference 
    for i in range(0, len(residue_names)):
        for val in ref_perres_data[i]:
            if abs(val) > cutoff:
                count_over_cutoff[residue_nums[i]] += 1
            abs_total_res_energies[residue_nums[i]] += abs(val)
        
    ordered_count_over_cutoff = sorted(count_over_cutoff.items(), key=lambda kv: kv[1], reverse=True)
    ordered_total_res_energies = sorted(abs_total_res_energies.items(), key=lambda kv: kv[1], reverse=True)

    # select residues to plot based on the number of times they
    # surpassed the cutoffs across all molecules
    for pair in ordered_count_over_cutoff:
        if pair[1] > 0:
            resnums_to_plot.append(pair[0])
        if len(resnums_to_plot) == desired_num_to_plot: break
    
    # possible that this will not yield enough residues to plot, so
    # if the list isn't long enough then look at energy magnitude
    if len(resnums_to_plot) < desired_num_to_plot:
        for pair in ordered_total_res_energies:
            if pair[0] not in resnums_to_plot:
                resnums_to_plot.append(pair[0])
            if len(resnums_to_plot) == desired_num_to_plot: break

# otherwise, the user has kindly provided a residue file to use
else:
    if verbose: print("Reading in residues to plot from file")
    # read in residues from provided file
    with open(args.input_residues, 'r') as f:
        for line in f.readlines():
            ll = line.split()
            if len(ll) == 1:
                if int(line) in residue_nums:
                    resnums_to_plot.append(int(line.strip()))
                else:
                    print(f"Skipping user-requested residue {line} as it is not in the residues of the footprints provided")

    if len(resnums_to_plot) > len(residue_names):
        print("ERROR: user-provided residue list is LONGER than the number of residues present! Please revise.")
        quit()


resnums_to_plot = sorted(resnums_to_plot)

# populate "remainder" list now that we have a final list of residues to plot
for num in residue_nums:
    if num not in resnums_to_plot:
        resnums_remainder.append(num)


#######################################
#      Calculate REMAIN values        #
#######################################

# gather remainder values for each molecule
if verbose: print("Calculating REMAIN values")
for mol in mol_names:
    remainders[mol] = [0,0]
    for num in residue_nums:
        abs_idx = residue_nums.index(num)
        if num not in resnums_to_plot:
            remainders[mol][0] += all_perres_data[mol][abs_idx][0] # vdw
            remainders[mol][1] += all_perres_data[mol][abs_idx][1] # es

# gather remainder values for reference
remainders["REF"] = [0,0]
for num in residue_nums:
    abs_idx = residue_nums.index(num)
    if num not in resnums_to_plot:
        remainders["REF"][0] += ref_perres_data[abs_idx][0] # vdw
        remainders["REF"][1] += ref_perres_data[abs_idx][1] # es

#######################################
#         Plot the footprint          #
#######################################

ref_vdw_values = []
ref_es_values = []
xticks = []
xlabels = []

# prepare reference values and axis parameters for plotting
i = 0
for num in resnums_to_plot:
    abs_idx = residue_nums.index(num)
    xlabels.append(residue_names[abs_idx]+str(num))
    xticks.append(i)
    ref_vdw_values.append(ref_perres_data[abs_idx][0])
    ref_es_values.append(ref_perres_data[abs_idx][1])
    i += 1
xticks.append(i)
xlabels.append("REMAIN")
ref_vdw_values.append(remainders["REF"][0])
ref_es_values.append(remainders["REF"][1])

# create figure and subplot axes for VDW and ES plots
if verbose: print("Creating plot and formatting")
fig, axs = plt.subplots(2, 1, figsize=[10, 12])
fig.suptitle(f"Footprint Plots", size=suptitLS)
plt.rcParams["axes.prop_cycle"] = plt.cycler("color",plt.cm.tab10.colors)

# set up VDW plot area
axs[0].set_ylim(-10, 5)
axs[0].set_xlim(0, len(resnums_to_plot)+1)
axs[0].xaxis.set_major_locator(MultipleLocator(1))
axs[0].xaxis.set_major_formatter(FormatStrFormatter('%s'))
axs[0].set_xticks(xticks)
axs[0].set_xticklabels(xlabels, rotation=90)
axs[0].set_ylabel("VDW Energy (kcal/mol)")
axs[0].xaxis.grid(which='major', color='black', linestyle='solid')

# set up ES plot area
axs[1].set_ylim(-10, 5)
axs[1].set_xlim(0, len(resnums_to_plot)+1)
axs[1].xaxis.set_major_locator(MultipleLocator(1))
axs[1].xaxis.set_major_formatter(FormatStrFormatter('%s'))
axs[1].set_xticks(xticks)
axs[1].set_xticklabels(xlabels, rotation=90)
axs[1].set_ylabel("ES Energy (kcal/mol)")
axs[1].xaxis.grid(which='major', color='black', linestyle='solid')


first_loop = True
if verbose: print("Plotting footprints")
for mol in mol_names:
    # plot the reference first if the plot is blank
    if individual or first_loop:
        for line in axs[0].lines:
            line.remove()
        for line in axs[1].lines:
            line.remove()
        axs[0].plot(xticks, ref_vdw_values, linewidth=lineweight, color="black")
        axs[1].plot(xticks, ref_es_values, linewidth=lineweight, color="black")    
        first_loop = False
    
    vdw_values = []
    es_values = []

    # gather data for plotting for VDW and ES from per-residue data
    i = 0
    for num in resnums_to_plot:
        abs_idx = residue_nums.index(num)
        vdw_values.append(all_perres_data[mol][abs_idx][0])
        es_values.append(all_perres_data[mol][abs_idx][1])
        i += 1
    vdw_values.append(remainders[mol][0])
    es_values.append(remainders[mol][1])

    # if we are plotting each individually, then generate the legend now,
    # use a set color, and save the figure as a PNG
    if individual:
        fig.suptitle(f"{mol} Footprint Plots", size=suptitLS)
        axs[0].plot(xticks, vdw_values, linewidth=lineweight, color="red")
        axs[1].plot(xticks, es_values, linewidth=lineweight, color="red")
        axs[0].legend(['Reference'] + [mol])
        axs[1].legend(['Reference'] + [mol])
        fig.tight_layout()
        plt.savefig(f"{mol}.{out_suff}.png",dpi=300)
    else: # otherwise, let the color cycler auto-generate
        axs[0].plot(xticks, vdw_values, linewidth=lineweight)
        axs[1].plot(xticks, es_values, linewidth=lineweight)
            

if not individual:
    axs[0].legend(['Reference'] + mol_names)
    axs[1].legend(['Reference'] + mol_names)
    fig.tight_layout()
    plt.savefig(f"{out_suff}.png",dpi=300)

if args.output_residues != "no":
    if verbose: print("Writing out residues plotted")
    with open(args.output_residues, 'r') as f:
        output = []
        for num in resnums_to_plot:
            abs_idx = residue_nums.index(num)
            output.append(residue_names[abs_idx] + str(residue_nums[abs_idx]))
        f.write("\n".join(output))
            
