import csv
import pandas as pd

configfile: "mmetsp.yml"
basename = config['basename']
out_dir = config['output_dir']
logs_dir = os.path.join(out_dir, 'logs')

# check params are in the right format, build alpha-ksize combos
alpha_ksize_scaled=[]
dna_param_strs=""
prot_param_strs=""
for alpha, info in config["alphabet_info"].items():
    scaled = info["scaled"]
    ksize = info["ksize"]
    if not isinstance(scaled, list):
        scaled = [scaled]
        config["alphabet_info"][alpha]["scaled"] = scaled
    min_scaled = min(scaled)
    if not isinstance(ksize, list):
        ksize=[ksize]
        config["alphabet_info"][alpha]["ksize"] = ksize
    # build a parameter for the right combinations
    #config["alphabet_info"][alpha]["select_params"] = expand("{alpha}-k{ksize}-scaled{scaled}", alpha=alpha, scaled=scaled, ksize=ksize)
    #these_alpha_ksize = expand("{alpha}-k{{ksize}}", ksize = ksize)
    
    # make param strings
    ks = [ f'k={k}' for k in ksize]
    ks = ",".join(ks)
    #scaled = min(scaled) #take minimum value of scaled list (just downsample for larger scaled vals)
    if alpha == 'nucleotide':
        alpha = "dna"
        this_pstr = f" -p {ks},scaled={min_scaled},abund,{alpha} "
        dna_param_strs += this_pstr
    else:
        this_pstr = f" -p {ks},scaled={min_scaled},abund,{alpha} "
        prot_param_strs += this_pstr
    alpha_ksize_scaled += expand(f"{alpha}-k{{ksize}}-sc{min_scaled}", ksize = ksize)#, scaled=scaled)


# some easier vars
alphabet_info = config['alphabet_info']
prot_ksizes = alphabet_info['protein']['ksize']
nucl_ksizes = alphabet_info['nucleotide']['ksize']

print(alpha_ksize_scaled)

rule all:
    input:
        expand(f"{out_dir}/databases/{basename}.{{aks}}.zip", aks =alpha_ksize_scaled),

all_param_str = {"nucleotide": dna_param_strs, "protein": prot_param_strs}
rule sketch_fromfile:
    input: 
        fromfile_csv = config['fromfile_csv'] 
    output: 
        f"{out_dir}/databases/{basename}.{{alphabet}}-k{{ksize}}-sc{{scaled}}.zip"
    log: 
        f"{logs_dir}/databases/{basename}.{{alphabet}}-k{{ksize}}-sc{{scaled}}.log"
    benchmark: 
        f"{logs_dir}/sketch-fromfile/{basename}.{{alphabet}}-k{{ksize}}-sc{{scaled}}.benchmark"
    threads: 1
    resources:
        mem_mb=3000,
        time=20000,
    conda: "conf/env/sourmash-4.4.yml" 
    shell:
        """
        sourmash sketch fromfile {input.fromfile_csv} -p k={wildcards.ksize},{wildcards.alphabet},scaled={wildcards.scaled},abund -o {output} --force 2> {log}
        """
        #sourmash sketch fromfile output.rank-compare/gtdb-rs207.fromfile.csv -o gtdb-rs207.no-prodigal.zip -p k=10,scaled=200,protein

