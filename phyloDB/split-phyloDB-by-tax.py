#! /usr/bin/env python
import os
import sys
import argparse
import screed

import pandas as pd
import sourmash
from sourmash.logging import notify
from sourmash.lca.lca_utils import make_lineage

def make_outdir(output_dirname):
    if not os.path.exists(output_dirname):
        try:
            os.makedirs(output_dirname)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

def split_lineages(row):
    lins = row['taxonomy']
    # do it this way in case not all ranks included for each taxonomy entry
    smashlin = make_lineage(lins)
    for lin_tup in smashlin[::-1]:
        rank = lin_tup.rank
        lin = lin_tup.name
        row[rank] = lin
    return row


def write_taxfasta(this_ident, tax_records, alltax, info_out, n, n_ff, output_dir):
    # first, check if we've seen this name before - should not happen
    if this_ident in alltax:
        notify(f'Error: found non-adjacent duplicated name, {this_ident}')
        sys.exit(-1)
    # add to alltax for bookkeeping
    alltax.add(this_ident)

    # make filename and write fasta
    short_id = this_ident.split(' ')[0]
    this_filename = f"{output_dir}/{short_id}.fa"

    if n_ff % 5000 == 0:
        print(f"writing file for {n_ff}'th strain, '{this_ident}' after processing {str(n+1)}th fasta entry\n")

    with open(this_filename, 'w') as out_fa:
        for rc in tax_records:
            out_fa.write(f">{rc.name}\n{rc.sequence}\n")

    # write fromfile entry
    info_out.write(f"{this_ident},,{this_filename}\n")

def main(args):
    # strain_name     peptide_count   taxonomy
    taxinfo = pd.read_csv(args.taxonomy_tsv, sep='\t')
    # make ourselves a better ident (rm difficult characters!)
    taxinfo['filefriendlyname'] = taxinfo['strain_name'].str.replace("\"", "").str.replace("\'", "").str.replace(" ", "_").str.replace('/', "__").str.replace(',', "-")
    #taxinfo['ident'] = 'phylodb_' + taxinfo.index.astype(str) + ' ' + taxinfo['strain_name'] # would prefer this, but **commas**, quotes, slashes break things
    taxinfo['ident'] = 'phylodb_' + taxinfo.index.astype(str) + ' ' + taxinfo['filefriendlyname']
    # make lineages spreadsheet
    notify(f"Making lineages csv: '{args.output_lineages}' \n")
    taxinfo = taxinfo.apply(split_lineages, axis=1)
    lineage_cols = ['ident'] + list(sourmash.lca.lca_utils.taxlist(include_strain=False))
    lineages = taxinfo[lineage_cols]
    lineages.to_csv(args.output_lineages, index=False)
    # now use taxinfo below
    taxinfo.set_index('filefriendlyname', inplace=True)

    notify(f"Splitting {args.fasta} by strain. Writing files to '{args.output_dir}' \n")
    output_dir = args.output_dir
    make_outdir(args.output_dir)

    # split fasta by strain:
    tax_records = []
    this_ident = None
    alltax = set()
    n_ff = 0
    with open(args.output_csv, 'w') as info_out:
        # make sourmash sketch fromfile csv
        info_out.write("name,genome_filename,protein_filename\n") # header needs these three, even if we only have proteins
        # loop through fasta file
        for n, record in enumerate(screed.open(args.fasta)):
            tax_name = record.name.rsplit("\t")[-1]
            tn = tax_name.replace("\"", "").replace("\'", "").replace(" ", "_").replace('/', "__").replace(',', "-")

            ident_name = taxinfo.at[tn, 'ident']
            # remove weird spaces and quotes

            #first entry: set this_ident
            if this_ident == None:
                this_ident = ident_name

            # same org: append record
            if ident_name == this_ident:
                tax_records.append(record)
            else:
                # we're onto a diff organism - write!
                write_taxfasta(this_ident, tax_records, alltax, info_out, n, n_ff, args.output_dir)

                # reset tax_records and this_ident
                tax_records = [record]
                this_ident = ident_name
                n_ff+=1

        # catch final org
        write_taxfasta(this_ident, tax_records, alltax, info_out, n, n_ff, args.output_dir)

    print(f"{str(n+1)} total entries written to {n_ff} individual fasta files\n")


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", default= "phylodb_1.076.pep.fa.gz")
    p.add_argument("--output-dir", default= "phylodb_1.076_taxsplit.ff")
    p.add_argument("--output-csv", default= "phylodb_1.076.fromfile.ff.csv")
    p.add_argument("--taxonomy-tsv", default="phylodb_1.076.taxonomy.txt")
    p.add_argument("--output-lineages", default="phylodb_1.076.lineages.ff.csv")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
