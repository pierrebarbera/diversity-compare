#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    Read all fasta files and build a sorted OTU contingency
    table. Usage: python OTU_contingency_table.py [input files]
"""

from __future__ import print_function

__author__ = "Frédéric Mahé <frederic.mahe@cirad.fr>"
__date__ = "2016/03/07"
__version__ = "$Revision: 5.0"

import os
import re
import sys
import operator
import argparse

#*****************************************************************************#
#                                                                             #
#                                  CLI                                        #
#                                                                             #
#*****************************************************************************#

def fail( msg="Unspecified Failure" ):
    print( "ERROR: {}".format( msg ) )
    sys.exit( 1 )

def command_line_args_parser():
    """
    Return an instance of argparse that can be used to process command line arguemnts.
    """

    def unit_interval(x):
        x = float(x)
        if x < 0.0 or x > 1.0:
            raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
        return x

    def existing_file( path ):
        if not os.path.isfile( path ) or not os.access( path, os.R_OK ):
            fail( "File doesn't exist or not readable: {}".format( path ) )
        return path


    # Init an args parser, with a group of required named arguments. It is just nicer to use named
    # arguments than having to rely on their order (i.e., use positional arguments instead).
    parser = argparse.ArgumentParser(
        description="Builds a contingency OTU table."
    )
    parser_required = parser.add_argument_group('required arguments')

    parser_required.add_argument(
        'per_sample_fasta_files',
        help="Per sample fasta files.",
        action='store',
        type=str,
        nargs='+'
    )

    parser_required.add_argument(
        "--stats",
        help="Clustering statistics file.",
        action="store",
        type=existing_file,
        dest="stats"
    )

    parser_required.add_argument(
        "--swarms",
        help="Swarms file.",
        action="store",
        type=existing_file,
        dest="swarms"
    )

    parser.add_argument(
        "--reps",
        help="Cluster representatives (fasta). If specified, enables 'length' and 'sequence' columns.",
        action="store",
        type=existing_file,
        dest="reps"
    )

    parser.add_argument(
        "--uchime",
        help="UCHIME file. If specified, enables 'chimera' column.",
        action="store",
        type=existing_file,
        dest="uchime"
    )

    parser.add_argument(
        "--quality",
        help="Quality file. If specified, enables 'quality' column.",
        action="store",
        type=existing_file,
        dest="quality"
    )

    parser.add_argument(
        "--assignments",
        help="Assignments file. If specified, enables 'identity', 'taxonomy' and 'references' columns",
        action="store",
        type=existing_file,
        dest="assignments"
    )

    parser.add_argument(
        '--threads',
        help="Number of threads to use.",
        action='store',
        dest='threads',
        default=1,
        type=int
    )

    parser.add_argument(
        "--verbose",
        help="Increase output verbosity.",
        action="store_true"
    )

    return parser

def command_line_args_postprocessor( args ):
    # Make sure that all paths are fully resolved and dirs have no trailing slashes.
    # args.jplace_file = os.path.abspath( os.path.realpath( args.jplace_file ))

    return args

def command_line_args():
    """
    Return a parsed and processed list of the command line arguments that were provided when
    running this script.
    """

    # Parse the given arguments from the command line, post-process them, return the result.
    parser = command_line_args_parser()
    args = parser.parse_args()
    args = command_line_args_postprocessor( args )
    return args

#*****************************************************************************#
#                                                                             #
#                                  Functions                                  #
#                                                                             #
#*****************************************************************************#


def representatives_parse( args ):
    """
    Get seed sequences.
    """
    if args.reps:
        separator = ";size="
        representatives_file = args.reps
        representatives = dict()
        with open(representatives_file, "rU") as representatives_file:
            for line in representatives_file:
                if line.startswith(">"):
                    amplicon = line.strip(">;\n").split(separator)[0]
                else:
                    representatives[amplicon] = line.strip()

        return representatives
    else:
        return None


def stats_parse( args ):
    """
    Map OTU seeds and stats.
    """
    separator = "\t"
    stats_file = args.stats
    stats = dict()
    seeds = dict()
    with open(stats_file, "rU") as stats_file:
        for line in stats_file:
            cloud, mass, seed, seed_abundance = line.strip().split(separator)[0:4]
            stats[seed] = int(mass)
            seeds[seed] = (int(seed_abundance), int(cloud))
    # Sort OTUs by decreasing mass
    sorted_stats = sorted(stats.iteritems(),
                          key=operator.itemgetter(1, 0))
    sorted_stats.reverse()

    return stats, sorted_stats, seeds


def swarms_parse( args ):
    """
    Map OTUs.
    """
    separator = "_[0-9]+|;size=[0-9]+;?| "  # parsing of abundance annotations
    swarms_file = args.swarms
    swarms = dict()
    with open(swarms_file, "rU") as swarms_file:
        for line in swarms_file:
            line = line.strip()
            amplicons = re.split(separator, line)[0::2]
            seed = amplicons[0]
            swarms[seed] = [amplicons]

    return swarms


def uchime_parse( args ):
    """
    Map OTU's chimera status.
    """
    if args.uchime:
        separator = " "
        uchime_file = args.uchime
        uchime = dict()
        with open(uchime_file, "rU") as uchime_file:
            for line in uchime_file:
                OTU = line.strip().split("\t")
                try:
                    seed = OTU[1].split(";")[0]
                except IndexError:  # deal with partial line (missing seed)
                    continue
                try:
                    status = OTU[17]
                except IndexError:  # deal with unfinished chimera detection runs
                    status = "NA"
                uchime[seed] = status

        return uchime
    else:
        return None


def quality_parse( args ):
    """
    List good amplicons.
    """
    if args.quality:
        quality_file = args.quality
        quality = dict()
        with open(quality_file, "rU") as quality_file:
            for line in quality_file:
                sha1, qual, length = line.strip().split()
                quality[sha1] = float(qual) / int(length)

        return quality
    else:
        return None


def stampa_parse( args ):
    """
    Map amplicon ids and taxonomic assignments.
    """
    if args.assignments:
        separator = "\t"
        stampa_file = args.assignments
        stampa = dict()
        with open(stampa_file, "rU") as stampa_file:
            for line in stampa_file:
                amplicon, abundance, identity, taxonomy, references = line.strip().split(separator)
                stampa[amplicon] = (identity, taxonomy, references)

        return stampa
    else:
        return None


def fasta_parse( args ):
    """
    Map amplicon ids, abundances and samples.
    """
    separator = ";size="
    fasta_files = args.per_sample_fasta_files
    samples = dict()
    amplicons2samples = dict()
    for fasta_file in fasta_files:
        sample = os.path.basename(fasta_file)
        sample = sample.split(".")[0]
        samples[sample] = samples.get(sample, 0) + 1
        with open(fasta_file, "rU") as fasta_file:
            for line in fasta_file:
                if line.startswith(">"):
                    amplicon, abundance = line.strip(">;\n").split(separator)
                    abundance = int(abundance)
                    if amplicon not in amplicons2samples:
                        amplicons2samples[amplicon] = {sample: abundance}
                    else:
                        # deal with duplicated samples
                        amplicons2samples[amplicon][sample] = amplicons2samples[amplicon].get(sample, 0) + abundance
    # deal with duplicated samples
    duplicates = [sample for sample in samples if samples[sample] > 1]
    if duplicates:
        print("Warning: some samples are duplicated", file=sys.stderr)
        print("\n".join(duplicates), file=sys.stderr)
    samples = sorted(samples.keys())

    return amplicons2samples, samples


def print_table(reps, stats, sorted_stats,
                swarms, uchime, amplicons2samples,
                samples, quality, seeds, stampa):
    """
    Export results.
    """
    header = ["OTU", "total", "cloud",
          "amplicon", "length", "abundance",
          "chimera", "spread", "quality",
          "sequence", "identity", "taxonomy", "references"]

    if not reps:
        header.remove( "length" )
        header.remove( "sequence" )

    if not uchime:
        header.remove( "chimera" )

    if not quality:
        header.remove( "quality" )

    if not stampa:
        header.remove( "identity" )
        header.remove( "taxonomy" )
        header.remove( "references" )

    # Print table header
    print("\t".join(header),
          "\t".join(samples),
          sep="\t", file=sys.stdout)

    # Print table content
    i = 1
    for seed, abundance in sorted_stats:
        line=[ i, abundance ]

        occurrences = dict([(sample, 0) for sample in samples])
        for amplicons in swarms[seed]:
            for amplicon in amplicons:
                for sample in samples:
                    occurrences[sample] += amplicons2samples[amplicon].get(sample, 0)
        spread = len([occurrences[sample] for sample in samples if occurrences[sample] > 0])
        sequence_abundance, cloud = seeds[seed]

        line.extend([ cloud, seed ])
        if reps:
            sequence = reps[seed]
            line.extend([ len(sequence) ])
        line.extend([ sequence_abundance ])

        # Chimera checking (deal with incomplete cases. Is it useful?)
        if uchime:
            if seed in uchime:
                chimera_status = uchime[seed]
            else:
                chimera_status = "NA"
            line.extend([ chimera_status ])

        line.extend([ spread ])

        # Quality
        if quality:
            if seed in quality:
                high_quality = quality[seed]
            else:
                high_quality = "NA"
            line.extend([ high_quality ])

        if reps:
            sequence = reps[seed]
            line.extend([ sequence ])

        # Chimera checking (deal with incomplete cases. Is it useful?)
        if stampa:
            if seed in stampa:
                identity, taxonomy, references = stampa[seed]
            else:
                identity, taxonomy, references = "NA", "NA", "NA"
            line.extend([ identity, taxonomy, references ])

        # output
        print("\t".join(map(str, line)),
            "\t".join([str(occurrences[sample]) for sample in samples]),
            sep="\t", file=sys.stdout)

        i += 1

    return


#*****************************************************************************#
#                                                                             #
#                                     Body                                    #
#                                                                             #
#*****************************************************************************#

if __name__ == '__main__':
    args = command_line_args()

    representatives = representatives_parse( args )

    # Parse stats
    stats, sorted_stats, seeds = stats_parse( args )

    # Parse swarms
    swarms = swarms_parse( args )

    # Parse uchime
    uchime = uchime_parse( args )

    # Parse quality
    quality = quality_parse( args )

    # Parse taxonomic assignment results
    stampa = stampa_parse( args )

    # Parse fasta files
    amplicons2samples, samples = fasta_parse( args )

    # Print table header
    print_table(representatives, stats, sorted_stats, swarms,
                uchime, amplicons2samples, samples, quality,
                seeds, stampa)
