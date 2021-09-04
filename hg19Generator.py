#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#
# This is a simple script to generate GRCh37/hg19 genomic coordinate for MIM genes.
#
# You will need ref_GRCh37.p5_top_level.gff3.gz which can downloaded from NCBI:
#
#   https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/H_sapiens/ARCHIVE/BUILD.37.3/GFF/ref_GRCh37.p5_top_level.gff3.gz
#
# And genemap2.txt which canb be downloaded from OMIM:
#
#   https://omim.org/downloads
#
# (registration required)
#


# Imports
import re


# Gene dict, K - 'geneID:chromosone', V - genomic coordinate
geneDict = dict()

# Process the GRCh37 data
with open('./ref_GRCh37.p5_top_level.gff3') as fileHandle:
    for line in fileHandle:

        # Skip comments
        if line.startswith('#'):
            continue

        # Strip trailing new line
        line = line.strip('\n')

        # Get the values
        valueList = line.split('\t')

        # Get the fields
        accessionNumber = valueList[0]
        sequenceType = valueList[2]
        genomicPositionStart = valueList[3]
        genomicPositionEnd = valueList[4]
        identifiers = valueList[8]

        # Skip non-genes
        if sequenceType != 'gene':
            continue

        # Skip non-genes
        if not accessionNumber.startswith('NC_'):
            continue

        # Extract the Entrez Gene ID
        matcher = re.search(r'GeneID:(\d+)', identifiers)
        if not matcher:
            continue
        entrezGeneID = matcher.group(1)

        # Extract the chromosome
        matcher = re.match(r'^NC_(\d{6})\.\d+', accessionNumber)
        if not matcher:
            continue
        chromosome = int(matcher.group(1))
        if chromosome == 23:
            chromosome = 'X'
        elif chromosome == 24:
            chromosome = 'Y'

        # Create the genomic coordinate
        genomicCoordinate = 'chr{0}:{1}-{2}'.format(chromosome, valueList[3], valueList[4])

        # Create the gene dict key and add the genomic coordinate to the gene dict
        geneDictKey = '{0}:{1}'.format(entrezGeneID, chromosome)
        geneDict[geneDictKey] = genomicCoordinate



# Process the OMIM genenmap2
with open('./genemap2.txt') as fileHandle:
    for line in fileHandle:

        # Skip comments
        if line.startswith('#'):
            continue

        # Strip trailing new line
        line = line.strip('\n')

        # Get the values
        valueList = line.split('\t')

        # Get the fields
        chromosome = valueList[0]
        mimNumber = valueList[5]
        approvedGeneSymbol = valueList[8]
        entrezGeneID = valueList[9]

        # Skip entries that don't have an Entrez Gene ID
        if not entrezGeneID:
            continue

        # Create the gene dict key and write out the data
        geneDictKey = '{0}:{1}'.format(entrezGeneID, chromosome[3:])
        if geneDictKey in geneDict:
            print('{0}\t{1}\t{2}\t{3}'.format(entrezGeneID, approvedGeneSymbol, mimNumber, geneDict.get(geneDictKey, '')))



