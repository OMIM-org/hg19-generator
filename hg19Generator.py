#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#
# This is a simple script to generate GRCh37/hg19 genomic coordinate for MIM genes.
#
# You will need the GCF_000001405.25_GRCh37.p13_genomic.gff.gz file which can downloaded from NCBI:
#
#   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz
#
# And mim2gene.txt which can be downloaded from OMIM:
#
#   https://omim.org/downloads
#


# Imports
import re


# Gene dict, K - 'geneID', V - genomic coordinate
geneDict = dict()

# Process the GRCh37 data
with open('./GCF_000001405.25_GRCh37.p13_genomic.gff') as fileHandle:
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

        # Add the genomic coordinate to the gene dict
        if entrezGeneID not in geneDict:
            geneDict[entrezGeneID] = set()
        geneDict[entrezGeneID].add(genomicCoordinate)



# Process the OMIM mim2gene file
with open('./mim2gene.txt') as fileHandle:
    for line in fileHandle:

        # Skip comments
        if line.startswith('#'):
            continue

        # Strip trailing new line
        line = line.strip('\n')

        # Get the values
        valueList = line.split('\t')

        # Get the fields
        mimNumber = valueList[0]
        entryType = valueList[1]
        entrezGeneID = ''
        if len(valueList) > 2:
            entrezGeneID = valueList[2]
        approvedGeneSymbol = ''
        if len(valueList) > 2:
            approvedGeneSymbol = valueList[3]
        ensemblGeneID = ''
        if len(valueList) > 2:
            ensemblGeneID = valueList[4]

        # Write out the data
        if entrezGeneID in geneDict:
            print('{0}\t{1}\t{2}\t{3}'.format(entrezGeneID, approvedGeneSymbol, mimNumber, '|'.join(geneDict[entrezGeneID])))
