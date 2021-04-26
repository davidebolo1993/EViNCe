#!/usr/bin/python3 env

import os
import sys
from datetime import date

import pyfaidx
import pybedtools
import pysam


def main():

	'''
	Convert simple BED with SVs from HACk to VCF compatible with truvari
	'''

	FASTA=os.path.abspath(sys.argv[1])
	BED=os.path.abspath(sys.argv[2])
	OUT=os.path.abspath(sys.argv[3])
	FA=pyfaidx.Fasta(FASTA)
	VCF=''
	VCF+='##fileformat=VCFv4.2\n'
	VCF+='##source=VISOR\n'
	VCF+='##filedate=' + ''.join(str(date.today()).split('-')) + '\n'
	VCF+='\n'.join(["##contig=<ID="+str(x)+",length="+str(len(FA[x]))+">" for x in FA.keys()])
	VCF+='\n##ALT=<ID=INS,Description="Insertion of novel sequence relative to the reference">\n'
	VCF+='##ALT=<ID=DEL,Description="Deletion relative to the reference">\n'
	VCF+='##ALT=<ID=DUP,Description="Region of elevated copy number relative to the reference">\n'
	VCF+='##ALT=<ID=INV,Description="Inversion of reference sequence">\n'
	VCF+='##ALT=<ID=BND,Description="Breakend of translocation">\n'
	VCF+='##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
	VCF+='##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n'
	VCF+='##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n'
	VCF+='##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
	VCF+='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSI00001\n'

	bed=pybedtools.BedTool(BED)
	bedsrtd=bed.sort()
	variants=[]


	for i,x in enumerate(bedsrtd):

		CHROM=x.chrom if x[3]!= 'translocation cut-paste' else x[4].split(':')[1]
		ID='var' + str(i+1)
		REF='N'
		FILTER='PASS'
		QUAL='.'
		GT='0/1'

		if x[3] == 'deletion':

			POS=str(x.start)
			ALT='<DEL>'
			INFO='SVTYPE=DEL;END=' + str(x.end) + ';SVLEN=-' + str(x.end-x.start)

		elif x[3] == 'insertion':

			POS=str(x.end)
			ALT='<INS>'
			INFO='SVTYPE=INS;END=' + str(x.end) + ';SVLEN=' + str(len(x[4]))

		elif x[3] == 'tandem duplication':

			POS=str(x.start)
			ALT='<DUP>'
			INFO='SVTYPE=DUP;END=' + str(x.end) + ';SVLEN=' + str(x.end-x.start)

		elif x[3] == 'inversion':

			POS=str(x.start)
			ALT='<INV>'
			INFO='SVTYPE=INV;END=' + str(x.end) + ';SVLEN=' + str(x.end-x.start)

		else:

			POS=x[4].split(':')[2]
			p=x.chrom+':'+str(x.start) if x[4].split(':')[3] == 'forward' else x.chrom+':'+str(x.end)
			ALT='N['+ p + '[' if x[4].split(':')[3] == 'forward' else 'N]'+ p + ']'
			INFO='SVTYPE=BND'

		variants.append('\t'.join([CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, 'GT', GT]))


	VCF+='\n'.join(variants)

	with open(OUT,'w') as vcfout:

		vcfout.write(VCF+'\n')


if __name__ == '__main__':

	main()

