#!/usr/bin/env python
# Copyright (C) 2019 Umesh Ghoshdastider, All rights reserved, version 1.0
# Installation Instructions

# Ubuntu, Mint and other Debian based distributions
# sudo pip install weblogo biopython
# sudo apt-get install blast cd-hit fasttree hmmer mafft muscle python-pip raxml 
# download usearch from www.drive5.com

# Mac OS

# sudo brew install blast cd-hit fasttree hmmer mafft muscle python python-pip raxml wxpython
# sudo pip install weblogo biopython matplotlib

# Windows: Download and install each programs separately

#from argparse import ArgumentParser
from Bio.SeqIO import convert,parse
import glob
from gooey import Gooey, GooeyParser
from matplotlib.pylab import *
import os
import shutil
from subprocess import call,Popen

# shlex.split('command') says how to parse shell arguments

# help='All input files need extension', file path must not include '.' in directory name, no gap expected in fasta header, output will be in input file folder
# variables 
# a = input_args, p = parser, sp = subparser, a1 = sp1_add_argument

# include full PATH name like /usr/bin/blastp if it is not in your path
# otherwise use only executable name like blastp

blastp_path='blastp' #'/usr/bin/blastp'
blastn_path='blastn' #'/usr/bin/blastn'
makeblastdb_path='makeblastdb' #'/usr/bin/makeblastdb'

cdhit_path='cdhit' #'/home/umesh/soft/bin/cdhit'
chitest_path='cdhit-est' #'/home/umesh/soft/bin/cdhit-est'

fasttree_path='fasttree' #'/home/umesh/soft/bin/fasttree'
hmmbuild_path='hmmbuild' #'/home/umesh/soft/bin/'
hmmsearch_path='hmmsearch' #'/home/umesh/soft/bin/'

jmodeltest_path='/home/umesh/soft/phylo/jmodeltest2/dist/jModelTest.jar'
mafft_path='/home/umesh/soft/msa/mafft-7.266-with-extensions/build/bin/mafft'
muscle_path='muscle' #'/home/umesh/soft/bin/muscle'
raxml_path='raxml' #'/home/umesh/soft/bin/raxml'
usearch_path='usearch' #'/home/umesh/soft/bin/usearch'
weblogo_path='weblogo' #'/usr/local/bin/weblogo'


format='clustal embl fasta genbank nexus phylip stockholm tab'.split()
	
def plot_table(a):
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_plot_table.png'
	f=loadtxt(a.ifile)
	if a.xcol1:
		for i in range(1,f.shape[1]):
			plot(f[:,0], f[:,i])
	else:
		plot(f)
	savefig(out,dpi=150)
	show()
	
def dna_complement(a):
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_compement.fa'

	with open(out,'w') as f:
		for rec in parse(a.ifile,'fasta'):
			f.write('>'+rec.id+'\n')
			f.write(str(rec.seq.complement())+'\n')
			
def dna_reverse_complement(a):
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_rev_complement.fa'

	with open(out,'w') as f:
		for rec in parse(a.ifile,'fasta'):
			f.write('>'+rec.id+'\n')
			f.write(str(rec.seq.reverse_complement())+'\n')
	
def dna_translate(a):
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_translate.fa'

	with open(out,'w') as f:
		for rec in parse(a.ifile,'fasta'):
			f.write('>'+rec.id+'\n')
			f.write(str(rec.seq.translate())+'\n')
	

def seq_format_conversion(a):
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'.'+a.oformat

	convert(a.ifile, a.iformat, out, a.oformat)


def remove_alignment_gap(a):
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_nogap.fa'
	with open(out,'w') as f:
		for rec in parse(a.ifile,'fasta'):
			f.write('>'+rec.id+'\n')
			f.write(str(rec.seq).replace('-','')+'\n')
		
		
def extract_geneid_fasta(a):
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_'+a.geneid+'.fa'
		
	with open(out,'w') as f:
		for rec in parse(a.ifile,'fasta'):
			if a.geneid in rec.id:
				f.write('>'+rec.id+'\n')
				f.write(str(rec.seq)+'\n')


def extract_genelist_fasta(a):
	if a.ofile:
		out=a.ofile
	else:
		out=os.path.basename(a.ifile).split('.')[0]+'_'+os.path.basename(a.genelist).split('.')[0]+'.fa'
	
	l=open(a.genelist).read().split('\n')
	
	with open(out,'w') as f:
		for rec in parse(a.ifile,'fasta'):
			for i in l:
				if i in rec.id:
					f.write('>'+rec.id+'\n')
					f.write(str(rec.seq)+'\n')
			
				
def alignment_muscle(a):
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_muscle.fa'
	
	cmd=[muscle_path,'-out', out,'-in', a.ifile]
	if a.maxiters:
		cmd.append('-maxiters')
		cmd.append(a.maxiters)
		call(cmd)
	else:
		call(cmd)


def alignment_mafft(a):
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_mafft.fa'
	
	if a.maxiterate:
		iter=a.maxiterate
	else:
		iter='0'
	
	if a.thread:
		th=a.thread
	else:
		th='1'
		
	with open(out,'w') as f:
		Popen([mafft_path, '--thread', th, '--maxiterate', iter, a.ifile], stdout=f)


def tree_fasttree(a):
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_ftree'
		
	if a.protein:
		call([fasttree_path,'-mlacc', '2', '-slownni', '-spr', '4', '-wag', '-gamma', '-out', out, a.ifile])
	else:
		call([fasttree_path,'-mlacc', '2', '-slownni', '-spr', '4', '-nt', '-gtr', '-gamma', '-out', out, a.ifile])		
		
		
def tree_raxml(a):
# -N 2 : two independent runs
	f=os.path.basename(a.ifile).split('.')[0]
	
	if a.protein:
		subs_matrix='PROTGAMMAILG'
	
	else:
		subs_matrix='GTRGAMMAI'
	
	if a.threads:
		th=a.threads
	else:
		th='1'
	

	if a.bootstrap:
		call([raxml_path, '-m', subs_matrix, '-p', '12345', '-x', '67890', '-N', a.bootstrap,'-f','a', '-s', a.ifile, '-n', f, '-T', th ])
	else:
		call([raxml_path, '-m', subs_matrix, '-p', '12345', '-N', '2', '-s', a.ifile, '-n', f, '-T', th ])
				
	dir='RAxML_'+f
	os.mkdir(dir)
	for fi in glob.glob('RAxML*'+f+'*'):
		shutil.move(fi, dir)


def tree_prep_ali(a): 
	'''Remove illegal characters ( ) [ ] ; : from alignment header for tree building'''
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_4tree.fa'
	with open(out,'w') as f:
		for rec in parse(a.ifile,'fasta'):
			for i in '( ) [ ] ; :'.split():
				rec.id=rec.id.replace(i,'')
			f.write('>'+rec.id  +'\n')
			f.write(str(rec.seq)+'\n')
			
			
def fasta_statistics(a):
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_stat.png'
	l=[]
	for rec in parse(a.ifile,'fasta'):
		l.append(len(rec))
	l=array(l)
	hist(l)
	title('Length- min: '+str(l.min())+', max: '+str(l.max())+', avg: '+str(l.mean())[:-7])
	xlabel('Sequence Length')
	ylabel('Frequency')
	
	savefig(out,dpi=150)
	show()
	
def remove_seq_from_ali(a):
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_50.fa'

	f1=open(a.ifile.split('.')[0]+'_short_seq','w')
	
	with open(out,'w') as f2:
		for rec in parse(a.ifile,'fasta'):
			gap=str(rec.seq).count('-')*100./len(rec.seq)
			if a.pctgap:
				pct=a.pctgap
			else:
				pct=50
			if gap<pct:
				f2.write('>'+rec.id  +'\n')
				f2.write(str(rec.seq)+'\n')
			else:
				f1.write(rec.id + '\t' + str(round(pct,2))+'\n')


def remove_empty_lines(a):
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_no_empty_line.fa'
	with open(out,'w') as f:
		for line in open(a.ifile):
			if line.strip()!='':
				f.write(line)


def rename_fasta2numbers(a):
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_rename2numbers.fa'
	n=0
	with open(out, 'w') as f:
		for line in open(a.ifile):
			if line[0]=='>':
				f.write('>'+str(n)+'\n')
				n+=1
			else:
				f.write(line)
				
				
def fasta_random_sample(a):
	if a.num:
		num=a.num
	else:
		num='100'
		
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_random'+num+'.fa'
		
	fa=list(parse(a.ifile,'fasta'))
	rand=random_integers(0,len(fa),int(num))
	
	n=0
	with open(out,'w') as f:
		for rec in fa:
			if n in rand:
				f.write('>'+rec.id  +'\n')
				f.write(str(rec.seq)+'\n')
			n+=1
				
	
def jmodeltest(a):
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_jmodeltest'
	if a.threads:
		threads=a.threads
	else:
		threads='1'
		
	if a.models:
		models=a.models
	else:
		models='7'
	
		
# java -jar ~/soft/phylo/jmodeltest2/dist/jModelTest.jar -d h1n1.fa		
# -i : invariant sites, f: unequal base freq, -g no of rate categoris	
	call([ 'java', '-jar', jmodeltest_path, '-d', a.ifile, '-o', out, '-tr', threads, '-s', models, '-i', '-f', '-AIC', '-g', '4' ])


def blastdb_make(a):
	if a.prot:
		call([ makeblastdb_path, '-dbtype', 'prot', '-in', a.ifile ])
	else:
		call([ makeblastdb_path, '-dbtype', 'nucl', '-in', a.ifile ])
		
		
def blast(a):
	out=a.ifile.split('.')[0]+'_blast'
#
	if a.prot:
		blast=blastp_path
	else:
		blast=blastn_path
		
#blastn -num_threads 1 -num_descriptions 5  -num_alignments 5  -db vp1ali99_1.fa -query test.fa -out test_blast	
	if a.threads:
		th=a.threads
	else:
		th='1'
		
	if a.num_descriptions:
		num_descriptions=a.num_descriptions
	else:
		num_descriptions='100'
		
	if a.num_alignments:
		num_alignments=a.num_alignments
	else:
		num_alignments='100'
	
	if a.max_target_seqs:
		max_target_seqs=a.max_target_seqs
	else:
		max_target_seqs='1'
	
	if a.tab:
		call([ blast, '-num_threads', th, '-max_target_seqs', max_target_seqs, '-db', a.db, '-query', a.ifile, '-out', out, '-outfmt', '7 qseqid sseqid pident slen sstart send evalue' ])
	else:
		call([ blast, '-num_threads', th, '-num_descriptions', num_descriptions, '-num_alignments', num_alignments, '-db', a.db, '-query', a.ifile, '-out', out])
		
		
def seqlogo_weblogo(a):
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_weblogo.png'
	
	call([weblogo_path, '-F', 'png', '-c', 'chemistry', '-n', '50', '--errorbars', 'NO', '--resolution', '300', '-f', a.ifile, '-o', out])
	

def cluster_cdhit(a):
	if a.pct:
		pct=a.pct
	else:
		pct='1'
		
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_cdhit_'+pct+'.fa'
	
	if a.protein:
		call([ cdhit_path, '-i', a.ifile, '-o', out, '-c', pct  ])
	else:
		call([ cdhit-est_path, '-i', a.ifile, '-o', out, '-c', pct  ])
		
		
def cluster_usearch(a):
	# usearch -cluster_fast $n.fa -id 0.9 -centroids $n\90.fa -uc $n\90.uc
	if a.pct:
		pct=a.pct
	else:
		pct='1'
		
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_usearch_'+pct+'.fa'
	
	call([ usearch_path, '-cluster_fast', a.ifile, '-centroids', out, '-id', pct, '-uc', out+'.clstr' ])
	
	
def hmmbuild(a):
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'.hmm'
		
	call([hmmbuild_path, out, a.ifile])
	
	
def hmmsearch(a):
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_hmmsearch'
		
	call([hmmsearch_path,'--tblout',out+'.tsv', '-o', out, a.ihmm, a.ifile ])
	
			
@Gooey(default_size=(800, 900),advanced=True) # (width,height)
def main():
	p=GooeyParser(description='Bioinformatics tools for sequence analysis and phylogenetics ')
#	p=ArgumentParser()
	sp=p.add_subparsers(help='commands', dest='command')
	

	sp1=sp.add_parser('seq_format_conversion')
	a1=sp1.add_argument
	a1('iformat',help='input format', default='fasta', choices=format)
	a1('oformat',help='output format', default='phylip', choices=format)
	a1('ifile',help='input file',  widget='FileChooser')
	a1('-ofile',help='output file, default: ifile.oformat', widget='FileChooser')
	sp1.set_defaults(func=seq_format_conversion)

	
	sp2=sp.add_parser('remove_alignment_gap')
	a2=sp2.add_argument
	a2('ifile',help='input fasta',  widget='FileChooser')
	a2('-ofile',help='output fasta, default: ifile_nogap.fa', widget='FileChooser')
	sp2.set_defaults(func=remove_alignment_gap)


	sp3=sp.add_parser('tree_fasttree')
	a3=sp3.add_argument
	a3('-protein', action='store_true', help='protein sequence')
	a3('ifile', help='fasta/phylip input file', widget='FileChooser')
	a3('-ofile', help='tree output file, default: ifile_ftree', widget='FileChooser')
	sp3.set_defaults(func=tree_fasttree) 
	

	sp4=sp.add_parser('alignment_muscle')
	a4=sp4.add_argument
	a4('-maxiters', help='maximum iterations, default: 16',choices='1 2 4 8 12 16'.split())
	a4('ifile', help='fasta input file', widget='FileChooser')
	a4('-ofile', help='fasta output file, default: ifile_muscle.fa', widget='FileChooser')
	sp4.set_defaults(func=alignment_muscle) 

	
	sp5=sp.add_parser('extract_geneid_fasta')
	a5=sp5.add_argument
	a5('ifile',help='input fasta',  widget='FileChooser')
	a5('geneid',help='geneid in fasta header')
	a5('-ofile',help='output fasta, default: ifile_geneid.fa', widget='FileChooser')
	sp5.set_defaults(func=extract_geneid_fasta)


	sp6=sp.add_parser('tree_raxml')
	a6=sp6.add_argument
	a6('-protein', action='store_true', help='protein sequence')
	a6('-bootstrap', help='# of bootstraps, default:no')
	a6('-threads',help='# of CPU threads, default: 1')
	a6('ifile', help='fasta/phylip input file', widget='FileChooser')
	sp6.set_defaults(func=tree_raxml) 
	

	sp7=sp.add_parser('tree_prep_ali')
	a7=sp7.add_argument
	a7('ifile',help='input fasta, Remove illegal characters ( ) [ ] ; : from alignment header for tree building',  widget='FileChooser')
	a7('-ofile',help='output fasta, default: ifile_4tree.fa', widget='FileChooser')
	sp7.set_defaults(func=tree_prep_ali)
	
	
	sp8=sp.add_parser('extract_genelist_fasta')
	a8=sp8.add_argument
	a8('ifile',help='input fasta',  widget='FileChooser')
	a8('genelist',help='geneids in fasta header',  widget='FileChooser')
	a8('-ofile',help='output fasta, default: ifile_genelist.fa', widget='FileChooser')
	sp8.set_defaults(func=extract_genelist_fasta)

	
	sp9=sp.add_parser('fasta_statistics')
	a9=sp9.add_argument
	a9('ifile',help='input fasta',  widget='FileChooser')
	a9('-ofile',help='output fasta, default: ifile_stat.png', widget='FileChooser')
	sp9.set_defaults(func=fasta_statistics)

	
	sp10=sp.add_parser('remove_seq_from_ali')
	a10=sp10.add_argument
	a10('ifile',help='input fasta, remove empty seq from alignment',  widget='FileChooser')
	a10('-ofile',help='output fasta, default: ifile_percent.png', widget='FileChooser')
	a10('-pctgap',help='remove seq over pecent of gap, default: 50')
	sp10.set_defaults(func=remove_seq_from_ali)

	
	sp11=sp.add_parser('remove_empty_lines')
	a11=sp11.add_argument
	a11('ifile',help='input fasta',  widget='FileChooser')
	a11('-ofile',help='output fasta, default: ifile_stat.png', widget='FileChooser')
	sp11.set_defaults(func=remove_empty_lines)

	
	sp12=sp.add_parser('jmodeltest')
	a12=sp12.add_argument
	a12('ifile',help='input fasta, rename fasta2numbers before using',  widget='FileChooser')
	a12('-ofile',help='output, default: ifile_jmodeltest', widget='FileChooser')
	a12('-threads',help='# of CPU threads, default: 1')
	a12('-models',help='# of substitution models to test, default: 7', choices='3 5 7 11 203'.split(),default='7')
	sp12.set_defaults(func=jmodeltest)

	
	sp13=sp.add_parser('rename_fasta2numbers')
	a13=sp13.add_argument
	a13('ifile',help='input fasta, renames fasta headers to numers',  widget='FileChooser')
	a13('-ofile',help='output fasta, default: ifile_rename2numbers.fa', widget='FileChooser')
	sp13.set_defaults(func=rename_fasta2numbers)
	
	
	sp14=sp.add_parser('fasta_random_sample')
	a14=sp14.add_argument
	a14('ifile',help='input fasta',  widget='FileChooser')
	a14('-ofile',help='output fasta, default: ifile_random_num.png', widget='FileChooser')
	a14('-num',help='# of random seq to sample, default: 100')
	sp14.set_defaults(func=fasta_random_sample)
	
	
	sp15=sp.add_parser('alignment_mafft')
	a15=sp15.add_argument
	a15('--maxiterate', help='maximum iterations, default: 0',choices='0 1 2 4 8 12 16'.split())
	a15('ifile', help='fasta input file', widget='FileChooser')
	a15('--thread',help='# of CPU threads, default: 1')
	a15('--ofile', help='fasta output file, default: ifile_mafft.fa', widget='FileChooser')
	a15('--path',help='mafft location',  widget='FileChooser')
	sp15.set_defaults(func=alignment_mafft) 

	
	sp16=sp.add_parser('blastdb_make')
	a16=sp16.add_argument
	a16('ifile',help='input fasta',  widget='FileChooser')
	a16('-prot', action='store_true', help='protein sequence')
	sp16.set_defaults(func=blastdb_make)

	
	sp17=sp.add_parser('blast')
	a17=sp17.add_argument
	a17('ifile',help='input fasta',  widget='FileChooser')
	a17('db', help='blast database',  widget='FileChooser')
	a17('-prot', action='store_true', help='protein sequence')
	a17('-num_alignments', default='100', help='# output alignments')
	a17('-num_descriptions', default='100', help='# output descriptions')
	a17('-threads', default='1', help='# threads')
	a17('-tab', action='store_true', help='tabular output: default:no')
	a17('-max_target_seqs', help='for tabular output only')
	sp17.set_defaults(func=blast)

	
	sp18=sp.add_parser('seqlogo_weblogo')
	a18=sp18.add_argument
	a18('ifile', help='fasta/other input file', widget='FileChooser')
	a18('-ofile', help='png/eps/svg output file, default: ifile_ftree', widget='FileChooser')
	sp18.set_defaults(func=seqlogo_weblogo) 

	
	sp19=sp.add_parser('cluster_cdhit')
	a19=sp19.add_argument
	a19('ifile', help='fasta input file', widget='FileChooser')
	a19('-ofile', help='fasta output file, default: ifile_cdhit_id', widget='FileChooser')
	a19('-protein', action='store_true', help='protein sequence')
	a19('-pct', help='sequence identity: 0-1, default: 1')
	sp19.set_defaults(func=cluster_cdhit) 

	
	sp20=sp.add_parser('cluster_usearch')
	a20=sp20.add_argument
	a20('ifile', help='fasta input file', widget='FileChooser')
	a20('-ofile', help='fasta output file, default: ifile_usearch_id', widget='FileChooser')
	a20('-protein', action='store_true', help='protein sequence')
	a20('-pct', help='sequence identity: 0-1, default: 1')
	sp20.set_defaults(func=cluster_usearch) 
	
	
	sp21=sp.add_parser('hmmbuild')
	a21=sp21.add_argument
	a21('ifile', help='fasta/other input alignment', widget='FileChooser')
	a21('-ofile', help='hmm output file, default: ifile.hmm', widget='FileChooser')
	sp21.set_defaults(func=hmmbuild) 
	
	
	sp22=sp.add_parser('hmmsearch')
	a22=sp22.add_argument
	a22('ifile', help='seq file to search', widget='FileChooser')
	a22('ihmm', help='input HMM obtained from hmmbuild', widget='FileChooser')
	a22('-ofile', help='hmm output file, default: ifile.hmm', widget='FileChooser')
	sp22.set_defaults(func=hmmsearch) 
	
	
	sp23=sp.add_parser('dna_translate')
	a23=sp23.add_argument
	a23('ifile', help='fasta input file', widget='FileChooser')
	a23('-ofile', help='translated output file, default: ifile_translate.fa', widget='FileChooser')
	sp23.set_defaults(func=dna_translate) 
	
	
	sp24=sp.add_parser('dna_complement')
	a24=sp24.add_argument
	a24('ifile', help='fasta input file', widget='FileChooser')
	a24('-ofile', help='output file, default: ifile_translate.fa', widget='FileChooser')
	sp24.set_defaults(func=dna_complement) 
	
	
	sp25=sp.add_parser('dna_reverse_complement')
	a25=sp25.add_argument
	a25('ifile', help='fasta input file', widget='FileChooser')
	a25('-ofile', help='output file, default: ifile_translate.fa', widget='FileChooser')
	sp25.set_defaults(func=dna_reverse_complement) 
	
	
	sp26=sp.add_parser('plot_table')
	a26=sp26.add_argument
	a26('ifile', help='table input file (columns separated by spaces or tab)', widget='FileChooser')
	a26('-xcol1', action='store_true', help='column 1 is x-axis')
	a26('-ofile', help='output file, default: ifile_plot_table.png', widget='FileChooser')
	sp26.set_defaults(func=plot_table) 
	
	

	
	a=p.parse_args()
		
	a.func(a)

	
	
if __name__ == '__main__':
	main()
