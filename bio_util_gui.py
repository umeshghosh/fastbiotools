#!/usr/bin/env python
# Copyright: Umesh Ghoshdastider
#from argparse import ArgumentParser
from Bio.SeqIO import convert,parse
import glob
from gooey import Gooey, GooeyParser
from pylab import array, hist, random_integers, savefig, title, xlabel, ylabel
import os
import shutil
from subprocess import call,Popen

# help='All input files need extension', file path must not include '.' in directory name, no gap expected in fasta file
# variables 
# a = input_args, p = parser, sp = subparser, a1 = sp1_add_argument

format='clustal embl fasta genbank nexus phylip stockholm tab'.split()
	
def seq_format_conv(a):
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
				
def muscle_alignment(a):
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_muscle.fa'
	
	cmd=['muscle','-out', out,'-in', a.ifile]
	if a.maxiters:
		cmd.append('-maxiters')
		cmd.append(a.maxiters)
		call(cmd)
	else:
		call(cmd)

def mafft_alignment(a):
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_mafft.fa'
	
	if a.maxiterate:
		iter=a.maxiterate
	else:
		iter='0'
	
	if a.path:
		path=a.path
	else:
		path='/home/umesh/soft/msa/mafft-7.266-with-extensions/build/bin/mafft'
	if a.thread:
		th=a.thread
	else:
		th='1'
		
	with open(out,'w') as f:
		Popen([path, '--thread', th, '--maxiterate', iter, a.ifile], stdout=f)

def fasttree(a):
	if a.ofile:
		out=a.ofile
	else:
		out=a.ifile.split('.')[0]+'_ftree'
		
	if a.protein:
		call(['fasttree','-wag', '-gamma', '-out', out, a.ifile])
	else:
		call(['fasttree','-nt', '-gtr', '-gamma', '-out', out, a.ifile])		
		
def raxml(a):
# -N 2 : two independent runs
	f=os.path.basename(a.ifile).split('.')[0]
	if a.protein:
		call(['raxml', '-m', 'PROTGAMMAILG',  '-p', '12345', '-N', '2', '-s', a.ifile, '-n', f ])
	else:
		call(['raxml', '-m', 'GTRGAMMAI', '-p', '12345', '-N', '2', '-s', a.ifile, '-n', f ])
				
	dir='RAxML_'+f
	os.mkdir(dir)
	for fi in glob.glob('RAxML*'+f+'*'):
		shutil.move(fi, dir)

def prep_ali4tree(a): 
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
			
def fasta_stat(a):
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

	savefig(out)
	
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
	
	if a.path:
		path=a.path
	else:
		path='/home/umesh/soft/phylo/jmodeltest2/dist/jModelTest.jar'
			
# java -jar ~/soft/phylo/jmodeltest2/dist/jModelTest.jar -d h1n1.fa		
# -i : invariant sites, f: unequal base freq, -g no of rate categoris	
	call([ 'java', '-jar', '/home/umesh/soft/phylo/jmodeltest2/dist/jModelTest.jar', '-d', a.ifile, '-o', out, '-tr', threads, '-s', models, '-i', '-f', '-AIC', '-g', '4' ])


def make_blastdb(a):
	if a.prot:
		call([ 'makeblastdb', '-dbtype', 'prot', '-in', a.ifile ])
	else:
		call([ 'makeblastdb', '-dbtype', 'nucl', '-in', a.ifile ])
		
		
def blast(a):
	out=a.ifile.split('.')[0]+'_blast'
#
	if a.prot:
		blast='blastp'
	else:
		blast='blastn'
		
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
	
	call([ blast, '-num_threads', th, '-num_descriptions', num_descriptions, '-num_alignments', num_alignments, '-db', a.db, '-query', a.ifile, '-out', out ])
		
			
@Gooey(default_size=(800, 700)) # (width,height)
def main():
	p=GooeyParser()
#	p=ArgumentParser()
	sp=p.add_subparsers(help='commands', dest='command')
	
	sp1=sp.add_parser('seq_format_conv')
	a1=sp1.add_argument
	a1('iformat',help='input format', default='fasta', choices=format)
	a1('oformat',help='output format', default='phylip', choices=format)
	a1('ifile',help='input file',  widget='FileChooser')
	a1('-ofile',help='output file, default: ifile.oformat', widget='FileChooser')
	sp1.set_defaults(func=seq_format_conv)
	
	sp2=sp.add_parser('remove_alignment_gap')
	a2=sp2.add_argument
	a2('ifile',help='input fasta',  widget='FileChooser')
	a2('-ofile',help='output fasta, default: ifile_nogap.fa', widget='FileChooser')
	sp2.set_defaults(func=remove_alignment_gap)

	sp3=sp.add_parser('fasttree')
	a3=sp3.add_argument
	a3('-protein', action='store_true', help='protein sequence')
	a3('ifile', help='fasta/phylip input file', widget='FileChooser')
	a3('-ofile', help='tree output file, default: ifile_ftree', widget='FileChooser')
	sp3.set_defaults(func=fasttree) 
	
	sp4=sp.add_parser('muscle_alignment')
	a4=sp4.add_argument
	a4('-maxiters', help='maximum iterations, default: 16',choices='1 2 4 8 12 16'.split())
	a4('ifile', help='fasta input file', widget='FileChooser')
	a4('-ofile', help='fasta output file, default: ifile_muscle.fa', widget='FileChooser')
	sp4.set_defaults(func=muscle_alignment) 
	
	sp5=sp.add_parser('extract_geneid_fasta')
	a5=sp5.add_argument
	a5('ifile',help='input fasta',  widget='FileChooser')
	a5('geneid',help='geneid in fasta header')
	a5('-ofile',help='output fasta, default: ifile_geneid.fa', widget='FileChooser')
	sp5.set_defaults(func=extract_geneid_fasta)

	sp6=sp.add_parser('raxml')
	a6=sp6.add_argument
	a6('-protein', action='store_true', help='protein sequence')
	a6('ifile', help='fasta/phylip input file', widget='FileChooser')
	sp6.set_defaults(func=raxml) 
	
	sp7=sp.add_parser('prep_alignment_4tree')
	a7=sp7.add_argument
	a7('ifile',help='input fasta, Remove illegal characters ( ) [ ] ; : from alignment header for tree building',  widget='FileChooser')
	a7('-ofile',help='output fasta, default: ifile_4tree.fa', widget='FileChooser')
	sp7.set_defaults(func=prep_ali4tree)
	
	
	sp8=sp.add_parser('extract_genelist_fasta')
	a8=sp8.add_argument
	a8('ifile',help='input fasta',  widget='FileChooser')
	a8('genelist',help='geneids in fasta header',  widget='FileChooser')
	a8('-ofile',help='output fasta, default: ifile_genelist.fa', widget='FileChooser')
	sp8.set_defaults(func=extract_genelist_fasta)
	
	sp9=sp.add_parser('fasta_stat')
	a9=sp9.add_argument
	a9('ifile',help='input fasta',  widget='FileChooser')
	a9('-ofile',help='output fasta, default: ifile_stat.png', widget='FileChooser')
	sp9.set_defaults(func=fasta_stat)
	
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
	a12('ifile',help='input fasta, rename fasta2numers before using',  widget='FileChooser')
	a12('-path',help='jModeltest.jar location',  widget='FileChooser')
	a12('-ofile',help='output, default: ifile_jmodeltest', widget='FileChooser')
	a12('-threads',help='# of CPU threads, default: 1')
	a12('-models',help='# of substitution models to test, default: 7', choices='3 5 7 11 203'.split(),default='7')
	sp12.set_defaults(func=jmodeltest)
	
	sp13=sp.add_parser('rename_fasta2numbers')
	a13=sp13.add_argument
	a13('ifile',help='input fasta',  widget='FileChooser')
	a13('-ofile',help='output fasta, default: ifile_rename2numbers.fa', widget='FileChooser')
	sp13.set_defaults(func=rename_fasta2numbers)
	
	
	sp14=sp.add_parser('fasta_random_sample')
	a14=sp14.add_argument
	a14('ifile',help='input fasta',  widget='FileChooser')
	a14('-ofile',help='output fasta, default: ifile_random_num.png', widget='FileChooser')
	a14('-num',help='# of random seq to sample, default: 100')
	sp14.set_defaults(func=fasta_random_sample)
	
	sp15=sp.add_parser('mafft_alignment')
	a15=sp15.add_argument
	a15('--maxiterate', help='maximum iterations, default: 0',choices='0 1 2 4 8 12 16'.split())
	a15('ifile', help='fasta input file', widget='FileChooser')
	a15('--thread',help='# of CPU threads, default: 1')
	a15('--ofile', help='fasta output file, default: ifile_mafft.fa', widget='FileChooser')
	a15('--path',help='mafft location',  widget='FileChooser')
	sp15.set_defaults(func=mafft_alignment) 
	
	sp16=sp.add_parser('make_blastdb')
	a16=sp16.add_argument
	a16('ifile',help='input fasta',  widget='FileChooser')
	a16('-prot', action='store_true', help='protein sequence')
	sp16.set_defaults(func=make_blastdb)
	
	sp17=sp.add_parser('blast')
	a17=sp17.add_argument
	a17('ifile',help='input fasta',  widget='FileChooser')
	a17('db', help='blast database',  widget='FileChooser')
	a17('-prot', action='store_true', help='protein sequence')
	a17('-num_alignments', default='100', help='# output alignments')
	a17('-num_descriptions', default='100', help='# output alignments')
	a17('-threads', default='1', help='# threads')
	sp17.set_defaults(func=blast)
	
	
	a=p.parse_args()
	
	a.func(a)

	
if __name__ == '__main__':
	main()