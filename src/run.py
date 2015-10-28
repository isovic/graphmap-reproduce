#! /usr/bin/python

import re;
import os;
import sys;
import subprocess;
import multiprocessing;
import fnmatch
import numpy as np;

SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));
tools_path = '%s/../tools' % (SCRIPT_PATH);
SAMSCRIPTS = '%s/../tools/samscripts/src/' % (SCRIPT_PATH);

# sys.path.append('%s/../tools/samscripts/src/' % (SCRIPT_PATH));
SAMSCRIPTS = '%s/../../samscripts/src/' % (SCRIPT_PATH);
sys.path.append('%s/../../samscripts/src/' % (SCRIPT_PATH));
try:
	import vcffilter;
	import fastqparser;
	import utility_sam;
	MODULE_VCFFILTER = True;
	MODULE_FASTQPARSER = True;
	MODULE_UTILITYSAM = True;
except:
	MODULE_VCFFILTER = False;
	MODULE_FASTQPARSER = False;
	MODULE_UTILITYSAM = False;

	# sys.stderr.write('Please run "%s setup" first! Exiting.\n\n' % (sys.argv[0]));
	# exit(1);

BAMSURGEON_PATH = '%s/../tools/bamsurgeon' % (SCRIPT_PATH);

def setup_tools():
	execute_command('sudo pip install pyvcf');
	execute_command('sudo apt-get install tabix');
	execute_command('sudo apt-get installa vcftools');
	execute_command('cd %s/../tools/; wget http://netassist.dl.sourceforge.net/project/lofreq/lofreq_star-2.1.2_linux-x86-64.tgz; tar xvf lofreq_star-2.1.2_linux-x86-64.tgz' % (SCRIPT_PATH));

	if (not os.path.exists('%s/../tools' % (SCRIPT_PATH))):
		execute_command('mkdir %s/../tools' % (SCRIPT_PATH));

	if (not os.path.exists('%s/../tools/samscripts' % (SCRIPT_PATH))):
		execute_command('cd %s/../tools; git clone https://github.com/isovic/samscripts.git' % (SCRIPT_PATH));
	
	if (not os.path.exists('%s/../tools/aligneval' % (SCRIPT_PATH))):
		execute_command('cd %s/../tools; git clone https://github.com/isovic/aligneval.git; cd aligneval; git checkout devel; ./setup.py aligners; ./setup.py tools' % (SCRIPT_PATH));

	if (not os.path.exists('%s/../tools/mutatrix/' % (SCRIPT_PATH))):
		# execute_command('cd %s/../packs; tar -xvf mutatrix.tar.gz' % (SCRIPT_PATH));
		# execute_command('mv packs/mutatrix tools/' % (SCRIPT_PATH));
		execute_command('cd %s/../tools; git clone --recursive https://github.com/ekg/mutatrix.git; cd mutatrix; git checkout 25d12240ec59b345da7f73d3fb5cdd3422a5ed46; git submodule update --init --recursive; make' % (SCRIPT_PATH));

	if (not os.path.exists('%s/../tools/bamsurgeon/' % (SCRIPT_PATH))):
		execute_command('cd %s/../tools; git clone --recursive https://github.com/adamewing/bamsurgeon.git' % (SCRIPT_PATH));




def setup_data():
	sys.stderr.write('Generating simulated data.\n');
	execute_command('cd %s/../tools; cd aligneval; ./setup.py simdata' % (SCRIPT_PATH));

	sys.stderr.write('Creating downloads folder in path: "%s/../data/downloads".\n' % (SCRIPT_PATH));
	os.makedirs('%s/../data/downloads' % (SCRIPT_PATH));

	sys.stderr.write('Downloading raw nanopore data.\n');
	
	sys.stderr.write('\tFetching the R7 dataset.\n');
	execute_command('cd %s/../data/downloads/; mkdir R7' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd R7; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA362/ERA362836/oxfordnanopore_native/Ecoli_R7_NONI.tgz' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd R7; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA362/ERA362836/oxfordnanopore_native/Ecoli_R7_ONI_flowcell_18.tar.gz' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd R7; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA362/ERA362836/oxfordnanopore_native/Ecoli_R7_ONI_flowcell_17.tar.gz' % (SCRIPT_PATH));

	sys.stderr.write('\tFetching the R7.3 data.\n');
	execute_command('cd %s/../data/downloads/; mkdir R7.3' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd R7.3; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA362/ERA362836/oxfordnanopore_native/Ecoli_R73.tgz' % (SCRIPT_PATH));

	sys.stderr.write('\tFetching the E. Coli UTI89 dataset, generated in-house.\n');	
	execute_command('cd %s/../data/downloads/; mkdir ecoli-uti89' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd ecoli-uti89; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA444/ERA444730/oxfordnanopore_native/reads.tar' % (SCRIPT_PATH));

	sys.stderr.write('\tFetching the Salmonella Typhi dataset.\n');	
	execute_command('cd %s/../data/downloads/; mkdir salmonella-typhi' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd salmonella-typhi; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA375/ERA375685/oxfordnanopore_native/H566_ON_inc.tar.gz' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd salmonella-typhi; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA375/ERA375987/oxfordnanopore_native/H566_30_min_inc.tar.gz' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd salmonella-typhi; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA376/ERA376255/oxfordnanopore_native/raw_2_rabsch_R7.tar.gz' % (SCRIPT_PATH));

	sys.stderr.write('\tFetching the amplicon sequencing dataset.\n');	
	execute_command('cd %s/../data/downloads/; mkdir amplicons' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd amplicons; wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR174/004/SRR1747434/SRR1747434.fastq.gz' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd amplicons; wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR174/000/SRR1748410/SRR1748410.fastq.gz' % (SCRIPT_PATH));

	sys.stderr.write('\tFetching the BE1 dataset.\n');	
	execute_command('cd %s/../data/downloads/; mkdir be1' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd be1; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA463/ERA463589/oxfordnanopore_native/FAA37759_GB2974_MAP005_20150423__2D_basecalling_v1.14_2D.tar.gz' % (SCRIPT_PATH));

	sys.stderr.write('\tFetching the ADP1 dataset.\n');	
	execute_command('cd %s/../data/downloads/; mkdir adp1' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd adp1; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA416/ERA416080/fastq/AWK_ONT_MN2064525.fastq.gz' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd adp1; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA416/ERA416080/fastq/AWK_ONT_MN2064006.fastq.gz' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd adp1; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA416/ERA416080/fastq/AWK_ONT_FAA43210.fastq.gz' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd adp1; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA416/ERA416080/fastq/AWK_ONT_FAA43204.fastq.gz' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd adp1; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA416/ERA416080/fastq/AWK_ONT_FAA17573.fastq.gz' % (SCRIPT_PATH));
	
	sys.stderr.write('\tFetching the Loman et al. dataset used for nanopore-only assembly.\n');
	execute_command('cd %s/../data/downloads/; mkdir ecoli-nmeth' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd ecoli-nmeth; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_20_LomanLabz_PC_Ecoli_K12_R7.3.tar' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd ecoli-nmeth; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_32_LomanLabz_K12_His_tag.tar' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd ecoli-nmeth; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_33_LomanLabz_PC_K12_0.4SPRI_Histag.tar' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd ecoli-nmeth; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_39.tar' % (SCRIPT_PATH));

	sys.stderr.write('\tFetching the MAP006-1 dataset.\n');
	execute_command('cd %s/../data/downloads/; mkdir ecoli-map006' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd ecoli-map006; wget http://nanopore.climb-radosgw01.bham.ac.uk/MAP006-1.basecalled.tar' % (SCRIPT_PATH));



def run_simulated_data():
	pass;

def run_real_data():
	RUN_CONSENSUS_TEST_LAMBDA();
	RUN_CONSENSUS_TEST_ECOLIR7();
	RUN_CONSENSUS_TEST_ECOLIR73();
	RUN_CONSENSUS_TEST_UTI89();
	RUN_CONSENSUS_TEST_TYPHI();
	RUN_CONSENSUS_TEST_ADP1();
	RUN_CONSENSUS_TEST_BE1();
	RUN_CONSENSUS_TEST_ECOLINMETH();
	RUN_CONSENSUS_TEST_MAP006();
	RUN_CONSENSUS_TEST_ECOLIR73_1d2d();

	RUN_MUTATED_REFERENCE_ADDITIONAL_TESTS(True);
	
	RUN_SV_TEST();
	RUN_AMPLICON_TEST();
	RUN_AMPLICON_TEST_1d2d();
	RUN_AMPLICON_TEST_GRAPHMAP_VS_BWAMEM();

	RUN_DRAFT_ASSEMBLY_REFERENCE_TESTS();

	RUN_VENN_COMPARISON_TEST('data/out/consensus-ecoliR7.3/GraphMap-ecoliR7.3.sam',
							'data/out/consensus-ecoliR7.3/LAST-ecoliR7.3.sam',
							'data/consensus-ecoliR7.3/reference/escherichia_coli.fa',
							'data/out/consensus-ecoliR7.3/analysis-venn/venn');

#########################################################
#########################################################
##### Consensus tests #####
#########################################################
#########################################################

def RUN_CONSENSUS_TEST_LAMBDA():
	reference_path = ('%s/../data/consensus-lambdaR6/reference/NC_001416.fa' % SCRIPT_PATH);
	reads_path = ('%s/../data/consensus-lambdaR6/reads/lambda_reads-enumerated.fastq' % SCRIPT_PATH)
	dataset_name = 'lambdaR6';
	out_path = '%s/../data/out/consensus-%s/' % (SCRIPT_PATH, dataset_name);
	mappers_to_run = ['daligner', 'graphmap', 'graphmap-anchor', 'last', 'bwamem', 'blasr', 'marginalign', 'marginaligngraphmap', 'marginaligngraphmap-anchor']

	run_all_mappers_only(reference_path, reads_path, dataset_name, out_path, 'nanopore', do_not_recalc=True, is_circular=True, select_mappers=mappers_to_run);
	evaluate_alignments(reference_path, reads_path, dataset_name, out_path, cov_threshold=20);
	# collect_alignments(reference_path, reads_path, dataset_name, out_path, cov_thresholds=[20]);

def RUN_CONSENSUS_TEST_ECOLIR7():
	reference_path = ('%s/../data/consensus-ecoliR7/reference/escherichia_coli.fa' % SCRIPT_PATH);
	reads_path = ('%s/../data/consensus-ecoliR7/reads/Ecoli_R7_CombinedFasta.fastq' % SCRIPT_PATH);
	dataset_name = 'ecoliR7';
	out_path = '%s/../data/out/consensus-%s/' % (SCRIPT_PATH, dataset_name);
	mappers_to_run = ['daligner', 'graphmap', 'graphmap-anchor', 'last', 'bwamem', 'blasr', 'marginalign', 'marginaligngraphmap', 'marginaligngraphmap-anchor']

	run_all_mappers_only(reference_path, reads_path, dataset_name, out_path, 'nanopore', do_not_recalc=True, is_circular=True, select_mappers=mappers_to_run);
	evaluate_alignments(reference_path, reads_path, dataset_name, out_path, cov_threshold=20);
	# collect_alignments(reference_path, reads_path, dataset_name, out_path, cov_thresholds=[20]);

def RUN_CONSENSUS_TEST_ECOLIR73():
	reference_path = ('%s/../data/consensus-ecoliR7.3/reference/escherichia_coli.fa' % SCRIPT_PATH);
	reads_path = ('%s/../data/consensus-ecoliR7.3/reads/ecoliR7.3.fastq' % SCRIPT_PATH)
	dataset_name = 'ecoliR7.3';
	out_path = '%s/../data/out/consensus-%s/' % (SCRIPT_PATH, dataset_name);
	mappers_to_run = ['daligner', 'graphmap', 'graphmap-anchor', 'last', 'bwamem', 'blasr', 'marginalign', 'marginaligngraphmap', 'marginaligngraphmap-anchor']

	run_all_mappers_only(reference_path, reads_path, dataset_name, out_path, 'nanopore', do_not_recalc=True, is_circular=True, select_mappers=mappers_to_run);
	evaluate_alignments(reference_path, reads_path, dataset_name, out_path, cov_threshold=20);
	# collect_alignments(reference_path, reads_path, dataset_name, out_path, cov_thresholds=[20]);

def RUN_CONSENSUS_TEST_UTI89():
	reference_path = ('%s/../data/consensus-uti89/reference/NC_007946.1.fa' % SCRIPT_PATH);
	reads_path = ('%s/../data/consensus-uti89/reads/reads_ecoli_uti89.fastq' % SCRIPT_PATH)
	dataset_name = 'uti89';
	out_path = '%s/../data/out/consensus-%s/' % (SCRIPT_PATH, dataset_name);
	mappers_to_run = ['daligner', 'graphmap', 'graphmap-anchor', 'last', 'bwamem', 'blasr', 'marginalign', 'marginaligngraphmap', 'marginaligngraphmap-anchor']

	run_all_mappers_only(reference_path, reads_path, dataset_name, out_path, 'nanopore', do_not_recalc=True, is_circular=True, select_mappers=mappers_to_run);
	evaluate_alignments(reference_path, reads_path, dataset_name, out_path, cov_threshold=10);
	# collect_alignments(reference_path, reads_path, dataset_name, out_path, cov_thresholds=[20]);

def RUN_CONSENSUS_TEST_TYPHI():
	reference_path = ('%s/../data/consensus-salmonellatyphi/reference/salmonella_typhi_Ty2.fa' % SCRIPT_PATH);
	reads_path = ('%s/../data/consensus-salmonellatyphi/reads/reads-typhi.fastq' % SCRIPT_PATH)
	dataset_name = 'salmonellatyphi';
	out_path = '%s/../data/out/consensus-%s/' % (SCRIPT_PATH, dataset_name);
	mappers_to_run = ['daligner', 'graphmap', 'graphmap-anchor', 'last', 'bwamem', 'blasr', 'marginalign', 'marginaligngraphmap', 'marginaligngraphmap-anchor']

	run_all_mappers_only(reference_path, reads_path, dataset_name, out_path, 'nanopore', do_not_recalc=True, is_circular=True, select_mappers=mappers_to_run);
	evaluate_alignments(reference_path, reads_path, dataset_name, out_path, cov_threshold=20);
	# collect_alignments(reference_path, reads_path, dataset_name, out_path, cov_thresholds=[20]);

def RUN_CONSENSUS_TEST_ADP1():
	reference_path = ('%s/../data/consensus-ADP1/reference/acinetobacter_baylyi_ADP1_NC_005966.1.fa' % SCRIPT_PATH);
	reads_path = ('%s/../data/consensus-ADP1/reads/all_1d2d.fastq' % SCRIPT_PATH)
	dataset_name = 'ADP1';
	out_path = '%s/../data/out/consensus-%s/' % (SCRIPT_PATH, dataset_name);
	mappers_to_run = ['daligner', 'graphmap', 'graphmap-anchor', 'last', 'bwamem', 'blasr', 'marginalign', 'marginaligngraphmap', 'marginaligngraphmap-anchor']

	run_all_mappers_only(reference_path, reads_path, dataset_name, out_path, 'nanopore', do_not_recalc=True, is_circular=True, select_mappers=mappers_to_run);
	evaluate_alignments(reference_path, reads_path, dataset_name, out_path, cov_threshold=20);
	# collect_alignments(reference_path, reads_path, dataset_name, out_path, cov_thresholds=[20]);

def RUN_CONSENSUS_TEST_BE1():
	reference_path = ('%s/../data/consensus-BE1/reference/bacteroides_fragilis-BFBE1.fa' % SCRIPT_PATH);
	reads_path = ('%s/../data/consensus-BE1/reads/reads-BE1.fastq' % SCRIPT_PATH)
	dataset_name = 'BE1';
	out_path = '%s/../data/out/consensus-%s/' % (SCRIPT_PATH, dataset_name);
	mappers_to_run = ['daligner', 'graphmap', 'graphmap-anchor', 'last', 'bwamem', 'blasr', 'marginalign', 'marginaligngraphmap', 'marginaligngraphmap-anchor']

	run_all_mappers_only(reference_path, reads_path, dataset_name, out_path, 'nanopore', do_not_recalc=True, is_circular=True, select_mappers=mappers_to_run);
	evaluate_alignments(reference_path, reads_path, dataset_name, out_path, cov_threshold=20);
	# collect_alignments(reference_path, reads_path, dataset_name, out_path, cov_thresholds=[20]);

def RUN_CONSENSUS_TEST_ECOLINMETH():
	reference_path = ('%s/../data/consensus-ecolinmeth/reference/escherichia_coli.fa' % SCRIPT_PATH);
	reads_path = ('%s/../data/consensus-ecolinmeth/reads/reads-nmeth-all_2d.fastq' % SCRIPT_PATH)
	dataset_name = 'ecolinmeth';
	out_path = '%s/../data/out/consensus-%s/' % (SCRIPT_PATH, dataset_name);
	mappers_to_run = ['daligner', 'graphmap', 'graphmap-anchor', 'last', 'bwamem', 'blasr', 'marginalign', 'marginaligngraphmap', 'marginaligngraphmap-anchor']

	run_all_mappers_only(reference_path, reads_path, dataset_name, out_path, 'nanopore', do_not_recalc=False, is_circular=True, select_mappers=mappers_to_run);
	evaluate_alignments(reference_path, reads_path, dataset_name, out_path, cov_threshold=20);
	# collect_alignments(reference_path, reads_path, dataset_name, out_path, cov_thresholds=[20]);

def RUN_CONSENSUS_TEST_MAP006():
	reference_path = ('%s/../data/consensus-ecoliMAP006/reference/escherichia_coli.fa' % SCRIPT_PATH);
	reads_path = ('%s/../data/consensus-ecoliMAP006/reads/reads-MAP006-1.fastq' % SCRIPT_PATH)
	dataset_name = 'ecoliMAP006-1';
	out_path = '%s/../data/out/consensus-%s/' % (SCRIPT_PATH, dataset_name);
	mappers_to_run = ['daligner', 'graphmap', 'graphmap-anchor', 'last', 'bwamem', 'blasr', 'marginalign', 'marginaligngraphmap', 'marginaligngraphmap-anchor']

	run_all_mappers_only(reference_path, reads_path, dataset_name, out_path, 'nanopore', do_not_recalc=False, is_circular=True, select_mappers=mappers_to_run);
	evaluate_alignments(reference_path, reads_path, dataset_name, out_path, cov_threshold=20);
	# collect_alignments(reference_path, reads_path, dataset_name, out_path, cov_thresholds=[20]);

### This is used to estimate the error rates in 1d and 2d data separately.
def RUN_CONSENSUS_TEST_ECOLIR73_1d2d():
	reference_path = ('%s/../data/consensus-ecoliR7.3-1d/reference/escherichia_coli.fa' % SCRIPT_PATH);
	reads_path = ('%s/../data/consensus-ecoliR7.3-1d/reads/ecoliR7.3-1d.fastq' % SCRIPT_PATH)
	dataset_name = 'ecoliR7.3-1d';
	out_path = '%s/../data/out/consensus-%s/' % (SCRIPT_PATH, dataset_name);
	mappers_to_run = ['daligner', 'graphmap', 'graphmap-anchor', 'last', 'bwamem', 'blasr', 'marginalign', 'marginaligngraphmap', 'marginaligngraphmap-anchor']
	run_all_mappers_only(reference_path, reads_path, dataset_name, out_path, 'nanopore', do_not_recalc=True, is_circular=True, select_mappers=mappers_to_run);
	evaluate_alignments(reference_path, reads_path, dataset_name, out_path, cov_threshold=20);

	reference_path = ('%s/../data/consensus-ecoliR7.3-2d/reference/escherichia_coli.fa' % SCRIPT_PATH);
	reads_path = ('%s/../data/consensus-ecoliR7.3-2d/reads/ecoliR7.3-2d.fastq' % SCRIPT_PATH)
	dataset_name = 'ecoliR7.3-2d';
	out_path = '%s/../data/out/consensus-%s/' % (SCRIPT_PATH, dataset_name);
	mappers_to_run = ['daligner', 'graphmap-anchor', 'last', 'bwamem', 'blasr', 'marginalign', 'marginaligngraphmap', 'marginaligngraphmap-anchor']
	run_all_mappers_only(reference_path, reads_path, dataset_name, out_path, 'nanopore', do_not_recalc=True, is_circular=True, select_mappers=mappers_to_run);
	evaluate_alignments(reference_path, reads_path, dataset_name, out_path, cov_threshold=20);



#########################################################
#########################################################
##### Mutated reference test #####
#########################################################
#########################################################

### This function mutates the E. Coli reference to contain an approximate number of SNPs as the draft assembly in Loman et al. "A complete bacterial genome assembled de novo using only nanopore sequencing data".
def RUN_MUTATED_REFERENCE_ADDITIONAL_TESTS(run_test=False):
	if (run_test == False):
		sys.stderr.write('*NOT* running the RUN_MUTATED_REFERENCE_ADDITIONAL_TESTS.\n');
		return;

	### This mutates the reference to include similar number of SNPs and indels as Loman/Simpson nanopore-only assembly pipeline (~3750 SNPs and ~42500 indels).
	generate_mutated_reference(('%s/../data/reference/escherichia_coli.fa' % SCRIPT_PATH), 0.0006, 0.0067, 'data/mutated-refs/draftlike');

	### This should run the alignment and collection of results on the entire E. Coli R7.3 dataset (1d+2d reads), on the artificially mutated reference.
	reference_path = ('%s/../data/reference/escherichia_coli.fa' % SCRIPT_PATH);
	assembly_draft = ('%s/../data/mutated-refs/draftlike_for_R7.3/mutated_escherichia_coli_snp0.000600_indel0.006700.fa' % SCRIPT_PATH);
	reads_path = ('%s/../data/reads-ecoliR7.3/ecoliR7.3.fastq' % SCRIPT_PATH);
	dataset_name = 'mutated_ref_draftlike_ecoliR7.3';
	out_path = '%s/../data/out/%s/' % (SCRIPT_PATH, dataset_name);
 
	run_all_mappers_only(assembly_draft, reads_path, dataset_name, out_path, 'nanopore', do_not_recalc=True, is_circular=True,
						select_mappers=['daligner', 'graphmap', 'graphmap-anchor', 'last', 'bwamem', 'blasr', 'marginalign', 'marginaligngraphmap', 'marginaligngraphmap-anchor']);
	# evaluate_alignments(assembly_draft, reads_path, dataset_name, out_path);
	# evaluate_unique_alignments(assembly_draft, reads_path, dataset_name, out_path);
	evaluate_consensus_sequences(reference_path, assembly_draft, dataset_name, out_path, 20);

	### This should run the alignment and collection of results on the Nature Methods assembly reads from the Loman et. al paper, on the artificially mutated reference.
	reference_path = ('%s/../data/reference/escherichia_coli.fa' % SCRIPT_PATH);
	assembly_draft = ('%s/../data/mutated-refs/draftlike_for_nmeth/mutated_escherichia_coli_snp0.000600_indel0.006700.fa' % SCRIPT_PATH);
	reads_path = ('%s/../data/nmeth/reads/reads-nmeth-all_2d.fastq' % SCRIPT_PATH);
	dataset_name = 'mutated_ref_draftlike_ecolinmeth';
	out_path = '%s/../data/out/%s/' % (SCRIPT_PATH, dataset_name);
 
	run_all_mappers_only(assembly_draft, reads_path, dataset_name, out_path, 'nanopore', do_not_recalc=True, is_circular=True,
						select_mappers=['daligner', 'graphmap', 'graphmap-anchor', 'last', 'bwamem', 'blasr', 'marginalign', 'marginaligngraphmap', 'marginaligngraphmap-anchor']);
	# evaluate_alignments(assembly_draft, reads_path, dataset_name, out_path);
	# evaluate_unique_alignments(assembly_draft, reads_path, dataset_name, out_path);
	evaluate_consensus_sequences(reference_path, assembly_draft, dataset_name, out_path, 20);



#########################################################
#########################################################
##### Structural Variant tests #####
#########################################################
#########################################################

def RUN_SV_TEST():
	original_reference = ('%s/../data/sv/original_reference/escherichia_coli.fa' % SCRIPT_PATH);
	mutated_reference = ('%s/../data/sv/mutated_reference/escherichia_coli-indel_events.fa' % SCRIPT_PATH);
	reads_path = ('%s/../data/sv/reads/all_2d_for_sv.fastq' % SCRIPT_PATH);
	
	## This is commented out just so it doesn't overwrite the files, because the evaluation is currently being tested.
	### First run the mappers on the original reference, to detect the differences that normaly exist and need to be omitted from further comparisons.
	run_all_mappers_only(original_reference, reads_path, 'all_2d_for_sv', '%s/../data/out/sv-normal_ref/' % (SCRIPT_PATH),
				'nanopore', do_not_recalc=True, is_circular=True,
				select_mappers=['daligner', 'graphmap', 'graphmap-anchor', 'last', 'bwamem', 'blasr', 'marginalign', 'marginaligngraphmap', 'marginaligngraphmap-anchor']);
	### Second, run the mappers on the modified reference.
	run_all_mappers_only(mutated_reference, reads_path, 'all_2d_for_sv', '%s/../data/out/sv-indel_ref/' % (SCRIPT_PATH),
				'nanopore', do_not_recalc=True, is_circular=True,
				select_mappers=['daligner', 'graphmap', 'graphmap-anchor', 'last', 'bwamem', 'blasr', 'marginalign', 'marginaligngraphmap', 'marginaligngraphmap-anchor']);

	### Detect the SV events.
	evaluate_sv_consensus(original_reference, 'all_2d_for_sv', '%s/../data/out/sv-normal_ref/' % (SCRIPT_PATH));
	evaluate_sv_consensus(mutated_reference, 'all_2d_for_sv', '%s/../data/out/sv-indel_ref/' % (SCRIPT_PATH));

	#############################################
	##### This part handles the evaluation. #####
	#############################################
	### First, load up the truth events. These are the indel events that were introduced to the reference,
	### and should be picked up by the mappers.
	fp = open('data/sv/sv_truth_events/truth_events-full.csv', 'r');
	truth_lines = [line.strip() for line in fp.readlines()];
	fp.close();

	### Create a lookup table for the offsets. The blacklisted events will be offset compared to the events
	### that are detected on the modified reference (because blacklisted events are detected on the original
	### reference that might be longer/shorter than the modified one).
	real_events = [];
	for line in truth_lines[1:]:
		line = line.strip();
		split_line = line.split('\t');
		### Create a list of events, where each event is described with: [end_pos, length, type, offset].
		real_events.append([int(split_line[2]), int(split_line[3]), split_line[0], 0]);
	sorted_real_events = sorted(real_events, key=lambda x: x[0]);
	current_offset = 0;
	for event in sorted_real_events:
		[pos, length, event_type, offset] = event;
		print event;
		current_offset += (length if (event_type == 'deletion_in_read') else (-length));
		# current_offset -= length if (event_type == 'deletion_in_read') else (-length);
		event[3] = current_offset;
		event[0] += current_offset;
	# print sorted_real_events;
	print '';
	for event in sorted_real_events:
		print event;
	print '';

	### Load up all the 'blacklisted' events.
	files_for_blacklist = find_files('data/out/sv-normal_ref/analysis-final/', '*chained.csv', -1);
	# files_for_blacklist = find_files('data/out/sv-normal_ref/analysis-final/', '*structvars.csv', -1);

	blacklist_lines = [];
	for file_for_blacklist in files_for_blacklist:
		fp = open(file_for_blacklist, 'r');
		lines = fp.readlines();
		fp.close();
		# if (len(blacklist_lines) == 0):
		# 	blacklist_lines.append(lines[0]);
		for line in lines[1:]:
			line = line.strip();
			line = line.split('#')[0];
			if (len(line) == 0):
				continue;
			split_line = line.split('\t');
			offset = find_event_offset(sorted_real_events, int(split_line[1]));
			split_line[1] = str(int(split_line[1]) + offset);
			split_line[2] = str(int(split_line[2]) + offset);

			split_line[-1] = 'blacklist';
			split_line.append('# %s' % (os.path.basename(file_for_blacklist).split('-')[0]));
			line = '\t'.join(split_line);
			blacklist_lines.append(line);

	### Create a compiled CSV table of the truth and blacklisted events.
	full_truth_events_path = 'data/out/sv-indel_ref/truth_events/truth_events_and_blacklist.csv';
	if (not os.path.exists(os.path.dirname(full_truth_events_path))):
		sys.stderr.write('Creating folder: "%s".\n' % (os.path.dirname(full_truth_events_path)));
		os.makedirs(os.path.dirname(full_truth_events_path));
	sys.stderr.write('Writing the truth events with blacklisted events to file: "%s".\n' % (full_truth_events_path));
	fp = open(full_truth_events_path, 'w');
	fp.write('\n'.join(truth_lines) + '\n');
	fp.write('\n'.join(blacklist_lines) + '\n');
	fp.close();

	files_for_eval = find_files('data/out/sv-indel_ref/analysis-final/', '*chained.csv', -1);
	# files_for_eval = find_files('data/out/sv-indel_ref/analysis-final/', '*structvars.csv', -1);
	command = '%s/evaluate_events.py %s %s' % (SCRIPT_PATH, full_truth_events_path, ' '.join(files_for_eval));
	execute_command(command);

def find_event_offset(sorted_real_events, position):
	i = 0;
	while (i < len(sorted_real_events)):
	# for event in sorted_real_events:
		event = sorted_real_events[i];
		if (event[0] > position):
			if (i == 0):
				return 0;
			else:
				print 'Found event for position %d:' % (position);
				print event;
				return (sorted_real_events[i-1][-1]);
		i += 1;
	return 0;





def calc_erroneous_windows(sam_file, window_length, min_error_rate):
	fp_in = None;
	
	try:
		fp_in = open(sam_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, sam_file));
		return [-1.0, 0, 0];
	
	num_windows = 0;
	num_erroneous_windows = 0;
	all_window_err_rates = [];

	i = 0;
	for line in fp_in:
		i += 1;
		if (len(line.strip()) == 0 or line[0] == '@'):
			continue;
		# sys.stderr.write('\rLine %d, num_accepted: %d, num_rejected: %d' % (i, num_accepted, num_rejected));

		# if (i < 100):
		# 	continue;
		# if (i > 200):
		# 	break;
		
		sam_line = utility_sam.SAMLine(line.rstrip());
		if (sam_line.IsMapped() == False):
			continue;
		if (len(sam_line.seq) < window_length):
			continue;

		# num_windows += len(sam_line.seq) - window_length + 1;
		[current_ratio, current_num_over_threshold, current_num_windows] = sam_line.CountErroneousWindows(window_length, min_error_rate);
		num_erroneous_windows += current_num_over_threshold;
		num_windows += current_num_windows;
		if (current_ratio >= 0.0):
			all_window_err_rates.append(current_ratio);

		rate = (-1.0) if (num_windows == 0) else (float(num_erroneous_windows) / float(num_windows));
		sys.stdout.write('\r[%d] Window rate: %.2f' % (i, rate));

	fp_in.close();

	sys.stdout.write('\n');
	# mean_window_err_rate = (-1.0) if (num_windows == 0) else (float(num_erroneous_windows) / float(num_windows));
	mean_window_err_rate = np.mean(all_window_err_rates);
	std_window_err_rate = np.std(all_window_err_rates);
	min_window_err_rate = np.min(all_window_err_rates);
	max_window_err_rate = np.max(all_window_err_rates);

	# print all_window_err_rates;
	# print rate;
	return [mean_window_err_rate, std_window_err_rate, min_window_err_rate, max_window_err_rate, num_erroneous_windows, num_windows];

### Takes two SAM files and compares the aligned qnames. It takes the qnames mapped by one that are not present in the other, and analyzes their error rates.
### sam_file1_path is the SAM file from which the alignments that are not mapped in sam_file2_path will be taken and analyzed.
def RUN_VENN_COMPARISON_TEST(sam_file1_path, sam_file2_path, reference_path, out_venn_prefix):
	sam_file1_path = os.path.abspath(sam_file1_path);
	sam_file2_path = os.path.abspath(sam_file2_path);
	reference_path = os.path.abspath(reference_path);
	out_venn_prefix = os.path.abspath(out_venn_prefix);

	sam_file1_basename = os.path.splitext(os.path.basename(sam_file1_path))[0];
	sam_file2_basename = os.path.splitext(os.path.basename(sam_file2_path))[0];

	# ### Find alignments that GraphMap mapped and LAST didn't, and vice-versa.
	execute_command('%s/venn_of_sams.py %s %s %s' % (SAMSCRIPTS, sam_file1_path, sam_file2_path, out_venn_prefix));

	# # ### Analyze the error rates in those alignments that were mapped only by the first mapper.
	alignments_only_in_first = '%s_qnames_only_in_%s.sam' % (out_venn_prefix, sam_file1_basename);
	execute_command('%s/errorrates.py event %s %s' % (SAMSCRIPTS, reference_path, alignments_only_in_first));

	# # ### Analyze the error rates in those alignments that were mapped only by the second mapper.
	alignments_only_in_second = '%s_qnames_only_in_%s.sam' % (out_venn_prefix, sam_file2_basename);
	execute_command('%s/errorrates.py event %s %s' % (SAMSCRIPTS, reference_path, alignments_only_in_second));

	# # ### Analyze the error rates in those alignments that were mapped by both, but aligned by the first mapper.
	injunct_alignments_path1 = '%s_qnames_in_both-alignments_from_%s.sam' % (out_venn_prefix, sam_file1_basename);
	execute_command('%s/errorrates.py event %s %s' % (SAMSCRIPTS, reference_path, injunct_alignments_path1));

	# # ### Analyze the error rates in those alignments that were mapped by both, but aligned by the first mapper.
	injunct_alignments_path2 = '%s_qnames_in_both-alignments_from_%s.sam' % (out_venn_prefix, sam_file2_basename);
	execute_command('%s/errorrates.py event %s %s' % (SAMSCRIPTS, reference_path, injunct_alignments_path2));


	# fp_out = sys.stdout;

	try:
		fp_out = open('%s_analysis.summary.txt' % (out_venn_prefix), 'w');
	except Exception, e:
		sys.stderr.write('Could not open file: "%s_analysis.summary.txt" for writing!\n' % (out_venn_prefix));
		sys.stderr.write('Using STDOUT instead.\n');
		fp_out = sys.stdout;
		sys.stderr.write('\n');

	files_for_window_count = [
								alignments_only_in_first,
								alignments_only_in_second,
								injunct_alignments_path1,
								injunct_alignments_path2
							];
	for file_path in files_for_window_count:
		lines = '';
		lines += 'Analysis of file: "%s"\n' % (file_path);
		fp_out.write(lines);
		fp_out.flush();
		[mean_window_err_rate, std_window_err_rate, min_window_err_rate, max_window_err_rate, num_err_windows, num_windows] = calc_erroneous_windows(file_path, 100, 0.30);
		lines = '';
		lines += 'Mean rate of erroneous windows: %.2f\n' % (mean_window_err_rate);
		lines += 'STD of rates of erroneous windows: %.2f\n' % (std_window_err_rate);
		lines += 'Min rate of erroneous windows: %.2f\n' % (min_window_err_rate);
		lines += 'Max rate of erroneous windows: %.2f\n' % (max_window_err_rate);
		lines += 'Number of erroneous windows: %d\n' % (num_err_windows);
		lines += 'Number of windows: %d\n' % (num_windows);
		lines += '\n';
		lines += '\n';
		fp_out.write(lines);
		fp_out.flush();
		if (fp_out != sys.stdout):
			sys.stdout.write(lines);
			sys.stdout.flush();

	if (fp_out != sys.stdout):
		fp_out.close();

	# injunct_alignments_path1 = '%s_qnames_in_both-from_sam1.sam' % (out_venn_prefix);
	# fp_out.write('Analysis of file: "%s"\n' % (injunct_alignments_path1));
	# [mean_window_err_rate, std_window_err_rate, min_window_err_rate, max_window_err_rate, num_err_windows, num_windows] = calc_erroneous_windows(injunct_alignments_path1, 100, 0.30);
	# fp_out.write('Mean rate of erroneous windows: %.2f\n' % (mean_window_err_rate));
	# fp_out.write('STD of rates of erroneous windows: %.2f\n' % (std_window_err_rate));
	# fp_out.write('Min rate of erroneous windows: %.2f\n' % (min_window_err_rate));
	# fp_out.write('Max rate of erroneous windows: %.2f\n' % (max_window_err_rate));
	# fp_out.write('Number of erroneous windows: %d\n' % (num_err_windows));
	# fp_out.write('Number of windows: %d\n' % (num_windows));
	# fp_out.write('\n');
	# fp_out.write('\n');
	# fp_out.flush();

	# injunct_alignments_path2 = '%s_qnames_in_both-from_sam2.sam' % (out_venn_prefix);
	# fp_out.write('Analysis of file: "%s"\n' % (injunct_alignments_path2));
	# [mean_window_err_rate, std_window_err_rate, min_window_err_rate, max_window_err_rate, num_err_windows, num_windows] = calc_erroneous_windows(injunct_alignments_path2, 100, 0.30);
	# fp_out.write('Mean rate of erroneous windows: %.2f\n' % (mean_window_err_rate));
	# fp_out.write('STD of rates of erroneous windows: %.2f\n' % (std_window_err_rate));
	# fp_out.write('Min rate of erroneous windows: %.2f\n' % (min_window_err_rate));
	# fp_out.write('Max rate of erroneous windows: %.2f\n' % (max_window_err_rate));
	# fp_out.write('Number of erroneous windows: %d\n' % (num_err_windows));
	# fp_out.write('Number of windows: %d\n' % (num_windows));
	# fp_out.write('\n');
	# fp_out.write('\n');
	# fp_out.flush();

	# if (fp_out != sys.stdout):
	# 	fp_out.close();




### This function takes an input SAM file, filters only unique alignment, filters only alignments withing the amplicon regions, and then runs variant calling. In case aligner was marginAlign, marginCaller is used.
def RUN_AMPLICON_TEST():
	reference = '%s/../data/amplicons-f1000/reference/ref_chr6_chr22-hg19_v38.fa' % (SCRIPT_PATH);
	reads = '%s/../data/amplicons-f1000/reads/reads_2d.fastq' % (SCRIPT_PATH);
	dataset_name = 'amplicons-f1000-2d';
	out_path = '%s/../data/out/%s/' % (SCRIPT_PATH, dataset_name);
	### Map all the amplicon reads to the chr6 and chr22 references.

	REGION_CYP2D6 = ['gi|224589814|ref|NC_000022.10|:42522077-42527144', 'CYP2D6'];
	REGION_HLAA = ['gi|224589818|ref|NC_000006.11|:29909854-29913805', 'HLA-A'];
	REGION_HLAB = ['gi|224589818|ref|NC_000006.11|:31321279-31325303', 'HLA-B'];

	regions = [REGION_CYP2D6, REGION_HLAA, REGION_HLAB];

#	dryrun = False;
	dryrun = True;

	sam_files = [

					'%s/GraphMap-%s.sam' % (out_path, dataset_name),
					'%s/GraphMap-anchor-%s.sam' % (out_path, dataset_name),
					'%s/LAST-%s.sam' % (out_path, dataset_name),
					'%s/BWAMEM-%s.sam' % (out_path, dataset_name),
					'%s/BLASR-%s.sam' % (out_path, dataset_name),
					'%s/marginAlign-%s.sam' % (out_path, dataset_name),
					'%s/marginAlignGraphMap-%s.sam' % (out_path, dataset_name),
					'%s/marginAlignGraphMap-anchor-%s.sam' % (out_path, dataset_name),
					'%s/DALIGNER-%s.sam' % (out_path, dataset_name),
				];

	fp_out = sys.stdout;
	fp_out_path = '%s/collected_variants.csv' % out_path;
	try:
		fp_out = open(fp_out_path, 'w');
		sys.stderr.write('Collecting variant calling results to path: "%s".\n' % (fp_out_path));
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for writing! Outputting to STDOUT instead.\n');

	fp_out.write('Mapper\tRegion\tTP\tFP\tTotal\tTruth\tTruth-PASS\tTruth-NonPASS\tFP-PASS\tFP-NonPASS\tFN-PASS\tFN-NonPASS\n');
	empty_line = '\t'*12;

	### Process all given SAM files.
	for sam_path in sam_files:
		sam_basename = os.path.basename(os.path.splitext(sam_path)[0]);
		sam_out_folder = '%s/inregion-%s/' % (out_path, sam_basename);
		current_region = 0;
		while (current_region < len(regions)):
		# for region in regions:
			region = regions[current_region];
			region_name = region[1];
			sys.stderr.write('Running region %s:' % (region[1]));

			reference_file_for_filtering = None;
			region_to_use = region;

			### First prepare the alignments for variant calling. This includes filtering the uniquely aligned reads, taking only 2d reads, and taking only reads that fully span the region.
			[bam_all_reads_in_region, bam_1d_reads_in_region, bam_2d_reads_in_region] = filter_spanning_reads(dryrun, region, reads, sam_path, sam_out_folder, reference, leftalign=False);
			sys.stderr.write('Return: "%s".\n' % (str([bam_all_reads_in_region, bam_1d_reads_in_region, bam_2d_reads_in_region])));

			sam_2d_reads_in_region = bam_2d_reads_in_region.replace('-sorted.bam', '.sam');
			sam_2d_basename = os.path.splitext(os.path.basename(sam_2d_reads_in_region))[0];
			out_margincaller_vcf = '%s/%s-margincaller.vcf' % (sam_out_folder, os.path.splitext(os.path.basename(sam_2d_reads_in_region))[0]);
			out_lofreq_vcf = '%s/%s-lofreq.vcf' % (sam_out_folder, os.path.splitext(os.path.basename(sam_2d_reads_in_region))[0]);

			run_marginCaller(dryrun, sam_2d_reads_in_region, reference, out_margincaller_vcf);
			### Evaluate the found variants.
			variant_caller = 'marginCaller';
			vcf_known_mutations_path = '%s/../data/amplicons-f1000/truth_variants/sorted-variants-dbSNP_and_NA12878-%s_amplicon-splitsnps.vcf' % (SCRIPT_PATH, region[1]);
			results = evaluate_vcf(sam_2d_basename, region[1], out_margincaller_vcf, vcf_known_mutations_path);
			fp_out.write('%s\t%s\t%s\t%s\n' % (sam_2d_basename, region_name, '\t'.join([str(value) for value in results]), variant_caller));
			fp_out.flush();

			run_lofreq(dryrun, bam_2d_reads_in_region, region_to_use, reference, out_lofreq_vcf);
			### Evaluate the found variants.
			variant_caller = 'lofreq';
			vcf_known_mutations_path = '%s/../data/amplicons-f1000/truth_variants/sorted-variants-dbSNP_and_NA12878-%s_amplicon-splitsnps.vcf' % (SCRIPT_PATH, region[1]);
			sam_2d_basename = os.path.splitext(os.path.basename(sam_2d_reads_in_region))[0];
			results = evaluate_vcf(sam_2d_basename, region[1], out_lofreq_vcf, vcf_known_mutations_path);
			fp_out.write('%s\t%s\t%s\t%s\n' % (sam_2d_basename, region_name, '\t'.join([str(value) for value in results]), variant_caller));
			fp_out.flush();

			current_region += 1;

		fp_out.write(empty_line + '\n');

	if (fp_out != sys.stdout):
		fp_out.close();
		sys.stderr.write('Collected variant calling results in path: "%s".\n' % (fp_out_path));

def RUN_AMPLICON_TEST_1d2d():
	reference = '%s/../data/amplicons-f1000-1d2d/reference/ref_chr6_chr22-hg19_v38.fa' % (SCRIPT_PATH);
	reads = '%s/../data/amplicons-f1000-1d2d/reads/reads_all-f1000.fastq' % (SCRIPT_PATH);
	dataset_name = 'amplicons-f1000-1d2d';
	out_path = '%s/../data/out/%s/' % (SCRIPT_PATH, dataset_name);
	### Map all the amplicon reads to the chr6 and chr22 references.

	REGION_CYP2D6 = ['gi|224589814|ref|NC_000022.10|:42522077-42527144', 'CYP2D6'];
	REGION_HLAA = ['gi|224589818|ref|NC_000006.11|:29909854-29913805', 'HLA-A'];
	REGION_HLAB = ['gi|224589818|ref|NC_000006.11|:31321279-31325303', 'HLA-B'];

	regions = [REGION_CYP2D6, REGION_HLAA, REGION_HLAB];

#	dryrun = False;
	dryrun = True;

	# sam_path = '%s/marginAlign-nanopore-nospecialchars-with_AS.sam' % (out_path);
	sam_files = [
					'%s/GraphMap-%s.sam' % (out_path, dataset_name),
					'%s/GraphMap-anchor-%s.sam' % (out_path, dataset_name),
					'%s/LAST-%s.sam' % (out_path, dataset_name),
					'%s/BWAMEM-%s.sam' % (out_path, dataset_name),
					'%s/BLASR-%s.sam' % (out_path, dataset_name),
					'%s/marginAlign-%s.sam' % (out_path, dataset_name),
					'%s/marginAlignGraphMap-%s.sam' % (out_path, dataset_name),
					'%s/marginAlignGraphMap-anchor-%s.sam' % (out_path, dataset_name),
					'%s/DALIGNER-%s.sam' % (out_path, dataset_name),
					];

	fp_out = sys.stdout;
	fp_out_path = '%s/collected_variants.csv' % out_path;
	try:
		fp_out = open(fp_out_path, 'w');
		sys.stderr.write('Collecting variant calling results to path: "%s".\n' % (fp_out_path));
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for writing! Outputting to STDOUT instead.\n');

	fp_out.write('Mapper\tRegion\tTP\tFP\tTotal\tTruth\tTruth-PASS\tTruth-NonPASS\tFP-PASS\tFP-NonPASS\tFN-PASS\tFN-NonPASS\n');
	empty_line = '\t'*12;

	### Process all given SAM files.
	for sam_path in sam_files:
		sam_basename = os.path.basename(os.path.splitext(sam_path)[0]);
		sam_out_folder = '%s/inregion-%s/' % (out_path, sam_basename);
		current_region = 0;
		while (current_region < len(regions)):
		# for region in regions:
			region = regions[current_region];
			region_name = region[1];
			sys.stderr.write('Running region %s:' % (region[1]));

			reference_file_for_filtering = None;
			region_to_use = region;

			filter_partial_reads(dryrun, region, reads, sam_path, sam_out_folder, reference, leftalign=False);

			current_region += 1;

		fp_out.write(empty_line + '\n');

	if (fp_out != sys.stdout):
		fp_out.close();
		sys.stderr.write('Collected variant calling results in path: "%s".\n' % (fp_out_path));

### This function tests the difference in number of variants between GraphMap and BWA-MEM when the coverage is the same.
### Since GraphMap produces a higher coverage of spanning reads, alignments were subsampled to the same coverage of BWA-MEM.
### Reads that were mapped by both mappers were used in the subsampled set of GraphMap, while the rest needed to fill the coverage
### was randomly selected.
def RUN_AMPLICON_TEST_GRAPHMAP_VS_BWAMEM():
	reference = '%s/../data/amplicons-f1000/reference/ref_chr6_chr22-hg19_v38.fa' % (SCRIPT_PATH);
	reads = '%s/../data/amplicons-f1000/reads/reads_2d.fastq' % (SCRIPT_PATH);
        dataset_name = 'amplicons-f1000-2d';
        out_path = '%s/../data/out/%s/bwamem_graphmap_comparison' % (SCRIPT_PATH, dataset_name);
	### Map all the amplicon reads to the chr6 and chr22 references.

	REGION_CYP2D6 = ['gi|224589814|ref|NC_000022.10|:42522077-42527144', 'CYP2D6'];
	REGION_HLAA = ['gi|224589818|ref|NC_000006.11|:29909854-29913805', 'HLA-A'];
	REGION_HLAB = ['gi|224589818|ref|NC_000006.11|:31321279-31325303', 'HLA-B'];

	regions = [REGION_CYP2D6, REGION_HLAA, REGION_HLAB];

	dryrun = False;
#	dryrun = True;

	sam_files = [
					'%s/subsampled-GraphMap-CYP2D6.sam' % (out_path),
					'%s/subsampled-GraphMap-HLA-A.sam' % (out_path),
					'%s/subsampled-GraphMap-HLA-B.sam' % (out_path),
					];

	fp_out = sys.stdout;
	fp_out_path = '%s/collected_variants.csv' % out_path;
	try:
		fp_out = open(fp_out_path, 'w');
		sys.stderr.write('Collecting variant calling results to path: "%s".\n' % (fp_out_path));
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for writing! Outputting to STDOUT instead.\n');

	fp_out.write('Mapper\tRegion\tTP\tFP\tTotal\tTruth\tTruth-PASS\tTruth-NonPASS\tFP-PASS\tFP-NonPASS\tFN-PASS\tFN-NonPASS\n');
	empty_line = '\t'*12;

	### Process all given SAM files.
	for sam_path in sam_files:
		sam_basename = os.path.basename(os.path.splitext(sam_path)[0]);
		sam_out_folder = '%s/inregion-%s/' % (out_path, sam_basename);
		current_region = 0;
		while (current_region < len(regions)):
		# for region in regions:
			region = regions[current_region];
			region_name = region[1];
			sys.stderr.write('Running region %s:' % (region[1]));

			reference_file_for_filtering = None;
			region_to_use = region;

			### First prepare the alignments for variant calling. This includes filtering the uniquely aligned reads, taking only 2d reads, and taking only reads that fully span the region.
			[bam_all_reads_in_region, bam_1d_reads_in_region, bam_2d_reads_in_region] = filter_spanning_reads(dryrun, region, reads, sam_path, sam_out_folder, reference, leftalign=False);
			sys.stderr.write('Return: "%s".\n' % (str([bam_all_reads_in_region, bam_1d_reads_in_region, bam_2d_reads_in_region])));

			sam_2d_reads_in_region = bam_2d_reads_in_region.replace('-sorted.bam', '.sam');
			sam_2d_basename = os.path.splitext(os.path.basename(sam_2d_reads_in_region))[0];
			out_lofreq_vcf = '%s/%s-lofreq.vcf' % (sam_out_folder, os.path.splitext(os.path.basename(sam_2d_reads_in_region))[0]);

			run_lofreq(dryrun, bam_2d_reads_in_region, region_to_use, reference, out_lofreq_vcf);
			### Evaluate the found variants.
			variant_caller = 'lofreq';
			vcf_known_mutations_path = '%s/../data/amplicons-f1000/truth_variants/sorted-variants-dbSNP_and_NA12878-%s_amplicon-splitsnps.vcf' % (SCRIPT_PATH, region[1]);
			sam_2d_basename = os.path.splitext(os.path.basename(sam_2d_reads_in_region))[0];
			results = evaluate_vcf(sam_2d_basename, region[1], out_lofreq_vcf, vcf_known_mutations_path);
			fp_out.write('%s\t%s\t%s\t%s\n' % (sam_2d_basename, region_name, '\t'.join([str(value) for value in results]), variant_caller));
			fp_out.flush();

			current_region += 1;

		fp_out.write(empty_line + '\n');

	if (fp_out != sys.stdout):
		fp_out.close();
		sys.stderr.write('Collected variant calling results in path: "%s".\n' % (fp_out_path));




def run_marginCaller(dry_run, sam_file, reference_path, out_vcf):
	if (not os.path.exists(os.path.dirname(out_vcf))):
		sys.stderr.write('Creating output path: "%s".\n' % os.path.dirname(out_vcf));
		os.makedirs(os.path.dirname(out_vcf));

	# fastq_header_hash = 
	# fastqparser.keep_gi_header
	### Adjust the reference file to conform to marginAlign/marginCaller bugs.
	sys.stderr.write('[] Keeping only the gi part of reference headers.\n');
	reference_giheader = '%s/%s-giheader.fa' % (os.path.dirname(sam_file), os.path.splitext(os.path.basename(reference_path))[0]);
	execute_command_w_dryrun(dry_run, '%s/fastqfilter.py giheader %s %s' % (SAMSCRIPTS, reference_path, reference_giheader));
	sys.stderr.write('\n');
	reference_marginalign = '%s/%s-giheader-marginalign.fa' % (os.path.dirname(sam_file), os.path.splitext(os.path.basename(reference_path))[0]);


	sys.stderr.write('[] Filtering FASTQ headers for marginAlign bugs. (reference_giheader: "%s", reference_marginalign: "%s").\n' % (reference_giheader, reference_marginalign));
	ref_header_hash = filter_fastq_for_marginalign(reference_giheader, reference_marginalign, None);

	sam_file_marginalign = '%s/%s-fixed_for_marginalign.sam' % (os.path.dirname(sam_file), os.path.splitext(os.path.basename(sam_file))[0]);
	# [qname_hash, rname_hash] = filter_sam_for_marginalign(sam_file, sam_file_marginalign);
	execute_command_w_dryrun(dry_run, '%s/samfilter.py marginalign %s %s' % (SAMSCRIPTS, sam_file, sam_file_marginalign));

	# sys.stderr.write('Return: "%s".\n' % (str([bam_all_reads_in_region, bam_1d_reads_in_region, bam_2d_reads_in_region])));
	jobtree = 'jobTree';
	if (os.path.exists(jobtree)):
		execute_command('rm -r %s' % (jobtree));

	### Call variants using marginCaller.
	# execute_command('%s/aligneval/aligners/marginAlign/marginCaller %s %s %s --jobTree %s' % (tools_path, sam_2d_reads_in_region, marginAlign_reference_file, out_vcf, jobtree));
	execute_command_w_dryrun(dry_run, '%s/aligneval/aligners/marginAlign/marginCaller %s %s %s.temp --jobTree %s' % (tools_path, sam_file_marginalign, reference_marginalign, out_vcf, jobtree));

	# fix_sam_qnames_after_marginAlign(sam, ref_header_hash, read_header_hash, out_sam_path);
	fix_vcf_rnames_after_margincaller('%s.temp' % (out_vcf), ref_header_hash, out_vcf);
	# execute_command_w_dryrun(dry_run, 'rm %s.temp' % (out_vcf));

def run_lofreq(dry_run, bam_file, region, reference_path, out_vcf):
	if (not os.path.exists(os.path.dirname(out_vcf))):
		sys.stderr.write('Creating output path: "%s".\n' % os.path.dirname(out_vcf));
		os.makedirs(os.path.dirname(out_vcf));

	if (os.path.exists(out_vcf)):
		sys.stderr.write('Output VCF file already exists. Need to move it because LoFreq won\'t overwrite. Moving "%s" to "%s.bak".\n' % (out_vcf, out_vcf));
		os.rename(out_vcf, out_vcf + '.bak');

	a = 0.01;
	q = 0;
	Q = 0;
	execute_command_w_dryrun(False, '%s/lofreq_star-2.1.2/bin/lofreq call -a %f -q %d -Q %d --no-default-filter --verbose -B -f %s -o %s %s -r "%s"' % (tools_path, a, q, Q, reference_path, out_vcf, bam_file, region[0]));

def fix_vcf_rnames_after_margincaller(in_vcf, ref_header_hash, out_vcf):
	try:
		fp_in = open(in_vcf, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Skipping.\n' % in_vcf);
		return;

	if (in_vcf == out_vcf):
		sys.stderr.write('ERROR: Output and input files are the same! Skipping.\n');
		return;
	try:
		fp_out = open(out_vcf, 'w');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for writing! Skipping.\n' % out_vcf);
		return;

	for line in fp_in:
		line = line.strip();
		if (len(line) == 0 or line[0] == '#'):
			fp_out.write(line + '\n');
			continue;
		split_line = line.split();
		rname = split_line[0];
		try:
			new_rname = ref_header_hash[rname];
		except:
			new_rname = rname;
		new_line = line.replace(rname, new_rname);
		fp_out.write(new_line + '\n');

	fp_in.close();
	fp_out.close();




def evaluate_sv_consensus(reference, dataset_name, out_path):
	# ### Missing part: evaluation of the SV consensus!
	# ../samscripts/src/samfilter.py uniquebest data/out/$sam_path/$samname.sam data/out/$sam_path/${samname}-uniquebest.sam
	# ../samscripts/src/svconsensus.py data/reference/escherichia_coli.fa data/out/$sam_path/$samname2.sam data/out/$sam_path/analysis-final/$samname2.structvars.csv
	mapper_names = ['GraphMap', 'GraphMap-anchor', 'LAST', 'BWAMEM', 'BLASR', 'DALIGNER', 'marginAlignGraphMap', 'marginAlign', 'marginAlignGraphMap-anchor'];
	for mapper_name in mapper_names:
		sam = '%s/%s-%s.sam' % (out_path, mapper_name, dataset_name);
		sam_with_AS = '%s/%s-%s-with_AS.sam' % (out_path, mapper_name, dataset_name);
		am_uniquebest = '%s/%s-%s-with_AS-uniquebest.sam' % (out_path, mapper_name, dataset_name);
		out_sv_path = '%s/analysis-final/%s.csv' % (out_path, mapper_name)
		execute_command('%s/samscripts/src/samfilter.py generateAS %s %s %s' % (tools_path, reference, out_sam, out_sam_with_AS));
		execute_command('%s/samscripts/src/samfilter.py uniquebest %s %s' % (tools_path, out_sam_with_AS, out_sam_uniquebest));
		execute_command('%s/svconsensus.py %s %s %s' % (tools_path, reference, out_sam_uniquebest, out_sv_path));

def evaluate_alignments(reference, reads, dataset_name, out_path, cov_threshold=20):
	out_collect_file = '%s/collected.csv' % (out_path);

	out_sam = '%s/DALIGNER-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file hcalc %s %s %s %d >> %s' % (tools_path, out_sam, reference, reads, cov_threshold, out_collect_file));

	out_sam = '%s/GraphMap-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s %d >> %s' % (tools_path, out_sam, reference, reads, cov_threshold, out_collect_file));

	out_sam = '%s/GraphMap-anchor-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s %d >> %s' % (tools_path, out_sam, reference, reads, cov_threshold, out_collect_file));

	out_sam = '%s/LAST-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s %d >> %s' % (tools_path, out_sam, reference, reads, cov_threshold, out_collect_file));

	out_sam = '%s/BWAMEM-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s %d >> %s' % (tools_path, out_sam, reference, reads, cov_threshold, out_collect_file));

	out_sam = '%s/BLASR-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s %d >> %s' % (tools_path, out_sam, reference, reads, cov_threshold, out_collect_file));

	out_sam = '%s/marginAlignGraphMap-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s %d >> %s' % (tools_path, out_sam, reference, reads, cov_threshold, out_collect_file));

	out_sam = '%s/marginAlign-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s %d >> %s' % (tools_path, out_sam, reference, reads, cov_threshold, out_collect_file));

	out_sam = '%s/marginAlignGraphMap-anchor-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s %d >> %s' % (tools_path, out_sam, reference, reads, cov_threshold, out_collect_file));

### Filters only one unique (best) alignment for each read, and then evaluates the results.
def evaluate_unique_alignments(reference, reads, dataset_name, out_path, cov_thresholds=[20]):
	out_collect_file = '%s/collected.csv' % (out_path);

	for cov_threshold in cov_thresholds:
		out_sam = '%s/LAST-%s.sam' % (out_path, dataset_name);
		out_sam_uniquebest = '%s/LAST-%s-uniquebest.sam' % (out_path, dataset_name);
		execute_command('%s/samscripts/src/samfilter.py uniquebest %s %s' % (tools_path, out_sam, out_sam_uniquebest));
		execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s %d >> %s' % (tools_path, out_sam_uniquebest, reference, reads, cov_threshold, out_collect_file));

		out_sam = '%s/DALIGNER-%s.sam' % (out_path, dataset_name);
		out_sam_uniquebest = '%s/DALIGNER-%s-uniquebest.sam' % (out_path, dataset_name);
		execute_command('%s/samscripts/src/samfilter.py uniquebest %s %s' % (tools_path, out_sam, out_sam_uniquebest));
		execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s %d >> %s' % (tools_path, out_sam_uniquebest, reference, reads, cov_threshold, out_collect_file));

		out_sam = '%s/marginAlign-%s.sam' % (out_path, dataset_name);
		out_sam_with_AS = '%s/marginAlign-%s-with_AS.sam' % (out_path, dataset_name);
		out_sam_uniquebest = '%s/marginAlign-%s-with_AS-uniquebest.sam' % (out_path, dataset_name);
		execute_command('%s/samscripts/src/samfilter.py generateAS %s %s %s' % (tools_path, reference, out_sam, out_sam_with_AS));
		execute_command('%s/samscripts/src/samfilter.py uniquebest %s %s' % (tools_path, out_sam_with_AS, out_sam_uniquebest));
		execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s %d >> %s' % (tools_path, out_sam_uniquebest, reference, reads, cov_threshold, out_collect_file));

		out_sam = '%s/marginAlignGraphMap-%s.sam' % (out_path, dataset_name);
		out_sam_with_AS = '%s/marginAlignGraphMap-%s-with_AS.sam' % (out_path, dataset_name);
		out_sam_uniquebest = '%s/marginAlignGraphMap-%s-with_AS-uniquebest.sam' % (out_path, dataset_name);
		execute_command('%s/samscripts/src/samfilter.py generateAS %s %s %s' % (tools_path, reference, out_sam, out_sam_with_AS));
		execute_command('%s/samscripts/src/samfilter.py uniquebest %s %s' % (tools_path, out_sam_with_AS, out_sam_uniquebest));
		execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s %d >> %s' % (tools_path, out_sam_uniquebest, reference, reads, cov_threshold, out_collect_file));

		out_sam = '%s/marginAlignGraphMap-anchor-%s.sam' % (out_path, dataset_name);
		out_sam_with_AS = '%s/marginAlignGraphMap-anchor-%s-with_AS.sam' % (out_path, dataset_name);
		out_sam_uniquebest = '%s/marginAlignGraphMap-anchor-%s-with_AS-uniquebest.sam' % (out_path, dataset_name);
		execute_command('%s/samscripts/src/samfilter.py generateAS %s %s %s' % (tools_path, reference, out_sam, out_sam_with_AS));
		execute_command('%s/samscripts/src/samfilter.py uniquebest %s %s' % (tools_path, out_sam_with_AS, out_sam_uniquebest));
		execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s %d >> %s' % (tools_path, out_sam_uniquebest, reference, reads, cov_threshold, out_collect_file));

def collect_alignments(reference, reads, dataset_name, out_path, cov_thresholds=[20]):
	out_collect_file = '%s/collected.csv' % (out_path);
	
	is_first_line = True;
	for cov_threshold in cov_thresholds:
		out_sam = '%s/DALIGNER-%s.sam' % (out_path, dataset_name);
		if (is_first_line == True):
			execute_command('%s/samscripts/src/alignmentstats.py file hcollect %s %s %s %d >> %s' % (tools_path, out_sam, reference, reads, cov_threshold, out_collect_file));
		else:
			execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s %d >> %s' % (tools_path, out_sam, reference, reads, cov_threshold, out_collect_file));

		out_sam = '%s/GraphMap-%s.sam' % (out_path, dataset_name);
		execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s %d >> %s' % (tools_path, out_sam, reference, reads, cov_threshold, out_collect_file));

		out_sam = '%s/GraphMap-anchor-%s.sam' % (out_path, dataset_name);
		execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s %d >> %s' % (tools_path, out_sam, reference, reads, cov_threshold, out_collect_file));

		out_sam = '%s/LAST-%s.sam' % (out_path, dataset_name);
		execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s %d >> %s' % (tools_path, out_sam, reference, reads, cov_threshold, out_collect_file));

		out_sam = '%s/BWAMEM-%s.sam' % (out_path, dataset_name);
		execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s %d >> %s' % (tools_path, out_sam, reference, reads, cov_threshold, out_collect_file));

		out_sam = '%s/BLASR-%s.sam' % (out_path, dataset_name);
		execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s %d >> %s' % (tools_path, out_sam, reference, reads, cov_threshold, out_collect_file));

		out_sam = '%s/marginAlignGraphMap-%s.sam' % (out_path, dataset_name);
		execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s %d >> %s' % (tools_path, out_sam, reference, reads, cov_threshold, out_collect_file));

		out_sam = '%s/marginAlign-%s.sam' % (out_path, dataset_name);
		execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s %d >> %s' % (tools_path, out_sam, reference, reads, cov_threshold, out_collect_file));

		out_sam = '%s/marginAlignGraphMap-anchor-%s.sam' % (out_path, dataset_name);
		execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s %d >> %s' % (tools_path, out_sam, reference, reads, cov_threshold, out_collect_file));
		
		is_first_line = False;

def collect_unique_alignments(reference, reads, dataset_name, out_path):
	out_collect_file = '%s/collected.csv' % (out_path);

	out_sam = '%s/LAST-%s-uniquebest.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file hcollect %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

	out_sam = '%s/DALIGNER-%s-uniquebest.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

	out_sam = '%s/marginAlign-%s-with_AS-uniquebest.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

### Params:
### 	orig_reference is the reference to which to map to.
###		orig_reads is the file containing the reads to map to.
###		dataset_name is the name that will be added to the suffix of the output SAM file.
###		out_path is the path to the folder
def run_all_mappers_only(orig_reference, orig_reads, dataset_name, out_path, machine_name, do_not_recalc=True, is_circular=True, select_mappers=['daligner', 'graphmap', 'graphmap-anchor', 'last', 'bwamem', 'blasr', 'marginalign', 'marginaligngraphmap', 'marginaligngraphmap-anchor']):
	if (not os.path.exists(out_path)):
		sys.stderr.write('Creating output path: "%s".\n' % out_path);
		execute_command('mkdir -p %s' % out_path);
	num_threads = multiprocessing.cpu_count() / 2;
	reads_basename = os.path.basename(os.path.splitext(orig_reads)[0]);
	out_collect_file = '%s/collected.csv' % (out_path);

	# dataset_name = os.path.splitext(os.path.basename(orig_reads))[0];
	sys.stderr.write('Dataset name: "%s".\n' % (dataset_name));

	for mapper in select_mappers:
		if (mapper == 'graphmap'):
			out_sam = '%s/GraphMap-%s.sam' % (out_path, dataset_name);
			if (do_not_recalc == False or (not os.path.exists(out_sam))):
				execute_command('%s/aligneval/wrappers/wrapper_graphmap.py run %s %s nanopore%s %s %s' % (tools_path, orig_reads, orig_reference, ('circ' if (is_circular == True) else ''), out_path, dataset_name));
			else:
				sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));

		elif (mapper == 'graphmap-anchor'):
			out_sam = '%s/GraphMap-anchor-%s.sam' % (out_path, dataset_name);
			if (do_not_recalc == False or (not os.path.exists(out_sam))):
				execute_command('%s/aligneval/wrappers/wrapper_graphmap.py run %s %s anchor%s %s anchor-%s' % (tools_path, orig_reads, orig_reference, ('circ' if (is_circular == True) else ''), out_path, dataset_name));
			else:
				sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));

		elif (mapper == 'last'):
			out_sam = '%s/LAST-%s.sam' % (out_path, dataset_name);
			if (do_not_recalc == False or (not os.path.exists(out_sam))):
				execute_command('%s/aligneval/wrappers/wrapper_lastal.py run %s %s %s %s %s' % (tools_path, orig_reads, orig_reference, machine_name, out_path, dataset_name));
			else:
				sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));

		elif (mapper == 'bwamem'):
			out_sam = '%s/BWAMEM-%s.sam' % (out_path, dataset_name);
			if (do_not_recalc == False or (not os.path.exists(out_sam))):
				execute_command('%s/aligneval/wrappers/wrapper_bwamem.py run %s %s %s %s %s' % (tools_path, orig_reads, orig_reference, machine_name, out_path, dataset_name));
			else:
				sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));

		elif (mapper == 'blasr'):
			out_sam = '%s/BLASR-%s.sam' % (out_path, dataset_name);
			if (do_not_recalc == False or (not os.path.exists(out_sam))):
				execute_command('%s/aligneval/wrappers/wrapper_blasr.py run %s %s %s %s %s' % (tools_path, orig_reads, orig_reference, machine_name, out_path, dataset_name));
			else:
				sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));

		elif (mapper == 'daligner'):
			out_sam = '%s/DALIGNER-%s.sam' % (out_path, dataset_name);
			if (do_not_recalc == False or (not os.path.exists(out_sam))):
				execute_command('%s/aligneval/wrappers/wrapper_daligner.py run %s %s %s %s %s' % (tools_path, orig_reads, orig_reference, machine_name, out_path, dataset_name));
			else:
				sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));

		elif (mapper == 'marginalign'):
			out_sam = '%s/marginAlign-%s.sam' % (out_path, dataset_name);
			if (do_not_recalc == False or (not os.path.exists(out_sam))):
				execute_command('%s/aligneval/wrappers/wrapper_marginalign.py run %s %s %s %s %s' % (tools_path, orig_reads, orig_reference, machine_name, out_path, dataset_name));
			else:
				sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));

		elif (mapper == 'marginaligngraphmap'):
			out_sam = '%s/marginAlignGraphMap-%s.sam' % (out_path, dataset_name);
			if (do_not_recalc == False or (not os.path.exists(out_sam))):
				execute_command('%s/aligneval/wrappers/wrapper_marginaligngraphmap.py run %s %s %s %s %s' % (tools_path, orig_reads, orig_reference, machine_name, out_path, dataset_name));
			else:
				sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));

		elif (mapper == 'marginaligngraphmap-anchor'):
			out_sam = '%s/marginAlignGraphMap-anchor-%s.sam' % (out_path, dataset_name);
			if (do_not_recalc == False or (not os.path.exists(out_sam))):
				execute_command('%s/aligneval/wrappers/wrapper_marginaligngraphmap.py run %s %s anchor %s %s' % (tools_path, orig_reads, orig_reference, out_path, 'anchor-' + dataset_name));
			else:
				sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));






### Mutates the given reference, and writes the mutations in a vcf file.
def generate_mutated_reference(reference_path, snp_rate, indel_rate, out_path):
	reference_path = os.path.abspath(reference_path);
	out_path = os.path.abspath(out_path);
	if (not os.path.exists(out_path)):
		os.makedirs(out_path);

	out_prefix = '%s/mutated_%s_snp%f_indel%f' % (out_path, os.path.splitext(os.path.basename(reference_path))[0], snp_rate, indel_rate);
	out_vcf = os.path.abspath('%s.vcf' % (out_prefix));
	out_rev_vcf = '%s/rev_%s.vcf' % (out_path, os.path.basename(out_prefix));
	ref_ext = os.path.splitext(reference_path)[-1];
	out_ref_file = '%s%s' % (out_prefix, ref_ext);

	sys.stderr.write('Mutating the reference using Mutatrix, output VCF file: "%s".\n' % (out_vcf));
	if (indel_rate != 0):
		execute_command('cd %s; %s/mutatrix/mutatrix --snp-rate %f --population-size 1 --microsat-min-len 0 --mnp-ratio 0 --indel-rate %f --indel-max 10 %s > %s' % (out_path, tools_path, snp_rate, indel_rate, reference_path, out_vcf));
	else:
		execute_command('cd %s; %s/mutatrix/mutatrix --snp-rate %f --population-size 1 --microsat-min-len 0 --mnp-ratio 0 --indel-rate 0 --indel-max 0 %s > %s' % (out_path, tools_path, snp_rate, reference_path, out_vcf));

	sys.stderr.write('Reversing the SNP bases in the VCF file, output VCF file: "%s".\n' % (out_rev_vcf));
	execute_command(r"cat %s | awk -F '\t' 'BEGIN {OFS = FS} {if ($0 == /^#.*/) print ; else {a=$4; $4=$5; $5=a; print } }' > %s" % (out_vcf, out_rev_vcf));

	sys.stderr.write('Compressing and indexing the VCF file.\n');
	execute_command('bgzip -c %s > %s.gz' % (out_rev_vcf, out_rev_vcf));
	execute_command('tabix -p vcf %s.gz' % (out_rev_vcf));

	### Mutatrix splits all reference sequences into separate files. This part of code joins them back into one FASTA file.
	[headers, lengths] = fastqparser.get_headers_and_lengths(reference_path);
	print headers;
	all_files = ['%s/1:%s:0%s' % (out_path, header.split(' ')[0], ref_ext) for header in headers];
	if (os.path.exists(out_ref_file)):
		os.rename(out_ref_file, '%s.bak' % (out_ref_file));
	for ref_file in all_files:
		### Take care of the special characters.
		escaped_ref_file = ref_file.replace('|', '\|');
		execute_command('cat %s >> %s' % (escaped_ref_file, out_ref_file));
		if (len(ref_file) > 0 and ('*' in ref_file) == False):
			print 'Removing file: "%s".' % (ref_file);
			os.remove(ref_file);

def make_consensus_reference_from_vcf(reference_file, vcf_file, out_consensus_sequence_file):
	if (not os.path.exists(os.path.dirname(out_consensus_sequence_file))):
		sys.stderr.write('Creating a folder on path: "%s".\n' % (os.path.dirname(out_consensus_sequence_file)));
		os.makedirs(os.path.dirname(out_consensus_sequence_file));

	sys.stderr.write('%s\n' % (reference_file));
	sys.stderr.write('%s\n' % (vcf_file));
	sys.stderr.write('%s\n' % (out_consensus_sequence_file));
	sys.stderr.write('\n');

	fp = open(vcf_file, 'r');
	fp_temp = open('.temp.vcf', 'w');
	for line in fp:
		line = line.strip();
		if (len(line) == 0 or (len(line) > 0 and line[0] == '#')):
			fp_temp.write(line + '\n');
			continue;
		split_line = line.split('\t');
		if (split_line[3] == 'N' and split_line[4] == 'N'):
			continue;
		fp_temp.write(line + '\n');
	fp.close();
	fp_temp.close();
	os.rename(vcf_file, '%s.bak' % (vcf_file));
	os.rename('.temp.vcf', vcf_file);

	sys.stderr.write('Making a Picard dictionary of the reference.\n');
	execute_command('java -jar %s/picard-tools-1.138/picard.jar CreateSequenceDictionary R= %s O= %s.dict' % (tools_path, reference_file, os.path.splitext(reference_file)[0]));
	execute_command('samtools faidx %s' % (reference_file));

	sys.stderr.write('Applying the VCF file to the reference to generate the alternate (consensus) sequence.\n');
	execute_command('java -jar %s/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R %s -o %s -V %s' % (tools_path, reference_file, out_consensus_sequence_file, vcf_file));

def run_dnadiff(reference_file, consensus_sequence_file, out_prefix):
	if (not os.path.exists(os.path.dirname(out_prefix))):
		sys.stderr.write('Creating a folder on path: "%s".\n' % (os.path.dirname(out_prefix)));
		os.makedirs(os.path.dirname(out_prefix));

	sys.stderr.write('Comparing the consensus sequence to the reference using MUMmer.\n');
	execute_command('dnadiff %s %s -p %s' % (reference_file, consensus_sequence_file, out_prefix));

def collect_consensus_sequence_comparison_results(dnadiff_path, out_collect_file):
	if (dnadiff_path.endswith('/') == False):
		dnadiff_path += '/';

	# print dnadiff_path;
	# print out_collect_file;
	sys.stderr.write('Collecting results in file: "%s"\n' % (out_collect_file));

	fp_out = open(out_collect_file, 'w');
	fp_out.write('Mapper\tIndels\tSNPs\tUncalled bases\n');

	report_files = find_files(dnadiff_path, '*.report');
	for report_file in report_files:
		mapper_name = report_file.split(dnadiff_path)[-1].split('/')[0];
		if ('BLASR' in mapper_name):
			mapper_name = 'BLASR';
		elif (('BWAMEM' in mapper_name)):
			mapper_name = 'BWA-MEM';
		elif (('ref_vs_draft' in mapper_name)):
			mapper_name = 'Original ref vs. Mutated ref';
		elif (('DALIGNER' in mapper_name) and ('uniquebest' in mapper_name)):
			mapper_name = 'DALIGNER unique best';
		elif (('DALIGNER' in mapper_name) and ('uniquebest' in mapper_name) == False):
			mapper_name = 'DALIGNER';

		elif (('LAST' in mapper_name) and ('marginAlign' in mapper_name) == False and ('uniquebest' in mapper_name)):
			mapper_name = 'LAST unique best';
		elif (('LAST' in mapper_name) and ('marginAlign' in mapper_name) == False and ('uniquebest' in mapper_name) == False):
			mapper_name = 'LAST';

		elif (('GraphMap' in mapper_name) and ('marginAlign' in mapper_name) == False and ('anchor' in mapper_name)):
			mapper_name = 'GraphMap anchored';
		elif (('GraphMap' in mapper_name) and ('marginAlign' in mapper_name) == False and ('anchor' in mapper_name) == False):
			mapper_name = 'GraphMap';

		elif (('marginAlign' in mapper_name) == True and (('LAST' in mapper_name) == True or (('LAST' in mapper_name) == False and ('GraphMap' in mapper_name) == False)) and ('uniquebest' in mapper_name)):
			mapper_name = 'marginAlign (LAST unique best)';
		elif (('marginAlign' in mapper_name) == True and (('LAST' in mapper_name) == True or (('LAST' in mapper_name) == False and ('GraphMap' in mapper_name) == False)) and ('uniquebest' in mapper_name) == False):
			mapper_name = 'marginAlign (LAST)';

		elif (('GraphMap' in mapper_name) and ('marginAlign' in mapper_name) == True and ('anchor' in mapper_name)):
			mapper_name = 'marginAlign (GraphMap anchored)';
		elif (('GraphMap' in mapper_name) and ('marginAlign' in mapper_name) == True and ('anchor' in mapper_name) == False):
			mapper_name = 'marginAlign (GraphMap)';

		fp = open(report_file, 'r');
		lines = fp.readlines();
		fp.close();
		snps = 0;
		indels = 0;
		insertions = 0;
		deletions = 0;
		unaligned_bases = 0;
		total_bases = 0;
		for line in lines:
			if (line.startswith('TotalSNPs')):
				snps = line.strip().split()[-1];
			elif (line.startswith('TotalIndels')):
				indels = line.strip().split()[-1];
			elif (line.startswith('UnalignedBases')):
				unaligned_bases = line.strip().split()[-1].split('(')[0];
			elif (line.startswith('TotalBases')):
				total_bases = line.strip().split()[-1];
		out_line = '%s\t%s\t%s\t%s' % (mapper_name, indels, snps, unaligned_bases);
		# sys.stdout.write(out_line + '\n');
		fp_out.write(out_line + '\n');
		# print mapper_name;
		# print report_file;
	fp_out.close();

def evaluate_consensus_sequences(reference_file, assembly_sequence_file, dataset_name, alignments_path, cov_threshold=20):
	dnadiff_prefix = '%s/dnadiff_cov_%d/%s/%s' % (alignments_path, cov_threshold, 'ref_vs_draft', 'ref_vs_draft');
	run_dnadiff(reference_file, assembly_sequence_file, dnadiff_prefix);

	vcf_file = '%s/analysis-intermediate/consensus-DALIGNER-%s-cov_%d.variant.vcf' % (alignments_path, dataset_name, cov_threshold);
	consensus_sequence_file = '%s/consensus_sequence/%s.fasta' % (alignments_path, os.path.basename(os.path.splitext(vcf_file)[0]));
	dnadiff_prefix = '%s/dnadiff_cov_%d/%s/%s' % (alignments_path, cov_threshold, os.path.basename(os.path.splitext(vcf_file)[0]), os.path.basename(os.path.splitext(vcf_file)[0]));
	make_consensus_reference_from_vcf(assembly_sequence_file, vcf_file, consensus_sequence_file);
	run_dnadiff(reference_file, consensus_sequence_file, dnadiff_prefix);

	vcf_file = '%s/analysis-intermediate/consensus-GraphMap-%s-cov_%d.variant.vcf' % (alignments_path, dataset_name, cov_threshold);
	consensus_sequence_file = '%s/consensus_sequence/%s.fasta' % (alignments_path, os.path.basename(os.path.splitext(vcf_file)[0]));
	dnadiff_prefix = '%s/dnadiff_cov_%d/%s/%s' % (alignments_path, cov_threshold, os.path.basename(os.path.splitext(vcf_file)[0]), os.path.basename(os.path.splitext(vcf_file)[0]));
	make_consensus_reference_from_vcf(assembly_sequence_file, vcf_file, consensus_sequence_file);
	run_dnadiff(reference_file, consensus_sequence_file, dnadiff_prefix);

	vcf_file = '%s/analysis-intermediate/consensus-GraphMap-anchor-%s-cov_%d.variant.vcf' % (alignments_path, dataset_name, cov_threshold);
	consensus_sequence_file = '%s/consensus_sequence/%s.fasta' % (alignments_path, os.path.basename(os.path.splitext(vcf_file)[0]));
	dnadiff_prefix = '%s/dnadiff_cov_%d/%s/%s' % (alignments_path, cov_threshold, os.path.basename(os.path.splitext(vcf_file)[0]), os.path.basename(os.path.splitext(vcf_file)[0]));
	make_consensus_reference_from_vcf(assembly_sequence_file, vcf_file, consensus_sequence_file);
	run_dnadiff(reference_file, consensus_sequence_file, dnadiff_prefix);

	vcf_file = '%s/analysis-intermediate/consensus-LAST-%s-cov_%d.variant.vcf' % (alignments_path, dataset_name, cov_threshold);
	consensus_sequence_file = '%s/consensus_sequence/%s.fasta' % (alignments_path, os.path.basename(os.path.splitext(vcf_file)[0]));
	dnadiff_prefix = '%s/dnadiff_cov_%d/%s/%s' % (alignments_path, cov_threshold, os.path.basename(os.path.splitext(vcf_file)[0]), os.path.basename(os.path.splitext(vcf_file)[0]));
	make_consensus_reference_from_vcf(assembly_sequence_file, vcf_file, consensus_sequence_file);
	run_dnadiff(reference_file, consensus_sequence_file, dnadiff_prefix);

	vcf_file = '%s/analysis-intermediate/consensus-BWAMEM-%s-cov_%d.variant.vcf' % (alignments_path, dataset_name, cov_threshold);
	consensus_sequence_file = '%s/consensus_sequence/%s.fasta' % (alignments_path, os.path.basename(os.path.splitext(vcf_file)[0]));
	dnadiff_prefix = '%s/dnadiff_cov_%d/%s/%s' % (alignments_path, cov_threshold, os.path.basename(os.path.splitext(vcf_file)[0]), os.path.basename(os.path.splitext(vcf_file)[0]));
	make_consensus_reference_from_vcf(assembly_sequence_file, vcf_file, consensus_sequence_file);
	run_dnadiff(reference_file, consensus_sequence_file, dnadiff_prefix);

	vcf_file = '%s/analysis-intermediate/consensus-BLASR-%s-cov_%d.variant.vcf' % (alignments_path, dataset_name, cov_threshold);
	consensus_sequence_file = '%s/consensus_sequence/%s.fasta' % (alignments_path, os.path.basename(os.path.splitext(vcf_file)[0]));
	dnadiff_prefix = '%s/dnadiff_cov_%d/%s/%s' % (alignments_path, cov_threshold, os.path.basename(os.path.splitext(vcf_file)[0]), os.path.basename(os.path.splitext(vcf_file)[0]));
	make_consensus_reference_from_vcf(assembly_sequence_file, vcf_file, consensus_sequence_file);
	run_dnadiff(reference_file, consensus_sequence_file, dnadiff_prefix);

	vcf_file = '%s/analysis-intermediate/consensus-marginAlign-%s-cov_%d.variant.vcf' % (alignments_path, dataset_name, cov_threshold);
	consensus_sequence_file = '%s/consensus_sequence/%s.fasta' % (alignments_path, os.path.basename(os.path.splitext(vcf_file)[0]));
	dnadiff_prefix = '%s/dnadiff_cov_%d/%s/%s' % (alignments_path, cov_threshold, os.path.basename(os.path.splitext(vcf_file)[0]), os.path.basename(os.path.splitext(vcf_file)[0]));
	make_consensus_reference_from_vcf(assembly_sequence_file, vcf_file, consensus_sequence_file);
	run_dnadiff(reference_file, consensus_sequence_file, dnadiff_prefix);

	vcf_file = '%s/analysis-intermediate/consensus-marginAlignGraphMap-%s-cov_%d.variant.vcf' % (alignments_path, dataset_name, cov_threshold);
	consensus_sequence_file = '%s/consensus_sequence/%s.fasta' % (alignments_path, os.path.basename(os.path.splitext(vcf_file)[0]));
	dnadiff_prefix = '%s/dnadiff_cov_%d/%s/%s' % (alignments_path, cov_threshold, os.path.basename(os.path.splitext(vcf_file)[0]), os.path.basename(os.path.splitext(vcf_file)[0]));
	make_consensus_reference_from_vcf(assembly_sequence_file, vcf_file, consensus_sequence_file);
	run_dnadiff(reference_file, consensus_sequence_file, dnadiff_prefix);

	vcf_file = '%s/analysis-intermediate/consensus-marginAlignGraphMap-anchor-%s-cov_%d.variant.vcf' % (alignments_path, dataset_name, cov_threshold);
	consensus_sequence_file = '%s/consensus_sequence/%s.fasta' % (alignments_path, os.path.basename(os.path.splitext(vcf_file)[0]));
	dnadiff_prefix = '%s/dnadiff_cov_%d/%s/%s' % (alignments_path, cov_threshold, os.path.basename(os.path.splitext(vcf_file)[0]), os.path.basename(os.path.splitext(vcf_file)[0]));
	make_consensus_reference_from_vcf(assembly_sequence_file, vcf_file, consensus_sequence_file);
	run_dnadiff(reference_file, consensus_sequence_file, dnadiff_prefix);



	vcf_file = '%s/analysis-intermediate/consensus-LAST-%s-uniquebest-cov_%d.variant.vcf' % (alignments_path, dataset_name, cov_threshold);
	consensus_sequence_file = '%s/consensus_sequence/%s.fasta' % (alignments_path, os.path.basename(os.path.splitext(vcf_file)[0]));
	dnadiff_prefix = '%s/dnadiff_cov_%d/%s/%s' % (alignments_path, cov_threshold, os.path.basename(os.path.splitext(vcf_file)[0]), os.path.basename(os.path.splitext(vcf_file)[0]));
	make_consensus_reference_from_vcf(assembly_sequence_file, vcf_file, consensus_sequence_file);
	run_dnadiff(reference_file, consensus_sequence_file, dnadiff_prefix);

	vcf_file = '%s/analysis-intermediate/consensus-DALIGNER-%s-uniquebest-cov_%d.variant.vcf' % (alignments_path, dataset_name, cov_threshold);
	consensus_sequence_file = '%s/consensus_sequence/%s.fasta' % (alignments_path, os.path.basename(os.path.splitext(vcf_file)[0]));
	dnadiff_prefix = '%s/dnadiff_cov_%d/%s/%s' % (alignments_path, cov_threshold, os.path.basename(os.path.splitext(vcf_file)[0]), os.path.basename(os.path.splitext(vcf_file)[0]));
	make_consensus_reference_from_vcf(assembly_sequence_file, vcf_file, consensus_sequence_file);
	run_dnadiff(reference_file, consensus_sequence_file, dnadiff_prefix);

	vcf_file = '%s/analysis-intermediate/consensus-marginAlign-%s-with_AS-uniquebest-cov_%d.variant.vcf' % (alignments_path, dataset_name, cov_threshold);
	consensus_sequence_file = '%s/consensus_sequence/%s.fasta' % (alignments_path, os.path.basename(os.path.splitext(vcf_file)[0]));
	dnadiff_prefix = '%s/dnadiff_cov_%d/%s/%s' % (alignments_path, cov_threshold, os.path.basename(os.path.splitext(vcf_file)[0]), os.path.basename(os.path.splitext(vcf_file)[0]));
	make_consensus_reference_from_vcf(assembly_sequence_file, vcf_file, consensus_sequence_file);
	run_dnadiff(reference_file, consensus_sequence_file, dnadiff_prefix);



	dnadiff_prefix = '%s/dnadiff_cov_%d/' % (alignments_path, cov_threshold);
	collect_consensus_sequence_comparison_results(dnadiff_prefix, '%s/collected-dnadiff-cov_%d.csv' % (dnadiff_prefix, cov_threshold));













#######################
### These are additional functions for the amplicon test (filtering and such stuff).
#######################
def filter_spanning_reads(dry_run, region, reads_path, sam_in_path, sam_out_folder, reference_path, leftalign=False):
	if not os.path.exists(sam_out_folder):
		sys.stderr.write('Creating folder "%s".\n' % (sam_out_folder));
		os.makedirs(sam_out_folder);

	sam_basename = os.path.basename(os.path.splitext(sam_in_path)[0]);

	# if (('marginalign' in os.path.basename(sam_in_path).lower()) or ('graphmap-params_20150525-all_reads-anchor_full_ref' in os.path.basename(sam_in_path).lower())):
	# 	sys.stderr.write('[-2] Removing special characters from the SAM qnames and rnames.\n');
	# 	execute_command_w_dryrun(dry_run, '%s/samfilter.py marginalign %s %s/%s-nospecialchars.sam' % (SAMSCRIPTS, sam_in_path, sam_out_folder, sam_basename));
	# 	sam_in_path = '%s/%s-nospecialchars.sam' % (sam_out_folder, os.path.basename(os.path.splitext(sam_in_path)[0]));

	# 	if ('marginalign' in os.path.basename(sam_in_path).lower()):
	# 		sys.stderr.write('[-1] Generating the alignment score so that alignments can be comparable.\n');
	# 		execute_command_w_dryrun(dry_run, '%s/samfilter.py generateAS %s %s/%s-nospecialchars.sam %s/%s-nospecialchars-with_AS.sam' % (SAMSCRIPTS, reference_path, sam_out_folder, sam_basename, sam_out_folder, sam_basename));
	# 		sam_in_path = '%s/%s-nospecialchars-with_AS.sam' % (sam_out_folder, os.path.basename(os.path.splitext(sam_in_path)[0]));

	sys.stderr.write('[-3] Fixing qnames and rnames in SAM files in case an aligner doesn\'t trim after whitespaces.\n');
	execute_command_w_dryrun(dry_run, '%s/samfilter.py fixhnames %s %s/%s-fxdhnames.sam' % (SAMSCRIPTS, sam_in_path, sam_out_folder, sam_basename));

	sys.stderr.write('[-2] Generating the alignment score in case it\'s not provided so that alignments can be comparable.\n');
	execute_command_w_dryrun(dry_run, '%s/samfilter.py generateAS %s %s/%s-fxdhnames.sam %s/%s-fxdhnames-with_AS.sam' % (SAMSCRIPTS, reference_path, sam_out_folder, sam_basename, sam_out_folder, sam_basename));

	sys.stderr.write('[-1] Check if base quality values are present, and add them if they are not. This is needed for LoFreq.\n');
	execute_command_w_dryrun(dry_run, '%s/samfilter.py addqv %s/%s-fxdhnames-with_AS.sam %s %s/%s-fxdhnames-with_AS-with_qv.sam' % (SAMSCRIPTS, sam_out_folder, sam_basename, reads_path, sam_out_folder, sam_basename));

	sam_in_path = '%s/%s-fxdhnames-with_AS-with_qv.sam' % (sam_out_folder, os.path.basename(os.path.splitext(sam_in_path)[0]));



	# if (('last' in os.path.basename(sam_in_path).lower()) or ('marginalign' in os.path.basename(sam_in_path).lower())):
	# 	sys.stderr.write('[0] Adding quality values to alignments...\n');
	# 	execute_command_w_dryrun(dry_run, '%s/samfilter.py addqv %s %s %s-addedqv.sam' % (SAMSCRIPTS, sam_in_path, reads_path, os.path.splitext(sam_in_path)[0]));
	# 	sys.stderr.write('\n');
	# 	sam_in_path = '%s-addedqv.sam' % (os.path.splitext(sam_in_path)[0]);

#	# temp_sam_region = '%s/%s-%s.sam' % (sam_out_folder, os.path.basename(os.path.splitext(sam_in_path)[0]), region[1]);
#	# temp_region_uniquebest = '%s/%s-uniquebest.sam' % (sam_out_folder, os.path.basename(os.path.splitext(temp_sam_region)[0]));
#	# temp_region_uniquebest_evalue = '%s/%s-evalue-1.sam' % (sam_out_folder, os.path.basename(os.path.splitext(temp_uniquebest_region)[0]));
#	# temp_before_2d = temp_uniquebest_region;

	temp_sam_uniquebest = '%s/%s-uniquebest.sam' % (sam_out_folder, os.path.basename(os.path.splitext(sam_in_path)[0]));
	temp_uniquebest_region = '%s/%s-%s.sam' % (sam_out_folder, os.path.basename(os.path.splitext(temp_sam_uniquebest)[0]), region[1]);
	temp_region_uniquebest_evalue = '%s/%s-evalue-1.sam' % (sam_out_folder, os.path.basename(os.path.splitext(temp_uniquebest_region)[0]));
	temp_before_2d = temp_uniquebest_region;

	### Extract only unique reads.
	sys.stderr.write('[1] Extracting only the unique reads with the best score "%s"...\n' % (temp_sam_uniquebest));
	execute_command_w_dryrun(dry_run, '%s/samfilter.py uniquebest %s %s' % (SAMSCRIPTS, sam_in_path, temp_sam_uniquebest));
	sys.stderr.write('\n');

	### Extract only the reads which fully cover the region:
	sys.stderr.write('[2] Extracting only reads that fully span region "%s" to file "%s"...\n' % (region[0], temp_uniquebest_region));
	execute_command_w_dryrun(dry_run, '%s/samfilter.py regionfull "%s" %s %s' % (SAMSCRIPTS, region[0], temp_sam_uniquebest, temp_uniquebest_region));
	# execute_command_w_dryrun('%s/samfilter.py regionpartial "%s" %s %s' % (SAMSCRIPTS, region[0], temp_sam_uniquebest, temp_uniquebest_region));
	sys.stderr.write('\n');

	if ('graphmap' in os.path.basename(sam_in_path).lower()):
		sys.stderr.write('[2.1] Filtering by E-value (GraphMap specific) to file "%s"...\n' % (temp_region_uniquebest_evalue));
		execute_command_w_dryrun(dry_run, '%s/samfilter.py evalue 1 %s %s' % (SAMSCRIPTS, temp_uniquebest_region, temp_region_uniquebest_evalue));
		sys.stderr.write('\n');
		temp_before_2d = temp_region_uniquebest_evalue;

	temp_1d = '%s/%s-1d.sam' % (sam_out_folder, os.path.basename(os.path.splitext(temp_before_2d)[0]));
	sys.stderr.write('[3] Extracting only the 1d reads to file "%s"...\n' % (temp_1d));
	execute_command_w_dryrun(dry_run, '%s/samfilter.py 1d %s %s' % (SAMSCRIPTS, temp_before_2d, temp_1d));
	sys.stderr.write('\n');

	temp_2d = '%s/%s-2d.sam' % (sam_out_folder, os.path.basename(os.path.splitext(temp_before_2d)[0]));
	sys.stderr.write('[4] Extracting only the 2d reads to file "%s"...\n' % (temp_2d));
	execute_command_w_dryrun(dry_run, '%s/samfilter.py 2d %s %s' % (SAMSCRIPTS, temp_before_2d, temp_2d));
	sys.stderr.write('\n');

	### Convert the SAM file to BAM format:
	sys.stderr.write('[5.1] Converting all aligned reads in the region to BAM format from file "%s"...' % (temp_before_2d));
	execute_command_w_dryrun(dry_run, '%s/convert_to_bam.sh %s' % (SAMSCRIPTS, os.path.splitext(temp_before_2d)[0]));
	bam_all = '%s-sorted.bam' % os.path.splitext(temp_before_2d)[0];
	sys.stderr.write('\n');
	sys.stderr.write('[5.2] Converting the 1d aligned reads to BAM format from file "%s"...' % (temp_1d));
	execute_command_w_dryrun(dry_run, '%s/convert_to_bam.sh %s' % (SAMSCRIPTS, os.path.splitext(temp_1d)[0]));
	bam_1d = '%s-sorted.bam' % os.path.splitext(temp_1d)[0];
	sys.stderr.write('\n');
	sys.stderr.write('[5.3] Converting the 2d aligned reads to BAM format from file "%s"...' % (temp_2d));
	execute_command_w_dryrun(dry_run, '%s/convert_to_bam.sh %s' % (SAMSCRIPTS, os.path.splitext(temp_2d)[0]));
	bam_2d = '%s-sorted.bam' % os.path.splitext(temp_2d)[0];
	sys.stderr.write('\n');

	bam_all_reads_in_region = bam_all;
	bam_1d_reads_in_region = bam_1d
	bam_2d_reads_in_region = bam_2d

	if (leftalign == True):
		if (reference_path == None):
			sys.stderr.write('ERROR: Reference not specified, cannot leftalign the BAM file.\n');
		else:
			sys.stderr.write('[6.1] Left aligning all aligned reads in the region from file "%s"...' % (bam_all));
			bam_all_leftaligned = '%s-leftaligned.bam' % (os.path.splitext(bam_all)[0]);
			execute_command_w_dryrun(dry_run, 'cat %s | %s/bamleftalign -f %s > %s' % (bam_all, TOOLS_PATH, reference_path, bam_all_leftaligned));
			execute_command_w_dryrun(dry_run, 'samtools index %s' % (bam_all_leftaligned));
			sys.stderr.write('\n');

			sys.stderr.write('[6.2] Left aligning all aligned reads in the region from file "%s"...' % (bam_1d));
			bam_1d_leftaligned = '%s-leftaligned.bam' % (os.path.splitext(bam_1d)[0]);
			execute_command_w_dryrun(dry_run, 'cat %s | %s/bamleftalign -f %s > %s' % (bam_1d, TOOLS_PATH, reference_path, bam_1d_leftaligned));
			execute_command_w_dryrun(dry_run, 'samtools index %s' % (bam_1d_leftaligned));
			sys.stderr.write('\n');

			sys.stderr.write('[6.3] Left aligning all aligned reads in the region from file "%s"...' % (bam_2d));
			bam_2d_leftaligned = '%s-leftaligned.bam' % (os.path.splitext(bam_2d)[0]);
			execute_command_w_dryrun(dry_run, 'cat %s | %s/bamleftalign -f %s > %s' % (bam_2d, TOOLS_PATH, reference_path, bam_2d_leftaligned));
			execute_command_w_dryrun(dry_run, 'samtools index %s' % (bam_2d_leftaligned));
			sys.stderr.write('\n');

			bam_all_reads_in_region = bam_all_leftaligned;
			bam_1d_reads_in_region = bam_1d_leftaligned
			bam_2d_reads_in_region = bam_2d_leftaligned

	return [bam_all_reads_in_region, bam_1d_reads_in_region, bam_2d_reads_in_region];

def filter_partial_reads(dry_run, region, reads_path, sam_in_path, sam_out_folder, reference_path, leftalign=False):
	if not os.path.exists(sam_out_folder):
		sys.stderr.write('Creating folder "%s".\n' % (sam_out_folder));
		os.makedirs(sam_out_folder);

	sam_basename = os.path.basename(os.path.splitext(sam_in_path)[0]);

	sys.stderr.write('[-3] Fixing qnames and rnames in SAM files in case an aligner doesn\'t trim after whitespaces.\n');
	execute_command_w_dryrun(dry_run, '%s/samfilter.py fixhnames %s %s/%s-fxdhnames.sam' % (SAMSCRIPTS, sam_in_path, sam_out_folder, sam_basename));

	sys.stderr.write('[-2] Generating the alignment score in case it\'s not provided so that alignments can be comparable.\n');
	execute_command_w_dryrun(dry_run, '%s/samfilter.py generateAS %s %s/%s-fxdhnames.sam %s/%s-fxdhnames-with_AS.sam' % (SAMSCRIPTS, reference_path, sam_out_folder, sam_basename, sam_out_folder, sam_basename));

	sys.stderr.write('[-1] Check if base quality values are present, and add them if they are not. This is needed for LoFreq.\n');
	execute_command_w_dryrun(dry_run, '%s/samfilter.py addqv %s/%s-fxdhnames-with_AS.sam %s %s/%s-fxdhnames-with_AS-with_qv.sam' % (SAMSCRIPTS, sam_out_folder, sam_basename, reads_path, sam_out_folder, sam_basename));

	sam_in_path = '%s/%s-fxdhnames-with_AS-with_qv.sam' % (sam_out_folder, os.path.basename(os.path.splitext(sam_in_path)[0]));

	temp_sam_uniquebest = '%s/%s-uniquebest.sam' % (sam_out_folder, os.path.basename(os.path.splitext(sam_in_path)[0]));
	temp_uniquebest_region = '%s/%s-%s.sam' % (sam_out_folder, os.path.basename(os.path.splitext(temp_sam_uniquebest)[0]), region[1]);
	temp_region_uniquebest_evalue = '%s/%s-evalue-1.sam' % (sam_out_folder, os.path.basename(os.path.splitext(temp_uniquebest_region)[0]));
	temp_before_2d = temp_uniquebest_region;

	### Extract only unique reads.
	sys.stderr.write('[1] Extracting only the unique reads with the best score "%s"...\n' % (temp_sam_uniquebest));
	execute_command_w_dryrun(dry_run, '%s/samfilter.py uniquebest %s %s' % (SAMSCRIPTS, sam_in_path, temp_sam_uniquebest));
	sys.stderr.write('\n');

	### Extract only the reads which partially cover the region:
	sys.stderr.write('[2] Extracting only reads that partially span region "%s" to file "%s"...\n' % (region[0], temp_uniquebest_region));
	execute_command_w_dryrun(dry_run, '%s/samfilter.py regionpartial "%s" %s %s' % (SAMSCRIPTS, region[0], temp_sam_uniquebest, temp_uniquebest_region));
	# execute_command_w_dryrun('%s/samfilter.py regionpartial "%s" %s %s' % (SAMSCRIPTS, region[0], temp_sam_uniquebest, temp_uniquebest_region));
	sys.stderr.write('\n');

	if ('graphmap' in os.path.basename(sam_in_path).lower()):
		execute_command_w_dryrun(False, 'grep -v "^@" %s | wc -l' % (temp_before_2d));
		sys.stderr.write('[2.1] Filtering by E-value (GraphMap specific) to file "%s"...\n' % (temp_region_uniquebest_evalue));
		execute_command_w_dryrun(dry_run, '%s/samfilter.py evalue 1 %s %s' % (SAMSCRIPTS, temp_uniquebest_region, temp_region_uniquebest_evalue));
		sys.stderr.write('\n');
		temp_before_2d = temp_region_uniquebest_evalue;

	execute_command_w_dryrun(False, 'grep -v "^@" %s | wc -l' % (temp_before_2d));

	return;
#	return [bam_all_reads_in_region, bam_1d_reads_in_region, bam_2d_reads_in_region];

def evaluate_vcf(mapper_name, region_name, vcf_file, vcf_known_mutations_path):
	if (not os.path.exists(vcf_file)):
		sys.stderr.write('ERROR: File "%s" does not exist! Continuing.\n' % (vcf_file));
		return [0]*10;

	# if (vcf_mutations_path.endswith('.gz') == False and os.path.exists(vcf_mutations_path + '.gz') == False):
	execute_command_w_dryrun(False, 'bgzip -c %s > %s.gz' % (vcf_known_mutations_path, vcf_known_mutations_path));
	execute_command_w_dryrun(False, 'tabix -f -p vcf %s.gz' % (vcf_known_mutations_path));
	if (vcf_known_mutations_path.endswith('.gz') == False):
		vcf_known_mutations_path_gz = vcf_known_mutations_path + '.gz';
		# print vcf_mutations_path;
	else:
		sys.stderr.write('Truth variants not in .VCF format, but gzipped! Extract them and give some text files!\n');
		return [0]*10;

	execute_command_w_dryrun(False, 'bgzip -c %s > %s.gz' % (vcf_file, vcf_file));
	execute_command_w_dryrun(False, 'tabix -f -p vcf %s.gz' % (vcf_file));
	# execute_command('bamsurgeon/etc/evaluator.py -v %s.gz -t %s -m SNV' % (vcf_file, vcf_mutations_path));
	fn_file = '%s-fn.vcf' % (os.path.splitext(vcf_file)[0]);
	fp_file = '%s-fp.vcf' % (os.path.splitext(vcf_file)[0]);

	command = '%s/etc/evaluator.py -v %s.gz -t %s -m SNV --fn %s --fp %s' % (BAMSURGEON_PATH, vcf_file, vcf_known_mutations_path_gz, fn_file, fp_file);
	[rc, output, err] = execute_command_with_ret(False, command);

	if (rc != 0):
		sys.stderr.write(str(err));
		exit(1);
	lines = output.split('\n');

	results = [0, 0, 0, 0];
	if (len(lines) > 1 and len(lines[1]) > 0):
		counts_line_id = -1;
		i = 0;
		for line in lines:
			i += 1;
			if (line.startswith('tpcount, fpcount, subrecs, trurecs:')):
				counts_line_id = i;
				break;

		# sys.stderr.write('Lines:\n%s\n\n' % ('\n'.join(lines)));
		# sys.stderr.write('lines[%d].split = %s\n' % (counts_line_id, str(lines[counts_line_id].split())));
		if (counts_line_id > 0):
			results = [int(value) for value in lines[counts_line_id].split()];

	sys.stderr.write('Truth pass counting...\n');
	sys.stderr.write(vcf_known_mutations_path + '\n');
	[num_pass_snps, num_nonpass_snps] = vcffilter.count_nonpass_variants(vcf_known_mutations_path, verbose=False);
	results.append(num_pass_snps);
	results.append(num_nonpass_snps);

	sys.stderr.write('FP pass counting...\n');
	[num_pass_snps, num_nonpass_snps] = vcffilter.count_nonpass_variants(fp_file, verbose=False);
	results.append(num_pass_snps);
	results.append(num_nonpass_snps);

	sys.stderr.write('FN pass counting...\n');
	[num_pass_snps, num_nonpass_snps] = vcffilter.count_nonpass_variants(fn_file, verbose=False);
	results.append(num_pass_snps);
	results.append(num_nonpass_snps);

	return results;













def filter_fastq_for_marginalign(input_fastq_path, out_fastq_path, fp_out):
	try:
		fp_in = open(input_fastq_path, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Skipping.\n' % input_fastq_path);
		return {};
	if (fp_out == None):
		if (input_fastq_path == out_fastq_path):
			sys.stderr.write('ERROR: Output and input files are the same! Skipping.\n');
			return {};
		try:
			fp_out = open(out_fastq_path, 'w');
		except:
			sys.stderr.write('ERROR: Could not open file "%s" for writing! Skipping.\n' % out_fastq_path);
			return {};

	num_matches = 0;
	header_hash = {};

	i = 0;
	while True:
		i += 1;
		[header, read] = fastqparser.get_single_read(fp_in);
		
		if (len(read) == 0):
			sys.stderr.write('Returning from filter_fastq_for_marginalign. i = %d\n' % (i));
			break;

		# read[0] = read[0][0] + read[0][1:].replace();
		# if (len(read[1]) <= 50000):
			# read[0] = read[0][0] + re.sub('[^0-9a-zA-Z]', '_', read[0][1:]); # re.sub("[|:", "_", read[0][1:]);
		new_header = read[0][0] + re.sub('[^0-9a-zA-Z]', '_', read[0][1:]); # re.sub("[|:", "_", read[0][1:]);
		header_hash[new_header[1:]] = read[0][1:];
		read[0] = new_header;
		fp_out.write('\n'.join(read) + '\n');

	sys.stderr.write('\n');
	fp_in.close();
	return header_hash;

def filter_sam_for_marginalign(sam_file, out_filtered_sam_file):
	fp_in = None;
	fp_out = None;
	
	rname_hash = {};
	qname_hash = {};

	try:
		fp_in = open(sam_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, sam_file));
		exit(1);
	
	try:
		fp_out = open(out_filtered_sam_file, 'w');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!' % (__name__, out_filtered_sam_file));
		exit(1);

	num_accepted = 0;
	num_rejected = 0;
	num_unmapped = 0;
	
	i = 0;
	for line in fp_in:
		line = line.strip();
		if (len(line) == 0 or line[0] == '@'):
			if (line.startswith('@SQ')):
				split_line = line.split('\t');
				found_hname = False;
				for param in split_line:
					if (param.startswith('SN:')):
						hname = param.split(':')[-1];
						new_hname = re.sub('[^0-9a-zA-Z]', '_', hname);
						sys.stderr.write('Found hname: "%s", replacing with: "%s".\n' % (hname, new_hname));
						new_line = line.replace(hname, new_hname);
						rname_hash[new_hname] = hname;
						fp_out.write(new_line + '\n');
						found_hname = True;
						break;
				if (found_hname == False):
					fp_out.write(line + '\n');
			else:
				fp_out.write(line + '\n');
			continue;

		i += 1;
		sys.stderr.write('\rLine %d, num_accepted: %d, num_rejected: %d' % (i, num_accepted, num_rejected));

		sam_line = utility_sam.SAMLine(line.rstrip());
		# split_line = line.split('\t');

		if (sam_line.IsMapped() == False):
			fp_out.write(line + '\n');
			num_unmapped += 1;

		qname = sam_line.qname;
		rname = sam_line.rname;

		new_line = line + '';

		if (qname != '*'):
			new_qname = re.sub('[^0-9a-zA-Z]', '_', qname);
			new_line = new_line.replace(qname, new_qname);
			qname_hash[new_qname] = qname;

		if (rname != '*'):
			new_rname = re.sub('[^0-9a-zA-Z]', '_', rname);
			new_line = new_line.replace(rname, new_rname);
			rname_hash[new_rname] = rname;

		fp_out.write(new_line + '\n');
		num_accepted += 1;

	fp_in.close();
	fp_out.close();
	
	sys.stderr.write('\n');
	sys.stderr.write('Done!\n');
	try:
		sys.stderr.write('num_accepted = %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
		sys.stderr.write('num_rejected = %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));
		sys.stderr.write('num_unmapped = %d (%.2f%%)\n' % (num_unmapped, (float(num_unmapped) / float(num_accepted + num_rejected)) * 100.0));
	except:
		pass;

	return [qname_hash, rname_hash];


def fix_sam_qnames_after_marginAlign(input_sam_path, ref_header_hash, read_header_hash, out_sam_path):
	if (input_sam_path == out_sam_path):
		sys.stderr.write('ERROR: Input and output SAM files are the same! Skipping the update of qname and rname values.\n');
		return;
	try:
		fp_in = open(input_sam_path, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_sam_path);
		exit(0);
	try:
		fp_out = open(out_sam_path, 'w');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for writing! Exiting.\n' % out_sam_path);
		exit(0);
	i = 0;
	for line in fp_in:
		i += 1;
		line = line.strip();
		if (len(line) == 0 or line[0] == '@'):
			if (line.startswith('@SQ')):
				split_line = line.split('\t');
				found_hname = False;
				for param in split_line:
					if (param.startswith('SN:')):
						hname = param.split(':')[-1];
						try:
							original_hname = ref_header_hash[hname];
							sys.stderr.write('Found hname: "%s".\n' % (hname));
						except:
							original_hname = hname;
							sys.stderr.write('Could not find hname "%s".\n' % (hname));

						new_line = line.replace(hname, original_hname.split()[0]);	### Split on whitespaces to report only the gi part of the header.
						fp_out.write(new_line + '\n');
						found_hname = True;
						break;
				if (found_hname == False):
					fp_out.write(line + '\n');

			else:
				fp_out.write(line + '\n');
			continue;
		sys.stderr.write('\rLine %d' % (i));

		split_line = line.split('\t');
		qname = split_line[0];
		rname = split_line[2];

		if (qname == '*' or rname == '*'):
			fp_out.write(line + '\n');
			continue;

		try:
			original_qname = read_header_hash[qname];
		except:
			original_qname = qname;
		try:
			original_rname = ref_header_hash[rname];
		except:
			original_rname = qname;

		new_line = line.replace(qname, original_qname.split()[0]);
		new_line = new_line.replace(rname, original_rname.split()[0]);
		fp_out.write(new_line + '\n');
	sys.stderr.write('\n');




















def execute_command(command):
	sys.stderr.write('[Executing] "%s"\n' % (command));
	subprocess.call(command, shell=True);
	sys.stderr.write('\n');

def execute_command_get_stdout(command):
	sys.stderr.write('[Executing] "%s"\n' % (command));
	p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
	[out, err] = p.communicate()
	sys.stderr.write('\n');
	return [out, err];

def execute_command_w_dryrun(dry_run, command):
	sys.stderr.write('[Executing] "%s"\n' % command);
	if (dry_run == False):
		subprocess.call(command, shell=True);
	sys.stderr.write('\n');

def execute_command_with_ret(dry_run, command):
	sys.stderr.write('Executing command: "%s"\n' % command);
	if (dry_run == False):
		p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE);
	[output, err] = p.communicate()
	rc = p.returncode
	sys.stderr.write('\n');
	return [rc, output, err];

# Finds all files with a given filename pattern starting from the given path, recursivelly.
def find_files(start_path, file_pattern, max_depth=-1):
	matches = []
	for root, dirnames, filenames in os.walk(start_path):
		for filename in fnmatch.filter(filenames, file_pattern):
			depth = len(root.split('/'));
			if (max_depth < 0 or (max_depth >= 0 and depth <= max_depth)):
				matches.append(os.path.join(root, filename))

	matches.sort();
	return matches;

def measure_command_wrapper(out_filename):
#	if (USE_BASICDEFINES_ == True):
#		return basicdefines.measure_command(out_filename);
#	else:
	return '/usr/bin/time --format "Command line: %%C\\nReal time: %%e s\\nCPU time: -1.0 s\\nUser time: %%U s\\nSystem time: %%S s\\nMaximum RSS: %%M kB\\nExit status: %%x" --quiet -o %s ' % out_filename;



def verbose_usage_and_exit():
	sys.stderr.write('Runs all analyses used in the GraphMap paper.\n\n');
	sys.stderr.write('Usage:\n');
	sys.stderr.write('\t%s mode\n' % sys.argv[0]);
	sys.stderr.write('\n');
	sys.stderr.write('\t- mode - either "run-simdata", "run-realdata", "setup-tools" or "setup-data"\n');
	sys.stderr.write('\n');

	exit(0);

if __name__ == "__main__":
	if (len(sys.argv) < 2 or len(sys.argv) > 7):
		verbose_usage_and_exit();

	if (sys.argv[1] == 'setup-tools'):
		setup_tools();
	elif (sys.argv[1] == 'setup-data'):
		setup_data();
	elif (sys.argv[1] == 'run-simdata'):
		if (MODULE_VCFFILTER == False or MODULE_FASTQPARSER == False or MODULE_UTILITYSAM == False)
			sys.stderr.write('Please run "%s setup-tools" first! Exiting.\n\n' % (sys.argv[0]));
			exit(1);
		run_simulated_data();
	elif (sys.argv[1] == 'run-realdata'):
		if (MODULE_VCFFILTER == False or MODULE_FASTQPARSER == False or MODULE_UTILITYSAM == False)
			sys.stderr.write('Please run "%s setup-tools" first! Exiting.\n\n' % (sys.argv[0]));
			exit(1);
		run_real_data();
	else:
		verbose_usage_and_exit();
