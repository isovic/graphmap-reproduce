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
# SAMSCRIPTS = '%s/../../samscripts/src/' % (SCRIPT_PATH);
sys.path.append('%s/../../samscripts/src/' % (SCRIPT_PATH));

def move_data_to_proper_folders(downloads_folder):
	sys.stderr.write('Copying the downloaded data to the data/ folder structure used in the tests.\n');

	sys.stderr.write('\tFinalizing the Mikheyev&Tin Lambda R6 dataset.\n');	
	folder_name = 'lambdaR6';
	fastq_file = 'lambda_reads-enumerated.fastq';
	reads_folder = '%s/../data/consensus-lambdaR6/reads/' % (SCRIPT_PATH);
	execute_command('mkdir -p %s; cp %s/%s/%s %s/'  % (reads_folder, downloads_folder, folder_name, fastq_file, reads_folder));

	sys.stderr.write('\tFinalizing the E. Coli R7 dataset.\n');	
	folder_name = 'R7';
	fastq_file = 'Ecoli_R7_CombinedFasta.fastq';
	fastq_file_2d = 'Ecoli_R7_CombinedFasta-2d.fastq';
	reads_folder = '%s/../data/consensus-ecoliR7/reads' % (SCRIPT_PATH);
	### All reads.
	execute_command('mkdir -p %s; cp %s/%s/%s %s/'  % (reads_folder, downloads_folder, folder_name, fastq_file, reads_folder));
	### 2d reads only.
	execute_command('mkdir -p %s; cp %s/%s/%s %s/'  % (reads_folder, downloads_folder, folder_name, fastq_file_2d, reads_folder));

	sys.stderr.write('\tFinalizing the E. Coli R7.3 dataset.\n');	
	folder_name = 'R7.3';
	reads_folder = '%s/../data/consensus-ecoliR7.3-1d/reads/' % (SCRIPT_PATH);
	fastq_file_1d = 'ecoliR7.3-1d.fastq';
	fastq_file_2d = 'ecoliR7.3-2d.fastq';
	### All reads.
	execute_command('mkdir -p %s; cp %s/%s/%s %s/'  % (reads_folder, downloads_folder, folder_name, fastq_file_1d, reads_folder));
	### 2d reads only.
	execute_command('mkdir -p %s; cp %s/%s/%s %s/'  % (reads_folder, downloads_folder, folder_name, fastq_file_2d, reads_folder));

	### Combining the 2d reads from the R7 and R7.3 datasets for SV tests.
	sys.stderr.write('\tFinalizing the SV dataset (2d reads from R7 and R7.3).\n');	
	folder_name = 'all_2d_for_sv';
	reads_folder = '%s/../data/sv/reads/' % (SCRIPT_PATH);
	fastq_file_2d = 'all_2d_for_sv.fastq';
	execute_command('mkdir -p %s; cp %s/%s/%s %s/'  % (reads_folder, downloads_folder, folder_name, fastq_file_2d, reads_folder));

	sys.stderr.write('\tFinalizing the UTI89 dataset.\n');	
	folder_name = 'ecoli-uti89';
	fastq_file = 'reads_ecoli_uti89.fastq';
	reads_folder = '%s/../data/consensus-uti89/reads/' % (SCRIPT_PATH);
	execute_command('mkdir -p %s; cp %s/%s/%s %s/'  % (reads_folder, downloads_folder, folder_name, fastq_file, reads_folder));

	sys.stderr.write('\tFinalizing the Salmonella dataset.\n');	
	folder_name = 'salmonella-typhi';
	fastq_file = 'reads-typhi.fastq';
	reads_folder = '%s/../data/consensus-salmonellatyphi/reads/' % (SCRIPT_PATH);
	execute_command('mkdir -p %s; cp %s/%s/%s %s/'  % (reads_folder, downloads_folder, folder_name, fastq_file, reads_folder));

	sys.stderr.write('\tFinalizing the amplicon sequencing dataset.\n');	
	folder_name = 'amplicons';
	fastq_file = 'reads_all-f1000.fastq';
	reads_folder = '%s/../data/amplicons-f1000-1d2d/reads/' % (SCRIPT_PATH);
	execute_command('mkdir -p %s/../data/amplicons-f1000-1d2d/reads/; cp %s/amplicons/reads_all-f1000.fastq %s/../data/amplicons-f1000-1d2d/reads/'  % (SCRIPT_PATH, downloads_folder, SCRIPT_PATH));
	execute_command('mkdir -p %s/../data/amplicons-f1000-1d2d/reads/; cp %s/amplicons/reads_2d.fastq %s/../data/amplicons-f1000-1d2d/reads/'  % (SCRIPT_PATH, downloads_folder, SCRIPT_PATH));
	execute_command('mkdir -p %s/../data/amplicons-f1000/reads/; cp %s/amplicons/reads_all-f1000.fastq %s/../data/amplicons-f1000/reads/'  % (SCRIPT_PATH, downloads_folder, SCRIPT_PATH));
	execute_command('mkdir -p %s/../data/amplicons-f1000/reads/; cp %s/amplicons/reads_2d.fastq %s/../data/amplicons-f1000/reads/'  % (SCRIPT_PATH, downloads_folder, SCRIPT_PATH));

	sys.stderr.write('\tFinalizing the BE1 dataset.\n');	
	folder_name = 'be1';
	fastq_file = 'reads-BE1.fastq';
	reads_folder = '%s/../data/consensus-BE1/reads/' % (SCRIPT_PATH);
	execute_command('mkdir -p %s; cp %s/%s/%s %s/'  % (reads_folder, downloads_folder, folder_name, fastq_file, reads_folder));

	sys.stderr.write('\tFinalizing the ADP1 dataset.\n');	
	folder_name = 'adp1';
	fastq_file = 'all_1d2d.fastq';
	reads_folder = '%s/../data/consensus-ADP1/reads/' % (SCRIPT_PATH);
	execute_command('mkdir -p %s; cp %s/%s/%s %s/'  % (reads_folder, downloads_folder, folder_name, fastq_file, reads_folder));

	sys.stderr.write('\tFinalizing the Loman et al. Nature Methods nanopore assembly dataset.\n');	
	folder_name = 'ecoli-nmeth';
	fastq_file = 'reads-nmeth-all_2d.fastq';
	reads_folder = '%s/../data/consensus-ecolinmeth/reads/' % (SCRIPT_PATH);
	execute_command('mkdir -p %s; cp %s/%s/%s %s/'  % (reads_folder, downloads_folder, folder_name, fastq_file, reads_folder));

	sys.stderr.write('\tFinalizing the MAP006 dataset.\n');	
	folder_name = 'ecoli-map006';
	fastq_file = 'reads-MAP006-1.fastq';
	reads_folder = '%s/../data/consensus-ecoliMAP006/reads/' % (SCRIPT_PATH);
	execute_command('mkdir -p %s; cp %s/%s/%s %s/'  % (reads_folder, downloads_folder, folder_name, fastq_file, reads_folder));



def get_real_data():
	yes_no = raw_input("Downloading and unpacking raw nanopore data will consume approximately ~800GB of disk space. Continue? [y/n] " % (MAPPER_NAME));
	if (yes_no != 'y'):
		return;

	sys.stderr.write('Downloading raw nanopore data and unpacking.\n');

	sys.stderr.write('Creating downloads folder in path: "%s/../data/downloads".\n' % (SCRIPT_PATH));
	if (not os.path.exists('%s/../data/downloads' % (SCRIPT_PATH))):
		os.makedirs('%s/../data/downloads' % (SCRIPT_PATH));

	sys.stderr.write('\tFetching the Mikheyev&Tin Lambda R6 dataset.\n');	
	execute_command('cd %s/../data/downloads/; mkdir lambdaR6; cd lambdaR6; git clone https://github.com/mikheyev/MinION-review.git' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/lambdaR6/MinION-review/data/; gunzip --keep lambda_reads_1d.1.fastq.gz' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/lambdaR6/MinION-review/data/; gunzip --keep lambda_reads_1d.2.fastq.gz' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/lambdaR6/MinION-review/data/; gunzip --keep lambda_reads_2d.fastq.gz' % (SCRIPT_PATH));
	sys.stderr.write('\tExtracting the Mikheyev&Tin Lambda R6 dataset.\n');	
	execute_command('cd %s/../data/downloads/lambdaR6/; cat MinION-review/data/lambda_reads_1d.1.fastq > lambda_reads.fastq' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/lambdaR6/; cat MinION-review/data/lambda_reads_1d.2.fastq >> lambda_reads.fastq' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/lambdaR6/; cat MinION-review/data/lambda_reads_2d.fastq >> lambda_reads.fastq' % (SCRIPT_PATH));
	execute_command('%s/fastqfilter.py enumerate %s/../data/downloads/lambdaR6/lambda_reads.fastq %s/../data/downloads/lambdaR6/lambda_reads-enumerated.fastq' % (SAMSCRIPTS, SCRIPT_PATH, SCRIPT_PATH));



	sys.stderr.write('\tFetching the R7 dataset.\n');
	execute_command('cd %s/../data/downloads/; mkdir R7' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd R7; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA362/ERA362836/oxfordnanopore_native/Ecoli_R7_NONI.tgz' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd R7; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA362/ERA362836/oxfordnanopore_native/Ecoli_R7_ONI_flowcell_18.tar.gz' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd R7; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA362/ERA362836/oxfordnanopore_native/Ecoli_R7_ONI_flowcell_17.tar.gz' % (SCRIPT_PATH));
	sys.stderr.write('\tExtracting the R7 dataset.\n');
	folder_name = 'R7';
	archive_basename = 'Ecoli_R7_NONI';					execute_command('cd %s/../data/downloads/; cd %s; mkdir %s; tar -xzvf %s.tgz -C %s/' % (SCRIPT_PATH, folder_name, archive_basename, archive_basename, archive_basename));
	archive_basename = 'Ecoli_R7_ONI_flowcell_18';		execute_command('cd %s/../data/downloads/; cd %s; mkdir %s; tar -xzvf %s.tar.gz -C %s/' % (SCRIPT_PATH, folder_name, archive_basename, archive_basename, archive_basename));
	archive_basename = 'Ecoli_R7_ONI_flowcell_17';		execute_command('cd %s/../data/downloads/; cd %s; mkdir %s; tar -xzvf %s.tar.gz -C %s/' % (SCRIPT_PATH, folder_name, archive_basename, archive_basename, archive_basename));
	### Prepare all reads.
	folder_name = 'R7';
	fastq_file = 'Ecoli_R7_CombinedFasta.fastq';
	reads_folder = '%s/../data/consensus-ecoliR7/reads' % (SCRIPT_PATH);
	archive_basename = 'Ecoli_R7_NONI'; execute_command('cd %s/../data/downloads/%s; poretools fastq %s/downloads/ > %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file));
	archive_basename = 'Ecoli_R7_ONI_flowcell_18'; execute_command('cd %s/../data/downloads/%s; poretools fastq %s/ >> %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file));
	archive_basename = 'Ecoli_R7_ONI_flowcell_17'; execute_command('cd %s/../data/downloads/%s; poretools fastq %s/ >> %s'  % (SCRIPT_PATH, folder_name, archive_basename, fastq_file));
	### Prepare only the 2d reads.
	folder_name = 'R7';
	reads_folder = '%s/../data/consensus-ecoliR7/reads' % (SCRIPT_PATH);
	fastq_file_2d = 'Ecoli_R7_CombinedFasta-2d.fastq';
	execute_command('mkdir -p %s/../data/downloads/%s' % (SCRIPT_PATH, folder_name));
	archive_basename = 'Ecoli_R7_NONI'; execute_command('cd %s/../data/downloads/%s; poretools fastq --type 2D %s/downloads/ > %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file_2d));
	archive_basename = 'Ecoli_R7_ONI_flowcell_18'; execute_command('cd %s/../data/downloads/%s; poretools fastq --type 2D %s/ >> %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file_2d));
	archive_basename = 'Ecoli_R7_ONI_flowcell_17'; execute_command('cd %s/../data/downloads/%s; poretools fastq --type 2D %s/ >> %s'  % (SCRIPT_PATH, folder_name, archive_basename, fastq_file_2d));

	sys.stderr.write('\tFetching the R7.3 data.\n');
	execute_command('cd %s/../data/downloads/; mkdir R7.3' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd R7.3; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA362/ERA362836/oxfordnanopore_native/Ecoli_R73.tgz' % (SCRIPT_PATH));
	sys.stderr.write('\tExtracting the R7.3 data.\n');
	folder_name = 'R7.3';
	archive_basename = 'Ecoli_R73';					execute_command('cd %s/../data/downloads/; cd %s; mkdir %s; tar -xzvf %s.tgz -C %s/' % (SCRIPT_PATH, folder_name, archive_basename, archive_basename, archive_basename));
	folder_name = 'R7.3';
	fastq_file = 'ecoliR7.3.fastq';
	reads_folder = '%s/../data/consensus-ecoliR7.3/reads/' % (SCRIPT_PATH);
	archive_basename = 'Ecoli_R73'; execute_command('cd %s/../data/downloads/%s; poretools fastq %s/downloads/ > %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file));
	execute_command('mkdir -p %s; cp %s/../data/downloads/%s/%s %s/'  % (reads_folder, SCRIPT_PATH, folder_name, fastq_file, reads_folder));
	### Prepare all reads.
	reads_folder = '%s/../data/consensus-ecoliR7.3-1d/reads/' % (SCRIPT_PATH);
	fastq_file_1d = 'ecoliR7.3-1d.fastq';
	archive_basename = 'Ecoli_R73'; execute_command('cd %s/../data/downloads/%s; poretools fastq --type fwd,rev %s/downloads/ > %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file_1d));
	### Prepare only the 2d reads.
	folder_name = 'R7.3';
	reads_folder = '%s/../data/consensus-ecoliR7.3-2d/reads/' % (SCRIPT_PATH);
	fastq_file_2d = 'ecoliR7.3-2d.fastq';
	execute_command('mkdir -p %s/../data/downloads/%s' % (SCRIPT_PATH, folder_name));
	archive_basename = 'Ecoli_R73'; execute_command('cd %s/../data/downloads/%s; poretools fastq --type 2D %s/downloads/ > %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file_2d));

	### Extracting the 2d reads from the R7 and R7.3 datasets.
	folder_name = 'all_2d_for_sv';
	reads_folder = '%s/../data/sv/reads/' % (SCRIPT_PATH);
	fastq_file_2d = 'all_2d_for_sv.fastq';
	execute_command('mkdir -p %s/../data/downloads/%s' % (SCRIPT_PATH, folder_name));
	execute_command('cd %s/../data/downloads/%s; cat %s/../data/consensus-ecoliR7/reads/Ecoli_R7_CombinedFasta-2d.fastq > %s'  % (SCRIPT_PATH, folder_name, SCRIPT_PATH, fastq_file_2d));
	execute_command('cd %s/../data/downloads/%s; cat %s/../data/consensus-ecoliR7.3-2d/reads/ecoliR7.3-2d.fastq >> %s'  % (SCRIPT_PATH, folder_name, SCRIPT_PATH, fastq_file_2d));

	sys.stderr.write('\tFetching the E. Coli UTI89 dataset, generated in-house.\n');	
	execute_command('cd %s/../data/downloads/; mkdir ecoli-uti89' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd ecoli-uti89; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA444/ERA444730/oxfordnanopore_native/reads.tar' % (SCRIPT_PATH));
	sys.stderr.write('\tExtracting the E. Coli UTI89 dataset, generated in-house.\n');	
	folder_name = 'ecoli-uti89';
	archive_basename = 'reads';					execute_command('cd %s/../data/downloads/; cd %s; mkdir %s; tar -xvf %s.tar -C %s/' % (SCRIPT_PATH, folder_name, archive_basename, archive_basename, archive_basename));
	folder_name = 'ecoli-uti89';
	fastq_file = 'reads_ecoli_uti89.fastq';
	reads_folder = '%s/../data/consensus-uti89/reads/' % (SCRIPT_PATH);
	archive_basename = 'reads'; execute_command('cd %s/../data/downloads/%s; poretools fastq %s/fast5/ > %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file));

	sys.stderr.write('\tFetching the Salmonella Typhi dataset.\n');	
	execute_command('cd %s/../data/downloads/; mkdir salmonella-typhi' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd salmonella-typhi; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA375/ERA375685/oxfordnanopore_native/H566_ON_inc.tar.gz' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd salmonella-typhi; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA375/ERA375987/oxfordnanopore_native/H566_30_min_inc.tar.gz' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd salmonella-typhi; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA376/ERA376255/oxfordnanopore_native/raw_2_rabsch_R7.tar.gz' % (SCRIPT_PATH));
	sys.stderr.write('\tExtracting the Salmonella Typhi dataset.\n');	
	folder_name = 'salmonella-typhi';
	archive_basename = 'H566_ON_inc';			execute_command('cd %s/../data/downloads/; cd %s; mkdir %s; tar -xzvf %s.tar.gz -C %s/' % (SCRIPT_PATH, folder_name, archive_basename, archive_basename, archive_basename));
	archive_basename = 'H566_30_min_inc';		execute_command('cd %s/../data/downloads/; cd %s; mkdir %s; tar -xzvf %s.tar.gz -C %s/' % (SCRIPT_PATH, folder_name, archive_basename, archive_basename, archive_basename));
	archive_basename = 'raw_2_rabsch_R7';		execute_command('cd %s/../data/downloads/; cd %s; mkdir %s; tar -xzvf %s.tar.gz -C %s/' % (SCRIPT_PATH, folder_name, archive_basename, archive_basename, archive_basename));
	folder_name = 'salmonella-typhi';
	fastq_file = 'reads-typhi.fastq';
	reads_folder = '%s/../data/consensus-salmonellatyphi/reads/' % (SCRIPT_PATH);
	archive_basename = 'H566_ON_inc'; execute_command('cd %s/../data/downloads/%s; poretools fastq %s/H566_ON_inc/ > %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file));
	archive_basename = 'H566_30_min_inc'; execute_command('cd %s/../data/downloads/%s/%s; find . -name "*.fast5" -exec mv {} ./ \;' % (SCRIPT_PATH, folder_name, archive_basename));
	archive_basename = 'H566_30_min_inc'; execute_command('cd %s/../data/downloads/%s; poretools fastq %s/ >> %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file));
	archive_basename = 'raw_2_rabsch_R7'; execute_command('cd %s/../data/downloads/%s; poretools fastq %s/raw_2_rabsch_R7/downloads/ >> %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file));

	sys.stderr.write('\tFetching the amplicon sequencing dataset.\n');	
	execute_command('cd %s/../data/downloads/; mkdir amplicons' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd amplicons; wget http://files.figshare.com/1866620/R7_3_pgx_metrichor2_23_1D_complement.fastq.bz2' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd amplicons; wget http://files.figshare.com/1866625/R7_3_pgx_metrichor2_23_1D.fastq.bz2' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd amplicons; wget http://files.figshare.com/1866623/R7_3_pgx_metrichor2_23_2D.fastq.bz2' % (SCRIPT_PATH));
	sys.stderr.write('\tExtracting the amplicon sequencing dataset.\n');	
	folder_name = 'amplicons';
	archive_basename = 'R7_3_pgx_metrichor2_23_1D';						execute_command('cd %s/../data/downloads/; cd %s; bzip2 -d --keep %s.fastq.bz2' % (SCRIPT_PATH, folder_name, archive_basename));
	archive_basename = 'R7_3_pgx_metrichor2_23_1D_complement';			execute_command('cd %s/../data/downloads/; cd %s; bzip2 -d --keep %s.fastq.bz2' % (SCRIPT_PATH, folder_name, archive_basename));
	archive_basename = 'R7_3_pgx_metrichor2_23_2D';						execute_command('cd %s/../data/downloads/; cd %s; bzip2 -d --keep %s.fastq.bz2' % (SCRIPT_PATH, folder_name, archive_basename));
	folder_name = 'amplicons';
	fastq_file = 'reads_all-f1000.fastq';
	reads_folder = '%s/../data/amplicons-f1000-1d2d/reads/' % (SCRIPT_PATH);
	archive_basename = 'R7_3_pgx_metrichor2_23_1D'; execute_command('cd %s/../data/downloads/%s; cat %s.fastq > reads_all-f1000.fastq' % (SCRIPT_PATH, folder_name, archive_basename));
	archive_basename = 'R7_3_pgx_metrichor2_23_1D_complement'; execute_command('cd %s/../data/downloads/%s; cat %s.fastq >> reads_all-f1000.fastq' % (SCRIPT_PATH, folder_name, archive_basename));
	archive_basename = 'R7_3_pgx_metrichor2_23_2D'; execute_command('cd %s/../data/downloads/%s; cat %s.fastq >> reads_all-f1000.fastq' % (SCRIPT_PATH, folder_name, archive_basename));
	archive_basename = 'R7_3_pgx_metrichor2_23_2D'; execute_command('cd %s/../data/downloads/%s; cat %s.fastq > reads_2d.fastq' % (SCRIPT_PATH, folder_name, archive_basename));

	sys.stderr.write('\tFetching the BE1 dataset.\n');	
	execute_command('cd %s/../data/downloads/; mkdir be1' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd be1; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA463/ERA463589/oxfordnanopore_native/FAA37759_GB2974_MAP005_20150423__2D_basecalling_v1.14_2D.tar.gz' % (SCRIPT_PATH));
	sys.stderr.write('\tExtracting the BE1 dataset.\n');	
	folder_name = 'be1';
	archive_basename = 'FAA37759_GB2974_MAP005_20150423__2D_basecalling_v1.14_2D';			execute_command('cd %s/../data/downloads/; cd %s; mkdir %s; tar -xzvf %s.tar.gz -C %s/' % (SCRIPT_PATH, folder_name, archive_basename, archive_basename, archive_basename));
	folder_name = 'be1';
	fastq_file = 'reads-BE1.fastq';
	reads_folder = '%s/../data/consensus-BE1/reads/' % (SCRIPT_PATH);
	archive_basename = 'FAA37759_GB2974_MAP005_20150423__2D_basecalling_v1.14_2D'; execute_command('cd %s/../data/downloads/%s; poretools fastq %s/ > %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file));

	sys.stderr.write('\tFetching the ADP1 dataset.\n');	
	execute_command('cd %s/../data/downloads/; mkdir adp1' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd adp1; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA416/ERA416080/fastq/AWK_ONT_MN2064525.fastq.gz' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd adp1; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA416/ERA416080/fastq/AWK_ONT_MN2064006.fastq.gz' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd adp1; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA416/ERA416080/fastq/AWK_ONT_FAA43210.fastq.gz' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd adp1; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA416/ERA416080/fastq/AWK_ONT_FAA43204.fastq.gz' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd adp1; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA416/ERA416080/fastq/AWK_ONT_FAA17573.fastq.gz' % (SCRIPT_PATH));
	sys.stderr.write('\tExtracting the ADP1 dataset.\n');	
	folder_name = 'adp1';
	archive_basename = 'AWK_ONT_MN2064525';			execute_command('cd %s/../data/downloads/; cd %s; gunzip --keep %s.fastq.gz' % (SCRIPT_PATH, folder_name, archive_basename));
	archive_basename = 'AWK_ONT_MN2064006';			execute_command('cd %s/../data/downloads/; cd %s; gunzip --keep %s.fastq.gz' % (SCRIPT_PATH, folder_name, archive_basename));
	archive_basename = 'AWK_ONT_FAA43210';			execute_command('cd %s/../data/downloads/; cd %s; gunzip --keep %s.fastq.gz' % (SCRIPT_PATH, folder_name, archive_basename));
	archive_basename = 'AWK_ONT_FAA43204';			execute_command('cd %s/../data/downloads/; cd %s; gunzip --keep %s.fastq.gz' % (SCRIPT_PATH, folder_name, archive_basename));
	archive_basename = 'AWK_ONT_FAA17573';			execute_command('cd %s/../data/downloads/; cd %s; gunzip --keep %s.fastq.gz' % (SCRIPT_PATH, folder_name, archive_basename));
	folder_name = 'adp1';
	fastq_file = 'all_1d2d.fastq';
	reads_folder = '%s/../data/consensus-ADP1/reads/' % (SCRIPT_PATH);
	archive_basename = 'AWK_ONT_MN2064525'; execute_command('cd %s/../data/downloads/%s; cat %s.fastq > %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file));
	archive_basename = 'AWK_ONT_MN2064006'; execute_command('cd %s/../data/downloads/%s; cat %s.fastq >> %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file));
	archive_basename = 'AWK_ONT_FAA43210'; execute_command('cd %s/../data/downloads/%s; cat %s.fastq >> %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file));
	archive_basename = 'AWK_ONT_FAA43204'; execute_command('cd %s/../data/downloads/%s; cat %s.fastq >> %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file));
	archive_basename = 'AWK_ONT_FAA17573'; execute_command('cd %s/../data/downloads/%s; cat %s.fastq >> %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file));

	sys.stderr.write('\tFetching the Loman et al. dataset used for nanopore-only assembly.\n');
	execute_command('cd %s/../data/downloads/; mkdir ecoli-nmeth' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd ecoli-nmeth; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_20_LomanLabz_PC_Ecoli_K12_R7.3.tar' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd ecoli-nmeth; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_32_LomanLabz_K12_His_tag.tar' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd ecoli-nmeth; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_33_LomanLabz_PC_K12_0.4SPRI_Histag.tar' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd ecoli-nmeth; wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_39.tar' % (SCRIPT_PATH));
	sys.stderr.write('\tExtracting the Loman et al. dataset used for nanopore-only assembly.\n');
	folder_name = 'ecoli-nmeth';
	archive_basename = 'flowcell_20_LomanLabz_PC_Ecoli_K12_R7.3';			execute_command('cd %s/../data/downloads/; cd %s; mkdir %s; tar -xvf %s.tar -C %s/' % (SCRIPT_PATH, folder_name, archive_basename, archive_basename, archive_basename));
	archive_basename = 'flowcell_32_LomanLabz_K12_His_tag';					execute_command('cd %s/../data/downloads/; cd %s; mkdir %s; tar -xvf %s.tar -C %s/' % (SCRIPT_PATH, folder_name, archive_basename, archive_basename, archive_basename));
	archive_basename = 'flowcell_33_LomanLabz_PC_K12_0.4SPRI_Histag';		execute_command('cd %s/../data/downloads/; cd %s; mkdir %s; tar -xvf %s.tar -C %s/' % (SCRIPT_PATH, folder_name, archive_basename, archive_basename, archive_basename));
	archive_basename = 'flowcell_39';										execute_command('cd %s/../data/downloads/; cd %s; mkdir %s; tar -xvf %s.tar -C %s/' % (SCRIPT_PATH, folder_name, archive_basename, archive_basename, archive_basename));
	folder_name = 'ecoli-nmeth';
	fastq_file = 'reads-nmeth-all_2d.fastq';
	reads_folder = '%s/../data/consensus-ecolinmeth/reads/' % (SCRIPT_PATH);
	archive_basename = 'flowcell_20_LomanLabz_PC_Ecoli_K12_R7.3'; execute_command('cd %s/../data/downloads/%s; poretools fastq --type 2D %s/flowcell_20/1.9/downloads/pass/ > %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file));
	archive_basename = 'flowcell_32_LomanLabz_K12_His_tag'; execute_command('cd %s/../data/downloads/%s; poretools fastq --type 2D %s/flowcell_32/downloads/pass/ >> %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file));
	archive_basename = 'flowcell_33_LomanLabz_PC_K12_0.4SPRI_Histag'; execute_command('cd %s/../data/downloads/%s; poretools fastq --type 2D %s/flowcell_33/downloads/pass/ >> %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file));
	archive_basename = 'flowcell_39'; execute_command('cd %s/../data/downloads/%s; poretools fastq --type 2D %s/flowcell_39_K12_Histag/downloads/pass/ >> %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file));

	sys.stderr.write('\tFetching the MAP006-1 dataset.\n');
	execute_command('cd %s/../data/downloads/; mkdir ecoli-map006' % (SCRIPT_PATH));
	execute_command('cd %s/../data/downloads/; cd ecoli-map006; wget http://nanopore.climb-radosgw01.bham.ac.uk/MAP006-1.basecalled.tar' % (SCRIPT_PATH));
	sys.stderr.write('\tExtracting the MAP006-1 dataset.\n');
	folder_name = 'ecoli-map006';
	archive_basename = 'MAP006-1.basecalled';			execute_command('cd %s/../data/downloads/; cd %s; mkdir %s; tar -xvf %s.tar -C %s/' % (SCRIPT_PATH, folder_name, archive_basename, archive_basename, archive_basename));
	folder_name = 'ecoli-map006';
	fastq_file = 'reads-MAP006-1.fastq';
	reads_folder = '%s/../data/consensus-ecoliMAP006/reads/' % (SCRIPT_PATH);
	archive_basename = 'MAP006-1.basecalled'; execute_command('cd %s/../data/downloads/%s; poretools fastq %s/MAP006-1_downloads/pass/ > %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file));
	archive_basename = 'MAP006-1.basecalled'; execute_command('cd %s/../data/downloads/%s; poretools fastq %s/MAP006-1_downloads/fail/ >> %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file));

def get_references(downloads_folder):

### Data that needs to be placed on the repo:
### - indel mutated reference for SV, in data/sv/mutated_reference/escherichia_coli-indel_events.fa
### - amplicons-f1000/truth_variants
### - 

	# ###
	# folder_name = 'R7';
	# reads_folder = '%s/../data/consensus-ecoliR7/reads' % (SCRIPT_PATH);
	# fastq_file_2d = 'Ecoli_R7_CombinedFasta-2d.fastq';
	# execute_command('mkdir -p %s/../data/downloads/%s' % (SCRIPT_PATH, folder_name));
	# archive_basename = 'Ecoli_R7_NONI'; execute_command('cd %s/../data/downloads/%s; poretools fastq --type 2D %s/downloads/ > %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file_2d));
	# archive_basename = 'Ecoli_R7_ONI_flowcell_18'; execute_command('cd %s/../data/downloads/%s; poretools fastq --type 2D %s/ >> %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file_2d));
	# archive_basename = 'Ecoli_R7_ONI_flowcell_17'; execute_command('cd %s/../data/downloads/%s; poretools fastq --type 2D %s/ >> %s'  % (SCRIPT_PATH, folder_name, archive_basename, fastq_file_2d));
	# execute_command('mkdir -p %s; cp %s/../data/downloads/%s/%s %s/'  % (reads_folder, SCRIPT_PATH, folder_name, fastq_file_2d, reads_folder));
	###
	# folder_name = 'R7.3';
	# reads_folder = '%s/../data/consensus-ecoliR7.3-2d/reads/' % (SCRIPT_PATH);
	# fastq_file_2d = 'ecoliR7.3-2d.fastq';
	# execute_command('mkdir -p %s/../data/downloads/%s' % (SCRIPT_PATH, folder_name));
	# archive_basename = 'Ecoli_R73'; execute_command('cd %s/../data/downloads/%s; poretools fastq --type fwd,rev %s/downloads/ > %s' % (SCRIPT_PATH, folder_name, archive_basename, fastq_file_2d));
	# execute_command('mkdir -p %s; cp %s/../data/downloads/%s/%s %s/'  % (reads_folder, SCRIPT_PATH, folder_name, fastq_file_2d, reads_folder));

	# folder_name = 'all_2d_for_sv';
	# reads_folder = '%s/../data/sv/reads/' % (SCRIPT_PATH);
	# fastq_file_2d = 'all_2d_for_sv.fastq';
	# execute_command('mkdir -p %s/../data/downloads/%s' % (SCRIPT_PATH, folder_name));
	# execute_command('cd %s/../data/downloads/%s; cat %s/../data/consensus-ecoliR7/reads/Ecoli_R7_CombinedFasta-2d.fastq > %s'  % (SCRIPT_PATH, folder_name, SCRIPT_PATH, fastq_file_2d));
	# execute_command('cd %s/../data/downloads/%s; cat %s/../data/consensus-ecoliR7.3-2d/reads/ecoliR7.3-2d.fastq >> %s'  % (SCRIPT_PATH, folder_name, SCRIPT_PATH, fastq_file_2d));
	# execute_command('mkdir -p %s; cp %s/../data/downloads/%s/%s %s/'  % (reads_folder, SCRIPT_PATH, folder_name, fastq_file_2d, reads_folder));



	### Downloading the E. Coli gi|48994873|gb|U00096.2|
	reference_dest = "escherichia_coli.fa";
	reference_url = "http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&sendto=on&log$=seqview&db=nuccore&dopt=fasta&sort=&val=50083297&from=begin&to=end&extrafeat=976&fmt_mask=294912";
	execute_command('mkdir -p %s/../data/downloads/references/; cd %s/../data/downloads/references/; wget %s -o %s'  % (SCRIPT_PATH, SCRIPT_PATH, reference_url, reference_dest));
	# ###
	# reference_folder = '%s/../data/reference/' % (SCRIPT_PATH);
	# execute_command('mkdir -p %s; cp %s/../data/downloads/references/%s %s'  % (reference_folder, SCRIPT_PATH, reference_dest, reference_folder));
	# ###
	# reference_folder = '%s/../data/consensus-ecoliR7/reference/' % (SCRIPT_PATH);
	# execute_command('mkdir -p %s; cp %s/../data/downloads/references/%s %s'  % (reference_folder, SCRIPT_PATH, reference_dest, reference_folder));
	# ###
	# reference_folder = '%s/../data/consensus-ecoliR7.3/reference/' % (SCRIPT_PATH);
	# execute_command('mkdir -p %s; cp %s/../data/downloads/references/%s %s'  % (reference_folder, SCRIPT_PATH, reference_dest, reference_folder));
	# ###
	# reference_folder = '%s/../data/nmeth/reference/' % (SCRIPT_PATH);
	# execute_command('mkdir -p %s; cp %s/../data/downloads/references/%s %s'  % (reference_folder, SCRIPT_PATH, reference_dest, reference_folder));
	# ###
	# reference_folder = '%s/../data/sv/original_reference/' % (SCRIPT_PATH);
	# execute_command('mkdir -p %s; cp %s/../data/downloads/references/%s %s'  % (reference_folder, SCRIPT_PATH, reference_dest, reference_folder));
	# ###
	# reference_folder = '%s/../data/consensus-ecoliMAP006/reference/' % (SCRIPT_PATH);
	# execute_command('mkdir -p %s; cp %s/../data/downloads/references/%s %s'  % (reference_folder, SCRIPT_PATH, reference_dest, reference_folder));
	# ###
	# reference_folder = '%s/../data/consensus-ecoliR7.3-1d/reference/' % (SCRIPT_PATH);
	# execute_command('mkdir -p %s; cp %s/../data/downloads/references/%s %s'  % (reference_folder, SCRIPT_PATH, reference_dest, reference_folder));
	# ###
	# reference_folder = '%s/../data/consensus-ecoliR7.3-2d/reference/' % (SCRIPT_PATH);
	# execute_command('mkdir -p %s; cp %s/../data/downloads/references/%s %s'  % (reference_folder, SCRIPT_PATH, reference_dest, reference_folder));



	# reference_folder = '%s/../data/consensus-lambdaR6/reference/' % (SCRIPT_PATH);
	# reference_dest = "NC_001416.fa";
	# reference_url = "http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&sendto=on&log$=seqview&db=nuccore&dopt=fasta&sort=&val=9626243&from=begin&to=end&extrafeat=976";
	# execute_command('mkdir -p %s/../data/downloads/references/; cd %s/../data/downloads/references/; wget %s -o %s'  % (SCRIPT_PATH, SCRIPT_PATH, reference_url, reference_dest));
	# execute_command('mkdir -p %s; cp %s/../data/downloads/references/%s %s'  % (reference_folder, SCRIPT_PATH, reference_dest, reference_folder));

	# reference_folder = '%s/../data/consensus-uti89/reference/' % (SCRIPT_PATH);
	# reference_dest = "NC_007946.1.fa";
	# reference_url = "http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&sendto=on&log$=seqview&db=nuccore&dopt=fasta&sort=&val=91209055&from=begin&to=end&extrafeat=976&fmt_mask=295336";
	# execute_command('mkdir -p %s/../data/downloads/references/; cd %s/../data/downloads/references/; wget %s -o %s'  % (SCRIPT_PATH, SCRIPT_PATH, reference_url, reference_dest));
	# execute_command('mkdir -p %s; cp %s/../data/downloads/references/%s %s'  % (reference_folder, SCRIPT_PATH, reference_dest, reference_folder));

	# reference_folder = '%s/../data/consensus-salmonellatyphi/reference/' % (SCRIPT_PATH);
	# reference_dest = "salmonella_typhi_Ty2.fa";
	# reference_url = "http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&sendto=on&log$=seqview&db=nuccore&dopt=fasta&sort=&val=29140543&from=begin&to=end&extrafeat=976&fmt_mask=294912";
	# execute_command('mkdir -p %s/../data/downloads/references/; cd %s/../data/downloads/references/; wget %s -o %s'  % (SCRIPT_PATH, SCRIPT_PATH, reference_url, reference_dest));
	# execute_command('mkdir -p %s; cp %s/../data/downloads/references/%s %s'  % (reference_folder, SCRIPT_PATH, reference_dest, reference_folder));

	# reference_dest = "ref_chr6_chr22-GRCh37.p13.fa";
	# chr6_url = "ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37.p13/Primary_Assembly/assembled_chromosomes/FASTA/chr6.fa.gz";
	# execute_command('mkdir -p %s/../data/downloads/references/; cd %s/../data/downloads/references/; wget %s; gunzip --keep chr6.fa.gz;'  % (SCRIPT_PATH, SCRIPT_PATH, chr6_url));
	# chr22_url = "ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37.p13/Primary_Assembly/assembled_chromosomes/FASTA/chr22.fa.gz";
	# execute_command('mkdir -p %s/../data/downloads/references/; cd %s/../data/downloads/references/; wget %s; gunzip --keep chr22.fa.gz;'  % (SCRIPT_PATH, SCRIPT_PATH, chr22_url));
	# execute_command('cd %s/../data/downloads/references/; cat chr6.fa > %s; cat chr22.fa >> %s'  % (SCRIPT_PATH, reference_dest, reference_dest));
	# ###
	# reference_folder = '%s/../data/amplicons-f1000/reference/' % (SCRIPT_PATH);
	# execute_command('mkdir -p %s; cp %s/../data/downloads/references/%s %s'  % (reference_folder, SCRIPT_PATH, reference_dest, reference_folder));
	# ###
	# reference_folder = '%s/../data/amplicons-f1000-1d2d/reference/' % (SCRIPT_PATH);
	# execute_command('mkdir -p %s; cp %s/../data/downloads/references/%s %s'  % (reference_folder, SCRIPT_PATH, reference_dest, reference_folder));

	# # reference_folder = '%s/../data/' % (SCRIPT_PATH);
	# # reference_dest = "";
	# # reference_url = "";
	# # execute_command('mkdir -p %s/../data/downloads/references/; cd %s/../data/downloads/references/; wget %s -o %s'  % (SCRIPT_PATH, SCRIPT_PATH, reference_url, reference_dest));
	# # execute_command('mkdir -p %s; cp %s/../data/downloads/references/%s %s'  % (reference_folder, SCRIPT_PATH, reference_dest, reference_folder));




	# reference_folder = '%s/../data/consensus-ADP1/reference/' % (SCRIPT_PATH);
	# reference_dest = "acinetobacter_baylyi_ADP1_NC_005966.1.fa";
	# reference_url = "http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&sendto=on&log$=seqview&db=nuccore&dopt=fasta&sort=&val=50083297&from=begin&to=end&extrafeat=976&fmt_mask=294912";
	# execute_command('mkdir -p %s/../data/downloads/references/; cd %s/../data/downloads/references/; wget %s -o %s'  % (SCRIPT_PATH, SCRIPT_PATH, reference_url, reference_dest));
	# execute_command('mkdir -p %s; cp %s/../data/downloads/references/%s %s'  % (reference_folder, SCRIPT_PATH, reference_dest, reference_folder));

	# reference_folder = '%s/../data/consensus-BE1/reference/' % (SCRIPT_PATH);
	# reference_dest = "bacteroides_fragilis-BFBE1.fa";
	# reference_url = "http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&sendto=on&log$=seqview&db=nuccore&dopt=fasta&sort=&val=928804415&from=begin&to=end&extrafeat=976&fmt_mask=294912";
	# execute_command('mkdir -p %s/../data/downloads/references/; cd %s/../data/downloads/references/; wget %s -o %s'  % (SCRIPT_PATH, SCRIPT_PATH, reference_url, reference_dest));
	# execute_command('mkdir -p %s; cp %s/../data/downloads/references/%s %s'  % (reference_folder, SCRIPT_PATH, reference_dest, reference_folder));

	return;

























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

if __name__ == "__main__":
	downloads_folder = '%s/../data/downloads' % (SCRIPT_PATH);

#	get_real_data();
#	move_data_to_proper_folders(downloads_folder);
	get_references(downloads_folder);
