# Scripts to reproduce the results from the GraphMap paper.

All tests on real data are implemented in a single script ```run.py```.  
Tests are implemented as individual functions that are called sequentially one after the other.  
Simulations are implemented in our other repository ```aligneval``` which will be downloaded and setup automatically.  


### Installation
The tools that are required for testing can be set-up using the following command:  
```  
run.py setup-tools  
```
This will automatically download and install: ```PyVCF```, ```Tabix```, ```VCFtools```, ```LoFreq```, ```samscripts```, ```Mutatrix``` and ```Bamsurgeon```.  
Sudo is required to install the first three on the list.  

The above command will also install ```aligneval``` which is our repo for evaluating alignments on simulated data.  
```Aligneval``` includes the wrapper scripts for each mapper/aligner required to run the tests. It will automatically download and install all aligners, but will not start generating simulated data right away. For this, the ```setup-simdata``` is used.  

**Manual setup** required for GATK and PicardTools because of their download/licencing process. Please obtain ```GenomeAnalysisTK-3.4-46``` and ```picard-tools-1.138``` and place them in the ```tools``` folder.  
Make sure that these are accessible: ```graphmap-reproduce/tools/picard-tools-1.138/picard.jar``` and ```graphmap-reproduce/tools/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar```.  

### Running tests on simulations  
```  
src/run-simulated.sh  
```  
Everything should setup and run automatically.  
This script will clone and install ```Aligneval```, unpack the references, unpack the pre-simulated data, generate some simulated data which wasn't packed, and download the hg19 GRCh37 reference.  
The entire process might take some time, especially if you consider the slowness of some mappers.  
The user will be prompted several times during the installation of ```Aligneval```, and twice at the very beginning of the alignment process. The prompts at the beginning of alignment are to test/skip the slower mappers.  
Final results will be placed in the ```graphmap-reproduce/results``` folder.  

### Running tests on real data
This is a bit more trickier, because the real data is huge. Most of the data is provided as raw nanopore ```.fast5``` files, and unpacking them might take over ```800GB``` of disk space.  
The entire process may take more than a day on a computer with a good Internet connection.  
To initiate this process, run:  
```  
src/setup-realdata.py  
```  

**Warning** - at this point, automatic setup of reference sequence has yet to be implemented. One can download them manually from NCBI and place them in the correct paths (specified in ```src/run-realdata.py```).  

To run the alignment and the evaluation of the results, run:  
```  
src/run-realdata.py  
```  
Please make sure that all reads and references are downloaded and correctly placed.  

Sudo is required to run the tests, as the memory/time measurements require access to the root folders.  
