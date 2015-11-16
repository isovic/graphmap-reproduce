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

### Data setup
To simulate datasets, run the script using:  
```run.py setup-simdata```  
This can take some time, but the process is performed only once.

To setup the real datasets, run the script using:  
```run.py setup-realdata```  
This will download and unpack **a lot** of data, ~700GB of disk space will be required.
Most of the public datasetsa are available only in the raw format, which consumes much more space than the final FASTQ files.
```setup-realdata``` will also setup the corresponding data folders which will be used to run the tests.  
The entire process may take more than a day on a computer with a good Internet connection.  

***Please note*** that ```setup-realdata``` is still not fully finished and **is under active development**, but it should be finished within a few days. All data links are given in the ```setup_data()``` function within the script, the slow part is to test the consistency of the setup commands because of the size of the data. We are currently working on this part, and should have the final setup function ready shortly.  
In the meantime, the user can set the data up manually by following the paths described within the script for each separate test.  



### Usage
Running the script with no parameters will list out the usage of the script to stderr.  
In short, to run all the tests on real datasets, use:  
```  
run.py run-realdata
```  
To run the simulations, use the following:  
```  
run.py run-simdata  
```  

Sudo is required to run the tests, as the memory/time measurements require access to the root folder.  


Automatic setup of datasets is still in development.  
Otherwise, most of the tests can be run if the data folder is set manually, by following the path specifications in each of the tests.  
