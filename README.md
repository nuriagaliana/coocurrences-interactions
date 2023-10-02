Power-laws in species' biotic interaction networks can be inferred from co-occurrence data

This repository contains all raw data and code used in the above manuscript from authors

Nuria Galiana, Jean-François Arnoldi, Frederico Mestre, Alejandro Rozenfeld, Miguel B. Araujo.

Date Created: 05-08-2022s

Email: galiana.nuria@gmail.com

In this README file we describe the contents of the archive (i.e. files and folders contained here) 
and provide a brief description of how the source code scripts should be executed. 

DATA FILES:
The archive is organised into a series of folders each containing the raw data and source code appropriate
for analysing that data set. The folder have been named according to the names of each data set. Thus, for example, the folder title 'Salix? contains the raw data for the raw data in the Salix dataset as well as the script used to process those data in order to construct networks of realised interactions and co-occurrence interactions. Within each folder named after the corresponding dataset, a subfolder named 'raw-data' contains the actual data. Please notice that the format of each raw dataset is different, and therefore, the code needed to read each dataset is also different. 
Specification of data format per dataset:
Garraf PP, Garraf PP2, Garraf HP, Olot, Montseny: raw data formatted in spreadsheets where each sheet has the information of the network of biotic interactions of each spatial unit of the dataset.
Gottin PP: raw data in the file ?plant_poll_all_interactions.csv?.Spatial units are defined in the column ?Site?; plant species are defined in the column ?P.Genus.Species?; pollinator species are defined in the column ?Genus.Species?.
Gottin HP: raw data in the file ?host_para_all_interactions.csv?. Spatial units are defined in the column ?Site?; host species are defined in column ?Genus.Species?; parasitoid species are defined in column ?P1.Species.Genus?.
Nahuel: raw data formatted as .txt files where each file describes the matrix of biotic interactions of each spatial unit.
Quercus: raw data formatted as spreadsheet file where each file describes the matrix of biotic interactions of each spatial unit.
Galpar: the format of this data is explained in detail in: Kopelke, J. P. et al. Food-web structure of willow-galling sawflies and their natural enemies across Europe. Ecology 98, 1730 (2017). 

SOURCE CODE FILES:
1.- 'data_processing_occurrence.r?: Within each dataset folder, a script called 'data_processing_occurrencer' is provided. 
This script contains the code necessary to read and process the raw data for each dataset and generate the realised and the co-occurrence networks. 
2.- 'data-analysis_occurrence.r?: This script, located within the parent folder, reads the output files from each of the 
datasets and produces the figures and performs the statistical analyses presented in the paper. 
3.- 'utils.r': An additional script containing a series of food web property measurement functions is provided
to aid with the calculations of network properties across datasets. This script is sourced from the 
'data_processing_occurrence.r? files above.


Additionally, we provide .csv files with the outputs of the already processed data in case you want to directly produce 
the figures or perform the analyses without having to generate the realised and the co-occurrence networks. Therefore, in each data folder you will find, for instance, ?network_real_pred_garraf-PP.csv? for the realised networks and ?network_coocurrence_pred_garraf-PP.csv? for the co-occurrence networks, both from the predators perspective. These .csv files are read in the ?data-analysis_occurrence.r? source code to perform the analyses.


If you have any questions / run into any issues, please contact galiana.nuria@gmail.com