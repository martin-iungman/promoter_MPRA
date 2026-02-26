## From FASTQ to count tables  

Download FASTQ data and place it in Data_processing/FASTQ.  
Run Data_processing/scripts/complete_analysis.sh from its location. It will call sequentially to the other scripts in the folder for each step of the processing of the data.  
Data_processing/scripts/data_processing.qmd contains the whole code and some explanations, along with code for a few plots to visualize the process.  

## External data

External_data/download_data.sh will download most of the external data used. External_data/download_fantom5.sh is specific for FANTOM5 necessary files.
