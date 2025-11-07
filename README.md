# 251030_SeahorseXLSX_to_ATPprod
Python script to convert seahorse assay output OCR and ECAR measures into ATP production rates

written by Brennan Wadsworth

Script or app will ask you to select a file, choose the desired seahorse output xlsx file, a desired directory for saving the results, and specification on which injection the oligomycin and Rotenone/Antimycin A was given.

The script requires these injections to be present and requires the "Assay Configuration", "Rate", and "Raw" sheets from the excel file. 

The script will scan for wells labelled "acid", "buffer", or H2SO4 to find wells to run the buffering power calculation. If it does so it will output a data file and graph to confirm the buffering power calculations were run. Otherwise it will assume a standard value (that calculated for Brennan's run on 31/10/2025).

The results file will be a CSV with the calculated values of OXPHOS and glycolytic ATP production rate in pmol ATP per min for each well.
