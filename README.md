# CRISPR-PairGuide_Design work flow in terminal

python Create_Alignment_files.py chr# start end                                                                             
python Sguide_calculation.py chr# start end                   
python Filter_specificity_threshold.py          
python Percent_chromatin_free.py               
python Guide_pairing.py           
python Oligonucleotide_generator.py           


#only need to run this once will create 214 negative controls for hg37        
python Negative_controls_hg37.py        
