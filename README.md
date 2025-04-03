# proab2
1.Download AbRSA from the official repository (http://aligncdr.labshare.cn/aligncdr/download.html) and extract the archive.
2.Modify the ABRSA_PATH variable in the auto_generate_config-final.py script to point to your AbRSA installation directory.
3.Obtain antibody-protein complex structures from the chai-1 server and place the prediction files (pred.rank_0.cif through pred.rank_4.cif) in the execution directory.
4.Run auto_generate_config-final.py to generate the configuration file input_config_cif.txt.
5. python fix-integrated-gbsa-analyzer-2.py -c input_config_cif.txt -o results -s 2000000 --stride 100
This procedure will perform molecular dynamics simulation analysis with 2,000,000 steps and sampling every 100 steps, outputting the results to the "results" directory.
