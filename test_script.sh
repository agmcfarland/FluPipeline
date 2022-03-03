
conda activate testenv

cd /home/agmcfarland/flu_project/FluPipeline

python FluPipeline.py -h

python FluPipeline.py \
--runtest --threads 6 --strain_sample_depth 3000

python FluPipeline.py \
--base_directory /home/agmcfarland/flu_project/shared_data/test_1_sample2 \
--sequence_directory /home/agmcfarland/flu_project/shared_data/test_data_1_sample \
--force \
--force_base_directory \
--threads 6 \
--cleanup 



python FluPipeline.py \
--base_directory /change/this \
--sequence_directory /home/agmcfarland/flu_project/shared_data/test_data_6_samples \
--force \
--force_base_directory \
--threads 6 \
--cleanup 


