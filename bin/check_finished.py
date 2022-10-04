

import os
from os.path import join as pjoin
import sys
import glob

runDir = sys.argv[1]

# runDir = '/data/Pipeline_output/061722_Output'
# runDir = '/data/flu_project/benchmarking_project/compare_caller/small/lofreq'

sampleLogsDir = pjoin(runDir,'sampleLogs')

store_results = {'finished':[],'unfinished':[],'error':[]}
for f in glob.glob(sampleLogsDir+'/*.txt'):
	if f.find('_ivar') > -1:
		continue
	else:
		with open(f,'r') as infile:
			finished = False
			error = False
			for l in infile:
				if l.find('FluPipeline Time:') > -1:
					finished = True
				if l.find('FluPipeline Error') > -1:
					error = True
			if finished == True:
				store_results['finished'].append(os.path.basename(f).replace('.txt',''))
			if error == True:
				store_results['error'].append(os.path.basename(f).replace('.txt',''))
			if finished == False:
				store_results['unfinished'].append(os.path.basename(f).replace('.txt',''))

print('\n')
print('Major variant stage results:')
for k, v in store_results.items():
	print('-'*10 + k + '-'*10)
	for i in v:
		print(i)
	print('Total:',len(v),'\n')




store_results = {'finished':[],'unfinished':[],'error':[]}
for f in glob.glob(sampleLogsDir+'/*.txt'):
	if f.find('_ivar') > -1:
		with open(f,'r') as infile:
			finished = False
			error = False
			for l in infile:
				if l.find('Finished running') > -1:
					finished = True
				if l.find('FluPipeline Error') > -1:
					error = True
			if finished == True:
				store_results['finished'].append(os.path.basename(f).replace('.txt',''))
			if error == True:
				store_results['error'].append(os.path.basename(f).replace('.txt',''))
			if finished == False:
				store_results['unfinished'].append(os.path.basename(f).replace('.txt',''))


print('Intrahost variant stage results:')
for k, v in store_results.items():
	print('-'*10 + k + '-'*10)
	for i in v:
		print(i)
	print('\n')
	print('Total:',len(v),'\n')





