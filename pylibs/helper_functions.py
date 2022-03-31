import os
import psutil
import sys
import subprocess
from .processing_classes import RunLogger

def call_Command(cmd, logger_, shell_=False):
	'''
	Runs a shell command using subprocess.run. Recods output to a logger object that has already been created.
	Default is to use subprocess.run. If shell is necessary then subprocess.Popen is used.
	'''
	try:
		if shell_ == False:
			capture = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True).stdout
			logger_.logger.info(capture)
		else:
			capture = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			logger_.logger.info(capture.communicate()[0].decode('utf-8'))
	except:
		logger_.logger.exception('\nFluPipeline Error:', exc_info=True)


def automatic_ThreadUsage(max_mem_per_thread, safety_factor=1, total_memory=psutil.virtual_memory()[0]/(1024.0**3), cpu_count=os.cpu_count()):
	'''
	Returns the number of threads to use based on the estimated memory usage per thread
	safety_factor: reduce the number of threads requested by this amount to allow for additional memory overflow.
	'''
	if max_mem_per_thread < 0:
		raise ValueError('--max_mem_per_thread must be a positive value')

	if max_mem_per_thread != None:
		threads_to_use = (total_memory/max_mem_per_thread) # calculate number of threads to use

		if str(threads_to_use).find('.')>-1:
			threads_to_use = int(str(threads_to_use)[:str(threads_to_use).find('.')]) # round down the calculated number of threads to use

		if threads_to_use == 0:
			raise ValueError('Not enough total memory availble for the estimated memory per thread --\nTotal memory available: {}Gb\nEstimated memory per thread: {}Gb'.format(total_memory,  max_mem_per_thread))

		if threads_to_use > cpu_count:
			threads_to_use = cpu_count # ensure number of threads is at least equal to cpu count

		if threads_to_use > safety_factor:
			threads_to_use =  threads_to_use - safety_factor # reduce number of threads by safety_factor to allow memory overflow

		if threads_to_use == cpu_count & cpu_count != 1:
			threads_to_use =  threads_to_use - safety_factor # if the cpu_count is equal to threads to use and there's more than 1, reduce threads calculated by safety_factor

		if cpu_count != 1:
			if (threads_to_use + safety_factor) * max_mem_per_thread < total_memory - max_mem_per_thread: 
				threads_to_use = threads_to_use + safety_factor # if the total memory used from all threads available is lower than the total memory minus one safety_factor, use all threads available

	return(threads_to_use)
