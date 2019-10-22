from __future__ import division
import re
import sys
from pyfaidx import Fasta
import getopt
import threading
import time
from queue import Queue

# 参数默认值
base_len = 1
base_num = 1
n_compute_threads=5
opts, args = getopt.getopt(sys.argv[1:], "hi:o:r:m:Y:X:Z:t:",["help", "input=", "output=", "screening position=", "screening length=", "the least number of screening changes=", "ref=", "filter mode=","thread="])

for opt_name, opt_value in opts:
	if opt_name in ('-i', '--input'):
		infile = opt_value
	if opt_name in ('-o', '--output'):
		outfile = opt_value
	if opt_name in ('-r', '--ref'):
		ref = opt_value  # 参考基因组，.fa格式
	if opt_name in ('-m', '--filter mode'):
		mode = opt_value  # 用来选择是做脱嘌呤还是脱氨基的过滤（depurination、deamination）
	if opt_name in ('-Z', '--screening position'):
		alter = opt_value  # 用来选择是5/3end两端都要符合还是任意一段（and、or）
	if opt_name in ('-Y', '--screening length'):
		base_len = int(opt_value)  # 用来选择脱氨基特征过滤的时候，前后长度read中存在脱氨基特征
	if opt_name in ('-X', '--the least number of screening changes'):
		base_num = int(opt_value)  # 用来选择脱氨基特征过滤的时候，前后特定长长度read中至少存在多少个脱氨基特征
	if opt_name in ('-t','--thread'):
		n_compute_threads=int(opt_value)

	if opt_name in ('-h', '--help'):
		print("required parameters: -i、-o、-m、-r")
		print("optional parameters: -h、-X、-Y、-Z、-t")
		print(
			"-i,--input:    input file(sam format)\n"
			"-o,--output:   output file(sam format)\n"
			"-r,--ref:	reference file(.fa format) to get the mapped seq in the reference\n"
			"-m,--mode:	mode(depurination/deamination),based on choosed mode to filter the comtamination\n"
			"-Z,--screening position:    alternative(and/or),\"and\" means that the program will filter the comtamination when the end of 3 and 5 both fit the conditions\n"
			"-Y,--screening length: it means that we set up condition that in X bp in 5/3 end there are some damage patterns(for --mode deamination)\n"
			"-X,--the least number of screening changes: it means that we set up condition that in fixed bp in 5/3 end there are X damage patterns(for --mode deamination)\n"
			"-h,--help: "
		)
##############
genes = Fasta(ref)
# Same with n_compute_threads
READ_BATCH_SIZE = 25
WRITE_BATCH_SIZE = 25
print("The program is using %d threads" % n_compute_threads)
# read_queue
r_q = Queue(READ_BATCH_SIZE)
# write_queue
w_q = Queue(WRITE_BATCH_SIZE)



class ProducerThread(threading.Thread):

	def __init__(self, infile, batch_size=10000):
		super(ProducerThread, self).__init__()
		self.infile = infile
		self.batch_size = batch_size

	def run(self):

		# compute batch_size
		num = 0
		lines = []
		f = open(self.infile, "r")

		for line in f:
			line = line.strip()
			if re.match('^@', line):
				continue
			global amount_nofilter
			elements = line.strip().split("\t")
			chr = elements[2]
			pos = int(elements[3])
			fragment = elements[9]

			if elements[1] == "4":  # flag=4 代表这个序列没有mapping到参考序列上
				continue
			if int(elements[4]) < 30:  # 跳过比对质量值小于30的reads
				continue
			if re.findall(r'^chrM|M', elements[2]):  # 线粒体是环状的，提取的代码比较麻烦，考虑环状。考虑到线粒体基因组小所以忽略
				continue
			if re.findall(r'[IDNSHF]', elements[5]):
				continue
			num += 1
			lines.append(line)
			if num % self.batch_size == 0:
				# Block until r_q is not full
				r_q.put(lines)
				print('Read %d,Start time: %s' % (num, time.ctime(int(time.time()))))
				sys.stdout.flush()
				lines = []

				time.sleep(0.1)
		f.close()

		if len(lines) > 0:
			r_q.put(lines)
			print('Read %d,Start time: %s' % (num, time.ctime(int(time.time()))))
			sys.stdout.flush()
			lines = []

		# Terminate Signal
		# for i in range(10):
		for i in range(READ_BATCH_SIZE):
			r_q.put(None)


class ComputeConsumerThread(threading.Thread):
	def __init__(self, genes):
		super(ComputeConsumerThread, self).__init__()
		self.genes = genes

	def run(self):
		while True:
			# Block until get
			queue_item = r_q.get()
			if queue_item == None:
				w_q.put(None)
				break
			lines = queue_item
			lines_filtered = []
			if mode=="depurination":
				for line in lines:
					if judge_depurination(line, tiqu_depurination(line)):
						lines_filtered.append(line)
				w_q.put(lines_filtered)
				time.sleep(0.1)
			elif mode=="deamination":
				for line in lines:
					elements = line.split("\t")
					fragment = elements[9]
					MD_tag = ""
					for element in elements:
						MD_tags = re.findall("^MD:Z:(.*)", element)
						if MD_tags:
							MD_tag = MD_tags[0]
					ref_seq = tiqu_deamination(MD_tag, fragment)
					if judge_deamination(fragment, ref_seq, base_len, base_num):
						lines_filtered.append(line)
				w_q.put(lines_filtered)
				time.sleep(0.1)


class WriteConsumerThread(threading.Thread):

	def __init__(self, outfile, n_compute_threads):
		super(WriteConsumerThread, self).__init__()
		self.outfile = outfile
		self.n_compute_threads = n_compute_threads

	def run(self):
		num = 0
		flush_count = 0

		fw = open(self.outfile, "w")
		while True:
			out_lines = w_q.get()

			if out_lines == None:
				num += 1
				if num >= int(self.n_compute_threads):
					break
				continue

			else:
				fw.write("\n".join(out_lines) + "\n")
				flush_count += 1
				if flush_count % 100 == 0:
					fw.fileno()
				#print('Write,Start time: %s' % time.ctime(int(time.time())))
				sys.stdout.flush()
				time.sleep(0.1)

		fw.flush()
		fw.close()


##############
def tiqu_depurination(read):
	try:
		genes = Fasta(ref)
		elements = re.split("\t", read)
		chr = elements[2]
		pos = int(elements[3])
		length = len(elements[9])
		# 提取reads
		ref_seq = genes[chr][pos - 1:pos + length - 1].seq
		if pos == 1:
			base_fommer = "NONE"
			base_latter = genes[chr][pos + length - 1:pos + length].seq
		elif length > len(ref_seq) or pos + length - 1 == len(genes[chr]):
			base_fommer = genes[chr][pos - 2:pos - 1].seq
			base_latter = "NONE"
		else:
			base_fommer = genes[chr][pos - 2:pos - 1].seq
			base_latter = genes[chr][pos + length - 1:pos + length].seq
		# 返回数组[比对上的参考序列，参考前一个base，参考后一个base]
		seq = [ref_seq, base_fommer, base_latter]
	except Exception as e:
		sys.stderr.write(str(e) + "\n")
		sys.stderr.write(read + "\n###############\n" + ref_seq)
	return seq


def judge_depurination(read, seq):
	# 基于脱嘌呤特征的判断
		# only analysis depurination in 5' end
		if re.match(r'[AGag]', seq[1]):
			return True
		else:
			return False

def tiqu_deamination(MD,fragment):
	alts=re.split("\d+",MD)
	nums=re.split("[a-zA-z]+",MD)
	seq=""
	count=0
	if alts[0]!="":
		for i in range(len(alts)):
			try:
				seq+=(alts[i]+fragment[count:count+int(nums[i+1])])
				count+=len(fragment[count:count+int(nums[i])]+alts[i+1])
			except IndexError:
				seq+=alts[i]

	if nums[0]!="":
		for i in range(len(nums)):
			try:
				seq+=(fragment[count:count+int(nums[i])]+alts[i+1])
				count+=len(fragment[count:count+int(nums[i])]+alts[i+1])
			except IndexError:
				seq+=fragment[count:count+int(nums[i])]
	return seq

def judge_deamination(raw_seq,ref_seq,base_len,base_num):
	count = 0 #计算有多少个脱氨基特征在这个长度的reads中
	five_end = three_end = False
	for base1, base2 in zip(raw_seq[:base_len], ref_seq[:base_len]):
		if re.findall('[Cc]', base2) and re.findall('[Tt]', base1):
			count += 1
	if count >= base_num:
		five_end = True
	count = 0
	for base1, base2 in zip(raw_seq[-base_len:], ref_seq[-base_len:]):
		if re.findall('[Gg]', base2) and re.findall('[Aa]', base1):
			count += 1
	if count >= base_num:
		three_end = True
	if alter == "and":
		if three_end and five_end:
			return True
	elif alter == "or":
		if three_end or five_end:
			return True


def thread_main(infile, outfile, genes, batch_size, n_compute_threads=5):
	p = ProducerThread(infile, batch_size)
	threads = []
	p.start()
	threads.append(p)

	for i in range(n_compute_threads):
		compute_c = ComputeConsumerThread(genes)
		compute_c.start()
		threads.append(compute_c)

	for i in range(1):
		write_c = WriteConsumerThread(outfile, n_compute_threads)
		write_c.start()
		threads.append(write_c)

	start_time = time.time()
	print("Compute Start: %s" % time.ctime(int(start_time)))
	for thread in threads:
		thread.join()

	end_time = time.time()
	print("Compute End: %s" % time.ctime(int(end_time)))

	last_time = end_time - start_time

	print("Last time: %.2f mins" % (last_time / 60))




if __name__ == '__main__':
	batch_size = 500000
	thread_main(infile, outfile, genes, batch_size, n_compute_threads)

sys.stderr.write("Program finished successfully\n")










