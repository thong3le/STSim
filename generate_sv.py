import random
import sys

class Node:
	def __init__(self, breakpoints, label=None, left=None, right=None):
		self.breakpoints = breakpoints
		self.label = label
		self.left = left
		self.right = right

def traverse(root, sv_dict, path = []):
	"""
		root: Node  -> the root of the tree
		sv_dict: dict -> map a leave label to a list of breakpoints
		path: list 	-> the path contains breakpoints so far
	"""
	if root:
		if root.left is None and root.right is None:
			path += root.breakpoints
			sv_dict[root.label] = path
			print "Leave ", root.label, ": ", len(path)
		else:
			traverse(root.left, sv_dict, path + root.breakpoints)
			traverse(root.right, sv_dict, path + root.breakpoints)

def reference(filename):
	"""
		filename: the name of reference genome fasta file
		only read the chr1
	"""
	chr1 = []
	ref = open(filename)
	for line in ref:
		if not line.startswith('>'):
			chr1.append(line.strip())
		elif 'chr1' not in line:
			return chr1
	return chr1

def generate_fa(ref, breakpoints):
	"""
		ref: list of String  -> ref is a list of strings of A,T,G,C
		SV: list of breakpoints
	"""
	size = len(ref)
	chosen1 = [True] * size
	chosen2 = [True] * size
	for l,r,t,label in breakpoints:
		for p in range(l,r):
			chosen1[p] = False
		if t == 0:
			for p in range(l,r):
				chosen2[p] = False
	newref1 = []
	newref2 = []
	for i in range(size):
		if chosen1[i]:
			newref1.append(ref[i])
		if chosen2[i]:
			newref2.append(ref[i])
	return newref1, newref2

def outfile_fa(filename, ref1, ref2):
	"""
		filename: String  -> the output fasta filename
		ref: list of String
	"""
	f = open(filename, 'w')
	f.write('>chr1_1\n')
	for line in ref1:
		f.write(line+'\n')
	f.write('>chr1_2\n')
	for line in ref2:
		f.write(line+'\n')
	f.close()

def outfile_pos(filename, breakpoints, ref_line_length):
	"""
		filename: String -> name of the output file for positions of SVs
		SV: list of breakpoints
		ref_line_segment: int -> number character per line in ref_genome fasta files
	"""
	breakpoints.sort()
	L = ref_line_length
	f = open(filename, 'w')
	for l,r,t,label in breakpoints:
		TYPE = "homozygous" if t == 0 else "heterozygous"
		f.write('chr1\t' + str(l*L)+'\t'+str(r*L)+'\t'+str(L*(r-l))+'\t'+TYPE+'\t'+label+'\n')
	f.close()

def readTree(filename):
	"""
		filename: String -> name of the input_tree file which contains the encoding of the somatic tree
	"""
	tree = {}
	f = open(filename, 'r')
	num_sv = 0
	for line in f:
		if line.startswith('#') or len(line.strip()) == 0:
			continue
		label, edges = line.split(':')
		if label.strip() == 'R':
			if len(edges.split()) != 2:
				print("Incorrect format: The first line should be 'R: x n'")
				sys.exit()
			rootlabel, w = edges.split()
			tree['R'] = [rootlabel, int(w)]
			num_sv += int(w)
		else:
			if len(edges.split()) != 4:
				print("Incorrect format: each internal node should have two children 'Bi: x1 n1 x2 n2'")
				sys.exit()
			label1, w1, label2, w2 = edges.split() 
			tree[label] = [label1, int(w1), label2, int(w2)]
			num_sv += int(w1) + int(w2)
	return tree, num_sv

def getTree(filename, ref_size):
	"""
		filename: String -> name of the input_tree file which contains the encoding of the somatic tree
		ref_size: int -> how many line in ref_genome fasta file (for chr1)
	"""
	deletion_min_size = 20  # 20 lines in FASTA file
	deletion_max_size = 200  # 100 lines in FASTA file
	tree, num_sv = readTree(filename)
	print "Total number of SVs is ", num_sv
	deletions = [random.randint(deletion_min_size, deletion_max_size) for _ in range(num_sv)]
	SV = []  # (l,r,type) where (l,r) is the position of SV and if type = 0, it is homozygous, else if type = 1 it is heterozygous
	#random.shuffle(deletions)
	bin = ref_size/num_sv
	offset = 10
	for i in range(num_sv):
		l = random.randint(i*bin+offset, (i+1)*bin - deletion_max_size - offset)
		r = l + deletions[i]
		t = 0 if random.random() < 0.5 else 1
		SV.append((l,r,t))

	random.shuffle(SV)

	if 'R' not in tree:
		print "Missing the root node."
		sys.exit()

	sv_index = 0
	label, w = tree['R']
	breakpoints = [(l,r,t,label) for l,r,t in SV[sv_index:sv_index+w]]
	root = Node(breakpoints, label)
	sv_index += w
	stack = [root]
	while stack:
		node = stack.pop()
		if node.label in tree:
			label1, w1, label2, w2 = tree[node.label]
			breakpoints = [(l,r,t,label1) for l,r,t in SV[sv_index:sv_index+w1]]
			node.left = Node(breakpoints, label1)
			sv_index += w1
			breakpoints = [(l,r,t,label2) for l,r,t in SV[sv_index:sv_index+w2]]
			node.right = Node(breakpoints, label2)
			sv_index += w2
			stack.append(node.left)
			stack.append(node.right)
	return root
	
def main():
	if len(sys.argv) != 3 and len(sys.argv) != 4:
		print sys.argv
		print "Usage: python generate_sv.py ref_genome input_tree [output_dir]"
		return
	else:
		ref_file = sys.argv[1]
		input_tree = sys.argv[2]
		output_prefix = "" if len(sys.argv) == 3 else sys.argv[3].strip('/') + '/'

	ref = reference(ref_file)
	#print len(ref)
	ref_size = len(ref)
	ref_line_length = len(ref[0])
	root = getTree(input_tree, ref_size)
	sv_dict = {}
	traverse(root, sv_dict)
	for leave in sv_dict:
		name = "subclone_at_leave_" + leave
		newref1, newref2 = generate_fa(ref, sv_dict[leave])
		outfile_fa(output_prefix + name + ".fa", newref1, newref2)
		outfile_pos(output_prefix + name + "_pos.txt", sv_dict[leave], ref_line_length)
	#print sv_dict


if __name__ == '__main__':
	main()
