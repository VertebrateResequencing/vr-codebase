#!/software/bin/python
# Aylwyn Scally 2008

import subprocess
import sys
import getopt
import re
import os

def usage():
	print 'usage: gcdepth.py depth_file [-i indiv_name] [-m max_depth] [-b data_bin | -s samp_freq] [-n n_gcbins] [-l Rfile] [-R path/to/R] [--sim] [--keepscript]'
	sys.exit(2)

def message(str):
	sys.stderr.write('%s: %s\n' % (os.path.basename(sys.argv[0]), str))


# defaults
label = ''
depfile = ''
sim = False
keepscript = False
rpath = '/software/R-2.6.0/bin/R'
sfile = 'gcdepth-script.R'
Ndatmax = 6e6
nbins = 30
#readlen = 35
Rlib = '/lustre/scratch1/dmc/g1k/bin/gcdepth.R'
#Rlib = '/software/solexa/bin/gcdepth/gcdepth.R'
datbin = 0
depmax = 'NULL'
indiv = ''
binned  = False
samp  = False

scrname = sys.argv[0]
try:
	opts, args = getopt.gnu_getopt(sys.argv[1:], 's:m:i:n:b:r:l:R:', ['sim', 'keepscript'])
except getopt.GetoptError:
	usage()

for (oflag, oarg) in opts:
	if oflag == '--sim':
		sim = True
	if oflag == '--keepscript':
		keepscript = True
	if oflag == '-m':
		depmax = oarg
	if oflag == '-r':
		reflink = oarg
	if oflag == '-l':
		Rib = oarg
	if oflag == '-i':
		indiv = oarg
	if oflag == '-R':
		rpath = oarg
	if oflag == '-b':
		datbin = int(oarg)
		binned = True
	if oflag == '-s':
		datbin = int(oarg)
		samp = True
	if oflag == '-n':
		nbins = int(oarg)

if len(args) > 0:
	depfile = args[0]
else:
	usage()

if not sim:
	fsf = open(sfile, 'w')
else:
	fsf = sys.stdout

lnum = 0
for line in open(depfile):
	lnum += 1
	if lnum == 1:
		pos0 = int(line.split()[1])
	if lnum == 2:
		datbin = int(line.split()[1]) - pos0
	if lnum == Ndatmax:
		break
Ndat = lnum

if not binned and not samp:
	bmatch = re.search(r'bindepth', depfile)
	if bmatch:
		binned = True
if not binned:
	samp = True
	message('Assuming sampled depth data')
	binstr = 'FALSE'
else:
	message('Assuming binned depth data')
	binstr = 'TRUE'
##depmax = int(1.5 * os.stat('aln.map').st_size / 3.5e9 * datbin / readlen)

if datbin == 0:
	message('error: Unable to determine data binsize from depth file')
	sys.exit(2)
else:
	message('Setting binsize = %d' % datbin)

fsf.write('source(\'%s\')\n' % Rlib)

fsf.write('depdat = read.depdat(\'%s\', Ndat = %d, bin = %d)\n' % (depfile, Ndat, datbin))
fsf.write('gcdepth(depdat, sname = \'%s\', depmax = %s, hc = TRUE, nbins = %d, binned = %s)\n' % (indiv, depmax, nbins, binstr))

if not sim:
	fsf.close()

cmd = '%s --no-save --slave < %s' % (rpath, sfile)
if sim:
	print(cmd)
else:
	subprocess.Popen(cmd, shell = True).wait()

if not sim and not keepscript:
	os.remove(sfile)
