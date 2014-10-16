#!/usr/bin/env python
#
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import csv,sys,argparse,os

def read_dat(fname,cn):
   csv.register_dialect('tab', delimiter='	', quoting=csv.QUOTE_NONE)
   if os.path.isdir(fname): fname = fname + '/dist.dat'
   with open(fname, 'rb') as f:
      reader = csv.reader(f, 'tab')
      for row in reader:
          if row[0][0]=='#': continue
          type = row[0]
          chr  = row[1]
          if type=='CN':
              cn[chr] = float(row[2])

def plot_copy_number(plot_fname,dat,titles):
    xlabs = {}
    for cn in dat:
        xlabs.update(cn)
    xlabels = sorted(xlabs.keys())
    symbols=['v','^','<','>']
    colors=['red','green','blue','orange','magenta','cyan']
    ymax = 0
    ymin = 1e9 
    fig, ax = plt.subplots(1, 1, figsize=(7,5))
    for j in range(len(dat)):
        cn = dat[j]
        title = titles[j]
        symbol = symbols[j % len(symbols)]
        color = colors[j % len(colors)]
        yvals = []
        xvals = []
        for i in range(len(xlabels)):
            if xlabels[i] in cn:
                if cn[xlabels[i]]==-1:
                    ax.annotate('?', xy=(i+j*0.05,0.5+j*0.01),va='center',ha='center',color=color,fontweight='bold')
                else:
                    if ymin>cn[xlabels[i]]: ymin = cn[xlabels[i]]
                xvals.append(i)
                yvals.append(cn[xlabels[i]])
                if ymax<cn[xlabels[i]]: ymax = cn[xlabels[i]]
        ax.plot(xvals,yvals,symbol,color=color,mec=color,label=title)
    ax.legend(loc='best',prop={'size':7},frameon=False,numpoints=1)
    ax.tick_params(axis='both', which='major', labelsize=9)
    ax.set_xticks(range(len(xlabels)))
    ax.set_xticklabels(xlabels,rotation=45)
    ax.set_xlim(-1,len(xlabels))
    ax.set_ylim(0,5.0)
    ax.set_yticks([1.0,2.0,3.0,4.0])
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Copy Number')
    plt.savefig(plot_fname)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output-file', help='Ouput file name')
    parser.add_argument('files',metavar='title:dist.dat', nargs='+')
    args = parser.parse_args()
    if args.output_file==None: parser.print_help(); sys.exit(1);
    dat = []
    titles = []
    for arg in args.files:
        cn = {}
        arg = arg.split('@')
        if len(arg)==1: arg = (arg,arg)
        read_dat(arg[1],cn)
        dat.append(cn)
        titles.append(arg[0])
    plot_copy_number(args.output_file,dat,titles)

if __name__ == '__main__':
   main()
