// g++ -Wall -O3 -o ~/bin/stats stats.cpp fastq.cpp fasta.cpp mh12_utils.cpp sam.cpp cigar.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <vector>

#include "fasta.h"
#include "fastq.h"
#include "mh12_utils.h"

using namespace std;


struct Stats
{
    double mean;
    unsigned long n50[9];
    unsigned long n50n[9];
    unsigned long longest;
    unsigned long shortest;
    unsigned long number;
    unsigned long totalLength;
    unsigned long nCount;
};


struct CmdLineOps
{
    unsigned long minLength;
    int infileStartIndex;
    bool ignoreNs;
    string outType;
};

void parseOptions(int argc, char** argv, CmdLineOps& ops);
Stats file2stats(string filename, CmdLineOps& ops);

void print_long(string fname, Stats& s);
void print_short(string fname, Stats& s);

int main(int argc, char* argv[])
{
    CmdLineOps ops;
    parseOptions(argc, argv, ops);
    bool first = true;

    for (int i = ops.infileStartIndex; i < argc; i++)
    {
        Stats s = file2stats(argv[i], ops);
        if (ops.outType.compare("long") == 0)
        {
            print_long(argv[i], s);
        }
        else
        {
            if (!first)
            {
                cout << "-------------------------------" << endl;
            }
            print_short(argv[i], s);
            first = false;
        }

    }

    return 0;
}



void parseOptions(int argc, char** argv, CmdLineOps& ops)
{
    string usage;
    ops.minLength = 1;
    ops.infileStartIndex = 1;
    ops.ignoreNs = false;
    ops.outType = "short";

    usage = "usage: stats [options] list of fasta/q files\n\n\
options:\n\
-l <int>\n\tMinimum length cutoff for each sequence [1]\n\
-n\n\tIgnore Ns in the sequences\n\
-s\n\tPrint shorter stats (this is the default, but\n\
\toption still here for backwards compatibility)\n\
-S\n\tPrint longer stats\n";

    if (argc < 2)
    {
        cerr << usage;
        exit(1);
    }

    while (argv[ops.infileStartIndex][0] == '-')
    {
        if (strcmp(argv[ops.infileStartIndex], "-s") == 0)
        {
            ops.outType = "short";
            ops.infileStartIndex++;
        }
        else if (strcmp(argv[ops.infileStartIndex], "-S") == 0)
        {
            ops.outType = "long";
            ops.infileStartIndex++;
        }
        else if (strcmp(argv[ops.infileStartIndex], "-l") == 0)
        {
            ops.minLength = atoi(argv[ops.infileStartIndex + 1]);
            ops.infileStartIndex += 2;
        }
        else if (strcmp(argv[ops.infileStartIndex], "-n") == 0)
        {
            ops.ignoreNs = true;
            ops.infileStartIndex++;
        }
        else
        {
            cerr << "error parsing options, somewhere around this: " << argv[ops.infileStartIndex] << endl;
            exit(1);
        }
    }
}


Stats file2stats(string filename, CmdLineOps& ops)
{
    Stats s;
    vector<int> seqLengths;
    ifstream inStream;
    int cumulativeLength = 0;
    short filetype = fastaOrFastq(filename);

    inStream.open(filename.c_str());

    if (! inStream.is_open())
    {
        cerr << "Error opening file " << filename << endl;
        exit(1);
    }

    s.totalLength = 0;
    s.nCount = 0;

    while(inStream.good())
    {
        Fasta* seq;

        if (filetype == mh12::FASTQ_FILE)
        {
            Fastq* p;
            p = new Fastq;
            seq = (Fastq*) p;
        }
        else if (filetype == mh12::FASTA_FILE)
        {
            seq = new Fasta();
        }
        else
        {
            cerr << "Input file type not recognised as fasta or fastq.  Aborting" << endl;
            exit(1);
        }

        if ( !(seq->fillFromFile(inStream)) )
        {
            break;
        }
        else if (seq->length() >= ops.minLength)
        {
            int l;
            int n = seq->nCount();

            if (ops.ignoreNs)
            {
                l = seq->length() - n;
            }
            else
            {
                l = seq->length();
            }

            seqLengths.push_back(l);
            s.totalLength += l;
            s.nCount += n;
        }

        delete seq;
    }

    inStream.close();

    for (int i = 0; i < 9; i++)
    {
        s.n50[i] = 0;
        s.n50n[i] = 0;
    }


    if (seqLengths.size() == 0)
    {
        s.longest = 0;
        s.shortest = 0;
        s.number = 0;
        s.mean = 0;
        s.totalLength = 0;
        return s;
    }

    sort(seqLengths.begin(), seqLengths.end());
    s.longest = seqLengths.back();
    s.shortest = seqLengths.front();
    s.number = seqLengths.size();
    s.mean = 1.0 * s.totalLength / s.number;


    int k = 0;

    for (int i = seqLengths.size() - 1; 0 <= i; i--)
    {
        cumulativeLength += seqLengths[i];

        while (k < 9 && cumulativeLength >= s.totalLength * (k + 1) / 10.0)
        {
            s.n50[k] = seqLengths[i];
            s.n50n[k] = seqLengths.size() - i;
            k++;
        }
    }

    return s;
}


void print_long(string fname, Stats& s)
{
    cout.precision(2);

    cout << fname << " " << "total_length " << s.totalLength << endl
         << fname << " " << "number       " << s.number << endl
         << fname << " " << "mean_length  " << fixed << s.mean << endl
         << fname << " " << "longest      " << s.longest << endl
         << fname << " " << "shortest     " << s.shortest << endl
         << fname << " " << "N_count      " << s.nCount << endl;

    for (int j = 0; j < 9; j++)
    {
        cout << fname << " " << "n" << j + 1 << "0  " << s.n50[j] << endl
             << fname << " " << "n" << j + 1 << "0n " << s.n50n[j] << endl;
    }
}

void print_short(string fname, Stats& s)
{
    cout.precision(2);

    cout << "stats for " << fname << endl
         << "sum = " << s.totalLength
         << ", n = " << s.number
         << ", ave = " << fixed << s.mean
         << ", largest = " << s.longest << endl
         << "N50 = " << s.n50[4] << ", n = " << s.n50n[4] << endl
         << "N60 = " << s.n50[5] << ", n = " << s.n50n[5] << endl
         << "N70 = " << s.n50[6] << ", n = " << s.n50n[6] << endl
         << "N80 = " << s.n50[7] << ", n = " << s.n50n[7] << endl
         << "N90 = " << s.n50[8] << ", n = " << s.n50n[8] << endl
         << "N100 = " << s.shortest << ", n = " << s.number << endl
         << "N_count = " << s.nCount << endl;
}
