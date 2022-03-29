import argparse
import subprocess
import hicstraw
import pandas as pd
from tools import run_command

def parseargs():
    parser = argparse.ArgumentParser(description='Download and dump HiC data')
    parser.add_argument('--hic_file', required=True, help="Path or url to .hic file.")
    parser.add_argument('--resolution', default=5000, help="Resolution of HiC to download. In units of bp.")
    parser.add_argument('--outdir', default=".")
    parser.add_argument('--include_raw', action="store_true", help="Download raw matrix in addtion to KR")
    parser.add_argument('--chromosomes', default="all", help="comma delimited list of chromosomes to download")
    parser.add_argument('--skip_gzip', action="store_true", help="dont gzip hic files")

    return parser.parse_args()

def save_to_dataframe(result, output):
    binX=[]
    binY=[]
    counts=[]
    for i in range(len(result)):
        binX.append(result[i].binX)
        binY.append(result[i].binY)
        counts.append(result[i].counts)
    df = pd.DataFrame()
    df[0]= binX
    df[1]= binY
    df[2]= counts
    df.to_csv("{}".format(output), sep="\t", index=False, header=False)

def main(args):

    if args.chromosomes == "all":
        chromosomes = list(range(1,23)) + ['X']
    else:
        chromosomes = args.chromosomes.split(",")

    for chromosome in chromosomes:
        print("Starting chr" + str(chromosome) + " ... ")
        outdir = "{0}/chr{1}/".format(args.outdir, chromosome)
        command = "mkdir -p " + outdir
        out = subprocess.getoutput(command)
        ## Download observed matrix with KR normalization
        result = hicstraw.straw('observed', 'SCALE', '{}'.format(args.hic_file), 'chr{}'.format(chromosome), 'chr{}'.format(chromosome), 'BP', int(args.resolution))
        save_to_dataframe(result, "{}/chr{}.INTERSCALEobserved.gz".format(outdir, chromosome))

        ## Download KR norm file
        result = hicstraw.straw('norm', 'SCALE', '{}'.format(args.hic_file), 'chr{}'.format(chromosome), 'chr{}'.format(chromosome), 'BP', int(args.resolution))
        binX=[]
        binY=[]
        counts=[]
        for i in range(len(result)):
            binX.append(result[i].binX)
            binY.append(result[i].binY)
            counts.append(result[i].counts)
        df = pd.DataFrame()
        df[0] = counts
        df.to_csv("{}/chr{}.INTERSCALEnorm.gz".format(outdir, chromosome), sep="\t", index=False, header=False)
        
        if args.include_raw:
            ## Download raw observed matrix
            result = hicstraw.straw('observed', 'NONE', '{}'.format(args.hic_file), 'chr{}'.format(chromosome), 'chr{}'.format(chromosome), 'BP', int(args.resolution))
            save_to_dataframe(result, "{}/chr{}.RAWobserved.gz".format(outdir, chromosome))
        

if __name__ == '__main__':
    args = parseargs()
    main(args)
