import pandas as pd


# Converts a list of probes to coordinates in bed format and outputs to file
def probe2bed(probes, file):
    id2coord = pd.read_table("../metadata/id2coord.txt", header=0, index_col=0)
    id2coord = id2coord.iloc[[i for i in range(len(id2coord)) if is_number(id2coord.CHR[i])]]
    id2coord = id2coord.astype(int)
    probesCoord = id2coord.iloc[[i for i in range(len(id2coord)) if id2coord.index[i] in probes]]
    chr = ["chr" + str(probesCoord['CHR'][i]) for i in range(len(probesCoord))]
    start = probesCoord['MAPINFO'].tolist()
    end = (probesCoord['MAPINFO'] + 1).tolist()
    bed = pd.DataFrame(dict(chr=chr, coord1=start, coord2=end, probe=probesCoord.index))
    bed = bed.sort_values(['chr', 'coord1'])
    bed.to_csv(file, index=False, header=False, sep='\t')
    
