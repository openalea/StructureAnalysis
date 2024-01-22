# Read .seq file and transform to an R file.
import os

def Setup(chdir=False, virtual=True):
    """
    Default setup and paths to data
    """

    print(os.getcwd())

    if not(virtual):
        home_dir = "D:\\"
    else:
        home_dir = "/media/sf_transfert"
    base_path = os.path.join(home_dir, "devlp_shared", "AppleCultivarsTreatments")
    data_path = base_path + os.sep + "Data"
    
    if chdir:
        os.chdir(base_path)

    output_file = base_path + os.sep + "Data" + os.sep + "multiseq_N_APPLE_R.csv"
    
    return base_path, data_path, output_file

def ReadSequence(data_path="", data_file="multiseq_N_APPLE_unix.seq"):
    """
    Read sequence from file and return Sequences object + python list
    """
    if data_path == "":
        base_path, data_path, output_file = Setup(chdir=False)
    
    # Python sequence 
    from openalea.sequence_analysis import sequences
    pyseq1 = sequences.Sequences(data_path + os.sep + data_file)
    pylist = []
    dim_seq = len(pyseq1[0][0])

    for i in range(len(pyseq1)):
        seq_i = []
        for l in range(len(pyseq1[i])):
            seq_i += [[pyseq1[i][l][d] for d in range(dim_seq)]]
        pylist += [seq_i]
        
    return pyseq1, pylist

def DefaultHeader():
    """
    Return default header for .csv
    """
    h = "######################################################################### \n"
    h += "# \n"
    h += "#  N-APPLE shoot \n"
    h += "# \n"
    h += "#  VARIABLE 1: treatment code (1: control, 2: treatment with 20 g N/tree dose, 3: treatment with 30 g N/h tree dose), \n"
    h += "#  VARIABLE 2: cultivar (1: Rubinola, 2: Topaz, 3: Golden Delicious), \n"
    h += "#  VARIABLE 3: axillary shoot type (0: no shoot, 1: short shoot, 2: medium shoot, 3: long shoot, 4: sylleptic shoot), \n"
    h += "#  VARIABLE 4: lateral flowering (0: no lateral flowering, 1: LF in median zone, 2: LF in distal zone, 3: LF in median and distal zone = Long LF zone), \n"
    h += "#  VARIABLE 5: terminal flowering (0: no flower cluster, 1: presence of a flower cluster). \n"
    h += "# \n"
    h += "#  Data : Martin Meszaros, Evelyne Costes, Jean-Baptiste Durand \n"
    h += "# \n"
    h += '#  read with read.csv(multiseq_N_APPLE_R.csv, header=TRUE, row.names=1, comment.char="#") \n'
    h += "# \n"
    h += "######################################################################### \n"
    h += '"", treatment code, '
    h += "cultivar, "
    h += "axillary shoot type, " 
    h += "lateral flowering, "
    h += "terminal flowering, "
    h += "sequence N. \n"    
    return h
    
def RestoredStatesHeader():
    """
    Return header for restored states
    """
    h = "######################################################################### \n"
    h += "# \n"
    h += "#  N-APPLE shoot \n"
    h += "# \n"
    h += "#  VARIABLE 1: restored state. \n"
    h += "#  VARIABLE 2: axillary shoot type (0: no shoot, 1: short shoot, 2: medium shoot, 3: long shoot, 4: sylleptic shoot), \n"
    h += "#  VARIABLE 3: lateral flowering (0: no lateral flowering, 1: LF in median zone, 2: LF in distal zone, 3: LF in median and distal zone = Long LF zone), \n"
    h += "#  VARIABLE 4: terminal flowering (0: no flower cluster, 1: presence of a flower cluster). \n"
    h += "# \n"
    h += "#  Data : Martin Meszaros, Evelyne Costes, Jean-Baptiste Durand \n"
    h += "# \n"
    h += '#  read with read.csv(multiseq_N_APPLE_R.csv, header=TRUE, row.names=1, comment.char="#") \n'
    h += "# \n"
    h += "######################################################################### \n"
    h += '"", restored state, '
    h += "axillary shoot type, " 
    h += "lateral flowering, "
    h += "terminal flowering, "
    h += "sequence N. \n"    
    return h
    
def WriteRSequence(pyseq1, output_file="", header=None):
    """
    Write sequence in a .txt file, one line per sequence index
    """
    if output_file == "":
        base_path, data_path, output_file = Setup(chdir=False)

    f = open(output_file, "w")
    
    if header is None:
        h = DefaultHeader()
    else:
        h = str(header)
    f.write(h)
    dim_seq = len(pyseq1[0][0])
    l = 1 # line
    for i in range(len(pyseq1)):
        # i: sequence
        for p in range(len(pyseq1[i])):
            # p: position
            f.write('"' + str(l) + '", ')
            for v in range(dim_seq):
                # write variables
                f.write(str(pyseq1[i][p][v]) + ", ")
            # write sequence id
            f.write(str(i) + "\n")
            l += 1
        
    f.close()
