# Export python dictionaries of sufficients HMSC statistics for 
# GLMMs to an R file.

import os


def ComputeGLMStatistics(pyseq1, seg_pyseq):
    """
    Compute statistics for GLM(M) estimation
    """
    dwell = {} # indexed by state, cultivar, tree
    transition = {} # indexed by past_state, treatment, cultivar, tree
    observation = {} # indexed by state, treatment, cultivar, tree
    for i in range(len(pyseq1)):
        d = 0
        s_past = -1 # past state

        for l in range(len(pyseq1[i])):
            va = pyseq1[i][l]
            t = va[0] # treatment
            c = va[1] # cultivar
            v = seg_pyseq[i][l]
            s = v[0]
            if s != s_past:
                if s_past != -1:
                    if dwell.has_key((s_past, t, c, i+1)):
                        dwell[s_past, t, c, i+1] += [[d, 1]]
                        transition[s_past, t, c, i+1] += [s]
                    else:
                        dwell[s_past, t, c, i+1] = [[d, 1]]
                        transition[s_past, t, c, i+1] = [s]
                d = 1 # dwell time 
            else:
                d += 1
            s_past = s
            if observation.has_key((s, t, c, i+1)):
                observation[s, t, c, i+1] += [v[1:]]
            else:
                observation[s, t, c, i+1] = [v[1:]]
        if dwell.has_key((s_past, t, c, i+1)):
            dwell[s_past, t, c, i+1] += [[d, 0]]
        else:
            dwell[s_past, t, c, i+1] = [[d, 0]]

    return dwell, transition, observation

def GLMStatisticsHeader():

    """
    Return header for GLM statistics
    """
    h = "######################################################################### \n"
    h += "# \n"
    h += "#  Statistics for GLM(M) estimation issued from HSMC restoration\n"
    h += "# \n"
    h += '#  VARIABLE 1: type of statistics ("d"well time, "t"ransition, "o"bservation). \n'
    h += '#  VARIABLE 2: restored state (current state for "d" and "o", previous state for "t"). \n'
    h += "#  VARIABLE 3: treatment. \n"
    h += "#  VARIABLE 4: cultivar. \n"
    h += "#  VARIABLE 5: value of statistics. \n"
    h += '#  VARIABLE 6: for "d", 0 if last dwell time (censored), 1 otherwise; for "t", NA \n'
    h += '#              for "o", id of variable. \n'
    h += "# \n"
    h += '#  read with read.csv(file_name.csv, header=TRUE, row.names=1, comment.char="#") \n'
    h += "# \n"
    h += "######################################################################### \n"
    h += '"", type, '
    h += "restored state, " 
    h += "treatment, "
    h += "cultivar, "
    h += "sequence N., "
    h += "value, "
    h += "info \n"
    return h
    
    
def WriteGLMStatistics(dwell, transition, observation, output_file):
    """
    Write statistics for GLM(M) estimation
    """
    f = open(output_file, "w")
    f.write(GLMStatisticsHeader())
    l = 1 # line 
    for (k, v) in dwell.items():
        f.write(str(l) + ", ")
        f.write("d, ")
        for o in k:
            f.write(str(o) + ", ")
        s = ""
        for o in v:
            for i in o:
                s += str(i) + ", "
        s = s[0:-2] + "\n"
        l += 1
        f.write(s)

    for (k, v) in transition.items():
        f.write(str(l) + ", ")
        f.write("t, ")
        for o in k:
            f.write(str(o) + ", ")
        s = ""
        for o in v:
            s += str(o) + ", "
        s = s[0:-2] + ", NA \n"
        l += 1
        f.write(s)

    for (k, v) in observation.items():
        # common part of all elements in current 
        # list of objects and associated string
        s_base = "o, " 
        for o in k:
            s_base += str(o) + ", "
        for o in v:
            for i in range(len(o)):
                # Variables
                s = s_base + str(o[i]) + ", " + str(i+1) + "\n"
                f.write(str(l) + ", " + s)
                l += 1

    f.close()

