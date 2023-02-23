from comparison_process import *
from util import readfile, writefile
import os




def compare(ref_mtg_name,tested_mtg_name,cacheprefix = None):
    print('read mtgs')
    ref_mtg = readfile(ref_mtg_name)
    tested_mtg = readfile(tested_mtg_name)
    # if cacheprefix and not os.path.exists(cacheprefix) : os.makedirs(cacheprefix)
    return comparison_process(ref_mtg,tested_mtg,cacheprefix)

def create_args(prefix) : return (prefix+'_validated_piperadius.bmtg',prefix+'_reconstruction_piperadius.bmtg',prefix)

def main():
    # compare('mais_validated_piperadius.bmtg','mais_reconstruction_piperadius.bmtg','puu4')
    #newcscore, tscore = compare(*create_args('puu1'))
    newcscore, tscore = compare('puu1_validated.bmtg','puu1_reconstruction.bmtg') #,'puu1')
    print(newcscore, tscore)
    
if __name__ == '__main__':
    main()