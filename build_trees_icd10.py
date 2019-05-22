""""""

import sys
import copy
import cPickle as pickle
import numpy as np
import pandas as pd

# ===== Function Definition =====

def store(obj, path):
    """Stores an object with the given path"""
    with open(path, "wb") as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


# ===== Execution flow =====

if __name__ == '__main__':
    infile = sys.argv[1]
    seqFile = sys.argv[2]
    typeFile = sys.argv[3]
    outFile = sys.argv[4]

    infd = open(infile, 'r')
    _ = infd.readline()

    seqs = pickle.load(open(seqFile, 'rb'))
    types = pickle.load(open(typeFile, 'rb'))

    startSet = set(types.keys())
    hitList = []
    missList = []
    # (n level) counts
    cat1count = 0
    cat2count = 0
    
    for line in infd:
        tokens = line.strip().split(',')
        icd9 = tokens[0][1:-1].strip()
        cat1 = tokens[1][1:-1].strip()
        desc1 = 'A_' + tokens[2][1:-1].strip()
        cat2 = tokens[3][1:-1].strip()
        desc2 = 'A_' + tokens[4][1:-1].strip()
        
        if icd9.startswith('E'):
            if len(icd9) > 4: icd9 = icd9[:4] + '.' + icd9[4:]
        else:
            if len(icd9) > 3: icd9 = icd9[:3] + '.' + icd9[3:]
        icd9 = 'D_' + icd9

        if icd9 not in types: 
            missList.append(icd9)
        else: 
            hitList.append(icd9)

        if desc1 not in types: 
            cat1count += 1
            types[desc1] = len(types)

        if len(cat2) > 0:
            if desc2 not in types: 
                cat2count += 1
                types[desc2] = len(types)

    infd.close()

    rootCode = len(types)
    types['A_ROOT'] = rootCode
    print rootCode

    set_hitlist = set(hitList)
    set_misslist = startSet - set_hitlist

    print("""
        ------
        Cat1 codes' count: {}
        Cat2 codes' count: {}
        All codes' count: {}
        Hitlist: {}
        Misslist: {}
        ------""".format(cat1count,
                         cat2count,
                         cat1count + cat2count + 1,
                         len(set_hitlist),
                         len(set_misslist)))

    # (n levels + 1) maps
    # map1 holds missed / singleton codes
    # map2 holds level-1 codes ('root' codes)
    # map3 holds level-2 codes
    maps = {
        "map3": {},
        "map2": {},
        "map1": { types[icd]: ([types[icd], rootCode]) for icd in set_misslist }
    }

    ccs_tree = pd.read_csv(infile)
    print(ccs_tree.columns)

    infd = open(infile, 'r')
    infd.readline()

    with open(outFile + ".log", "a") as f:
        for line in infd:
            tokens = line.strip().split(',')
            icd9 = tokens[0][1:-1].strip()
            cat1 = tokens[1][1:-1].strip()
            desc1 = 'A_' + tokens[2][1:-1].strip()
            cat2 = tokens[3][1:-1].strip()
            desc2 = 'A_' + tokens[4][1:-1].strip()

            if icd9.startswith('E'):
                if len(icd9) > 4: icd9 = icd9[:4] + '.' + icd9[4:]
            else:
                if len(icd9) > 3: icd9 = icd9[:3] + '.' + icd9[3:]
            icd9 = 'D_' + icd9


            if icd9 not in types:
                # print("Comparing original Diagonsis codes with CCS codes")
                # print("If not included in Diagnosis, it's not useful to populate the tree with them.")
                f.write("ICD not in Diagnosis types. Level {}; included in set_misslist: True\n".format(
                    len(cat2.split('.'))
                ))
                # jumps the loop to the next element
                continue
            
            icdCode = types[icd9]

            codeVec = []
            

            if len(cat2) > 0:
                code2 = types[desc2]
                code1 = types[desc1]
                maps["map3"][icdCode] = [icdCode, rootCode, code1, code2]
            elif len(cat2) <= 0:
                print "Code1: " + code1
                code1 = types[desc1]
                maps["map2"][icdCode] = [icdCode, rootCode, code1]
    
    print("map 2 values length: " + str(len(maps["map2"].values())))

    # Now we re-map the integers to all medical codes.
    # (n levels + 1) re-maps
    # store maps in a dict
    remap = {
        "level3": {},
        "level2": {},
        "level1": {}
    }

    code_types = {}
    rtypes = { v: k for k, v in types.iteritems() }

    # print("rtypes", rtypes)

    codeCount = 0

    for nm, vs in maps.iteritems():
        for icd_code, ancestors in vs.iteritems():
            # do a re-count of all code types no matter its level
            code_types[rtypes[icd_code]] = codeCount
            
            # re-map all codes in each level by its read order
            # the first value of each key is the key itself
            # also drop the first ancestor (the value itself)
            remap["level" + nm[-1]][codeCount] = [codeCount] + ancestors[1:]
            
            # update the count
            codeCount += 1

    # re-map the patient visits by its code
    seqs_by_code = []
    for patient in seqs:
        new_patient = [
            [code_types[rtypes[code]] for code in visit]
            for visit in patient
        ]
        seqs_by_code.append(new_patient)

    # Store objects in a pickle way
    """
    for nm, dt in remap.iteritems():
        # Print state before storing
        print(nm, len(dt.values()))
        store(dt, "{}.{}.pk".format(outFile, nm))
    """
    
    # As the code parsing doesn't find any rootCode belonging to 
    # both input and ccs tree, overwrite this level
    # remap["level2"] = None
    store(remap["level3"], "{}.{}.pk".format(outFile, "level2"))
    store(remap["level1"], "{}.{}.pk".format(outFile, "level1"))

    store(seqs_by_code, outFile + ".seqs")
    store(code_types, outFile + ".types")
