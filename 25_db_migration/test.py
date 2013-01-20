#!/usr/bin/python

from alignment import * # C implementation of smith waterman

#########################################
def main():
    

    human_seq      =  "B001--MEEPQSDPSVEPPL-SQETFSDLWKZB002-LPE-NNVL---Z"
    pep_seq_to_fix =  "B001VITDERDPGTSPPVPTGHSFATPLALLLAGTWFSPISLZ"
    #human_seq      =  "--MEEPQSDPSVEPPL-SQETFSDLWKZ-LPE-NNVL---"
    #pep_seq_to_fix =  "VITDERDPGTSPPVPTGHSFATPLALLLAGTWFSPISL"
    [aligned_human, aligned_ortho] \
        = smith_waterman_context (human_seq, pep_seq_to_fix)

    print human_seq
    print pep_seq_to_fix
    print
    print ">human"
    print aligned_human
    print ">ortho"
    print aligned_ortho

    return True


#########################################
if __name__ == '__main__':
    main()
