import numpy as np

def set_xor_2d(matrice1,matrice2):
    
    # if len(matrice1[0,:]) != len(matrice2[0,:]) or len(matrice1[0,:]) != 2:
        # print "ERROR: Different number of columns! in 'setxor2d'"
        # return
    
    list1 = []
    for i in range(len(matrice1[:,0])):
        list1.append((matrice1[i,0],matrice1[i,1]))
    list1 = set(list1)
    
    list2 = []
    for i in range(len(matrice2[:,0])):
        list2.append((matrice2[i,0],matrice2[i,1]))
    list2 = set(list2)
    
    output = list1.symmetric_difference(list2)    
    
    output = list(output)
    output = np.array(output)
    
    return output
