import numpy as np

def harmonize_topo(elems,v,MatList):

    finished_topo_flag = False
    harmonized_elems = []
    topoListNums = np.unique(MatList)
    num_topologies = len(topoListNums)
    elems = np.array(elems)
    
    #VARIABLE: elems_one_topo_group
    #INDEX 1:  topology group (max = num_topologies)
    #INDEX 2: instance of elems in topo group
    elems_one_topo_group = [[] for i in range(num_topologies)]
        
    #initialize as a list of lists
    for region in range(num_topologies):
        elems_one_topo_group[region] = np.array(MatList == topoListNums[region]).nonzero()[0].tolist()
        elems_one_topo_group[region] = [val for val in elems_one_topo_group[region]]

    #now for target harmonious topologies, construct a for loop with an inner while loop
    for region_number in range(num_topologies):
        
        harmonized_elems = []
        print("Harmonizing Region %i (%i elements)..." %(region_number,len(elems_one_topo_group[region_number])))
        
        while finished_topo_flag == False:
            if np.shape(elems)[1] == 8:
                elems, harmonized_elems, finished_topo_flag = HarmonizeOneRegion3D(
                                                elems,
                                                v,
                                                elems_one_topo_group[region_number],
                                                harmonized_elems,
                                                region_number)

            else: #np.shape(elems)[1] should be '4'
                elems, harmonized_elems, finished_topo_flag = HarmonizeOneRegion2D(
                                                elems,
                                                v,
                                                elems_one_topo_group[region_number],
                                                harmonized_elems,
                                                region_number)

            #~ print "harmonized_elems = ", sorted(harmonized_elems)
            #~ print "length harmonized = ", len(harmonized_elems)
        
        finished_topo_flag = False
    
    elems = elems.tolist()
    
    return elems
    
def HarmonizeOneRegion3D(elems, v, curr_elem_list, harmonized_elems, topo_number):

    finished_topo_flag = False
    
    #~ print " "
    #~ print "EXECUTING ONE ITERATION OF harmonizeOneRegion()"
    #~ print "curr_elem_list in topology = ", curr_elem_list
    #~ print "harmonized elems on entering = ", harmonized_elems
    
    #if none of curr elems belong to list  harmonized elements, look for one whose topology is regular
    if [val for val in curr_elem_list if val in harmonized_elems] == []:
        
        #use a "seed" element whose topology will be enforced on all surrounding elements
        harmonized_elems.append(curr_elem_list[0])
        
        #make sure the seed element is right-handed
        vertex0 = elems[curr_elem_list[0]][0] #origin
        vertex1 = elems[curr_elem_list[0]][1] #u-direction
        vertex2 = elems[curr_elem_list[0]][2] #v-direction
        vertex4 = elems[curr_elem_list[0]][4] #w-direction
        
        u_vec = np.array([ v[vertex1][0] - v[vertex0][0]  , v[vertex1][1] - v[vertex0][1]  , v[vertex1][2] - v[vertex0][2] ])
        v_vec = np.array([ v[vertex2][0] - v[vertex0][0]  , v[vertex2][1] - v[vertex0][1]  , v[vertex2][2] - v[vertex0][2] ])
        w_vec = np.array([ v[vertex4][0] - v[vertex0][0]  , v[vertex4][1] - v[vertex0][1]  , v[vertex4][2] - v[vertex0][2] ])
        
        #note this does not execute when triple product = 0, as is the case when I am making lazy test cases
        if np.dot( u_vec , np.cross(v_vec,w_vec) ) < 0:
            ce = curr_elem_list[0] #ce = "current element"
            elems[ce,:] = [ elems[ce][1], elems[ce][0], elems[ce][3], elems[ce][2], elems[ce][5], elems[ce][4], elems[ce][7], elems[ce][6] ]
        
    #~ print "harmonized elems after calling generator = ", harmonized_elems
    
    #now, create a set of seed elems
    harm_elems_in_topo = [ne for ne in curr_elem_list if ne in harmonized_elems]
    unharmonized_elems_in_region = [ne for ne in curr_elem_list if ne not in harmonized_elems]
    
    #~ print "unharmonized elems in region = ", unharmonized_elems_in_region

    elems_to_harmonize = []
    ref_elems_to_harm = []
    
    #query the instances of other elements in the topology group sharing exactly 4 vertices (i.e. a face if no duplicate nodes)
    setElems = elems.tolist()
    setElems = [set(hex) for hex in setElems]
    for unharm_ne in unharmonized_elems_in_region:
        #get the current nodes of the unharmonized element "ne"
        
        #compare nodes of the unharmonized element "ne" with each element already harmonized
        for harm_ne in harm_elems_in_topo:
            
            if setElems[unharm_ne] & setElems[harm_ne]:
                if len(setElems[unharm_ne] & setElems[harm_ne]) == 4:
                    elems_to_harmonize.append(unharm_ne)
                    ref_elems_to_harm.append(harm_ne)
                    #break out once you find a single instance satisfying this criterion
                    break
    
    #~ print "elems to harmonize (sharing face) = ", elems_to_harmonize
    #~ print "ref elems sharing faces with current elems = ", ref_elems_to_harm
    
    for ind, harm_ne in enumerate(elems_to_harmonize):
        curr_ref_elem = ref_elems_to_harm[ind]
        
        elems2harm_nodes = elems[harm_ne,:]
        ref_elems_nodes = elems[curr_ref_elem,:]
        
        elems[harm_ne,:] = Harmonize2Elems3D(ref_elems_nodes, elems2harm_nodes)
        harmonized_elems.append(harm_ne)
        
    if len(harmonized_elems) == len(curr_elem_list):
        finished_topo_flag = True

    return elems, harmonized_elems, finished_topo_flag
        
def Harmonize2Elems3D(ref_elem, other_elem):
    
    ref_elem = list(ref_elem)
    other_elem = list(other_elem)

    ref_elem_shared_indices = [ind for ind,val in enumerate(ref_elem) if val in other_elem]
    ref_elem_shared_verts = [val for ind,val in enumerate(ref_elem) if val in other_elem]
    other_elem_shared_indices = [other_elem.index(val) for val in ref_elem_shared_verts]
        
    #find the plane of interest for "other_elem" and make a list of tuples for the conjugate values

    #~ print "ref elem shared indices = ", ref_elem_shared_indices
    #~ print "other elem shared indices = ", other_elem_shared_indices
    
    #case 1: "other elem" plane xi1 = 0 is shared
    if sorted(other_elem_shared_indices) == [0,2,4,6]:
        other_elem_tuples = [(val,val+1) for val in other_elem_shared_indices]
    #case2: plane xi1 = 1 is shared
    elif sorted(other_elem_shared_indices) == [1,3,5,7]:                
        other_elem_tuples = [(val,val-1) for val in other_elem_shared_indices]
    #case3: plane xi2 = 0 is shared
    elif sorted(other_elem_shared_indices) == [0,1,4,5]:
        other_elem_tuples = [(val,val+2) for val in other_elem_shared_indices]
    #case4: plane xi2 = 1 is shared
    elif sorted(other_elem_shared_indices) == [2,3,6,7]:
        other_elem_tuples = [(val,val-2) for val in other_elem_shared_indices]
    #case5: plane xi3 = 0 is shared
    elif sorted(other_elem_shared_indices) == [0,1,2,3]:
        other_elem_tuples = [(val,val+4) for val in other_elem_shared_indices]
    #case6: plane xi3 = 1 is shared
    elif sorted(other_elem_shared_indices) == [4,5,6,7]:
        other_elem_tuples = [(val,val-4) for val in other_elem_shared_indices]
    
    #~ print "other elem tuples = ", other_elem_tuples
    
    #default is no change
    correct_index_order = [0,1,2,3,4,5,6,7]
    
    #case 1: plane xi1 = 0 is shared
    if sorted(ref_elem_shared_indices) == [0,2,4,6]:
        
        correct_index_order[0] = other_elem_tuples[0][1]
        correct_index_order[1] = other_elem_shared_indices[0]
        correct_index_order[2] = other_elem_tuples[1][1]
        correct_index_order[3] = other_elem_shared_indices[1]
        correct_index_order[4] = other_elem_tuples[2][1]
        correct_index_order[5] = other_elem_shared_indices[2]
        correct_index_order[6] = other_elem_tuples[3][1]
        correct_index_order[7] = other_elem_shared_indices[3]
        
    #case2: plane xi1 = 1 is shared
    elif sorted(ref_elem_shared_indices) == [1,3,5,7]:
        
        #can pull the correct harmonized nodes from the ordered ref_elem_shared_indices   
        correct_index_order[0] = other_elem_shared_indices[0]
        correct_index_order[1] = other_elem_tuples[0][1]
        correct_index_order[2] = other_elem_shared_indices[1]
        correct_index_order[3] = other_elem_tuples[1][1]
        correct_index_order[4] = other_elem_shared_indices[2]
        correct_index_order[5] = other_elem_tuples[2][1]
        correct_index_order[6] = other_elem_shared_indices[3]
        correct_index_order[7] = other_elem_tuples[3][1]
        
    #case3: plane xi2 = 0 is shared
    elif sorted(ref_elem_shared_indices) == [0,1,4,5]:
        
        correct_index_order[0] = other_elem_tuples[0][1]
        correct_index_order[1] = other_elem_tuples[1][1]
        correct_index_order[2] = other_elem_shared_indices[0]
        correct_index_order[3] = other_elem_shared_indices[1]
        correct_index_order[4] = other_elem_tuples[2][1]
        correct_index_order[5] = other_elem_tuples[3][1]
        correct_index_order[6] = other_elem_shared_indices[2]
        correct_index_order[7] = other_elem_shared_indices[3]
        
    #case4: plane xi2 = 1 is shared
    elif sorted(ref_elem_shared_indices) == [2,3,6,7]:
        
        correct_index_order[0] = other_elem_shared_indices[0]
        correct_index_order[1] = other_elem_shared_indices[1]
        correct_index_order[2] = other_elem_tuples[0][1]
        correct_index_order[3] = other_elem_tuples[1][1]
        correct_index_order[4] = other_elem_shared_indices[2]
        correct_index_order[5] = other_elem_shared_indices[3]
        correct_index_order[6] = other_elem_tuples[2][1]
        correct_index_order[7] = other_elem_tuples[3][1]
        
    #case5: plane xi3 = 0 is shared
    elif sorted(ref_elem_shared_indices) == [0,1,2,3]:
        
        correct_index_order[0] = other_elem_tuples[0][1]
        correct_index_order[1] = other_elem_tuples[1][1]
        correct_index_order[2] = other_elem_tuples[2][1]
        correct_index_order[3] = other_elem_tuples[3][1]
        correct_index_order[4] = other_elem_shared_indices[0]
        correct_index_order[5] = other_elem_shared_indices[1]
        correct_index_order[6] = other_elem_shared_indices[2]
        correct_index_order[7] = other_elem_shared_indices[3]
        
    #case6: plane xi3 = 1 is shared
    elif sorted(ref_elem_shared_indices) == [4,5,6,7]:
        
        correct_index_order[0] = other_elem_shared_indices[0]
        correct_index_order[1] = other_elem_shared_indices[1]
        correct_index_order[2] = other_elem_shared_indices[2]
        correct_index_order[3] = other_elem_shared_indices[3]
        correct_index_order[4] = other_elem_tuples[0][1]
        correct_index_order[5] = other_elem_tuples[1][1]
        correct_index_order[6] = other_elem_tuples[2][1]
        correct_index_order[7] = other_elem_tuples[3][1]
        
    #~ print "correct index order = ", correct_index_order
    
    other_elem = np.array(other_elem)
    harmonized_elem = other_elem[correct_index_order]
        
    return harmonized_elem
    
def HarmonizeOneRegion2D(elems, v, curr_elem_list, harmonized_elems, topo_number):

    finished_topo_flag = False
    
    #~ print " "
    #~ print "EXECUTING ONE ITERATION OF harmonizeOneRegion()"
    #~ print "curr_elem_list in topology = ", curr_elem_list
    #~ print "harmonized elems on entering = ", harmonized_elems
    
    #if none of curr elems belong to list  harmonized elements, look for one whose topology is regular
    if [val for val in curr_elem_list if val in harmonized_elems] == []:
        
        #use a "seed" element whose topology will be enforced on all surrounding elements
        harmonized_elems.append(curr_elem_list[0])
        
    #~ print "harmonized elems after calling generator = ", harmonized_elems
    
    #now, create a set of seed elems
    harm_elems_in_topo = [ne for ne in curr_elem_list if ne in harmonized_elems]
    unharmonized_elems_in_region = [ne for ne in curr_elem_list if ne not in harmonized_elems]
    
    #~ print "unharmonized elems in region = ", unharmonized_elems_in_region

    elems_to_harmonize = []
    ref_elems_to_harm = []
    
    #query the instances of other elements in the topology group sharing exactly 4 vertices (i.e. a face if no duplicate nodes)
    for unharm_ne in unharmonized_elems_in_region:
        #get the current nodes of the unharmonized element "ne"
        unharmonized_nodes = set(elems[unharm_ne,:])
        
        #compare nodes of the unharmonized element "ne" with the 
        for harm_ne in harm_elems_in_topo:
            harmonized_nodes = set(elems[harm_ne,:])
            
            if len(unharmonized_nodes.intersection(harmonized_nodes)) == 2:
                elems_to_harmonize.append(unharm_ne)
                ref_elems_to_harm.append(harm_ne)
                #break out once you find a single instance satisfying this criterion
                break
    
    #~ print "elems to harmonize (sharing face) = ", elems_to_harmonize
    #~ print "ref elems sharing faces with current elems = ", ref_elems_to_harm
    
    for ind, harm_ne in enumerate(elems_to_harmonize):
        curr_ref_elem = ref_elems_to_harm[ind]
        
        elems2harm_nodes = elems[harm_ne,:]
        ref_elems_nodes = elems[curr_ref_elem,:]
        
        elems[harm_ne,:] = Harmonize2Elems2D(ref_elems_nodes, elems2harm_nodes)
        harmonized_elems.append(harm_ne)
        
    if len(harmonized_elems) == len(curr_elem_list):
        finished_topo_flag = True

    return elems, harmonized_elems, finished_topo_flag
    
def Harmonize2Elems2D(ref_elem, other_elem):
    
    ref_elem = list(ref_elem)
    other_elem = list(other_elem)
    
    ref_elem_shared_indices = [ind for ind,val in enumerate(ref_elem) if val in other_elem]
    ref_elem_shared_verts = [val for ind,val in enumerate(ref_elem) if val in other_elem]
    other_elem_shared_indices = [other_elem.index(val) for val in ref_elem_shared_verts]
        
    #find the plane of interest for "other_elem" and make a list of tuples for the conjugate values

    #~ print "ref elem shared indices = ", ref_elem_shared_indices
    #~ print "other elem shared indices = ", other_elem_shared_indices
    
    #case 1: "other elem" line xi1 = 0 is shared
    if sorted(other_elem_shared_indices) == [0,2]:
        other_elem_tuples = [(val,val+1) for val in other_elem_shared_indices]
    #case2: line xi1 = 1 is shared
    elif sorted(other_elem_shared_indices) == [1,3]:                
        other_elem_tuples = [(val,val-1) for val in other_elem_shared_indices]
    #case3: line xi2 = 0 is shared
    elif sorted(other_elem_shared_indices) == [0,1]:
        other_elem_tuples = [(val,val+2) for val in other_elem_shared_indices]
    #case4: line xi2 = 1 is shared
    elif sorted(other_elem_shared_indices) == [2,3]:
        other_elem_tuples = [(val,val-2) for val in other_elem_shared_indices]
    
    #~ print "other elem tuples = ", other_elem_tuples
    
    #default is no change
    correct_index_order = [0,1,2,3]
    
    #case 1: line xi1 = 0 is shared
    if sorted(ref_elem_shared_indices) == [0,2]:
        
        correct_index_order[0] = other_elem_tuples[0][1]
        correct_index_order[1] = other_elem_shared_indices[0]
        correct_index_order[2] = other_elem_tuples[1][1]
        correct_index_order[3] = other_elem_shared_indices[1]
        
    #case2: line xi1 = 1 is shared
    elif sorted(ref_elem_shared_indices) == [1,3]:
        
        #can pull the correct harmonized nodes from the ordered ref_elem_shared_indices   
        correct_index_order[0] = other_elem_shared_indices[0]
        correct_index_order[1] = other_elem_tuples[0][1]
        correct_index_order[2] = other_elem_shared_indices[1]
        correct_index_order[3] = other_elem_tuples[1][1]
        
    #case3: line xi2 = 0 is shared
    elif sorted(ref_elem_shared_indices) == [0,1]:
        
        correct_index_order[0] = other_elem_tuples[0][1]
        correct_index_order[1] = other_elem_tuples[1][1]
        correct_index_order[2] = other_elem_shared_indices[0]
        correct_index_order[3] = other_elem_shared_indices[1]
        
    #case4: line xi2 = 1 is shared
    elif sorted(ref_elem_shared_indices) == [2,3]:
        
        correct_index_order[0] = other_elem_shared_indices[0]
        correct_index_order[1] = other_elem_shared_indices[1]
        correct_index_order[2] = other_elem_tuples[0][1]
        correct_index_order[3] = other_elem_tuples[1][1]
        
    #~ print "correct index order = ", correct_index_order
    
    other_elem = np.array(other_elem)
    harmonized_elem = other_elem[correct_index_order]
        
    return harmonized_elem
