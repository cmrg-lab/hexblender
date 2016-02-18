import bpy
import time


def find_hex(selectionOnly = False):
    # created by Matt Gonzales 2011
    
    #obj = sce.objects.active
    #context.active_object.data

    try:
        obj = bpy.context.active_object
        if obj.type != 'MESH':
            raise TypeError("Active object is not a Mesh")
            return
    except Exception as msg:
        print("Error getting mesh: %s" % msg)
        return

    # Get editmode changes -- are we sure about this?
    # the original command was mesh.update_from_editmode()
    #obj.update_from_editmode()
    #mesh = obj.data


    if len(mesh.polygons) < 1:
        raise ValueError("Mesh has no faces")

    if selectionOnly:
        verts = [vert for vert in mesh.vertices if vert.select]
        if len(verts) == 0:
            print("ERROR: nothing selected!")
            return
        edges = [edge for edge in mesh.edges if edge.select]
        faces = mesh.polygons[mesh.polygons.active]

        nboolSelected = [vert in verts for vert in range(len(mesh.vertices))]
        vertSel_map = list(np.cumsum(np.array(nboolSelected))-1)
        edgeboolSelected = [edge in edges for edge in range(len(mesh.edges))]
        edgeInverse = list(np.cumsum(np.array(edgeboolSelected))-1)
        faceBoolSelected = [face in faces for face in range(len(mesh.polygons))]
        faceInverse = list(np.cumsum(np.array(faceBoolSelected))-1)
    else:
        verts = mesh.vertices
        edges = mesh.edges
        faces = mesh.polygons

    print("begin finding hexes...")
    now = time.time()
    # now = time.time()-now
    #~ print "time: %.5f" %(now)
    now = time.time()
    #get a list of (4) edges that correspond to each face
    faces_edges = [[] for a in range(len(faces))]
    #~ print "faces' edges"
    for ind,face in enumerate(faces):
        for edge in range(4):
            faces_edges[ind].append(mesh.findEdges(*face.edge_keys[edge]))
            #~ print faces_edges[face]

    now = time.time()-now
    #~ print "time: %.5f" %(now)
    now = time.time()
    edges_vertices = [[] for i in range(len(edges))]
    #~ print "edges vertices = "
    for ind,edge in enumerate(edges):
        curr_verts_for_edges = [edge.v1.index, edge.v2.index]
        edges_vertices[ind].extend(curr_verts_for_edges)
        #~ print edges_vertices[edge]

    now = time.time()-now
    print("finished finding lists of vertices associated with each edge in %.1f seconds" %(now))
    now = time.time()
    # print "build a structure that returns the edge numbers associated with one vertex"
    verts_edges = [ [] for i in range(len(verts))]
    #~ print "verts edges = "
    # for i in range(len(verts)):
        # for j in range(len(verts)):
            # if mesh.findEdges(i,j) != None:
                # verts_edges[i].append(mesh.findEdges(i,j))
        #~ print verts_edges[i]

    for ind,vert in enumerate(verts):
        for j in range(len(edges_vertices)):
            if vert.index in edges_vertices[j]:
                if selectionOnly:
                    verts_edges[ind].append(edges[j].index)
                else:
                    verts_edges[ind].append(j)
    #~ print "verts_edges = ", verts_edges

    now = time.time()-now
    #~ print "time: %.5f" %(now)
    now = time.time()
    verts_OneNeighbor = [ [] for i in range(len(verts))]
    #~ print "verts_OneNeighbor = "
    for ind,vert in enumerate(verts):
        #grab edge from a list of edges for vertex "ind"
        for edge in verts_edges[ind]:
            if selectionOnly:
                # edge is the numbered edge in the entire mesh, needs to 
                # be mapped to Selected edge index
                verts_OneNeighbor[ind].extend(edges_vertices[mesh.edges.selected().index(edge)])
            else:
                verts_OneNeighbor[ind].extend(edges_vertices[edge])
        verts_OneNeighbor[ind] = set(verts_OneNeighbor[ind])
        verts_OneNeighbor[ind] = [y for y in verts_OneNeighbor[ind] if y != vert.index]
        #~ print "vert.index = ", vert.index
        #~ print verts_OneNeighbor[ind]

    now = time.time()-now
    print("finished building structure finding vertices' vertex neighbors in %.1f seconds" %(now))
    # now = time.time()
    # print "cycle through faces_edges and return 'faces_edge_vert' which returns the two vertices belonging to any face"
    #index 1 - face
    #index 2 - edge
    #index 3 - length 2 for vertex1 or vertex2
    faces_edge_vert = [ [ [] for j in range(4) ] for i in range(len(faces)) ]
    faces_edge_vert_unique = [ [] for i in range(len(faces)) ]

    for ind,face in enumerate(faces):
        for edge in range(4):
            edge_number = faces_edges[ind][edge]
            if selectionOnly:
                edge_number = mesh.edges.selected().index(edge_number)
            faces_edge_vert[ind][edge].extend([edges[edge_number].v1.index, edges[edge_number].v2.index])

        #make a copy with the unique vertex members
        [faces_edge_vert_unique[ind].extend(x) for x in faces_edge_vert[ind]]
        faces_edge_vert_unique[ind] = list(set(faces_edge_vert_unique[ind]))

    #~ print "faces edge vert = "
    #~ print faces_edge_vert
    #~ print "faces edge vert unique = "
    #~ print faces_edge_vert_unique

    #~ print "faces neighboring faces"

    now = time.time()-now
    #~ print "time: %.5f" %(now)
    now = time.time()
    faces_neighboring_faces = [[] for a in range(len(faces))]
    # print "faces_neighboring_faces\n", faces_neighboring_faces
    # for face in range(len(faces)):
        # curr_face_edge_list = faces_edges[face]
        # print "curr_face_edge_list\n", curr_face_edge_list

        ## inner loop - compare current face to all other faces except itself
        # for face2 in range(len(faces)):
            # if not face2 == face:
                # curr_face_edge_list2 = set(faces_edges[face2])
                ## query to see if any edges are shared between face and face2
                # if [val for val in curr_face_edge_list if val in curr_face_edge_list2] != []:
                    # faces_neighboring_faces[face].append(face2)
    
    faces_edges = [set(myList) for myList in faces_edges]
    for ind1,face1 in enumerate(faces):
        for ind2,face2 in enumerate(faces):
            # only need to compare ind1 for ind1 < ind2 if you are assigning 2
            # values after the 'if'te(faces):
            if ind1==ind2: 
                break 
            #i.e. if the intersection of these sets is not empty
            if faces_edges[ind1] & faces_edges[ind2]: 
                faces_neighboring_faces[ind1].append(face2.index)
                faces_neighboring_faces[ind2].append(face1.index)

    #~ print faces_neighboring_faces

    now = time.time()-now
    print("finished building face neighbor list in %.1f seconds" %(now))
    now = time.time()
    # now populate list of face's neighboring faces, plus the neighboring 
    # faces of those neighboring faces
    #~ print "exclude the 'reference' face"
    #return a unique list
    faces_neighbor2 = [[] for a in range(len(faces))]
    one_neighbor_list = [[] for a in range(len(faces))]
        
    for ind,face in enumerate(faces):

        curr_face_neighbor_list = faces_neighboring_faces[ind]
        for instance in curr_face_neighbor_list:
            if selectionOnly:
                faces_neighbor2[ind].extend(faces_neighboring_faces[mesh.faces.selected().index(instance)])
            else:
                faces_neighbor2[ind].extend(faces_neighboring_faces[instance])

        # first cull out the current face, then make a copy for next routine, 
        # and then find the unique set of face neighbors
        faces_neighbor2[ind] = [val for val in faces_neighbor2[ind] if val != face.index]
        one_neighbor_list[ind] = faces_neighbor2[ind][:] #maybe needs a copy?
        faces_neighbor2[ind] = set(faces_neighbor2[ind])
        faces_neighbor2[ind] = [val for val in faces_neighbor2[ind]]

        #~ print "face #", face
        #~ print faces_neighbor2[ind]

    now = time.time()-now
    #~ print "time: %.5f" %(now)
    now = time.time()
    #~ print "now take one_neighbor_list, which is an uncompressed list of all the neighboring elements' 1-neighbors"
    # (i.e. the original elements' 2-neighbors)
    # if there is any face # which has incidence of 4, there is a hex 
    # associated with the ref. element, the found face #, and the face 
    # identities of the 4 elements which articulated with it
    sixth_faces = [[] for a in range(len(faces))]
    for ind,face in enumerate(faces):
        temp_count_list = []
        for member in one_neighbor_list[ind]:
            temp_count_list.append(one_neighbor_list[ind].count(member))

        #~ print "temp list count, face # ", face
        #~ print temp_count_list

        for member2 in range(len(temp_count_list)):
            if temp_count_list[member2] == 4:
                sixth_faces[ind].append(one_neighbor_list[ind][member2])

        sixth_faces[ind] =set(sixth_faces[ind])
        sixth_faces[ind] = [val for val in sixth_faces[ind]]

    #~ print "sixth faces = "
    # print "sixth_faces", sixth_faces

    now = time.time()-now
    #~ print "time: %.5f" %(now)
    now = time.time()
    # list must be further culled because this alone doesn't filter out 
    # all the hits with 4 neighbors who should be a hex
    #~ print "assemble groups of 6 faces based on the previous routine"
    # first index: global face
    # second index: instance of a hex candidate
    # third index: those 6 faces that MIGHT be a hex
    hex_candidates = [ [ [] for instance in range(len(sixth_faces[a])) ] for a in range(len(faces)) ]
    #~ print "hex candidates = "
    for ind, face in enumerate(faces):
        for instance in range(len(hex_candidates[ind])): 
            # first entry
            hex_candidates[ind][instance].append(face.index)
            # sixth face detected is second entry
            hex_candidates[ind][instance].append(sixth_faces[ind][instance]) 

            # now loop through faces_neighboring_faces for instances of 
            # the sixth face in the first face's neighbor list
            for neigh in faces_neighboring_faces[ind]:
                if selectionOnly:
                    if sixth_faces[ind][instance] in faces_neighboring_faces[mesh.faces.selected().index(neigh)]:
                        hex_candidates[ind][instance].append(neigh)
                else:
                    if sixth_faces[ind][instance] in faces_neighboring_faces[neigh]:
                        hex_candidates[ind][instance].append(neigh)

    #~ for a in range(len(faces)):
        #~ print hex_candidates[a]

    # hex_candidates holds FACE NUMBERS, len8_hex_candidates holds VERTEX NUMBERS

    # this guy will hold the hex candidates that pass the "8 unique vertex" test
    len8_hex_candidates = []

    #~ print "curr vert list = "

    # variable flag_hexes store a list of tuples that contains the
    # index1 - face number
    # index2 - candidate hex instance number
    flag_hexes = []

    now = time.time()-now
    #~ print "time: %.5f" %(now)
    now = time.time()
    #~ print "now go through and check to see if each of these candidates has exactly 8 unique vertices between all the faces"
    for ind,face in enumerate(faces):

        # instance is a FACE instance
        for candidate_hex_instance in range(len(hex_candidates[ind])):

            curr_vert_list = []
            # should always be 6
            for face_inst in range(
                len(hex_candidates[ind][candidate_hex_instance])): 
                # get corresponding face structure for this instance
                if selectionOnly:
                    correspond_face = faces[mesh.faces.selected().index(hex_candidates[ind][candidate_hex_instance][face_inst])]
                else:
                    correspond_face = faces[hex_candidates[ind][candidate_hex_instance][face_inst]]
                for vert in correspond_face.verts:
                    curr_vert_list.append(vert.index)

        #~ print "face = ", face
        #~ print curr_vert_list
        #~ print set(curr_vert_list)

            if len(set(curr_vert_list)) == 8:
                len8_hex_candidates.append(list(set(curr_vert_list)))
                if selectionOnly:
                    #append a tuple
                    flag_hexes.append((mesh.faces.selected().index(face.index),
                        candidate_hex_instance)) 
                else:
                    #append a tuple
                    flag_hexes.append((face.index,candidate_hex_instance)) 

    #~ print "len8 hex candidates = "
    #~ print len8_hex_candidates

    #~ print "flag hexes = "
    #~ print flag_hexes


    unique_hexes = []
    unique_hexes_sorted = []
    unique_face_list = []
    #~ print "unique hexes = "
    for instance in range(len(len8_hex_candidates)):
        vertlist_instance = len8_hex_candidates[instance]
        if sorted(vertlist_instance) not in unique_hexes_sorted:
            unique_hexes.append(vertlist_instance)
            unique_hexes_sorted.append(sorted(vertlist_instance))
            unique_face_list.append(hex_candidates[flag_hexes[instance][0]][flag_hexes[instance][1]])

    # for instance in range(len(unique_hexes)):
        # print unique_hexes[instance]
    # print "unique face list = "
    # for instance in range(len(unique_hexes)):
        # print unique_face_list[instance]

    # get material-number of cubes:
    # primary criterion: highest number of faces for material
    # if equal: lower material-number
    HexMat = []
    for instance in range(len(unique_hexes)):
        HexFacesMat = []
        SumMat = []
        for face in unique_face_list[instance]:
            if selectionOnly:
                HexFacesMat.append(faces[mesh.faces.selected().index(face)].mat)
            else:
                HexFacesMat.append(faces[face].mat)
        setHexFacesMat = set(HexFacesMat)
        setHexFacesMat = list(setHexFacesMat)
        for mat in setHexFacesMat:
            SumMat.append(sum([val == mat for val in HexFacesMat]))
        HexMat.append(setHexFacesMat[SumMat.index(max(SumMat))])
    # print "HexMat", HexMat

    # NOW we have our hexes!  now need to get the proper ordering so 
    # Continuity understands!
    # get the first face and use its vertices as a reference
    # need to call the native face.verts because then you know the vertex 
    # ordering of the face is CW or CCW

    now = time.time()-now
    #~ print "time: %.5f" %(now)
    now = time.time()
    # WARNING we toss out the ordering in variable len8_hex_candidates and 
    # number the first 4 according to
    #~ print "first-numbered surface who we pluck from"
    #~ print "ordered hex final = "
    ordered_hex_vertices = [[] for b in range(len(unique_hexes))]
    for loop_num in range(len(unique_hexes)):
        curr_verts = unique_hexes[loop_num]
        curr_faces = unique_face_list[loop_num]
        curr_first_face = curr_faces[0]

        # get four verts belonging to the first-numbered surface, then 
        # rearrange them so that it is not CW/CCW
        first_four_verts = []
        if selectionOnly:
            for vert in faces[mesh.faces.selected().index(curr_first_face)].verts:
                first_four_verts.append(vert.index)            
        else:
            for vert in faces[curr_first_face].verts:
                first_four_verts.append(vert.index)
        first_four_verts = [first_four_verts[0], first_four_verts[1], first_four_verts[3], first_four_verts[2]]

        remaining_four_verts = []
        for val in curr_verts:
            if val not in first_four_verts:
                remaining_four_verts.append(val)

        #~ print "first four verts = ", first_four_verts
        #~ print "remaining four verts = ", remaining_four_verts

        #find position #5 in eventual hex
        vertex_position1 = first_four_verts[0]
        vertex_position2 = first_four_verts[1]
        vertex_position3 = first_four_verts[2]
        vertex_position4 = first_four_verts[3]
        
        if selectionOnly:
            pos5_vertex = [x for x in remaining_four_verts if x in verts_OneNeighbor[mesh.verts.selected().index(vertex_position1)]]
            pos6_vertex = [y for y in remaining_four_verts if y in verts_OneNeighbor[mesh.verts.selected().index(vertex_position2)]]
            pos7_vertex = [z for z in remaining_four_verts if z in verts_OneNeighbor[mesh.verts.selected().index(vertex_position3)]]
            pos8_vertex = [a for a in remaining_four_verts if a in verts_OneNeighbor[mesh.verts.selected().index(vertex_position4)]]
        else:
            pos5_vertex = [x for x in remaining_four_verts if x in verts_OneNeighbor[vertex_position1]]
            pos6_vertex = [y for y in remaining_four_verts if y in verts_OneNeighbor[vertex_position2]]
            pos7_vertex = [z for z in remaining_four_verts if z in verts_OneNeighbor[vertex_position3]]
            pos8_vertex = [a for a in remaining_four_verts if a in verts_OneNeighbor[vertex_position4]]
            
        #~ print "pos5 vert = ", pos5_vertex
        #~ print "pos6 vert = ", pos6_vertex
        #~ print "pos7 vert = ", pos7_vertex
        #~ print "pos8 vert = ", pos8_vertex

        final_hex = first_four_verts
        final_hex.extend(pos5_vertex)
        final_hex.extend(pos6_vertex)
        final_hex.extend(pos7_vertex)
        final_hex.extend(pos8_vertex)
        ordered_hex_vertices[loop_num].extend(final_hex)

    now = time.time()-now
    #~ print "time: %.5f" %(now)
    print("Found", len(ordered_hex_vertices), "hexes! \n\nFinding contiguous material regions...\n")
    
    # find out whether there are any non-intesecting regions that belong 
    # to the same labeled topology region and separate them
    HexMatNew = contiguous_regions(list(ordered_hex_vertices), HexMat)
    
    print("Harmonizing topology regions...")
    
    # find verts and put them into an array
    verts_array = np.zeros([len(mesh.verts),3])
    for ind, vert in enumerate(mesh.verts):
        verts_array[ind,:] = vert.co
    
    # harmonize hexes
    ordered_hex_vertices = harmonize_topo(ordered_hex_vertices,
                                          verts_array,
                                          HexMatNew)
    
    # save the elems and mats so you don't have to find every time that 
    # you do something
    if not selectionOnly:
        cache_data(ordered_hex_vertices, HexMatNew)
    
    print("Completed hex_find")
    
    return ordered_hex_vertices, HexMat
