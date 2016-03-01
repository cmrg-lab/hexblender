import bpy
import time
import json
import bmesh
import numpy as np
import pickle
from hexblender.harmonize_hex import harmonize_topo

def get_active_mesh():
    ''' Simply returns the current active mesh '''
    orig_mode = 'OBJECT'
    if hasattr(bpy.context.object, "mode") and bpy.context.object.mode != 'OBJECT' :
        orig_mode = bpy.context.object.mode
        print("NOTICE: Changing to OBJECT mode...")
        bpy.ops.object.mode_set(mode='OBJECT')

    if bpy.context.active_object is not None and\
       bpy.context.active_object.type == "MESH":
        print("Getting mesh with name: %s" % bpy.context.active_object.data.name)
        return bpy.context.active_object.data, orig_mode
    else:
        print("No active mesh found!")
        return None, orig_mode


def get_active_mesh_and_bmesh():
    ''' Simply returns the current active mesh in a bmesh'''

    mesh, orig_mode = get_active_mesh()

    if mesh is not None:
        if bpy.context.object.mode != 'EDIT':
            print("NOTICE: We need to be in EDIT MODE for this command!")
            print("NOTICE: Changing to EDIT mode...")
            bpy.ops.object.mode_set(mode='EDIT')

        return mesh, bmesh.from_edit_mesh(mesh)
    else:
        print("Error: Unable to create bmesh!")


def reset_mode(orig_mode):
    if orig_mode != 'OBJECT':
        print("Changing back into mode: %s" % orig_mode)
        bpy.ops.object.mode_set(mode=orig_mode)


def get_material_data(mesh = None):

    try:
        faces = mesh.faces
    except Exception as msg:
        mesh, orig_mode = get_active_mesh()
        faces = mesh.polygons

    matList = [[] for val in faces]

    for face in faces:
        matList[face.index] = face.material_index

    return matList


def cache_data(ordered_hex_vertices, hex_mat_new):
    """ Retain information of original mesh.

    When a mesh comes from Continutiy we do not want to lose that information, 
    so we will store it here.  That information is:
    number of elements
    which nodes are connected to which elements
    number of materials
    """
        
    mesh, orig_mode = get_active_mesh()

    # Create a property to store the original mesh information
    # There might be a better way to do this, but this seems to work
    # for now.
    bpy.types.Mesh.cached_data_json = bpy.props.StringProperty()

    # see if we have already cached the data
    try:
        cached_data = bpy.types.Mesh.cached_data[mesh.name]
    except Exception as msg:
        cached_data = {}

    for ind, hex_vert in enumerate(ordered_hex_vertices):
        cached_data[str(ind)] = hex_vert

    cached_data['num_elems'] = len(ordered_hex_vertices)
    cached_data['mats'] = hex_mat_new

    # Store the cached data in our property, which seems to be the only?
    # way to have the data serialized in a blend file
    mesh.cached_data_json = json.dumps(cached_data)
    
    # XXX Leaving for reference should we run into problems with larger meshes
    '''
    # the mesh properties class has problems when the length of an array 
    # is longer than about 10,000 elements, so I break up the material list 
    # into lists of 10,000 and concatenate them in ReloadHex()
    numMatLists = (int(len(hex_mat_new)) / 10000) + 1     # floor( length / 10000 ) + 1
    for matListNumber in range(numMatLists):
        stringKey = 'mats' + str(matListNumber)
        mesh.properties[stringKey] = hex_mat_new[matListNumber*10000 : (matListNumber+1)*10000]
        #python smart enough to stop slicing at the last entry (since modulo generally won't be zero)
    mesh.properties['numMatLists'] = numMatLists
    '''
    print("cached data for %d hexes" % len(ordered_hex_vertices))


def set_mat_data_hex(matList, new_cube_faces, mesh = None, subdivided=True):

    if mesh is None:
        print("Mesh is None, getting new mesh!!!!")
        mesh, orig_mode = get_active_mesh()

    faces = mesh.polygons

    print("#### working on mesh with name: %s" % mesh.name)
    nMatList = []
    for i in range(len(matList)):
        if subdivided:
            for j in range(8):
                nMatList.append(matList[i])
        else:
            nMatList.append(matList[i])

    facesInVerts = []
    for face in faces:
        c_face = []
        c_face.append(int(face.vertices[0]))
        c_face.append(int(face.vertices[1]))
        c_face.append(int(face.vertices[2]))
        c_face.append(int(face.vertices[3]))
        cSet_face = sorted(c_face)
        facesInVerts.append(cSet_face)

    start = time.time()
    for i,face in enumerate(new_cube_faces):
        c_face = sorted(face)
        # if c_face in facesInVerts:
        c_Mat = nMatList[int(i/6)]
        # else:
        if faces[facesInVerts.index(c_face)].material_index == 0 or (faces[facesInVerts.index(c_face)].material_index > c_Mat and c_Mat != 0):
            faces[facesInVerts.index(c_face)].material_index = c_Mat
    print("cube_faces loop took: %s" % (time.time() - start))

    # XXX: define the material instances so the user doesn't have to create
    # the materials manually, although they'll still need to specify their
    # colors, unless we also want to do that for them.
    #bpy.ops.material.new()
    #bpy.ops.object.material_slot_add()


def add_new_data(verts_new, edges_keys_index, faces_keys_index, cubes_verts_index, newObj = False):

    new_cube_faces = []
    # takes the clockwise order of verts
    for cube in cubes_verts_index:
        new_cube_faces.append([int(value) for value in cube[:4]])
        new_cube_faces.append([int(value) for value in cube[4:]])

        for i in range(4):
            new_cube_faces.append([int(cube[i]),
                                   int(cube[(i+1)%4]),
                                   int(cube[4+(i+1)%4]),
                                   int(cube[(i+4)])])
    # add all new faces
    for face_key in faces_keys_index:
        new_cube_faces.append(face_key)
   
    if newObj:
        mesh = bpy.data.meshes.new('pickledObject')
        mesh.from_pydata(verts_new,
                         [], #edges
                         new_cube_faces, #faces
                        )

        mesh.update()
        print("Validate mesh!")
        mesh.validate(verbose=True)

        # create the new mesh and make it active 
        obj = bpy.data.objects.new('pickledObject', mesh)
        bpy.context.scene.objects.link(obj)
        bpy.context.scene.objects.active = obj

    else:
        ''' old way
        sce = bpy.data.scenes.active
        ob = sce.objects.active
        mesh = ob.getData(mesh=1)
        '''
        mesh, orig_mode = get_active_mesh()

        faces = b_mesh.faces
        edges = b_mesh.edges
        verts = b_mesh.verts

        now = time.time()

        # add all new verts
        # not using bmesh

        # find how many verts we need to add then do
        vert.new(num_new_nodes)

        # now enter each coordinat
        count = 0
        for index in range(num_old_verts, (num_old_verts + num_new_verts)):
            verts[index].co = new_verts[counter]
            count += 0

        counter = 0
        orig_vert_num = len(verts)
        verts.ensure_lookup_table()
        counter = 0
        face_counter = 0
        for vert in verts_new:
            verts.new(vert)
            verts.ensure_lookup_table()
            verts[-1].index = orig_vert_num + counter
            counter += 1
            face_counter += 1
            # add the last three verts to create a new face
            if face_counter == 3:
                faces.ensure_lookup_table()
                faces.new((verts[i] for i in range(-3,0))) 
                faces.ensure_lookup_table()
                face_counter = 0
        verts.ensure_lookup_table()
        faces.ensure_lookup_table()

        # add all new edges
        for i in range(len(edges_keys_index)):
            edges.new(edges_keys_index[i])

        # add all new faces
        # do we really need to explicity add faces
        # even after addind verts and edges?
        for i in range(len(faces_keys_index)):
            faces.new(faces_keys_index[i])
        #faces = mesh.faces
        #edges = mesh.edges

        #for f in new_cube_faces:
        #    bm.faces.new(new_cube_faces)

        print("in add new data, # faces: %s" % len(faces))

    return new_cube_faces

def reload_hex(cached_data):
    """ Will retrieve the orginal element data from the Continuity pickle file """

    num_elems = cached_data['num_elems']
    cubes = []
    for index in range(num_elems):
        cubes.append(cached_data[str(index)])

    mat_list = cached_data['mats']

    return cubes, mat_list

def get_boundary_data(cubes, mesh = None):

    try:
        faces = mesh.faces
        verts = mesh.verts
    except Exception as msg:
        mesh, orig_mode = get_active_mesh()
        faces = mesh.polygons
        verts = mesh.vertices

    edges = mesh.edges

    allVerts = []
    for vert in verts:
        v_co = []
        v_co.append(vert.co.x)
        v_co.append(vert.co.y)
        v_co.append(vert.co.z)
        allVerts.append(v_co)

    allFaces = []
    for face in faces:
        face_in_verts = []
        face_in_verts.append(int(face.vertices[0]))
        face_in_verts.append(int(face.vertices[1]))
        face_in_verts.append(int(face.vertices[2]))
        face_in_verts.append(int(face.vertices[3]))
        allFaces.append(face_in_verts)

    cubes_trans = []
    for cube in cubes:
        cube_new = []
        cube_new.append(cube[0])
        cube_new.append(cube[1])
        cube_new.append(cube[3])
        cube_new.append(cube[2])
        cube_new.append(cube[4])
        cube_new.append(cube[5])
        cube_new.append(cube[7])
        cube_new.append(cube[6])
        cubes_trans.append(cube_new)

    return cubes_trans, allFaces, allVerts


def contiguous_regions(cubes, HexMat):
    
    topologies = np.unique(HexMat)
    tot_topologies = len(topologies)
    HexMatNew = [None for val in HexMat]
    final_index_lists = []    
    
    for topo_num in topologies:
        # get indices of hexes that are in the current region
        curr_topo_hex_indices = [ind for ind, val in enumerate(HexMat) if val == topo_num]
        curr_topo_node_lists = [cubes[ind][:] for ind, val in enumerate(HexMat) if val == topo_num ]
        
        # initialize remaining topo indices to curr topo hex indices
        remaining_topo_indices = curr_topo_hex_indices[:]
        one_topos_groups = []
        curr_topo_nodes = []
        
        #~ print "curr_topo_hex_indices = ", curr_topo_hex_indices
        #~ print "curr topo node lists = ", curr_topo_node_lists
        
        curr_topo_it = 0
        continueLoop1 = True
        while continueLoop1:
            
            continueLoop2 = True
            initialized_elem = remaining_topo_indices.pop(0)
            one_topos_groups.append([])
            one_topos_groups[curr_topo_it].append(initialized_elem)
            curr_topo_nodes.append(list(cubes[initialized_elem][:]))
            
            while continueLoop2:

                #~ print "one_topos_groups = ", one_topos_groups
                
                init_length = len(one_topos_groups[curr_topo_it])
                
                #~ print "init length = ", init_length
                #~ print "curr topo it = ", curr_topo_it
                
                #~ print "remaining topo indices = ", remaining_topo_indices
                indices_to_remove = []
                
                for ind in remaining_topo_indices:
                    #~ print "ind = ", ind
                    
                    curr_elem_nodes = cubes[ind][:]
                    #~ print "curr elem nodes = ", curr_elem_nodes
                    #~ print "curr topo nodes = ", curr_topo_nodes[curr_topo_it]
                    #~ print [val in curr_topo_nodes[curr_topo_it] for val in curr_elem_nodes]
                    
                    if any([val in curr_topo_nodes[curr_topo_it] for val in curr_elem_nodes]):
                        if ind not in one_topos_groups[curr_topo_it]:
                            #~ print "ind in add= ", ind
                            one_topos_groups[curr_topo_it].append(ind)
                            curr_topo_nodes[curr_topo_it].extend(curr_elem_nodes)
                            curr_topo_nodes[curr_topo_it] = list(set(curr_topo_nodes[curr_topo_it]))
                            indices_to_remove.append(ind)
                        
                    #~ print "one topo groups = ", one_topos_groups
                    #~ print "curr topo nodes = ", curr_topo_nodes
                
                for val in indices_to_remove:
                    remaining_topo_indices.remove(val)
                
                if init_length == len(one_topos_groups[curr_topo_it]): #no change
                    #~ print "final one topo group = ", one_topos_groups
                    continueLoop2 = False
            
            curr_topo_it += 1
            if remaining_topo_indices == []:
                continueLoop1 = False
                
        final_index_lists.append(one_topos_groups)
    
    # find the total number of new topology regions
    counter = 0
    for instance in final_index_lists:
        counter += len(instance)
    num_regions = counter
    counter2 = 0
    for instance in final_index_lists:
        for region in instance:
            counter2 += 1
            for ele_ind in region:
                HexMatNew[ele_ind] = counter2-1
    # currently this does not preserve old topology region numbers... 
    # maybe add this later if we need it
    return HexMatNew


def contiguous_regions(cubes, HexMat):
    
    topologies = np.unique(HexMat)
    tot_topologies = len(topologies)
    HexMatNew = [None for val in HexMat]
    final_index_lists = []    
    
    for topo_num in topologies:
        #get indices of hexes that are in the current region
        curr_topo_hex_indices = [ind for ind, val in enumerate(HexMat) if val == topo_num]
        curr_topo_node_lists = [cubes[ind][:] for ind, val in enumerate(HexMat) if val == topo_num ]
        
        #initialize remaining topo indices to curr topo hex indices
        remaining_topo_indices = curr_topo_hex_indices[:]
        one_topos_groups = []
        curr_topo_nodes = []
        
        #~ print "curr_topo_hex_indices = ", curr_topo_hex_indices
        #~ print "curr topo node lists = ", curr_topo_node_lists
        
        curr_topo_it = 0
        continueLoop1 = True
        while continueLoop1:
            
            continueLoop2 = True
            initialized_elem = remaining_topo_indices.pop(0)
            one_topos_groups.append([])
            one_topos_groups[curr_topo_it].append(initialized_elem)
            curr_topo_nodes.append(list(cubes[initialized_elem][:]))
            
            while continueLoop2:

                #~ print "one_topos_groups = ", one_topos_groups
                
                init_length = len(one_topos_groups[curr_topo_it])
                
                #~ print "init length = ", init_length
                #~ print "curr topo it = ", curr_topo_it
                
                #~ print "remaining topo indices = ", remaining_topo_indices
                indices_to_remove = []
                
                for ind in remaining_topo_indices:
                    #~ print "ind = ", ind
                    
                    curr_elem_nodes = cubes[ind][:]
                    #~ print "curr elem nodes = ", curr_elem_nodes
                    #~ print "curr topo nodes = ", curr_topo_nodes[curr_topo_it]
                    #~ print [val in curr_topo_nodes[curr_topo_it] for val in curr_elem_nodes]
                    
                    if any([val in curr_topo_nodes[curr_topo_it] for val in curr_elem_nodes]):
                        if ind not in one_topos_groups[curr_topo_it]:
                            #~ print "ind in add= ", ind
                            one_topos_groups[curr_topo_it].append(ind)
                            curr_topo_nodes[curr_topo_it].extend(curr_elem_nodes)
                            curr_topo_nodes[curr_topo_it] = list(set(curr_topo_nodes[curr_topo_it]))
                            indices_to_remove.append(ind)
                        
                    #~ print "one topo groups = ", one_topos_groups
                    #~ print "curr topo nodes = ", curr_topo_nodes
                
                for val in indices_to_remove:
                    remaining_topo_indices.remove(val)
                
                if init_length == len(one_topos_groups[curr_topo_it]): #no change
                    #~ print "final one topo group = ", one_topos_groups
                    continueLoop2 = False
            
            curr_topo_it += 1
            if remaining_topo_indices == []:
                continueLoop1 = False
                
        final_index_lists.append(one_topos_groups)
    
    #find the total number of new topology regions
    counter = 0
    for instance in final_index_lists:
        counter += len(instance)
    num_regions = counter
    counter2 = 0
    for instance in final_index_lists:
        for region in instance:
            counter2 += 1
            for ele_ind in region:
                HexMatNew[ele_ind] = counter2-1
    #currently this does not preserve old topology region numbers... maybe add this later if we need it
    return HexMatNew


def find_hex(selectionOnly = False, verts = None):

    # created by Matt Gonzales 2011
    
    mesh, orig_mode = get_active_mesh()

    ###################################################
    ############### DEBUG FIND HEX
    ###################################################
    do_sort = False

    if selectionOnly:
        if verts is None:
            print("!!!! verts is None!!!!")
            #verts = [v for v in mesh.vertices if v.select]
            verts = [v for v in b_mesh.verts if v.select]
            if len(verts) == 0:
                print("ERROR: nothing is selected!")
                return

        # Get selected edge and polygon objects
        edges = [e for e in mesh.edges if e.select]
        faces = [p for p in mesh.polygons if p.select]
            
        # The old way to get things
        #verts = [mesh.verts[ind] for ind in mesh.verts.selected()]
        #edges = [mesh.edges[ind] for ind in mesh.edges.selected()]
        #faces = [mesh.faces[ind] for ind in mesh.faces.selected()]

        # Seems like we're getting ALL objects, not just selected ones
        nboolSelected = [vert in verts for vert in range(len(mesh.vertices))]
        vertSel_map = list(np.cumsum(np.array(nboolSelected))-1)
        edgeboolSelected = [edge in edges for edge in range(len(mesh.edges))]
        edgeInverse = list(np.cumsum(np.array(edgeboolSelected))-1)
        faceBoolSelected = [face in faces for face in range(len(mesh.polygons))]
        faceInverse = list(np.cumsum(np.array(faceBoolSelected))-1)

        # For ease of compatibilty with the existing code, we need a list
        # containg the selected edge numbers
        selected_verts = [v.index for v in verts]
        selected_edges = [e.index for e in edges]
        selected_faces = [f.index for f in faces]
    else:
        verts = mesh.vertices
        edges = mesh.edges
        faces = mesh.polygons

    print("begin finding hexes...\n")
    now = time.time()
    # now = time.time()-now
    #~ print "time: %.5f" %(now)
    now = time.time()
    #get a list of (4) edges that correspond to each face
    faces_edges = [[] for a in range(len(faces))]
    #~ print "faces' edges"
    for ind,face in enumerate(faces):
        for edge in range(4):
            faces_edges[ind].append(mesh.edge_keys.index(face.edge_keys[edge]))
            # old way
            #faces_edges[ind].append(mesh.findEdges(*face.edge_keys[edge]))
            #~ print faces_edges[face]

    with open('/Users/jeffvandorn/faces_edges_2_72.pickle', mode='wb') as f:
        pickle.dump(faces_edges, f, 2)

    now = time.time()-now
    #~ print "time: %.5f" %(now)
    now = time.time()
    edges_vertices = [[] for i in range(len(edges))]
    #~ print "edges vertices = "
    for ind,edge in enumerate(edges):
        # old way
        #curr_verts_for_edges = [edge.v1.index, edge.v2.index]
        curr_verts_for_edges = [edge.vertices[0], edge.vertices[1]]
        edges_vertices[ind].extend(curr_verts_for_edges)
        #~ print edges_vertices[edge]

    with open('/Users/jeffvandorn/edges_vertices_2_72.pickle', mode='wb') as f:
        pickle.dump(edges_vertices, f, 2)

    now = time.time()-now
    print("finished finding lists of vertices associated with each edge in %.1f seconds" %(now))
    now = time.time()
    #~ print "build a structure that returns the edge numbers associated with one vertex"
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

    with open('/Users/jeffvandorn/verts_edges_2_72.pickle', mode='wb') as f:
        pickle.dump(verts_edges, f, 2)

    now = time.time()-now
    #~ print "time: %.5f" %(now)
    now = time.time()
    verts_OneNeighbor = [ [] for i in range(len(verts))]
    #~ print "verts_OneNeighbor = "
    for ind,vert in enumerate(verts):
        #grab edge from a list of edges for vertex "ind"
        for edge in verts_edges[ind]:
            if selectionOnly:
                #edge is the numbered edge in the entire mesh, needs to be mapped to Selected edge index
                # old way
                #verts_OneNeighbor[ind].extend(edges_vertices[mesh.edges.selected().index(edge)])
                verts_OneNeighbor[ind].extend(edges_vertices[selected_edges.index(edge)])
            else:
                verts_OneNeighbor[ind].extend(edges_vertices[edge])
        if ind < 10:
            print("\nverts_OneNeighbor[%s] before set:\n%s" % (ind, verts_OneNeighbor[ind]))
        if do_sort:
            verts_OneNeighbor[ind] = sorted(set(verts_OneNeighbor[ind]))
        else:
            verts_OneNeighbor[ind] = set(verts_OneNeighbor[ind])
        if ind < 10:
            print("verts_OneNeighbor after set:\n%s" % verts_OneNeighbor[ind])
        verts_OneNeighbor[ind] = [y for y in verts_OneNeighbor[ind] if y != vert.index]
        if ind < 10:
            print("verts_OneNeighbor after assign:\n%s" % verts_OneNeighbor[ind])

        #~ print "vert.index = ", vert.index
        #~ print verts_OneNeighbor[ind]

    with open('/Users/jeffvandorn/verts_one_neighbor_2_72.pickle', mode='wb') as f:
        pickle.dump(verts_OneNeighbor, f, 2)

    now = time.time()-now
    print("finished building structure finding vertices' vertex neighbors in %.1f seconds" %(now))
    #~ now = time.time()
    #~ print "cycle through faces_edges and return 'faces_edge_vert' which returns the two vertices belonging to any face"
    #index 1 - face
    #index 2 - edge
    #index 3 - length 2 for vertex1 or vertex2
    faces_edge_vert = [ [ [] for j in range(4) ] for i in range(len(faces)) ]
    #faces_edge_vert_unique = [ [] for i in range(len(faces)) ]

    for ind,face in enumerate(faces):
        for edge in range(4):
            edge_number = faces_edges[ind][edge]
            if selectionOnly:
                edge_number = selected_edges.index(edge_number)
            faces_edge_vert[ind][edge].extend([edges[edge_number].vertices[0], edges[edge_number].vertices[1]])

        #make a copy with the unique vertex members
        # THIS (faces_edges_vert_unique) IS NOT USED ANYWHERE ELSE
        # commenting out for now
        #[faces_edge_vert_unique[ind].extend(x) for x in faces_edge_vert[ind]]
        #faces_edge_vert_unique[ind] = list(set(faces_edge_vert_unique[ind]))

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
            if ind1==ind2: break #only need to compare ind1 for ind1 < ind2 if you are assigning 2 values after the 'if'
            if faces_edges[ind1] & faces_edges[ind2]: #i.e. if the intersection of these sets is not empty
                faces_neighboring_faces[ind1].append(face2.index)
                faces_neighboring_faces[ind2].append(face1.index)

    #~ print faces_neighboring_faces
    with open('/Users/jeffvandorn/faces_neighboring_faces_272.pickle', mode='wb') as f:
        pickle.dump(faces_neighboring_faces, f, 2)

    now = time.time()-now
    print("finished building face neighbor list in %.1f seconds" %(now))
    now = time.time()
    #now populate list of face's neighboring faces, plus the neighboring faces of those neighboring faces
    #~ print "exclude the 'reference' face"
    #return a unique list
    faces_neighbor2 = [[] for a in range(len(faces))]
    one_neighbor_list = [[] for a in range(len(faces))]
        
    for ind,face in enumerate(faces):

        curr_face_neighbor_list = faces_neighboring_faces[ind]
        for instance in curr_face_neighbor_list:
            if selectionOnly:
                # old way
                #faces_neighbor2[ind].extend(faces_neighboring_faces[mesh.faces.selected().index(instance)])
                faces_neighbor2[ind].extend(faces_neighboring_faces[selected_faces.index(instance)])
            else:
                faces_neighbor2[ind].extend(faces_neighboring_faces[instance])

        #first cull out the current face, then make a copy for next routine, and then find the unique set of face neighbors
        faces_neighbor2[ind] = [val for val in faces_neighbor2[ind] if val != face.index]
        one_neighbor_list[ind] = faces_neighbor2[ind][:] #maybe needs a copy?
        if do_sort:
            faces_neighbor2[ind] = sorted(set(faces_neighbor2[ind]))
        else:
            faces_neighbor2[ind] = set(faces_neighbor2[ind])
        faces_neighbor2[ind] = [val for val in faces_neighbor2[ind]]

        #~ print "face #", face
        #~ print faces_neighbor2[ind]

    with open('/Users/jeffvandorn/one_neighbor_list_2_72.pickle', mode='wb') as f:
        pickle.dump(one_neighbor_list, f, 2)

    now = time.time()-now
    #~ print "time: %.5f" %(now)
    now = time.time()
    #~ print "now take one_neighbor_list, which is an uncompressed list of all the neighboring elements' 1-neighbors"
    #(i.e. the original elements' 2-neighbors)
    #if there is any face # which has incidence of 4, there is a hex associated with the ref. element, the found face #, and the
    #face identities of the 4 elements which articulated with it
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

        if ind < 6:
            print("sixth_faces before set:\n%s" % sixth_faces)
        if do_sort:
            sixth_faces[ind] = sorted(set(sixth_faces[ind]))
        else:
            sixth_faces[ind] = set(sixth_faces[ind])
        if ind < 6:
            print("\nsixth_faces after set:\n%s" % sixth_faces)
        sixth_faces[ind] = [val for val in sixth_faces[ind]]
        if ind < 6:
            print("\nsixth_faces after reassign:\n%s" % sixth_faces)

    with open('/Users/jeffvandorn/sixth_faces_2_72.pickle', mode='wb') as f:
        pickle.dump(sixth_faces, f, 2)

    #~ print "sixth faces = "
    # print "sixth_faces", sixth_faces

    now = time.time()-now
    #~ print "time: %.5f" %(now)
    now = time.time()
    #list must be further culled because this alone doesn't filter out all the hits with 4 neighbors who should be a hex
    #~ print "assemble groups of 6 faces based on the previous routine"
    #first index: global face
    #second index: instance of a hex candidate
    #third index: those 6 faces that MIGHT be a hex
    hex_candidates = [ [ [] for instance in range(len(sixth_faces[a])) ] for a in range(len(faces)) ]
    #~ print "hex candidates = "
    for ind, face in enumerate(faces):
        for instance in range(len(hex_candidates[ind])):
            hex_candidates[ind][instance].append(face.index) #first entry
            hex_candidates[ind][instance].append(sixth_faces[ind][instance]) #sixth face detected is second entry

            #now loop through faces_neighboring_faces for instances of the sixth face in the first face's neighbor list
            for neigh in faces_neighboring_faces[ind]:
                if selectionOnly:
                    if sixth_faces[ind][instance] in faces_neighboring_faces[selected_faces.index(neigh)]:
                        hex_candidates[ind][instance].append(neigh)
                else:
                    if sixth_faces[ind][instance] in faces_neighboring_faces[neigh]:
                        hex_candidates[ind][instance].append(neigh)

    #~ for a in range(len(faces)):
        #~ print hex_candidates[a]

    with open('/Users/jeffvandorn/hex_candidates_2_72.pickle', mode='wb') as f:
        pickle.dump(hex_candidates, f, 2)

    #hex_candidates = pickle.loads(open("/Users/jeffvandorn/hex_candidates.pickle", "rb").read())

    #with open('/Users/jeffvandorn/hex_candidates_2_72_from_249b.pickle', mode='wb') as f:
    #    pickle.dump(hex_candidates, f, 2)

    # hex_candidates holds FACE NUMBERS, len8_hex_candidates holds VERTEX NUMBERS

    #this guy will hold the hex candidates that pass the "8 unique vertex" test
    len8_hex_candidates = []

    #~ print "curr vert list = "

    #variable flag_hexes store a list of tuples that contains the
    #index1 - face number
    #index2 - candidate hex instance number
    flag_hexes = []

    now = time.time()-now
    #~ print "time: %.5f" %(now)
    now = time.time()
    #~ print "now go through and check to see if each of these candidates has exactly 8 unique vertices between all the faces"
    for ind,face in enumerate(faces):

        #instance is a FACE instance
        for candidate_hex_instance in range(len(hex_candidates[ind])):

            curr_vert_list = []
            for face_inst in range(len(hex_candidates[ind][candidate_hex_instance])): #should always be 6
                #get corresponding face structure for this instance
                if selectionOnly:
                    correspond_face = faces[selected_faces.index(hex_candidates[ind][candidate_hex_instance][face_inst])]
                else:
                    correspond_face = faces[hex_candidates[ind][candidate_hex_instance][face_inst]]
                for vert in correspond_face.vertices:
                    curr_vert_list.append(vert)

        #~ print "face = ", face
        #~ print curr_vert_list
        #~ print set(curr_vert_list)

            if len(set(curr_vert_list)) == 8:
                if do_sort:
                    len8_hex_candidates.append(sorted(set(curr_vert_list)))
                else:
                    len8_hex_candidates.append(list(set(curr_vert_list)))
                if selectionOnly:
                    flag_hexes.append((selected_faces.index(face.index),candidate_hex_instance)) #append a tuple
                else:
                    flag_hexes.append((face.index,candidate_hex_instance)) #append a tuple

    #~ print "len8 hex candidates = "
    #~ print len8_hex_candidates

    #~ print "flag hexes = "
    #~ print flag_hexes

    with open('/Users/jeffvandorn/len8_hex_candidates_272.pickle', mode='wb') as f:
        pickle.dump(len8_hex_candidates, f, 2)

    with open('/Users/jeffvandorn/flag_hexes_2_72.pickle', mode='wb') as f:
        pickle.dump(flag_hexes, f, 2)

    #flag_hexes = pickle.loads(open("/Users/jeffvandorn/flag_hexes.pickle", "rb").read())
    #len8_hex_candidates = pickle.loads(open("/Users/jeffvandorn/len8_hex_candidates.pickle", "rb").read())

    unique_hexes = []
    unique_hexes_sorted = []
    unique_face_list = []
    #~ print "unique hexes = "
    counter = 0
    for instance in range(len(len8_hex_candidates)):
        vertlist_instance = len8_hex_candidates[instance]
        if sorted(vertlist_instance) not in unique_hexes_sorted:
            unique_hexes.append(vertlist_instance)
            unique_hexes_sorted.append(sorted(vertlist_instance))
            unique_face_list.append(hex_candidates[flag_hexes[instance][0]][flag_hexes[instance][1]])
            print("\n## unique_face_list[%s] = %s" % (counter, unique_face_list[counter]))
            print("   obtained from: flag_hexes[0]: %s   and flag_hexes[1]: %s" % (flag_hexes[instance][0], flag_hexes[instance][1]))
            counter += 1

    # for instance in range(len(unique_hexes)):
        # print unique_hexes[instance]
    # print "unique face list = "
    # for instance in range(len(unique_hexes)):
        # print unique_face_list[instance]

    with open('/Users/jeffvandorn/unique_face_list_2_72.pickle', mode='wb') as f:
        pickle.dump(unique_face_list, f, 2)

    unique_face_list = pickle.loads(open("/Users/jeffvandorn/unique_face_list.pickle", "rb").read())
    print("unique_face_list[16]: %s" % unique_face_list[16])
    print("unique_face_list[17]: %s" % unique_face_list[17])

    # get material-number of cubes:
    # primary criterion: highest number of faces for material
    # if equal: lower material-number
    HexMat = []
    counter = 0
    for instance in range(len(unique_hexes)):
        HexFacesMat = []
        SumMat = []
        for face in unique_face_list[instance]:
            if selectionOnly:
                HexFacesMat.append(faces[selected_faces.index(face)].material_index)
            else:
                HexFacesMat.append(faces[face].material_index)
        if counter >= 16 and counter <= 17:
            print("\nHexFacesMat[%s] before set\n%s" % (counter, HexFacesMat))
        if do_sort:
            setHexFacesMat = sorted(set(HexFacesMat))
        else:
            setHexFacesMat = set(HexFacesMat)
            setHexFacesMat = list(setHexFacesMat)
        if counter >= 16 and counter <= 17:
            print("HexFacesMat after set\n%s" % setHexFacesMat)

        for mat in setHexFacesMat:
            SumMat.append(sum([val == mat for val in HexFacesMat]))

        HexMat.append(setHexFacesMat[SumMat.index(max(SumMat))])
        counter += 1
    # print "HexMat", HexMat

    with open('/Users/jeffvandorn/hex_mat_2_72.pickle', mode='wb') as f:
        pickle.dump(HexMat, f, 2)
        
    unique_face_list = pickle.loads(open("/Users/jeffvandorn/unique_face_list.pickle", "rb").read())

    #NOW we have our hexes!  now need to get the proper ordering so Continuity understands!
    #get the first face and use its vertices as a reference
    #need to call the native face.verts because then you know the vertex ordering of the face is CW or CCW

    now = time.time()-now
    #~ print "time: %.5f" %(now)
    now = time.time()
    #WARNING we toss out the ordering in variable len8_hex_candidates and number the first 4 according to
    #~ print "first-numbered surface who we pluck from"
    #~ print "ordered hex final = "
    ordered_hex_vertices = [[] for b in range(len(unique_hexes))]
    for loop_num in range(len(unique_hexes)):
        curr_verts = unique_hexes[loop_num]
        curr_faces = unique_face_list[loop_num]
        curr_first_face = curr_faces[0]

        #get four verts belonging to the first-numbered surface, then rearrange them so that it is not CW/CCW
        first_four_verts = []
        if selectionOnly:
            for vert in faces[selected_faces.index(curr_first_face)].vertices:
                first_four_verts.append(vert)            
        else:
            for vert in faces[curr_first_face].vertices:
                first_four_verts.append(vert)
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
            pos5_vertex = [x for x in remaining_four_verts if x in verts_OneNeighbor[selected_verts.index(vertex_position1)]]
            pos6_vertex = [y for y in remaining_four_verts if y in verts_OneNeighbor[selected_verts.index(vertex_position2)]]
            pos7_vertex = [z for z in remaining_four_verts if z in verts_OneNeighbor[selected_verts.index(vertex_position3)]]
            pos8_vertex = [a for a in remaining_four_verts if a in verts_OneNeighbor[selected_verts.index(vertex_position4)]]
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

    with open('/Users/jeffvandorn/ordered_hex_2_72.pickle', mode='wb') as f:
        pickle.dump(ordered_hex_vertices, f, 2)

    now = time.time()-now
    #~ print "time: %.5f" %(now)
    print ("Found %s hexes! \n\nFinding contiguous material regions..." % len(ordered_hex_vertices))
    
    #find out whether there are any non-intesecting regions that belong to the same labeled topology region
    #and separate them
    HexMatNew = contiguous_regions(list(ordered_hex_vertices), HexMat)
    
    print("Harmonizing topology regions...")
    
    #find verts and put them into an array
    verts_array = np.zeros([len(mesh.vertices),3])
    for ind, vert in enumerate(mesh.vertices):
        verts_array[ind,:] = vert.co
    
    #harmonize hexes
    ordered_hex_vertices = harmonize_topo(ordered_hex_vertices,verts_array,HexMatNew)
    
    #save the elems and mats so you don't have to find every time that you do something
    if not selectionOnly:
        cache_data(ordered_hex_vertices,HexMatNew)
    
    print("\nFINISHED FIND HEX!")

    return ordered_hex_vertices, HexMat

def delete_all(mesh_name):
    ''' Deletes the mesh with the given name'''
    bpy.ops.object.mode_set(mode='OBJECT')
    bpy.ops.object.select_by_type(type='MESH')
    bpy.ops.object.delete(use_global=False)
    for mesh in bpy.data.meshes:
        if mesh.name == mesh_name:
            bpy.data.meshes.remove(mesh)

def write_pickle(file_id, **data):
    return pickle.dump(data, file_id)
