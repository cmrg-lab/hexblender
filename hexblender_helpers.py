import bpy
import blf
import bgl
import time
import json
import bmesh
import math
import numpy as np
import pickle
from mathutils import Vector
from random import random
from bpy_extras.view3d_utils import location_3d_to_region_2d as loc3d2d
from hexblender.harmonize_hex import harmonize_topo
from hexblender.hermite_element_derivatives import hermite_element_derivatives
from hexblender.hex_interp_subdiv import hex_interp_subdiv
from hexblender.hex_interp_subdiv_no_priorities import hex_interp_subdiv_no_priorities
from hexblender.regularize_elements_pullback_vectorized import regularize_elements_vect
from hexblender.export_16pt_hermite import subdivide_quad_16pt_hermite

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
        bpy.ops.object.mode_set(mode='EDIT')
        # XXX: I have to do this command twice for it to actually
        # set the mode to EDIT, why?
        bpy.ops.object.mode_set(mode='EDIT')
        b_mesh = bmesh.from_edit_mesh(mesh)
        return mesh, b_mesh
    else:
        print("Error: Unable to create bmesh!")


def reset_mode(orig_mode):
    if orig_mode != 'OBJECT':
        print("Changing back into mode: %s" % orig_mode)
        bpy.ops.object.mode_set(mode=orig_mode)


def get_material_data(mesh = None):
    orig_mode = None
    
    try:
        faces = mesh.faces
    except Exception as msg:
        mesh, orig_mode = get_active_mesh()
        if mesh == None:
            return None
        else:
            faces = mesh.polygons

    matList = [[] for val in faces]

    for face in faces:
        matList[face.index] = face.material_index

    if orig_mode is not None:
       reset_mode(orig_mode)

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
    if len(mesh.cached_data_json) > 0:
        try:
            cached_data = json.loads(mesh.cached_data_json)
        except Exception as msg:
            print("Error reloading hex data: %s" % msg)
            cached_data = {}
    else:
        print("No cached data found.")
        cached_data = {}

    for ind, hex_vert in enumerate(ordered_hex_vertices):
        cached_data[str(ind)] = hex_vert

    cached_data['num_elems'] = len(ordered_hex_vertices)
    cached_data['mats'] = hex_mat_new

    # Store the cached data in our property, which seems to be the only?
    # way to have the data serialized in a blend file
    mesh.cached_data_json = json.dumps(cached_data)
    mesh.update()
    
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


def set_mat_data_hex(matList, new_cube_faces, mesh = None, subdivided=True, 
                     existing_mats = None):

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
    faces_in_verts_dict = {}
    counter = 0
    for face in faces:
        c_face = [0]*4
        c_face[0] = int(face.vertices[0])
        c_face[1] = int(face.vertices[1])
        c_face[2] = int(face.vertices[2])
        c_face[3] = int(face.vertices[3])
        cSet_face = sorted(c_face)
        facesInVerts.append(cSet_face)
        faces_in_verts_dict["%s" % cSet_face] = counter
        counter += 1

    start = time.time()
    for i,face in enumerate(new_cube_faces):
        c_face = sorted(face)
        # if c_face in facesInVerts:
        c_Mat = nMatList[int(i/6)]
        # else:
        if(faces[faces_in_verts_dict["%s" % c_face]].material_index == 0 or 
            (faces[faces_in_verts_dict["%s" % c_face]].material_index > c_Mat and c_Mat != 0)):
            faces[faces_in_verts_dict["%s" % c_face]].material_index = c_Mat

    print("cube_faces loop took: %s" % (time.time() - start))

    num_mats = len(set(matList))

    # Create new materials and colors if they don't already exist
    # otherwise, use the exisitng materials
    if existing_mats == None:
        for i in range(num_mats):
            mat = bpy.data.materials.new("Material_%i" % i)
            mat.diffuse_color = random(), random(), random()
            mesh.materials.append(mat)
    else:
        # This will most often be empty, but just being careful
        mat_name_list = [obj.name for obj in mesh.materials]
        for mat in existing_mats:
            # We don't want duplicates
            if mat.name not in mat_name_list:
                mesh.materials.append(mat)


def add_new_data(verts_new, edges_keys_index, faces_keys_index, cubes_verts_index, newObj = False, mesh_name = None):

    def check_for_dups(new_cube_faces, unique_cubes):
        sorted_cube = sorted(new_cube_faces[-1])
        try:
            unique_cubes["%s" % sorted_cube] == 1
            new_cube_faces.pop()
        except KeyError as msg:
            unique_cubes["%s" % sorted_cube] = 1

    full_new_cube_faces = []
    new_cube_faces = []
    # takes the clockwise order of verts
    unique_cubes = {}
    for cube in cubes_verts_index:
        new_cube_faces.append([int(value) for value in cube[:4]])
        full_new_cube_faces.append(new_cube_faces[-1])
        check_for_dups(new_cube_faces, unique_cubes)
        
        new_cube_faces.append([int(value) for value in cube[4:]])
        full_new_cube_faces.append(new_cube_faces[-1])
        check_for_dups(new_cube_faces, unique_cubes)

        for i in range(4):
            new_cube_faces.append([int(cube[i]),
                                   int(cube[(i+1)%4]),
                                   int(cube[4+(i+1)%4]),
                                   int(cube[(i+4)])])
            full_new_cube_faces.append(new_cube_faces[-1])
            check_for_dups(new_cube_faces, unique_cubes)

    # add all new faces
    for face_key in faces_keys_index:
        new_cube_faces.append(face_key)

    if newObj:
        if mesh_name == None:
            mesh_name = 'pickledObject'
        mesh = bpy.data.meshes.new(mesh_name)
        mesh.from_pydata(verts_new,
                         [], #edges
                         new_cube_faces, #faces
                        )

        mesh.update()

        # create the new mesh and make it active 
        obj = bpy.data.objects.new(mesh_name, mesh)
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

    return full_new_cube_faces

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
    start_time = time.time() 
    mesh, orig_mode = get_active_mesh()

    if selectionOnly:
        if verts is None:
            print("!!!! verts is None!!!!")
            verts = [v for v in mesh.vertices if v.select]
            if len(verts) == 0:
                print("ERROR: nothing is selected!")
                return

        # Get selected edge and polygon objects
        edges = [e for e in mesh.edges if e.select]
        faces = [p for p in mesh.polygons if p.select]
            
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

    print("got mesh in %f" % (time.time() - start_time))

    print("begin finding hexes...\n")
    now = time.time()
    # now = time.time()-now
    #~ print "time: %.5f" %(now)
    now = time.time()
    #get a list of (4) edges that correspond to each face
    faces_edges = [[] for a in range(len(faces))]
    #~ print "faces' edges"
    now = time.time()
    
    # make a lookup for the mesh edge_keys
    mesh_edge_keys_dict = {}
    for ind, keys in enumerate(mesh.edge_keys):
        mesh_edge_keys_dict[keys] = ind

    for ind,face in enumerate(faces):
        for edge in range(4):
            faces_edges[ind].append(mesh_edge_keys_dict[face.edge_keys[edge]])
    print("got face_edges in %f" % (time.time() - now))

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

    print("got edge_vertices in %f" % (time.time() - now))
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

    print("got verts_edges in %f" % (time.time() - now))

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
        verts_OneNeighbor[ind] = set(verts_OneNeighbor[ind])
        verts_OneNeighbor[ind] = [y for y in verts_OneNeighbor[ind] if y != vert.index]

        #~ print "vert.index = ", vert.index
        #~ print verts_OneNeighbor[ind]

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
        faces_neighbor2[ind] = set(faces_neighbor2[ind])
        faces_neighbor2[ind] = [val for val in faces_neighbor2[ind]]

        #~ print "face #", face
        #~ print faces_neighbor2[ind]

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

        sixth_faces[ind] = set(sixth_faces[ind])
        sixth_faces[ind] = [val for val in sixth_faces[ind]]

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
                len8_hex_candidates.append(list(set(curr_vert_list)))
                if selectionOnly:
                    flag_hexes.append((selected_faces.index(face.index),candidate_hex_instance)) #append a tuple
                else:
                    flag_hexes.append((face.index,candidate_hex_instance)) #append a tuple

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
                HexFacesMat.append(faces[selected_faces.index(face)].material_index)
            else:
                HexFacesMat.append(faces[face].material_index)
        setHexFacesMat = set(HexFacesMat)
        setHexFacesMat = list(setHexFacesMat)

        for mat in setHexFacesMat:
            SumMat.append(sum([val == mat for val in HexFacesMat]))

        HexMat.append(setHexFacesMat[SumMat.index(max(SumMat))])
    # print "HexMat", HexMat

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
        cache_data(ordered_hex_vertices, HexMatNew)
    
    print("\nFINISHED FIND HEX!")

    reset_mode(orig_mode)
    print("find hex took: %f" % (time.time() - start_time))
    return ordered_hex_vertices, HexMat


def delete_all(mesh_name):
    ''' Deletes the mesh object with the given name'''

    print("!! Deleting mesh object: %s" % mesh_name)
    bpy.ops.object.mode_set(mode='OBJECT')
    for ob in bpy.context.scene.objects:
        if ob.type == 'MESH' and ob.name == mesh_name:
            ob.select = True
        else:
            ob.select = False
            
    bpy.ops.object.delete()

    try:
        bpy.data.meshes.remove(bpy.data.meshes[mesh_name])
    except Exception as msg:
        print("WARNING removing mesh: %s" % msg)


def delete_mesh():
    ''' Deletes the elements of the mesh'''
    
    mesh, b_mesh = get_active_mesh_and_bmesh()
    bmesh.ops.delete(b_mesh, geom=b_mesh.faces, context=6)
    bmesh.ops.delete(b_mesh, geom=b_mesh.verts, context=6)

    # Show the updates in the viewport
    bmesh.update_edit_mesh(mesh, True)     
    

def write_pickle(file_id, **data):
    return pickle.dump(data, file_id)


def get_priorities(data_string):
    '''This will turn the given string into a list'''

    tmp = eval(data_string)

    data = []
    if type(tmp) is tuple:
        for item in tmp:
            data.append(item)
    else:
        data = [tmp]

    return data


def get_cached_data(mesh, contiguous = False):
    ''' Retrieve the cached data '''
    bpy.types.Mesh.cached_data_json = bpy.props.StringProperty()

    # See if we have the original mesh element and material data
    if len(mesh.cached_data_json) > 0:
        try:
            cubes, mat_list = reload_hex(json.loads(mesh.cached_data_json))
            print("Reloaded cached hex data!")
        except Exception as msg:
            print("Error reloading hex data: %s" % msg)
    else:
        print("No cached data found, finding cubes...")
        cubes, mat_list = find_hex(selectionOnly=False)

        if contiguous:
            # Preseving comment from Blender 2.49b HexBlender
            # there might be a little bug here not properly using
            # contiguousRegions() that can be fixed by tabbing in and out
            # of edit mode. I didn't want to add the normal command Window.
            # EditMode(1) here to fix this but someone can try it :)
            hex_mat_new = contiguous_regions(cubes, mat_list)
            return cubes, hex_mat_new

    return cubes, mat_list

def write_out_vertices(filename = None):
    try:
        mesh, orig_mode = get_active_mesh()
        verts = mesh.vertices
    except Exception as msg:
        print("Error getting vertices: %" % msg)
        return {'CANCELLED'}, None

    # === Header ===
    if hex_prop.export_vert_weights:
        # export verts with weights
        vert_string = 'Coords_1_val\tCoords_2_val\tCoords_3_val'

        # from panels, get names of selected vertex groups/weight paint layer(s) to export,
        # and add field variable header(s)
        group_name = hex_prop.group_name
        field_name = hex_prop.field_name
        # need to be able to specify which field user wants to assign from weight paint layer
        # nodal_fields = ['FieldVec1_Var1','FieldVec1_Var2','FieldVec1_Var3',
        #                'FieldVec2_Var4','FieldVec2_Var5','FieldVec2_Var6',
        #                'FieldVec3_Var7','FieldVec3_Var8','FieldVec3_Var9',
        #                'FieldVec4_Var10','FieldVec4_Var11','FieldVec4_Var12',
        #                'FieldVec5_Var13','FieldVec5_Var14','FieldVec5_Var15']
        # field_names  = ['Field 1','Field 2','Field 3',
        #                'Field 4','Field 5','Field 6',
        #                'Field 7','Field 8','Field 9',
        #                'Field 10','Field 11','Field 12']
        field_dict = {'Fiber angle':'FibAng_Fiber','Transverse angle':'FibAng_Trans','Sheet angle':'FibAng_Sheet',
                      'Field 1':'FieldVec1_Var1','Field 2':'FieldVec1_Var2','Field 3':'FieldVec1_Var3',
                      'Field 4':'FieldVec2_Var4','Field 5':'FieldVec2_Var5','Field 6':'FieldVec2_Var6',
                      'Field 7':'FieldVec3_Var7','Field 8':'FieldVec3_Var8','Field 9':'FieldVec3_Var9',
                      'Field 10'':FieldVec4_Var10','Field 11':'FieldVec4_Var11','Field 12':'FieldVec4_Var12',
                      'Field 13'':FieldVec5_Var13','Field 14':'FieldVec5_Var14','Field 15':'FieldVec5_Var15'}
        # 
        for k in field_name:
            vert_string = ('\t%s_val' % field_dict[k])
        vert_string += ('\tLabel\tNodes\n')
    else:
        # export verts only
        vert_string = 'Coords_1_val\tCoords_2_val\tCoords_3_val\tLabel\tNodes\n'

    # === Vertex List ===
    if hex_prop.export_vert_weights:
        # CV: get name of selected vertex group to export
        for i, v in enumerate(verts):
            vert_string += ('%.6f\t%.6f\t%.6f' % (tuple(v.co)))
            # CV: Need to figure a way to get weights of selected vertex group (weight paint layer).
            # When weight painting, Blender only includes painted vertices in the vertex group.
            # If a vertex does not belong to the weight paint group, it's default value should 
            # be set to 0 when written to the output string.
            # An alternative is to create a vertex group (weight paint layer) that includes all
            # vertices; all vertices need to be painted/initialized with zero weight so that
            # when we call a vertex group by name, we can be confident that all vertices are
            # included in the vertex group, and the output string will be easier to write in
            # vert_string += ('%.6f\t%d\t%d\n' % (v[vertex_group_names].weight,i+1,i+1))

            vert_string += ('%.6f\t%d\t%d\n' % (v.weight,i+1,i+1))
    else:
        for i, v in enumerate(verts):
            vert_string += ('%.6f\t%.6f\t%.6f' % (tuple(v.co)))
            vert_string += ('\t%d\t%d\n' % (i+1,i+1))

    if filename is not None:
        myfile = open(filename, 'w')
        myfile.write(vert_string)
        myfile.close()
        print("Successfully exported %s" % filename) 

    reset_mode(orig_mode)
    return {'FINISHED'}, vert_string

########### MD 1/12/2016 ###########
def compute_write_out_vertex_weights(context, filename = None):
    try:
        # identifies vertex and weights in the selected group
        mesh, orig_mode = get_active_mesh()
        verts = mesh.vertices
        obj = bpy.context.active_object
        bpy.ops.object.mode_set(mode='OBJECT',toggle=True)
        mesh = obj.data
        selVerts = [v for v in mesh.vertices]
        indexVal = obj.vertex_groups.active_index
        weights = []
        vert_weights = []
        for v in selVerts:
            for n in v.groups:
                if n.group == indexVal:
                    weights.append([v,n.weight])
                    vert_weights.append([weights[-1][0].index,weights[-1][1]])
    except Exception as msg:
        print("Error getting vertex weights: %" % msg)
        return {'CANCELLED'}, None
    
    try:
        hex_scene = context.scene.hexblender
    except Exception as msg:
        print("ERROR: Unable to get hexblender scene!  %s" % msg)
        return {'CANCELLED'}, None

    scale_f = hex_scene.hexblender_properties.vert_weight_scalar
    vert_weights_scaled = []
    if scale_f != 1.00:
        for v in vert_weights:
            vert_weights_scaled.append([v[0], scale_f*v[1]])
    else:
        vert_weights_scaled = vert_weights

    # === Header ===
    # vertweight_string = 'Vertex\tWeight\n'
    # vert_string = 'Coords_1_val\tCoords_2_val\tCoords_3_val\tLabel\tNodes\n'

    # === Vertex List ===
    for i in vert_weights_scaled:
        vertweight_string += ('%d\t%.6f\n' % (i[0],i[1]))
        vert_string += ('%.6f\t%.6f\t%.6f' % (tuple(v.co)))
        vertweight_string += ('%.6f\t%.6f\t%.6f' % (tuple(v.co)))
        vertweight_string += ('\t%d\t%d\n' % (i+1,i+1))

    if filename is not None:
        myfile = open(filename, 'w')
        myfile.write(vertweight_string)
        myfile.close()
        print("Successfully exported %s" % filename) 

    reset_mode(orig_mode)
    return {'FINISHED'}, vertweight_string
####################################

def write_out_3d_elements(filename = None):
    try:
        mesh, orig_mode = get_active_mesh()
        cubes, hex_mat_new = get_cached_data(mesh, contiguous = True)
    except Exception as msg:
        print("Error getting elements: %s" % msg)
        return {'CANCELLED'}, None

    # Write headers
    elem_string = 'Node_0_Val\tNode_1_Val\tNode_2_Val\tNode_3_Val\tNode_4_Val\tNode_5_val\tNode_6_Val\tNode_7_Val\tLabel\tElement\n'
    
    # Write elements
    for num_elem, elem_verts in enumerate(cubes):
        for v in elem_verts:
            elem_string += '%i\t' % (v+1)
        elem_string += 'region%d\t%d\n' %(hex_mat_new[num_elem]+1,num_elem+1)

    if filename is not None:
        myfile = open(filename, 'w')
        myfile.write(elem_string)
        myfile.close()
        print("Successfully exported %s" % filename) 

    reset_mode(orig_mode)
    return {'FINISHED'}, elem_string


def compute_write_hermite_3d_derivs(context, filename = None):

    def write_hermite_deriv(filename, Bx, By, Bz):
        #~ print ordered_hex_vertices[loop_num]

        # write headers
        deriv_string = 'u\tdu_dxi1\tdu_dxi2\td2u_dxi1xi2\tdu_dxi3\td2u_dxi2dxi3\td2u_dxi1dxi3\td3u_dxi\tCoord\tNode\tElement\n'

        # write elements
        for i in range(len(Bx)):
            for j in range(8):
                for k in range(8):
                    deriv_string += '%f\t' % Bx[i,j,k]
                deriv_string += '%i\t%i\t%i\n' %(1,j+1,i+1)
                for k in range(8):
                    deriv_string += '%f\t' % By[i,j,k]
                deriv_string += '%i\t%i\t%i\n' %(2,j+1,i+1)
                for k in range(8):
                    deriv_string += '%f\t' % Bz[i,j,k]
                deriv_string += '%i\t%i\t%i\n' %(3,j+1,i+1)

        if filename is not None:
            myfile = open(filename, 'w')
            myfile.write(deriv_string)
            myfile.close()
            print("Successfully exported %s" % filename) 
            return None

        return deriv_string

    try:
        hex_scene = context.scene.hexblender
    except Exception as msg:
        print("ERROR: Unable to get hexblender scene!  %s" % msg)
        return

    itt = hex_scene.hexblender_properties.tri_cubic_iters
    write_out_sfm = True
    loadAdjustedNodes = False

    # We have to do subdivision if filename is None, as that indicates
    # that we are exporting a bundle file and that requires doing the 
    # subdivision
    subdiv_type = int(hex_scene.hexblender_properties.export_subdiv_types)
    if subdiv_type > 0 or filename is None:
        useExistingNodes = False
    else:
        useExistingNodes = True
        print("Using existing mesh for tricubic derivs.")

    if loadAdjustedNodes:
        dirName = 'C:\Continuity\Models\Simulation92/'
        cubes_normal,dummy = ReloadHex()
        eNN = []
        for hex in cubes_normal:
            eNN.append([hex[0],hex[1],hex[3],hex[2],hex[4],hex[5],hex[7],hex[6]])
        eNN = np.array(eNN)
        nNN = np.load(dirName + 'adjustedNodes.npy')
    
    if not loadAdjustedNodes and not useExistingNodes:
        mesh, orig_mode = get_active_mesh()
        cubes, matList = get_cached_data(mesh)

        if not write_out_sfm:
            if len(mesh.cached_data_json) > 0:
                hexMatNew = ContiguousRegions(cubes, matList)
                orig_matlist = copy.copy(hexMatNew)
            else:
                orig_matlist = copy.copy(matList)
            orig_cubes = copy.copy(cubes)

        cubes,tmp0,n = get_boundary_data(cubes)
        cubes = np.array(cubes)
        n = np.array(n)

        # by default, default is the totality of all material numbers, 
        # and the subdivision nor regularization are done piecewise
        # this here is an internal switch

        # Now getting priorities from the GUI, but of course
        # people can still overwrite/change them here should they
        # want.  
        interpPriorities = get_priorities(
            hex_scene.hexblender_properties.interp_priorities)

        # Example settings
        #interpPriorities = [np.unique(matList).tolist()] # default
        #interpPriorities = [[0,1,2,3,4,5,6,7,8,9,12,13,14,15],[10,11],[16]]
        #interpPriorities = [[0,1,2,3,4,5,6,7,8,9,10]]
        
        # Priority groups for the element regularization. 
        # By default, uses all material numbers simultaneously
        regPriorities = get_priorities(
            hex_scene.hexblender_properties.reg_priorities)

        # Example settings
        #regPriorities = [np.unique(matList).tolist()] # default
        #regPriorities = [[0,1,2,3,4,5,6,7,8,9,12,13,14,15],[10,11],[16]]
        #regPriorities = [[0,1,2,3,4,5,6,7,8,9,10]]

        splinePriorities = get_priorities(
            hex_scene.hexblender_properties.spline_priorities)

        # Example settings
        #splinePriorities = [[0],[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12],[13],[14],[15],[16]]
        #splinePriorities = [[0],[1],[2],[3],[4],[5],[6],[7],[8],[9],[10]]
        #splinePriorities = [np.unique(matList).tolist()] # default
        
        MatListSubdivOnce = []
        for val in matList:
            MatListSubdivOnce.extend([val, val, val, val, val, val, val, val])
        MatListSubdivTwice = []
        for val in MatListSubdivOnce:
            MatListSubdivTwice.extend([val, val, val, val, val, val, val, val])
        
        # options for HexInterpSubdiv:
        # 1) break into topology regions, 
        # 2) thin plate spline mapping
        print("starting first subdivision!")

        if subdiv_type == 1:
            print("Using the following priorities for subdiv/deriv calcs")
            print("interpPriorities: %s" % interpPriorities)
            print("regPriorities: %s" % regPriorities)
            print("splinePriorities: %s" % splinePriorities)

            eN, nN = hex_interp_subdiv(cubes,
                                       n,
                                       MatList=matList,
                                       priorities=interpPriorities, 
                                       thinPlateMapping=True, 
                                       thinPlateRegions=splinePriorities)
            print("starting second subdivision!")
            eNN, nNN = hex_interp_subdiv(eN,
                                         nN,
                                         MatList=MatListSubdivOnce,
                                         priorities=interpPriorities)
        elif subdiv_type == 2:
            print("Not using priorities for subdiv/deriv calcs")
            eN, nN = hex_interp_subdiv_no_priorities(cubes,
                                       n,
                                       MatList=matList,
                                       priorities=interpPriorities, 
                                       thinPlateMapping=True, 
                                       thinPlateRegions=splinePriorities)
            print("starting second subdivision!")
            eNN, nNN = hex_interp_subdiv_no_priorities(eN,
                                         nN,
                                         MatList=MatListSubdivOnce,
                                         priorities=interpPriorities)
        else:
            print("!! Invalid subdivision type: %s" % subdiv_type)
            return
        
        #cache result (eNN and MatListSubdivTwice)
        cubes_normal, newfaces, newverts = get_boundary_data(eNN)
        # I am perplexed as to why, but I have to turn cubes_normal into 
        # an array (then I turn it back into a list) 
        # or else the IDProperty type chokes 
        #cache_data(np.array(cubes_normal).tolist(),MatListSubdivTwice)

        if itt > 0:
            # regularize the mesh piecewise. If 
            # priorities = [np.unique(HexMat).tolist()] (the default), 
            # this executes once flip the priorities list so that the 
            # last list executes first, and the first list executes last
            regPriorities.reverse()
            
            for priority_group in regPriorities:
                
                #pick out the hexes for the base mesh
                currHexesInMatlGrp = np.array([ind for ind,val in enumerate(matList) if val in priority_group])
                cubes_group = np.array(cubes).copy()
                cubes_group = cubes_group[np.array(currHexesInMatlGrp),:].tolist()
                
                allVertInds = np.unique(np.array(cubes_group[:]).flatten())
                allVertIndsSet = set(allVertInds)
                isVertIncluded = np.array([vert in allVertIndsSet for vert in range(n.shape[0])])
                vertSelectionMap = np.cumsum(isVertIncluded)-1
                
                verts_array = np.zeros([np.shape(allVertInds)[0],3])
                for ind,val in enumerate(verts_array):
                    verts_array[ind,0] = n[allVertInds[ind],0]
                    verts_array[ind,1] = n[allVertInds[ind],1]
                    verts_array[ind,2] = n[allVertInds[ind],2]
                
                cubes_reduced = vertSelectionMap[np.array(cubes_group)]
                cubes_trans = np.zeros(np.shape(cubes_reduced),'int32')
                for ind,cube in enumerate(cubes_reduced):
                    cubes_trans[ind,0] = cube[0]
                    cubes_trans[ind,1] = cube[1]
                    cubes_trans[ind,2] = cube[3]
                    cubes_trans[ind,3] = cube[2]
                    cubes_trans[ind,4] = cube[4]
                    cubes_trans[ind,5] = cube[5]
                    cubes_trans[ind,6] = cube[7]
                    cubes_trans[ind,7] = cube[6]
                
                #pick out the hexes for the twice-subdivided mesh
                currHexesInMatlGrpTwSub = np.array([ind for ind,val in enumerate(MatListSubdivTwice) if val in priority_group])
                cubes_groupTwSub = eNN.copy()
                cubes_groupTwSub = cubes_groupTwSub[np.array(currHexesInMatlGrpTwSub),:].tolist()
                
                allVertIndsTwSub = np.unique(np.array(cubes_groupTwSub[:]).flatten())
                allVertIndsTwSubSet = set(allVertIndsTwSub)
                isVertIncludedTwSub = np.array([vert in allVertIndsTwSubSet for vert in range(np.shape(nNN)[0])])
                vertSelectionMapTwSub = np.cumsum(isVertIncludedTwSub)-1
                
                verts_arrayTwSub = np.zeros([np.shape(allVertIndsTwSub)[0],3])
                for ind,val in enumerate(verts_arrayTwSub):
                    verts_arrayTwSub[ind,0] = nNN[allVertIndsTwSub[ind],0]
                    verts_arrayTwSub[ind,1] = nNN[allVertIndsTwSub[ind],1]
                    verts_arrayTwSub[ind,2] = nNN[allVertIndsTwSub[ind],2]
                
                cubes_reducedTwSub = vertSelectionMapTwSub[np.array(cubes_groupTwSub)]
                
                nNNReg_reduced =\
                   regularize_elements_vect(cubes_reduced,
                                            verts_array,
                                            cubes_reducedTwSub,
                                            verts_arrayTwSub,
                                            itt,
                                            immobilizeRidges=False)
    
                for ind,vert_num in enumerate(allVertIndsTwSub):
                    nNN[vert_num,0] = nNNReg_reduced[ind,0]
                    nNN[vert_num,1] = nNNReg_reduced[ind,1]
                    nNN[vert_num,2] = nNNReg_reduced[ind,2]
        
        # Do we keep coarse mesh?
        existing_mats = get_mats_used_unused_sorted(mesh)
        mesh_name = bpy.context.active_object.name
        if hex_scene.hexblender_properties.delete_orig_mesh:
            delete_all(mesh_name)
            print("Deleted coarse mesh. \nCreating subdivided mesh.")

        new_cube_faces = add_new_data(nNN, [], [], eNN, newObj=True, mesh_name = mesh_name)
        set_mat_data_hex(MatListSubdivOnce, new_cube_faces, existing_mats=existing_mats)
        #set_mat_data_hex(MatListSubdivOnce, new_cube_faces)

        # show the new mesh and cache it
        print("Updating view.")
        context.scene.objects.active = context.scene.objects.active 
        cache_data(np.array(cubes_normal).tolist(),MatListSubdivTwice)

    if (useExistingNodes):
        mesh, orig_mode = get_active_mesh()

        cubes_normal, dummy = get_cached_data(mesh)
        eNN = []
        for hex in cubes_normal:
            eNN.append([hex[0],hex[1],hex[3],hex[2],hex[4],hex[5],hex[7],hex[6]])
        eNN = np.array(eNN)

        verts_export = []
        for vert in mesh.vertices:
            verts_export.append([vert.co.x, vert.co.y, vert.co.z])
        nNN = np.array(verts_export)

    Bx, By, Bz = hermite_element_derivatives(eNN,nNN)
    # print "Bx\n", Bx[0]
    # print "Bx3\n", Bx[0,:,3]
    # print "Bx4\n", Bx[0,:,4]

    Bx = np.swapaxes(Bx,2,0)
    Bx = np.swapaxes(Bx,2,1)
    By = np.swapaxes(By,2,0)
    By = np.swapaxes(By,2,1)
    Bz = np.swapaxes(Bz,2,0)
    Bz = np.swapaxes(Bz,2,1)

    for i in range(len(Bx)):
        tmp = Bx[i,:,3].copy()
        Bx[i,:,3] = Bx[i,:,4].copy()
        Bx[i,:,4] = tmp.copy()
        tmp = Bx[i,:,5].copy()
        Bx[i,:,5] = Bx[i,:,6].copy()
        Bx[i,:,6] = tmp.copy()

        tmp = Bx[i,2,:].copy()
        Bx[i,2,:] = Bx[i,3,:].copy()
        Bx[i,3,:] = tmp.copy()
        tmp = Bx[i,6,:].copy()
        Bx[i,6,:] = Bx[i,7,:].copy()
        Bx[i,7,:] = tmp.copy()

    for i in range(len(By)):
        tmp = By[i,:,3].copy()
        By[i,:,3] = By[i,:,4].copy()
        By[i,:,4] = tmp.copy()
        tmp = By[i,:,5].copy()
        By[i,:,5] = By[i,:,6].copy()
        By[i,:,6] = tmp.copy()

        tmp = By[i,2,:].copy()
        By[i,2,:] = By[i,3,:].copy()
        By[i,3,:] = tmp.copy()
        tmp = By[i,6,:].copy()
        By[i,6,:] = By[i,7,:].copy()
        By[i,7,:] = tmp.copy()

    for i in range(len(Bz)):
        tmp = Bz[i,:,3].copy()
        Bz[i,:,3] = Bz[i,:,4].copy()
        Bz[i,:,4] = tmp.copy()
        tmp = Bz[i,:,5].copy()
        Bz[i,:,5] = Bz[i,:,6].copy()
        Bz[i,:,6] = tmp.copy()

        tmp = Bz[i,2,:].copy()
        Bz[i,2,:] = Bz[i,3,:].copy()
        Bz[i,3,:] = tmp.copy()
        tmp = Bz[i,6,:].copy()
        Bz[i,6,:] = Bz[i,7,:].copy()
        Bz[i,7,:] = tmp.copy()

    deriv_string = write_hermite_deriv(filename, Bx, By, Bz)
    reset_mode(orig_mode)
    return {'FINISHED'}, deriv_string

def write_out_2d_elements(filename = None):
    # WARNING - HARMONIZE QUAD IS NOT CALLED HERE, AND CONTIGUOUS 
    # REGIONS WILL NOT, IN GENERAL, SHARE TOPOLOGIES
    # I DO THIS BECAUSE I WANT TO BE CONSISTENT WITH THE NODAL 
    # POSITIONS READ BY BICUBIC HERMITE EXPORT RIGHT NOW
    
    mesh, orig_mode = get_active_mesh()
    faces = mesh.polygons
    verts = mesh.vertices

    matGroup = []
    faces_verts = []
    for face in faces:
        face_verts = []
        for i in range(4):
            face_verts.append(face.vertices[i])
        faces_verts.append(face_verts)
        matGroup.append(face.material_index)
    
    matGroup = contiguous_regions(faces_verts, matGroup)

    #write headers
    elem_string = 'Node_0_Val\tNode_1_Val\tNode_2_Val\tNode_3_Val\tLabel\tElement\n'

    #write elements
    for i in range(len(faces_verts)):
        elem_string += '%i\t%i\t%i\t%i\tregion%i\t%i\n' % (faces_verts[i][0]+1, faces_verts[i][1]+1, faces_verts[i][3]+1, faces_verts[i][2]+1, matGroup[i]+1, i+1)

    if filename is not None:
        fid = open(filename, 'w')
        fid.write(elem_string)
        fid.close()
        print("Successfully exported %s" % filename)

    reset_mode(orig_mode)
    return {'FINISHED'}, elem_string

def compute_write_hermite_2d_derivs(context, filename = None):

    def write_hermite(filename, Bx, By, Bz):
        #~ print ordered_hex_vertices[loop_num]

        #write headers
        deriv_string = 'u\tdu_dxi1\tdu_dxi2\td2u_dxi1xi2\tCoord\tNode\tElement\n'

        Hmat_pts = [(0,0),(2,0),(0,2),(2,2),(1,0),(3,0),(1,2),(3,2),(0,1),(2,1),(0,3),(2,3),(1,1),(3,1),(1,3),(3,3)]

        #write elements
        for i in range(len(Bx)):    #element loop
            for j in range(4):     #local node loop

                for k in range(4):  # 4 x 4 hermite matrix loop
                    deriv_string += '%f\t' % Bx[i, Hmat_pts[4*j+k][0], Hmat_pts[4*j+k][1]]
                deriv_string += '%i\t%i\t%i\n' %(1,j+1,i+1)

                for k in range(4):
                    deriv_string += '%f\t' % By[i,Hmat_pts[4*j+k][0],Hmat_pts[4*j+k][1]]
                deriv_string += '%i\t%i\t%i\n' %(2,j+1,i+1)

                for k in range(4):
                    deriv_string += '%f\t' % Bz[i,Hmat_pts[4*j+k][0],Hmat_pts[4*j+k][1]]
                deriv_string += '%i\t%i\t%i\n' %(3,j+1,i+1)

        if filename is not None:
            dfile = open(filename, 'w')
            dfile.write(deriv_string)
            dfile.close()
            print("Successfully exported %s" % filename) 
            return None
        else:
            return deriv_string

    try:
        hex_scene = context.scene.hexblender
    except Exception as msg:
        print("ERROR: Unable to get hexblender scene!  %s" % msg)
        return

    mesh, orig_mode = get_active_mesh()

    tmp,f,n = get_boundary_data([])
    f = np.array(f)
    n = np.array(n)

    Bx, By, Bz = subdivide_quad_16pt_hermite(f,n)
    deriv_string = write_hermite(filename, Bx, By, Bz)
    reset_mode(orig_mode)
    return {'FINISHED'}, deriv_string

def get_object_data(object_name):
    if object_name in bpy.data.objects:
        object_data = bpy.data.objects[object_name].data
    else:
        print("ERROR: Unable to find object named: %s" % object_name)
        print("Known objects are:")
        print(bpy.data.objects.keys())
        return None

    # OK, we have a known object, we want to create a numpy array
    # of it's coordinates and an additional column with ones in it
    object_verts = [[item.co.x, item.co.y, item.co.z, 1] for item in object_data.vertices]

    return np.array(object_verts)


def eig_s(data):
    l, v = np.linalg.eig(data)

    # add sort function to make conventionally sorted vectors and
    # values (Continuity's convention); 
    # smallest (l(1)) to largest (l(3)) in eigenvalues
    # need to sort negative eigenvalues in descending order and positive in ascending order
    
    i_s = l.argsort()
    l_s = np.sort(abs(l))
    ld = np.diag(l)
        
    l_s = np.diag(ld[i_s])
    v_s = v[:,i_s]

    return v_s, l_s


def read_tensor_data_file(filename):
    data = np.loadtxt(filename)

    dt_data = np.zeros([3,3,len(data)]);
    coords = np.zeros([len(data),3]);
    es = np.zeros([len(data),3])
    vs = np.zeros([3,3,len(data)])

    # Image resolution (cm) (used in Continuity for anatomical fitting mask))
    xres = 0.937504;
    yres = 0.937504;
    zres = 0.937504;

    size_x = 100;
    size_y = 100;
    size_z = 100;

    dt_data[0,0,:] = data[:,4] 
    dt_data[0,1,:] = data[:,5] 
    dt_data[0,2,:] = data[:,7] 
    dt_data[1,0,:] = data[:,5] 
    dt_data[1,1,:] = data[:,6]
    dt_data[1,2,:] = data[:,8] 
    dt_data[2,0,:] = data[:,7] 
    dt_data[2,1,:] = data[:,8] 
    dt_data[2,2,:] = data[:,9]

    coords[:,0] = data[:,1] * xres; 
    coords[:,1] = data[:,2] * yres; 
    coords[:,2] = data[:,3] * zres;

    es[:,0] = data[:,12] 
    es[:,1] = data[:,11] 
    es[:,2] = data[:,10]

    vs[0,2,:] = data[:,13] 
    vs[1,2,:] = data[:,14]
    vs[2,2,:] = data[:,15]
    vs[0,1,:] = data[:,16] 
    vs[1,1,:] = data[:,17]
    vs[2,1,:] = data[:,18]
    vs[0,0,:] = data[:,19] 
    vs[1,0,:] = data[:,20]
    vs[2,0,:] = data[:,21]

    return dt_data, coords, es, vs

def data_dt_form(filename, coords, data):
    '''
    dataDTform(FILE,COORDS,D) writes an output text file named FILE 
    with the diffusion tensors D and their spatial coordinates into a 
    file suitably formatted for rendering in Continuity.
    Input arguments:
       FILE is the path to the output file.
       COORDS is a nx3 array of the x,y,z Cartesian coordinates of 
       n diffusion tensors.
       D is a 3x3xn array of the diffusion tensor data to be rendered.

    The output FILE is to be imported into Continuity using 
    Mesh->Render->Raw Diffusion Tensors.
    In the 'Render Raw Tensors Form' popup menu:
       1. Set the path to FILE in 'File Name:'.
       2. Select 'Data from: FILE' radio button.
       3. Ignore 'X Slice(s)', 'Y Slice(s)', 'Z slice(s)', 
          and 'Tensor thinning factor'.
       4. Set eigenvalue scaling to appropriate value. 
          For COORDS in cm, use ~0.1; for COORDS in mm, use ~0.01; 
          for COORDS in um, use ~0.001.
       5. 'Superquadric squareness factor' controls the sharpness 
          of the edges of the rendered glyphs; lower values produce 
          glyphs with soft/round edges; higher values produce glyphs 
          with sharper and squarer edges.
       6. Check 'Normalize eigenvalues'.
       7. Set 'Coloring' to 'Fractional anisotropy'. This setting 
          actually colors the tensors by the x-component 
          (red Continuity axis) of the primary eigenvector to 
           illustrate the longitudinal/horizontal orientation of the glyph.
       8. Click 'OK' to render.
    '''
    print("Number of data point: %s" % len(coords.squeeze()))

    labels = [num+1 for num in range(len(coords.squeeze()))]
    np_labels = np.array(labels, int)

    data_hdr = 'coord_1_val\tcoord_2_val\tcoord_3_val\tdxx_val\tdyy_val\tdzz_val\tdxy_val\tdxz_val\tdyz_val\tData'
    fmt_str = '%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%d'

    # was test.out
    np.savetxt(filename, (np.transpose([coords.squeeze()[:,0], coords.squeeze()[:,1], coords.squeeze()[:,2],
              data[0,0,:].squeeze(), data[1,1,:].squeeze(), data[2,2,:].squeeze(),
              data[0,1,:].squeeze(), data[0,2,:].squeeze(), data[1,2,:].squeeze(), np_labels])), 
              fmt=fmt_str, delimiter='\t', header=data_hdr, comments='')

    print("Wrote out data to: %s" % filename)
    
def adjust_list(in_list, x, y):
    return [[old_x + x, old_y + y] for (old_x, old_y) in in_list]


def generate_points(width, height):
    amp = 5  # radius fillet

    width += 2
    height += 4
    width = ((width/2) - amp) + 2
    height -= (2*amp)

    pos_list, final_list = [], []

    n_points = 12
    seg_angle = 2 * math.pi / n_points
    for i in range(n_points + 1):
        angle = i * seg_angle
        x = math.cos(angle) * amp
        y = math.sin(angle) * amp
        pos_list.append([x, -y])

    w_list, h_list = [1, -1, -1, 1], [-1, -1, 1, 1]
    slice_list = [[i, i+4] for i in range(0, n_points, 3)]

    for idx, (start, end) in enumerate(slice_list):
        point_array = pos_list[start:end]
        w = width * w_list[idx]
        h = height * h_list[idx]
        final_list += adjust_list(point_array, w, h)

    return final_list


def get_points(index, point_dict):
    '''
    index:   string representation of the index number
    returns: rounded rect point_list used for background.

    the neat thing about this is if a width has been calculated once, it
    is stored in a dict and used if another polygon is saught with that width.
    '''
    width, height = blf.dimensions(0, index)
    if not (width in point_dict):
        point_dict[width] = generate_points(width, height)

    return point_dict[width]


# calculate locations and store them as ID property in the mesh
def draw_callback_px(self, context, point_dict, mesh):
    # polling
    if context.mode != "EDIT_MESH":
        return

    # get screen information
    region = context.region
    rv3d = context.space_data.region_3d
    this_object = context.active_object
    matrix_world = this_object.matrix_world

    text_height = 13
    blf.size(0, text_height, 72)

    def draw_index(rgb, index, coord):

        vector3d = matrix_world * coord
        x, y = loc3d2d(region, rv3d, vector3d)

        index = str(index)
        polyline = get_points(index, point_dict)

        ''' draw polygon '''
        bgl.glColor4f(0.103, 0.2, 0.2, 0.2)
        bgl.glBegin(bgl.GL_POLYGON)
        for pointx, pointy in polyline:
            bgl.glVertex2f(pointx+x, pointy+y)
        bgl.glEnd()

        ''' draw text '''
        txt_width, txt_height = blf.dimensions(0, index)
        bgl.glColor3f(*rgb)
        blf.position(0, x - (txt_width / 2), y - (txt_height / 2), 0)
        blf.draw(0, index)

    vert_idx_color = (1.0, 1.0, 1.0)
    edge_idx_color = (1.0, 1.0, 0.0)
    face_idx_color = (1.0, 0.8, 0.8)
    hex_idx_color = (1.0, 0.6, 0.0)

    bm = bmesh.from_edit_mesh(mesh)

    try:
        hex_scene = context.scene.hexblender
    except Exception as msg:
        print("ERROR: Unable to get hexblender scene!  %s" % msg)
        return {'FINISHED'}

    hex_prop = hex_scene.hexblender_properties

    if hex_prop.live_mode:
        mesh.update()

    if hex_prop.display_vert_index:
        for v in bm.verts:
            if not v.hide and (v.select or not hex_prop.display_sel_only):
                ## CoDEmanx: bm.verts.index_update()?
                draw_index(vert_idx_color, v.index+1, v.co.to_4d())

    if hex_prop.display_edge_index:
        for e in bm.edges:
            if not e.hide and (e.select or not hex_prop.display_sel_only):
                v1 = e.verts[0].co
                v2 = e.verts[1].co
                loc = v1 + ((v2 - v1) / 2)
                draw_index(edge_idx_color, e.index+1, loc.to_4d())

    if hex_prop.display_face_index:
        for f in bm.faces:
            if not f.hide and (f.select or not hex_prop.display_sel_only):
                draw_index(
                    face_idx_color, f.index+1, f.calc_center_median().to_4d())

    if hex_prop.display_hex_index:
        bpy.types.Mesh.cached_data_json = bpy.props.StringProperty()

        # see if we have already cached the data
        if len(mesh.cached_data_json) > 0:
            try:
                cubes, mat_list = reload_hex(json.loads(mesh.cached_data_json))
            except Exception as msg:
                print("Error reloading hex data: %s" % msg)
        else:
            print("finding cubes...")
            cubes, mat_list = find_hex(selectionOnly=False)
            bm = bmesh.from_edit_mesh(mesh)
            bm.verts.ensure_lookup_table()

        # get the average x,y,z for the verts of this hex
        for index, elem in enumerate(cubes):
            plot = []
            x_loc = 0
            y_loc = 0
            z_loc = 0
            for vert in elem:
                x_loc += bm.verts[vert].co.x
                y_loc += bm.verts[vert].co.y
                z_loc += bm.verts[vert].co.z

                # handle the selection only case, based on vertices
                # if all 8 verts were selected, show the hex index
                if hex_prop.display_sel_only:
                    if bm.verts[vert].select:
                        plot.append(True)
                    else:
                        plot.append(False)
                else:
                    plot = [True]*8

            if False not in plot and len(plot) == 8:
                vect = Vector((x_loc/len(elem), y_loc/len(elem), z_loc/len(elem), 1))

                draw_index(hex_idx_color, index+1, vect)

def parse_list(listStr):   #utility function parses string to listStr delimiting',' and '-'
    ParseListError = "Error: Undefined element flag!"
    sc = listStr.split(',')
    for (i,s) in enumerate(sc):
        sc[i] = s.strip()

    sd = []
    for i in sc:
        if i.find('-') != -1:
            temp = i.split('-')
            min = int(temp[0])
            max = int(temp[1])
            for i in range(min, max + 1):
                sd.append(i)
        else:
            try:
                sd.append(int(i))
            except: 
                raise ParseListError
    return sd

def get_mats_defined(mesh):
    '''
    Simple getter that will return the currently defined material OBJECTs
    in the supplied mesh
    NOTE: This does not mean that all the returns materials are actually
    used in the object (i.e., associated with any faces/etc), just that it
    has been defined.
    '''

    return [obj for obj in mesh.materials]

def get_mats_used_unused_sorted(mesh):
    '''
    This will return all the currently defined material OBJECTs, but will
    sort them by usage.  Meaning, all the matls at the beginning of the list
    are actually used by elements within the mesh, whereas the ones at 
    the end are not and have simply been defined.
    '''

    # Go through all the faces and get the index of the material used
    # I'm certainly open to finding other ways to know which elements
    # of a mesh belong to a given material, but this is the best I 
    # could find at the moment
    now = time.time()
    used_mat_index = [face.material_index for face in mesh.polygons]
    used_mat_index_unique = list(set(used_mat_index))

    used_mat_names = []
    # Get the names of the used materials
    for index in used_mat_index_unique:
        used_mat_names.append(mesh.materials[index].name)

    # Get all the defined materials
    defined_mats = get_mats_defined(mesh)

    # Make the lists of used/unused matls and then return the joined lists
    used_mat_objs = []
    unused_mat_objs = []
    for mat in defined_mats:
        if mat.name in used_mat_names:
            used_mat_objs.append(mat)
        else:
            unused_mat_objs.append(mat) 

    used_mat_objs.extend(unused_mat_objs)
    print("Found used materials in: %f" % ((time.time() - now)))

    return used_mat_objs
