# get the data
# I'm guessing we need to create a numpy array of the stl data:
import numpy as np
from scipy.linalg import logm

def render_dt_form(filename, coords, es, vs):
    '''
    renderDTform Write out a diffusion tensor rendering form for Continuity.
    renderDTform(FILE,COORDS,ES,VS) writes a text file named FILE
    with the data point coordinates (COORDS), sorted eigenvalues (ES),
    and sorted eigenvectors (VS). It assumes that the eigenvalues are
    sorted in descending order (i.e. es = [es_largest es_mid es_smallest];
    vs = [v_primary; v_secondary; vs_tertiary]; the Continuity data
    renderer reverses this order and instead expects es to be [smallest
    largest]. Simply change the order of the eigenvalues in the output to match,
    keeping the original order of the eigenvectors; this is the
    explanation for the rearrangement of the eigenvectors for the rendering
    form.
    CVillongco
    October 2011
    '''

    print("Number of data points:%s" % len(coords))
    ones_array = np.ones((len(coords.squeeze())))

    labels = [num+1 for num in range(len(coords.squeeze()))]
    np_labels = np.array(labels, int)

    fmt_str = '%.6e\t%d\t%.6e\t%d\t%.6e\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d'

    data_hdr = 'coord1_val\tcoord1_weight\tcoord2_val\tcoord2_weight\tcoord3_val\tcoord3_weight\tevec11_val\tevec11_weight\tevec12_val\tevec12_weight\tevec13_val\tevec13_weight\tevec21_val\tevec21_weight\tevec22_val\tevec22_weight\tevec23_val\tevec23_weight\tevec31_val\tevec31_weight\tevec32_val\tevec32_weight\tevec33_val\tevec33_weight\teval1_val\teval1_weight\teval2_val\teval2_weight\teval3_val\teval3_weight\tData\n'

    import pdb;pdb.set_trace()
    np.savetxt('test_render.out', 
        (np.transpose([coords[:,0].squeeze(), ones_array, coords[:,1].squeeze(), ones_array, coords[:,2].squeeze(),
              vs[0,2,:].squeeze(), ones_array, vs[0,1,:].squeeze(), ones_array, vs[0,0,:].squeeze(),
              vs[1,2,:].squeeze(), ones_array, vs[1,1,:].squeeze(), ones_array, vs[1,0,:].squeeze(),
              vs[2,2,:].squeeze(), ones_array, vs[2,1,:].squeeze(), ones_array, vs[2,0,:].squeeze(), 
              es[:,2].squeeze(), ones_array, es[:,1].squeeze(), ones_array, es[:,0].squeeze(), ones_array, np_labels])), 
              fmt=fmt_str, delimiter='\t', header=data_hdr, comments='')


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

    np.savetxt('test.out', (np.transpose([coords.squeeze()[:,0], coords.squeeze()[:,1], coords.squeeze()[:,2],
              data[0,0,:].squeeze(), data[1,1,:].squeeze(), data[2,2,:].squeeze(),
              #data[0,1,:].squeeze(), data[0,2,:].squeeze(), data[1,2,:].squeeze(), np.array(labels)])), 
              data[0,1,:].squeeze(), data[0,2,:].squeeze(), data[1,2,:].squeeze(), np_labels])), 
              fmt=fmt_str, delimiter='\t', header=data_hdr, comments='')

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

def make_test_data():
    coords = np.zeros([5,3])
    dt = np.zeros([3,3,5])
    es = np.zeros([5,3])
    vs = np.zeros([3,3,5])

    dt[0,0,0] = 0.001286000000000   
    dt[0,1,0] = 0.000116000000000   
    dt[0,2,0] = 0.000136000000000
    dt[1,0,0] = 0.000116000000000   
    dt[1,1,0] = 0.001525000000000  
    dt[1,2,0] = -0.000371000000000
    dt[2,0,0] = 0.000136000000000  
    dt[2,1,0] = -0.000371000000000   
    dt[2,2,0] = 0.001563000000000

    dt[0,0,1] = 0.848000000000000e-03  
    dt[0,1,1] = -0.116000000000000e-03  
    dt[0,2,1] = -0.015000000000000e-03
    dt[1,0,1] = -0.116000000000000e-03   
    dt[1,1,1] = 0.370000000000000e-03  
    dt[1,2,1] = -0.065000000000000e-03
    dt[2,0,1] = -0.015000000000000e-03  
    dt[2,1,1] = -0.065000000000000e-03   
    dt[2,2,1] = 0.927000000000000e-03

    dt[0,0,2] = 0.000574000000000   
    dt[0,1,2] = 0.000032000000000  
    dt[0,2,2] = -0.000002000000000
    dt[1,0,2] = 0.000032000000000   
    dt[1,1,2] = 0.001011000000000   
    dt[1,2,2] = 0.000029000000000
    dt[2,0,2] = -0.000002000000000   
    dt[2,1,2] = 0.000029000000000   
    dt[2,2,2] = 0.001393000000000

    dt[0,0,3] = 0.000634000000000  
    dt[0,1,3] = -0.000304000000000  
    dt[0,2,3] = -0.000140000000000
    dt[1,0,3] = -0.000304000000000   
    dt[1,1,3] = 0.000264000000000   
    dt[1,2,3] = 0.000158000000000
    dt[2,0,3] = -0.000140000000000   
    dt[2,1,3] = 0.000158000000000   
    dt[2,2,3] = 0.001543000000000

    dt[0,0,4] = 0.000960000000000   
    dt[0,1,4] = 0.000006000000000  
    dt[0,2,4] = -0.000099000000000
    dt[1,0,4] = 0.000006000000000   
    dt[1,1,4] = 0.000710000000000   
    dt[1,2,4] = 0.000009000000000
    dt[2,0,4] = -0.000099000000000   
    dt[2,1,4] = 0.000009000000000   
    dt[2,2,4] = 0.001038000000000

    coords[0,:] = [-0.180000000000000,   0.820000000000000,   0.840000000000000]
    coords[1,:] = [-0.180000000000000,   0.880000000000000,   0.840000000000000]
    coords[2,:] = [-0.180000000000000,   0.820000000000000,   0.860000000000000]
    coords[3,:] = [-0.180000000000000,   0.840000000000000,   0.860000000000000]
    coords[4,:] = [-0.180000000000000,   0.860000000000000,   0.860000000000000]

    es[0,:] = [0.001043000000000,   0.001415000000000,   0.001916000000000]
    es[1,:] = [0.000336000000000,   0.000875000000000,   0.000935000000000]
    es[2,:] = [0.000571000000000,   0.001012000000000,   0.001396000000000]
    es[3,:] = [0,   0.000755000000000,   0.001597000000000]
    es[4,:] = [0.000709000000000,   0.000894000000000,   0.001106000000000]

    vs[0,0,0] = -0.590430000000000   
    vs[0,1,0] =0.806428000000000  
    vs[0,2,0] = -0.032655000000000
    vs[1,0,0] =0.576200000000000   
    vs[1,1,0] = 0.449508000000000   
    vs[1,2,0] = 0.682595000000000
    vs[2,0,0] = 0.565143000000000   
    vs[2,1,0] = 0.384209000000000  
    vs[2,2,0] =-0.730067000000000
   
    vs[0,0,1] = -0.223094000000000   
    vs[0,1,1] =0.974576000000000  
    vs[0,2,1] =-0.020764000000000
    vs[1,0,1] =-0.968260000000000  
    vs[1,1,1] =-0.224011000000000  
    vs[1,2,1] =-0.110870000000000
    vs[2,0,1] =-0.112702000000000  
    vs[2,1,1] =-0.004629000000000   
    vs[2,2,1] =0.993618000000000
  
    vs[0,0,2] =  -0.997271000000000   
    vs[0,1,2] = 0.073827000000000  
    vs[0,2,2] = -0.000571000000000
    vs[1,0,2] = 0.073657000000000   
    vs[1,1,2] = 0.994384000000000  
    vs[1,2,2] = -0.075998000000000
    vs[2,0,2] = -0.005043000000000  
    vs[2,1,2] = -0.075833000000000  
    vs[2,2,2] = -0.997108000000000

    vs[0,0,3] = -0.477006000000000   
    vs[0,1,3] =0.857901000000000  
    vs[0,2,3] =-0.190971000000000
    vs[1,0,3] =-0.877494000000000  
    vs[1,1,3] =-0.452580000000000   
    vs[1,2,3] =0.158668000000000
    vs[2,0,3] =0.049692000000000   
    vs[2,1,3] =0.243262000000000   
    vs[2,2,3] =0.968687000000000

    vs[0,0,4] =-0.037370000000000   
    vs[0,1,4] =0.826254000000000   
    vs[0,2,4] =0.562057000000000
    vs[1,0,4] =0.998597000000000   
    vs[1,1,4] =0.051988000000000  
    vs[1,2,4] =-0.010031000000000
    vs[2,0,4] =-0.037508000000000   
    vs[2,1,4] =0.560893000000000  
    vs[2,2,4] =-0.827038000000000

    return dt, coords, es, vs

def read_data_file(filename):
    data = np.loadtxt(filename)

    dt_data = np.zeros([3,3,len(data)]);
    coords = np.zeros([len(data),3]);
    es = np.zeros([len(data),3])
    vf = np.zeros([3,3,len(data)])

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


def read_off_file(filename):
    with open(filename, 'r') as infile:
        linecount = 0
        array_count = 0
        for line in infile:
            if linecount == 0:
                pass
            if linecount == 1:
                if "\t" in line:
                    spacer = "\t"
                else:
                    spacer = " "
                numlines = int(line.split(spacer)[0])
                contents = np.zeros([numlines, 4], float)
            if linecount > 1:
                contents[array_count,0] = float(line.split(spacer)[0])
                contents[array_count,1] = float(line.split(spacer)[1])
                contents[array_count,2] = float(line.split(spacer)[2])
                contents[array_count,3] = 1
                array_count += 1

            linecount += 1

        #return np.matrix(contents)
        return contents


def compute_affine_xform(aligned_filename, unaligned_filename, data_filename):
    '''
    stl_unaligned = bpy.data.objects[stl_unaligned_name]
    #stl_unaligned_np =  np.zeros([len(stl_unaligned.vertices),4])

    # this is going to hurt
    stl_list = [item.co for item in stl_unaligned.data.vertices]
    '''

    aligned_data = read_off_file(aligned_filename)
    unaligned_data = read_off_file(unaligned_filename)
    #dt, coords, es, vs = read_data_file(data_filename)
    dt, coords, es, vs = make_test_data()
    print("data read...")

    # Defines the coordinates to the transformed data that fits to
    # Continuity standards
    #coords_seg = coords_from_data_file
    coords_seg = unaligned_data

    # number of data point
    #num_data = len(aligned_data)
    num_data = 5

    # compute 4x4 transformation
    transformation1 = np.dot(unaligned_data.conj().T, unaligned_data)
    transformation2 = np.dot(unaligned_data.conj().T, aligned_data)
    transformation = np.linalg.solve(transformation1, transformation2)
    print("computed transformation:")
    print(transformation)

    C = transformation[0:3, 0:3]

    # Apply transformation to scanner coordinates
    coords_temp = coords_seg
    coords_temp[:,3] = 1
    coords_r = np.dot(coords_temp, transformation)
    coords_r = np.delete(coords_r, 3, axis = 1)

    print('applied transformation')

    ####################################################################
    # rotated diffusion tensors
    dt_r = np.zeros(dt.shape)
    # rotated diffusion tensors, log transforms
    dt_rl = np.zeros(dt.shape)
    # roated diffusion tensors eigenvectors
    vpr = np.zeros(dt.shape)

    indsrem = []
    k = 0

    # Apply transformation to DT data
    for index in range(num_data):
        if index % 10000 == 0:
            print("%d tensors processed" % index)

        # DT rotation
        a = np.dot(C.conj().T, dt[:,:,index])
        dt_r[:,:,index] = np.dot(a, C)

        # finds the eigenvectors (v) and eigenvalues (e) of the DT data
        v, e = eig_s(dt_r[:,:,index])

        # Eigenvectors
        vpr[:,:,index] = v

        # Eigenvalues
        # don't need epr(j,:)  = es_scale(j,:);
        # don't need DT_rscale(:,:,j) = vpr(:,:,j)*diag(epr(j,:))*vpr(:,:,j)';

        # Compute matrix logarithm of DT
        dt_rl[:,:,index] = logm(dt_r[:,:,index])
        if np.isnan(dt_rl).any() or np.isinf(dt_rl).any():
            # mark data indices with nonreal tensors
            indsrem.append(index);

    #################################################################### 
    # Finds the index of the original coordinates at a certain slice in 
    # the model.  Because of signifcant figures, must first round the 
    # coordinates before finding slice at which to render 

    # find median image slice of DT data in scanner reference frametensors.
    median_image = np.median(coords[:,1])
    ind_slice = np.where(coords[:,1] == median_image)

    # initialize matrix to track DT data indices after removing bad voxels
    inds_slice2 = np.zeros((num_data)) 

    # DT data indices belonging to median slice get set to 1
    inds_slice2[ind_slice] = 1 

    # DT data indices belonging to bad voxels get reset back to 0
    inds_slice2[indsrem] = 0  

    # find DT data indices that are in median slice AND are not bad voxels
    ind_slice3 = np.where(inds_slice2==1)

    inds_dat = np.ones((num_data))
    inds_dat[indsrem] = 0
    inds_dat2 = np.where(inds_dat==1)

    #################################################################### 
    # save all DT data to Continuity data form
    data_dt_form('Sham_20140605_DTform_test_jvd.txt',
                coords_r[inds_dat2,:],
                dt_rl[:,:,inds_dat2])

    # Write out Continuity data file to check rotated positions
    #filename = 'TAC_20140617_rendertensors-med_test_jvd';
    #print('Data file written to %s with %d points' % (filename,len(coords_r(inds_dat2,:)))
    #renderDTform([file '-Cont.txt'],coords_r(ind_slice3,:),epr(ind_slice3,:),vpr(:,:,ind_slice3))
    #renderDTform([file '-Scanner.txt'],coords(ind_slice3,:),es(ind_slice3,:),vs(:,:,ind_slice3))
    filename = 'render_sample'
    import pdb;pdb.set_trace()
    render_dt_form(filename+'_cont.txt',
        coords_r[ind_slice3,:],
        es[ind_slice3,:],
        vpr[:,:,ind_slice3])

    render_dt_form(filename+'_scanner.txt',
        coords[ind_slice3,:],
        es[ind_slice3,:],
        vs[:,:,ind_slice3])



unaligned_filename = '/Users/jeffvandorn/work/my_stuff/NBCR_summer_schoo_2015/DTfitting/TAC_20140617_DTdata_coords.off'
aligned_filename = '/Users/jeffvandorn/work/my_stuff/NBCR_summer_schoo_2015/DTfitting/TAC_20140617_2_Transrot_DTdata_coords.off'
data_filename = '/Users/jeffvandorn/work/my_stuff/NBCR_summer_schoo_2015/6-10_fin_dump_masked.txt'
compute_affine_xform(aligned_filename, unaligned_filename, data_filename)
