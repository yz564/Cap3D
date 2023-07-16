import os
import numpy as np
from netCDF4 import Dataset
import glob
import datetime

# Global variable MatRotY
MatRotY = None

def ConvertCubit_g_2_elx_V8(nstr='*', RotYDeg=0):
    if nstr == '*':
        nstr = '*'
        RotYDeg = 0

    if RotYDeg is None:
        RotYDeg = 0

    global MatRotY
    beta = RotYDeg * np.pi / 180
    cosb = np.cos(beta)
    sinb = np.sin(beta)
    MatRotY = np.array([[cosb, 0, sinb], [0, 1, 0], [-sinb, 0, cosb]])

    gtype = '.g'
    etype = '.elx'
    gv = glob.glob(nstr + gtype)
    elxv = glob.glob(nstr + etype)

    M = len(gv)
    for m in range(M):
        gname = os.path.splitext(os.path.basename(gv[m]))[0]
        gv[m] = gname

    N = len(elxv)
    for m in range(N):
        elxname = os.path.splitext(os.path.basename(elxv[m]))[0]
        elxv[m] = elxname

    nelem_nnod_natt = []
    for m in range(M):
        gname = gv[m]
        gexist = gname in elxv
        if not gexist:
            ncg2elxV3(gname + gtype, gname + etype, False)
            print('\t{} {}: converted successfully!'.format(gname, etype))
        else:
            print('\t{} {}: NO need to be converted!'.format(gname, etype))

def get_NodXYZ(ncid):
    var_names = ncid.variables.keys()
    xyz_vec = None

    # Look for the variable containing coordinate data
    if 'coord' in var_names:
        xyz_vec = ncid.variables['coord'][:]
    elif 'coordx' in var_names and 'coordy' in var_names and 'coordz' in var_names:
        x_vec = ncid.variables['coordx'][:]
        y_vec = ncid.variables['coordy'][:]
        z_vec = ncid.variables['coordz'][:]
        xyz_vec = np.column_stack((x_vec, y_vec, z_vec))

    return xyz_vec


def ncg2elxV3(ncg_file, file_out, pf=True):
    global MatRotY  # Global variable should be defined here

    ncid = Dataset(ncg_file, 'r')

    # Get dimensions
    num_node = ncid.dimensions['num_nodes'].size
    num_elem = ncid.dimensions['num_elem'].size
    num_blk = ncid.dimensions['num_el_blk'].size

    # Get element information
    nelem_nnod_natt = []
    for cnt in range(1, num_blk + 1):
        blk_el_str = 'num_el_in_blk' + str(cnt)
        blk_nn_str = 'num_nod_per_el' + str(cnt)
        blk_at_str = 'num_att_in_blk' + str(cnt)

        nelem = ncid.dimensions[blk_el_str].size
        nnod = ncid.dimensions[blk_nn_str].size
        natt = ncid.dimensions[blk_at_str].size

        nelem_nnod_natt.append([nelem, nnod, natt])

    if sum([elem[0] for elem in nelem_nnod_natt]) != num_elem or any([elem[2] != 1 for elem in nelem_nnod_natt]):
        print('ERROR: in subdomain element numbers or attributes !!!')
        return

    # Read attribute and connect data
    att_connect_eltype = []
    elem_type_vec = []
    for cnt in range(1, num_blk + 1):
        cstr = 'attrib' + str(cnt)
        att_vec = ncid.variables[cstr][:]

        nelem = nelem_nnod_natt[cnt - 1][0]
        natt = nelem_nnod_natt[cnt - 1][2]

        #ConnectAttrib = np.transpose(np.reshape(att_vec, (natt, nelem)))
        ConnectAttrib = np.reshape(att_vec, (nelem,natt))
        cstr = 'connect' + str(cnt)
        con_vec = ncid.variables[cstr][:]
        nnod = nelem_nnod_natt[cnt - 1][1]
        #Connect = np.transpose(np.reshape(con_vec, (nnod, nelem)))
        Connect = np.reshape(con_vec, (nelem,nnod))
        att_connect_eltype.append([Connect, ConnectAttrib])

        elem_typ_str = ncid.variables[cstr].elem_type
        L = len(elem_typ_str)
        while elem_typ_str[L - 1].isdigit():
            elem_typ_str = elem_typ_str[:L - 1]
            L = len(elem_typ_str)

        if elem_typ_str == 'SHELL':
            if len(ConnectAttrib[0]) == 2:
                shell_attrib = ConnectAttrib[0, 0] * 1000 + ConnectAttrib[0, 1]
            else:
                shell_attrib = ConnectAttrib[0, 0]
            elem_type_vec.append('QUAD4')
        elif elem_typ_str == 'TETRA':
            elem_type_vec.append('TETRA10')
        elif elem_typ_str == 'SHELL':
            elem_type_vec.append('QUAD9')
        elif elem_typ_str == 'TRI':
            elem_type_vec.append('TRI3')
        else:
            print('ERROR: Element Type needs to be modified !!!')
            return

    # Read coordinate data
    xyz_vec = get_NodXYZ(ncid)
    if xyz_vec is None:
        print('ERROR: Could not find coordinate data in the NetCDF file!')
        return

    num_node = len(xyz_vec)
    NodXYZ = np.zeros((num_node, 3))

    # Assign the coordinate data
    if xyz_vec.shape[1] == 3:
        NodXYZ[:, :] = xyz_vec
    else:
        print('ERROR: Incompatible coordinate data shape in the NetCDF file!')
        return

    # Rotate the coordinates
    #if MatRotY is not None:
    #    NodXYZ = np.dot(NodXYZ, MatRotY.T)

    # Write data to the output file
    affix = '.elx'
    if file_out.endswith(affix):
        file_out = file_out[:-len(affix)]
    with open(file_out + affix, 'w') as fwid:
        fwid.write(f'#Date: {str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))}, by Yang Zhong\n')
        fwid.write(f'{num_blk:4d}{num_elem:8d}{num_node:8d} \t#num_blk, num_elem, num_node\n')
        fwid.write(f'{nelem_nnod_natt[0][1]:4d}{elem_type_vec[0]:>12s} \t#num_nod_per_elem, elem_type\n')

        wfmt = ''
        for m in range(nelem_nnod_natt[0][1]):
            wfmt += '%0d '
        wfmt += '\r\n'

        for cnt in range(1, num_blk + 1):
            fwid.write(f'{nelem_nnod_natt[cnt - 1][0]:8d}{cnt:4d}{int(att_connect_eltype[cnt - 1][1][0, 0]):8d} \t#num_elem_blk, blk_id, atrib_blk\n')
            for m in range(nelem_nnod_natt[cnt - 1][0]):
                fwid.write(wfmt % tuple(att_connect_eltype[cnt - 1][0][m, :]))

        for m in range(num_node):
            v = NodXYZ[m, :]
            w = np.dot(MatRotY, v.T)
            fwid.write('%.16G %.16G %.16G \n' % (w[0], w[1], w[2]))

    ncid.close()

    if pf:
        print(f'Mesh converted successfully: {file_out}{affix} !!!')


ConvertCubit_g_2_elx_V8()

