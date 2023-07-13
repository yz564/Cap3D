import os
import numpy as np
from netCDF4 import Dataset

def ConvertCubit_g_2_elx_V8(nstr='*', RotYDeg=0):
    if nstr == '*' and RotYDeg == 0:
        nstr = '*'
        RotYDeg = 0
    
    global MatRotY
    beta = np.deg2rad(RotYDeg)
    cosb = np.cos(beta)
    sinb = np.sin(beta)
    MatRotY = np.array([[cosb, 0, sinb], [0, 1, 0], [-sinb, 0, cosb]])

    gtype = '.g'
    etype = '.elx'
    gv = [file for file in os.listdir() if file.endswith(gtype) and nstr in file]
    elxv = [file for file in os.listdir() if file.endswith(etype) and nstr in file]

    M = len(gv)
    for m in range(M):
        gname = gv[m]
        L = len(gname) - len(gtype)
        gnam = gname[:L]
        gv[m] = gnam

    N = len(elxv)

    print('Convert *.g to *.elx:')
    for m in range(M):
        gname = gv[m]
        gexist = False
        for n in range(N):
            gexist = (gname + etype) == elxv[n]
            if gexist:
                break
        if not gexist:
            ncg2elxV3(gname + gtype, gname + etype, False)
            print(f'\t{gname + etype}: converted successfully!')
        else:
            print(f'\t{gname + etype}: NO need to be converted!')

def ncg2elxV3(ncg_file, file_out, pf=True):
    global MatRotY

    finfo = Dataset(ncg_file, 'r')
    dim_num = len(finfo.dimensions)
    for dim in finfo.dimensions.values():
        tch = dim.name
        if tch == 'num_nodes':
            num_node = dim.size
        elif tch == 'num_elem':
            num_elem = dim.size
        elif tch == 'num_el_blk':
            num_blk = dim.size

    nelem_nnod_natt = np.zeros((num_blk, 3))
    for cnt in range(1, num_blk + 1):
        blk_el_str = f'num_el_in_blk{cnt}'
        blk_nn_str = f'num_nod_per_el{cnt}'
        blk_at_str = f'num_att_in_blk{cnt}'

        for dim in finfo.dimensions.values():
            tch = dim.name
            if tch == blk_el_str:
                nelem = dim.size
            elif tch == blk_nn_str:
                nnod = dim.size
            elif tch == blk_at_str:
                natt = dim.size

        nelem_nnod_natt[cnt-1, :] = [nelem, nnod, natt]

    if (np.sum(nelem_nnod_natt[:, 0]) != num_elem) or (not np.all(nelem_nnod_natt[:, 2] == 1)):
        print('ERROR: in subdomain element numbers or attributes !!!')
        return

    att_connect_eltype = [[f'attrib{cnt}', f'connect{cnt}'] for cnt in range(1, num_blk + 1)]
    elem_type_vec = [None] * num_blk
    for cnt in range(1, num_blk + 1):
        att_connect_eltype[cnt-1][0] = f'attrib{cnt}'
        att_connect_eltype[cnt-1][1] = f'connect{cnt}'

        att_var = finfo.variables[att_connect_eltype[cnt-1][1]]
        typstr = att_var.elem_type

        if typstr == 'SHELL4':
            att_connect_eltype[cnt-1][2] = 'QUAD4'
        elif typstr == 'TETRA10':
            att_connect_eltype[cnt-1][2] = typstr
        elif typstr == 'SHELL9':
            att_connect_eltype[cnt-1][2] = 'QUAD9'
        elif typstr == 'TRI3':
            att_connect_eltype[cnt-1][2] = typstr
        else:
            print('ERROR: Element Type needs to be modified !!!')
            return
        elem_type_vec[cnt-1] = att_connect_eltype[cnt-1][2]

    ConnectAttrib = [[None, None] for _ in range(num_blk)]
    for cnt in range(1, num_blk + 1):
        cstr = att_connect_eltype[cnt-1][0]
        att_vec = finfo.variables[cstr][:]
        nelem = nelem_nnod_natt[cnt-1, 0]
        natt = nelem_nnod_natt[cnt-1, 2]
        ConnectAttrib[cnt-1][1] = att_vec.reshape((natt, nelem)).T

        cstr = att_connect_eltype[cnt-1][1]
        con_vec = finfo.variables[cstr][:]
        nnod = nelem_nnod_natt[cnt-1, 1]
        ConnectAttrib[cnt-1][0] = con_vec.reshape((nnod, nelem)).T

    NodXYZ = np.zeros((num_node, 3))
    for var in finfo.variables.values():
        tch = var.name
        if tch == 'coord':
            xyz_vec = var[:]
            NodXYZ = xyz_vec.reshape((num_node, 3))
            break
        elif tch == 'coordx':
            xyz_vec = var[:]
            NodXYZ[:, 0] = xyz_vec
        elif tch == 'coordy':
            xyz_vec = var[:]
            NodXYZ[:, 1] = xyz_vec
        elif tch == 'coordz':
            xyz_vec = var[:]
            NodXYZ[:, 2] = xyz_vec

    avgxyz = np.mean(np.abs(NodXYZ))
    nodf0 = np.abs(NodXYZ) < avgxyz * 1.0E-12
    NodXYZ[nodf0] = 0

    L = len(file_out)
    affix = '.elx'
    if file_out[(L-3):L] == affix:
        file_out = file_out[:L-4]

    with open(file_out + affix, 'w') as fwid:
        fwid.write(f'#Date: {str(datetime.now())}, by Breeze\n')
        fwid.write(f'{num_blk:4d}{num_elem:8d}{num_node:8d} \t#num_blk, num_elem, num_node\n')
        fwid.write(f'{nelem_nnod_natt[0, 1]:4d}{att_connect_eltype[0][2]:12s} \t#num_nod_per_elem, elem_type\n')

        wfmt = ' '.join(['%0d'] * nelem_nnod_natt[0, 1]) + '\n'

        for cnt in range(1, num_blk + 1):
            elem_typ_str = att_connect_eltype[0][2]
            L = len(elem_typ_str)
            while elem_typ_str[L-1].isdigit():
                elem_typ_str = elem_typ_str[:L-1]
                L = len(elem_typ_str)
            if elem_typ_str == 'SHELL':
                if len(ConnectAttrib[cnt-1][1][0]) == 2:
                    shell_attrib = ConnectAttrib[cnt-1][1][0, 0] * 1000 + ConnectAttrib[cnt-1][1][0, 1]
                else:
                    shell_attrib = ConnectAttrib[cnt-1][1][0]
                fwid.write(f'{nelem_nnod_natt[cnt-1, 0]:8d}{cnt:4d}{shell_attrib:8d} \t#num_elem_blk, blk_id, atrib_blk\n')
            else:
                fwid.write(f'{nelem_nnod_natt[cnt-1, 0]:8d}{cnt:4d}{ConnectAttrib[cnt-1][1][0]:8d} \t#num_elem_blk, blk_id, atrib_blk\n')

            for m in range(nelem_nnod_natt[cnt-1, 0]):
                fwid.write(wfmt % tuple(ConnectAttrib[cnt-1][0][m]))

        for m in range(num_node):
            v = NodXYZ[m]
            w = np.dot(MatRotY, v)
            fwid.write(f'%.16G %.16G %.16G \n' % tuple(w))

    if pf:
        print(f'Mesh converted successfully: {file_out + affix} !!!')


