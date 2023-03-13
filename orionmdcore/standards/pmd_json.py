# (C) 2022 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

import json

import parmed

from orionmdcore.mdengine import json_numpy

import numpy

import ast


def parmed_to_json(pmd_structure):

    class ResidueEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, parmed.Residue):
                return obj.__getstate__()

            return json.JSONEncoder.default(self, obj)

    class AtomEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, parmed.Atom):
                return obj.__getstate__()

            return json.JSONEncoder.default(self, obj)

    class AtomTypeEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, parmed.AtomType):
                dic = dict(name=obj.name, number=obj.number, mass=obj.mass,
                           atomic_number=obj.atomic_number, bond_type=obj._bond_type,
                           charge=obj.charge, epsilon=obj.epsilon, rmin=obj.rmin,
                           epsilon_14=obj.epsilon_14, rmin_14=obj.rmin_14, nbfix=obj.nbfix.copy())
                return dic
                # return str(obj.__str__())

            return json.JSONEncoder.default(self, obj)

    class BondTypeEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, parmed.BondType):
                return str(obj.__getstate__())

            return json.JSONEncoder.default(self, obj)

    class AngletypeEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, parmed.AngleType):
                return str(obj.__getstate__())

            return json.JSONEncoder.default(self, obj)

    class DihedraltypeEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, parmed.DihedralType):
                return str(obj.__getstate__())

            return json.JSONEncoder.default(self, obj)

    class ndarraytypeEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, numpy.ndarray):
                return json_numpy.dumps(obj)

            return json.JSONEncoder.default(self, obj)

    class FortranFormatEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, parmed.amber.FortranFormat):
                return obj.__getstate__()

            return json.JSONEncoder.default(self, obj)

    class CMapTypeEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, parmed.CmapType):
                return obj.__getstate__()

            return json.JSONEncoder.default(self, obj)

    class CmapGridEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, parmed.topologyobjects._CmapGrid):
                return {'resolution': obj.resolution, '_data': obj._data}

            return json.JSONEncoder.default(self, obj)

    class NonbondedExceptionTypeEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, parmed.NonbondedExceptionType):
                return obj.__getstate__()

            return json.JSONEncoder.default(self, obj)

    class MultipleJsonEncoders():

        def __init__(self, *encoders):
            self.encoders = encoders
            self.args = ()
            self.kwargs = {}

        def default(self, obj):
            for encoder in self.encoders:
                try:
                    return encoder(*self.args, **self.kwargs).default(obj)
                except TypeError:
                    pass
            raise TypeError(f'Object of type {obj.__class__.__name__} is not JSON serializable')

        def __call__(self, *args, **kwargs):
            self.args = args
            self.kwargs = kwargs
            enc = json.JSONEncoder(*args, **kwargs)
            enc.default = self.default
            return enc

    #################

    encoders = MultipleJsonEncoders(ResidueEncoder,
                                    AtomEncoder,
                                    AtomTypeEncoder,
                                    BondTypeEncoder,
                                    AngletypeEncoder,
                                    DihedraltypeEncoder,
                                    ndarraytypeEncoder,
                                    FortranFormatEncoder,
                                    CMapTypeEncoder,
                                    CmapGridEncoder,
                                    NonbondedExceptionTypeEncoder)

    pmd_dic_dic = pmd_structure.__getstate__()

    pmd_dic_dic['links'] = pmd_structure.links

    json_str = json.dumps(pmd_dic_dic, cls=encoders)

    return json_str


def json_to_parmed(json_str):

    pmd_dic = json.loads(json_str)

    for attr in pmd_dic:
        if attr == 'residues':
            residues_list = pmd_dic['residues']
            pmd_residues_list = parmed.ResidueList()

            if residues_list:
                for res_dic in residues_list:
                    atoms_list = res_dic['atoms']
                    new_pmd_atom_list = list()
                    for at_dic in atoms_list:
                        at_type = at_dic['atom_type']
                        pmd_atom_type = parmed.AtomType(at_type['name'],
                                                        at_type['number'],
                                                        at_type['mass'],
                                                        atomic_number=at_type['atomic_number'],
                                                        bond_type=at_type['bond_type'],
                                                        charge=at_type['charge']
                                                        )
                        pmd_atom_type.set_lj_params(at_type['epsilon'],
                                                    at_type['rmin'],
                                                    eps14=at_type['epsilon_14'],
                                                    rmin14=at_type['rmin_14'])
                        pmd_atom = parmed.Atom()
                        pmd_atom.__setstate__(at_dic)
                        pmd_atom.atom_type = pmd_atom_type
                        new_pmd_atom_list.append(pmd_atom)

                    res_dic['atoms'] = new_pmd_atom_list

                    pmd_residue = parmed.Residue("UNK")
                    pmd_residue.__setstate__(res_dic)
                    pmd_residues_list.append(pmd_residue)

            pmd_dic['residues'] = pmd_residues_list

        if attr == 'bond_types':

            bond_types_list = pmd_dic['bond_types']
            pmd_bond_types_list = parmed.TrackedList()

            if bond_types_list:

                for bn_type_str in bond_types_list:
                    bn_type_dic = ast.literal_eval(bn_type_str)

                    pmd_bn_type = parmed.BondType(k=bn_type_dic['k'], req=bn_type_dic['req'])
                    pmd_bn_type.__dict__.update(bn_type_dic)
                    pmd_bond_types_list.append(pmd_bn_type)

            pmd_dic['bond_types'] = pmd_bond_types_list

        if attr == 'angle_types':

            angle_types_list = pmd_dic['angle_types']
            pmd_angle_types_list = parmed.TrackedList()

            if angle_types_list:
                for an_type_str in angle_types_list:
                    an_type_dic = ast.literal_eval(an_type_str)
                    pmd_an_type = parmed.AngleType(k=an_type_dic['k'], theteq=an_type_dic['theteq'])
                    pmd_an_type.__dict__.update(an_type_dic)
                    pmd_angle_types_list.append(pmd_an_type)

            pmd_dic['angle_types'] = pmd_angle_types_list

        if attr == 'dihedral_types':

            dihedral_types_list = pmd_dic['dihedral_types']
            pmd_dihedral_types_list = parmed.TrackedList()

            if dihedral_types_list:
                for dh_type_str in dihedral_types_list:
                    dh_type_dic = ast.literal_eval(dh_type_str)
                    pmd_dh_type = parmed.DihedralType(phi_k=dh_type_dic['phi_k'], per=dh_type_dic['per'], phase=dh_type_dic['phase'])
                    pmd_dh_type.__dict__.update(dh_type_dic)
                    pmd_dihedral_types_list.append(pmd_dh_type)

            pmd_dic['dihedral_types'] = pmd_dihedral_types_list

        if attr == 'urey_bradley_types':

            urey_bradley_types_list = pmd_dic['urey_bradley_types']
            pmd_urey_bradley_types_list = parmed.TrackedList()

            # TODO MISSING IMPLEMENTATION
            if urey_bradley_types_list:
                for ub_type_str in urey_bradley_types_list:
                    raise ValueError("NOT Implemented")

            pmd_dic['urey_bradley_types'] = pmd_urey_bradley_types_list

        if attr == 'improper_types':

            improper_types_list = pmd_dic['improper_types']
            pmd_improper_types = parmed.TrackedList()

            # TODO MISSING IMPLEMENTATION
            if improper_types_list:
                for im_type_str in improper_types_list:
                    raise ValueError("NOT Implemented")

            pmd_dic['improper_types'] = pmd_improper_types

        if attr == 'rb_torsion_types':

            rb_torsion_types_list = pmd_dic['rb_torsion_types']
            pmd_rb_torsion_types = parmed.TrackedList()

            # TODO MISSING IMPLEMENTATION
            if rb_torsion_types_list:
                for rb_type_str in rb_torsion_types_list:
                    raise ValueError("NOT Implemented")

            pmd_dic['rb_torsion_types'] = pmd_rb_torsion_types

        if attr == 'cmap_types':

            cmap_types_list = pmd_dic['cmap_types']
            pmd_cmap_types_list = parmed.TrackedList()

            if cmap_types_list:
                for cmap_dic in cmap_types_list:
                    pmd_cmap = parmed.CmapType(cmap_dic['resolution'], cmap_dic['grid']['_data'])
                    pmd_cmap_types_list.append(pmd_cmap)

            pmd_dic['cmap_types'] = pmd_cmap_types_list

        if attr == 'trigonal_angle_types':

            trigonal_angle_types_list = pmd_dic['trigonal_angle_types']
            pmd_trigonal_angle_types = parmed.TrackedList()

            # TODO MISSING IMPLEMENTATION
            if trigonal_angle_types_list:
                for ta_type_str in trigonal_angle_types_list:
                    raise ValueError("NOT Implemented")

            pmd_dic['trigonal_angle_types'] = pmd_trigonal_angle_types

        if attr == 'adjust_types':

            adjust_types_list = pmd_dic['adjust_types']
            pmd_adjust_types_list = parmed.TrackedList()

            if adjust_types_list:
                for adj_type_dic in adjust_types_list:
                    pmd_adj_type = parmed.NonbondedExceptionType(rmin=adj_type_dic['rmin'],
                                                                 epsilon=adj_type_dic['epsilon'],
                                                                 chgscale=adj_type_dic['chgscale'])
                    pmd_adj_type.__dict__.update(adj_type_dic)
                    pmd_adjust_types_list.append(pmd_adj_type)

            pmd_dic['adjust_types'] = pmd_adjust_types_list

        if attr == 'out_of_plane_bend_types':

            out_of_plane_bend_types_list = pmd_dic['out_of_plane_bend_types']
            pmd_out_of_plane_bend_types = parmed.TrackedList()

            # TODO MISSING IMPLEMENTATION
            if out_of_plane_bend_types_list:
                for out_plane_type_str in out_of_plane_bend_types_list:
                    raise ValueError("NOT Implemented")

            pmd_dic['out_of_plane_bend_types'] = pmd_out_of_plane_bend_types

        if attr == 'pi_torsion_types':

            pi_torsion_types_list = pmd_dic['pi_torsion_types']
            pmd_pi_torsion_types = parmed.TrackedList()

            # TODO MISSING IMPLEMENTATION
            if pi_torsion_types_list:
                for pi_type_str in pi_torsion_types_list:
                    raise ValueError("NOT Implemented")

            pmd_dic['pi_torsion_types'] = pmd_pi_torsion_types

        if attr == 'groups':

            groups_list = pmd_dic['groups']
            pmd_groups = parmed.TrackedList()

            # TODO MISSING IMPLEMENTATION
            if groups_list:
                for groups_str in groups_list:
                    raise ValueError("NOT Implemented")

            pmd_dic['groups'] = pmd_groups

        if attr == 'stretch_bend_types':

            stretch_bend_types_list = pmd_dic['stretch_bend_types']
            pmd_stretch_bend_types = parmed.TrackedList()

            # TODO MISSING IMPLEMENTATION
            if stretch_bend_types_list:
                for stretch_str in stretch_bend_types_list:
                    raise ValueError("NOT Implemented")

            pmd_dic['stretch_bend_types'] = pmd_stretch_bend_types

        if attr == 'torsion_torsion_types':

            torsion_torsion_types_list = pmd_dic['torsion_torsion_types']
            pmd_torsion_torsion_types = parmed.TrackedList()

            # TODO MISSING IMPLEMENTATION
            if torsion_torsion_types_list:
                for tt_str in torsion_torsion_types_list:
                    raise ValueError("NOT Implemented")

            pmd_dic['torsion_torsion_types'] = pmd_torsion_torsion_types

        if attr == 'links':
            links_types_list = pmd_dic['links']
            pmd_links_types = parmed.TrackedList()

            # TODO MISSING IMPLEMENTATION
            if links_types_list:
                for tt_str in links_types_list:
                    raise ValueError("NOT Implemented")

            pmd_dic['links'] = pmd_links_types

        if attr == '_coordinates':
            if pmd_dic['_coordinates']:
                coords = json_numpy.loads(pmd_dic['_coordinates'])
                pmd_dic['_coordinates'] = coords
            else:
                pmd_dic['_coordinates'] = None

        if attr == '_box':
            if pmd_dic['_box']:
                box = json_numpy.loads(pmd_dic['_box'])
                pmd_dic['_box'] = box
            else:
                pmd_dic['_box'] = None

    new_structure = parmed.structure.Structure()

    new_structure.__setstate__(pmd_dic)

    return new_structure
