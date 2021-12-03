# (C) 2020 OpenEye Scientific Software Inc. All rights reserved.
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

try:
    from openeye import oechem

    from datarecord import CustomHandler

    from orionclient.session import in_orion, OrionSession, get_session

    from orionclient.types import File

    from orionclient.types import Shard, ShardCollection

    from orionclient.helpers.collections import (
        try_hard_to_create_shard,
        try_hard_to_download_shard,
    )
except ImportError:
    from orionmdcore import __installation__error__
    raise ImportError(__installation__error__)

import pickle

import parmed

from orionmdcore.mdengine.utils import MDState

from orionmdcore.forcefield import MDComponents

import copy

from os import environ

import os

import json


class ParmedData(CustomHandler):
    @staticmethod
    def get_name():
        return "Parmed"

    @classmethod
    def validate(cls, value):
        return isinstance(value, parmed.structure.Structure)

    @classmethod
    def copy(cls, value):
        return parmed.structure.copy(value)

    @staticmethod
    def serialize(structure):
        struct_dict = structure.__getstate__()
        pkl_obj = pickle.dumps(struct_dict)
        return bytes(pkl_obj)

    @staticmethod
    def deserialize(data):
        new_structure = parmed.structure.Structure()
        new_structure.__setstate__(pickle.loads(bytes(data)))
        return new_structure


class MDStateData(CustomHandler):
    @staticmethod
    def get_name():
        return "MDState"

    @classmethod
    def validate(cls, value):
        return isinstance(value, MDState)

    @classmethod
    def copy(cls, value):
        return copy.deepcopy(value)

    @staticmethod
    def serialize(state):
        pkl_obj = pickle.dumps(state)
        return bytes(pkl_obj)

    @staticmethod
    def deserialize(data):
        new_state = pickle.loads(bytes(data))
        return new_state


class MDComponentData(CustomHandler):
    @staticmethod
    def get_name():
        return "MDComponents"

    @classmethod
    def validate(cls, value):
        return isinstance(value, MDComponents)

    @classmethod
    def copy(cls, components):
        return copy.deepcopy(components)

    @staticmethod
    def serialize(components):
        comp_dic = components.__getstate__()
        json_str = json.dumps(comp_dic)
        json_str.encode(encoding='utf8')
        return json_str.encode(encoding='utf8')

    @staticmethod
    def deserialize(json_str_bytes):
        json_str = json_str_bytes.decode()
        comp_dic = json.loads(json_str)
        new_md_components = MDComponents()
        new_md_components.__setstate__(comp_dic)
        return new_md_components


class DesignUnit(CustomHandler):
    @staticmethod
    def get_name():
        return "DesignUnit"

    @classmethod
    def validate(cls, value):
        return isinstance(value, oechem.OEDesignUnit)

    @classmethod
    def copy(cls, value):
        return copy.deepcopy(value)

    @staticmethod
    def serialize(du):
        return oechem.OEWriteDesignUnitToBytes(du)

    @staticmethod
    def deserialize(du_bytes):

        design_unit = oechem.OEDesignUnit()

        if not oechem.OEReadDesignUnitFromBytes(design_unit, du_bytes):
            raise ValueError("It was not possible to deserialize the Design Unit")

        return design_unit


def check_filename(filename):
    special_shell_characters = [
        "~",
        "`",
        "#",
        "$",
        "&",
        "*",
        "(",
        ")",
        "\\",
        "|",
        "[",
        "]",
        "{",
        "}",
        ";",
        "'",
        '"',
        ">",
        "<",
        "/",
        "?",
        "!",
    ]

    if any(c in special_shell_characters for c in filename):
        new_filename = ""
        for c in filename:
            if c in special_shell_characters:
                new_filename += "_"
            else:
                new_filename += c

        return new_filename
    else:
        return filename


def upload_file(filename, orion_ui_name="OrionFile"):

    if in_orion():

        # session = APISession

        session = OrionSession(
            requests_session=get_session(
                retry_dict={
                    403: 5,
                    404: 20,
                    409: 45,
                    460: 15,
                    500: 2,
                    502: 45,
                    503: 45,
                    504: 45,
                }
            )
        )

        file_upload = File.upload(session, orion_ui_name, filename)

        session.tag_resource(file_upload, "Trajectory")

        job_id = environ.get("ORION_JOB_ID")

        if job_id:
            session.tag_resource(file_upload, "Job {}".format(job_id))

        file_id = file_upload.id

    else:
        file_id = filename

    return file_id


def download_file(file_id, filename):

    if in_orion():

        # session = APISession

        session = OrionSession(
            requests_session=get_session(
                retry_dict={
                    403: 5,
                    404: 20,
                    409: 45,
                    460: 15,
                    500: 2,
                    502: 45,
                    503: 45,
                    504: 45,
                }
            )
        )

        resource = session.get_resource(File, file_id)

        resource.download_to_file(filename)

        fn_local = filename

    else:
        fn_local = file_id

    if not os.path.isfile(fn_local):
        raise IOError("The trajectory file has not been found: {}".format(fn_local))

    return fn_local


def delete_file(file_id):

    if in_orion():

        # session = APISession

        session = OrionSession(
            requests_session=get_session(
                retry_dict={
                    403: 5,
                    404: 20,
                    409: 45,
                    460: 15,
                    500: 2,
                    502: 45,
                    503: 45,
                    504: 45,
                }
            )
        )

        resource = session.get_resource(File, file_id)

        session.delete_resource(resource)

    else:
        os.remove(file_id)

    return True


def upload_data(filename, collection_id=None, shard_name=""):

    if in_orion():

        if collection_id is None:
            raise ValueError("The Collection ID is None")

        # session = APISession

        session = OrionSession(
            requests_session=get_session(
                retry_dict={
                    403: 5,
                    404: 20,
                    409: 45,
                    460: 15,
                    500: 2,
                    502: 45,
                    503: 45,
                    504: 45,
                }
            )
        )

        collection = session.get_resource(ShardCollection, collection_id)

        if collection.state == "open":
            pass
        elif collection.state == "ready":
            collection.open()
        else:
            raise ValueError("Collection is not in an Open State: {}".format(collection.state))

        shard = try_hard_to_create_shard(collection, filename, name=shard_name)

        file_id = shard.id

        shard.close()

    else:
        file_id = filename

    return file_id


def download_data(file_id, path, collection_id=None):

    if in_orion():

        if collection_id is None:
            raise ValueError("The Collection ID is None")

        # session = APISession

        session = OrionSession(
            requests_session=get_session(
                retry_dict={
                    403: 5,
                    404: 20,
                    409: 45,
                    460: 15,
                    500: 2,
                    502: 45,
                    503: 45,
                    504: 45,
                }
            )
        )

        collection = session.get_resource(ShardCollection, collection_id)

        shard = session.get_resource(Shard(collection=collection), file_id)

        from orionmdcore.standards import MDFileNames

        fn_local = os.path.join(path, MDFileNames.mddata)

        try_hard_to_download_shard(shard, fn_local)

        shard.close()

    else:
        fn_local = file_id

    if not os.path.isfile(fn_local):
        raise IOError("The MD data file has not been found: {}".format(fn_local))

    return fn_local


def delete_data(file_id, collection_id=None):

    if in_orion():

        if collection_id is None:
            raise ValueError("The Collection ID is None")

        # session = APISession

        session = OrionSession(
            requests_session=get_session(
                retry_dict={
                    403: 5,
                    404: 20,
                    409: 45,
                    460: 15,
                    500: 2,
                    502: 45,
                    503: 45,
                    504: 45,
                }
            )
        )

        collection = session.get_resource(ShardCollection, collection_id)

        if collection.state == "open":
            pass
        elif collection.state == "ready":
            collection.open()
        else:
            raise ValueError("Collection is not in an Open State: {}".format(collection.state))

        session.delete_resource(Shard(collection=collection, id=file_id))

    else:
        os.remove(file_id)

    return True
