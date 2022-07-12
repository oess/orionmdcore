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


def orion_session(f):
    def wrap(*args, **kwargs):

        if in_orion():
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
        else:
            session = None

        return f(*args, session=session, **kwargs)

    return wrap


# class ParmedData(CustomHandler):
#     @staticmethod
#     def get_name():
#         return "Parmed"
#
#     @classmethod
#     def validate(cls, value):
#         return isinstance(value, parmed.structure.Structure)
#
#     @classmethod
#     def copy(cls, value):
#         return parmed.structure.copy(value)
#
#     @staticmethod
#     def serialize(structure):
#         struct_dict = structure.__getstate__()
#         pkl_obj = pickle.dumps(struct_dict)
#         return bytes(pkl_obj)
#
#     @staticmethod
#     def deserialize(data):
#         new_structure = parmed.structure.Structure()
#         new_structure.__setstate__(pickle.loads(bytes(data)))
#         return new_structure


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
        pmd_struct_dict = structure.__getstate__()
        json_str = json.dumps(pmd_struct_dict)
        return json_str.encode(encoding='utf8')

    @staticmethod
    def deserialize(json_str_bytes):
        try:
            json_str = json_str_bytes.decode()
            pmd_dic = json.loads(json_str)
            new_structure = parmed.structure.Structure()
            new_structure.__setstate__(pmd_dic)
            return new_structure
        except:
            new_structure = parmed.structure.Structure()
            new_structure.__setstate__(pickle.loads(bytes(json_str_bytes)))
            return new_structure


# class MDStateData(CustomHandler):
#     @staticmethod
#     def get_name():
#         return "MDState"
#
#     @classmethod
#     def validate(cls, value):
#         return isinstance(value, MDState)
#
#     @classmethod
#     def copy(cls, value):
#         return copy.deepcopy(value)
#
#     @staticmethod
#     def serialize(state):
#         pkl_obj = pickle.dumps(state)
#         return bytes(pkl_obj)
#
#     @staticmethod
#     def deserialize(data):
#         new_state = pickle.loads(bytes(data))
#         return new_state

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
    def serialize(mdstate):
        mdstate_dic = mdstate.__getstate__()
        json_str = json.dumps(mdstate_dic)
        return json_str.encode(encoding='utf8')

    @staticmethod
    def deserialize(json_str_bytes):
        try:
            json_str = json_str_bytes.decode()
            mdstate_dic = json.loads(json_str)
            new_md_state = MDState()
            new_md_state.__setstate_new_api__(mdstate_dic)
            return new_md_state
        except:
            new_state = pickle.loads(bytes(json_str_bytes))
            return new_state

# class MDComponentData(CustomHandler):
#     @staticmethod
#     def get_name():
#         return "MDComponents"
#
#     @classmethod
#     def validate(cls, value):
#         return isinstance(value, MDComponents)
#
#     @classmethod
#     def copy(cls, components):
#         return copy.deepcopy(components)
#
#     @staticmethod
#     def serialize(components):
#         pkl_obj = pickle.dumps(components)
#         return bytes(pkl_obj)
#
#     @staticmethod
#     def deserialize(components):
#         new_components = pickle.loads(bytes(components))
#         return new_components


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
        return json_str.encode(encoding='utf8')

    @staticmethod
    def deserialize(json_str_bytes):
        try:
            json_str = json_str_bytes.decode()
            comp_dic = json.loads(json_str)
            new_md_components = MDComponents()
            new_md_components.__setstate_new_api__(comp_dic)
            return new_md_components
        except:
            print("Deprecation Warning: Old API in use. This will not be supported in further releases at any point")
            new_components = pickle.loads(bytes(json_str_bytes))
            return new_components


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


@orion_session
def upload_file(filename, session=None, orion_ui_name="OrionFile"):

    # Check if file is empty
    file_size = os.path.getsize(filename)
    if file_size == 0:
        raise ValueError("The file is empty")

    if in_orion():

        if session is None:
            raise ValueError("Session cannot be None")

        file_upload = File.upload(session, orion_ui_name, filename)

        session.tag_resource(file_upload, "Trajectory")

        job_id = environ.get("ORION_JOB_ID")

        if job_id:
            session.tag_resource(file_upload, "Job {}".format(job_id))

        file_id = file_upload.id

        # # Delete the local Orion file after uploading it
        # os.remove(filename)

    else:
        file_id = filename

    return file_id


@orion_session
def download_file(file_id, filename, session=None):

    if in_orion():

        if session is None:
            raise ValueError("Session cannot be None")

        resource = session.get_resource(File, file_id)

        resource.download_to_file(filename)

        fn_local = filename

    else:
        fn_local = file_id

    if not os.path.isfile(fn_local):
        raise IOError("The trajectory file has not been found: {}".format(fn_local))

    return fn_local


@orion_session
def delete_file(file_id, session=None):

    if in_orion():

        if session is None:
            raise ValueError("Session cannot be None")

        resource = session.get_resource(File, file_id)

        session.delete_resource(resource)

    else:
        os.remove(file_id)

    return True


@orion_session
def upload_data(filename, collection_id=None, shard_name="", session=None):

    if in_orion():

        if session is None:
            raise ValueError("Session cannot be None")

        if collection_id is None:
            raise ValueError("The Collection ID is None")

        collection = session.get_resource(ShardCollection, collection_id)

        if collection.state == "open":
            pass
        elif collection.state == "ready":
            collection.open()
        else:
            raise ValueError("Collection is not in an Open State: {}".format(collection.state))

        shard = Shard.create(collection, name=shard_name)
        shard.upload_file(filename)

        file_id = shard.id

        shard.close()

    else:
        file_id = filename

    return file_id


@orion_session
def download_data(file_id, path, session=None, collection_id=None):

    if in_orion():

        if session is None:
            raise ValueError("Session cannot be None")

        if collection_id is None:
            raise ValueError("The Collection ID is None")

        collection = session.get_resource(ShardCollection, collection_id)

        shard = session.get_resource(Shard(collection=collection), file_id)

        from orionmdcore.standards import MDFileNames

        fn_local = os.path.join(path, MDFileNames.mddata)

        shard.download_to_file(fn_local)

        shard.close()

    else:
        fn_local = file_id

    if not os.path.isfile(fn_local):
        raise IOError("The MD data file has not been found: {}".format(fn_local))

    return fn_local


@orion_session
def delete_data(file_id, session=None, collection_id=None):

    if in_orion():

        if session is None:
            raise ValueError("Session cannot be None")

        if collection_id is None:
            raise ValueError("The Collection ID is None")

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
