# (C) 2022 Cadence Design Systems, Inc. (Cadence) 
# All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of Cadence products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. Cadence claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable Cadence offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall Cadence be
# liable for any damages or liability in connection with the Sample Code
# or its use.
import json
from base64 import b64decode, b64encode

import numpy as np
from numpy.lib.format import dtype_to_descr, descr_to_dtype


def default(obj):
    if isinstance(obj, (np.ndarray, np.generic)):
        return {
            '__numpy__': b64encode(obj.data if obj.flags.c_contiguous else obj.tobytes()).decode('ascii'),
            'dtype': dtype_to_descr(obj.dtype),
            'shape': obj.shape
        }
    raise TypeError(f'Object of type {type(obj)} is not JSON serializable')


def object_hook(dct):
    if '__numpy__' in dct:
        np_obj = np.frombuffer(b64decode(dct['__numpy__']), descr_to_dtype(dct['dtype']))
        shape = dct['shape']
        return np_obj.reshape(shape) if shape else np_obj[0]  # Scalar test
    return dct


_dumps = json.dumps
_loads = json.loads
_dump = json.dump
_load = json.load


def dumps(*args, **kwargs):
    kwargs.setdefault('default', default)
    return _dumps(*args, **kwargs)


def loads(*args, **kwargs):
    kwargs.setdefault('object_hook', object_hook)
    return _loads(*args, **kwargs)


def dump(*args, **kwargs):
    kwargs.setdefault('default', default)
    return _dump(*args, **kwargs)


def load(*args, **kwargs):
    kwargs.setdefault('object_hook', object_hook)
    return _load(*args, **kwargs)


def patch():
    """Monkey patches the json module in order to support serialization/deserialization of Numpy arrays and scalars."""
    json.dumps = dumps
    json.loads = loads
    json.dump = dump
    json.load = load
