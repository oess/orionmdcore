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

import importlib.util
import pickle

OPENMM_SPEC = importlib.util.find_spec("openmm")


class MDOrionUnpickler(pickle.Unpickler):
    """Unpickler that supports deserialization of old-format MD records."""

    def find_class(self, module, name):
        if module.startswith("oemdtoolbox") and name == "MDComponents":
            module = "orionmdcore.forcefield.md_components"
        elif module.startswith("MDOrion") and name == "MDState":
            module = "orionmdcore.mdengine.utils"
        elif OPENMM_SPEC is not None and module.startswith("simtk"):
            words = module.split(".")
            if len(words) > 1 and words[1] == "openmm":
                words = words[1:]  # pop "simtk" prefix
            else:
                words[0] = "openmm"  # replace "simtk" with "openmm"
            module = ".".join(words)
        return super().find_class(module, name)
