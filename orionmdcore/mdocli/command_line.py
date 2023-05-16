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

try:

    from openeye import oechem

    from datarecord import read_records, OEField

    from datarecord import (Types,
                            OEWriteRecord,
                            OERecord,
                            Meta, OEFieldMeta)

    from orionmdcore.standards import Fields

    from orionmdcore.mdrecord import MDDataRecord

    from orionmdcore.standards import CollectionsNames

    import orionclient

    from orionclient.types import File, Shard, ShardCollection

    from orionclient.session import get_session
except ImportError:
    from orionmdcore import __installation__error__

    raise ImportError(__installation__error__)

import os

from tempfile import TemporaryDirectory

from tqdm import tqdm

import click

import gzip

import json

from orionmdcore.standards.utils import json_to_parmed


@click.group(
    context_settings={
        "help_option_names": ("-h", "--help")
    }
)
@click.pass_context
@click.option("--profile", help="OCLI profile name", default="default")
def main(ctx, profile):
    ctx.obj = dict()
    ctx.obj['session'] = orionclient.session.OrionSession(orionclient.session.get_profile_config(profile=profile),
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
                                                          ))


@main.group()
@click.argument('filename', type=click.Path(exists=True))
@click.option("--id", help="Record ID number", multiple=True, default=['all'])
@click.pass_context
def dataset(ctx, filename, id):
    """Records Extraction"""

    ctx.obj['filename'] = filename

    ifs = oechem.oeifstream(filename)

    records = []

    for rec in read_records(ifs):
        records.append(rec)
    ifs.close()

    if id[0] == 'all':
        ctx.obj['records'] = records
    else:

        list_rec = []
        for idx in id:

            if int(idx) < len(records):
                list_rec.append(records[int(idx)])
            else:
                raise ValueError("Wrong record number selection: {} > max = {}".format(id, len(records)))

        ctx.obj['records'] = list_rec


@dataset.command("makelocal")
@click.option("--name", help="Edit the trajectory file name", default="local.oedb")
@click.pass_context
def data_trajectory_extraction(ctx, name):

    def recursive(rec_in, rec_out):

        mdrecord_out = MDDataRecord(rec_out)
        mdrecord_in = MDDataRecord(rec_in)

        for fd in rec_in.get_fields(include_meta=True):

            # print(">>>>>>>>>>>", fd.get_name())

            if fd.get_name() == Fields.collections.get_name():
                continue

            elif fd.get_name() == Fields.md_stages.get_name():
                new_stages = []

                stages = mdrecord_in.get_stages

                system_title = mdrecord_in.get_title
                sys_id = mdrecord_in.get_flask_id

                for stage in stages:

                    stg_type = stage.get_value(Fields.stage_type)
                    new_stage = OERecord(stage)

                    with TemporaryDirectory() as output_directory:
                        data_fn = os.path.basename(output_directory) + '_' + system_title + '_' + str(
                            sys_id) + '-' + stg_type + '.tar.gz'
                        shard_id = stage.get_value(oechem.OEField("MDData_OPLMD", Types.Int))

                        protein_shard = session.get_resource(Shard(collection=collection), shard_id)
                        protein_shard.download_to_file(data_fn)

                        new_stage.delete_field(oechem.OEField("MDData_OPLMD", Types.Int))
                        new_stage.set_value(Fields.mddata, data_fn)

                        if stage.has_field(oechem.OEField("Trajectory_OPLMD", Types.Int)):
                            trj_field = stage.get_field("Trajectory_OPLMD")

                            trj_meta = trj_field.get_meta()
                            md_engine = trj_meta.get_attribute(Meta.Annotation.Description)

                            trj_id = stage.get_value(trj_field)
                            trj_fn = os.path.basename(output_directory) + '_' + system_title + '_' + str(
                                sys_id) + '-' + stg_type + '_traj' + '.tar.gz'

                            resource = session.get_resource(File, trj_id)
                            resource.download_to_file(trj_fn)

                            trj_meta = OEFieldMeta()
                            trj_meta.set_attribute(Meta.Annotation.Description, md_engine)
                            new_trj_field = oechem.OEField(Fields.trajectory.get_name(), Fields.trajectory.get_type(),
                                                           meta=trj_meta)

                            new_stage.delete_field(oechem.OEField("Trajectory_OPLMD", Types.Int))
                            new_stage.set_value(new_trj_field, trj_fn)

                    new_stages.append(new_stage)

                mdrecord_out.get_record.set_value(Fields.md_stages, new_stages)

            elif fd.get_name() == Fields.pmd_structure.get_name():
                # print("PARMED")

                pmd_id = rec_in.get_value(OEField('Structure_Parmed_OPLMD', Types.Int))
                protein_shard = session.get_resource(Shard(collection=collection), pmd_id)

                with TemporaryDirectory() as output_directory:

                    try:
                        parmed_fn = os.path.join(output_directory, "parmed.json")

                        protein_shard.download_to_file(parmed_fn)

                        with gzip.open(parmed_fn, 'rt', encoding='utf8') as f:
                            pmd_json = json.load(f)

                        pmd_structure = json_to_parmed(pmd_json)

                        mdrecord_out.set_parmed(pmd_structure)

                    except:
                        raise ValueError("Error Creating the Parmed Structure")

            elif fd.get_type() == oechem.Types.Record:

                rec_tmp_in = rec_in.get_value(fd)
                rec_tmp_out = oechem.OERecord()

                gac_tmp = recursive(rec_tmp_in, rec_tmp_out)

                rec_out.set_value(fd, gac_tmp)

            elif fd.get_type() == oechem.Types.RecordVec:
                rec_list = rec_in.get_value(fd)
                out_rec_list = list()
                for rec_tmp in rec_list:
                    new_rec = oechem.OERecord()
                    o_tmp_rec = recursive(rec_tmp, new_rec)
                    out_rec_list.append(o_tmp_rec)
                rec_out.set_value(fd, out_rec_list)

            else:
                rec_out.set_value(fd, rec_in.get_value(fd))

        return rec_out

    session = ctx.obj['session']

    ofs = oechem.oeofstream(name)

    for record in tqdm(ctx.obj['records']):

        new_record = OERecord(record)

        collections = record.get_value(Fields.collections)

        collection_id = collections[CollectionsNames.md]

        collection = session.get_resource(ShardCollection, collection_id)

        exit_rec = recursive(new_record, OERecord())

        OEWriteRecord(ofs, exit_rec, fmt='binary')

    ofs.close()


@dataset.command("logs")
@click.option("--stgn", help="MD Stage name", default="last")
@click.pass_context
def logs_extraction(ctx, stgn):

    for record in ctx.obj['records']:

        mdrecord = MDDataRecord(record)

        info = mdrecord.get_stage_info(stg_name=stgn)

        print(info)


@dataset.command("protein")
@click.pass_context
def protein_extraction(ctx):

    for record in ctx.obj['records']:

        mdrecord = MDDataRecord(record)

        if not record.has_value(Fields.protein):
            print("No protein have been found in the selected record")
            return
        else:
            title = mdrecord.get_title
            fn = title.split('_')[0] + ".oeb"
            with oechem.oemolostream(fn) as ofs:
                oechem.OEWriteConstMolecule(ofs, record.get_value(Fields.protein))
        print("Protein file generated: {}".format(fn))


@dataset.command("ligand")
@click.pass_context
def ligand_extraction(ctx):

    for record in ctx.obj['records']:

        mdrecord = MDDataRecord(record)

        if not record.has_value(Fields.ligand):
            print("No ligand have been found in the selected record")
            return
        else:
            title = mdrecord.get_title.split("_")[1:]
            title = "_".join(title)
            id = mdrecord.get_flask_id
            fn = title + "_" + str(id)+".oeb"
            with oechem.oemolostream(fn) as ofs:
                oechem.OEWriteConstMolecule(ofs, record.get_value(Fields.ligand))
        print("Ligand file generated: {}".format(fn))


@dataset.command("info")
@click.pass_context
def info_extraction(ctx):

    rec_size = 0

    def GetHumanReadable(size, precision=2):

        suffixes = ['B', 'KB', 'MB', 'GB', 'TB']
        suffixIndex = 0

        while size > 1024 and suffixIndex < 4:
            suffixIndex += 1  # increment the index of the suffix
            size = size / 1024.0  # apply the division

        color = ""
        if suffixIndex > 1:
            color = "\033[33m"

        return "%s %.*f %s \033[00m" % (color, precision, size, suffixes[suffixIndex])

    def recursive_record(record, level=0):

        nonlocal rec_size

        for field in record.get_fields():

            field_type = field.get_type()

            blank = "       "
            print("{} |".format(blank * (level + 1)))
            print("{} |".format(blank * (level + 1)))
            dis = "______"

            if not field_type == Types.Record and not field_type == Types.RecordVec:

                rec_size += len(record.get_bytes(field))

                if (field.get_type() == Types.String or
                        field.get_type() == Types.Int or
                        field.get_type() == Types.Float):
                    print("{} {} name = {}\n        "
                          "{}type = {}\n        "
                          "{}value = {}\n        "
                          "{}size = {}".format(blank * (level + 1),
                                               dis,
                                               field.get_name(),
                                               blank * (level + 1),
                                               field.get_type_name(),
                                               blank * (level + 1),
                                               str(record.get_value(field))[0:30],
                                               blank * (level + 1),
                                               GetHumanReadable(len(record.get_bytes(field)))
                                               ))
                else:
                    try:
                        size = GetHumanReadable(len(record.get_bytes(field)))
                    except:
                        size = "None"

                    print("{} {} name = {}\n        "
                          "{}type = {}\n        "
                          "{}size = {}".format(blank * (level + 1),
                                               dis,
                                               field.get_name(),
                                               blank * (level + 1),
                                               field.get_type_name(),
                                               blank * (level + 1),
                                               size
                                               ))

            elif field_type == Types.Record:
                print("{} {} RECORD: {}".format(blank * (level + 1), dis, field.get_name()))
                recursive_record(record.get_value(field), level + 1)

            elif field_type == Types.RecordVec:
                vec = record.get_value(field)
                print("{} {} RECORD VECTOR: {} containing {} records".format(blank * (level + 1),
                                                                             dis,
                                                                             field.get_name(),
                                                                             len(vec)))
                print("{} |".format(blank * (level + 2)))
                print("{} |".format(blank * (level + 2)))

                for idx in range(0, len(vec)):
                    print("{} {} RECORD # {}".format(blank * (level + 2),
                                                     dis,
                                                     idx))

                    recursive_record(vec[idx], level + 2)
                    if idx != len(vec) - 1:
                        print("{} |".format(blank * (level + 2)))
                        print("{} |".format(blank * (level + 2)))
            else:
                raise ValueError("Field type error: {}".format(field_type))

    for idx in range(0, len(ctx.obj['records'])):
        print(30 * "*" + " RECORD {}/{} ".format(idx + 1, len(ctx.obj['records'])) + 30 * "*")
        rec_size = 0
        recursive_record(ctx.obj['records'][idx], 0)
        print("\n" + 22 * "*" + " END RECORD - SIZE {} ".format(GetHumanReadable(rec_size)) + 23 * "*" + "\n")


@main.group()
@click.argument('filename', type=click.Path(exists=True))
@click.option("--id", help="Record ID number", default="all")
@click.pass_context
def analysis(ctx, filename, id):
    """Records Extraction"""

    ctx.obj['filename'] = filename

    ifs = oechem.oeifstream(filename)

    records = []

    for rec in read_records(ifs):
        records.append(rec)
    ifs.close()

    if id == 'all':
        ctx.obj['records'] = records
    else:
        if int(id) < len(records):
            ctx.obj['records'] = [records[int(id)]]
        else:
            raise ValueError("Wrong record number selection: {} > max = {}".format(int(id), len(records)))


def check_sys_id(record):
    if not record.has_value(Fields.title):
        raise ValueError("Flask title field is not present on the record")

    title = record.get_value(Fields.title)

    if not record.has_value(Fields.flaskid):
        raise ValueError("Flask ID not present on the record")

    flaskid = record.get_value(Fields.flaskid)

    sys_id = title + "_" + str(flaskid)

    return sys_id


@analysis.command("energy")
@click.pass_context
def energy_extraction(ctx):

    for record in ctx.obj['records']:

        sys_id = check_sys_id(record)

        if not record.has_field(Fields.Analysis.oeintE_rec):
            raise ValueError("Interaction Energy Record field is missing")

        oeintE_rec = record.get_value(Fields.Analysis.oeintE_rec)

        for fd in oeintE_rec.get_fields():
            name = fd.get_name()
            vec_list = oeintE_rec.get_value(fd)
            f = open(sys_id+"_"+name + '.txt', 'w')
            for val in vec_list:
                f.write(str(val) + "\n")
            f.close()


@analysis.command("mmpbsa")
@click.pass_context
def mmpbsa_extraction(ctx):

    for record in ctx.obj['records']:

        sys_id = check_sys_id(record)

        if not record.has_field(Fields.Analysis.oepbsa_rec):
            raise ValueError("PBSA record field is missing")

        oepbsa_rec = record.get_value(Fields.Analysis.oepbsa_rec)

        for fd in oepbsa_rec.get_fields():
            name = fd.get_name()
            vec_list = oepbsa_rec.get_value(fd)
            f = open(sys_id+"_"+name + '.txt', 'w')
            for val in vec_list:
                f.write(str(val) + "\n")
            f.close()


@analysis.command("clusters")
@click.pass_context
def cluster_extraction(ctx):

    for record in ctx.obj['records']:

        sys_id = check_sys_id(record)

        if not record.has_field(Fields.Analysis.oeclus_rec):
            raise ValueError("Cluster record field is missing")

        oeclus_rec = record.get_value(Fields.Analysis.oeclus_rec)

        clust_names = ["ClusLigAvgMol",
                       "ClusProtAvgMol",
                       "ClusLigMedMol",
                       "ClusProtMedMol"]

        ofs = oechem.oemolostream(sys_id+"_clusters.oeb")

        for fd in oeclus_rec.get_fields():

            name = fd.get_name()

            if name in clust_names:
                mol = oeclus_rec.get_value(fd)
                mol.SetTitle(name)
                oechem.OEWriteConstMolecule(ofs, mol)

        ofs.close()


@analysis.command("traj_confs")
@click.pass_context
def traj_conf_extraction(ctx):

    for record in ctx.obj['records']:

        sys_id = check_sys_id(record)

        if not record.has_field(Fields.Analysis.oetraj_rec):
            raise ValueError("Multi Conf Trajectory record field is missing")

        oetraj_rec = record.get_value(Fields.Analysis.oetraj_rec)

        traj_names = ["LigTraj",
                      "ProtTraj_OPLMD",
                      "WatTraj"]

        ofs = oechem.oemolostream(sys_id+"_traj_confs.oeb")

        for fd in oetraj_rec.get_fields():

            name = fd.get_name()

            if name in traj_names:
                mol = oetraj_rec.get_value(fd)
                mol.SetTitle(name)
                oechem.OEWriteConstMolecule(ofs, mol)

        ofs.close()
