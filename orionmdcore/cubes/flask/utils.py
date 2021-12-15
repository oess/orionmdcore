# (C) 2021 OpenEye Scientific Software Inc. All rights reserved.
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

from datarecord import OEField, Types
from orionmdcore.mdrecord import MDDataRecord
from orionmdcore.mdengine.utils import MDState
from orionmdcore.forcefield import MDComponents
from openeye import oechem


def get_human_readable(size, precision=2):

    suffixes = ["B", "KB", "MB", "GB", "TB"]
    suffixIndex = 0

    while size > 1024 and suffixIndex < 4:
        suffixIndex += 1  # increment the index of the suffix
        size = size / 1024.0  # apply the division

    return "%.*f %s" % (precision, size, suffixes[suffixIndex])


def parser_log_old_api(logs):
    lines = logs.split("\n")

    filtered_lines = list()
    for ln in lines:
        if '=' in ln:
            filtered_lines.append(ln)
        else:
            continue

    info_dic = dict()
    for ln in filtered_lines:
        ln = ln.split()
        tmp_key = ""
        for idx, word in enumerate(ln):
            if word == '=':
                try:
                    if '.' in ln[idx + 1]:
                        info_dic[tmp_key] = float(ln[idx + 1])
                    else:
                        info_dic[tmp_key] = int(ln[idx + 1])
                except Exception as e:
                    if ln[idx + 1] == 'None':
                        info_dic[tmp_key] = None
                    else:
                        info_dic[tmp_key] = ln[idx + 1]
            else:
                if idx == 0:
                    tmp_key = word
                else:
                    tmp_key += ' ' + word

    return info_dic


def convert(record):
    oetraj_rec_field = OEField("OETraj", Types.Record)
    lig_traj_field = OEField("LigTraj", Types.Chem.Mol)
    water_traj_field = OEField("WatTraj", Types.Chem.Mol)

    md_record = MDDataRecord(record)

    if md_record.has_stages:

        md_comp = md_record.get_md_components

        if str(type(md_comp)) == '<class \'oemdtoolbox.ForceField.md_components.MDComponents\'>':

            new_md_comp = MDComponents()

            for comp_name, comp in md_comp.get_components.items():
                if comp_name == 'ligand':
                    oechem.OEAssignAromaticFlags(comp, oechem.OEAroModelOpenEye)
                new_md_comp.set_component_by_name(comp_name, comp)

            new_md_comp.set_title(md_comp.get_title)

            if md_comp.has_box_vectors:
                new_md_comp.set_box_vectors(md_comp.get_box_vectors)

            md_record.set_md_components(new_md_comp)

            if record.has_field(oetraj_rec_field):

                oetraj_rec = record.get_value(oetraj_rec_field)
                md_oetraj_rec = MDDataRecord(oetraj_rec)

                if oetraj_rec.has_field(lig_traj_field):
                    lig_traj = oetraj_rec.get_value(lig_traj_field)
                    lig_traj.SetTitle("LigandTraj")
                    md_oetraj_rec.set_ligand_traj(lig_traj, "LigTrajConf_")

                if oetraj_rec.has_field(water_traj_field):
                    water_traj = oetraj_rec.get_value(water_traj_field)
                    water_traj.SetTitle("WaterTraj")
                    md_oetraj_rec.set_water_traj(water_traj, "WatTrajConf_")

                record.set_value(oetraj_rec_field, oetraj_rec)

            stage_names = md_record.get_stages_names

            for stgn in stage_names:

                if md_record.has_stage_info(stg_name=stgn):
                    pass
                else:
                    logs = md_record.get_stage_logs(stg_name=stgn)

                    if logs is not None:
                        info_dic = parser_log_old_api(logs)
                        md_record.set_stage_info(info_dic, stg_name=stgn)

                pmd = md_record.get_parmed(sync_stage_name=stgn)
                new_state = MDState(pmd)
                md_record.set_stage_state(new_state, stg_name=stgn)

        else:
            print("No Conversion is Needed. Bypass record")

    else:  # Analysis Output

        if md_record.has_md_components:

            md_comp = md_record.get_md_components

            if str(type(md_comp)) == '<class \'oemdtoolbox.ForceField.md_components.MDComponents\'>':

                new_md_comp = MDComponents()

                for comp_name, comp in md_comp.get_components.items():
                    if comp_name == 'ligand':
                        oechem.OEAssignAromaticFlags(comp, oechem.OEAroModelOpenEye)
                    new_md_comp.set_component_by_name(comp_name, comp)

                new_md_comp.set_title(md_comp.get_title)

                if md_comp.has_box_vectors:
                    new_md_comp.set_box_vectors(md_comp.get_box_vectors)

                md_record.set_md_components(new_md_comp)
            else:
                print("No Conversion is Needed. Bypass record")

        if md_record.has_field(oetraj_rec_field):
            convert_oetraj = False
            oetraj_rec = md_record.get_value(oetraj_rec_field)
            md_oetraj_rec = MDDataRecord(oetraj_rec)

            if oetraj_rec.has_field(lig_traj_field):
                lig_traj = oetraj_rec.get_value(lig_traj_field)
                lig_traj.SetTitle("LigandTraj")
                md_oetraj_rec.set_ligand_traj(lig_traj, "LigTrajConf_")
                oetraj_rec.delete_field(lig_traj_field)
                convert_oetraj = True

            if oetraj_rec.has_field(water_traj_field):
                water_traj = oetraj_rec.get_value(water_traj_field)
                water_traj.SetTitle("WaterTraj")
                md_oetraj_rec.set_water_traj(water_traj, "WatTrajConf_")
                oetraj_rec.delete_field(water_traj_field)
                convert_oetraj = True

            if convert_oetraj:
                md_record.set_value(oetraj_rec_field, oetraj_rec)

        #######

        oetrajconf_rec_fd = OEField("Lig_Conf_Data", Types.RecordVec)

        if record.has_value(oetrajconf_rec_fd):

            oetrajconf_records = record.get_value(oetrajconf_rec_fd)
            new_conf_records = list()

            for conf_rec in oetrajconf_records:

                md_conf_rec = MDDataRecord(conf_rec)

                if md_conf_rec.has_stages:

                    md_comp = md_conf_rec.get_md_components

                    if str(type(md_comp)) == '<class \'oemdtoolbox.ForceField.md_components.MDComponents\'>':

                        new_md_comp = MDComponents()

                        for comp_name, comp in md_comp.get_components.items():
                            if comp_name == 'ligand':
                                oechem.OEAssignAromaticFlags(comp, oechem.OEAroModelOpenEye)
                            new_md_comp.set_component_by_name(comp_name, comp)

                        new_md_comp.set_title(md_comp.get_title)

                        if md_comp.has_box_vectors:
                            new_md_comp.set_box_vectors(md_comp.get_box_vectors)

                        md_conf_rec.set_md_components(new_md_comp)

                        if md_conf_rec.has_field(oetraj_rec_field):

                            oetraj_rec = md_conf_rec.get_value(oetraj_rec_field)
                            md_oetraj_rec = MDDataRecord(oetraj_rec)

                            if oetraj_rec.has_field(lig_traj_field):
                                lig_traj = oetraj_rec.get_value(lig_traj_field)
                                lig_traj.SetTitle("LigandTraj")
                                md_oetraj_rec.set_ligand_traj(lig_traj, "LigTrajConf_")
                                oetraj_rec.delete_field(lig_traj_field)

                            if oetraj_rec.has_field(water_traj_field):
                                water_traj = oetraj_rec.get_value(water_traj_field)
                                water_traj.SetTitle("WaterTraj")
                                md_oetraj_rec.set_water_traj(water_traj, "WatTrajConf_")
                                oetraj_rec.delete_field(water_traj_field)

                            md_conf_rec.set_value(oetraj_rec_field, oetraj_rec)

                        stage_names = md_conf_rec.get_stages_names

                        for stgn in stage_names:

                            if md_conf_rec.has_stage_info(stg_name=stgn):
                                pass
                            else:
                                logs = md_conf_rec.get_stage_logs(stg_name=stgn)

                                if logs is not None:
                                    info_dic = parser_log_old_api(logs)
                                    md_conf_rec.set_stage_info(info_dic, stg_name=stgn)

                            pmd = md_conf_rec.get_parmed(sync_stage_name=stgn)
                            new_state = MDState(pmd)
                            md_conf_rec.set_stage_state(new_state, stg_name=stgn)

                    else:
                        print("No Conversion is Needed. Bypass record")

                new_conf_records.append(md_conf_rec)

            if new_conf_records:
                record.set_value(oetrajconf_rec_fd, new_conf_records)

    return record
