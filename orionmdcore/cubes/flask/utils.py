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
