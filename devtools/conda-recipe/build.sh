#!/bin/bash

username=$(python -c """import os
home = os.path.expanduser('~')

with open(os.path.join(home, '.pypirc'), 'r') as f:
    lines = f.readlines()

for ln in lines:
    if 'username' in ln:
        username = ln.split(':')[-1].strip()
print(username)
""")

password=$(python -c """import os
home = os.path.expanduser('~')

with open(os.path.join(home, '.pypirc'), 'r') as f:
    lines = f.readlines()

for ln in lines:
    if 'password' in ln:
        password = ln.split(':')[-1].strip()
print(password)
""")

export username
export password


python -m pip install https://$username:$password@magpie.eyesopen.com/download/503b/OpenEye-orionplatform-4.0.0.tar.gz
python -m pip install https://$username:$password@magpie.eyesopen.com/download/8053/OpenEye-Snowball-0.21.0.tar.gz
python -m pip install https://$username:$password@magpie.eyesopen.com/download/5626/OpenEye-brood-0.1.11.tar.gz
python -m pip install https://$username:$password@magpie.eyesopen.com/download/b929/OpenEye-brood-osx-x64-0.1.11.tar.gz
python -m pip install https://$username:$password@magpie.eyesopen.com/download/c682/OpenEye-floereport-0.1.11.tar.gz
python -m pip install https://$username:$password@magpie.eyesopen.com/download/6928/OpenEye-spruce-0.15.3.tar.gz
python -m pip install https://$username:$password@magpie.eyesopen.com/download/1991/OpenEye-mmds-client-1.1.5.tar.gz
python -m pip install https://$username:$password@magpie.eyesopen.com/download/1db0/OpenEye-szybki-0.1.2.tar.gz
python -m pip install https://$username:$password@magpie.eyesopen.com/download/7ef9/OpenEye-szybki-osx-x64-0.1.1.tar.gz
python -m pip install https://$username:$password@magpie.eyesopen.com/download/1e50/OpenEye-szmap-0.1.10.tar.gz
python -m pip install https://$username:$password@magpie.eyesopen.com/download/0be5/OpenEye-szmap-osx-x64-0.1.10.tar.gz
python -m pip install https://$username:$password@magpie.eyesopen.com/download/2675/OpenEye-sitehopper-0.15.3.tar.gz
python -m pip install https://$username:$password@magpie.eyesopen.com/download/a0d8/OpenEye-client-utils-0.5.5.tar.gz

python -m pip install https://$username:$password@magpie.eyesopen.com/download/0590/pygtop-2.1.4.tar.gz
python -m pip install https://$username:$password@magpie.eyesopen.com/download/40f8/docopt-0.6.2.tar.gz

python -m pip install https://$username:$password@magpie.eyesopen.com/download/80e4/molecupy-1.1.0-py3-none-any.whl


$PYTHON setup.py install
