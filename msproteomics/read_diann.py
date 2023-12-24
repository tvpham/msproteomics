# coding=utf-8
# Copyright 2023 Thang V Pham
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import sys
import sysconfig
import subprocess
import site

def main():

    cmd = os.path.join(sysconfig.get_path('platlib'), 'msproteomics', 'siteloc')

    print('Activate siteloc (use -h for help) ...')

    subprocess.run([cmd] + sys.argv[1:])