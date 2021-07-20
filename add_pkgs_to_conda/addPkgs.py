"""
Add module to anaconda3
Add new module path to '/home/yingh/anaconda3/lib/python3.6/site-packages/myModule.pth'

@author: yinghuang
"""
import os


default_conf_file_path = '/home/yingh/anaconda3/lib/python3.6/site-packages/myModule.pth'

def addpkgs(*ipaths, conf_file_path=default_conf_file_path):
    for p in ipaths:
        assert isinstance(p, str) is True, \
        '({}) is not str, please input a "str" format value.'
        assert os.path.exists(p) is True, \
        '({}) is not exists, please input an available path.'

    with open(conf_file_path, 'rt') as f:
        pkgs_lst = set(f.read().split('\n'))
        pkgs_lst.update(set(ipaths))

    with open(conf_file_path, 'wt') as f:
        f.write('\n'.join(list(pkgs_lst)))
