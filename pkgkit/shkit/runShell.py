# -*- coding: utf-8 -*-
"""
Created on 2018-07-10 17:37:49
Last Modified on 2018-07-10 17:37:49

Used as a module, to help run shell in python

@Author: Ying Huang
"""
import subprocess
from collections import namedtuple

out = namedtuple('out',['stdout', 'stderr'])

def runshell(cmd):
    stdout, stderr = subprocess.Popen(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE).communicate()

    sh_out = '[STDOUT]:\n' + stdout.decode('utf-8') + \
        '[STDERR]:\n' + stderr.decode('utf-8')

    return sh_out

def sh(cmd):
    stdout, stderr = subprocess.Popen(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE).communicate()

    out.stdout = stdout.decode('utf-8')
    out.stderr = stderr.decode('utf-8')

    return out
