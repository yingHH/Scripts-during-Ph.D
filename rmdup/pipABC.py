# -*- coding: utf-8 -*-
"""
Created on Sun Mar 18 16:32:04 2018

@author: yinghuang
"""
import os


PWD = os.path.split(os.path.realpath(__file__))[0]
shPath = os.path.join(PWD, 'piplib')


class RunPip:

    def __init__(self, shName, ifile, odir):
        """
        pipABC:function
        shName:str
        ifile:str
        odir:str
        """
        self.shName = shName
        self.ifile = ifile
        self.odir = odir

    def pipABC(self):
        script = os.path.join(shPath, self.shName)
        cmd = 'bash {} {} {}'
        signal = os.system(
            cmd.format(script, self.ifile, self.odir)
        )
        return signal

    def __repr__(self):
        self.signal = self.pipABC()
        return "Shell return signal is {}".format(str(self.signal))
