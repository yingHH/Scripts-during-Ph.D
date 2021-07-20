# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 16:03:35 2018

@author: yinghuang
"""
from collections import namedtuple
import fire
import inspect
import pandas as pd


def main(ifile):
    print(get_log_info(ifile))


def get_log_info(ifile):
    logInfo = StatLog(ifile)
    res = pd.DataFrame(list(logInfo.Result))
    return res


class StatLog:

    Obj_RegisterCls = None
    Result = tuple()

    def __init__(self, ifile):
        for line in self._read_log(ifile):
            #print(line)
            self._choose_RegisterCls(line)
            if self.Obj_RegisterCls is None:
                continue
            res = self._handle_log(line)
            self._update_Result(res)

    def _read_log(self, ifile):
        with open(ifile, 'rt') as f:
            for line in f:
                yield line.strip('\n')

    def _choose_RegisterCls(self, line):
        try:
            # If this line is start marker of registered class,
            # refresh "self.Obj_RegisterCls" to new Register Class
            self._refresh_RegisterCls(RegisterClsIdx[line])
            print('> get marker line:', line)
        except:
            # If this line is end marker of registered class,
            # refresh "self.Obj_RegisterCls" to NoneType
            if self.Obj_RegisterCls is not None and self.Obj_RegisterCls.edLine == line:
                self._refresh_RegisterCls(None)
                print('> get marker line:', line)

    def _refresh_RegisterCls(self, regCls):
        if regCls is not None:
            self.Obj_RegisterCls = regCls
        else:
            self.Obj_RegisterCls = None

    def _handle_log(self, line):
        if self.Obj_RegisterCls is not None:
            obj = self.Obj_RegisterCls.cls()
            return obj(line)

    def _update_Result(self, ituples):
        if isinstance(ituples, tuple):
            self.Result += ituples


"""The follow Classes are executants to extract information from log file."""

class CutadaptLogStat:

    Result = None

    # These methods for PE type
    def _get_TotalReadsPairs(self, line):
        headLine = "Total read pairs processed"
        if line.startswith(headLine):
            key = "Total read pairs processed"
            txt = line.strip(' ')
            value = int(txt.split(' ')[-1].replace(',', ''))
            self.Result = ((key, value),)
            print('> get info:', ((key, value),))

    def _get_R1WithAdapter(self, line):
        headLine = "  Read 1 with adapter"
        if line.startswith(headLine):
            key = "Read 1 with adapter"
            txt = line.strip(' ')
            value = int(txt.split(' ')[-2].replace(',', ''))
            self.Result = ((key, value),)
            print('> get info:', ((key, value),))

    def _get_R2WithAdapter(self, line):
        headLine = "  Read 2 with adapter"
        if line.startswith(headLine):
            key = "Read 2 with adapter"
            txt = line.strip(' ')
            value = int(txt.split(' ')[-2].replace(',', ''))
            self.Result = ((key, value),)
            print('> get info:', ((key, value),))

    def _get_PairTooShort(self, line):
        headLine = "Pairs that were too short"
        if line.startswith(headLine):
            key = "Pairs that were too short"
            txt = line.strip(' ')
            value = int(txt.split(' ')[-2].replace(',', ''))
            self.Result = ((key, value),)
            print('> get info:', ((key, value),))

    def _get_PairWrite(self, line):
        headLine = "Pairs written (passing filters)"
        if line.startswith(headLine):
            key = "Pairs written (passing filters)"
            txt = line.strip(' ')
            value = int(txt.split(' ')[-2].replace(',', ''))
            self.Result = ((key, value),)
            print('> get info:', ((key, value),))

    # These methods for SE type
    def _get_TotalReads(self, line):
        headLine = "Total reads processed:"
        if line.startswith(headLine):
            key = "Total reads processed"
            txt = line.strip(' ')
            value = int(txt.split(' ')[-1].replace(',', ''))
            self.Result = ((key, value),)
            print('> get info:', ((key, value),))

    def _get_ReadsWithAdapter(self, line):
        headLine = "Reads with adapters:"
        if line.startswith(headLine):
            key = "Reads with adapters"
            txt = line.strip(' ')
            value = int(txt.split(' ')[-2].replace(',', ''))
            self.Result = ((key, value),)
            print('> get info:', ((key, value),))

    def _get_ReadsTooShort(self, line):
        headLine = "Reads that were too short:"
        if line.startswith(headLine):
            key = "Reads that were too short"
            txt = line.strip(' ')
            value = int(txt.split(' ')[-2].replace(',', ''))
            self.Result = ((key, value),)
            print('> get info:', ((key, value),))

    def _get_ReadsWrite(self, line):
        headLine = "Reads written (passing filters):"
        if line.startswith(headLine):
            key = "Reads written (passing filters)"
            txt = line.strip(' ')
            value = int(txt.split(' ')[-2].replace(',', ''))
            self.Result = ((key, value),)
            print('> get info:', ((key, value),))

    def _lst_handle_func(self):
        lst = [(name, value)
         for name, value in inspect.getmembers(CutadaptLogStat)
         if name.startswith('_get_') and inspect.isfunction(value)]
        #print(lst)
        return tuple(lst)

    def __call__(self, line):
        for name, value in self._lst_handle_func():
            getattr(self, name)(line)
            #print('*** method : ', name)
            if self.Result is not None:
                #print("... cutadapt result:", self.Result)
                return self.Result



class Hisat2LogStat:

    Result = None

    # These methods for PE type
    def _get_TotalPairs(self, line):
        tailLine = "were paired; of these:"
        if line.endswith(tailLine):
            key = "Total Pairs"
            txt = line.strip(' ')
            value = int(txt.split(' ')[0])
            self.Result = ((key, value),)
            print('> get info:', ((key, value),))

    def _get_0Concordantly(self, line):
        tailLine = "aligned concordantly 0 times"
        if line.endswith(tailLine):
            key = "Aligned Concordantly 0 Times"
            txt = line.strip(' ')
            value = int(txt.split(' ')[0])
            self.Result = ((key, value),)
            print('> get info:', ((key, value),))

    def _get_1Concordantly(self, line):
        tailLine = "aligned concordantly exactly 1 time"
        if line.endswith(tailLine):
            key = "Aligned Concordantly 1 Times"
            txt = line.strip(' ')
            value = int(txt.split(' ')[0])
            self.Result = ((key, value),)
            print('> get info:', ((key, value),))

    def _get_m1Concordantly(self, line):
        tailLine = "aligned concordantly >1 times"
        if line.endswith(tailLine):
            key = "Aligned Concordantly >1 Times"
            txt = line.strip(' ')
            value = int(txt.split(' ')[0])
            self.Result = ((key, value),)
            print('> get info:', ((key, value),))

    # These methods for SE type
    def _get_TotalReads(self, line):
        tailLine = "were unpaired; of these:"
        if line.endswith(tailLine):
            key = "Total Reads"
            txt = line.strip(' ')
            value = int(txt.split(' ')[0])
            self.Result = ((key, value),)
            print('> get info:', ((key, value),))

    def _get_Aligned0Times(self, line):
        tailLine = "aligned 0 times"
        if line.endswith(tailLine):
            key = "Aligned 0 Times"
            txt = line.strip(' ')
            value = int(txt.split(' ')[0])
            self.Result = ((key, value),)
            print('> get info:', ((key, value),))

    def _get_Aligned1Times(self, line):
        tailLine = "aligned exactly 1 time"
        if line.endswith(tailLine):
            key = "Aligned 1 Times"
            txt = line.strip(' ')
            value = int(txt.split(' ')[0])
            self.Result = ((key, value),)
            print('> get info:', ((key, value),))

    def _get_AlignedM1Times(self, line):
        tailLine = "aligned >1 times"
        if line.endswith(tailLine):
            key = "Aligned >1 Times"
            txt = line.strip(' ')
            value = int(txt.split(' ')[0])
            self.Result = ((key, value),)
            print('> get info:', ((key, value),))

    # The common method
    def _get_mappingRate(self, line):
        tailLine = "overall alignment rate"
        if line.endswith(tailLine):
            key = "Overall alignment rate"
            txt = line.strip(' ')
            value = txt.split(' ')[0]
            self.Result = ((key, value),)
            print('> get info:', ((key, value),))

    def _lst_handle_func(self):
        lst = [(name, value)
         for name, value in inspect.getmembers(Hisat2LogStat)
         if name.startswith('_get_') and inspect.isfunction(value)]
        #print(lst)
        return tuple(lst)

    def __call__(self, line):
        for name, value in self._lst_handle_func():
            getattr(self, name)(line)
            #print('*** method : ', name)
            if self.Result is not None:
                #print("... cutadapt result:", self.Result)
                return self.Result


class RMDupLogStat:

    Result = None

    def _get_ReadWriteReads(self, line):
        startLine = "READ"
        if line.startswith(startLine):
            txt = line.strip(' ').split(' ')
            rkey = "READ Reads"
            rvalue = int(txt[1])
            wkey = "WRITTEN Reads"
            wvalue = int(txt[3])
            self.Result = ((rkey, rvalue), (wkey, wvalue))
            print('> get info:', ((rkey, rvalue), (wkey, wvalue)))

    """
    def _get_ExaminedReads(self, line):
        headLine = "EXCLUDED"
        if line.startswith(headLine):
            key = "EXAMINED Reads"
            txt = line.strip(' ')
            value = int(txt.split(' ')[3])
            self.Result = ((key, value),)
            print('> get info:', ((key, value),))

    def _get_WriteReads(self, line):
        startLine = "READ"
        if line.startswith(startLine):
            txt = line.strip(' ').split(' ')
            key = "WRITTEN Reads"
            value = int(txt[3])
            self.Result = ((key, value),)
            print('> get info:', ((key, value),))

    def _get_DuplicateTotaldReads(self, line):
        headLine = "DUPLICATE TOTAL"
        if line.startswith(headLine):
            key = "DUPLICATE TOTAL Reads"
            txt = line.strip(' ')
            value = int(txt.split(' ')[2])
            self.Result = ((key, value),)
            print('> get info:', ((key, value),))
    """

    def _lst_handle_func(self):
        lst = [(name, value)
         for name, value in inspect.getmembers(RMDupLogStat)
         if name.startswith('_get_') and inspect.isfunction(value)]
        #print(lst)
        return tuple(lst)

    def __call__(self, line):
        for name, value in self._lst_handle_func():
            getattr(self, name)(line)
            #print('*** method : ', name)
            if self.Result is not None:
                #print("... cutadapt result:", self.Result)
                return self.Result


# Register new class when a new class is used to analysis log file
RegisterCls = namedtuple('RegisterCls', ['stLine','edLine','cls'])
    # "stLine" : the start part of log file which the new class will use
    # "edLine" : the end part of log file which the new class will use
    # "cls" : the reference to new class

# Use "RegisterCls" to register new class ("Cutadapt, Hisat2, RMDup, ...")
Cutadapt = RegisterCls("=== Running cutadapt ... ===", "=== cutadapt Done ===", CutadaptLogStat)
Hisat2 = RegisterCls("=== Running hisat2 ... ===", "=== hisat2 Done ===", Hisat2LogStat)
RMDup = RegisterCls("=== mark and remove duplication ===", "=== index bam ===", RMDupLogStat)
# After new class registered, registered class must be added to "RegisterClsIdx" dict for searching
RegisterClsIdx = {
    Cutadapt.stLine : Cutadapt,
    Hisat2.stLine : Hisat2,
    RMDup.stLine : RMDup,
}


if __name__ == '__main__':
    fire.Fire(main)
