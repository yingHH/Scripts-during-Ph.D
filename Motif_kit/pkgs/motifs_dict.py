# -*- coding: utf-8 -*-
"""
Created on 2018-10-19 21:45:44
Last Modified on 2018-10-19 21:45:44

Classes for reading in motifs

@Author: Ying Huang
"""
from collections import deque
import re

from .motif_format import MEMEmotifFormat


class MEMEmotifsDict():

    """
    Store a Motif Database in MEMEmotifsDict() object.
    Version, alphabet, strands is stored in string format respectivly.
    Motifs part is stored in dict format, which each values is a single motif
     part in MEMEmotifFormat() class format.

    Example:
    # create MEMEmotifsDict() object.
	>>> motif_db = MEMEmotifsDict('/usr/local/bin/meme/db/motif_databases/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme')

    # Version part
	>>> motif_db.version
    
    # Frequencies part
	>>> motif_db.frequencies

    # Alphabet part
	>>> motif_db.alphabet

    # Strands part
	>>> motif_db.strands
	
    # Search motif wanted and output
    >>> motif_db.search_motifs(list(motif_db.motif.keys())[:5])
    >>> motif_db.motif_searched
	>>> motif_db.searched_motif_output('./motif.test.meme')
    """

    Version_Marker = 'MEME version'

    Alphabet_Marker = 'ALPHABET='

    Strands_Marker = 'strands:'

    Frequencies_Start_Marker = 'Background letter frequencies'
    Frequencies_End_Marker = r'A [01]\.[0-9]+ C [01]\.[0-9]+ G [01]\.[0-9]+ T [01]\.[0-9]+'

    Motif_start_Marker = 'MOTIF'
    Motif_end_Marker = 'URL' 
    

    def __init__(self, meme_ifile):

        """meme_ifile: <str> path to meme motif database."""

        self._ifile = meme_ifile

        # Five main compartments of MEME format Motif.
        self.version = None  # Required by MEME format Motif.
        self.alphabet = None  # Optional by MEME format Motif.
        self.strands = None  # Optional by MEME format Motif.
        self.frequencies = None  # Optional by MEME format Motif.
        self.motif = dict()  # Required by MEME format Motif.

        self.tmp_motif = deque()
        self._is_read_motif = False

        # After running 'self._make_motif_dict()', Five main compartments of
        # MEME format Motif are refreshed by input motif database.
        self._make_motif_dict()

        self.motif_searched = dict()

    def _make_motif_dict(self):

        with open(self._ifile, 'rt') as f:
            print("> Reading file: {}".format(self._ifile))
            for line in f:
                if line.strip('\n') is None:
                    #print('NoneType line')
                    # Note: 'URL' row of Motif lines is optional, that means some motif
                    # lines don't end with 'URL' row, instead of 'matrix' row.
                    # Thus, if '_is_read_motif' is True and current line is NoneType,
                    # we determine current line the end of a motif.
                    if self._is_read_motif is True:
                        #print('*** Motif line end with None')
                        # Current line is in the URL part (ending row) of motif.
                        # And this motif don't have 'url' part.
                        # 1. Now 'tmp_motif' has recieved a complete motif, and
                        # should be stored in 'self.motif' (a dict).
                        a_motif = MEMEmotifFormat(''.join(self.tmp_motif))
                        # If motif ID is reduplicative, raise error.
                        if a_motif.motif_id in self.motif.keys():
                            raise Exception(
                                "ERR: Get reduplicative motif ID '{}', which can not be a key of dict.".format(
                                    a_motif.motif_id)
                            )
                        self.motif[a_motif.motif_id] = a_motif
                        # Clear tmp_motif for getting new line.
                        self.tmp_motif.clear()
                        # 3. 'self._is_read_motif' should change to False, to 
                        # stop read in line.
                        self._is_read_motif = False
                    # Skip empty line
                    continue
                elif line.startswith(self.Version_Marker):
                    #print('Version line')
                    self.version = line
                elif line.startswith(self.Alphabet_Marker):
                    #print('Alphabet line')
                    self.alphabet = line
                elif line.startswith(self.Strands_Marker):
                    #print('Strand line')
                    self.strands = line
                # Frequencies contains two rows.
                elif line.startswith(self.Frequencies_Start_Marker):
                    #print('freq line start')
                    # Current line is at first row of Frequencies.
                    self.frequencies = line
                elif re.match(self.Frequencies_End_Marker, line):
                    #print('freq line end')
                    # Current line is at second row of Frequencies.
                    self.frequencies += line
                # Motif consist of three parts: name (required), matri
                # (required), url (optional). 
                elif line.startswith(self.Motif_start_Marker):
                    #print('Motif line start')
                    # When meeting a line starts with 'Motif_start_Marker', it 
                    # means current line is the start line of motif.
                    # 1. Current line can be read in as part of 'tmp_motif'.
                    self.tmp_motif.append(line)
                    # 2. 'self._is_read_motif' should change to True, to start
                    # read in line.
                    self._is_read_motif = True
                elif self._is_read_motif is True:
                    #print('Motif line matrix (or URL)')
                    # Current line is in the matrix (or URL) part of motif.
                    # Current line can be read in as part of 'tmp_motif'.
                    self.tmp_motif.append(line)

                    if line.startswith(self.Motif_end_Marker):
                        #print('Motif line URL')
                        # Current line is in the URL part (ending row) of motif.
                        # 1. Now 'tmp_motif' has recieved a complete motif, and 
                        # should be stored in 'self.motif' (a dict).
                        
                        a_motif = MEMEmotifFormat(''.join(self.tmp_motif))
            
                        # If motif ID is reduplicative, raise error.
                        if a_motif.motif_id in self.motif.keys():
                            raise Exception(
                                "ERR: Get reduplicative motif ID '{}', which can not be a key of dict.".format(a_motif.motif_id)
                            )
                        self.motif[a_motif.motif_id] = a_motif
                        # Clear tmp_motif for getting new line.
                        self.tmp_motif.clear() 
                        # 3. 'self._is_read_motif' should change to False, to stop
                        # read in line.
                        self._is_read_motif = False

    def search_motifs(self, motif_ids):
        """motif_ids: <iter> a list of motif IDs."""

        self.motif_searched.clear()
        
        for motif_id in motif_ids:
            self.motif_searched[motif_id] = self.motif[motif_id]
        
        return self

    def searched_motif_output(self, ofile):
        
        header = '\n'.join([self.version, self.alphabet, self.strands, self.frequencies])

        with open(ofile, 'wt') as f:
            f.write(header)
        
        out_motif = deque()

        for motif in self.motif_searched.values():
            out_motif.append(str(motif))

        with open(ofile, 'at') as f:
            out = '\n' + '\n'.join(out_motif)
            f.write(out)

