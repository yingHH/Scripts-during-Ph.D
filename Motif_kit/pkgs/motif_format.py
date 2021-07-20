# -*- coding: utf-8 -*-
"""
Created on 2018-10-19 19:52:15
Last Modified on 2018-10-19 19:52:15

Classes for all kinds of motif format

@Author: Ying Huang
"""
import re


class MEMEmotifFormat():
    """
    Object to format motif in MEME format.

    A motif record in MEME format:
        MOTIF MA0004.1 Arnt
        letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0
         0.200000  0.800000  0.000000  0.000000
         0.950000  0.000000  0.050000  0.000000
         0.000000  1.000000  0.000000  0.000000
         0.000000  0.000000  1.000000  0.000000
         0.000000  0.000000  0.000000  1.000000
         0.000000  0.000000  1.000000  0.000000
        URL http://jaspar2018.genereg.net/matrix/MA0004.1

    Example:
        >>> a_motif = \
        '''OTIF MA0004.1 Arnt
        letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0
         0.200000  0.800000  0.000000  0.000000
         0.950000  0.000000  0.050000  0.000000
         0.000000  1.000000  0.000000  0.000000
         0.000000  0.000000  1.000000  0.000000
         0.000000  0.000000  0.000000  1.000000
         0.000000  0.000000  1.000000  0.000000
        URL http://jaspar2018.genereg.net/matrix/MA0004.1'''

        >>> motif = MEMEmotifFormat(a_motif)

        >>> motif
            MOTIF MA0004.1 Arnt
            letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0
             0.200000  0.800000  0.000000  0.000000
             0.950000  0.000000  0.050000  0.000000
             0.000000  1.000000  0.000000  0.000000
             0.000000  0.000000  1.000000  0.000000
             0.000000  0.000000  0.000000  1.000000
             0.000000  0.000000  1.000000  0.000000
            URL http://jaspar2018.genereg.net/matrix/MA0004.1

        >>> motif.motif_id
            'MA0004.1'

        >>> motif.motif_name
            'Arnt'

        >>> motif.name
            'MOTIF MA0004.1 Arnt'

        >>> motif.matrix
            'letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0\n 0.200000  0.800000  0.000000  0.000000\n 0.950000  0.000000  0.050000  0.000000\n 0.000000  1.000000  0.000000  0.000000\n 0.000000  0.000000  1.000000  0.000000\n 0.000000  0.000000  0.000000  1.000000\n 0.000000  0.000000  1.000000  0.000000\n'

        >>> motif.url
        'URL http://jaspar2018.genereg.net/matrix/MA0004.1'

    """

    def __init__(self, a_motif):
        """a_motif: <str> A motif record in MEME format."""
        self._motif = a_motif

        self._reg_motif_name = re.search(
            r'MOTIF ([A-Z]+[0-9]+\.[0-9]+) ([A-Za-z0-9:]+)',
            self._motif
        )

        self._reg_letter_matrix = re.search(
            r'(letter-probability matrix: alength= [0-9]+ w= [0-9]+ nsites= [0-9]+ E= [0-9]+)\n(((\s+[01]\.[0-9]+){4}\n)+)',
            self._motif
        )

        self._reg_URL = re.search(
            r'URL (http://.+)',
            self._motif
        )

    def __str__(self):
        return str(self._motif)

    def __repr__(self):
        return str(self._motif)

    @property
    def name(self):
        try:
            return self._reg_motif_name.group(0)
        except:
            raise Exception(
                "Please check the 'name' row of your motif record below:\n{}\n".format(self._motif))

    @property
    def motif_id(self):
        try:
            return self._reg_motif_name.group(1)
        except:
            raise Exception(
                "Please check 'motif id' in the 'name' row of your motif record below:\n{}\n".format(self._motif))

    @property
    def motif_name(self):
        try:
            return self._reg_motif_name.group(2)
        except:
            raise Exception(
                "Please check 'motif name' in the 'name' row of your motif record below:\n{}\n".format(self._motif))

    @property
    def matrix(self):
        try:
            return str(self._reg_letter_matrix.group(0))
        except:
            raise Exception(
                "Please check the 'matrix' rows of your motif record below:\n{}\n".format(self._motif))

    @property
    def url(self):
        try:
            return self._reg_URL.group(0)
        except:
            print("No URL in Motif record.")
