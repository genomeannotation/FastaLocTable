#!/usr/bin/env python

import sys

class Sequence:

    def __init__(self, header="", bases=""):
        self.header = header
        self.bases = bases

    def __str__(self):
        result = "Sequence " + self.header
        result += " of length " + str(len(self.bases))
        return result
