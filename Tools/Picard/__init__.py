#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
from Tools.Picard.SortVcf import SortVcf
from Tools.Picard.MarkDuplicates import MarkDuplicates
from Tools.Picard.AddOrReplaceReadGroups import AddOrReplaceReadGroups
from Tools.Picard.CreateSequenceDictionary import CreateSequenceDictionary

SortVcf = SortVcf()
MarkDuplicates = MarkDuplicates()
AddOrReplaceReadGroups = AddOrReplaceReadGroups()
CreateSequenceDictionary = CreateSequenceDictionary()
