__author__ = 'mahajrod'

from Data.Nucleotides import degenerated_nucleotides


class PrimerRoutines:
    
    def count_GC_AT_composition(self, primer_sequence):
        GC_count_min = 0
        GC_count_max = 0
        AT_count_min = 0
        AT_count_max = 0

        for nucleotide in str(primer_sequence).upper():
            if nucleotide == "G" or nucleotide == "C":
                GC_count_max += 1
                GC_count_min += 1
            elif nucleotide == "A" or nucleotide == "T" or nucleotide == "U":
                AT_count_max += 1
                AT_count_min += 1
            else:
                if nucleotide not in degenerated_nucleotides:
                    raise ValueError("Wrong nucleotide")
                else:
                    if nucleotide == "S":
                        GC_count_max += 1
                        GC_count_min += 1
                    elif nucleotide == "W":
                        AT_count_max += 1
                        AT_count_min += 1
                    else:
                        GC_count_max += 1
                        AT_count_max += 1
        return {"GC": [GC_count_min, GC_count_max],
                "AT": [AT_count_min, AT_count_max]}

    def melting_temperature_formula1(self, primer_sequence):
        # Tm = (wA+xT) * 2 + (yG+zC) * 4
        GC_AT_composition = self.count_GC_AT_composition(primer_sequence)

        max_temp = GC_AT_composition["GC"][1]*4 + GC_AT_composition["AT"][0]*2
        min_temp = GC_AT_composition["GC"][0]*4 + GC_AT_composition["AT"][1]*2

        return min_temp, max_temp

    def melting_temperature_formula2(self, primer_sequence):
        # Tm = 64.9 +41*(yG+zC-16.4)/(wA+xT+yG+zC)
        GC_AT_composition = self.count_GC_AT_composition(primer_sequence)

        max_temp = 64.9 + 41.0 * (float(GC_AT_composition["GC"][1]) - 16.4)/(float(len(primer_sequence)))
        min_temp = 64.9 + 41.0 * (float(GC_AT_composition["GC"][0]) - 16.4)/(float(len(primer_sequence)))

        return min_temp, max_temp
