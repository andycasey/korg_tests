import re
from grok.transitions.formula import Formula

roman_numerals = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"]

class Species(object):

    def __init__(self, representation):

        representation = representation.strip(" 0")
        tokens = re.split("[\s\._]", representation)
        if len(tokens) > 2:
            raise ValueError(f"Invalid species code '{representation}' ({type(representation)})")
        
        self.formula = Formula(tokens[0])

        if len(tokens) == 1 or len(tokens[1]) == 0:
            charge = 0
        else:
            if tokens[1] in roman_numerals:
                charge = roman_numerals.index(tokens[1])
            elif tokens[1].isdigit():
                # Take the first digit only.
                charge = int(tokens[1][0])
                # TODO: Do anything with the rest of the information? 
                #       MOOG encodes isotopic ratios in this way.
            else:
                raise ValueError(f"Invalid charge '{tokens[1]}' in '{representation}'")

        self.charge = charge

        return None

    def __repr__(self):
        return f"<{''.join(self.atoms)} {roman_numerals[self.charge]}>"


    @property
    def atoms(self):
        return self.formula.atoms