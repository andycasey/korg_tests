from grok.utils import periodic_table

class Formula(object):
    
    _Nz = 3 # support up to tri-atomic molecules

    def __init__(
            self, 
            representation
        ):

        # We assume that if a tuple is given, you (mostly) know what you're doing.
        if isinstance(representation, (tuple, list)):
            Zs = list(representation)
        else:
            # Try various things.
            if isinstance(representation, int):
                if 1 <= representation < len(periodic_table):
                    Zs = [representation]
                else:
                    raise ValueError("Must be some molecule, but you should use string representation.")
    
            elif isinstance(representation, str):
                if representation.isdigit():
                    # Numeric codes like 0801 -> OH
                    L = len(representation)
                    if L <= 2:
                        Zs = [int(representation)]
                    elif L <= 4:
                        code = f"{representation:0>4s}"
                        Z1 = int(code[0:2])
                        Z2 = int(code[2:4])
                        Zs = [Z1, Z2]
                    else:
                        raise ValueError(f"Numeric codes for molecules with more than 4 chars like "
                                         f"'{representation}' are not supported.")
                else:
                    # Should be something like "OH", "FeH", "Li", "C2"
                    indices = [i for i, char in enumerate(representation) if char.isdigit() or char.isupper()]
                    indices.append(len(representation))
                    codes = []
                    for j, index in enumerate(indices[:-1]):
                        codes.append(representation[index:indices[j+1]])
            
                    Zs = []
                    for sub_code in codes:
                        if sub_code.isnumeric():
                            previous_Z = Zs[-1]
                            for j in range(int(sub_code) - 1):
                                Zs.append(previous_Z)
                        else:
                            Zs.append(periodic_table.index(sub_code) + 1)
            else:
                raise TypeError(f"Invalid type representation ({type(representation)}) for '{representation}'")

        # Prepend with zeros.
        Zs = (([0] * self._Nz) + Zs)[-self._Nz:]

        # TODO: Should we sort here or just in time for outputting?
        self.Zs = tuple(Zs)
        return None


    def __repr__(self):
        return f"<{''.join(self.atoms)}>"


    @property
    def atoms(self):
        return tuple([periodic_table[Z-1] for Z in self.Zs if Z > 0])