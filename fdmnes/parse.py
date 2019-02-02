from os import linesep
from . import settings


def param2str(param):
    if isinstance(param, bool):
        return ""
    if isinstance(param, str):
        return param
    if isinstance(param, int):
        return "%i"%param
    if isinstance(param, float):
        return "%g"%param
    if isinstance(param, list):
        assert all(map(lambda x: not isinstance(x, list), param))
        return linesep.join(map(param2str, param))
    if isinstance(param, tuple):
        assert all(map(lambda x: isinstance(x, (bool,str,int,float,tuple)),
                       param))
        return " ".join(map(param2str, param))

def mkfloat(string):
    i = string.find("(")
    if i>=0:
        string = string[:i]
    return float(string)


def keyword_exists(key):
    key = key.lower()
    keyset = None
    if key in settings.synonyms:
        return settings.synonyms[key]
    for Group in settings.Defaults:
        keys = Group.keys()
        keysl = map(str.lower, keys)
        keysmask = [s==key for s in keysl]
        if any(keysmask):
            i = keysmask.index(True)
            keyset = list(keys)[i]
            break
    return keyset



def parse_input_file(fpath, outpath_only=False):
    with open(fpath, "r") as fh:
        content = fh.readlines()

    content = map(str.strip, content)
    content = filter(lambda s: not s.startswith("!"), content)
    content = map(lambda s: s.split("!")[0], content)
    content = map(str.strip, content)
    content = filter(lambda s: bool(s), content)
    content = list(content) #python3

    Param = dict()
    keyw = None
    sg = None
    conv = False
    while content:
        line = content.pop(0)
        lline = line.lower()
        if lline == "filout" or lline=="conv_out":
            path_out = content.pop(0)
            if outpath_only:
                return path_out
            if lline=="conv_out":
                conv = True
        elif lline == "spgroup":
            sg = content.pop(0)
        elif lline.startswith("crystal") or lline.startswith("molecule"):
            keyw = lline.capitalize()
            Param[keyw] = []
        elif lline=="end":
            break
        else:
            nextkeyw = keyword_exists(line)
            if nextkeyw:
                keyw = nextkeyw
                Param[keyw] = []
            elif keyw!=None:
                Param[keyw].append(line)

    structure_keyw = ["Crystal", "Crystal_t", "Crystal_p",
                      "Molecule", "Molecule_t"]

    found_structure = [(k in Param) for k in structure_keyw]

    if sum(found_structure) > 1:
        raise ConsistencyError(
            "Several structure definitions found in input-file. "\
            "Structure can be specified in input-file by the "\
            "keywords Crystal or Molecule, not both!")
    elif sum(found_structure) == 1:
        ind = found_structure.index(True)
        structure = dict(structure_keyw=Param.pop(structure_keyw[ind]))
    elif not any(found_structure):
        #raise IOError("No structure has been defined in input-file")
        structure = dict()

    if conv:
        bavfile = ""
    elif "Extract" in Param and Param["Extract"][0]:
        bavfile = Param["Extract"][0]
    else:
        bavfile = path_out + "_bav.txt"


    return dict(Param=Param, sg=sg, path_out=path_out, conv=conv,
                bavfile=bavfile, structure=structure)



def parse_bavfile(bavfile):
    bavinfo = {}
    with open(bavfile, "rb") as bf:
        bavcont = bf.read()
        bf.seek(-60, 2)
        tail = [s.decode().strip() for s in bf.readlines()]
        #tail = map(str.strip, bf.readlines())
    bavinfo["success"] = (tail[-1]=='Have a beautiful day !')
    bavinfo["num_absorber"] = bavcont.count(b"Index of the absorbing atom")

    return bavinfo
