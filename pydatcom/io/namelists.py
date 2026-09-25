"""
DATCOM's 23 input namelists: XNAM1 to XNAM23.

Each XNAM routine defines one namelist (its name, variables, dimensions and
storage positions), gathers the COMMON words it reads into a buffer, calls
NAMER (``IOP`` 0, read) or NAMEW (``IOP`` 1, echo) on the buffer, scatters
the buffer back, and for some namelists finishes with a few derived words.

COMMON blocks are passed as 1-based lists keyed by the COMMON name, with the
words numbered as the XNAM routines declare them (so ``IBW`` word 1 is the
block's leading flag word).  ``VTI`` must carry 324 words: XNAM21 copies
eight words past the 316 its COMMON declares, into whatever follows the
block, and back.

The definitions were extracted from the sources' DATA statements by
parsing.

Reference: datcom-legacy/datcom_2000/xnam1.f ... xnam23.f
"""

from typing import IO, Dict, List, MutableSequence, Optional, Tuple

from pydatcom.io.fortran_format import fortran_write
from pydatcom.io.namelist_io import namew
from pydatcom.io.namelist_reader import namer
from pydatcom.utils.constants import KAND, UNUSED

Blocks = Dict[str, MutableSequence]

NAMELISTS: Dict[int, dict] = {}
NAMELISTS.update(
    {1: {'name': 'FLTCON',
         'unit': 9,
         'buffer': 'A2A',
         'size': 161,
         'names': ['NMACH', 'MACH', 'NALPHA', 'ALSCHD', 'RNNUB', 'HYPERS',
                   'STMACH', 'TSMACH', 'TR', 'ALT', 'PINF', 'TINF', 'VINF',
                   'WT', 'GAMMA', 'NALT', 'LOOP', 'ALPHA'],
         'dims': [1, 20, 1, 20, 20, -1, 1, 1, 1, 20, 20, 20, 20, 1, 1, 1, 1,
                  20],
         'locs': [1, 3, 2, 23, 43, 161, 94, 95, 96, 97, 74, 117, 137, 157, 158,
                  159, 160, 23]},
     2: {'name': 'OPTINS',
         'unit': 9,
         'buffer': 'A3',
         'size': 4,
         'names': ['SREF', 'CBARR', 'ROUGFC', 'BLREF'],
         'dims': [1, 1, 1, 1],
         'locs': [1, 2, 3, 4]},
     3: {'name': 'BODY',
         'unit': 9,
         'buffer': 'A4',
         'size': 129,
         'names': ['NX', 'X', 'S', 'P', 'R', 'ZU', 'ZL', 'BNOSE', 'BTAIL',
                   'BLN', 'BLA', 'DS', 'ITYPE', 'METHOD', 'ELLIP'],
         'dims': [1, 20, 20, 20, 20, 20, 20, 1, 1, 1, 1, 1, 1, 1, 1],
         'locs': [1, 2, 22, 42, 62, 82, 102, 122, 123, 124, 125, 126, 127, 128,
                  129]},
     4: {'name': 'WGPLNF',
         'unit': 9,
         'buffer': 'A',
         'size': 15,
         'names': ['CHRDBP', 'CHRDR', 'CHRDTP', 'CHSTAT', 'SSPN', 'SSPNE',
                   'SSPNOP', 'SAVSI', 'SAVSO', 'SWAFP', 'TWISTA', 'TYPE',
                   'SSPNDD', 'DHDADI', 'DHDADO'],
         'dims': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
         'locs': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
         'iequ': [5, 6, 1, 9, 4, 3, 2, 7, 8, 10, 11, 15, 12, 13, 14]},
     5: {'name': 'WGSCHR',
         'unit': 9,
         'buffer': 'A6',
         'size': 253,
         'names': ['TOVC', 'DELTAY', 'XOVC', 'CLI', 'ALPHAI', 'CLALPA', 'CLMAX',
                   'CAMBER', 'CM0', 'XOVCO', 'CM0T', 'LERI', 'LERO', 'TOVCO',
                   'CMO', 'CMOT', 'TCEFF', 'KSHARP', 'CLMAXL', 'SLOPE', 'CLAMO',
                   'CLAM0', 'ARCL', 'XAC', 'DWASH', 'YCM', 'CLD', 'TYPEIN',
                   'NPTS', 'XCORD', 'YUPPER', 'YLOWER', 'MEAN', 'THICK',
                   'ALPHAO', 'ALPHA0'],
         'dims': [1, 1, 1, 1, 1, 20, 20, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 6,
                  1, 1, 1, 20, 1, 1, 1, 1, 1, 50, 50, 50, 50, 50, 1, 1],
         'locs': [16, 17, 18, 19, 20, 21, 41, 64, 61, 66, 67, 62, 63, 65, 61,
                  67, 70, 71, 68, 95, 69, 69, 92, 72, 101, 93, 94, 253, 102,
                  103, 153, 203, 153, 203, 10, 10]},
     6: {'name': 'SYNTHS',
         'unit': 9,
         'buffer': 'A7',
         'size': 19,
         'names': ['XCG', 'XW', 'ZW', 'ALIW', 'ZCG', 'XH', 'ZH', 'ALIH', 'XV',
                   'VERTUP', 'HINAX', 'XVF', 'SCALE', 'ZV', 'ZVF', 'YV', 'YF',
                   'PHIV', 'PHIF'],
         'dims': [1, 1, 1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
         'locs': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
                  19]},
     7: {'name': 'HTPLNF',
         'unit': 9,
         'buffer': 'A',
         'size': 75,
         'names': ['CHRDBP', 'CHRDR', 'CHRDTP', 'CHSTAT', 'SSPN', 'SSPNE',
                   'SSPNOP', 'SAVSI', 'SAVSO', 'SWAFP', 'TWISTA', 'TYPE',
                   'SSPNDD', 'DHDADI', 'DHDADO', 'RLPH', 'SHB', 'SEXT'],
         'dims': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 20, 20, 20],
         'locs': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 36,
                  56],
         'iequ': [5, 6, 1, 9, 4, 3, 2, 7, 8, 10, 11, 15, 12, 13, 14, 95, 115,
                  135],
         'xtype': ['STRA', 'DDUB', 'CRAN', 'CURV']},
     8: {'name': 'HTSCHR',
         'unit': 9,
         'buffer': 'HS',
         'size': 306,
         'names': ['TOVC', 'DELTAY', 'XOVC', 'CLI', 'ALPHAI', 'CLALPA', 'CLMAX',
                   'CAMBER', 'CM0', 'XOVCO', 'CM0T', 'LERI', 'LERO', 'TOVCO',
                   'CMO', 'CMOT', 'TCEFF', 'KSHARP', 'CLAM0', 'CLAMO', 'ARCL',
                   'CLMAXL', 'YCM', 'CLD', 'XAC', 'TYPEIN', 'NPTS', 'XCORD',
                   'YUPPER', 'YLOWER', 'MEAN', 'THICK', 'ALPHAO', 'ALPHA0'],
         'dims': [1, 1, 1, 1, 1, 20, 20, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  1, 1, 1, 1, 20, 1, 1, 50, 50, 50, 50, 50, 1, 1],
         'locs': [16, 17, 18, 19, 20, 21, 41, 64, 61, 66, 67, 62, 63, 65, 61,
                  67, 70, 71, 69, 69, 92, 68, 93, 94, 72, 306, 155, 156, 206,
                  256, 206, 256, 10, 10]},
     9: {'name': 'VTPLNF',
         'unit': 9,
         'buffer': 'VT',
         'size': 75,
         'names': ['CHRDTP', 'SSPNOP', 'SSPNE', 'SSPN', 'CHRDBP', 'CHRDR',
                   'SAVSI', 'CHSTAT', 'SWAFP', 'TWISTA', 'SSPNDD', 'DHDADI',
                   'DHDADO', 'TYPE', 'SAVSO', 'SVWB', 'SVB', 'SVHB'],
         'dims': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 20, 20, 20],
         'locs': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 36,
                  56],
         'iequ': [1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 8, 95, 115,
                  135],
         'xtype': ['STRA', 'DDUB', 'CRAN', 'CURV']},
     10: {'name': 'VTSCHR',
          'unit': 9,
          'buffer': 'VTS',
          'size': 314,
          'names': ['TOVC', 'DELTAY', 'XOVC', 'CLI', 'ALPHAI', 'CLALPA',
                    'CLMAX', 'CAMBER', 'CM0', 'XOVCO', 'CM0T', 'LERI', 'LERO',
                    'TOVCO', 'CMO', 'CMOT', 'TCEFF', 'KSHARP', 'CLAM0', 'CLAMO',
                    'ARCL', 'CLMAXL', 'YCM', 'CLD', 'XAC', 'TYPEIN', 'NPTS',
                    'XCORD', 'YUPPER', 'YLOWER', 'MEAN', 'THICK', 'ALPHAO',
                    'ALPHA0'],
          'dims': [1, 1, 1, 1, 1, 20, 20, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 20, 1, 1, 50, 50, 50, 50, 50, 1, 1],
          'locs': [16, 17, 18, 19, 20, 21, 41, 64, 61, 66, 67, 62, 63, 65, 61,
                   67, 70, 71, 69, 69, 92, 68, 155, 155, 72, 314, 163, 164, 214,
                   264, 214, 264, 10, 10]},
     11: {'name': 'PROPWR',
          'unit': 9,
          'buffer': 'AC',
          'size': 50,
          'names': ['AIETLP', 'NENGSP', 'THSTCP', 'PHALOC', 'PHVLOC', 'PRPRAD',
                    'ENGFCT', 'BWAPR3', 'BWAPR6', 'BWAPR9', 'NOPBPE', 'BAPR75',
                    'CROT', 'YP'],
          'dims': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, 1],
          'locs': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 29, 28]},
     12: {'name': 'JETPWR',
          'unit': 9,
          'buffer': 'AD',
          'size': 50,
          'names': ['AIETLJ', 'NENGSJ', 'THSTCJ', 'JIALOC', 'JEVLOC', 'JEALOC',
                    'JINLTA', 'JEANGL', 'JEVELO', 'AMBTMP', 'JESTMP', 'JELLOC',
                    'JETOTP', 'AMBSTP', 'JERAD'],
          'dims': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
          'locs': [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]},
     13: {'name': 'LARWB',
          'unit': 9,
          'buffer': 'AE',
          'size': 50,
          'names': ['ZB', 'SREF', 'DELTEP', 'SFRONT', 'AR', 'R3LEOB', 'DELTAL',
                    'L', 'SWET', 'PERBAS', 'SBASE', 'HB', 'BB', 'BLF', 'XCG',
                    'THETAD', 'ROUNDN', 'SBS', 'SBSLB', 'XCENSB', 'XCENW'],
          'dims': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, 1, 1, -1, 1, 1, 1,
                   1],
          'locs': [30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44,
                   45, 46, 47, 48, 49, 50]},
     14: {'name': 'GRNDEF',
          'unit': 9,
          'buffer': 'AF',
          'size': 141,
          'names': ['NGH', 'GRDHT'],
          'dims': [1, 10],
          'locs': [63, 64]},
     15: {'name': 'TVTPAN',
          'unit': 9,
          'buffer': 'AG',
          'size': 162,
          'names': ['BVP', 'BV', 'BDV', 'BH', 'SV', 'VPHITE', 'VLP', 'ZP'],
          'dims': [1, 1, 1, 1, 1, 1, 1, 1],
          'locs': [155, 156, 157, 158, 159, 160, 161, 162]},
     16: {'name': 'SYMFLP',
          'unit': 9,
          'buffer': 'AIA',
          'size': 117,
          'names': ['CHRDFI', 'CHRDFO', 'SPANFI', 'SPANFO', 'NDELTA', 'PHETEP',
                    'PHETE', 'FTYPE', 'NTYPE', 'SCHA', 'CB', 'TC', 'SCHD',
                    'DELTA', 'CPRMEI', 'CPRMEO', 'SCLD', 'SCMD', 'CMU',
                    'DELJET', 'JETFLP', 'EFFJET', 'CAPINB', 'CAPOUT', 'DOBDEF',
                    'DOBCIN', 'DOBCOT'],
          'dims': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 9, 9, 9, 9, 9, 1, 9,
                   1, 9, 9, 9, 9, 1, 1],
          'locs': [12, 13, 14, 15, 16, 61, 11, 17, 62, 18, 59, 60, 18, 1, 39,
                   49, 19, 29, 63, 64, 74, 75, 85, 95, 105, 115, 116]},
     17: {'name': 'ASYFLP',
          'unit': 9,
          'buffer': 'AJ',
          'size': 137,
          'names': ['DELTAL', 'DELTAR', 'DELTAD', 'DELTAS', 'XSOC', 'HSOC',
                    'STYPE', 'XSPRME', 'NDELTA', 'CHRDFI', 'CHRDFO', 'SPANFI',
                    'SPANFO', 'PHETE'],
          'dims': [9, 9, 9, 9, 9, 9, 1, 1, 1, 1, 1, 1, 1, 1],
          'locs': [19, 29, 1, 39, 49, 60, 18, 59, 16, 12, 13, 14, 15, 11]},
     18: {'name': 'HYPEFF',
          'unit': 9,
          'buffer': 'AK',
          'size': 137,
          'names': ['ALITD', 'XHL', 'TWOTI', 'CF', 'LAMNR', 'HNDLTA', 'HDELTA'],
          'dims': [1, 1, 1, 1, -1, 1, 10],
          'locs': [1, 2, 3, 4, 15, 16, 5]},
     19: {'name': 'TRNJET',
          'unit': 9,
          'buffer': 'AL',
          'size': 137,
          'names': ['TIME', 'FC', 'ALPHA', 'LAMNRJ', 'NT', 'ME', 'ISP', 'SPAN',
                    'PHE', 'GP', 'CC', 'LFP'],
          'dims': [10, 10, 10, -10, 1, 1, 1, 1, 1, 1, 1, 1],
          'locs': [1, 12, 22, 39, 11, 32, 33, 34, 35, 36, 37, 38]},
     20: {'name': 'VFPLNF',
          'unit': 9,
          'buffer': 'VFP',
          'size': 75,
          'names': ['CHRDTP', 'SSPNOP', 'SSPNE', 'SSPN', 'CHRDBP', 'CHRDR',
                    'SAVSI', 'CHSTAT', 'SWAFP', 'TWISTA', 'SSPNDD', 'DHDADI',
                    'DHDADO', 'TYPE', 'SAVSO', 'SVWB', 'SVB', 'SVHB'],
          'dims': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 20, 20, 20],
          'locs': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 36,
                   56],
          'iequ': [163, 164, 165, 166, 167, 168, 169, 171, 172, 173, 174, 175,
                   176, 177, 170, 257, 277, 297],
          'xtype': ['STRA', 'DDUB', 'CRAN', 'CURV']},
     21: {'name': 'VFSCHR',
          'unit': 9,
          'buffer': 'VFS',
          'size': 476,
          'names': ['TOVC', 'DELTAY', 'XOVC', 'CLI', 'ALPHAI', 'CLALPA',
                    'CLMAX', 'CAMBER', 'CM0', 'XOVCO', 'CM0T', 'LERI', 'LERO',
                    'TOVCO', 'CMO', 'CMOT', 'TCEFF', 'KSHARP', 'CLAM0', 'CLAMO',
                    'ARCL', 'CLMAXL', 'YCM', 'CLD', 'XAC', 'TYPEIN', 'NPTS',
                    'XCORD', 'YUPPER', 'YLOWER', 'MEAN', 'THICK', 'ALPHAO',
                    'ALPHA0'],
          'dims': [1, 1, 1, 1, 1, 20, 20, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 20, 1, 1, 50, 50, 50, 50, 50, 1, 1],
          'locs': [178, 179, 180, 181, 182, 183, 203, 226, 223, 228, 229, 224,
                   225, 227, 223, 229, 232, 233, 231, 231, 254, 230, 0, 0, 234,
                   476, 325, 326, 376, 426, 376, 426, 10, 10]},
     22: {'name': 'CONTAB',
          'unit': 9,
          'buffer': 'AM',
          'size': 137,
          'names': ['TTYPE', 'CFITC', 'CFOTC', 'CFITT', 'CFOTT', 'BITC', 'BOTC',
                    'BITT', 'BOTT', 'B1', 'B2', 'B3', 'B4', 'D1', 'D2', 'D3',
                    'GCMAX', 'KS', 'RL', 'BGR', 'DELR'],
          'dims': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1],
          'locs': [117, 118, 119, 122, 123, 120, 121, 124, 125, 126, 127, 128,
                   129, 130, 131, 132, 133, 134, 135, 136, 137]},
     23: {'name': 'EXPR',
          'unit': 10,
          'buffer': 'AH',
          'size': 469,
          'names': ['CDB', 'CLB', 'CMB', 'CLAB', 'CMAB', 'CDW', 'CLW', 'CMW',
                    'CLAW', 'CMAW', 'CDH', 'CLH', 'CMH', 'CLAH', 'CMAH', 'CDWB',
                    'CLWB', 'CMWB', 'CLAWB', 'CMAWB', 'QOQINF', 'EPSLON',
                    'DEODA', 'CDV', 'ALPOW', 'ALPLW', 'ALPOH', 'ALPLH', 'ACLMW',
                    'CLMW', 'ACLMH', 'CLMH'],
          'dims': [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
                   20, 20, 20, 20, 20, 20, 20, 20, 1, 1, 1, 1, 1, 1, 1, 1, 1],
          'locs': [1, 21, 41, 61, 81, 101, 121, 141, 161, 181, 201, 221, 241,
                   261, 281, 301, 321, 341, 361, 381, 401, 421, 441, 461, 462,
                   463, 464, 465, 466, 467, 468, 469]}}
)

# Namelists read straight into a COMMON block: number -> (block, words).
_DIRECT = {2: ('OPTION', 4), 3: ('BODYI', 129), 6: ('SYNTSS', 19),
           11: ('POWER', 50), 12: ('POWER', 50), 13: ('POWER', 50),
           14: ('FLGTCD', 141), 15: ('VTI', 162), 17: ('FLAPIN', 137),
           18: ('FLAPIN', 137), 19: ('FLAPIN', 137), 22: ('FLAPIN', 137)}
_TYPE_ERROR = '(42H ERROR TYPE CANNOT BE GREATER THAN 4 TYPE=,E12.5,' \
    '9H SET TO 1)'

Link = Tuple[int, str, int]


def _planform(block: str, iequ, skip: int) -> Tuple[List[Link],
                                                     List[Link]]:
    """XNAM7, 9 and 20: fifteen scalars through ``IEQU``, then three
    20-word arrays; the write-back skips the word ``skip``."""
    links = [(buffer_word, block, iequ[buffer_word - 1])
             for buffer_word in range(1, 16)]
    for array_base, base_offset in ((16, -1), (17, 18), (18, 37)):
        links += [(array_base + element + base_offset, block,
                   iequ[array_base - 1] + element - 1)
                  for element in range(1, 21)]
    back = [ln for ln in links if not (ln[0] <= 15 and ln[2] == skip)]
    return links, back


def _links(n: int) -> Tuple[List[Link], List[Link]]:
    """The buffer's links to COMMON words, as ``(buffer word, block,
    block word)``: those gathered, and those written back."""
    d = NAMELISTS[n]
    if n in _DIRECT:
        block, words = _DIRECT[n]
        links = [(buffer_word, block, buffer_word)
                 for buffer_word in range(1, words + 1)]
        return links, links
    if n == 1:
        links = [(buffer_word, 'FLGTCD', buffer_word)
                 for buffer_word in range(1, 161)] + \
            [(161, 'FLOLOG', 19)]
        return links, links
    if n == 4:
        links = [(buffer_word, 'WINGI', d['iequ'][buffer_word - 1])
                 for buffer_word in range(1, 16)]
        return links, [ln for ln in links if ln[2] != 10]
    if n in (5, 8, 10, 21):
        top, other, shift = {5: (101, 'IBW', 111), 8: (154, 'IBH', 58),
                             10: (162, 'IBV', 50),
                             21: (324, 'IVF', -112)}[n]
        size = d['size']
        block = {5: 'WINGI', 8: 'HTI', 10: 'VTI', 21: 'VTI'}[n]
        links = [(buffer_word, block, buffer_word)
                 for buffer_word in range(1, top + 1)] + \
            [(size, other, 132)] + \
            [(buffer_word, other, buffer_word + shift)
             for buffer_word in range(top + 1, size)]
        back = [ln for ln in links if not (n == 10 and ln[0] == 155)]
        return links, back
    if n == 7:
        return _planform('HTI', d['iequ'], 10)
    if n == 9:
        return _planform('VTI', d['iequ'], 10)
    if n == 20:
        return _planform('VTI', d['iequ'], 172)
    if n == 16:
        links = [(buffer_word, 'FLAPIN', buffer_word)
                 for buffer_word in range(1, 117)]
        return links, links
    if n == 23:
        links = []
        for buffer_word in range(1, 101):
            common_word = (buffer_word + 1 if buffer_word <= 60
                           else buffer_word + 41)
            for off, block in ((0, 'IBODY'), (100, 'IWING'), (200, 'IHT'),
                               (300, 'IBW')):
                links.append((buffer_word + off, block, common_word))
        links += [(buffer_word + 400, 'IDWASH', buffer_word + 1)
                  for buffer_word in range(1, 61)]
        links += [(461, 'IVT', 2)]
        links += [(buffer_word + 461, 'EXPER', buffer_word + 112)
                  for buffer_word in range(1, 5)]
        links += [(466, 'SBETA', 172), (467, 'SBETA', 173),
                  (468, 'SBETA', 280), (469, 'SBETA', 281)]
        return links, links
    raise ValueError(f'no namelist {n}')


def _check_type(blocks: Blocks, block: str, word: int,
                xtype) -> List[str]:
    """XNAM7, 9, 20: a planform TYPE given as a Hollerith word becomes its
    index; a TYPE of 5 or more is reported and set to 1."""
    value = blocks[block][word]
    for type_index, name in enumerate(xtype, 1):
        if value == name:
            value = float(type_index)
    lines = []
    if not (isinstance(value, (int, float)) and value < 5.):
        lines = fortran_write(_TYPE_ERROR, [value])
        value = 1.0
    blocks[block][word] = value
    return lines


def _finish(n: int, blocks: Blocks) -> List[str]:
    """The words each namelist derives after reading or echoing."""
    b = blocks
    if n == 3:
        t = b['BODYI'][127]
        if t == UNUSED:
            t = 2.0
        if t < 1.0:
            t = 1.0
        if t > 3.0:
            t = 3.0
        b['BODYI'][127] = t
    elif n == 4:
        b['WINGD'][106] = b['WINGI'][7]
        b['WINGD'][112] = b['WINGI'][8]
        b['WINGD'][138] = 0.0
    elif n == 6:
        for k, word in enumerate((33, 65, 74, 77, 82), 1):
            b['BDATA'][word] = b['SYNTSS'][k]
    elif n in (7, 9, 20):
        source, target, off, type_word = {
            7: ('HTI', 'HTDATA', 0, 15), 9: ('VTI', 'VTDATA', 0, 15),
            20: ('VTI', 'VTDATA', 195, 177)}[n]
        b[target][off + 106] = b[source][7]
        b[target][off + 112] = b[source][8]
        b[target][off + 138] = 0.0
        return _check_type(b, source, type_word, NAMELISTS[n]['xtype'])
    elif n == 13:
        b['OPTION'][1] = b['POWER'][31]
        b['OPTION'][2] = b['POWER'][37]
    return []


def xnam(n: int, iop: int, blocks: Blocks,
         units: Optional[Dict[int, IO[str]]] = None,
         state: Optional[dict] = None) -> dict:
    """Translate XNAMn: read (``iop`` 0) or echo (``iop`` 1) namelist
    ``n`` against the COMMON ``blocks`` (edited in place).

    Args:
        n: 1 to 23.  iop: 0 to read with NAMER, 1 to print with NAMEW.
        blocks: The COMMON blocks the namelist touches.
        units: The open card files by unit number (9, or 10 for XNAM23).
        state: Saved locals: ``buffer16`` (XNAM16's word 117, never
            copied out), and NAMER's and NAMEW's.

    Returns:
        ``lines`` (records printed) and ``ieof`` (read only).
    """
    d = NAMELISTS[n]
    st = state if state is not None else {}
    buf = [0.0] * (d['size'] + 1)
    gather, back = _links(n)
    for i, block, word in gather:
        buf[i] = blocks[block][word]
    if n == 16:
        buf[117] = st.get('buffer16', 0.0)
    vname = list(''.join(d['names']))
    lenvn = [len(v) for v in d['names']]
    result = {'lines': [], 'ieof': None}
    if iop == 0:
        r = namer(KAND, units[d['unit']], list(d['name']), vname, lenvn,
                  d['dims'], buf, d['locs'],
                  st.setdefault('namer', {}))
        result.update(r)
    elif iop == 1:
        lines, st['vtype'] = namew(KAND, list(d['name']), vname, lenvn,
                                   d['dims'], buf[1:], d['locs'],
                                   st.get('vtype', 0))
        result['lines'] = lines
    for i, block, word in back:
        blocks[block][word] = buf[i]
    if n == 16:
        st['buffer16'] = buf[117]
    result['lines'] = result['lines'] + _finish(n, blocks)
    return result


def exsubt(nf: int, mdata: bool, nnames: int, rewind, read_xnam23) -> int:
    """Translate EXSUBT: reread the experimental-data namelist, once per
    name, from the top of unit 10.  Returns the reads made."""
    if nf < 0 or not mdata:
        return 0
    rewind()
    for _ in range(nnames):
        read_xnam23()
    return nnames
