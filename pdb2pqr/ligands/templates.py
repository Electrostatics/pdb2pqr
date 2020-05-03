"""Templates for various types of molecular groups.

The classification of the acids and bases is from the following publication:
    L Xing, RC Glen, RD Clark, JChInfCS, 2003, 43, 87-879

Author:  Paul Czodrowski
"""


TEMPLATES = {}

"""(I) ACIDS

    (1) Aromatic Acids; Sulfonic and Sulfinic Acids
    ===============================================
    distinction of ortho/meta/para is of paramount importance

       C.ar-C.ar        O.co2
      /         \      /
    C.ar        C.ar-C.2
      \         /      \
       C.ar-C.ar        O.co2
"""

TEMPLATES['SimpleBenzAcid']={'C.2':['O.co2','O.co2','C.ar']}
TEMPLATES['SortedSimpleBenzAcid']={'C.2':['C.ar','O.co2','O.co2']}

"""(2) Phenols and Thiophenols
    ===============================================
    only 6 thiophenols in Glen's data set => too exotic
    again, detection of ortho/meta/para substitution needed!

    (3) Aliphatic and Alicyclic Carboxylic, Sulfonic, and Sulfinic Acids
    ===============================================

                      O.co2
                     /
     -C.ar-O.3-C.3-C.2
                     \
                      O.co2
"""

TEMPLATES['SimpleAliphaticAcid']={'C.2': ['C.3','O.co2','O.co2'],
                                 'C.3': ['C.2','H','H','O.3'],
                                 'O.3': ['C.3','C.ar'],
                                 'C.ar':['C.ar','C.ar','O.3'],
                                 'root_atoms': ['C.2','C.3','O.3','C.ar']}
TEMPLATES['Acid_TypesNames']={'O1':{'sybyl_type':'O.co2','neighbours':['C3']},
                             'O2':{'sybyl_type':'O.co2','neighbours':['C3']},
                             'C3':{'sybyl_type':'C.2','neighbours':['O1','O2','C4']},
                             'C4':{'sybyl_type':'C.3','neighbours':['C3','C5','H4A','H4B']},
                             'O5':{'sybyl_type':'O.3','neighbours':['C4','C6']},
                             'C6':{'sybyl_type':'C.ar','neighbours':['O5','C6A','C6B']}}


"""(4) Aliphatic and Alicyclic Alcohols and Thiols
    ===============================================
    only alcohols with proximal strong electron withdrawing groups thiols are
    more considerably more acidic than alcohols
"""


"""(5) Acidic Nitrogens and Carbons
    ===============================================
    When strong electron withdrawing groups (e.g. nitro, nitrile, carbonyl,
    etc.) are attached to nitrogen or carbon the proton on the nitrogen or
    carbon atom may become appreciably acidic. Multiple tautomers usually
    coexist for these compounds, which can make it difficult to determine which
    class it belongs to.
"""


"""(I) BASES
    ===========

    (1a) Pyridines
    six-membered ring with one nitrogen

       C.ar-C.ar
      /         \ 
    C.ar        N.ar
      \         / 
       C.ar-C.ar 

    (1b) other six-membered rings
    pyridazine (2Ns, directly bonded)
    pyrimdine (2N, 1C in between)
    pyrazines (2N, 2C in between)

    (2) Anilines

    (3) Imidazoles

    (4) Alkylamines 
"""