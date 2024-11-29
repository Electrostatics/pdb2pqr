"""Expected results for ligand tests"""

from io import StringIO

TORSION_RESULTS = {
    "acetylcholine.mol2": {
        ("CAA", "NAD", "CAE", "CAF"),
        ("CAC", "NAD", "CAE", "CAF"),
        ("CAE", "CAF", "OAG", "CAH"),
        ("CAF", "CAE", "NAD", "CAA"),
        ("CAF", "CAE", "NAD", "CAC"),
        ("CAF", "CAE", "NAD", "CAJ"),
        ("CAF", "OAG", "CAH", "CAI"),
        ("CAF", "OAG", "CAH", "OAB"),
        ("CAH", "OAG", "CAF", "CAE"),
        ("CAI", "CAH", "OAG", "CAF"),
        ("CAJ", "NAD", "CAE", "CAF"),
        ("NAD", "CAE", "CAF", "OAG"),
        ("OAB", "CAH", "OAG", "CAF"),
        ("OAG", "CAF", "CAE", "NAD"),
    },
    "ethanol.mol2": set(),
    "pyrrole.mol2": {
        ("CAA", "CAD", "CAE", "CAC"),
        ("CAA", "NAB", "CAC", "CAE"),
        ("CAC", "CAE", "CAD", "CAA"),
        ("CAC", "NAB", "CAA", "CAD"),
        ("CAD", "CAA", "NAB", "CAC"),
        ("CAD", "CAE", "CAC", "NAB"),
        ("CAE", "CAC", "NAB", "CAA"),
        ("CAE", "CAD", "CAA", "NAB"),
        ("NAB", "CAA", "CAD", "CAE"),
        ("NAB", "CAC", "CAE", "CAD"),
    },
    "tetramethylammonium.mol2": set(),
    "glycerol.mol2": {
        ("CAA", "CAC", "CAE", "OAF"),
        ("CAE", "CAC", "CAA", "OAB"),
        ("OAB", "CAA", "CAC", "CAE"),
        ("OAB", "CAA", "CAC", "OAD"),
        ("OAD", "CAC", "CAA", "OAB"),
        ("OAD", "CAC", "CAE", "OAF"),
        ("OAF", "CAE", "CAC", "CAA"),
        ("OAF", "CAE", "CAC", "OAD"),
    },
    "cyclohexane.mol2": {
        ("CAA", "CAB", "CAC", "CAF"),
        ("CAA", "CAD", "CAE", "CAF"),
        ("CAB", "CAA", "CAD", "CAE"),
        ("CAB", "CAC", "CAF", "CAE"),
        ("CAC", "CAB", "CAA", "CAD"),
        ("CAC", "CAF", "CAE", "CAD"),
        ("CAD", "CAA", "CAB", "CAC"),
        ("CAD", "CAE", "CAF", "CAC"),
        ("CAE", "CAD", "CAA", "CAB"),
        ("CAE", "CAF", "CAC", "CAB"),
        ("CAF", "CAC", "CAB", "CAA"),
        ("CAF", "CAE", "CAD", "CAA"),
    },
}


RING_RESULTS = {
    "ethanol.mol2": set(),
    "glycerol.mol2": set(),
    "tetramethylammonium.mol2": set(),
    "acetate.mol2": set(),
    "acetonitrile.mol2": set(),
    "acetylcholine.mol2": set(),
    "fatty-acid.mol2": set(),
    "pyrrole.mol2": {("CAA", "CAD", "CAE", "CAC", "NAB")},
    "adp.mol2": {
        ("C4", "C5", "N7", "C8", "N9"),
        ("C2", "N1", "C6", "C5", "C4", "N3"),
        ("C1'", "O4'", "C4'", "C3'", "C2'"),
    },
    "cyclohexane.mol2": {("CAA", "CAD", "CAE", "CAF", "CAC", "CAB")},
    "naphthalene.mol2": {
        ("CAA", "CAF", "CAG", "CAH", "CAC", "CAB"),
        ("CAC", "CAH", "CAI", "CAJ", "CAE", "CAD"),
    },
    "anthracene.mol2": {
        ("CAA", "CAH", "CAI", "CAJ", "CAC", "CAB"),
        ("CAC", "CAJ", "CAK", "CAL", "CAE", "CAD"),
        ("CAE", "CAL", "CAM", "CAN", "CAG", "CAF"),
    },
    "1HPX-ligand.mol2": {
        ("C1", "C2", "C3", "C4", "C5", "N1"),
        ("C18", "N4", "C20", "S2", "C19"),
        ("C3", "C6", "C7", "C8", "C9", "C4"),
        ("C28", "C33", "C32", "C31", "C30", "C29"),
    },
    "1QBS-ligand.mol2": {
        ("C1", "N7", "C6", "C5", "C4", "C3", "N2"),
        ("C31", "C36", "C35", "C34", "C33", "C32"),
        ("C61", "C66", "C65", "C64", "C63", "C62"),
        ("C71", "C76", "C75", "C74", "C73", "C72"),
        ("C21", "C26", "C25", "C24", "C23", "C22"),
    },
    "1US0-ligand.mol2": {
        ("C24", "C27", "C28", "C25", "C29", "C26"),
        ("C2", "C4", "C7", "C3", "C6", "C5"),
    },
}


FORMAL_CHARGE_RESULTS = {
    "1HPX-ligand.mol2": 87 * [0],
    "1QBS-ligand.mol2": 80 * [0],
    "1US0-ligand.mol2": 22 * [0] + [-0.5, -0.5] + 11 * [0],
    "acetate.mol2": [-0.5, 0, -0.5] + 4 * [0],
    "acetylcholine.mol2": 13 * [0] + [1] + 12 * [0],
    "adp.mol2": [-1] + 5 * [0] + [-1] + 32 * [0],
    "anthracene.mol2": 24 * [0],
    "acetonitrile.mol2": 6 * [0],
    "cyclohexane.mol2": 18 * [0],
    "ethanol.mol2": 9 * [0],
    "fatty-acid.mol2": [-0.5, 0, -0.5] + 26 * [0],
    "glycerol.mol2": 14 * [0],
    "naphthalene.mol2": 18 * [0],
    "crown-ether.mol2": 42 * [0],
    "pyrrole.mol2": 10 * [0],
    "tetramethylammonium.mol2": 4 * [0] + [1] + 12 * [0],
}

PARAMETER_RESULTS = {
    "1HPX-ligand.mol2": [
        {
            "name": "C1",
            "type": "C.ar",
            "charge": -0.043081773453678196,
            "radius": 1.87,
        },
        {
            "name": "N1",
            "type": "N.ar",
            "charge": -0.35678069534508244,
            "radius": 1.4,
        },
        {
            "name": "O1",
            "type": "O.3",
            "charge": -0.5711988943947174,
            "radius": 1.52,
        },
        {
            "name": "S1",
            "type": "S.3",
            "charge": -0.36230725415046,
            "radius": 2.15,
        },
        {
            "name": "C2",
            "type": "C.ar",
            "charge": -0.07546448002333427,
            "radius": 1.87,
        },
        {
            "name": "N2",
            "type": "N.am",
            "charge": -0.703336840791423,
            "radius": 1.4,
        },
        {
            "name": "O2",
            "type": "O.3",
            "charge": -0.5909676444057169,
            "radius": 1.52,
        },
        {
            "name": "S2",
            "type": "S.3",
            "charge": -0.31326568130416904,
            "radius": 2.15,
        },
        {
            "name": "C3",
            "type": "C.ar",
            "charge": 0.012184036151037315,
            "radius": 1.87,
        },
        {
            "name": "N3",
            "type": "N.am",
            "charge": -0.7139711103762134,
            "radius": 1.4,
        },
        {
            "name": "O3",
            "type": "O.2",
            "charge": -0.2919660371248733,
            "radius": 1.52,
        },
        {
            "name": "C4",
            "type": "C.ar",
            "charge": -0.024891112570323084,
            "radius": 1.87,
        },
        {
            "name": "N4",
            "type": "N.am",
            "charge": -0.7120419698971743,
            "radius": 1.4,
        },
        {
            "name": "O4",
            "type": "O.2",
            "charge": -0.2689251910425875,
            "radius": 1.52,
        },
        {
            "name": "C5",
            "type": "C.ar",
            "charge": -0.02721768449096584,
            "radius": 1.87,
        },
        {
            "name": "N5",
            "type": "N.am",
            "charge": -0.7179288700589408,
            "radius": 1.4,
        },
        {
            "name": "O5",
            "type": "O.2",
            "charge": -0.2884010361503686,
            "radius": 1.52,
        },
        {
            "name": "C6",
            "type": "C.ar",
            "charge": 0.0753432557874654,
            "radius": 1.87,
        },
        {
            "name": "O6",
            "type": "O.2",
            "charge": -0.2912545916584642,
            "radius": 1.52,
        },
        {
            "name": "C7",
            "type": "C.ar",
            "charge": -0.06673754373256466,
            "radius": 1.87,
        },
        {
            "name": "C8",
            "type": "C.ar",
            "charge": -0.11391059914410395,
            "radius": 1.87,
        },
        {
            "name": "C9",
            "type": "C.ar",
            "charge": -0.10324637080138911,
            "radius": 1.87,
        },
        {
            "name": "C10",
            "type": "C.3",
            "charge": 0.0955779782243598,
            "radius": 1.87,
        },
        {
            "name": "C11",
            "type": "C.2",
            "charge": 0.13325937240601735,
            "radius": 1.87,
        },
        {
            "name": "C12",
            "type": "C.3",
            "charge": -0.007454479991213576,
            "radius": 1.87,
        },
        {
            "name": "C13",
            "type": "C.2",
            "charge": 0.1300074666489752,
            "radius": 1.87,
        },
        {
            "name": "C14",
            "type": "C.3",
            "charge": 0.03316079200653884,
            "radius": 1.87,
        },
        {
            "name": "C15",
            "type": "C.3",
            "charge": 0.1162788489650129,
            "radius": 1.87,
        },
        {
            "name": "C16",
            "type": "C.3",
            "charge": -0.05012048448800629,
            "radius": 1.87,
        },
        {
            "name": "C17",
            "type": "C.2",
            "charge": 0.17629988149001455,
            "radius": 1.87,
        },
        {
            "name": "C18",
            "type": "C.3",
            "charge": 0.1273124481017021,
            "radius": 1.87,
        },
        {
            "name": "C19",
            "type": "C.3",
            "charge": 0.009718364614545697,
            "radius": 1.87,
        },
        {
            "name": "C20",
            "type": "C.3",
            "charge": 0.05373224349483613,
            "radius": 1.87,
        },
        {
            "name": "C21",
            "type": "C.2",
            "charge": 0.14236781888812639,
            "radius": 1.87,
        },
        {
            "name": "C22",
            "type": "C.3",
            "charge": -0.018858356195459093,
            "radius": 1.87,
        },
        {
            "name": "C23",
            "type": "C.3",
            "charge": -0.08872387572707603,
            "radius": 1.87,
        },
        {
            "name": "C24",
            "type": "C.3",
            "charge": -0.08872387572707603,
            "radius": 1.87,
        },
        {
            "name": "C25",
            "type": "C.3",
            "charge": -0.08872387572707603,
            "radius": 1.87,
        },
        {
            "name": "C26",
            "type": "C.3",
            "charge": 0.10593153852032608,
            "radius": 1.87,
        },
        {
            "name": "C27",
            "type": "C.3",
            "charge": -0.05782750620937734,
            "radius": 1.87,
        },
        {
            "name": "C28",
            "type": "C.ar",
            "charge": -0.09940049289330599,
            "radius": 1.87,
        },
        {
            "name": "C29",
            "type": "C.ar",
            "charge": -0.11788247675515531,
            "radius": 1.87,
        },
        {
            "name": "C30",
            "type": "C.ar",
            "charge": -0.12669321002628098,
            "radius": 1.87,
        },
        {
            "name": "C31",
            "type": "C.ar",
            "charge": -0.12840215383108972,
            "radius": 1.87,
        },
        {
            "name": "C32",
            "type": "C.ar",
            "charge": -0.12669321002628098,
            "radius": 1.87,
        },
        {
            "name": "C33",
            "type": "C.ar",
            "charge": -0.11788247675515531,
            "radius": 1.87,
        },
        {
            "name": "HN2",
            "type": "H",
            "charge": 0.47365062709293504,
            "radius": 1.1,
        },
        {
            "name": "HN3",
            "type": "H",
            "charge": 0.46813799183459115,
            "radius": 1.1,
        },
        {
            "name": "HN5",
            "type": "H",
            "charge": 0.4699450251195792,
            "radius": 1.1,
        },
        {
            "name": "HO2",
            "type": "H",
            "charge": 0.4135083335659837,
            "radius": 1.1,
        },
        {
            "name": "H1",
            "type": "H",
            "charge": 0.1930651515377243,
            "radius": 1.1,
        },
        {
            "name": "H2",
            "type": "H",
            "charge": 0.14930498376978826,
            "radius": 1.1,
        },
        {
            "name": "H5",
            "type": "H",
            "charge": 0.2004887586723136,
            "radius": 1.1,
        },
        {
            "name": "H7",
            "type": "H",
            "charge": 0.1568207646303434,
            "radius": 1.1,
        },
        {
            "name": "H8",
            "type": "H",
            "charge": 0.13234810374569977,
            "radius": 1.1,
        },
        {
            "name": "H9",
            "type": "H",
            "charge": 0.13832988708343302,
            "radius": 1.1,
        },
        {
            "name": "1H10",
            "type": "H",
            "charge": 0.2090610332096663,
            "radius": 1.1,
        },
        {
            "name": "2H10",
            "type": "H",
            "charge": 0.2090610332096663,
            "radius": 1.1,
        },
        {
            "name": "1H12",
            "type": "H",
            "charge": 0.12111462041395545,
            "radius": 1.1,
        },
        {
            "name": "2H12",
            "type": "H",
            "charge": 0.12111462041395545,
            "radius": 1.1,
        },
        {
            "name": "H14",
            "type": "H",
            "charge": 0.1675599590015787,
            "radius": 1.1,
        },
        {
            "name": "H15",
            "type": "H",
            "charge": 0.22760838606623204,
            "radius": 1.1,
        },
        {
            "name": "1H16",
            "type": "H",
            "charge": 0.0852867334715187,
            "radius": 1.1,
        },
        {
            "name": "2H16",
            "type": "H",
            "charge": 0.0852867334715187,
            "radius": 1.1,
        },
        {
            "name": "H18",
            "type": "H",
            "charge": 0.24317031751195953,
            "radius": 1.1,
        },
        {
            "name": "1H19",
            "type": "H",
            "charge": 0.12776560931027464,
            "radius": 1.1,
        },
        {
            "name": "2H19",
            "type": "H",
            "charge": 0.12776560931027464,
            "radius": 1.1,
        },
        {
            "name": "1H20",
            "type": "H",
            "charge": 0.19565330098691133,
            "radius": 1.1,
        },
        {
            "name": "2H20",
            "type": "H",
            "charge": 0.19565330098691133,
            "radius": 1.1,
        },
        {
            "name": "1H23",
            "type": "H",
            "charge": 0.05671902293773316,
            "radius": 1.1,
        },
        {
            "name": "2H23",
            "type": "H",
            "charge": 0.05671902293773316,
            "radius": 1.1,
        },
        {
            "name": "3H23",
            "type": "H",
            "charge": 0.05671902293773316,
            "radius": 1.1,
        },
        {
            "name": "1H24",
            "type": "H",
            "charge": 0.05671902293773316,
            "radius": 1.1,
        },
        {
            "name": "2H24",
            "type": "H",
            "charge": 0.05671902293773316,
            "radius": 1.1,
        },
        {
            "name": "3H24",
            "type": "H",
            "charge": 0.05671902293773316,
            "radius": 1.1,
        },
        {
            "name": "1H25",
            "type": "H",
            "charge": 0.05671902293773316,
            "radius": 1.1,
        },
        {
            "name": "2H25",
            "type": "H",
            "charge": 0.05671902293773316,
            "radius": 1.1,
        },
        {
            "name": "3H25",
            "type": "H",
            "charge": 0.05671902293773316,
            "radius": 1.1,
        },
        {
            "name": "H26",
            "type": "H",
            "charge": 0.21276929701184144,
            "radius": 1.1,
        },
        {
            "name": "1H27",
            "type": "H",
            "charge": 0.0838296767957033,
            "radius": 1.1,
        },
        {
            "name": "2H27",
            "type": "H",
            "charge": 0.0838296767957033,
            "radius": 1.1,
        },
        {
            "name": "3H27",
            "type": "H",
            "charge": 0.0838296767957033,
            "radius": 1.1,
        },
        {
            "name": "H29",
            "type": "H",
            "charge": 0.1344214537548114,
            "radius": 1.1,
        },
        {
            "name": "H30",
            "type": "H",
            "charge": 0.12944135663694809,
            "radius": 1.1,
        },
        {
            "name": "H31",
            "type": "H",
            "charge": 0.12895177093126134,
            "radius": 1.1,
        },
        {
            "name": "H32",
            "type": "H",
            "charge": 0.12944135663694809,
            "radius": 1.1,
        },
        {
            "name": "H33",
            "type": "H",
            "charge": 0.1344214537548114,
            "radius": 1.1,
        },
    ],
    "1QBS-ligand.mol2": [
        {
            "name": "C1",
            "type": "C.2",
            "charge": 0.2819290256115998,
            "radius": 1.87,
        },
        {
            "name": "C3",
            "type": "C.3",
            "charge": 0.044938629062860544,
            "radius": 1.87,
        },
        {
            "name": "C4",
            "type": "C.3",
            "charge": 0.029905546428555488,
            "radius": 1.87,
        },
        {
            "name": "C5",
            "type": "C.3",
            "charge": 0.029905546428555488,
            "radius": 1.87,
        },
        {
            "name": "C6",
            "type": "C.3",
            "charge": 0.044938629062860544,
            "radius": 1.87,
        },
        {
            "name": "C20",
            "type": "C.3",
            "charge": 0.02388740859601932,
            "radius": 1.87,
        },
        {
            "name": "C21",
            "type": "C.ar",
            "charge": -0.06923087776665914,
            "radius": 1.87,
        },
        {
            "name": "C22",
            "type": "C.ar",
            "charge": -0.10856468032765192,
            "radius": 1.87,
        },
        {
            "name": "C23",
            "type": "C.ar",
            "charge": -0.11036338907520785,
            "radius": 1.87,
        },
        {
            "name": "C24",
            "type": "C.ar",
            "charge": -0.08244776383718133,
            "radius": 1.87,
        },
        {
            "name": "C25",
            "type": "C.ar",
            "charge": -0.11036338907520785,
            "radius": 1.87,
        },
        {
            "name": "C26",
            "type": "C.ar",
            "charge": -0.10856468032765192,
            "radius": 1.87,
        },
        {
            "name": "C27",
            "type": "C.3",
            "charge": -0.021297911008685094,
            "radius": 1.87,
        },
        {
            "name": "C30",
            "type": "C.3",
            "charge": -0.04323885537049762,
            "radius": 1.87,
        },
        {
            "name": "C31",
            "type": "C.ar",
            "charge": -0.09772541979676261,
            "radius": 1.87,
        },
        {
            "name": "C32",
            "type": "C.ar",
            "charge": -0.11774626937995367,
            "radius": 1.87,
        },
        {
            "name": "C33",
            "type": "C.ar",
            "charge": -0.12668933876013122,
            "radius": 1.87,
        },
        {
            "name": "C34",
            "type": "C.ar",
            "charge": -0.12840215383108972,
            "radius": 1.87,
        },
        {
            "name": "C35",
            "type": "C.ar",
            "charge": -0.12668933876013122,
            "radius": 1.87,
        },
        {
            "name": "C36",
            "type": "C.ar",
            "charge": -0.11774626937995367,
            "radius": 1.87,
        },
        {
            "name": "C60",
            "type": "C.3",
            "charge": -0.04323885537049762,
            "radius": 1.87,
        },
        {
            "name": "C61",
            "type": "C.ar",
            "charge": -0.09772541979676261,
            "radius": 1.87,
        },
        {
            "name": "C62",
            "type": "C.ar",
            "charge": -0.11774626937995367,
            "radius": 1.87,
        },
        {
            "name": "C63",
            "type": "C.ar",
            "charge": -0.12668933876013122,
            "radius": 1.87,
        },
        {
            "name": "C64",
            "type": "C.ar",
            "charge": -0.12840215383108972,
            "radius": 1.87,
        },
        {
            "name": "C65",
            "type": "C.ar",
            "charge": -0.12668933876013122,
            "radius": 1.87,
        },
        {
            "name": "C66",
            "type": "C.ar",
            "charge": -0.11774626937995367,
            "radius": 1.87,
        },
        {
            "name": "C70",
            "type": "C.3",
            "charge": 0.02388740859601932,
            "radius": 1.87,
        },
        {
            "name": "C71",
            "type": "C.ar",
            "charge": -0.06923087776665914,
            "radius": 1.87,
        },
        {
            "name": "C72",
            "type": "C.ar",
            "charge": -0.10856468032765192,
            "radius": 1.87,
        },
        {
            "name": "C73",
            "type": "C.ar",
            "charge": -0.11036338907520785,
            "radius": 1.87,
        },
        {
            "name": "C74",
            "type": "C.ar",
            "charge": -0.08244776383718133,
            "radius": 1.87,
        },
        {
            "name": "C75",
            "type": "C.ar",
            "charge": -0.11036338907520785,
            "radius": 1.87,
        },
        {
            "name": "C76",
            "type": "C.ar",
            "charge": -0.10856468032765192,
            "radius": 1.87,
        },
        {
            "name": "C77",
            "type": "C.3",
            "charge": -0.021297911008685094,
            "radius": 1.87,
        },
        {
            "name": "N2",
            "type": "N.am",
            "charge": -0.7121260761527837,
            "radius": 1.4,
        },
        {
            "name": "N7",
            "type": "N.am",
            "charge": -0.7121260761527837,
            "radius": 1.4,
        },
        {
            "name": "O1",
            "type": "O.2",
            "charge": -0.19835509547371977,
            "radius": 1.52,
        },
        {
            "name": "O4",
            "type": "O.3",
            "charge": -0.6060823059727466,
            "radius": 1.52,
        },
        {
            "name": "O5",
            "type": "O.3",
            "charge": -0.6060823059727466,
            "radius": 1.52,
        },
        {
            "name": "O27",
            "type": "O.3",
            "charge": -0.6190964861641765,
            "radius": 1.52,
        },
        {
            "name": "O77",
            "type": "O.3",
            "charge": -0.6190964861641765,
            "radius": 1.52,
        },
        {
            "name": "H1",
            "type": "H",
            "charge": 0.1939484454048374,
            "radius": 1.1,
        },
        {
            "name": "H2",
            "type": "H",
            "charge": 0.1783837779801701,
            "radius": 1.1,
        },
        {
            "name": "H3",
            "type": "H",
            "charge": 0.1783837779801701,
            "radius": 1.1,
        },
        {
            "name": "H4",
            "type": "H",
            "charge": 0.1939484454048374,
            "radius": 1.1,
        },
        {
            "name": "H5",
            "type": "H",
            "charge": 0.1660607436631621,
            "radius": 1.1,
        },
        {
            "name": "H6",
            "type": "H",
            "charge": 0.1660607436631621,
            "radius": 1.1,
        },
        {
            "name": "H7",
            "type": "H",
            "charge": 0.13696512588510704,
            "radius": 1.1,
        },
        {
            "name": "H8",
            "type": "H",
            "charge": 0.13656070949143695,
            "radius": 1.1,
        },
        {
            "name": "H9",
            "type": "H",
            "charge": 0.13656070949143695,
            "radius": 1.1,
        },
        {
            "name": "H10",
            "type": "H",
            "charge": 0.13696512588510704,
            "radius": 1.1,
        },
        {
            "name": "H11",
            "type": "H",
            "charge": 0.1336796298532153,
            "radius": 1.1,
        },
        {
            "name": "H12",
            "type": "H",
            "charge": 0.1336796298532153,
            "radius": 1.1,
        },
        {
            "name": "H13",
            "type": "H",
            "charge": 0.08892832134853815,
            "radius": 1.1,
        },
        {
            "name": "H14",
            "type": "H",
            "charge": 0.08892832134853815,
            "radius": 1.1,
        },
        {
            "name": "H15",
            "type": "H",
            "charge": 0.13442892363549896,
            "radius": 1.1,
        },
        {
            "name": "H16",
            "type": "H",
            "charge": 0.12944135663694809,
            "radius": 1.1,
        },
        {
            "name": "H17",
            "type": "H",
            "charge": 0.12895177093126134,
            "radius": 1.1,
        },
        {
            "name": "H18",
            "type": "H",
            "charge": 0.12944135663694809,
            "radius": 1.1,
        },
        {
            "name": "H19",
            "type": "H",
            "charge": 0.13442892363549896,
            "radius": 1.1,
        },
        {
            "name": "H20",
            "type": "H",
            "charge": 0.08892832134853815,
            "radius": 1.1,
        },
        {
            "name": "H21",
            "type": "H",
            "charge": 0.08892832134853815,
            "radius": 1.1,
        },
        {
            "name": "H22",
            "type": "H",
            "charge": 0.13442892363549896,
            "radius": 1.1,
        },
        {
            "name": "H23",
            "type": "H",
            "charge": 0.12944135663694809,
            "radius": 1.1,
        },
        {
            "name": "H24",
            "type": "H",
            "charge": 0.12895177093126134,
            "radius": 1.1,
        },
        {
            "name": "H25",
            "type": "H",
            "charge": 0.12944135663694809,
            "radius": 1.1,
        },
        {
            "name": "H26",
            "type": "H",
            "charge": 0.13442892363549896,
            "radius": 1.1,
        },
        {
            "name": "H27",
            "type": "H",
            "charge": 0.1660607436631621,
            "radius": 1.1,
        },
        {
            "name": "H28",
            "type": "H",
            "charge": 0.1660607436631621,
            "radius": 1.1,
        },
        {
            "name": "H29",
            "type": "H",
            "charge": 0.13696512588510704,
            "radius": 1.1,
        },
        {
            "name": "H30",
            "type": "H",
            "charge": 0.13656070949143695,
            "radius": 1.1,
        },
        {
            "name": "H31",
            "type": "H",
            "charge": 0.13656070949143695,
            "radius": 1.1,
        },
        {
            "name": "H32",
            "type": "H",
            "charge": 0.13696512588510704,
            "radius": 1.1,
        },
        {
            "name": "H33",
            "type": "H",
            "charge": 0.1336796298532153,
            "radius": 1.1,
        },
        {
            "name": "H34",
            "type": "H",
            "charge": 0.1336796298532153,
            "radius": 1.1,
        },
        {
            "name": "H35",
            "type": "H",
            "charge": 0.4089032076940611,
            "radius": 1.1,
        },
        {
            "name": "H36",
            "type": "H",
            "charge": 0.4089032076940611,
            "radius": 1.1,
        },
        {
            "name": "H37",
            "type": "H",
            "charge": 0.40353983279195227,
            "radius": 1.1,
        },
        {
            "name": "H38",
            "type": "H",
            "charge": 0.40353983279195227,
            "radius": 1.1,
        },
    ],
    "1US0-ligand.mol2": [
        {
            "name": "C2",
            "type": "C.ar",
            "charge": -0.011692507946738866,
            "radius": 1.87,
        },
        {
            "name": "C3",
            "type": "C.ar",
            "charge": -0.07564532850595106,
            "radius": 1.87,
        },
        {
            "name": "C4",
            "type": "C.ar",
            "charge": 0.0995106215449818,
            "radius": 1.87,
        },
        {
            "name": "C5",
            "type": "C.ar",
            "charge": 0.087912610051304,
            "radius": 1.87,
        },
        {
            "name": "C6",
            "type": "C.ar",
            "charge": -0.06822502752899269,
            "radius": 1.87,
        },
        {
            "name": "C7",
            "type": "C.ar",
            "charge": 0.071206877934909,
            "radius": 1.87,
        },
        {
            "name": "BR8",
            "type": "Br",
            "charge": -0.1781415724110881,
            "radius": 1.85,
        },
        {
            "name": "F9",
            "type": "F",
            "charge": -0.24609828999108224,
            "radius": 2.4,
        },
        {
            "name": "C11",
            "type": "C.2",
            "charge": 0.01969593460035858,
            "radius": 1.87,
        },
        {
            "name": "C13",
            "type": "C.3",
            "charge": -0.012644223822613768,
            "radius": 1.87,
        },
        {
            "name": "F14",
            "type": "F",
            "charge": -0.24298404863346873,
            "radius": 2.4,
        },
        {
            "name": "O15",
            "type": "O.3",
            "charge": -0.564853763228142,
            "radius": 1.52,
        },
        {
            "name": "S16",
            "type": "S.3",
            "charge": -0.08677321834382315,
            "radius": 2.15,
        },
        {
            "name": "N17",
            "type": "N.pl3",
            "charge": -0.5406144014197267,
            "radius": 1.4,
        },
        {
            "name": "C20",
            "type": "C.3",
            "charge": 0.11379374529001238,
            "radius": 1.87,
        },
        {
            "name": "C24",
            "type": "C.ar",
            "charge": -0.025165867794008426,
            "radius": 1.87,
        },
        {
            "name": "C25",
            "type": "C.ar",
            "charge": 0.05247378966494793,
            "radius": 1.87,
        },
        {
            "name": "C26",
            "type": "C.ar",
            "charge": -0.09820934352177407,
            "radius": 1.87,
        },
        {
            "name": "C27",
            "type": "C.ar",
            "charge": 0.09943434454444673,
            "radius": 1.87,
        },
        {
            "name": "C28",
            "type": "C.ar",
            "charge": -0.030058244343191105,
            "radius": 1.87,
        },
        {
            "name": "C29",
            "type": "C.ar",
            "charge": -0.08190135488592247,
            "radius": 1.87,
        },
        {
            "name": "C32",
            "type": "C.2",
            "charge": 0.13523301177634597,
            "radius": 1.87,
        },
        {
            "name": "O33",
            "type": "O.co2",
            "charge": -0.6823845278785808,
            "radius": 1.76,
        },
        {
            "name": "O34",
            "type": "O.co2",
            "charge": -0.6823845278785808,
            "radius": 1.76,
        },
        {
            "name": "H1",
            "type": "H",
            "charge": 0.1785741592725812,
            "radius": 1.1,
        },
        {
            "name": "H2",
            "type": "H",
            "charge": 0.1460165676666413,
            "radius": 1.1,
        },
        {
            "name": "H3",
            "type": "H",
            "charge": 0.15178137252163487,
            "radius": 1.1,
        },
        {
            "name": "H4",
            "type": "H",
            "charge": 0.12305674443192938,
            "radius": 1.1,
        },
        {
            "name": "H5",
            "type": "H",
            "charge": 0.12305674443192938,
            "radius": 1.1,
        },
        {
            "name": "H6",
            "type": "H",
            "charge": 0.3361191041929824,
            "radius": 1.1,
        },
        {
            "name": "H7",
            "type": "H",
            "charge": 0.21824692354069947,
            "radius": 1.1,
        },
        {
            "name": "H8",
            "type": "H",
            "charge": 0.21824692354069947,
            "radius": 1.1,
        },
        {
            "name": "H9",
            "type": "H",
            "charge": 0.13902103205782942,
            "radius": 1.1,
        },
        {
            "name": "H10",
            "type": "H",
            "charge": 0.16785538090669733,
            "radius": 1.1,
        },
        {
            "name": "H11",
            "type": "H",
            "charge": 0.14654036016275432,
            "radius": 1.1,
        },
    ],
    "acetate.mol2": [
        {
            "name": "OAC",
            "type": "O.co2",
            "charge": -0.6920003959699945,
            "radius": 1.76,
        },
        {
            "name": "CAB",
            "type": "C.2",
            "charge": 0.08786702981704235,
            "radius": 1.87,
        },
        {
            "name": "OAD",
            "type": "O.co2",
            "charge": -0.6920003959699945,
            "radius": 1.76,
        },
        {
            "name": "CAA",
            "type": "C.3",
            "charge": -0.007013553584021696,
            "radius": 1.87,
        },
        {
            "name": "HAB",
            "type": "H",
            "charge": 0.10104910523565629,
            "radius": 1.1,
        },
        {
            "name": "HAC",
            "type": "H",
            "charge": 0.10104910523565629,
            "radius": 1.1,
        },
        {
            "name": "HAA",
            "type": "H",
            "charge": 0.10104910523565629,
            "radius": 1.1,
        },
    ],
    "acetonitrile.mol2": [
        {
            "name": "NAC",
            "type": "N.1",
            "charge": -0.3658493397941428,
            "radius": 1.4,
        },
        {
            "name": "CAB",
            "type": "C.1",
            "charge": 0.047166532439326746,
            "radius": 1.87,
        },
        {
            "name": "CAA",
            "type": "C.3",
            "charge": 0.0031209482290866543,
            "radius": 1.87,
        },
        {
            "name": "HAB",
            "type": "H",
            "charge": 0.10518728637524313,
            "radius": 1.1,
        },
        {
            "name": "HAC",
            "type": "H",
            "charge": 0.10518728637524313,
            "radius": 1.1,
        },
        {
            "name": "HAA",
            "type": "H",
            "charge": 0.10518728637524313,
            "radius": 1.1,
        },
    ],
    "acetylcholine.mol2": [
        {
            "name": "CAI",
            "type": "C.3",
            "charge": -0.008419213752965279,
            "radius": 1.87,
        },
        {
            "name": "HAK",
            "type": "H",
            "charge": 0.09728030818644891,
            "radius": 1.1,
        },
        {
            "name": "HAL",
            "type": "H",
            "charge": 0.09728030818644891,
            "radius": 1.1,
        },
        {
            "name": "HAM",
            "type": "H",
            "charge": 0.09728030818644891,
            "radius": 1.1,
        },
        {
            "name": "CAH",
            "type": "C.2",
            "charge": 0.1545682666391745,
            "radius": 1.87,
        },
        {
            "name": "OAB",
            "type": "O.2",
            "charge": -0.27967627950752316,
            "radius": 1.52,
        },
        {
            "name": "OAG",
            "type": "O.3",
            "charge": -0.5429502778137623,
            "radius": 1.52,
        },
        {
            "name": "CAF",
            "type": "C.3",
            "charge": 0.037889841506971396,
            "radius": 1.87,
        },
        {
            "name": "HAI",
            "type": "H",
            "charge": 0.1710371780188628,
            "radius": 1.1,
        },
        {
            "name": "HAJ",
            "type": "H",
            "charge": 0.1710371780188628,
            "radius": 1.1,
        },
        {
            "name": "CAE",
            "type": "C.3",
            "charge": 0.0572339093088662,
            "radius": 1.87,
        },
        {
            "name": "HAG",
            "type": "H",
            "charge": 0.16582144135089144,
            "radius": 1.1,
        },
        {
            "name": "HAH",
            "type": "H",
            "charge": 0.16582144135089144,
            "radius": 1.1,
        },
        {
            "name": "NAD",
            "type": "N.3",
            "charge": -0.6101582744216567,
            "radius": 1.4,
        },
        {
            "name": "CAJ",
            "type": "C.3",
            "charge": 0.016251330821632083,
            "radius": 1.87,
        },
        {
            "name": "HAO",
            "type": "H",
            "charge": 0.13079998580857152,
            "radius": 1.1,
        },
        {
            "name": "HAP",
            "type": "H",
            "charge": 0.13079998580857152,
            "radius": 1.1,
        },
        {
            "name": "HAN",
            "type": "H",
            "charge": 0.13079998580857152,
            "radius": 1.1,
        },
        {
            "name": "CAA",
            "type": "C.3",
            "charge": 0.016251330821632083,
            "radius": 1.87,
        },
        {
            "name": "HAB",
            "type": "H",
            "charge": 0.13079998580857152,
            "radius": 1.1,
        },
        {
            "name": "HAC",
            "type": "H",
            "charge": 0.13079998580857152,
            "radius": 1.1,
        },
        {
            "name": "HAA",
            "type": "H",
            "charge": 0.13079998580857152,
            "radius": 1.1,
        },
        {
            "name": "CAC",
            "type": "C.3",
            "charge": 0.016251330821632083,
            "radius": 1.87,
        },
        {
            "name": "HAE",
            "type": "H",
            "charge": 0.13079998580857152,
            "radius": 1.1,
        },
        {
            "name": "HAF",
            "type": "H",
            "charge": 0.13079998580857152,
            "radius": 1.1,
        },
        {
            "name": "HAD",
            "type": "H",
            "charge": 0.13079998580857152,
            "radius": 1.1,
        },
    ],
    "adp.mol2": [
        {
            "name": "O2B",
            "type": "O.3",
            "charge": -0.9292599002864984,
            "radius": 1.52,
        },
        {
            "name": "PB",
            "type": "P.3",
            "charge": 0.2901422925812394,
            "radius": 1.8,
        },
        {
            "name": "O3B",
            "type": "O.3",
            "charge": -0.15563495903544203,
            "radius": 1.52,
        },
        {
            "name": "O1B",
            "type": "O.2",
            "charge": -0.10031235497346548,
            "radius": 1.52,
        },
        {
            "name": "O3A",
            "type": "O.3",
            "charge": -0.22961271227927002,
            "radius": 1.52,
        },
        {
            "name": "PA",
            "type": "P.3",
            "charge": 0.20924046046691175,
            "radius": 1.8,
        },
        {
            "name": "O2A",
            "type": "O.3",
            "charge": -0.9415448332469695,
            "radius": 1.52,
        },
        {
            "name": "O1A",
            "type": "O.2",
            "charge": -0.11862974402094492,
            "radius": 1.52,
        },
        {
            "name": "O5'",
            "type": "O.3",
            "charge": -0.46621302263040854,
            "radius": 1.52,
        },
        {
            "name": "C5'",
            "type": "C.3",
            "charge": 0.0798878056064067,
            "radius": 1.87,
        },
        {
            "name": "H5'",
            "type": "H",
            "charge": 0.18839134309240957,
            "radius": 1.1,
        },
        {
            "name": "H5S",
            "type": "H",
            "charge": 0.18839134309240957,
            "radius": 1.1,
        },
        {
            "name": "C4'",
            "type": "C.3",
            "charge": 0.0695445224862847,
            "radius": 1.87,
        },
        {
            "name": "H4'",
            "type": "H",
            "charge": 0.20776211836320221,
            "radius": 1.1,
        },
        {
            "name": "O4'",
            "type": "O.3",
            "charge": -0.5870759018169631,
            "radius": 1.52,
        },
        {
            "name": "C3'",
            "type": "C.3",
            "charge": 0.03852622354001534,
            "radius": 1.87,
        },
        {
            "name": "H3'",
            "type": "H",
            "charge": 0.17992130819444715,
            "radius": 1.1,
        },
        {
            "name": "O3'",
            "type": "O.3",
            "charge": -0.6050957763635193,
            "radius": 1.52,
        },
        {
            "name": "H8L",
            "type": "H",
            "charge": 0.408917570005205,
            "radius": 1.1,
        },
        {
            "name": "C2'",
            "type": "C.3",
            "charge": 0.04785444900747468,
            "radius": 1.87,
        },
        {
            "name": "H2'",
            "type": "H",
            "charge": 0.18698431697045484,
            "radius": 1.1,
        },
        {
            "name": "O2'",
            "type": "O.3",
            "charge": -0.6030269458228865,
            "radius": 1.52,
        },
        {
            "name": "H8M",
            "type": "H",
            "charge": 0.40956200683924104,
            "radius": 1.1,
        },
        {
            "name": "C1'",
            "type": "C.3",
            "charge": 0.1072884172431409,
            "radius": 1.87,
        },
        {
            "name": "H1'",
            "type": "H",
            "charge": 0.25190824252954014,
            "radius": 1.1,
        },
        {
            "name": "N9",
            "type": "N.ar",
            "charge": -0.32267909962343944,
            "radius": 1.4,
        },
        {
            "name": "C8",
            "type": "C.ar",
            "charge": 0.03749022549830809,
            "radius": 1.87,
        },
        {
            "name": "H8",
            "type": "H",
            "charge": 0.24861896888732188,
            "radius": 1.1,
        },
        {
            "name": "N7",
            "type": "N.ar",
            "charge": -0.2902489369083555,
            "radius": 1.4,
        },
        {
            "name": "C5",
            "type": "C.ar",
            "charge": 0.15701292932781932,
            "radius": 1.87,
        },
        {
            "name": "C4",
            "type": "C.ar",
            "charge": 0.16341092847650182,
            "radius": 1.87,
        },
        {
            "name": "N3",
            "type": "N.ar",
            "charge": -0.28837259982632546,
            "radius": 1.4,
        },
        {
            "name": "C2",
            "type": "C.ar",
            "charge": 0.00966560470099504,
            "radius": 1.87,
        },
        {
            "name": "H2",
            "type": "H",
            "charge": 0.23340110948012485,
            "radius": 1.1,
        },
        {
            "name": "N1",
            "type": "N.ar",
            "charge": -0.3751603212355695,
            "radius": 1.4,
        },
        {
            "name": "H1",
            "type": "H",
            "charge": 0.2881016272229071,
            "radius": 1.1,
        },
        {
            "name": "C6",
            "type": "C.ar",
            "charge": 0.10850097968415028,
            "radius": 1.87,
        },
        {
            "name": "N6",
            "type": "N.2",
            "charge": -0.46957750571325696,
            "radius": 1.4,
        },
        {
            "name": "H6",
            "type": "H",
            "charge": 0.3719198204868037,
            "radius": 1.1,
        },
    ],
    "anthracene.mol2": [
        {
            "name": "CAB",
            "type": "C.ar",
            "charge": -0.10948195185220932,
            "radius": 1.87,
        },
        {
            "name": "HAB",
            "type": "H",
            "charge": 0.13718781330297145,
            "radius": 1.1,
        },
        {
            "name": "CAC",
            "type": "C.ar",
            "charge": -0.05951797220652635,
            "radius": 1.87,
        },
        {
            "name": "CAD",
            "type": "C.ar",
            "charge": -0.090739817012503,
            "radius": 1.87,
        },
        {
            "name": "HAD",
            "type": "H",
            "charge": 0.14537539214857495,
            "radius": 1.1,
        },
        {
            "name": "CAE",
            "type": "C.ar",
            "charge": -0.05951797220652635,
            "radius": 1.87,
        },
        {
            "name": "CAF",
            "type": "C.ar",
            "charge": -0.10948195185220932,
            "radius": 1.87,
        },
        {
            "name": "HAF",
            "type": "H",
            "charge": 0.13718781330297145,
            "radius": 1.1,
        },
        {
            "name": "CAG",
            "type": "C.ar",
            "charge": -0.1252004685570065,
            "radius": 1.87,
        },
        {
            "name": "HAG",
            "type": "H",
            "charge": 0.12969479174473472,
            "radius": 1.1,
        },
        {
            "name": "CAN",
            "type": "C.ar",
            "charge": -0.1252004685570065,
            "radius": 1.87,
        },
        {
            "name": "HAN",
            "type": "H",
            "charge": 0.12969479174473472,
            "radius": 1.1,
        },
        {
            "name": "CAM",
            "type": "C.ar",
            "charge": -0.10948195185220932,
            "radius": 1.87,
        },
        {
            "name": "HAM",
            "type": "H",
            "charge": 0.13718781330297145,
            "radius": 1.1,
        },
        {
            "name": "CAL",
            "type": "C.ar",
            "charge": -0.05951797220652635,
            "radius": 1.87,
        },
        {
            "name": "CAK",
            "type": "C.ar",
            "charge": -0.090739817012503,
            "radius": 1.87,
        },
        {
            "name": "HAK",
            "type": "H",
            "charge": 0.14537539214857495,
            "radius": 1.1,
        },
        {
            "name": "CAJ",
            "type": "C.ar",
            "charge": -0.05951797220652635,
            "radius": 1.87,
        },
        {
            "name": "CAI",
            "type": "C.ar",
            "charge": -0.10948195185220932,
            "radius": 1.87,
        },
        {
            "name": "HAI",
            "type": "H",
            "charge": 0.13718781330297145,
            "radius": 1.1,
        },
        {
            "name": "CAH",
            "type": "C.ar",
            "charge": -0.1252004685570065,
            "radius": 1.87,
        },
        {
            "name": "HAH",
            "type": "H",
            "charge": 0.12969479174473472,
            "radius": 1.1,
        },
        {
            "name": "CAA",
            "type": "C.ar",
            "charge": -0.1252004685570065,
            "radius": 1.87,
        },
        {
            "name": "HAA",
            "type": "H",
            "charge": 0.12969479174473472,
            "radius": 1.1,
        },
    ],
    "cyclohexane.mol2": [
        {
            "name": "CAA",
            "type": "C.3",
            "charge": -0.10295775116564801,
            "radius": 1.87,
        },
        {
            "name": "HAA",
            "type": "H",
            "charge": 0.051478875582824005,
            "radius": 1.1,
        },
        {
            "name": "HAB",
            "type": "H",
            "charge": 0.051478875582824005,
            "radius": 1.1,
        },
        {
            "name": "CAB",
            "type": "C.3",
            "charge": -0.10295775116564801,
            "radius": 1.87,
        },
        {
            "name": "HAC",
            "type": "H",
            "charge": 0.051478875582824005,
            "radius": 1.1,
        },
        {
            "name": "HAD",
            "type": "H",
            "charge": 0.051478875582824005,
            "radius": 1.1,
        },
        {
            "name": "CAC",
            "type": "C.3",
            "charge": -0.10295775116564801,
            "radius": 1.87,
        },
        {
            "name": "HAE",
            "type": "H",
            "charge": 0.051478875582824005,
            "radius": 1.1,
        },
        {
            "name": "HAF",
            "type": "H",
            "charge": 0.051478875582824005,
            "radius": 1.1,
        },
        {
            "name": "CAF",
            "type": "C.3",
            "charge": -0.10295775116564801,
            "radius": 1.87,
        },
        {
            "name": "HAK",
            "type": "H",
            "charge": 0.051478875582824005,
            "radius": 1.1,
        },
        {
            "name": "HAL",
            "type": "H",
            "charge": 0.051478875582824005,
            "radius": 1.1,
        },
        {
            "name": "CAE",
            "type": "C.3",
            "charge": -0.10295775116564801,
            "radius": 1.87,
        },
        {
            "name": "HAI",
            "type": "H",
            "charge": 0.051478875582824005,
            "radius": 1.1,
        },
        {
            "name": "HAJ",
            "type": "H",
            "charge": 0.051478875582824005,
            "radius": 1.1,
        },
        {
            "name": "CAD",
            "type": "C.3",
            "charge": -0.10295775116564801,
            "radius": 1.87,
        },
        {
            "name": "HAH",
            "type": "H",
            "charge": 0.051478875582824005,
            "radius": 1.1,
        },
        {
            "name": "HAG",
            "type": "H",
            "charge": 0.051478875582824005,
            "radius": 1.1,
        },
    ],
    "ethanol.mol2": [
        {
            "name": "CAA",
            "type": "C.3",
            "charge": -0.09830853015557071,
            "radius": 1.87,
        },
        {
            "name": "HAA",
            "type": "H",
            "charge": 0.04965533009730527,
            "radius": 1.1,
        },
        {
            "name": "HAB",
            "type": "H",
            "charge": 0.04965533009730527,
            "radius": 1.1,
        },
        {
            "name": "HAC",
            "type": "H",
            "charge": 0.04965533009730527,
            "radius": 1.1,
        },
        {
            "name": "CAB",
            "type": "C.3",
            "charge": -0.05333201263132795,
            "radius": 1.87,
        },
        {
            "name": "HAD",
            "type": "H",
            "charge": 0.1136006498229285,
            "radius": 1.1,
        },
        {
            "name": "HAE",
            "type": "H",
            "charge": 0.1136006498229285,
            "radius": 1.1,
        },
        {
            "name": "OAC",
            "type": "O.3",
            "charge": -0.6256712436819194,
            "radius": 1.52,
        },
        {
            "name": "HAF",
            "type": "H",
            "charge": 0.4011444965310452,
            "radius": 1.1,
        },
    ],
    "fatty-acid.mol2": [
        {
            "name": "OAA",
            "type": "O.co2",
            "charge": -0.6905415761361866,
            "radius": 1.76,
        },
        {
            "name": "CAB",
            "type": "C.2",
            "charge": 0.09491530638614049,
            "radius": 1.87,
        },
        {
            "name": "OAL",
            "type": "O.co2",
            "charge": -0.6905415761361866,
            "radius": 1.76,
        },
        {
            "name": "CAC",
            "type": "C.3",
            "charge": 0.012573264939255872,
            "radius": 1.87,
        },
        {
            "name": "HAC",
            "type": "H",
            "charge": 0.11848768797200306,
            "radius": 1.1,
        },
        {
            "name": "HAD",
            "type": "H",
            "charge": 0.11848768797200306,
            "radius": 1.1,
        },
        {
            "name": "CAD",
            "type": "C.3",
            "charge": -0.07742980670292797,
            "radius": 1.87,
        },
        {
            "name": "HAE",
            "type": "H",
            "charge": 0.06140059278092227,
            "radius": 1.1,
        },
        {
            "name": "HAF",
            "type": "H",
            "charge": 0.06140059278092227,
            "radius": 1.1,
        },
        {
            "name": "CAE",
            "type": "C.3",
            "charge": -0.076959643935905,
            "radius": 1.87,
        },
        {
            "name": "HAG",
            "type": "H",
            "charge": 0.07285165771612633,
            "radius": 1.1,
        },
        {
            "name": "HAH",
            "type": "H",
            "charge": 0.07285165771612633,
            "radius": 1.1,
        },
        {
            "name": "CAF",
            "type": "C.2",
            "charge": -0.24236491728150775,
            "radius": 1.87,
        },
        {
            "name": "HAI",
            "type": "H",
            "charge": 0.16485599597989423,
            "radius": 1.1,
        },
        {
            "name": "CAG",
            "type": "C.2",
            "charge": -0.2426026937449157,
            "radius": 1.87,
        },
        {
            "name": "HAJ",
            "type": "H",
            "charge": 0.1648388876026803,
            "radius": 1.1,
        },
        {
            "name": "CAH",
            "type": "C.3",
            "charge": -0.07965667787077134,
            "radius": 1.87,
        },
        {
            "name": "HAK",
            "type": "H",
            "charge": 0.07232360507910963,
            "radius": 1.1,
        },
        {
            "name": "HAL",
            "type": "H",
            "charge": 0.07232360507910963,
            "radius": 1.1,
        },
        {
            "name": "CAI",
            "type": "C.3",
            "charge": -0.09862724378017657,
            "radius": 1.87,
        },
        {
            "name": "HAM",
            "type": "H",
            "charge": 0.05403847507524275,
            "radius": 1.1,
        },
        {
            "name": "HAN",
            "type": "H",
            "charge": 0.05403847507524275,
            "radius": 1.1,
        },
        {
            "name": "CAJ",
            "type": "C.3",
            "charge": -0.10608184708836539,
            "radius": 1.87,
        },
        {
            "name": "HAO",
            "type": "H",
            "charge": 0.04966180988232633,
            "radius": 1.1,
        },
        {
            "name": "HAP",
            "type": "H",
            "charge": 0.04966180988232633,
            "radius": 1.1,
        },
        {
            "name": "CAK",
            "type": "C.3",
            "charge": -0.11173869603286389,
            "radius": 1.87,
        },
        {
            "name": "HAR",
            "type": "H",
            "charge": 0.04061118893012504,
            "radius": 1.1,
        },
        {
            "name": "HAS",
            "type": "H",
            "charge": 0.04061118893012504,
            "radius": 1.1,
        },
        {
            "name": "HAQ",
            "type": "H",
            "charge": 0.04061118893012504,
            "radius": 1.1,
        },
    ],
    "glycerol.mol2": [
        {
            "name": "OAB",
            "type": "O.3",
            "charge": -0.6193219106971865,
            "radius": 1.52,
        },
        {
            "name": "HAG",
            "type": "H",
            "charge": 0.4033320347133161,
            "radius": 1.1,
        },
        {
            "name": "CAA",
            "type": "C.3",
            "charge": -0.025227212803925057,
            "radius": 1.87,
        },
        {
            "name": "HAA",
            "type": "H",
            "charge": 0.1316040991547254,
            "radius": 1.1,
        },
        {
            "name": "HAB",
            "type": "H",
            "charge": 0.1316040991547254,
            "radius": 1.1,
        },
        {
            "name": "CAC",
            "type": "C.3",
            "charge": -0.00300604112002907,
            "radius": 1.87,
        },
        {
            "name": "HAC",
            "type": "H",
            "charge": 0.16288785945019288,
            "radius": 1.1,
        },
        {
            "name": "OAD",
            "type": "O.3",
            "charge": -0.6115908338990877,
            "radius": 1.52,
        },
        {
            "name": "HAD",
            "type": "H",
            "charge": 0.4077267965256133,
            "radius": 1.1,
        },
        {
            "name": "CAE",
            "type": "C.3",
            "charge": -0.025227212803925054,
            "radius": 1.87,
        },
        {
            "name": "HAE",
            "type": "H",
            "charge": 0.1316040991547254,
            "radius": 1.1,
        },
        {
            "name": "HAF",
            "type": "H",
            "charge": 0.1316040991547254,
            "radius": 1.1,
        },
        {
            "name": "OAF",
            "type": "O.3",
            "charge": -0.6193219106971865,
            "radius": 1.52,
        },
        {
            "name": "HAH",
            "type": "H",
            "charge": 0.4033320347133161,
            "radius": 1.1,
        },
    ],
    "naphthalene.mol2": [
        {
            "name": "CAB",
            "type": "C.ar",
            "charge": -0.10985826436427909,
            "radius": 1.87,
        },
        {
            "name": "HAB",
            "type": "H",
            "charge": 0.13714843913419308,
            "radius": 1.1,
        },
        {
            "name": "CAC",
            "type": "C.ar",
            "charge": -0.06352790423000582,
            "radius": 1.87,
        },
        {
            "name": "CAD",
            "type": "C.ar",
            "charge": -0.10985826436427909,
            "radius": 1.87,
        },
        {
            "name": "HAD",
            "type": "H",
            "charge": 0.13714843913419308,
            "radius": 1.1,
        },
        {
            "name": "CAE",
            "type": "C.ar",
            "charge": -0.12522022371734542,
            "radius": 1.87,
        },
        {
            "name": "HAE",
            "type": "H",
            "charge": 0.1296940010624343,
            "radius": 1.1,
        },
        {
            "name": "CAJ",
            "type": "C.ar",
            "charge": -0.12522022371734542,
            "radius": 1.87,
        },
        {
            "name": "HAJ",
            "type": "H",
            "charge": 0.1296940010624343,
            "radius": 1.1,
        },
        {
            "name": "CAI",
            "type": "C.ar",
            "charge": -0.10985826436427909,
            "radius": 1.87,
        },
        {
            "name": "HAI",
            "type": "H",
            "charge": 0.13714843913419308,
            "radius": 1.1,
        },
        {
            "name": "CAH",
            "type": "C.ar",
            "charge": -0.06352790423000582,
            "radius": 1.87,
        },
        {
            "name": "CAG",
            "type": "C.ar",
            "charge": -0.10985826436427909,
            "radius": 1.87,
        },
        {
            "name": "HAG",
            "type": "H",
            "charge": 0.13714843913419308,
            "radius": 1.1,
        },
        {
            "name": "CAF",
            "type": "C.ar",
            "charge": -0.12522022371734542,
            "radius": 1.87,
        },
        {
            "name": "HAF",
            "type": "H",
            "charge": 0.1296940010624343,
            "radius": 1.1,
        },
        {
            "name": "CAA",
            "type": "C.ar",
            "charge": -0.12522022371734542,
            "radius": 1.87,
        },
        {
            "name": "HAA",
            "type": "H",
            "charge": 0.1296940010624343,
            "radius": 1.1,
        },
    ],
    "pyrrole.mol2": [
        {
            "name": "CAE",
            "type": "C.3",
            "charge": -0.02725132001441914,
            "radius": 1.87,
        },
        {
            "name": "HAE",
            "type": "H",
            "charge": 0.104834463409739,
            "radius": 1.1,
        },
        {
            "name": "HAF",
            "type": "H",
            "charge": 0.104834463409739,
            "radius": 1.1,
        },
        {
            "name": "CAD",
            "type": "C.2",
            "charge": -0.19753848158167037,
            "radius": 1.87,
        },
        {
            "name": "HAD",
            "type": "H",
            "charge": 0.18534981863154645,
            "radius": 1.1,
        },
        {
            "name": "CAA",
            "type": "C.2",
            "charge": -0.1112561716335608,
            "radius": 1.87,
        },
        {
            "name": "HAA",
            "type": "H",
            "charge": 0.2699321111667989,
            "radius": 1.1,
        },
        {
            "name": "NAB",
            "type": "N.2",
            "charge": -0.44475473234270757,
            "radius": 1.4,
        },
        {
            "name": "CAC",
            "type": "C.2",
            "charge": -0.13358085588726296,
            "radius": 1.87,
        },
        {
            "name": "HAC",
            "type": "H",
            "charge": 0.24943070484179752,
            "radius": 1.1,
        },
    ],
    "tetramethylammonium.mol2": [
        {
            "name": "CAA",
            "type": "C.3",
            "charge": 0.013840200623221964,
            "radius": 1.87,
        },
        {
            "name": "HAA",
            "type": "H",
            "charge": 0.12997929051575582,
            "radius": 1.1,
        },
        {
            "name": "HAB",
            "type": "H",
            "charge": 0.12997929051575582,
            "radius": 1.1,
        },
        {
            "name": "HAC",
            "type": "H",
            "charge": 0.12997929051575582,
            "radius": 1.1,
        },
        {
            "name": "NAC",
            "type": "N.3",
            "charge": -0.6151122886819584,
            "radius": 1.4,
        },
        {
            "name": "CAB",
            "type": "C.3",
            "charge": 0.013840200623221978,
            "radius": 1.87,
        },
        {
            "name": "HAE",
            "type": "H",
            "charge": 0.12997929051575582,
            "radius": 1.1,
        },
        {
            "name": "HAF",
            "type": "H",
            "charge": 0.12997929051575582,
            "radius": 1.1,
        },
        {
            "name": "HAD",
            "type": "H",
            "charge": 0.12997929051575582,
            "radius": 1.1,
        },
        {
            "name": "CAE",
            "type": "C.3",
            "charge": 0.013840200623221978,
            "radius": 1.87,
        },
        {
            "name": "HAK",
            "type": "H",
            "charge": 0.12997929051575582,
            "radius": 1.1,
        },
        {
            "name": "HAL",
            "type": "H",
            "charge": 0.12997929051575582,
            "radius": 1.1,
        },
        {
            "name": "HAJ",
            "type": "H",
            "charge": 0.12997929051575582,
            "radius": 1.1,
        },
        {
            "name": "CAD",
            "type": "C.3",
            "charge": 0.013840200623221978,
            "radius": 1.87,
        },
        {
            "name": "HAH",
            "type": "H",
            "charge": 0.12997929051575582,
            "radius": 1.1,
        },
        {
            "name": "HAI",
            "type": "H",
            "charge": 0.12997929051575582,
            "radius": 1.1,
        },
        {
            "name": "HAG",
            "type": "H",
            "charge": 0.12997929051575582,
            "radius": 1.1,
        },
    ],
}

STRING_1HPX = """HETATM 1    C1    LIG  1    -0.0430817735  8.2957558296
HETATM 2    N1    LIG  1    -0.3567806953  8.8642031247
HETATM 3    O1    LIG  1    -0.5711988944  9.1244919327
HETATM 4    S1    LIG  1    -0.3623072542  8.4098272259
HETATM 5    C2    LIG  1    -0.0754644800  8.0825452574
HETATM 6    N2    LIG  1    -0.7033368408  9.2085255331
HETATM 7    O2    LIG  1    -0.5909676444  8.7116035261
HETATM 8    S2    LIG  1    -0.3132656813  8.7605435981
HETATM 9    C3    LIG  1     0.0121840362  8.6192189498
HETATM 10   N3    LIG  1    -0.7139711104  9.0429239921
HETATM 11   O3    LIG  1    -0.2919660371 11.7608669889
HETATM 12   C4    LIG  1    -0.0248911126  8.3975503043
HETATM 13   N4    LIG  1    -0.7120419699  9.1281779687
HETATM 14   O4    LIG  1    -0.2689251910 11.9726798365
HETATM 15   C5    LIG  1    -0.0272176845  8.3884310504
HETATM 16   N5    LIG  1    -0.7179288701  8.9894054162
HETATM 17   O5    LIG  1    -0.2884010362 11.7926192375
HETATM 18   C6    LIG  1     0.0753432558  9.0487754381
HETATM 19   O6    LIG  1    -0.2912545917 11.7702834794
HETATM 20   C7    LIG  1    -0.0667375437  8.1434847512
HETATM 21   C8    LIG  1    -0.1139105991  7.8613684842
HETATM 22   C9    LIG  1    -0.1032463708  7.9252346777
HETATM 23   C10   LIG  1     0.0955779782  8.5536222328
HETATM 24   C11   LIG  1     0.1332593724 10.1521664812
HETATM 25   C12   LIG  1    -0.0074544800  7.9366419563
HETATM 26   C13   LIG  1     0.1300074666 10.1317750510
HETATM 27   C14   LIG  1     0.0331607920  8.1612521219
HETATM 28   C15   LIG  1     0.1162788490  8.6734320906
HETATM 29   C16   LIG  1    -0.0501204845  7.6824020199
HETATM 30   C17   LIG  1     0.1762998815 10.4714335775
HETATM 31   C18   LIG  1     0.1273124481  8.7567099170
HETATM 32   C19   LIG  1     0.0097183646  8.0374376393
HETATM 33   C20   LIG  1     0.0537322435  8.3045888878
HETATM 34   C21   LIG  1     0.1423678189 10.2168467012
HETATM 35   C22   LIG  1    -0.0188583562  7.8769054795
HETATM 36   C23   LIG  1    -0.0887238757  7.4582556830
HETATM 37   C24   LIG  1    -0.0887238757  7.4582556830
HETATM 38   C25   LIG  1    -0.0887238757  7.4582556830
HETATM 39   C26   LIG  1     0.1059315385  8.5989689579
HETATM 40   C27   LIG  1    -0.0578275062  7.6467371471
HETATM 41   C28   LIG  1    -0.0994004929  7.9768121661
HETATM 42   C29   LIG  1    -0.1178824768  7.8502328918
HETATM 43   C30   LIG  1    -0.1266932100  7.7958113822
HETATM 44   C31   LIG  1    -0.1284021538  7.7861998095
HETATM 45   C32   LIG  1    -0.1266932100  7.7958113822
HETATM 46   C33   LIG  1    -0.1178824768  7.8502328918
HETATM 47   HN2   LIG  1     0.4736506271  9.3793810920
HETATM 48   HN3   LIG  1     0.4681379918  9.3514454934
HETATM 49   HN5   LIG  1     0.4699450251  9.3649880489
HETATM 50   HO2   LIG  1     0.4135083336  9.0347534947
HETATM 51   H1    LIG  1     0.1930651515  7.9310778118
HETATM 52   H2    LIG  1     0.1493049838  7.7494684898
HETATM 53   H5    LIG  1     0.2004887587  7.9579791263
HETATM 54   H7    LIG  1     0.1568207646  7.7783750214
HETATM 55   H8    LIG  1     0.1323481037  7.6911054255
HETATM 56   H9    LIG  1     0.1383298871  7.7127230099
HETATM 57   1H10  LIG  1     0.2090610332  7.9845419050
HETATM 58   2H10  LIG  1     0.2090610332  7.9845419050
HETATM 59   1H12  LIG  1     0.1211146204  7.6358478484
HETATM 60   2H12  LIG  1     0.1211146204  7.6358478484
HETATM 61   H14   LIG  1     0.1675599590  7.8253253808
HETATM 62   H15   LIG  1     0.2276083861  8.0619341012
HETATM 63   1H16  LIG  1     0.0852867335  7.4987774313
HETATM 64   2H16  LIG  1     0.0852867335  7.4987774313
HETATM 65   H18   LIG  1     0.2431703175  8.1293441537
HETATM 66   1H19  LIG  1     0.1277656093  7.6574434197
HETATM 67   2H19  LIG  1     0.1277656093  7.6574434197
HETATM 68   1H20  LIG  1     0.1956533010  7.9423134276
HETATM 69   2H20  LIG  1     0.1956533010  7.9423134276
HETATM 70   1H23  LIG  1     0.0567190229  7.3921792295
HETATM 71   2H23  LIG  1     0.0567190229  7.3921792295
HETATM 72   3H23  LIG  1     0.0567190229  7.3921792295
HETATM 73   1H24  LIG  1     0.0567190229  7.3921792295
HETATM 74   2H24  LIG  1     0.0567190229  7.3921792295
HETATM 75   3H24  LIG  1     0.0567190229  7.3921792295
HETATM 76   1H25  LIG  1     0.0567190229  7.3921792295
HETATM 77   2H25  LIG  1     0.0567190229  7.3921792295
HETATM 78   3H25  LIG  1     0.0567190229  7.3921792295
HETATM 79   H26   LIG  1     0.2127692970  7.9981857070
HETATM 80   1H27  LIG  1     0.0838296768  7.4951353599
HETATM 81   2H27  LIG  1     0.0838296768  7.4951353599
HETATM 82   3H27  LIG  1     0.0838296768  7.4951353599
HETATM 83   H29   LIG  1     0.1344214538  7.7011379372
HETATM 84   H30   LIG  1     0.1294413566  7.6832603969
HETATM 85   H31   LIG  1     0.1289517709  7.6818306603
HETATM 86   H32   LIG  1     0.1294413566  7.6832603969
HETATM 87   H33   LIG  1     0.1344214538  7.7011379372
"""
CHARGES_1HPX = []
for line in StringIO(STRING_1HPX):
    words = line.split()
    charge = float(words[5])
    CHARGES_1HPX.append(charge)
