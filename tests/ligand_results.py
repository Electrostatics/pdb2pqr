"""Expected results for ligand tests"""
TORSION_RESULTS = {
    "ethanol.mol2": {
        ('CAA', 'CAB', 'OAC', 'HAF'), ('HAF', 'OAC', 'CAB', 'CAA')},
    "glycerol.mol2": {
        ('CAA', 'CAB', 'CAC', 'OAF'), ('CAA', 'CAB', 'OAE', 'HAG'),
        ('CAB', 'CAA', 'OAD', 'HAF'), ('CAB', 'CAC', 'OAF', 'HAH'),
        ('CAC', 'CAB', 'CAA', 'OAD'), ('CAC', 'CAB', 'OAE', 'HAG'),
        ('HAF', 'OAD', 'CAA', 'CAB'), ('HAG', 'OAE', 'CAB', 'CAA'),
        ('HAG', 'OAE', 'CAB', 'CAC'), ('HAH', 'OAF', 'CAC', 'CAB'),
        ('OAD', 'CAA', 'CAB', 'CAC'), ('OAD', 'CAA', 'CAB', 'OAE'),
        ('OAE', 'CAB', 'CAA', 'OAD'), ('OAE', 'CAB', 'CAC', 'OAF'),
        ('OAF', 'CAC', 'CAB', 'CAA'), ('OAF', 'CAC', 'CAB', 'OAE')},
    "cyclohexane.mol2": {
        ('CAA', 'CAB', 'CAC', 'CAF'), ('CAA', 'CAD', 'CAE', 'CAF'),
        ('CAB', 'CAA', 'CAD', 'CAE'), ('CAB', 'CAC', 'CAF', 'CAE'),
        ('CAC', 'CAB', 'CAA', 'CAD'), ('CAC', 'CAF', 'CAE', 'CAD'),
        ('CAD', 'CAA', 'CAB', 'CAC'), ('CAD', 'CAE', 'CAF', 'CAC'),
        ('CAE', 'CAD', 'CAA', 'CAB'), ('CAE', 'CAF', 'CAC', 'CAB'),
        ('CAF', 'CAC', 'CAB', 'CAA'), ('CAF', 'CAE', 'CAD', 'CAA')}
}


RING_RESULTS = {
    "ethanol.mol2": set(),
    "glycerol.mol2": set(),
    "cyclohexane.mol2": {('CAA', 'CAD', 'CAE', 'CAF', 'CAC', 'CAB')},
    "naphthalene.mol2": {
        ('CAA', 'CAB', 'CAC', 'CAH', 'CAG', 'CAF'),
        ('CAC', 'CAH', 'CAI', 'CAJ', 'CAE', 'CAD')},
    "anthracene.mol2": {
        ('CAC', 'CAJ', 'CAK', 'CAL', 'CAE', 'CAD'),
        ('CAE', 'CAL', 'CAM', 'CAN', 'CAG', 'CAF'),
        ('CAA', 'CAB', 'CAC', 'CAJ', 'CAI', 'CAH')},
    "crown-ether.mol2": {
        ('CAA', 'CAE', 'CAF', 'CAG', 'CAC', 'CAB'),
        ('CAD', 'CAI', 'CAJ', 'CAK', 'CAF', 'CAE'),
        ('CAF', 'CAG', 'CAH', 'CAM', 'CAL', 'CAK')}
}


BOND_RESULTS = {
    "cyclohexane.mol2": 6 * ["single"],
    "ethanol.mol2": [
        "single", "single", None, "single", "single"],
    "glycerol.mol2": [
        None, "single", "single", "single", "single", None, "single", None],
    "acetylcholine.mol2": [
        "single", "double", "single", "single", "single", "single", "single",
        "single", "single"],
    "acetonitrile.mol2": [
        "triple", "single"],
    "pyrrole.mol2": [
        "aromatic", "aromatic", "aromatic", "aromatic", "aromatic", None],
    "fatty-acid.mol2": [
        "double", "double", "single", "single", "single", "single", "double",
        "single", "single", "single", "single"],
    "trimethylamine.mol2": ["single", None, "single", "single"],
    "naphthalene.mol2": 11 * ["aromatic"]
}


FORMAL_CHARGE_RESULTS = {
    "1HPX-ligand.mol2": 87*[0],
    "1QBS-ligand.mol2": 80*[0],
    "1US0-ligand.mol2": 22*[0] + [-0.5, -0.5] + 11*[0],
    "acetate.mol2": [-0.5, 0, -0.5] + 4*[0],
    "acetylcholine.mol2": 13*[0] + [1] + 12*[0],
    "adp.mol2": [-1] + 5*[0] + [-1] + 32*[0],
    "anthracene.mol2": 24*[0],
    "acetonitrile.mol2": 6*[0],
    "cyclohexane.mol2": 18*[0],
    "ethanol.mol2": 9*[0],
    "fatty-acid.mol2": [-0.5, 0, -0.5] + 26*[0],
    "glycerol.mol2": 14*[0],
    "naphthalene.mol2": 18*[0],
    "crown-ether.mol2": 42*[0],
    "pyrrole.mol2": 10*[0],
    "trimethylamine.mol2": 4*[0] + [1] + 12*[0]
}
