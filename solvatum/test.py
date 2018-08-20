
import ui

d = ui.Database()

def test_atoms():
    assert d.atoms_in_sol('Hexane') == 20

def test_acomatic():
    assert d.sol_is_aromatic('Anisole') == True



