import numpy as np
from pyscf import gto, scf, mcscf, fci
from pyscf.mcscf.addons import project_init_guess
import basis_set_exchange as bse
import pyscf
import pennylane as qml 
from pyscf.fci import direct_spin1, addons


import os
Hs = []
Es = []
hf_energies = []
casci_energies = []


#Building reference calculation

def build_u2_reference(basis_input, ncas, nelecas):
    mol_ref = gto.Mole()
    mol_ref.atom = """
    U 0 0 0
    U 0 0 2.5
    """
    mol_ref.basis = basis_input
    mol_ref.charge = 0
    mol_ref.spin = 0
    mol_ref.symmetry = False
    mol_ref.build()
    mol_ref.verbose = 0

    mf_ref = scf.ROHF(mol_ref).sfx2c1e()
    mf_ref.level_shift = 0.5
    mf_ref.diis_space = 12
    mf_ref.max_cycle = 100
    mf_ref.verbose = 0
    mf_ref.kernel()
    if not mf_ref.converged:
        mf_ref = scf.newton(mf_ref).run()

    mycas_ref = mcscf.CASSCF(mf_ref, ncas, nelecas)
    mycas_ref.fix_spin_(ss=0.0)
    mycas_ref.verbose = 0
    en = mycas_ref.kernel()

    mo_ref = mycas_ref.mo_coeff
    return mol_ref, mo_ref



def H_gen(basis_input, elements, geom, spin, charge, ncas, nelecas,
          mol_ref, mo_ref,
          save=True, savefile="H_data.npz", geom_id=None):
    # only CASCI scan part; remove the whole mol_ref / mf_ref / mycas_ref block

    mol = gto.Mole()
    mol.atom = geom
    mol.basis = basis_input
    mol.charge = charge
    mol.spin = spin
    mol.symmetry = False
    mol.build()
    mol.verbose = 0

    mf = scf.ROHF(mol).sfx2c1e()
    mf.level_shift = 0.5
    mf.diis_space = 12
    mf.max_cycle = 100
    mf.verbose = 0
    hf_energy = mf.kernel()
    if not mf.converged:
        mf = scf.newton(mf).run()
    hf_energy = mf.e_tot

    mycas = mcscf.CASCI(mf, ncas, nelecas)
    from pyscf.mcscf.addons import project_init_guess
    mo_proj = project_init_guess(mycas, mo_ref, prev_mol=mol_ref)
    fcis = direct_spin1.FCI(mol)
    fcis.spin = 0
    fcis.nroots = 1
    fcis = addons.fix_spin_(fcis, ss=0.0, shift=0.8)
    mycas.fcisolver = fcis

    # stick to singlet and follow same root
    mycas.fix_spin_(ss=0.0)
    
    h1ecas, ecore = mycas.get_h1eff(mo_proj)
    h2ecas = mycas.get_h2eff(mo_proj)
    mycas.verbose = 0
    casci_energy = mycas.kernel(mo_coeff=mo_proj)[0]

    hf_energies.append(hf_energy)
    casci_energies.append(casci_energy)
    mo_prev = mf.mo_coeff

    two_mo = pyscf.ao2mo.restore('1', h2ecas, norb=ncas)
    two_mo = np.swapaxes(two_mo, 1, 3)


    one_mo = h1ecas
    core_constant = np.array([ecore]) 

    H_fermionic = qml.qchem.fermionic_observable(core_constant, one_mo, two_mo, cutoff=1e-20)


    H = qml.jordan_wigner(H_fermionic)
    

    if save:
        if geom_id is None:
            label = geom
        else:
            label = geom_id

        if os.path.exists(savefile):
            data = np.load(savefile, allow_pickle=True)
            labels = list(data["labels"])
            Hs = list(data["Hs"])
            Es = list(data["casci_energies"])
        else:
            labels, Hs, Es = [], [], []

        labels.append(label)
        Hs.append(H)
        Es.append(casci_energy)

        np.savez_compressed(
            savefile,
            labels=np.array(labels, dtype=object),
            Hs=np.array(Hs, dtype=object),
            casci_energies=np.array(Es),
        )

    return H, casci_energy




