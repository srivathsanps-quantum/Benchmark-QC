import numpy as np
from pyscf import gto, scf, mcscf, fci
from pyscf.mcscf.addons import project_init_guess
import basis_set_exchange as bse
import pyscf
import pennylane as qml 

import os
Hs = []
Es = []


def H_gen(basis_input, elements, geom, spin, charge, ncas, nelecas,
          save=True, savefile="H_data.npz", geom_id=None):
    # Define bond distance range for potential energy surface
    hf_energies = []
    casci_energies = []

    print("Computing U2 potential energy surface...")
    print("Bond Distance (Å)     HF Energy (Hartree)      CASCI Energy (Hartree)")
    print("-" * 75)

    # ---------- Step 1: Reference CASSCF calculation ----------
    mol_ref = gto.Mole()
    mol_ref.atom = '''
    U 0 0 0
    U 0 0 2.5
    '''
    # Use a built-in basis (fix for BasisNotFoundError)
    mol_ref.basis = basis_input        # or try 'def2-TZVP' / 'ano-rcc-mb' if available
    mol_ref.charge = 0
    mol_ref.spin = 0
    mol_ref.symmetry = False
    mol_ref.build()

    # Perform ROHF + X2C (spin-free)
    mf_ref = scf.ROHF(mol_ref).sfx2c1e()
    mf_ref.level_shift = 0.5
    mf_ref.diis_space = 12
    mf_ref.max_cycle = 100
    #mf_ref.verbose = 0
    mf_ref.kernel()
    if not mf_ref.converged:
        mf_ref = scf.newton(mf_ref).run()

    # One-shot CASSCF to get proper active orbitals
    mycas_ref = mcscf.CASSCF(mf_ref, ncas, nelecas)
    mycas_ref.fix_spin_(ss=0.0)
    en = mycas_ref.kernel()
    print('Ref.CASSCF energy:', en[0])
    mo_ref = mycas_ref.mo_coeff
    print('\n')

    # ---------- Step 2: CASCI scan with projected orbitals ----------
    mo_prev = None
    print('-----------------Bond distances ---------', geom)
    mol = gto.Mole()
    mol.atom = geom
    mol.basis = basis_input          # keep same basis as ref
    mol.charge = charge
    mol.spin = spin
    mol.symmetry = False
    mol.build()

    # ROHF + X2C (spin-free)
    mf = scf.ROHF(mol).sfx2c1e()
    mf.level_shift = 0.5
    mf.diis_space = 12
    mf.max_cycle = 100
    if mo_prev is not None:
        #print('MO prev is used')
        mf.mo_coeff = mo_prev
    hf_energy = mf.kernel()
    if not mf.converged:
        mf = scf.newton(mf).run()
    hf_energy = mf.e_tot


    from pyscf.mcscf.addons import project_init_guess

    mycas = mcscf.CASCI(mf, ncas, nelecas)
    mo_proj = project_init_guess(mycas, mo_ref, prev_mol=mol_ref)
    #print('-------------------Going to get the mo_proj---------------------')
    #print('The bond length I am working on', u2_bond)
    


    from pyscf.fci import direct_spin1, addons
    fcis = direct_spin1.FCI(mol)
    fcis.spin = 0
    fcis.nroots = 1
    fcis = addons.fix_spin_(fcis, ss=0.0, shift=0.8)
    mycas.fcisolver = fcis

    # stick to singlet and follow same root
    mycas.fix_spin_(ss=0.0)
    
    h1ecas, ecore = mycas.get_h1eff(mo_proj)
    h2ecas = mycas.get_h2eff(mo_proj)





    casci_energy = mycas.kernel(mo_coeff=mo_proj)[0]
    print('casci_energy', casci_energy)

    hf_energies.append(hf_energy)
    casci_energies.append(casci_energy)
    mo_prev = mf.mo_coeff

    #print(f"{u2_bond:8.1f}        {hf_energy:12.8f}           {casci_energy:12.8f}")
    print('\n')

  



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