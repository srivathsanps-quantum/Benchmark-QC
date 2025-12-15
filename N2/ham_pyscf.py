import basis_set_exchange as bse
import pyscf
from pyscf import gto, scf, mcscf
import pennylane as qml
import numpy as np
import os

def H_gen(basis_input, elements, geom, spin, charge, ncas, nelecas,
          save=True, savefile="H_data.npz", geom_id=None):
    """
    Generate H and CASCI energy.
    If save=True, append this point to `savefile`.
    geom_id can be any label (e.g. bond length) stored with the data.
    """

    basis = bse.get_basis(basis_input, elements=elements, fmt="nwchem")

    mol = gto.Mole()
    mol.atom = geom
    mol.basis = basis
    mol.charge = charge
    mol.spin = spin
    mol.build()

    if mol.spin == 0:
        mf = scf.RHF(mol)
    else:
        mf = scf.ROHF(mol)
    mf.verbose = 0
    hf = mf.kernel()

    if not mf.converged:
        mf = mf.newton(mf).run()

    mycas = mcscf.CASCI(mf, ncas, nelecas)
    mycas.verbose = 0
    cas_energy = mycas.kernel()[0]

    one_mo, ecore = mycas.get_h1eff(mycas.mo_coeff)
    h2ecas = mycas.get_h2eff(mycas.mo_coeff)
    two_mo = pyscf.ao2mo.restore("1", h2ecas, norb=ncas)
    two_mo = np.swapaxes(two_mo, 1, 3)

    core_constant = np.array([ecore])
    H_fermionic = qml.qchem.fermionic_observable(
        core_constant, one_mo, two_mo, cutoff=1e-20
    )
    H = qml.jordan_wigner(H_fermionic)

    # -------- saving logic (single point or PES) ----------
    if save:
        # use geom_id if given, else store the raw geom string
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
        Es.append(cas_energy)

        np.savez_compressed(
            savefile,
            labels=np.array(labels, dtype=object),
            Hs=np.array(Hs, dtype=object),
            casci_energies=np.array(Es),
        )

    return H, cas_energy
