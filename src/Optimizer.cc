#include "Optimizer.hpp"

namespace fragdock {
  fltype Optimizer::optimize(Molecule &mol, const EnergyCalculator& ec) const {
    using namespace std;
    srand(0); // fix seed

    // hill climbing
    fltype opt = ec.getEnergy(mol, receptor);
    const fltype pi = acos(-1.0);
    const int NEAREST_NUM = 200;
    fltype trans_step = 0.5;
    fltype rotate_step = pi / 30.0;
    Molecule first_mol = mol;
    while(true) {
      fltype next_val = 1e10;
      Molecule next_mol = mol;
      Vector3d center = mol.getCenter();
      for(int i = 0; i < NEAREST_NUM; i++) {
        Vector3d dv = Vector3d(utils::randf(-trans_step, trans_step),
                               utils::randf(-trans_step, trans_step),
                               utils::randf(-trans_step, trans_step));
        //dv = dv.unit();
        fltype th  = utils::randf(-rotate_step, rotate_step);
        fltype phi = utils::randf(-rotate_step, rotate_step);
        fltype psi = utils::randf(-rotate_step, rotate_step);

        Molecule tmp_mol = mol;
        tmp_mol.translate(-center);
        tmp_mol.rotate(th, phi, psi);
        tmp_mol.translate(center + dv);
        // tmp_mol.translate(dv);

        if (first_mol.calcRMSD(tmp_mol) > max_rmsd) continue;

        fltype val = ec.getEnergy(tmp_mol, receptor);
        if(val < next_val) {
          next_val = val;
          next_mol = tmp_mol;
        }
      }
      //cout << "next_val: " << next_val << endl;
      if(next_val < opt) {
        opt = next_val;
        mol = next_mol;
      } else break;
    }

    return opt;
  }

  fltype Optimizer_Grid::optimize(Molecule &mol) const {
    using namespace std;
    srand(0); // fix seed

    // hill climbing
    fltype opt = calcInterEnergy(mol) + mol.getIntraEnergy();
    const fltype pi = acos(-1.0);
    const int NEAREST_NUM = 200;
    fltype trans_step = 0.5;
    fltype rotate_step = pi / 30.0;
    Molecule first_mol = mol;
    while(true) {
      fltype next_val = 1e10;
      Molecule next_mol = mol;
      Vector3d center = mol.getCenter();
      for(int i = 0; i < NEAREST_NUM; i++) {
        Vector3d dv = Vector3d(utils::randf(-trans_step, trans_step),
                               utils::randf(-trans_step, trans_step),
                               utils::randf(-trans_step, trans_step));
        //dv = dv.unit();
        fltype th  = utils::randf(-rotate_step, rotate_step);
        fltype phi = utils::randf(-rotate_step, rotate_step);
        fltype psi = utils::randf(-rotate_step, rotate_step);

        Molecule tmp_mol = mol;
        tmp_mol.translate(-center);
        tmp_mol.rotate(th, phi, psi);
        tmp_mol.translate(center + dv);
        // tmp_mol.translate(dv);

        if (first_mol.calcRMSD(tmp_mol) > max_rmsd) continue;

        fltype val = calcInterEnergy(tmp_mol) + mol.getIntraEnergy();
        if(val < next_val) {
          next_val = val;
          next_mol = tmp_mol;
        }
      }
      //cout << "next_val: " << next_val << endl;
      if(next_val < opt) {
        opt = next_val;
        mol = next_mol;
      } else break;
    }

    return opt;
  }
}
