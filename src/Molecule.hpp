#ifndef MOLECULE_H_
#define MOLECULE_H_

#include "common.hpp"
#include "Atom.hpp"
#include "Vector3d.hpp"
#include <set>
#include <map>
#include <iostream>

// #define FRAGMENTTEST

namespace fragdock {
  struct Bond {
    int atom_id1, atom_id2;
    bool is_rotor;
    friend std::ostream& operator<< (std::ostream& os, const Bond& bond);
    Bond(const int atom_id1, const int atom_id2, const bool is_rotor) : atom_id1(atom_id1), atom_id2(atom_id2), is_rotor(is_rotor) {}
  };
  class Molecule {
  protected:
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    // std::vector<std::vector<int> > bond_ids;
    // fltype radius;
    int heavy_num;
    std::string title, smiles, identifier;
    fltype intraEnergy;
  public:
    std::vector<std::vector<int> > bond_ids; //結合箇所を辺とみなしたグラフ?
    Molecule(): heavy_num(-1) {}
    Molecule(const std::vector<Atom> &atoms, const std::string title, const std::string smiles);
    // void calcRadius();
    fltype getRadius() const;
    void translate(const Vector3d &vec);
    void rotate(fltype theta, fltype phi, fltype psi);
    void rotate(const Vector3d& vec);
    // void axisRotate(const Vector3d& axis, fltype th);
    // void axisRotate(const Vector3d& base, const Vector3d& axis, fltype th, const std::vector<int>& id_set);
    void append(const Molecule &o);
    void append(const Atom &o);
    void append(const Bond &bond);
    Vector3d getCenter() const;
    int size() const { return atoms.size(); }
    int heavyNum();
    const Atom& getAtom(int i) const { return atoms[i]; }
    const std::vector<Atom>& getAtoms() const { return atoms; }
    const std::vector<Bond>& getBonds() const { return bonds; }
    friend std::ostream& operator<< (std::ostream& os, const Molecule& mol);

    void renumbering(int newsz, const std::vector<unsigned int>& vec);
    void deleteHydrogens();

    const std::string& gettitle() const { return title; }
    const std::string& getsmiles() const { return smiles; }
    const std::string& getIdentifier() const { return identifier; }

#ifdef FRAGMENTTEST
    void settitle(const std::string& _title) { title = _title; }
    void setsmiles(const std::string& _smiles) { smiles = _smiles; }
    void setIdentifier(const std::string& _identifier) { identifier = _identifier; }
#endif

    std::vector<std::vector<int> > getGraphDistances() const;
    void setIntraEnergy(fltype energy) { intraEnergy = energy; }
    fltype getIntraEnergy() const { return intraEnergy; }
    fltype getNrots() const;

    bool operator<(const Molecule& o) const { return identifier < o.identifier; }
    bool operator>(const Molecule& o) const { return identifier > o.identifier; }

  };
}
#endif
