#include "Atom.hpp"

namespace fragdock {
  std::ostream& operator<< (std::ostream& os, const Atom& atom){
    os << "atom_id: " << atom.getId() << " "
       << "atom_type: " << xs_name(atom.getXSType()) << " "
       << "coordination: (" << atom.x << "," << atom.y << "," << atom.z << ") ";
      return os;
  }
}
