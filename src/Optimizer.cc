#include "Optimizer.hpp"

namespace fragdock {
  fltype Optimizer::optimize(Molecule &mol, const EnergyCalculator& ec) const {
    using namespace std;

    // hill climbing
    fltype opt = ec.getEnergy(mol, receptor);
    const fltype pi = acos(-1.0);
    const int NEAREST_NUM = 200;
    fltype trans_step = 0.5;
    fltype rotate_step = pi / 30.0;
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
        tmp_mol.translate(-center); //rotateは原点中心の回転なので目的リガンドを原点に一度持ってきて回転させると楽ということ．
        tmp_mol.rotate(th, phi, psi);
        tmp_mol.translate(center + dv);
        // tmp_mol.translate(dv);

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
  
  Eigen::VectorXf Optimizer_Grid::calcGradient(Molecule &mol, const Eigen::VectorXf& X) {
    
    const fltype &x   = X[0];
    const fltype &y   = X[1];
    const fltype &z   = X[2];
    const fltype &phi = X[3];
    const fltype &the = X[4];
    const fltype &psi = X[5];
    
    Molecule tmp_mol = mol;
    Vector3d center = tmp_mol.getCenter();
    tmp_mol.translate(-center);
    tmp_mol.rotate(phi, the, psi); 
    tmp_mol.translate(Vector3d(x, y, z));
    
    const fltype cosH = cos(phi), sinH = sin(phi), cosT = cos(the), sinT = sin(the), cosS = cos(psi), sinS = sin(psi);
        Eigen::Matrix3f R(3, 3);
        R << cosS*cosH - cosT*sinS*sinH,  cosT*sinS*cosH + cosS*sinH, sinT*sinS,
             -cosT*cosS*sinH - sinS*cosH, cosT*cosS*cosH - sinS*sinH, sinT*cosS,    
                     sinT*sinH,                   -sinT*cosH,            cosT;
        
    Eigen::Matrix3f R_t = R.transpose();
    
    int dim = 6;
    
    Eigen::VectorXf grad = Eigen::VectorXf::Zero(dim);
    
    Vector3d center0 = mol.getCenter();
    
    bool out_of_grid = false;
    fltype energy = 0.0;
    
    for (int i = 0; i < mol.getAtoms().size(); i++) {
      
      auto &atoms = tmp_mol.getAtoms();
      auto &atoms0 = mol.getAtoms();
      auto &a = atoms[i];
      auto &a0 = atoms0[i];
      
      Vector3d Xi_dash = a0 - center0;
      fltype xi_dash = Xi_dash.x;
      fltype yi_dash = Xi_dash.y;
      fltype zi_dash = Xi_dash.z;
      
      Eigen::MatrixXf Qi_t(3, 3);
      Qi_t << zi_dash*sinT*cosS - yi_dash*cosT, xi_dash*cosT - zi_dash*sinT*sinS, yi_dash*sinT*sinS - xi_dash*sinT*cosS,
                        -zi_dash*sinS,                   -zi_dash*cosS,                xi_dash*sinS + yi_dash*cosS,
                          -yi_dash,                         xi_dash,                                0;
          
      const AtomEnergyGrid& ex = gradient_grids[a.getXSType()][0];
      const AtomEnergyGrid& ey = gradient_grids[a.getXSType()][1];
      const AtomEnergyGrid& ez = gradient_grids[a.getXSType()][2];
      
      fltype del_xi = LIMIT_ENERGY;
      fltype del_yi = LIMIT_ENERGY;
      fltype del_zi = LIMIT_ENERGY;
      
      //cout << del_xi << " " << del_yi << " " << del_zi << endl;
      
      int ax = ex.convertX(a);
      int ay = ey.convertY(a);
      int az = ez.convertZ(a);
      
      Point3d<int> num = gradient_grids[a.getXSType()][0].getNum();
      if (ax < 0 or ay < 0 or az < 0 or ax >= num.x or ay >= num.y or az >= num.z) {
        out_of_grid = true;
        
        int tx = ax < 0 ? 0 : ax >= num.x ? num.x - 1 : ax;
        int ty = ay < 0 ? 0 : ay >= num.y ? num.y - 1 : ay;
        int tz = az < 0 ? 0 : az >= num.z ? num.z - 1 : az;
        
        energy += atom_grids[a.getXSType()].getEnergy(tx, ty, tz);
        
        
        del_xi = 0;
        del_yi = 0;
        del_zi = 0;
            
      } else {
        del_xi = ex.getEnergy(a);
        del_yi = ey.getEnergy(a);
        del_zi = ez.getEnergy(a);
      }
          
      grad[0] += del_xi;
      grad[1] += del_yi;
      grad[2] += del_zi;
      
      fltype del_xi_dash = R_t(0, 0)*del_xi + R_t(1, 0)*del_yi + R_t(2, 0)*del_zi;
      fltype del_yi_dash = R_t(0, 1)*del_xi + R_t(1, 1)*del_yi + R_t(2, 1)*del_zi;
      fltype del_zi_dash = R_t(0, 2)*del_xi + R_t(1, 2)*del_yi + R_t(2, 2)*del_zi;
      
      grad[3] += Qi_t(0, 0)*del_xi_dash + Qi_t(0, 1)*del_yi_dash + Qi_t(0, 2)*del_zi_dash;
      grad[4] += Qi_t(1, 0)*del_xi_dash + Qi_t(1, 1)*del_yi_dash + Qi_t(1, 2)*del_zi_dash;
      grad[5] += Qi_t(2, 0)*del_xi_dash + Qi_t(2, 1)*del_yi_dash + Qi_t(2, 2)*del_zi_dash;
          
    }
    
    
    if (out_of_grid) {
      for(int i = 0; i < dim; i++) grad[i] = 0;
    }
    
    return grad;
  }



  fltype Optimizer_Grid::optimizeBySteepest(Molecule &mol, int max_iterations, fltype step_size) {
    using namespace std;
    
    Vector3d center = mol.getCenter();
    fltype opt = calcscore(mol);
    int dim = 6;
    Eigen::VectorXf X(dim); X << center.x, center.y, center.z, 0, 0, 0;
    Eigen::VectorXf preX(dim); preX << 1e9, 1e9, 1e9, 1e9, 1e9, 1e9;

    int iter;
    const fltype pi = acos(-1.0);
    fltype lim_trans_step = 0.5;
    fltype lim_rotate_step = pi / 1800.0;
    
    Eigen::VectorXf minX(dim); minX << center.x, center.y, center.z, 0, 0, 0;
    fltype min = calcscore(mol);
    for (iter = 0; iter < max_iterations; iter++) {
      Eigen::VectorXf grad = calcGradient(mol, X);
      
      
      Eigen::VectorXf dX = step_size * grad;
      
      float max_trans_step = max({abs(dX[0]), abs(dX[1]), abs(dX[2])});
      if (max_trans_step > lim_trans_step) {
        for (int j = 0; j < 3; j++)
          dX[j] *= lim_trans_step / max_trans_step;
      }
      
      float max_rotate_step = max({abs(dX[3]), abs(dX[4]), abs(dX[5])});
      if (max_rotate_step > lim_rotate_step) {
        for (int j = 3; j < 6; j++)
          dX[j] *= lim_rotate_step / max_rotate_step;
      }
      
      //dX[0]=dX[1]=dX[2]=0;//x,y,z完全固定
      dX[4]=dX[5]=0;//the,psi完全固定
      X -= dX;
      
      const fltype x   = X[0];
      const fltype y   = X[1];
      const fltype z   = X[2];
      const fltype phi = X[3];
      const fltype the = X[4];
      const fltype psi = X[5];
      
      Molecule tmpmol = mol;
      tmpmol.translate(-center);
      tmpmol.rotate(phi, the, psi);
      tmpmol.translate(Vector3d(x, y, z));
      fltype tmpscore = calcscore(tmpmol);
      
      cout<<"grad "<<iter<<':';
      for(int i=0;i<dim;i++) cout<<dX[i]<<",";
      cout<<"score:"<<tmpscore<<'\n';
      
      
      if (tmpscore < min) {
        min = tmpscore;
        minX = X;
      }


      const fltype EPS = 1e-2;
      if ((X-preX).norm() < EPS) {
        cout<<'\n';
        cout<<"preX:";
        for (int j=0;j<dim;j++) printf("%.6f,", preX[j]);cout<<'\n';
        cout<<"X:";
        for (int j=0;j<dim;j++) printf("%.6f,", X[j]);cout<<'\n';
        break;
      } 
      
      preX = X;
    }
    cout<<'\n';

    Eigen::VectorXf oriX(dim); oriX << center.x, center.y, center.z, 0, 0, 0;
    cout<<"orig:";
    for (int j=0;j<dim;j++) printf("%.6f,", oriX[j]);cout<<'\n';
    cout<<"last:";
    for (int j=0;j<dim;j++) printf("%.6f,", X[j]);cout<<'\n';
    cout<<"min :";
    for (int j=0;j<dim;j++) printf("%.6f,", minX[j]);cout<<"\n\n";
    cout<<"opt :"<<opt<<"\n\n";
    
    const fltype x   = minX[0];
    const fltype y   = minX[1];
    const fltype z   = minX[2];
    const fltype phi = minX[3];
    const fltype the = minX[4];
    const fltype psi = minX[5];
    
    mol.translate(-center);
    mol.rotate(phi, the, psi); 
    mol.translate(Vector3d(x, y, z));
    
    opt = calcscore(mol);
    
    return opt;
  }



  fltype Optimizer_Grid::optimize(Molecule &mol) const {
    using namespace std;

    // hill climbing
    fltype opt = calcscore(mol);//cout<<opt<<" ";
    const fltype pi = acos(-1.0);
    const int NEAREST_NUM = 200;
    fltype trans_step = 0.5;
    fltype rotate_step = pi / 30.0;
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
        tmp_mol.translate(-center); //rotateは原点中心の回転なので目的リガンドを原点に一度持ってきて回転させると楽ということ．
        tmp_mol.rotate(th, phi, psi); 
        tmp_mol.translate(center + dv);
        // tmp_mol.translate(dv);

        fltype val = calcscore(tmp_mol);
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
//cout<<opt<<"\n";
    return opt;
  }
  
}
