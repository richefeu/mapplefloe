#include "Driving.hpp"

// Remark: if the creation fails, nullptr is returned
Driving *Driving::create(const std::string &token) {
  if (token == "imposedVelocities") {
    return new imposedVelocities();
  }  else if (token == "imposedVelocityCycles") {
    return new imposedVelocityCycles();
  }

  return nullptr;
}

Driving::~Driving() {}

// ========== imposedVelocities ==========

imposedVelocities::imposedVelocities() {}

void imposedVelocities::read(std::istream &is) {
  is >> vx >> vy >> vz >> vrot;
}

void imposedVelocities::write(std::ostream &os) {
  os << "imposedVelocities " << vx << ' ' << vy << ' ' << vz << ' ' << vrot << std::endl;
}

void imposedVelocities::set(double) {}

// ========== imposedVelocityCycles ==========

imposedVelocityCycles::imposedVelocityCycles() {}

void imposedVelocityCycles::read(std::istream &is) {
  is >> T1 >> vx1 >> vy1 >> vz1 >> vrot1 >> T2 >> vx2 >> vy2 >> vz2 >> vrot2 >> T3 >> vx3 >> vy3 >> vz3 >> vrot3;
}

void imposedVelocityCycles::write(std::ostream &os) {
  os << "imposedVelocityCycles" << std::endl
     << T1 << ' ' << vx1 << ' ' << vy1 << ' ' << vz1 << ' ' << vrot1 << std::endl
     << T2 << ' ' << vx2 << ' ' << vy2 << ' ' << vz2 << ' ' << vrot2 << std::endl
     << T3 << ' ' << vx3 << ' ' << vy3 << ' ' << vz3 << ' ' << vrot3 << std::endl;
}

void imposedVelocityCycles::set(double t) {
  double sumT = T1 + T2 + T3;
  double ct = std::fmod(t, sumT);
  
  if (ct < T1) {
    vx = vx1;
    vy = vy1;
    vz = vz1;
    vrot = vrot1;    
  } else if (ct > T1+T2) {
    vx = vx3;
    vy = vy3;
    vz = vz3;
    vrot = vrot3; 
  } else {
    vx = vx2;
    vy = vy2;
    vz = vz2;
    vrot = vrot2;
  }
  
  
  
}
