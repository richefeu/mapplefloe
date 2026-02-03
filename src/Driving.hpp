#pragma once

#include <iostream>
#include <cmath>

// Mother class for any lodaing condition
class Driving {
public:
  double vx{0.0};
  double vy{0.0};
  double vz{0.0};
  double vrot{0.0};
  
  static Driving *create(const std::string &token);

  virtual void read(std::istream &is) = 0;
  virtual void write(std::ostream &os) = 0;
  virtual void set(double) {}

  //Drive() = delete; // deactivated Ctor
  virtual ~Driving(); // virtual Dtor
};

class imposedVelocities : public Driving {
public:
  
  imposedVelocities();

  virtual void read(std::istream &is);
  virtual void write(std::ostream &os);
  virtual void set(double);
};

class imposedVelocityCycles : public Driving {
public:
  
  double T1{0.0};
  double T2{0.0};
  double T3{0.0};
  
  double vx1{0.0};
  double vx2{0.0};
  double vx3{0.0};
  
  double vy1{0.0};
  double vy2{0.0};
  double vy3{0.0};
  
  double vz1{0.0};
  double vz2{0.0};
  double vz3{0.0};
  
  double vrot1{0.0};
  double vrot2{0.0};
  double vrot3{0.0};
  
  imposedVelocityCycles();

  virtual void read(std::istream &is);
  virtual void write(std::ostream &os);
  virtual void set(double);
};


