#pragma once

#ifndef VECTORMATH_H
#define VECTORMATH_H
#include "systems.hpp"
#ifdef SYSTEM_IS_SET

NAMESPACESTART
class Vec2D : public CoordinateSystem2D {
private:
  float abs;
  void calcAbs();

public:
  // Basic operrator overloading:

  Vec2D operator+(Vec2D &Mathvector);
  Vec2D operator-(Vec2D &Mathvector);
  float operator*(Vec2D &Mathvector);

#ifdef CARTESIAN_IS_2D_STANDARD
  Vec2D(float X, float Y);
  Vec2D(float X, float Y, float originX, float originY);
#endif

#ifdef POLAR_IS_STANDARD
  Vec2D(float radius, float phi);
#endif

  // Math with scalar results
  float getAbs();                     // Calculates the absolutvalue
  float lenght();                     // Calculates the length of a vector
  float dotProduct(Vec2D Mathvector); // Calculates Scalar product
  float calcAngle(Vec2D Mathvector);  // Calculates the angle between 2 vectors

  float polarSystemArea(); // Calculates the circle of a polar system
  float
  polarSystemCircumfrance(); // Calculates the circumfrance of a polar system
  float rectangleArea(Vec2D MathVector); // Calculates the area of a rectangle
                                         // with length of vector 1 and vector 2
  float rectangleCircumfrance(
      Vec2D MathVector); // Calculates the circumfrance of a rectangle with
                         // length of vector 1 and vector 2

  /* Setters:
   * If data is set/changed it will automaticly change the values for other
   * coordinate types as well.
   */

  virtual void setX(float value);
  virtual void setY(float value);
  virtual void setPhi(float value);
  virtual void setRadius(float radius);

  // You can use methods insteadof operators
  Vec2D add(Vec2D Mathvector);      // Adding 2 vectors
  Vec2D subtract(Vec2D Mathvector); // Subtract 2 vectors
  Vec2D rotateVector(float angle);  // Rotates a vector dependent op the angle

  // spezial Vectors
  static Vec2D polarVector(float radius,
                           float angle); // Creates a Vector in a polar system
  static Vec2D zeroVector();             // Creates a vector with lenght 0
};

class Vec3D : public CoordinateSystem3D {
private:
  float abs;
  void calcAbs();

public: // Operator overloading
  Vec3D operator+(Vec3D &Mathvector);
  Vec3D operator-(Vec3D &Mathvector);
  float operator*(Vec3D &Mathvector);

#ifdef CARTESIAN_IS_3D_STANDARD
  Vec3D(float X, float Y, float Z);
  Vec3D(float X, float Y, float Z, float originX, float originY);
#endif

#ifdef SPHERE_IS_STANDARD
  Vec3D(float radius, float phi, float theta);
#endif

#ifdef CYLINDER_IS_STANDARD
  Vec3D(float radius, float phi, float height);
#endif

  // Setters
  void setX(float X);
  void setY(float Y);
  void setZ(float Z);

  void setPhi(float phi); // For spheresystems (radius, angle Phi, angle Theta)
  void setHeight(float phi); // For cylinder systems (radius, angle, height);
  void
  setTheta(float phi); // For spheresystems (radius, angle Phi, angle Theta);

  // Calculations that return a scalar value
  float dotProduct(Vec3D Mathvector);
  float calcAngle(Vec3D Mathvector);

  // Calculations that return a vector
  Vec3D vecProduct(Vec3D Mathvector);
  Vec3D add(Vec3D Mathvector);
  Vec3D subtract(Vec3D Mathvector);

  // Return a spezial vector
  static Vec3D zeroVector(); // Returns a Vector with lenght 0
  static Vec3D sphereVector(float radius, float angleOne, float angleTwo);
  static Vec3D cylinderVector(float radius, float angle, float height);

  // Volumes and surfaces of sphere and cylinder systems
  float cylinderSystemVolume();
  float cylinderSystemSurface();
  float cylinderSystemLateralSurface();

  float sphereSystemVolume();
  float sphereSystemSurface();

  // Volumes and surface from a cuboid build with vectors
  float cuboidVolume(Vec3D MathVector, Vec3D MathVectorTwo);
  float cuboidSurface(Vec3D MathVector, Vec3D MathVectorTwo);

  float getAbs();
  float getSphereRadius();
  float getCylinderRadius();
};

#endif // SYSTEM_IS_SET
NAMESPACEEND
#endif // VECTORMATH_H