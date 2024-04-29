#pragma once

#ifndef VECTORMATH_H
#define VECTORMATH_H
#include "systems.hpp"
#ifdef SYSTEM_IS_SET

NAMESPACESTART
class Vec2D : public CoordinateSystem2D {
  private:
    double aX;                                               // Absolute X the coordinate of a vector that does not start at zero
    double aY;                                               // Absolute Y the coordinate of a vector that does not start at zero
    double abs;                                              // The absolute value or lenght of a vector
    double distanceToZero;   
    double vectorRotation;                                
    void calcDZ();                                          // Calculates the distance to zero
    void calcAbs();                                         // Calculates the absolute value
    void calcAbsXY();                                       // Calculates and sets the absolute x and y coordinates
    

  public:
    // Basic operrator overloading:
    Vec2D operator/(double scalarValue);                     //divides all components of the Vector with the scalar value and returns it
    Vec2D operator+(Vec2D& Mathvector);                     //Basic Vector addition
    Vec2D operator-(Vec2D& Mathvector);                     //Basic Vector subtractions
    double operator*(Vec2D& Mathvector);                     //Calculates the dotproduct of two vectors
    
#ifdef CARTESIAN_IS_2D_STANDARD
    Vec2D(double XY);                                        // creates a Vector with the same x and y value
    Vec2D(double X, double Y);                                // basic x and y vector
    Vec2D(double X, double Y, double originX, double originY);  // crates a vector with a specific non zero startpoint
#endif // CARTESIAN_IS_2D_STANDARD

#ifdef POLAR_IS_STANDARD
    Vec2D(double radius, double phi);                          // Creates a Vector in a polar coordinatesystem
#endif // POLAR_IS_STANDARD

    // Math with scalar results

    double getAX();                                       // Gets the absolute X coordinate from a vector with a non zero start
    double getAY();                                       // Gets the absolute X coordinate from a vector with a non zero start
    double getDistanceToZero();                           // Gets the distance to [ 0 ; 0 ] from a vector with a non zero start
    double getRotationAngle();                            // Gets the rotation angle for rotated vectors
    double getAbs();                                      // Calculates the absolutvalue
    double lenght();                                      // Calculates the length of a vector

    double dotProduct(Vec2D Mathvector);                  // Calculates Scalar product
    double calcAngle(Vec2D Mathvector);                   // Calculates the angle between 2 vectors
    double distance(Vec2D Mathvector);                    // Calculates the distance between vectors
                    
    double polarSystemArea();                             // Calculates the circle of a polar system
    double polarSystemCircumfrance();                     // Calculates the circumfrance of a polar system
    double rectangleArea(Vec2D MathVector);               // Calculates the area of a rectangle
    double rectangleCircumfrance(Vec2D MathVector);       // Calculates the circumfrance of the rectangle
    
    /* Setters:
     * If data is set/changed it will automaticly change the values for other
     * coordinate types as well.
     */

    void setX(double value);                              // Basic x setter
    void setY(double value);                              // Basic y setter
    void setPhi(double value);                            // Setter for the angle in a polar system
    void setRadius(double radius);                        // Setter for the radius in a polar system

    void normalize();

    void setOriginX(double newOriginX);
    void setOriginY(double newOriginY);
   
    void moveVectorX(double move);                        // Moves the vector in a specific x-direction
    void moveVectorY(double move);                        // Moves the vector in a specific y-direction
    void moveVector(double moveX, double moveY);           // Moves the vector in a specific x and y direction 
    void rotateThisVector(double value);

    // You can use methods instead of operators
    Vec2D add(Vec2D Mathvector);                         // Adding 2 vectors
    Vec2D subtract(Vec2D Mathvector);                    // Subtract 2 vectors
    Vec2D rotateVector(double radians);                   // Rotates a vector dependent op the angle

    // spezial Vectors
    static Vec2D polarVector(double radius, double angle); // Creates a Vector in a polar system
    static Vec2D zeroVector();                           // Creates a vector with lenght 0
};


//Most of the methods in the Vec3D class ar simular to the methods in the Vec2D class
class Vec3D : public CoordinateSystem3D {
  private:
    double aX;
    double aY;
    double aZ;
    double distanceToZero;
    double abs;
    void calcAbs();
    void calcDTZ();
    void calcAXYZ();

  public: // Operator overloading
    Vec3D operator+(Vec3D& Mathvector);
    Vec3D operator-(Vec3D& Mathvector);
    Vec3D operator/(double scalarValue);
    double operator*(Vec3D& Mathvector);


#ifdef CARTESIAN_IS_3D_STANDARD
    Vec3D(double XYZ);
    Vec3D(double X, double Y, double Z);
    Vec3D(double X, double Y, double Z, double originX, double originY, double originZ );
#endif

#ifdef SPHERE_IS_STANDARD
    Vec3D(double radius, double phi, double theta);
#endif

#ifdef CYLINDER_IS_STANDARD
    Vec3D(double radius, double phi, double height);
#endif

    // Setters
    void setX(double X);
    void setY(double Y);
    void setZ(double Z);

    void setOriginX(double value);
    void setOriginY(double value);
    void setOriginZ(double value);

    void normalize();


    //SetOrigin!!!!!!!

    void setPhi(double phi);       // For spheresystems (radius, angle Phi, angle Theta)
    void setHeight(double phi);    // For cylinder systems (radius, angle, height);
    void setTheta(double phi);     // For spheresystems (radius, angle Phi, angle Theta);

    void moveVectorX(double move); // Moves the vector in a specific x-direction
    void moveVectorY(double move);  
    void moveVectorZ(double Z);        
    void moveVector(double moveX, double moveY, double moveZ);           
    void rotateThisVector(double phi, double Theta);

    // Calculations that return a scalar value
    double dotProduct(Vec3D Mathvector);
    double calcAngle(Vec3D Mathvector);

    // Calculations that return a vector
    Vec3D vecProduct(Vec3D Mathvector);
    Vec3D add(Vec3D Mathvector);
    Vec3D subtract(Vec3D Mathvector);

    Dmath::Vec3D rotateVector(double radiansPhi, double radiansTheta);

    // Return a spezial vector
    static Vec3D zeroVector(); // Returns a Vector with lenght 0
    static Vec3D sphereVector(double radius, double angleOne, double angleTwo);
    static Vec3D cylinderVector(double radius, double angle, double height);

    // Volumes and surfaces of sphere and cylinder systems
    double cylinderSystemVolume();
    double cylinderSystemSurface();
    double cylinderSystemLateralSurface();

    double sphereSystemVolume();
    double sphereSystemSurface();

    // Volumes and surface from a cuboid build with vectors
    double cuboidVolume(Vec3D MathVector, Vec3D MathVectorTwo);
    double cuboidSurface(Vec3D MathVector, Vec3D MathVectorTwo);

    double getAbs();
    double getSphereRadius();
    double getCylinderRadius();
};

#endif // SYSTEM_IS_SET
NAMESPACEEND
#endif // VECTORMATH_H