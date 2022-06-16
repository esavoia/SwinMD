/** Matrix3.h -- definition and operations for a 3 by 3 matrix
 **
 ** [ m00 m01 m02 ]   
 ** [ m10 m11 m12 ]   representation of a 3 by 3 matrix
 ** [ m20 m21 m22 ]   
 **
 ** referenc: PROTOMOL 
 **/

#ifndef MATRIX3_H
#define MATRIX3_H

#include <iostream>
#include <math.h>
#include <stdio.h>
#include "../NEMD_defs.h"
#include "Vector3.h"

using std::endl;
using std::ostream;


class Matrix3 {

/**
 ** My data members
 **/
public:
  Double xx, xy, xz;
  Double yx, yy, yz;
  Double zx, zy, zz;

private:
  static const Matrix3 unitMatrix;

/**
 ** Constructors, destructors
 **/
public:
  Matrix3();
 
  Matrix3(Double x00,   Double x01,   Double x02,
      Double x10,   Double x11,   Double x12,
      Double x20,   Double x21,   Double x22);

  Matrix3(Double mat[9]);

  Matrix3(const Matrix3 &tm);

  // construct matrix by cross product of two Vector3
  Matrix3(const Vector3 &a, const Vector3 &b);

  // construct the matrix with three Vector3 as rows.
  Matrix3(const Vector3 &v1, const Vector3 &v2, const Vector3 &v3);
  
public:
  // Matrix3& operator=(const Matrix3&);  // Use default version.

public:
  // ~Matrix3();  // Use default version.
  
/**
 ** New methods of class Matrix3
 **/
public:
  const Matrix3& get_unit_matrix();  // Return the unit matrix.
  void unit_matrix();  // Set 'this' to the unit matrix.
  void zero_matrix();// Set 'this' to the zero matrix.
  void set_matrix(Double x00,   Double x01,   Double x02,
           Double x10,   Double x11,   Double x12,
           Double x20,   Double x21,   Double x22);
  bool is_zero() const; // Check if all elements equals zero.
  Double det();

  // access an element through index -- M(i,j) => value
  Double operator()(int i, int j) const;
  // Set an element.
  void operator()(int i, int j, Double x);
  
  int operator==(const Matrix3& tm) const;
  int operator!=(const Matrix3& tm) const;
 
  // Matrix3& operator*=(const Matrix3& tm);
  void operator*=(const Matrix3& tm);

  // Matrix multiplication.
  Matrix3 operator*(const Matrix3& tm) const;  

  // matrix multiply a Vector, return a Vector
  Vector3 operator*(const Vector3& tm) const;
  
  // Matrix multiply constant. 
  void operator*=(const Double tm);
  Matrix3 operator*(const Double tm) const;  
  
  int operator/=(const Double tm);
  Matrix3 operator/(const Double tm) const;  
  // Matrix div constant.

  void operator+=(const Matrix3& tm);
  Matrix3 operator+(const Matrix3& tm) const;  
  // Matrix addition.

  void operator-=(const Matrix3& tm);
  Matrix3 operator-(const Matrix3& tm) const;  
  // Matrix subtraction.

  Matrix3 operator-(void) const;  
  // negate operator that returns (- (*this)).

  // Matrix transpose
  void transpose_me(); // transpose 'this'. 
  void transpose(const Matrix3& tm); // Set 'this' to tm transposed, leave tm unchanged.
  // Return the transposed matrix, leave the original unchanged.
  Matrix3 transposed() const; 
  
  // Compute the inverse of the given matrix. 
  // Returns: false -> singular matrix
  //          true  -> success.
  bool invert();      
  
  // Concat scale matrix to the current transformation.
  void scale(Double sx, Double sy, Double sz);
  void scale(Double s);
  void scale(const Vector3& scaleFactor);
  
  void convert(const Matrix3& from, Double to[9]);
  
  // need to check correctness of this function
  friend Vector3 operator*(const Matrix3& tm, const Vector3& point)
  {
    return Vector3(tm.xx*point.x + tm.yx*point.y + tm.zx*point.z,
        tm.xy*point.x + tm.yy*point.y + tm.zy*point.z,
        tm.xz*point.x + tm.yz*point.y + tm.zz*point.z);
  }

  friend ostream& operator<<(ostream &os, const Matrix3 &tm)
  {
    os << "[[" << tm.xx << ", " << tm.xy << ", " <<tm.xz << "], ";
    os <<  "[" << tm.yx << ", " << tm.yy << ", " <<tm.yz << "], ";
    os <<  "[" << tm.zx << ", " << tm.zy << ", " <<tm.zz << "]]";
    return os;
  }
  
};


#endif // MATRIX3_H
