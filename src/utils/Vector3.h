/**  Vector3.h - definition of vector and vectorised operations for position,
 **     velocity, force ... 'inline' functions used for better performance.
 **
 **  Copyright (c) 2002 by
 **  Centre for Molecular Simulation, School of Information Technology
 **  Swinburne University of Technology
 **  
 **  Author: Zhongwu Zhou
 **
 **/

#ifndef VECTOR3_H
#define VECTOR3_H

#include <iostream>
#include <math.h>
#include <stdio.h>
#include "../NEMD_defs.h"
using namespace std;

class Vector3{
   public:
     Double x,y,z;
     
	 /** 
	  ** Constructors 
	  **/
     inline Vector3(void) : x(0.0), y(0.0), z(0.0) { ; }

     inline Vector3( Double newx, Double newy, Double newz){
            x = newx;
	        y = newy;
	        z = newz;
	 }

	 // construct Vector3 from another copy 
     inline Vector3(const Vector3 &v) { 
		    x = v.x;
	        y = v.y;
	        z = v.z;
	 }

	 /**
	  ** operators
	  **/
	 // Index access
     inline Double operator[](int i) const {
       if (i == 0)        return  x;
       else if (i == 1)   return  y;
       else if (i == 2)   return  z;
       // (Error("Vector reference out of bounds.");
       return 0.0;
     }

    // V = v1
    inline Vector3& operator=(const Vector3 &v1)
    {
        x = v1.x;  y = v1.y;  z = v1.z;
        return *this;
    }

     //  V = value 
     inline Vector3& operator=(const Double &value) {
       x = value;
       y = value;
       z = value;
       return *this;
     }

     //  V += v1 
     inline void operator+=(const Vector3 &v1) {
       x += v1.x;
       y += v1.y;
       z += v1.z;
     }

     // V -= v1
     inline void operator-=(const Vector3 &v1) {
       x -= v1.x;
       y -= v1.y;
       z -= v1.z;
     }

     // V *= const
     inline void operator*=(const Double &value) {
       x *= value;
       y *= value;
       z *= value;
     }

     // V /= const
     inline void operator/=(const Double& value) {
       x /= value;
       y /= value;
       z /= value;
     }

	 // comparing the equility of two double vectors could be of some problem in practive use
	 // we should avoid to use following two operators 
     inline friend int operator == (const Vector3& v1, const Vector3& v2) {
       return (v1.x == v2.x && v1.y == v2.y && v1.z == v2.z);
     }
     inline friend int operator != (const Vector3& v1, const Vector3& v2) {
       return (v1.x != v2.x || v1.y != v2.y || v1.z != v2.z);
     }

     // addition of two vectors -- V3 = v1 + v2
     inline friend Vector3 operator+(const Vector3& v1, const Vector3& v2) {
       return Vector3( v1.x+v2.x, v1.y+v2.y, v1.z+v2.z);
     }

     // negation -- V3 = -v1
     inline friend Vector3 operator-(const Vector3 &v1) {
       return Vector3( -v1.x, -v1.y, -v1.z);
     }

     // subtraction -- V3 = v1 - v2
     inline friend Vector3 operator-(const Vector3 &v1, const Vector3 &v2) {
       return Vector3( v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);
     }
     // inner ("dot") product -- value = v1*v2
     inline friend Double operator*(const Vector3 &v1, const Vector3 &v2) {
       return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
     }
     // scalar product -- V3 = f*v1
     inline friend Vector3 operator*(const Double &f, const Vector3 &v1) {
       return Vector3(f*v1.x, f*v1.y, f*v1.z);
     }
     // scalar product -- V3 = v1*f
     inline friend Vector3 operator*(const Vector3 &v1, const Double &f) {
       return Vector3(f*v1.x, f*v1.y, f*v1.z);
     }
     // division by a scalar -- V3 = v1/f
     inline friend Vector3 operator/(const Vector3 &v1, const Double &f) {
//       if (!f)
//         Error_exit("Division by 0 on a vector operation.");
       return Vector3(v1.x/f, v1.y/f, v1.z/f);
     }
     
     /**
	  ** Vector operations
	  **/
     // add a vector to this vector -- V.add(v1), equl. V+=v1
     inline void add(const Vector3 &v1) {
       x+=v1.x; 
	   y+=v1.y; 
	   z+=v1.z;
     }

	 // add a constant factor to each element of a vector -- V.add(c) => V
     // void add_const(Double c)
	 inline void add(Double c)
     {
	   x+=c;
	   y+=c;
	   z+=c;
     }

     // subtract the vector from this one -- V.sub(v1), equl. V-=v1
     inline void sub(const Vector3 &v1) {
       x-=v1.x; 
	   y-=v1.y; 
	   z-=v1.z;
     }

	 // rescale everything by a scalar --  V.mult(f), equl. V*=f
     inline void mult(Double f) {
       x*=f; 
	   y*=f; 
	   z*=f;
     }

	 // value = V.dot(v2), equl. value = V*v2
     inline Double dot(const Vector3 &v2) {
       return x*v2.x + y*v2.y + z*v2.z;
     }

     // divide each element by a scalar -- V.div(f), equl. V/=f
     void div(Double f) {
       x/=f; 
	   y/=f; 
	   z/=f;
     }

     // distance from the origin(0, 0, 0)
     inline Double length(void) const {
       return sqrt(x*x+y*y+z*z);
     }
     
     inline Double length2(void) const {
       return (x*x + y*y + z*z);
     }

     // return the unit vector in the same direction
     inline Vector3 unitVector(void) const {
	   Double len = length();
       
	   if (len == 0.0)
		   return Vector3(0.0, 0.0, 0.0);

       return Vector3(x/len, y/len, z/len);
     }
     
     
     // one cross product  v3 = cross(v1, v2)
     inline friend Vector3 cross(const Vector3 &v1, const Vector3 &v2) {
       return Vector3( v1.y*v2.z-v2.y*v1.z,
                      v2.x*v1.z-v1.x*v2.z,
                      v1.x*v2.y-v2.x*v1.y  );
     }

     // multiplying a cross product by a scalar is very common
     // one cross product  v3 = k*cross(v1, v2)
     inline friend Vector3 cross(const Double &k, const Vector3 &v1, const Vector3 &v2) {
       return Vector3( k*(v1.y*v2.z-v2.y*v1.z),
                      k*(v2.x*v1.z-v1.x*v2.z),
                      k*(v1.x*v2.y-v2.x*v1.y) );
     }

     // cross product between matrices  v3 = cross(double[][], double[][],vones)
     inline friend Vector3 cross(const Double (&k)[3][3], const Double (&m)[3][3], const Vector3 &v1) {
       return Vector3( k[1][0]*m[2][0] - k[2][0]*m[1][0] + k[1][1]*m[2][1] - k[2][1]*m[1][1] + k[1][2]*m[2][2] - k[2][2]*m[1][2] ,
                       k[2][0]*m[0][0] - k[0][0]*m[2][0] + k[2][1]*m[0][1] - k[0][1]*m[2][1] + k[2][2]*m[0][2] - k[0][2]*m[2][2] ,
                       k[0][0]*m[1][0] - k[1][0]*m[0][0] + k[0][1]*m[1][1] - k[1][1]*m[0][1] + k[0][2]*m[1][2] - k[1][2]*m[0][2] );
     }

     friend ostream &operator<<(ostream &output , const Vector3 &v1)
     {
         output << v1.x << '\t' << v1.y << '\t' << v1.z << '\n';
         return output;
     }
	   
};


#endif
