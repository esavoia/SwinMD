/** Matrix3.cpp
 **/

#include "Matrix3.h"


const Matrix3 Matrix3::unitMatrix = Matrix3(1,0,0, 0,1,0, 0,0,1);

/**
 ** Constructors
 **/
Matrix3::Matrix3()
{ 
  zero_matrix();
}

Matrix3::Matrix3(Double x00, Double x01, Double x02,
		       Double x10, Double x11, Double x12,
		       Double x20, Double x21, Double x22)
{
  set_matrix(x00, x01, x02, x10, x11, x12, x20, x21, x22);
}
  
Matrix3::Matrix3(Double mat[9])
{
  xx = *mat++; xy = *mat++; xz = *mat++;
  yx = *mat++; yy = *mat++; yz = *mat++;
  zx = *mat++; zy = *mat++; zz = *mat;
}

Matrix3::Matrix3(const Matrix3 &tm)
{
  xx = tm.xx; xy = tm.xy; xz = tm.xz;
  yx = tm.yx; yy = tm.yy; yz = tm.yz;
  zx = tm.zx; zy = tm.zy; zz = tm.zz;
}

//construct a matrix from the cross product of two Vectors
Matrix3::Matrix3(const Vector3 &a, const Vector3 &b)
{
  set_matrix(a.x*b.x, a.x*b.y, a.x*b.z,
             a.y*b.x, a.y*b.y, a.y*b.z,
             a.z*b.x, a.z*b.y, a.z*b.z);
}

Matrix3::Matrix3(const Vector3 &v1, const Vector3 &v2, 
		       const Vector3 &v3)
{
  set_matrix(v1.x,v1.y,v1.z, v2.x,v2.y,v2.z, v3.x,v3.y,v3.z);
}

/**
 ** Functions
 **/
void Matrix3::unit_matrix() 
{
  set_matrix(1,0,0, 0,1,0, 0,0,1);
}

void Matrix3::zero_matrix() 
{
  set_matrix(0,0,0, 0,0,0, 0,0,0); 
}

bool Matrix3::is_zero() const
{
  return *this == Matrix3(0,0,0, 0,0,0, 0,0,0);
}

const Matrix3& Matrix3::get_unit_matrix()
{
  return unitMatrix;
}

// access an element through index
Double Matrix3::operator()(int i, int j) const
{
  if (i == 0) {
    if (j == 0)       return xx;
    else if (j == 1)  return xy;
    else if (j == 2)  return xz;
  }
  else if (i == 1) {
    if (j == 0)       return yx;
    else if (j == 1)  return yy;
    else if (j == 2)  return yz;
  }
  else if (i == 2) {
    if (j == 0)       return zx;
    else if (j == 1)  return zy;
    else if (j == 2)  return zz;
  }
  // Error("[Matrix3::operator()] index out of range.");
  return 0;
}

// set [i][j] element to a value
void Matrix3::operator()(int i, int j, Double value)
{
  if (i == 0) {
    if (j == 0)      { xx = value; return;}
    else if (j == 1) { xy = value; return;}
    else if (j == 2) { xz = value; return;}
  }
  else if (i == 1) {
    if (j == 0)      { yx = value; return;}
    else if (j == 1) { yy = value; return;}
    else if (j == 2) { yz = value; return;}
  }
  else if (i == 2) {
    if (j == 0)      { zx = value; return;}
    else if (j == 1) { zy = value; return;}
    else if (j == 2) { zz = value; return;}
  }
  // Error("[Matrix3::operator()] index out of range.");
}


void Matrix3::set_matrix(Double x00, Double x01, Double x02,
			    Double x10, Double x11, Double x12,
			    Double x20, Double x21, Double x22)
{
  xx = x00; xy = x01; xz = x02;
  yx = x10; yy = x11; yz = x12;
  zx = x20; zy = x21; zz = x22;
}

int Matrix3::operator==(const Matrix3& tm) const
{
  return (xx==tm.xx && xy==tm.xy && xz==tm.xz &&
          yx==tm.yx && yy==tm.yy && yz==tm.yz &&
	      zx==tm.zx && zy==tm.zy && zz==tm.zz);
}

int Matrix3::operator!=(const Matrix3& tm) const
{
  return !(*this == tm);
}

// M *= value
void Matrix3::operator*=(const Double value)
{
  xx *= value; xy *= value; xz *= value;
  yx *= value; yy *= value; yz *= value;
  zx *= value; zy *= value; zz *= value;
}

// M /= value
int Matrix3::operator/=(const Double value)
{
  if (value == 0) 
  {
    //Error("devide by zero error"); -- leave the caller to handle this
    return -1;
  }
  xx /= value; xy /= value; xz /= value;
  yx /= value; yy /= value; yz /= value;
  zx /= value; zy /= value; zz /= value;
  return 0;
}

// M *= tm
void Matrix3::operator*=(const Matrix3& tm)
{
  Matrix3 temp(*this);

  xx = temp.xx*tm.xx + temp.xy*tm.yx + temp.xz*tm.zx;
  xy = temp.xx*tm.xy + temp.xy*tm.yy + temp.xz*tm.zy;
  xz = temp.xx*tm.xz + temp.xy*tm.yz + temp.xz*tm.zz;
  
  yx = temp.yx*tm.xx + temp.yy*tm.yx + temp.yz*tm.zx;
  yy = temp.yx*tm.xy + temp.yy*tm.yy + temp.yz*tm.zy;
  yz = temp.yx*tm.xz + temp.yy*tm.yz + temp.yz*tm.zz;
  
  zx = temp.zx*tm.xx + temp.zy*tm.yx + temp.zz*tm.zx;
  zy = temp.zx*tm.xy + temp.zy*tm.yy + temp.zz*tm.zy;
  zz = temp.zx*tm.xz + temp.zy*tm.yz + temp.zz*tm.zz;

}

// M * value => newM
Matrix3 Matrix3::operator*(const Double value) const 
{
   return Matrix3(xx*value, xy*value, xz*value, 
                  yx*value, yy*value, yz*value,
                  zx*value, zy*value, zz*value);
}

// M/value => newM
Matrix3 Matrix3::operator/(const Double value) const
{
  return Matrix3(xx/value, xy/value, xz/value, 
                 yx/value, yy/value, yz/value,
                 zx/value, zy/value, zz/value);
}

// M*tm => newM
Matrix3 Matrix3::operator*(const Matrix3& tm) const 
{
  Matrix3 res;
  res.xx = xx*tm.xx +  xy*tm.yx + xz*tm.zx;
  res.xy = xx*tm.xy +  xy*tm.yy + xz*tm.zy;
  res.xz = xx*tm.xz +  xy*tm.yz + xz*tm.zz;

  res.yx = yx*tm.xx +  yy*tm.yx + yz*tm.zx;
  res.yy = yx*tm.xy +  yy*tm.yy + yz*tm.zy;
  res.yz = yx*tm.xz +  yy*tm.yz + yz*tm.zz;

  res.zx = zx*tm.xx +  zy*tm.yx + zz*tm.zx;
  res.zy = zx*tm.xy +  zy*tm.yy + zz*tm.zy;
  res.zz = zx*tm.xz +  zy*tm.yz + zz*tm.zz;
  return res;
}

// M * vector => newVector
Vector3 Matrix3::operator*(const Vector3& tm) const 
{
  Vector3 res;
  res.x = xx*tm.x +  xy*tm.y + xz*tm.z;
  res.y = yx*tm.x +  yy*tm.y + yz*tm.z;
  res.z = zx*tm.x +  zy*tm.y + zz*tm.z;
  return res;
}

// M -= tm
void Matrix3::operator-=(const Matrix3& tm)
{
  xx -= tm.xx; 
  xy -= tm.xy;
  xz -= tm.xz;
  
  yx -= tm.yx;
  yy -= tm.yy;
  yz -= tm.yz;
  
  zx -= tm.zx;
  zy -= tm.zy;
  zz -= tm.zz;
}

// M - tm => newM
Matrix3 Matrix3::operator-(const Matrix3& tm) const 
{
  Matrix3 res;

  res.xx = xx - tm.xx; 
  res.xy = xy - tm.xy;
  res.xz = xz - tm.xz;
  
  res.yx = yx - tm.yx;
  res.yy = yy - tm.yy;
  res.yz = yz - tm.yz;
  
  res.zx = zx - tm.zx;
  res.zy = zy - tm.zy;
  res.zz = zz - tm.zz;
  return res;

}

// negative
Matrix3 Matrix3::operator-(void) const  
{
  Matrix3 res;
   
  res.xx = -xx;
  res.xy = -xy;
  res.xz = -xz;
  res.yx = -yx;
  res.yy = -yy;
  res.yz = -yz;
  res.zx = -zx;
  res.zy = -zy;
  res.zz = -zz;
  return res;
}

// transpose this matrix to a new one without changing this one
// M.transposed() => newM
Matrix3 Matrix3::transposed() const 
{
  Matrix3 res;
    
  res.xx = xx;
  res.xy = yx;
  res.xz = zx;
  
  res.yx = xy;
  res.yy = yy;
  res.yz = zy;
  
  res.zx = xz;
  res.zy = yz;
  res.zz = zz;
  return res;
}

// set this one as transposed tm, - M.transpose(tm) => M
void Matrix3::transpose(const Matrix3& tm) 
{

  // NB: make sure that *this == tm also works
  Double a = tm.yx;
     yx = tm.xy;
     xy = a;

  Double b  = tm.zx;
      zx = tm.xz;  
      xz = b;

  Double c  = tm.zy;
      zy = tm.yz;
      yz = c;

  xx = tm.xx;
  yy = tm.yy;
  zz = tm.zz;
}

void Matrix3::transpose_me() 
{
  Double a = yx;
     yx = xy;
     xy = a;

  Double b  = zx;
      zx = xz;  
      xz = b;

  Double c  = zy;
      zy = yz;
      yz = c;
}

// M += tm
void Matrix3::operator+=(const Matrix3& tm)
{
  xx += tm.xx; 
  xy += tm.xy;
  xz += tm.xz;
  
  yx += tm.yx;
  yy += tm.yy;
  yz += tm.yz;
  
  zx += tm.zx;
  zy += tm.zy;
  zz += tm.zz;
}

// M + tm => newM
Matrix3 Matrix3::operator+(const Matrix3& tm) const 
{
  Matrix3 res;

  res.xx = xx + tm.xx; 
  res.xy = xy + tm.xy;
  res.xz = xz + tm.xz;
  
  res.yx = yx + tm.yx;
  res.yy = yy + tm.yy;
  res.yz = yz + tm.yz;
  
  res.zx = zx + tm.zx;
  res.zy = zy + tm.zy;
  res.zz = zz + tm.zz;
  return res;
}


Double Matrix3::det()
{
  Double a00 = yy*zz - yz*zy;
  Double a10 = yx*zz - yz*zx;
  Double a20 = yx*zy - yy*zx;
  return xx*a00 - xy*a10 + xz*a20;
}

bool Matrix3::invert()
{
  Matrix3 tmp;
  Double det;

  //
  // Calculate the determinant of submatrix A (optimized version:
  // don,t just compute the determinant of A)
  //
  tmp.xx = yy*zz - yz*zy;
  tmp.yx = yx*zz - yz*zx;
  tmp.zx = yx*zy - yy*zx;

  tmp.xy = xy*zz - xz*zy;
  tmp.yy = xx*zz - xz*zx;
  tmp.zy = xx*zy - xy*zx;

  tmp.xz = xy*yz - xz*yy;
  tmp.yz = xx*yz - xz*yx;
  tmp.zz = xx*yy - xy*yx;

  det = xx*tmp.xx - xy*tmp.yx + xz*tmp.zx;

  //
  // singular matrix ?
  //
//  if (fabs(det) < EPSILON*EPSILON)
//    return false;

  det = 1/det;

  //
  // inverse(A) = adj(A)/det(A)
  //
  tmp.xx *= det;
  tmp.xz *= det;
  tmp.yy *= det;
  tmp.zx *= det;
  tmp.zz *= det;

  det = -det;

  tmp.xy *= det;
  tmp.yx *= det;
  tmp.yz *= det;
  tmp.zy *= det;

  *this = tmp;
  return true;
}

/*
 * Concat scale matrix to the current transformation.
 */

void Matrix3::scale(Double s)
{
  xx *= s; xy *= s; xz *= s;
  yx *= s; yy *= s; yz *= s;
  zx *= s; zy *= s; zz *= s;
}


void Matrix3::scale(Double sx, Double sy, Double sz)
{
  xx *= sx; xy *= sy; xz *= sz;
  yx *= sx; yy *= sy; yz *= sz;
  zx *= sx; zy *= sy; zz *= sz;
}

void Matrix3::scale(const Vector3& scaleFactor)
{
  xx *= scaleFactor.x; xy *= scaleFactor.y; xz *= scaleFactor.z;
  yx *= scaleFactor.x; yy *= scaleFactor.y; yz *= scaleFactor.z;
  zx *= scaleFactor.x; zy *= scaleFactor.y; zz *= scaleFactor.z;
}

// convert Matrix elements into an array
void convert(const Matrix3& from, Double to[9])
{
  to[0]  = from.xx; to[1]  = from.xy; to[2]  = from.xz;
  to[3]  = from.yx; to[4]  = from.yy; to[5]  = from.yz;
  to[6]  = from.zx; to[7]  = from.zy; to[8]  = from.zz;
}
