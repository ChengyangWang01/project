#ifndef BASE_HPP
#define BASE_HPP

#include <vector>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>

template <typename T>
struct Matrix3D
{
  public:
  typedef T           value_type;
  typedef T&          value_reference;
  typedef const T&    const_reference;
  typedef T*          iterator_type;
  typedef int         size_type;

  public:
  size_type nx, ny, nz;

  value_reference operator()(size_type);
  const_reference operator()(size_type) const;

  value_reference operator()(size_type, size_type, size_type);
  const_reference operator()(size_type, size_type, size_type) const;

  iterator_type address(size_type, size_type, size_type);

  public:
  Matrix3D();
  Matrix3D(size_type, size_type, size_type);

  template <typename U>
  friend std::ostream& operator<<(std::ostream&, const Matrix3D<U>&);

  void write_to_file(const std::string&) const;
  void swap(Matrix3D &);

  private:
  std::vector<value_type> __data;
  size_type __indexing(size_type, size_type, size_type) const;



}; // end of struct Matrix 


// =================================== CPP parts ====================================

template <typename T>
void Matrix3D<T>::write_to_file(const std::string& fn) const 
{
  std::ofstream ofs(fn);

  ofs << static_cast<std::size_t>(nx) << ' ' 
      << static_cast<std::size_t>(ny) << ' ' 
      << static_cast<std::size_t>(nz) << '\n';

  ofs << *this;

  ofs.close();
}

template <typename T>
void Matrix3D<T>::swap(Matrix3D<T> & other) 
  { __data.swap(other.__data); }

template <typename T>
std::ostream& operator<<(std::ostream& os, const Matrix3D<T>& in)
{
  for (std::size_t ix = 0; ix < in.nx; ++ix)
  {
    os << "";
    for (std::size_t iy = 0; iy < in.ny; ++iy)
    {
      os << "";
      for (std::size_t iz = 0; iz < in.nz; ++iz)
      {
        os << std::fixed << std::setw(9) << std::setprecision(3) << in(ix, iy, iz) << " ";
      }
      os << "" << std::endl;
    }
    os << "" << std::endl;
  }
  return os;
}

/// @brief 
/// @tparam T 
/// @param x_idx 
/// @param y_idx 
/// @param z_idx 
/// @return 
template <typename T>
  inline
  int Matrix3D<T>::__indexing(size_type x_idx, size_type y_idx, size_type z_idx) 
    const
    {
      assert(x_idx < nx);
      assert(y_idx < ny);
      assert(z_idx < nz);

      return x_idx * (nz * ny) + y_idx * nz + z_idx;
    }

template <typename T>
  inline
  T* Matrix3D<T>::address(size_type x_idx, size_type y_idx, size_type z_idx)
  { return __data.data() + __indexing(x_idx, y_idx, z_idx); }

/// @brief 
/// @tparam T 
/// @param nx 
/// @param ny 
/// @param nz 
template <typename T>
  inline
  Matrix3D<T>::Matrix3D(size_type nx, size_type ny, size_type nz)
    : nx{nx}, ny{ny}, nz{nz}, __data{std::vector<value_type>(nx*ny*nz, 0)}
    { }

template <typename T>
  inline 
  Matrix3D<T>::Matrix3D()
   : nx{0}, ny{0}, nz{0}, __data{std::vector<value_type>(0, 0)}
   { }


/// @brief 
/// @tparam T 
/// @param index 
/// @return 
template <typename T>
  inline
  T& Matrix3D<T>::operator()(size_type index)
    { return __data.at(index); }

/// @brief 
/// @tparam T 
/// @param x_idx 
/// @param y_idx 
/// @param z_idx 
/// @return 
template <typename T>
  inline
  T& Matrix3D<T>::operator()(size_type x_idx, size_type y_idx, size_type z_idx) 
  { return __data[ __indexing(x_idx, y_idx, z_idx) ]; }

/// @brief 
/// @tparam T 
/// @param index 
/// @return 
template <typename T>
  inline
  const 
  T& Matrix3D<T>::operator()(size_type index) 
  const
    { return __data.at(index); }


/// @brief 
/// @tparam T 
/// @param x_idx 
/// @param y_idx 
/// @param z_idx 
/// @return 
template <typename T>
  inline
  const T& Matrix3D<T>::operator()(size_type x_idx, size_type y_idx, size_type z_idx) 
  const
  { return __data[ __indexing(x_idx, y_idx, z_idx) ]; }


#endif // end def BASE_HPP