// //////////////////////////////////////////////////////////
// fnv.h
// Copyright (c) 2011 Stephan Brumme. All rights reserved.
// see http://create.stephan-brumme.com/disclaimer.html
//

// compute FNV1a hash
// originally developed by Fowler, Noll and Vo
// http://isthe.com/chongo/tech/comp/fnv/

#pragma once

#include <string>
#include <cassert>


// 32 bit integers
#ifdef _MSC_VER
typedef unsigned int uint32_t;
#else
#include <stdint.h>
#endif


namespace FNV
{
  // default values recommended by http://isthe.com/chongo/tech/comp/fnv/
  const uint32_t Prime = 0x01000193; //   16777619
  const uint32_t Seed  = 0x811C9DC5; // 2166136261

  /// hash a single byte
  inline uint32_t fnv1a(unsigned char oneByte, uint32_t hash = Seed)
  {
    return (oneByte ^ hash) * Prime;
    // FNV1: return (oneByte * Prime) ^ hash;
  }


  /// hash a short (two bytes)
  inline uint32_t fnv1a(unsigned short twoBytes, uint32_t hash = Seed)
  {
    const unsigned char* ptr = (const unsigned char*) &twoBytes;
    hash = fnv1a(*ptr++, hash);
    return fnv1a(*ptr  , hash);
  }


  /// hash a 32 bit integer (four bytes)
  inline uint32_t fnv1a(uint32_t fourBytes, uint32_t hash = Seed)
  {
    const unsigned char* ptr = (const unsigned char*) &fourBytes;
    hash = fnv1a(*ptr++, hash);
    hash = fnv1a(*ptr++, hash);
    hash = fnv1a(*ptr++, hash);
    return fnv1a(*ptr  , hash);
  }


  /// hash a block of memory
  inline uint32_t fnv1a(const void* data, size_t numBytes, uint32_t hash = Seed)
  {
    assert(data);
    const unsigned char* ptr = (const unsigned char*)data;
    while (numBytes--)
      hash = (*ptr++ ^ hash) * Prime;
      // same as hash = fnv1a(*ptr++, hash); but much faster in debug mode
    return hash;
  }


  /// hash a C-style string
  inline uint32_t fnv1a(const char* text, uint32_t hash = Seed)
  {
    assert(text);
    const unsigned char* ptr = (const unsigned char*)text;
    while (*ptr)
      hash = (*ptr++ ^ hash) * Prime;
      // same as hash = fnv1a(*ptr++, hash); but much faster in debug mode
    return hash;
  }


  /// hash an std::string
  inline uint32_t fnv1a(const std::string& text, uint32_t hash = Seed)
  {
    return fnv1a(text.c_str(), text.length(), hash);
    // or: fnv1a(text.c_str(), hash);
  }


  /// hash a float
  inline uint32_t fnv1a(float number, uint32_t hash = Seed)
  {
    return fnv1a(&number, sizeof(number), hash);
  }


  /// hash a double
  inline uint32_t fnv1a(double number, uint32_t hash = Seed)
  {
    return fnv1a(&number, sizeof(number), hash);
  }


  /// hash a block of memory
  template <unsigned int Unroll>
  uint32_t fnv1a_unrolled(const void* data, size_t numBytes, uint32_t hash = Seed)
  {
#ifndef NDEBUG
    // unrolling isn't performed when optimizations are disabled
    return fnv1a(data, numBytes, hash);
#endif
    assert(data);
    const unsigned char* ptr = (const unsigned char*)data;
    // unroll
    while (numBytes >= Unroll)
    {
      // Unroll is a constant and smart compilers like GCC and Visual C++ unroll properly
      hash = fnv1a(ptr, Unroll, hash);
      ptr += Unroll;
      numBytes -= Unroll;
    }
    // process remaining bytes
    return fnv1a(ptr, numBytes, hash);
  }
  /// catch invalid Unroll value
  template <>
  inline uint32_t fnv1a_unrolled<0>(const void* data, size_t numBytes, uint32_t hash /*= Seed*/)
  {
    return fnv1a(data, numBytes, hash);
  }
  /// not unrolled at all
  template <>
  inline uint32_t fnv1a_unrolled<1>(const void* data, size_t numBytes, uint32_t hash /*= Seed*/)
  {
    return fnv1a(data, numBytes, hash);
  }
}
