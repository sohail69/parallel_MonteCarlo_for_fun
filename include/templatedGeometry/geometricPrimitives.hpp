#pragma once
#include <array>
#include <vector>

//Facades used for making
//the code more clear
template<typename real, size_t sdim>
using Point = std::array<real,sdim>

template<typename real, size_t sdim>
using VecND = std::array<real,sdim>
