// Copyright (c) 2024 Thomas Klietsch, all rights reserved.
//
// Licensed under the GNU Lesser General Public License, version 3.0 or later
//
// This program is free software: you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation, either version 3 of
// the License, or ( at your option ) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General
// Public License along with this program.If not, see < https://www.gnu.org/licenses/>. 

#pragma once

#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>
#include <limits>
#include <numbers>
#include <sstream>
#include <vector>

// C++23
#if __STDCPP_FLOAT128_T__ == 1
#include <stdfloat>
using Real = std::float128_t;
#else
using Real = long double;
#endif

std::string RealToString(
	Real const& value,
	uint8_t const& decimals = 8)
{
	auto constexpr max_digits{ std::numeric_limits<Real>::digits10 + 1 };

	// Decimals = '0' outputs all available decimals,
	// unlike std::setprecision(0) which outputs none.
	// A prefix space is added if value is positive (or zero)
	std::stringstream stream{
		decimals ?
		std::format("{: .{}}", value, std::min<uint8_t>(decimals, max_digits)) :
		std::format("{: }", value)
	};
	return stream.str();
};

#if __STDCPP_FLOAT128_T__ == 1
// Currently no support in C++23 for std::cout of std::float128_t
std::ostream& operator<<(
	std::ostream& os,
	Real const& value)
{
	return os << RealToString(value, std::cout.precision());
};
#endif

Real constexpr pi = std::numbers::pi_v<Real>;

// Return 'Not a Number', without throwing an exception
Real constexpr NaN = std::numeric_limits<long double>::quiet_NaN();

// Numeric stability
// Smallest value such that 1+epsilon evaluates to 1
Real constexpr numeric_epsilon = std::numeric_limits<Real>::epsilon();

// The magnitude below which a number is considered to be zero
Real constexpr magnitude_zero = std::numeric_limits<float>::epsilon();

namespace GaussJordan
{

	// Only square or overdetermined matrices are supported
	class Matrix final
	{

	private:

		uint8_t n_rows{ 0 };
		uint8_t n_columns{ 0 };

		std::vector<std::vector<Real>> cell{ 0 };

		// Used if out of range row/column is used
		// for read/write to matrix cell
		Real dummy{ NaN };

	public:

		Matrix() {};

		Matrix(
			uint8_t const& rows,
			uint8_t const& columns)
		{
			if (rows < columns)
				return;

			n_rows = rows;
			n_columns = columns;

			cell = std::vector<std::vector<Real>>(rows, std::vector<Real>(columns, 0));
		};

		// Read/Write the value at [row,column]
		Real& operator () (
			uint8_t const& row,
			uint8_t const& column)
		{
			if ((row >= n_rows) || (column >= n_columns))
				return dummy;

			return cell[row][column];
		};

		// Get the value at [row,column]
		Real operator () (
			uint8_t const& row,
			uint8_t const& column) const
		{
			if ((row >= n_rows) || (column >= n_columns))
				return NaN;

			return cell[row][column];
		};

		void set_row(
			uint8_t const& index,
			std::vector<Real> const& data)
		{
			if ((!n_columns) || (data.size() != n_columns))
				return;

			for (uint16_t i{ 0 };i < n_columns;++i)
				if (!std::isfinite(data[i]))
					return;

			cell[index] = data;
		};

		void set_column(
			uint8_t const& index,
			std::vector<Real> const& data)
		{
			if ((!n_rows) || (data.size() != n_rows))
				return;

			for (uint16_t i{ 0 };i < n_rows;++i)
				if (!std::isfinite(data[i]))
					return;

			for (uint16_t i{ 0 };i < n_rows;++i)
				cell[i][index] = data[i];
		};

		// Returns number of rows and columns
		std::tuple< uint8_t, uint8_t> size() const
		{
			return std::tuple(n_rows, n_columns);
		};

		// Tests diagonal for zeroes
		bool is_diagonal_nonzero(
			Real const& epsilon = magnitude_zero) const
		{
			if (!is_valid())
				return false;

			for (uint16_t i{ 0 }; i < n_columns; ++i)
				if (std::abs(cell[i][i]) < epsilon)
					return false;

			return true;
		};

		// True if the matrix is square or overdetermined,
		// and neither number of rows or columns are zero
		bool is_valid() const
		{
			if (!n_rows || !n_columns || (n_rows < n_columns))
				return false;

			return true;
		};

		void print() const
		{
			std::cout << n_rows + 0 << "x" << n_columns + 0 << " matrix\n";
			for (uint16_t i = 0;i < n_rows;++i)
			{
				for (uint16_t j = 0;j < n_columns;++j)
					std::cout << cell[i][j] << " ";
				std::cout << "\n";
			}
		};

	};

	// Gauss Jordan elimination
	std::vector<Real> Solve(
		Matrix matrix,
		std::vector<Real> equal,
		Real const& epsilon = magnitude_zero)
	{
		if (!matrix.is_valid())
			return {};

		auto const [n_rows, n_columns] = matrix.size();

		if (n_rows != equal.size())
			return {};

		// Method uses the diagonal, so it must be non-zero
		if (!matrix.is_diagonal_nonzero(epsilon))
		{
			// TODO Mathness: craptastic code I threw together

			// The non-zero entries for each column
			std::vector<std::vector<int16_t>> nze_list(n_columns);
			for (uint16_t column{ 0 };column < n_columns;++column)
			{
				for (uint16_t row{ 0 };row < n_rows;++row)
					if (std::abs(matrix(row, column)) >= epsilon)
						nze_list[column].push_back(row);
				// Column has only zeroes
				if (!nze_list[column].size())
					return {};
			}

			// Set current used index in row_list
			std::vector<uint16_t> list_index(n_columns, 0);
			// Sort matrix so that the major diagonal is non-zero
			while (1)
			{
				std::vector<bool> row_used(n_rows, false);

				// Set row_used state
				for (uint16_t i{ 0 }; i < n_rows; ++i)
					row_used[list_index[i]] = true;

				uint16_t n_used{ 0 };
				for (uint16_t i{ 0 }; i < n_rows; ++i)
					if (row_used[i])
						++n_used;

				if (n_used == n_columns)
					break;
				else
				{
					// No valid major diagonal found, increment test_indices
					++list_index[0];
					for (uint16_t i{ 0 }; i < n_columns - 1; ++i)
						if (list_index[i] >= nze_list[i].size())
						{
							list_index[i] = 0;
							++list_index[i + 1];
						}

					// No combination gives a non-zero diagonal
					if (list_index[n_columns - 1] >= nze_list[n_columns - 1].size())
						return {};
				}
			};

			// Construct new (square) matrix and matching equal
			Matrix new_matrix{ n_columns, n_columns };
			std::vector<Real> new_equal(n_columns);
			for (uint16_t row{ 0 }; row < n_columns; ++row)
			{
				uint8_t const index = list_index[row];
				new_equal[row] = equal[index];
				for (uint16_t column{ 0 }; column < n_columns; ++column)
					new_matrix(row, column) = matrix(index, column);
			}
			std::swap(matrix, new_matrix);
			std::swap(equal, new_equal);
		}

		// Solve as a square matrix
		// Gauss Jordan elimination
		for (uint16_t column{ 0 }; column < n_columns; ++column)
			for (uint16_t row{ 0 }; row < n_columns; ++row)
				if (row != column)
				{
					Real const scalar = matrix(row, column) / matrix(column, column);
					if (!std::isfinite(scalar))
						return {};

					for (uint16_t k{ 0 }; k < n_columns; ++k)
						matrix(row, k) -= matrix(column, k) * scalar;

					equal[row] -= equal[column] * scalar;
				}

		// Rescale diagonal to 1, to get result
		std::vector<Real> result(n_columns);
		for (uint16_t i{ 0 }; i < n_columns; ++i)
		{
			Real value = equal[i] / matrix(i, i);
			if (!std::isfinite(value))
				return {};

			result[i] = value;
		}

		return result;
	};

	// For a square matrix this should be zero.
	// An overdetermined matrix has more equations than unknowns,
	// so the result may not be correct for all of them
	// TODO Mathness: does this have a proper name?
	Real ErrorEstimate(
		Matrix const& matrix,
		std::vector<Real> const& equal,
		std::vector<Real> const& result)
	{
		if (!matrix.is_valid())
			return NaN;

		auto const [n_rows, n_columns] = matrix.size();

		if ((n_rows != equal.size()) || (n_columns != result.size()))
			return NaN;

		Real error{ 0 };
		for (uint16_t i{ 0 }; i < n_rows; ++i)
		{
			Real sum{ 0 };
			for (uint16_t j{ 0 };j < n_columns;++j)
				sum += matrix(i, j) * result[j];

			error += std::abs(sum - equal[i]);
		}

		return error;
	};

};
