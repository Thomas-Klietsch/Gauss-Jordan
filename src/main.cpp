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

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "./gauss_jordan.hpp"

int main()
{
	// Set up equations to solve for Kronrod extension of Gauss Lobatto quadrature
	std::vector<Real> equal{ 2.q, 2.q / 3.q, 2.q / 5.q, 2.q / 7.q, 2.q / 9.q };

	Real x{ 1.q / 5.q };
	GaussJordan::Matrix matrix(5, 4);
	matrix.set_column(0, { 2, 2, 2, 2, 2 });
	matrix.set_column(1, { 2, 2 * x, 2 * x * x, 2 * x * x * x, 2 * x * x * x * x });
	matrix.set_column(3, { 1, 0, 0, 0, 0 });

	// Generate data for gnuplot
	std::ofstream file("plot_data.txt", std::ios::out);
	std::cout << std::setprecision(10);
	if (file.is_open())
	{
		for (uint16_t i{ 0 }; i < 500; ++i)
		{
			x = 0.002q * i;
			matrix.set_column(2, { 2, 2 * x, 2 * x * x, 2 * x * x * x, 2 * x * x * x * x });

			auto result = GaussJordan::Solve(matrix, equal);

			Real y = GaussJordan::ErrorEstimate(matrix, equal, result);

			if (std::isfinite(y))
				file << x << "  " << y << "\n";
		}
		file.close();

		std::system("gnuplot plot_config.txt");
	}

	// Bracket search for minimal "error"
	Real const start{ 1.q / 5.q };
	Real const end{ 1.q };
	x = (start + end) / 2;
	Real const span = std::abs(end - x);
	uint8_t loop{ 1 };
	do
	{
		int8_t index{ 0 };
		Real min_error{ 1e20 };
		Real const h = span * std::pow(4, -loop);
		for (int8_t i{ -3 };i <= 3;++i)
		{
			Real t = x + h * i;
			matrix.set_column(2, { 2, 2 * t, 2 * t * t, 2 * t * t * t, 2 * t * t * t * t });
			auto result = GaussJordan::Solve(matrix, equal);
			Real error = GaussJordan::ErrorEstimate(matrix, equal, result);
			if (error < min_error)
			{
				index = i;
				min_error = error;
			}
		}
		x += index * h;
		if ((min_error < numeric_epsilon) || (h < numeric_epsilon * 2))
			break;
	} while (++loop < 100);

	auto result = GaussJordan::Solve(matrix, equal);
	std::cout << std::setprecision(20) << "x = " << x << "\n";
	std::cout << "A = " << result[0] << "\n";
	std::cout << "B = " << result[1] << "\n";
	std::cout << "C = " << result[2] << "\n";
	std::cout << "D = " << result[3] << "\n";

	return 0;
};
