#pragma once
#include <boost/numeric/ublas/matrix.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;


namespace mcec
{

	class CellList
	{
	public:
		CellList(unsigned long int nums_particle, unsigned int lx, unsigned int ly);
		~CellList();


	private:
		matrix<float> m_cell_list;


	};

}
