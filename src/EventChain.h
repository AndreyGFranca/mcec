#pragma once
#include <pybind11/pybind11.h>

using namespace mcec
{

	class EventChain
	{
	public:
		EventChain();
		~EventChain();

		void evolveSystem(mcec::System system, unsigned int steps);

	private:

	};

}