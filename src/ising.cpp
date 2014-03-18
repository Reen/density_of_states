#include "simulation.h"


#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#ifndef __USE_GNU
#define __USE_GNU
#endif

#include <execinfo.h>

#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <cxxabi.h>

#include "config.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

void my_terminate(void);

namespace {
	// invoke set_terminate as part of global constant initialization
	static const bool SET_TERMINATE = std::set_terminate(my_terminate);
}

void my_terminate() {
	static bool tried_throw = false;

	try {
		// try once to re-throw currently active exception
		if (!tried_throw) {
			tried_throw = true;
			throw;
		}
	}
	catch (const std::exception &e) {
		std::cerr << __FUNCTION__ << " caught unhandled exception. what(): "
			<< e.what() << std::endl;
	}
	catch (...) {
		std::cerr << __FUNCTION__ << " caught unknown/unhandled exception."
			<< std::endl;
	}

	void * array[50];
	int size = backtrace(array, 50);

	std::cerr << __FUNCTION__ << " backtrace returned "
		<< size << " frames\n\n";

	char ** messages = backtrace_symbols(array, size);

	for (int i = 0; i < size && messages != NULL; ++i) {
		std::cerr << "[bt]: (" << i << ") " << messages[i] << std::endl;
	}
	std::cerr << std::endl;

	free(messages);

	abort();
}

int main(int argc, char *argv[])
{
#ifdef USE_MPI
	MPI_Init(&argc, &argv);
#endif
	Simulation s;
	int ret = s.exec(argc, argv);
#ifdef USE_MPI
	MPI_Finalize();
#endif
	return ret;
}
