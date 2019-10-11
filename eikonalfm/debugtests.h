#pragma once

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>

#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
// Replace _NORMAL_BLOCK with _CLIENT_BLOCK if you want the allocations to be of _CLIENT_BLOCK type
#define new DBG_NEW


/*
// add this 3 lines to the files for memory leak detection
#ifdef _DEBUG
	#include "debugtests.h"
#endif
*/