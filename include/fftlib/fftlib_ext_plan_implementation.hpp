// Copyright (c) 2016 Florian Wende (flwende@gmail.com)
//
// Distributed under the BSD 2-clause Software License 
// (See accompanying file LICENSE)

namespace fftlib_internal
{
	ext_plan::ext_plan()
		: nthreads(1),
		  composed_fft(false),
		  ball_cube(NONE),
		  offset_low(0),
		  offset_high(0),
		  in(NULL),
		  out(NULL),
		  p_1(NULL),
		  p_2(NULL),
		  p_3(NULL)
	{
	}

	ext_plan::~ext_plan()
	{
	}
}
