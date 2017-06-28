#include <cstring>

namespace fftlib_internal
{
	using namespace fftlib;

	template<backend B, transformation T>
	void composed_fft<B, T>::cleanup()
	{
		if (planner_buffer_1 != NULL)
			_mm_free(planner_buffer_1);

		planner_buffer_1 = NULL;
		planner_buffer_1_size = 0;

		if (planner_buffer_2 != NULL)
			_mm_free(planner_buffer_2);

		planner_buffer_2 = NULL;
		planner_buffer_2_size = 0;
	}

	template<backend B, transformation T>
	template<std::int32_t D>
	void composed_fft<B, T>::create(typename fftlib_internal::ext_plan& p, const typename data<B, INT>::type_t* n, typename trafo<B, T>::in_t* in, typename trafo<B, T>::out_t* out, const typename data<B, INT>::type_t sign, const typename data<B, UINT>::type_t flags)
	{
		// empty
	}

	template<>
	template<>
	void composed_fft<FFTW, C2C_64>::create<3>(typename fftlib_internal::ext_plan& p, const typename data<FFTW, INT>::type_t* n, typename trafo<FFTW, C2C_64>::in_t* in, typename trafo<FFTW, C2C_64>::out_t* out, const typename data<FFTW, INT>::type_t sign, const typename data<FFTW, UINT>::type_t flags)
	{
		#if defined(FFTLIB_BALL_CUBE_OPT)
			// default values.
			p.ball_cube = NONE;
			p.howmany = n[0];
			// use ball-cube optimization only if at least two dimensions are larger than 48.
			if (((n[0] > 32 ? 1 : 0) + (n[1] > 32 ? 1 : 0) + (n[2] > 32 ? 1 : 0)) >= 2)
				{
					// points in xy-plane used for scanning in z-direction.
					const std::int32_t base_pos[] = {0, n[2] / 2, n[2] - 1,
									(n[1] / 2) * n[2] + 0, (n[1] / 2) * n[2] + n[2] / 2, (n[1] / 2) * n[2] + n[2] - 1,
									(n[1] - 1) * n[2] + 0, (n[1] - 1) * n[2] + n[2] / 2, (n[1] - 1) * n[2] + n[2] - 1};
					// initial value for the number of layers containing zeros only: we will take the minimum afterwards (which can be 0).
					std::int32_t zero_layers = n[0];
					// scan from below (0 to n[0] / 2).
					for (std::int32_t i = 0; i < (sizeof(base_pos) / sizeof(std::int32_t)); ++i)
						{
							std::int32_t z;
							const typename data<FFTW, R_64>::type_t* ptr_in = reinterpret_cast<const typename data<FFTW, R_64>::type_t*>(&in[base_pos[i]]);
							// scan slab in z-direction.
							for (z = 0; z < (n[0] / 2); ++z)
								{
									// if values are non-zero, break the loop.
									if (ptr_in[0] != 0.0 || ptr_in[1] != 0.0)
										break;
									// otherwise, continue with the next layer.
									ptr_in += (2 * n[1] * n[2]);
								}
							// determine the new value of 'zero_layers'.
							zero_layers = std::min(zero_layers, z);
							// if it is zero, break the outer loop, as it cannot become larger than zero anymore.
							if (zero_layers == 0)
								break;
						}
					// there is a ball in the center.
					if (zero_layers > 0 && zero_layers != (n[0] / 2))
						{
							// use the ball-cube optimization only if the fraction of the zero-layers on n[0] satisfies the threshold criteria.
							if (((2.0 * zero_layers) / n[0]) > FFTLIB_BALL_CUBE_THRESHOLD)
								{
									p.ball_cube = CENTER;
									p.howmany = n[0] - 2 * (zero_layers - 1);
									p.howmany = ((p.howmany + 1) / 2) * 2;
									p.howmany = std::min(p.howmany, n[0] - zero_layers);
									p.offset_low = zero_layers;
								}
						}
					// if the ball is not in the center.
					else if (zero_layers == 0)
						{
							zero_layers = n[0];
							// scan from center (n[0] / 2 to 0).
							for (std::int32_t i = 0; i < (sizeof(base_pos) / sizeof(std::int32_t)); ++i)
								{
									std::int32_t z;
									const typename data<FFTW, R_64>::type_t* ptr_in = reinterpret_cast<const typename data<FFTW, R_64>::type_t*>(&in[(n[0] / 2) * n[1] * n[2] + base_pos[i]]);
									// scan slab in z-direction.
									for (z = (n[0] / 2); z > 0; --z)
										{
											// if values are non-zero, break the loop.
											if (ptr_in[0] != 0.0 || ptr_in[1] != 0.0)
												break;
											// otherwise, continue with the next layer.
											ptr_in -= (2 * n[1] * n[2]);
										}
									// determine the new value of 'zero_layers'.
									zero_layers = std::min(zero_layers, (n[0] / 2) - z);
									// if it is zero, break the outer loop, as it cannot become larger than zero anymore.
									if (zero_layers == 0)
										break;
								}
							// there is a ball at the bottom.
							if (zero_layers > 0 && zero_layers != (n[0] / 2))
								{
									// use the ball-cube optimization only if the fraction of the zero-layers on n[0] satisfies the threshold criteria.
									if (((2.0 * zero_layers) / n[0]) > FFTLIB_BALL_CUBE_THRESHOLD)
										{
											p.ball_cube = BORDER;
											p.howmany = (n[0] / 2) - (zero_layers - 1);
											p.howmany = ((p.howmany + 1) / 2) * 2;
											p.howmany = std::min(p.howmany, n[0] / 2);
											p.offset_low = 0;
											p.offset_high = n[0] - p.howmany;
										}
								}
						}
				}
		#else
			// if ball-cube optimization is not enabled, all n[0] layers in z-direction are transformed.
			p.ball_cube = NONE;
			p.howmany = n[0];

		#endif

		#if defined(FFTLIB_USE_MKL)
			typename trafo<FFTW, C2C_64>::in_t* ptr_in = in;
			typename trafo<FFTW, C2C_64>::out_t* ptr_out = out;
		#else
			typename trafo<FFTW, C2C_64>::in_t* ptr_in = NULL;
			typename trafo<FFTW, C2C_64>::out_t* ptr_out = NULL;

			std::size_t size = n[0] * n[1] * n[2] * std::max(sizeof(typename trafo<FFTW, C2C_64>::in_t), sizeof(typename trafo<FFTW, C2C_64>::out_t));
			if (planner_buffer_1_size < size)
				{
					_mm_free(planner_buffer_1);
					planner_buffer_1_size = static_cast<std::size_t>(1.2 * size);
					planner_buffer_1 = reinterpret_cast<void*>(_mm_malloc(planner_buffer_1_size, FFTLIB_ALIGNMENT));
				}

			ptr_in = reinterpret_cast<typename trafo<FFTW, C2C_64>::in_t*>(planner_buffer_1);
			
			if (reinterpret_cast<void*>(in) == reinterpret_cast<void*>(out))
				{
					ptr_out = reinterpret_cast<typename trafo<FFTW, C2C_64>::out_t*>(planner_buffer_1);
				}
			else
				{
					std::size_t size = n[0] * n[1] * n[2] * sizeof(typename trafo<FFTW, C2C_64>::out_t);
					if (planner_buffer_2_size < size)
						{
							_mm_free(planner_buffer_2);
							planner_buffer_2_size = static_cast<std::size_t>(1.2 * size);
							planner_buffer_2 = reinterpret_cast<void*>(_mm_malloc(planner_buffer_2_size, FFTLIB_ALIGNMENT));
						}
					ptr_out = reinterpret_cast<typename trafo<FFTW, C2C_64>::out_t*>(planner_buffer_2);
				}
		#endif

		// 2dfft
		fftlib_internal::configuration<FFTW, C2C_64, 2> c_2dfft = fftlib_internal::configuration<FFTW, C2C_64, 2>(&n[1], p.howmany, ptr_in, &n[1], 1, n[1] * n[2], ptr_out, &n[1], 1, n[1] * n[2], sign, FFTW_MEASURE, p.nthreads);
		fft<FFTW, C2C_64, 2>& my_2dfft = fft<FFTW, C2C_64, 2>::get_instance();
		if ((p.p_1 = reinterpret_cast<const void*>(my_2dfft.find_plan(c_2dfft))) == NULL)
			p.p_1 = reinterpret_cast<const void*>(my_2dfft.create_plan(c_2dfft, ptr_in, ptr_out));
		
		// 1dfft
		fftlib_internal::configuration<FFTW, C2C_64, 1> c_1dfft = fftlib_internal::configuration<FFTW, C2C_64, 1>(&n[0], n[1] * n[2], ptr_out, &n[0], n[1] * n[2], 1, ptr_out, &n[0], n[1] * n[2], 1, sign, FFTW_MEASURE, p.nthreads);
		fft<FFTW, C2C_64, 1>& my_1dfft = fft<FFTW, C2C_64, 1>::get_instance();
		if ((p.p_2 = reinterpret_cast<const void*>(my_1dfft.find_plan(c_1dfft))) == NULL)
			p.p_2 = reinterpret_cast<const void*>(my_1dfft.create_plan(c_1dfft, ptr_out, ptr_out));
	}

	template<backend B, transformation T>
	void composed_fft<B, T>::execute(typename fftlib_internal::ext_plan& p, typename trafo<B, T>::in_t* in, typename trafo<B, T>::out_t* out)
	{
		// empty
	}

	template<>
	void composed_fft<FFTW, C2C_64>::execute(typename fftlib_internal::ext_plan& p, typename trafo<FFTW, C2C_64>::in_t* in, typename trafo<FFTW, C2C_64>::out_t* out)
	{
		fftlib_internal::dynamic_lib<FFTW>& dl = fftlib_internal::dynamic_lib<FFTW>::get_instance();

		// 2dfft
		#if defined(FFTLIB_BALL_CUBE_OPT)
			switch (p.ball_cube)
				{
				case NONE:
					dl.template execute_plan<C2C_64>(*reinterpret_cast<const typename trafo<FFTW, C2C_64>::plan_t*>(p.p_1), in, out, p.nthreads);
					break;
					
				case CENTER:
					dl.template execute_plan<C2C_64>(*reinterpret_cast<const typename trafo<FFTW, C2C_64>::plan_t*>(p.p_1), &in[p.offset_low * p.n[1] * p.n[2]], &out[p.offset_low * p.n[1] * p.n[2]], p.nthreads);
					break;
					
				case BORDER:
					dl.template execute_plan<C2C_64>(*reinterpret_cast<const typename trafo<FFTW, C2C_64>::plan_t*>(p.p_1), &in[p.offset_low * p.n[1] * p.n[2]], &out[p.offset_low * p.n[1] * p.n[2]], p.nthreads);
					dl.template execute_plan<C2C_64>(*reinterpret_cast<const typename trafo<FFTW, C2C_64>::plan_t*>(p.p_1), &in[p.offset_high * p.n[1] * p.n[2]], &out[p.offset_high * p.n[1] * p.n[2]], p.nthreads);
					break;
				}
		#else
			dl.template execute_plan<C2C_64>(*reinterpret_cast<const typename trafo<FFTW, C2C_64>::plan_t*>(p.p_1), in, out, p.nthreads);
		#endif

		// 1dfft
		dl.template execute_plan<C2C_64>(*reinterpret_cast<const typename trafo<FFTW, C2C_64>::plan_t*>(p.p_2), out, out, p.nthreads);
	}
}
