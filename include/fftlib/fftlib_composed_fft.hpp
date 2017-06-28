// Copyright (c) 2016 Florian Wende (flwende@gmail.com)
//
// Distributed under the BSD 2-clause Software License 
// (See accompanying file LICENSE)

#if !defined(FFTLIB_COMPOSED_FFT_HPP)
#define FFTLIB_COMPOSED_FFT_HPP

namespace fftlib_internal
{
	template<backend B = FFTW, transformation T = C2C_64>
	class composed_fft
	{
	public:
		static void cleanup();

		template<std::int32_t D = 3>
		static void create(typename fftlib_internal::ext_plan& p, const typename data<B, INT>::type_t* n, typename trafo<B, T>::in_t* in, typename trafo<B, T>::out_t* out, const typename data<B, INT>::type_t sign, const typename data<B, UINT>::type_t flags);

		static void execute(typename fftlib_internal::ext_plan& p, typename trafo<B, T>::in_t* in, typename trafo<B, T>::out_t* out);
		
	private:
		static std::int32_t extend[32];

		static void* planner_buffer_1;

		static void* planner_buffer_2;

		static std::size_t planner_buffer_1_size;

		static std::size_t planner_buffer_2_size;
	};

	template<backend B, transformation T>
	void* composed_fft<B, T>::planner_buffer_1 = NULL;

	template<backend B, transformation T>
	void* composed_fft<B, T>::planner_buffer_2 = NULL;

	template<backend B, transformation T>
	std::size_t composed_fft<B, T>::planner_buffer_1_size = 0;

	template<backend B, transformation T>
	std::size_t composed_fft<B, T>::planner_buffer_2_size = 0;
}

#endif
