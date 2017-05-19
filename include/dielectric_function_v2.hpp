#ifndef TCM_DIELECTRIC_FUNCTION_HPP
#define TCM_DIELECTRIC_FUNCTION_HPP

#include <cmath>
#include <cassert>

#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <regex>
#include <algorithm>

#include <boost/core/demangle.hpp>

#include <benchmark.hpp>
#include <logging.hpp>

#include <constants.hpp>
#include <matrix.hpp>
#include <blas.hpp>




namespace tcm {

template <class T> struct huge_val;

template <> 
struct huge_val<float> { static constexpr float value = HUGE_VALF; };

template <> 
struct huge_val<double> { static constexpr double value = HUGE_VAL; };

template <> 
struct huge_val<long double> { static constexpr long double value = HUGE_VALL; };

constexpr float huge_val<float>::value;


///////////////////////////////////////////////////////////////////////////////
/// \brief Computes \f$ f(E) \f$

/// Calculates the Fermi-Dirac distribution:
/// \f[ f(E) := \frac{1}{\exp(\frac{E - \mu}{k_{\text{B}}T}) + 1}. \f]
///
/// In the calculation we use that p
/// \f[ \begin{aligned}
///         \lim_{E\to -\infty} f(E) &= 1 //
///         \lim_{E\to +\infty} f(E) &= 0
///     \end{algned}
/// \f]
/// If `std::exp` returns HUGE_VAL, HUGE_VALF or HUGE_VALL, 0 is returned.
/// \param     E       Energy, i.e. a real number. Thus `_F` must be floating
///                    point: `std::is_floating_point<_F>::value == true`.
/// \param     mu      Chemical potential \f$ \mu \f$. It has the dimensions
///                    of energy \f$\implies\f$ `_R` is also floating point.
/// \param     t       Temperature \f$ T \f$.
/// \param     kb      Boltzmann constant \f$ k_\text{B} \f$.
///
/// \return \f$ f(E) \f$.
///////////////////////////////////////////////////////////////////////////////
template<class _F, class _R>
auto fermi_dirac(_F const E, _R const t, _R const mu, _R const kb ) noexcept
{
	static_assert(std::is_floating_point<_F>::value, "Energy must be real.");
	static_assert(std::is_floating_point<_R>::value, "Chemical potential, " 
		"temperature and Boltzmann's constant must be real.");

	auto const x = std::exp((E - mu) / (kb * t));

	if (std::isinf(x)) return 0.0;
	if (x < std::numeric_limits<decltype(x)>::epsilon()) return 1.0;
	return 1.0 / (x + 1.0); 
}


///////////////////////////////////////////////////////////////////////////////
/// \brief Defines tools to calculate the \f$ G(\omega) \f$ matrix.
///////////////////////////////////////////////////////////////////////////////
namespace g_function {

///////////////////////////////////////////////////////////////////////////////
/// \brief Computes \f$G_{i, j}(\omega)=\frac{f_i-f_j}{E_i-E_j-\omega}\f$.

/// \tparam    _F      Represents a real/complex field.
/// \tparam    _R      Represents a real field.
/// \tparam    _Number Just some abstract field.
/// \param     i       Row of the matrix \f$G(\omega)\f$.
/// \param     j       Column of the matrix \f$G(\omega)\f$.
/// \param     E       Pointer into the array of energies. `_F` must be real,
///                    and size of \p E must be at least `max(i,j)`.
/// \param     f       Pointer into the array of occupational numbers. `_R`
///                    must be real, and size of \p f must be at least
///                    `max(i,j)`.
///
/// \return \f$G_{i,j}(\omega)\f$.
///////////////////////////////////////////////////////////////////////////////
template<class _F, class _R, class _Number>
auto at( std::size_t const i, std::size_t const j
       , _Number const omega
       , _F const* E
       , _R const* f ) noexcept
{
	static_assert(std::is_floating_point<_F>::value, "Energy must be real.");
	static_assert(std::is_floating_point<_R>::value, "Occupational numbers " 
		"must be real.");
	return (f[i] - f[j]) / (E[i] - E[j] - omega);
}


///////////////////////////////////////////////////////////////////////////////
/// \brief Computes \f$ G(\omega) \f$.

/// \f[ G_{i,j}(\omega) = \frac{f_i - f_j}{E_i - E_j - \omega}. \f]
///
/// \param omega    Frequency \f$\omega\f$ at which to calculate \f$ G \f$.
///                 It may be either real or complex.
/// \param E        Energies of the system: \f$1 \times N\f$ matrix (i.e. a 
///                 column vector). `_F` must be floating point.
/// \param cs       Constants map. This function requires the availability of
///                 \f$ T, \mu, k_\text{B}, \hbar \f$ to run. Checks are
///                 performed at runtime, i.e. this function <b>may throw</b>!
///                 `_R`, the type of constants in \p cs must also be floating
///                 point.
/// \param lg       The logger.
/// \return         \f$G(\omega)\f$.
/// \exception      May throw.
///////////////////////////////////////////////////////////////////////////////
template<class _Number, class _F, class _R, class _Logger>
auto make( _Number const omega
         , Matrix<_F> const& E
         , std::map<std::string, _R> const& cs 
         , _Logger & lg )
{
	static_assert(std::is_floating_point<_F>::value, "Energy must be real.");
	static_assert(std::is_floating_point<_R>::value, "Physical constants " 
		"such as chemical potential and temperature must be real.");
	using Real = decltype( fermi_dirac( std::declval<_F>()
	                                  , std::declval<_R>() 
	                                  , std::declval<_R>()
	                                  , std::declval<_R>() ));
	using Complex = decltype( at( std::declval<std::size_t>()
	                            , std::declval<std::size_t>()
	                            , std::declval<_Number>()
	                            , std::declval<_F const*>() 
	                            , std::declval<_R const*>() ));

	TCM_MEASURE( "g_function::make<" + boost::core::demangle(
		typeid(Complex).name()) + ">()" );
	LOG(lg, debug) << "Calculating G for omega = " << omega << "...";
	require(__PRETTY_FUNCTION__, cs, "temperature");
	require(__PRETTY_FUNCTION__, cs, "chemical-potential");
	require(__PRETTY_FUNCTION__, cs, "boltzmann-constant");
	assert( is_column(E) == 1 );
	const auto t  = cs.at("temperature");
	const auto mu = cs.at("chemical-potential");
	const auto kb = cs.at("boltzmann-constant");
	const auto N  = E.height();

	Matrix<Real> f{N, 1};
	std::transform( E.data(), E.data() + N, f.data()
	              , [t, mu, kb](auto Ei) noexcept
	                { return fermi_dirac(Ei, t, mu, kb); } );
	Matrix<Complex> G{N, N};
	for (std::size_t j = 0; j < N; ++j) {
		for (std::size_t i = 0; i < N; ++i) {
			G(i,j) = at(i, j, omega, E.data(), f.data());
		}
	}

	LOG(lg, debug) << "Successfully calculated G.";
	return G;
}


} // namespace g_function





///////////////////////////////////////////////////////////////////////////////
/// Defines tools to compute \f$ \chi(\omega) \f$ matrix.
///////////////////////////////////////////////////////////////////////////////
namespace chi_function {

///////////////////////////////////////////////////////////////////////////////
/// \brief Calculates \f$ \chi_{a,b}(\omega) \f$.

/// We compute \f$ \chi_{a,b}(\omega) \f$ as
/// \f[ \begin{aligned}
///         A &:= \psi_a \circ \psi_b^*, 
///         \text{ i.e. } A_i = \psi_{i,a}\cdot\psi{i,b}^* \\ \text{}
///         \chi_{a,b}(\omega) &= A^\dagger \times G(\omega) \times A. 
///     \end{aligned}
/// \f]
/// 
/// Hadamard product \f$ \circ \f$ is computed using `std::transform`, which
/// can easily be parallelized using GNU parallel mode of the `stdlibc++`.
/// Matrix-vector products are calculated by first calling `?GEMV` and then
/// `?DOTC`. Intel MKL has a highly parallel implementation of both of this
/// functions.
///
/// \param a    Row of the matrix.
/// \param b    Column of the matrix.
/// \param Psi  Eigenstates of the system.
/// \param G    G function calculated by calling g_function::make().
///
/// \returns    \f$ \chi_{a,b}(\omega) \f$.
/// \exception  May throw if memory allocations or LAPACK operations fail.
///////////////////////////////////////////////////////////////////////////////
template <class _F, class _C>
auto at( std::size_t const a, std::size_t const b
       , Matrix<_F> const& Psi, Matrix<_C> const& G )
{
	TCM_MEASURE( "chi_function::at<" + boost::core::demangle(
		typeid(_F).name()) + ", " + boost::core::demangle(
		typeid(_C).name()) + ">()" );
	using Complex = std::common_type_t<_F, _C>;

	const auto N = Psi.height();
	Matrix<Complex>    A{N, 1};
	Matrix<Complex> temp{N, 1};

	std::transform( Psi.cbegin_row(a), Psi.cend_row(a)
	              , Psi.cbegin_row(b)
	              , A.data()
	              , [](auto x, auto y) { return x * std::conj(y); } );
	blas::gemv( blas::Operator::T
	          , Complex{1}, G, A
	          , Complex{0}, temp );
	return Complex{2} * blas::dot(A, temp);
}

namespace {
// For the case that eigenstates are actually real.
template<class _Number, class _F, class _R, class _Logger>
auto make_impl( _Number const omega
              , Matrix<_F> const& E
              , Matrix<_F> const& Psi
              , std::map<std::string, _R> const& cs
              , _Logger & lg )
{
	static_assert(std::is_floating_point<_F>::value, "Energy must be real.");
	static_assert(std::is_floating_point<_R>::value, "Physical constants " 
		"such as chemical potential and temperature must be real.");
	TCM_MEASURE( "chi_function::make_impl<" + boost::core::demangle(
		typeid(_F).name()) + ">()" );
	auto const N = E.height();
	auto const G = g_function::make(omega, E, cs, lg);

	using T = decltype( at( std::declval<std::size_t>()
	                      , std::declval<std::size_t>()
	                      , std::declval<Matrix<_F>>()
	                      , std::declval<decltype(G)>() ));
	Matrix<T> Chi{N, N};

	auto const _total_points_ = static_cast<double>(N * (N - 1) / 2);
	auto _start_ = std::chrono::system_clock::now();
	// fill the diagonal
	for (std::size_t i = 0; i < N; ++i) {
			if ((_start_ - std::chrono::system_clock::now()) >
				std::chrono::minutes{5}) {
				
				LOG(lg, info) << "at " << std::round((i + 1) / _total_points_)
				              << "% ..." << std::flush;
				_start_ = std::chrono::system_clock::now();
			}
			Chi(i, i) = at(i, i, Psi, G);
	}
	// calculate upper triangle
	for (std::size_t j = 0; j < N; ++j) {
		for (std::size_t i = 0; i < j; ++i) {
			if ((_start_ - std::chrono::system_clock::now()) >
				std::chrono::minutes{5}) {
				
				LOG(lg, info) << "at " << std::round((N + (i + 1)* (j + 1)) 
					/ _total_points_)
				              << "% ..." << std::flush;
				_start_ = std::chrono::system_clock::now();
			}
			Chi(i, j) = at(i, j, Psi, G);
			Chi(j, i) = Chi(i, j);
		}
	}
	return Chi;
}

// For the case that eigenstates are complex.
template<class _Number, class _F, class _R, class _Logger>
auto make_impl( _Number const omega
              , Matrix<_F> const& E
              , Matrix<std::complex<_F>> const& Psi
              , std::map<std::string, _R> const& cs
              , _Logger & lg )
{
	static_assert(std::is_floating_point<_F>::value, "Energy must be real.");
	static_assert(std::is_floating_point<_R>::value, "Physical constants " 
		"such as chemical potential and temperature must be real.");
	TCM_MEASURE( "chi_function::make_impl<" + boost::core::demangle(
		typeid(std::complex<_F>).name()) + ">()" );
	auto const N = E.height();
	auto const G = g_function::make(omega, E, cs, lg);

	using T = decltype( at( std::declval<std::size_t>()
	                      , std::declval<std::size_t>()
	                      , std::declval<Matrix<std::complex<_F>>>()
	                      , std::declval<decltype(G)>() ));
	Matrix<T> Chi{N, N};

	auto const _total_points_ = static_cast<double>(N * N);
	auto _start_ = std::chrono::system_clock::now();
	for (std::size_t j = 0; j < N; ++j) {
		for (std::size_t i = 0; i < N; ++i) {
			if ((_start_ - std::chrono::system_clock::now()) >
				std::chrono::minutes{5}) {
				
				LOG(lg, info) << "at " << std::round(((i + 1) * (j + 1)) / _total_points_)
				              << "% ..." << std::flush;
				_start_ = std::chrono::system_clock::now();
			}
			Chi(i, j) = at(i, j, Psi, G);
		}
	}
	return Chi;
}
} // end unnamed namespace


///////////////////////////////////////////////////////////////////////////////
/// \brief Calculates \f$ \chi(\omega) \f$.

/// First, calculates \f$ G(\omega) \f$ by calling tcm::g_function::make and
/// then computes each matrix element by calling at().
/// \tparam _Number   Some abstract field.
/// \tparam _F        Complex or real field.
/// \tparam _C        Complex field.
/// \tparam _F        Real field.
/// \tparam _Logger   Type of the logger. 
/// \param omega      Frequency \f$ \omega \f$ at which to calculate 
///                   \f$\chi\f$.
/// \param E          Eigenenergies of the system.
/// \param Psi        Eigenstates of the system.
/// \param cs         Constants. 
/// \param lg         Logger object.
///
/// \returns \f$ \chi(\omega) \f$ as a Matrix<_C>.
/// \exception May throw.
///////////////////////////////////////////////////////////////////////////////
template<class _Number, class _F, class _C, class _R, class _Logger>
auto make( _Number const omega
         , Matrix<_F> const& E
         , Matrix<_C> const& Psi
         , std::map<std::string, _R> const& cs
         , _Logger & lg )
{
	TCM_MEASURE( "chi_function::make<" + boost::core::demangle(
		typeid(_C).name()) + ">()" );
	LOG(lg, debug) << "Calculating chi for omega = " << omega << "...";

	const auto N = E.height();
	assert( is_column(E) );
	assert( is_square(Psi) );
	assert( N == Psi.height() );

	auto const Chi = make_impl(omega, E, Psi, cs, lg);
	LOG(lg, debug) << "Successfully calculating chi.";
	return Chi;
}


} // namespace chi_function



///////////////////////////////////////////////////////////////////////////////
/// \brief Defines tools to calculate the Coulomb interaction potential.
///////////////////////////////////////////////////////////////////////////////
namespace coulomb {

namespace {
///////////////////////////////////////////////////////////////////////////////
/// \brief Computes the distance between two points in \f$ F^3 \f$.
///////////////////////////////////////////////////////////////////////////////
template<class _F>
auto distance(std::array<_F, 3> const& v, std::array<_F, 3> const& w)
{
	return std::sqrt( std::pow(std::abs(v[0] - w[0]), 2)
	                + std::pow(std::abs(v[1] - w[1]), 2)
	                + std::pow(std::abs(v[2] - w[2]), 2)
	                );
}
} // unnamed namespace


///////////////////////////////////////////////////////////////////////////////
/// \brief Calculates \f$ V_{i,j} \f$.

/// Potential is calculated as following
/// \f[ V_{i,j} = \left\{ 
///     \begin{aligned}
///         &\frac{e}{4\pi\varepsilon_0|\mathbf{r_i} - \mathbf{r_j}|}
///              , &\text{if } i \neq j, \\ \text{}
///         &V_0, &\text{if } i = j.
///     \end{aligned}
///     \right.
/// \f]
/// 
/// \tparam _F        Real field.
/// \tparam _R        Another real field.
/// \param i          Index of the first atom, i.e. row of the matrix.
/// \param j          Index of the second atom,.i.e. column of the matrix.
/// \param positions  Positions of the atoms, in __meters__.
/// \param e          Elementary charge, in Coulombs.
/// \param pi         \f$\pi\f$.
/// \param eps0       \f$\varepsilon_0\f$.
/// \param v0         Self-interaction potential \f$V_0\f$, in eV.
///
/// \returns Potential in eV.
/// \exceptions Should not throw.
///////////////////////////////////////////////////////////////////////////////
template<class _F, class _R>
auto at( std::size_t const i, std::size_t const j
       , std::vector<std::array<_F, 3>> const& positions
       , _R const e
       , _R const pi
       , _R const eps0
       , _R const v0 ) noexcept
{
	return (i == j)
		? v0
		: e / (_R{4.0} * pi * eps0 * distance(positions[i], positions[j]));
}



///////////////////////////////////////////////////////////////////////////////
/// \brief Construct an Elemental-like matrix \f$ V \f$, representing the 
/// Coulomb potential.

/// \tparam _T       Element type of the output Matrix.
/// \tparam _F       Real field.
/// \tparam _R       Another real field.
/// \tparam _Logger  The logger type.
/// \param positions Array with positions of atoms.
/// \param cs        Constants. `elementary-charge`, `vacuum-permittivity`, 
///                  `pi` and `self-interaction-potential` are needed.
/// \param lg        The logger.
///
/// \return Matrix describing the Coulomb interaction.
/// \exception May throw.
///////////////////////////////////////////////////////////////////////////////
template<class _T, class _F, class _R, class _Logger>
auto make( std::vector<std::array<_F, 3>> const& positions
         , std::map<std::string, _R> const& cs 
         , _Logger & lg ) -> Matrix<_T>
{
	TCM_MEASURE( "coulomb::make<" + boost::core::demangle(
		typeid(_T).name()) + ">()" );
	LOG(lg, debug) << "Calculating V...";

	require(__PRETTY_FUNCTION__, cs, "elementary-charge");
	require(__PRETTY_FUNCTION__, cs, "pi");
	require(__PRETTY_FUNCTION__, cs, "vacuum-permittivity");
	require(__PRETTY_FUNCTION__, cs, "self-interaction-potential");

	auto const N         = positions.size();
	auto const e         = cs.at("elementary-charge");
	auto const pi        = cs.at("pi");
	auto const eps0      = cs.at("vacuum-permittivity");
	auto const v0        = cs.at("self-interaction-potential");

	auto const V = build_matrix
		( N, N
		, [&positions, e, pi, eps0, v0] (auto i, auto j)
		  { return boost::numeric_cast<_T>(
		        at(i, j, positions, e, pi, eps0, v0)
			);
		  }
		);

	LOG(lg, debug) << "Successfully calculated V.";
	return V;
}


} // namespace coulomb







///////////////////////////////////////////////////////////////////////////////
/// \brief Defines tools to compute the dielectric function.
///////////////////////////////////////////////////////////////////////////////
namespace dielectric_function {


///////////////////////////////////////////////////////////////////////////////
/// \brief Calculates the dielectric function matrix \f$\epsilon(\omega)\f$.
///////////////////////////////////////////////////////////////////////////////
template< class _Number, class _F, class _C, class _R, class _T, class _Logger>
auto make( _Number const omega
         , Matrix<_F> const& E
         , Matrix<_C> const& Psi
         , Matrix<_T> const& V
         , std::map<std::string, _R> const& cs 
         , _Logger & lg )
{
	TCM_MEASURE( "dielectric_function::make<" + boost::core::demangle(
		typeid(_C).name()) + ">()" );
	LOG(lg, debug) << "Calculating epsilon for omega = " << omega << "...";

	const auto N = E.height();
	assert( is_column(E) );
	assert( is_square(Psi) );
	assert( is_square(V) );
	assert( Psi.height() == N );
	assert( V.height() == N );

	auto const Chi = chi_function::make(omega, E, Psi, cs, lg);

	static_assert(std::is_same<_T, typename decltype(Chi)::value_type>::value, "");
	Matrix<_T> epsilon{N, N};
	for (std::size_t j = 0; j < N; ++j) {
		for (std::size_t i = 0; i < N; ++i) {
			epsilon(i, j) = (i == j) ? 1.0 : 0.0;
		}
	}

	blas::gemm( blas::Operator::None, blas::Operator::None
	          , _T{-1.0}, V, Chi
	          , _T{ 1.0}, epsilon );

	LOG(lg, debug) << "Successfully calculating epsilon.";
	return epsilon;
}



} // namespace dielectric function



} // namespace tcm


#endif // TCM_DIELECTRIC_FUNCTION_HPP
