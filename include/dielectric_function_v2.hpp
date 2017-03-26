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


#include <config.hpp>
#include <logging.hpp>
#include <benchmark.hpp>

#include <constants.hpp>
#include <matrix.hpp>
#include <blas.hpp>




namespace tcm {


///////////////////////////////////////////////////////////////////////////////
/// \brief Computes \f$ f(E) \f$

/// Calculates the Fermi-Dirac distribution:
/// \f[ f(E) := \frac{1}{\exp(\frac{E - \mu}{k_{\text{B}}T}) + 1}. \f]
/// All parameters are templated to allow the use of floating point numbers 
/// with different precision. The result, however, is (implicitely) converted 
/// to \p _F, i.e. the type of \p E.
/// \tparam    _F      Represents a real/compex field.
/// \tparam    _R      Represents a real field.
/// \param     mu      \f$ \mu \f$ is the <em>chemical potential</em>
/// \param     t       \f$ T \f$ is the <em>temperature</em> 
/// \param     kb      \f$ k_\text{B} \f$ is the <em>Boltzmann constant</em>.
///
/// \return \f$ f(E) \f$.
/// \except Should not throw.
///////////////////////////////////////////////////////////////////////////////
template<class _F, class _R>
auto fermi_dirac( _F const E
                , _R const t, _R const mu, _R const kb ) noexcept -> _F
{
	auto const arg = (E - mu) / (kb * t);
	if(arg >   700.0) return 0.0;
	if(arg < - 700.0) return 1.0;
	return 1.0 / (std::exp(arg) + 1.0); 
}


///////////////////////////////////////////////////////////////////////////////
/// \brief Defines tools to calculate the \f$ G(\omega) \f$ matrix.
///////////////////////////////////////////////////////////////////////////////
namespace g_function {


namespace {
///////////////////////////////////////////////////////////////////////////////
/// \brief Computes \f$G_{i, j}(\omega)=\frac{f_i-f_j}{E_i-E_j-\omega}\f$.

/// \tparam    _F      Represents a real/complex field.
/// \tparam    _R      Represents a real field.
/// \tparam    _Number Just some abstract field.
/// \param     i       Row of the matrix \f$G(\omega)\f$.
/// \param     j       Column of the matrix \f$G(\omega)\f$.
/// \param     E       Pointer into the array of energies. `E[k]` must be
///                    equal to \f$ E_k \f$ for all `k`s.
/// \param     f       Pointer into the array of occupational numbers. `f[k]`
///                    must be equal to \f$ f(E_k) \f$ for all `k`s.
///
/// \return Element of the matrix \f$G(\omega)\f$ at row \p i and column \p j,
/// i.e. \f$G_{i,j}(\omega)\f$.
///////////////////////////////////////////////////////////////////////////////
template<class _F, class _R, class _Number>
auto at( std::size_t const i, std::size_t const j
       , _Number const omega
       , _F const* E
       , _R const* f ) noexcept
{
	return (f[i] - f[j]) / (E[i] - E[j] - omega);
}
} // unnamed namespace



///////////////////////////////////////////////////////////////////////////////
/// \brief Computes \f$ G(\omega) \f$.

/// \f[ G_{i,j}(\omega) = \frac{f_i - f_j}{E_i - E_j - \omega}. \f]
///
/// \tparam T       Element type of the output matrix. This allow for the
/// construction of complex matrices even if the result of calling at() is
/// real.
/// \tparam _Number Some abstract field.
/// \tparam _F      A real/complex field.
/// \tparam _R      A real field.
/// \tparam _Logger Type of the logger. Usually something like 
///                 `boost::log::sources::logger`.
/// \param omega    Frequency \f$\omega\f$ at which to calculate \f$ G \f$.
/// \param E        Energies of the system: \f$1 \times N\f$ matrix (i.e. a 
///                 column vector
/// \param cs       Constants map. This function requires the availability of
///                 \f$ T, \mu, k_\text{B}, \hbar \f$ to run. Checks are
///                 performed at runtime, i.e. this function <b>may throw</b>!
/// \param lg       The logger.
/// \return         \f$N \times N\f$ Matrix<T> \f$G(\omega)\f$.
/// \exception      May throw.
///////////////////////////////////////////////////////////////////////////////
template<class T, class _Number, class _F, class _R, class _Logger>
auto make( _Number const omega
         , Matrix<_F> const& E
         , std::map<std::string, _R> const& cs 
         , _Logger & lg ) -> Matrix<T>
{
	MEASURE;
	LOG(lg, debug) << "Calculating G for omega = " << omega << "...";

	require(__PRETTY_FUNCTION__, cs, "temperature");
	require(__PRETTY_FUNCTION__, cs, "chemical-potential");
	require(__PRETTY_FUNCTION__, cs, "boltzmann-constant");
	assert( is_column(E) == 1 );

	const auto t  = cs.at("temperature");
	const auto mu = cs.at("chemical-potential");
	const auto kb = cs.at("boltzmann-constant");
	const auto N  = E.height();

	Matrix<_F> f{N, 1};
	std::transform( E.data(), E.data() + N
	              , f.data()
	              , [t, mu, kb](auto Ei) noexcept
	                { return fermi_dirac(Ei, t, mu, kb); } 
                  );
	
	auto G = build_matrix
		( N, N
		, [omega, &E, &f](auto i, auto j)
		  { return T{ at(i, j, omega, E.data(), f.data()) }; }
		);

	LOG(lg, debug) << "Successfully calculated G.";
	return G;
}


} // namespace g_function





///////////////////////////////////////////////////////////////////////////////
/// Defines tools to compute \f$ \chi(\omega) \f$ matrix.
///////////////////////////////////////////////////////////////////////////////
namespace chi_function {

namespace {
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
/// \tparam _C  Complex field. Should be either `std::complex<float>` or
///             `std::comlex<double>`. Just `float` or `double` should work,
///             too. Other types are not supported by LAPACK.
/// \param a    Row of the matrix.
/// \param b    Column of the matrix.
/// \param Psi  Eigenstates of the system.
/// \param G    G function calculated by calling g_function::make().
///
/// \returns    \f$ \chi_{a,b}(\omega) \f$.
/// \exception  May throw if memory allocations or LAPACK operations fail.
///////////////////////////////////////////////////////////////////////////////
template<class _C>
auto at( std::size_t const a, std::size_t const b
       , Matrix<_C> const& Psi
       , Matrix<_C> const& G ) -> _C
{
	MEASURE;

	const auto N = Psi.height();
	Matrix<_C>    A{N, 1};
	Matrix<_C> temp{N, 1};

	// A := (Psi_a o Psi_b*)
	std::transform( Psi.cbegin_row(a), Psi.cend_row(a)
	              , Psi.cbegin_row(b)
	              , A.data()
	              , [](auto x, auto y) { return x * std::conj(y); }
	              );
	// temp := G . A
	blas::gemv( blas::Operator::T
	          , _C{1.0}, G, A
	          , _C{0.0}, temp
	          );
	// 2 A* . temp = 2 A* . G . A
	return _C{2.0} * blas::dot(A, temp);
}
} // unnamed namespace


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
         , _Logger & lg ) -> Matrix<_C>
{
	MEASURE;
	LOG(lg, debug) << "Calculating chi for omega = " << omega << "...";

	const auto N = E.height();
	assert( is_column(E) );
	assert( is_square(Psi) );
	assert( N == Psi.height() );

	auto const G = g_function::make<_C>(omega, E, cs, lg);
	auto const Chi = build_matrix
		( N, N
		, [&Psi, &G] (auto a, auto b) { return at(a, b, Psi, G); }
	    );

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


namespace {
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
} // unnamed namespace



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
	MEASURE;
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
template< class _Number, class _F, class _C, class _R, class _Logger>
auto make( _Number const omega
         , Matrix<_F> const& E
         , Matrix<_C> const& Psi
         , Matrix<_C> const& V
         , std::map<std::string, _R> const& cs 
         , _Logger & lg )
{
	MEASURE;
	LOG(lg, debug) << "Calculating epsilon for omega = " << omega << "...";

	const auto N = E.height();
	assert( is_column(E) );
	assert( is_square(Psi) );
	assert( is_square(V) );
	assert( Psi.height() == N );
	assert( V.height() == N );

	auto const Chi = chi_function::make(omega, E, Psi, cs, lg);
	Matrix<_C> epsilon{N, N};

	for(std::size_t i = 0; i < N; ++i)
		epsilon(i, i) = 1.0;

	blas::gemm( blas::Operator::None, blas::Operator::None
	          , _C{-1.0}, V, Chi
	          , _C{ 1.0}, epsilon );


	LOG(lg, debug) << "Successfully calculating epsilon.";
	return epsilon;
}



} // namespace dielectric function



} // namespace tcm


#endif // TCM_DIELECTRIC_FUNCTION_HPP
