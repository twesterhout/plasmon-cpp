#ifndef TCM_CONSTANTS_HPP
#define TCM_CONSTANTS_HPP

#include <functional>
#include <type_traits>
#include <map>
#include <unordered_map>
#include <cmath>

#include <boost/program_options.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>


///////////////////////////////////////////////////////////////////////////////
/// \file constants.hpp
/// \brief Defines the 'default' constants map.
///
/// \detail We implement some useful functions to manipulate _physical_
///         constants. These constants are represented as a symbol table
///         mapping #std::string to some numerical type _R representing real
///         field.
///
/// For the calculations, we're interested in the following constants:
/// |              Key             |     Notation     | Value(default) | Dimension |
/// | ---------------------------- | ---------------- | -------------- | --------- |
/// | `pi`                         | \f$\pi\f$        | 3.14159...     | 1         |
/// | `boltzmann-constant`         | \f$k_\text{B}\f$ | 8.61733...E-5  | eV        |
/// | `elementary-charge`          | \f$e\f$          | 1.60217...E-19 | C         |
/// | `planck-constant`            | \f$\hbar\f$      | 6.58212...E-16 | eV        |
/// | `self-interaction-potential` | \f$V_0\f$        | 15.78          | eV        |
/// | `temperature`                | \f$T\f$          | 300.0          | K         |
/// | `vacuum-permittivity`        | \f$\epsilon_0\f$ | 8.85419...E-12 | F/m       |
/// A table with these values can be easily obtained using the 
/// #default_constants() function.
///
/// For different simulations different values of, for example, temperature
/// may be desired. To account for this case, we implement two additional
/// functions: 
/// * #init_constants_options() simplifies the creation of 
///   boost::program_options::options_description for the user to be able
///   to specify constants.
/// * #load_constants() allows the creation of constants table from
///   boost::program_options::variables_map, i.e. from command line arguments.
///
/// Here is an example how these functions can be used together:
/// \include src/constants_example.cpp
///////////////////////////////////////////////////////////////////////////////


namespace tcm {

//////////////////////////////////////////////////////////////////////////////
/// \brief Returns a table with default values for relevalnt physical 
///        constants.

/// This function has the following behavior:
/// \snippet include/constants.hpp Default constants
///
/// \todo make this function return std::unordered_map in place of std::map.
/// Or even better, make the container type a template parameter.
//////////////////////////////////////////////////////////////////////////////
template <class _R>
auto default_constants() -> std::map<std::string, _R>
{
//! [Default constants]
	std::map<std::string, _R> constants =
		{ {"pi",                          M_PI}
		, {"boltzmann-constant",          8.6173303E-5}
		, {"chemical-potential",          0.4}
		, {"elementary-charge",           1.6021766208E-19}
		, {"planck-constant",             6.582119514E-16}
		, {"self-interaction-potential",  15.78}
		, {"temperature",                 300.0}
		, {"vacuum-permittivity",         8.854187817E-12}
		, {"tau",                         6.0E-3}
		};
	return constants;
//! [Default constants]
}


///////////////////////////////////////////////////////////////////////////////
/// \brief Simplifies the creation of #options_descriptions concerning
///        constants.

/// Creates a #boost::program_options::options_description called 
/// _Simulation Constants_ with options of the form "in.constants.KEY" where
/// KEYs are the names of the constants such as `pi` or `planck-constant`.
/// Each of these options has a default value provided by the
/// #default_constants() function. Template parameter `_R` specifies the type
/// of the constants. It is used to create the default constants table and to
/// extract values from `any` where program_options stores the values.
///////////////////////////////////////////////////////////////////////////////
template <class _R>
auto init_constants_options() -> boost::program_options::options_description
{
	using namespace boost::program_options;
	options_description description{"Simulation Constants"};
	auto const cs = tcm::default_constants<_R>();
	description.add_options()
		( "in.constants.pi"
		, value<_R>()->default_value(cs.at("pi"))
		, "PI." )
		( "in.constants.boltzmann-constant"
		, value<_R>()->default_value(cs.at("boltzmann-constant"))
		, "Boltzmann constant." )
		( "in.constants.chemical-potential"
		, value<_R>()->default_value(cs.at("chemical-potential"))
		, "Chemical potential." )
		( "in.constants.elementary-charge"
		, value<_R>()->default_value(cs.at("elementary-charge"))
		, "Elementary charge." )
		( "in.constants.planck-constant"
		, value<_R>()->default_value(cs.at("planck-constant"))
		, "Planck constant" )
		( "in.constants.self-interaction-potential"
		, value<_R>()->default_value(cs.at("self-interaction-potential"))
		, "Self interaction Coulomb potential." )
		( "in.constants.temperature"
		, value<_R>()->default_value(cs.at("temperature"))
		, "Temperature" )
		( "in.constants.vacuum-permittivity"
		, value<_R>()->default_value(cs.at("vacuum-permittivity"))
		, "Vacuum permittivity." )
		( "in.constants.tau"
		, value<_R>()->default_value(cs.at("tau"))
		, "Relaxation time tau." );
	return description;
}


///////////////////////////////////////////////////////////////////////////////
/// \brief Creates a constants table from a 
///        #boost::program_options::variables_map.

/// This function iterates over entries in \p vm, searches for all options
/// that are of the form "in.constants.*", converts the values from
/// boost::program_options::variable_value to _To and returns a table of 
/// type _Map consisting only of those options.
/// \snippet include/constants.hpp Load constants
///////////////////////////////////////////////////////////////////////////////
template< class _To
        , class _From = _To
        , class _Map = std::unordered_map<std::string, _To> 
        >
auto load_constants(boost::program_options::variables_map const& vm) -> _Map
{
	static_assert( std::is_same<typename _Map::key_type, std::string>::value
	             , "Wrong key_type: must be std::string." );
	static_assert( std::is_same<typename _Map::mapped_type, _To>::value
	             , "Wrong mapped_type: must be _To." );
//! [Load constants]
	auto const is_constant = [](auto const& p) { 
		return boost::algorithm::starts_with(p.first, "in.constants."); 
	};

	using value_type = boost::program_options::variables_map::value_type;
	auto const mk_constant = [](value_type const& p) {
		return std::make_pair(p.first, p.second.as<_From>());
	};

	auto rng = vm | boost::adaptors::filtered(std::cref(is_constant))
	              | boost::adaptors::transformed(std::cref(mk_constant));
	return _Map{rng.begin(), rng.end()};
//! [Load constants]
}


///////////////////////////////////////////////////////////////////////////////
/// \brief Checks whether \p key is in \p constants map and throws if it is
///        not.

/// \snippet include/constants.hpp Require constant
///////////////////////////////////////////////////////////////////////////////
template<class _Map>
auto require( std::string const& func_name
            , _Map const& constants_map
            , std::string const& key ) -> void
{
//! [Require constant]
	if (!constants_map.count(key)) {
		throw std::runtime_error{ "Constant `" + key + "` is required to run "
		                          "`" + func_name + "`!" };
	}
//! [Require constant]
}

} // namespace tcm

#endif // TCM_CONSTANTS_HPP
