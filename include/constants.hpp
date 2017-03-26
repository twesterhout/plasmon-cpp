#ifndef TCM_CONSTANTS_HPP
#define TCM_CONSTANTS_HPP

#include <functional>
#include <map>
#include <cmath>

#include <boost/program_options.hpp>


//! \file constants.hpp
//! \brief Defines the 'default' constants map.

namespace tcm {

//! \brief Defines some constants

//! * `pi` - \f$ \pi = 3.14159... \f$
//! * `boltzmann-constant` - \f$ k_\text{B} = 8.61733... \cdot 10^{-5} \text{ eV}\f$ 
//!    (see <a href="https://en.wikipedia.org/wiki/Boltzmann_constant">wikipedia</a>).
//! * `chemical-potential` - \f$ \mu = \epsilon_\text{F} \f$, i.e. Fermi energy
//!    in the article.
//! * `elementary-charge` - \f$ e = 1.60217... \text{ C}\f$ 
//!    (see <a href="https://en.wikipedia.org/wiki/Elementary_charge">wikipedia</a>).
//! * `planck-constant` - \f$ \hbar = 6.58212... \cdot 10^{-16} \text{ eV}\cdot\text{s}\f$ 
//!    (see <a href="https://en.wikipedia.org/wiki/Planck_constant">wikipedia</a>).
//! * `self-interaction-potential` - \f$ V_0 = 0.58 \text{ a.u.} \approx 15.78 \text{ eV}\f$,
//!    from the article
//! * `temperature` - \f$ T \f$
//! * `vacuum-permittivity` - \f$ \epsilon_0 = 8.85419... \cdot 10^{-12} \text{ F/m}\f$ 
//!    (see <a href="https://en.wikipedia.org/wiki/Vacuum_permittivity">wikipedia</a>).


template <class _R>
auto default_constants() -> std::map<std::string, _R>
{
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
}


template <class _R>
auto init_constants_options() -> boost::program_options::options_description
{
	using namespace boost::program_options;
	options_description description;
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


template<class _To, class _From = _To>
auto load_constants(boost::program_options::variables_map const& vm) -> std::map<std::string, _To>
{
	std::map<std::string, _To> cs = 
		{ {"pi",                          vm["in.constants.pi"].as<_From>()}
		, {"boltzmann-constant",          vm["in.constants.boltzmann-constant"].as<_From>()}
		, {"chemical-potential",          vm["in.constants.chemical-potential"].as<_From>()}
		, {"elementary-charge",           vm["in.constants.elementary-charge"].as<_From>()}
		, {"planck-constant",             vm["in.constants.planck-constant"].as<_From>()}
		, {"self-interaction-potential",  vm["in.constants.self-interaction-potential"].as<_From>()}
		, {"temperature",                 vm["in.constants.temperature"].as<_From>()}
		, {"vacuum-permittivity",         vm["in.constants.vacuum-permittivity"].as<_From>()}
		, {"tau",                         vm["in.constants.tau"].as<_From>()}
		};
	return cs;
}


template<class T>
auto require( std::string const& func_name
            , std::map<std::string, T> const& constants_map
            , std::string const& key
            ) -> void
{
	if(!constants_map.count(key)) {
		throw std::runtime_error
		{"Constant `" + key + "` is required to run "
		 "`" + func_name + "`!"};
	}
}


} // namespace tcm

#endif // TCM_CONSTANTS_HPP
